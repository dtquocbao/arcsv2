"""ML2 sample selection, coverage, blacklist fallback, and workflow markers."""

from __future__ import annotations

import gzip
import hashlib
import shutil
import subprocess
import urllib.error
from pathlib import Path

import pytest
from scripts.fetch_filter_beds import fetch_blacklist, write_blacklist_source
from scripts.select_ml2_samples import HPRC_COMMIT, choose, main

from arcsv2.coverage import (
    chrom_median_windows,
    coverage_rows,
    inside_windows,
    plan_regions,
    rows_from_totals,
)

ROOT = Path(__file__).resolve().parents[1]
SUPERS = ["AFR", "AMR", "EAS", "EUR", "SAS"]
POPS = {"AFR": "ACB", "AMR": "CLM", "EAS": "CHB", "EUR": "FIN", "SAS": "BEB"}


def test_long_span_uses_twenty_even_windows():
    windows = inside_windows(0, 100_000)
    assert len(windows) == 20
    assert windows[0] == (0, 1000)
    assert windows[-1] == (99_000, 100_000)
    assert all(end - start == 1000 for start, end in windows)
    assert inside_windows(1000, 2000) == [(1000, 2000)]


def test_zero_flank_depth_is_na():
    calls = [{"chrom": "chrSim", "minbp": "1000", "maxbp": "2000", "id": "c1"}]
    lengths = {"chrSim": 4000}
    _regions, meta = plan_regions(calls, lengths, {"chrSim": []}, seed=1)
    totals = {}
    for name, (kind, _key, length) in meta.items():
        totals[name] = 0 if kind == "flank" else 4 * length
    row = rows_from_totals(calls, meta, totals)[0]
    assert row["id"] == "c1"
    assert row["cov_inside"] == "4.000000"
    assert row["cov_flank"] == "0.000000"
    assert row["coverage_ratio"] == "NA"
    assert float(row["cov_chrom_median"]) == pytest.approx(4.0)
    assert float(row["coverage_ratio_chrom"]) == pytest.approx(1.0)


def test_chrom_median_windows_skip_gaps_and_repeat():
    gaps = [(0, 3000)]
    first = chrom_median_windows("chrSim", 5000, gaps, seed=7, n=30)
    second = chrom_median_windows("chrSim", 5000, gaps, seed=7, n=30)
    assert first == second
    assert len(first) == 30
    for start, end in first:
        assert end <= 5000
        assert not (start < 3000 and end > 0)


def _uniform_bam(tmp_path: Path) -> tuple[Path, Path, Path]:
    """Depth 4 across chrSim, plus MAPQ-10 reads that must be ignored."""
    sam = tmp_path / "reads.sam"
    lines = ["@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:chrSim\tLN:4000"]
    seq = "A" * 100
    qual = "I" * 100
    n = 0
    for start in range(1, 4001, 100):
        for _copy in range(4):
            n += 1
            lines.append(f"r{n}\t0\tchrSim\t{start}\t60\t100M\t*\t0\t0\t{seq}\t{qual}")
    for _extra in range(10):
        n += 1
        lines.append(f"r{n}\t0\tchrSim\t1001\t10\t100M\t*\t0\t0\t{seq}\t{qual}")
    sam.write_text("\n".join(lines) + "\n", encoding="utf-8")
    bam = tmp_path / "reads.bam"
    subprocess.run(["samtools", "view", "-bS", "-o", str(bam), str(sam)], check=True)
    subprocess.run(["samtools", "sort", "-o", str(bam), str(bam)], check=True)
    subprocess.run(["samtools", "index", str(bam)], check=True)
    calls = tmp_path / "arcsv_out.tab"
    calls.write_text(
        "chrom\tminbp\tmaxbp\tid\nchrSim\t1000\t2000\tc1\n",
        encoding="utf-8",
    )
    gaps = tmp_path / "gaps.bed"
    gaps.write_text("chrSim\t3500\t4000\n", encoding="utf-8")
    return bam, calls, gaps


@pytest.mark.skipif(shutil.which("samtools") is None, reason="samtools is not on PATH")
def test_coverage_on_label_smoke_depths(tmp_path: Path):
    bam, calls, gaps = _uniform_bam(tmp_path)
    rows = coverage_rows(
        [{"chrom": "chrSim", "minbp": "1000", "maxbp": "2000", "id": "c1"}],
        str(bam),
        gaps=gaps,
        lengths={"chrSim": 4000},
    )
    row = rows[0]
    assert float(row["cov_inside"]) == pytest.approx(4.0)
    assert float(row["cov_flank"]) == pytest.approx(4.0)
    assert float(row["coverage_ratio"]) == pytest.approx(1.0)
    assert float(row["cov_chrom_median"]) == pytest.approx(4.0)
    assert calls.name == "arcsv_out.tab"


def test_label_smk_temp_markers_and_unchanged_call_snakefile():
    label = (ROOT / "workflows" / "label.smk").read_text(encoding="utf-8")
    fetch = label.split("rule fetch_region:", 1)[1].split("\nrule ", 1)[0]
    assert "temp(" in fetch
    assert ".bam" in fetch and ".bai" in fetch
    index = label.split("rule index_haps:", 1)[1].split("\nrule ", 1)[0]
    assert "temp(" in index and ".mmi" in index
    coverage = label.split("rule coverage:", 1)[1].split("\nrule ", 1)[0]
    assert "coverage.tsv" in coverage
    assert "PYTHONPATH=python" in coverage
    call = label.split("rule call:", 1)[1].split("\nrule ", 1)[0]
    assert "PYTHONPATH" not in call
    merge = label.split("rule merge:", 1)[1].split("\nrule ", 1)[0]
    assert "coverage.tsv" in merge
    assert "{input.tabs:q}" in merge
    call_only = (ROOT / "workflows" / "Snakefile").read_text(encoding="utf-8")
    assert "temp(" not in call_only
    sbatch = (ROOT / "scripts" / "hpc" / "run_label.sbatch").read_text(encoding="utf-8")
    assert "--partition=long" in sbatch
    assert "--cpus-per-task=1" in sbatch
    assert "--mem=4G" in sbatch
    assert "--time=7-00:00:00" in sbatch
    assert "run_label.sh" in sbatch
    assert "--rerun-incomplete" in sbatch


def _catalog(data: Path) -> None:
    catalog = data / "catalog"
    ped = catalog / "pedigree"
    ped.mkdir(parents=True)
    kg_header = (
        "#ENA_FILE_PATH\tMD5SUM\tRUN_ID\tSTUDY_ID\tSTUDY_NAME\tCENTER_NAME\t"
        "SUBMISSION_ID\tSUBMISSION_DATE\tSAMPLE_ID\tSAMPLE_NAME\tPOPULATION"
    )
    unrelated = [kg_header]
    related = [kg_header]
    hprc = ["sample_id,haplotype,assembly,assembly_md5,assembly_fai,assembly_gzi"]
    ped_lines = []

    def add(sample: str, pop: str, family: str, related_sample: bool = False) -> None:
        run = "ERR" + sample
        bucket = (
            "s3://1000genomes/1000G_2504_high_coverage/additional_698_related/data/"
            if related_sample
            else "s3://1000genomes/1000G_2504_high_coverage/data/"
        )
        line = (
            f"{bucket}{run}/{sample}.final.cram\tabc\t{run}\tSTUDY\tSTUDY\tNYGC\t"
            f"SUB\t2020-01-01\t{sample}\t{sample}\t{pop}"
        )
        (related if related_sample else unrelated).append(line)
        for hap in ("1", "2"):
            hprc.append(
                f"{sample},{hap},s3://human-pangenomics/{sample}.hap{hap}.fa.gz,"
                f"s3://human-pangenomics/{sample}.hap{hap}.fa.gz.md5,"
                f"s3://human-pangenomics/{sample}.hap{hap}.fa.gz.fai,"
                f"s3://human-pangenomics/{sample}.hap{hap}.fa.gz.gzi"
            )
        ped_lines.append(f"{family} {sample} 0 0 1 0")

    add("ML0A", "ACB", "FAM_ML0")
    add("REL1", "ACB", "FAM_ML0", related_sample=True)
    add("DUP_A", "ACB", "FAM_DUP")
    add("DUP_B", "ACB", "FAM_DUP")
    for i in range(30):
        superpop = SUPERS[i % 5]
        add(f"S{i:02d}", POPS[superpop], f"FAM{i:02d}")

    (catalog / "1000G_2504_high_coverage.sequence.index").write_text(
        "\n".join(unrelated) + "\n", encoding="utf-8"
    )
    (catalog / "1000G_698_related_high_coverage.sequence.index").write_text(
        "\n".join(related) + "\n", encoding="utf-8"
    )
    commit12 = HPRC_COMMIT[:12]
    (catalog / f"assemblies_release2_v1.0.{commit12}.csv").write_text(
        "\n".join(hprc) + "\n", encoding="utf-8"
    )
    (ped / "20130606_g1k.ped").write_text("\n".join(ped_lines) + "\n", encoding="utf-8")


def _ml0(path: Path) -> None:
    path.write_text(
        "sample\tpopulation\tsuperpopulation\tsex\tfamily_id\t"
        "cram_url\tcrai_url\thap1_url\thap2_url\t"
        "cram_checksum\tcrai_checksum\thap1_checksum\thap2_checksum\n"
        "ML0A\tACB\tAFR\tfemale\tFAM_ML0\thttps://example/a.cram\thttps://example/a.cram.crai\t"
        "https://example/a.fa\thttps://example/b.fa\tmd5:abc\tnone\tmd5:a\tmd5:b\n",
        encoding="utf-8",
    )


def _samples(path: Path) -> list[str]:
    ids = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#") or line.startswith("sample\t"):
            continue
        ids.append(line.split("\t")[0])
    return ids


def test_selection_is_deterministic_and_drops_families(tmp_path: Path):
    data = tmp_path / "data"
    _catalog(data)
    ml0 = tmp_path / "ml0.tsv"
    _ml0(ml0)
    outs = []
    for name in ("a", "b"):
        out = tmp_path / name
        code = main(
            [
                "--data",
                str(data),
                "--ml0",
                str(ml0),
                "--out-dir",
                str(out),
                "--scratch",
                "/scratch/user/dtquocbao",
                "--seed",
                "204178949",
            ]
        )
        assert code == 0
        outs.append(out)
    files_a = sorted(p.name for p in outs[0].iterdir())
    files_b = sorted(p.name for p in outs[1].iterdir())
    assert files_a == files_b
    for name in files_a:
        assert (outs[0] / name).read_bytes() == (outs[1] / name).read_bytes()

    chosen = []
    for path in outs[0].glob("samples.ml2.*.tsv"):
        chosen.extend(_samples(path))
    assert "ML0A" not in chosen
    assert "REL1" not in chosen
    assert "DUP_A" in chosen
    assert "DUP_B" not in chosen
    assert len(_samples(outs[0] / "samples.ml2.test.tsv")) == 20
    assert len(_samples(outs[0] / "samples.ml2.batch01.tsv")) == 10
    test_rows = [
        line.split("\t")
        for line in (outs[0] / "samples.ml2.test.tsv").read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#") and not line.startswith("sample\t")
    ]
    header = [
        line.split("\t")
        for line in (outs[0] / "samples.ml2.test.tsv").read_text(encoding="utf-8").splitlines()
        if line.startswith("sample\t")
    ][0]
    super_i = header.index("superpopulation")
    split_i = header.index("split")
    counts = {superpop: 0 for superpop in SUPERS}
    for row in test_rows:
        assert row[split_i] == "test"
        counts[row[super_i]] += 1
    assert max(counts.values()) - min(counts.values()) <= 1

    yaml_text = (outs[0] / "label.ml2.batch01.yaml").read_text(encoding="utf-8")
    assert "call_min_per_mb: 2.0" in yaml_text
    assert "mem_mb_floor: 28000" in yaml_text
    assert "threads: 4" in yaml_text
    assert "label_mem_mb: 30000" in yaml_text
    assert "data_root: /scratch/user/dtquocbao/arcsv2/data" in yaml_text
    assert "https://1000genomes.s3.amazonaws.com/" in yaml_text
    assert "- chrX" in yaml_text
    assert "streamed" in yaml_text


def test_selection_caps_at_120():
    rows = []
    for i in range(200):
        superpop = SUPERS[i % 5]
        rows.append(
            {
                "sample": f"N{i:03d}",
                "population": POPS[superpop],
                "superpopulation": superpop,
                "sex": "",
                "family_id": f"F{i:03d}",
                "cram_url": "https://example/x.cram",
            }
        )
    test, train, dropped, excluded = choose(rows, ml0_ids=set(), ml0_families=set(), seed=3)
    assert excluded == 0
    assert dropped == []
    assert len(test) == 20
    assert len(train) == 100
    again, train2, _, _ = choose(rows, ml0_ids=set(), ml0_families=set(), seed=3)
    assert [row["sample"] for row in test] == [row["sample"] for row in again]
    assert [row["sample"] for row in train] == [row["sample"] for row in train2]


def test_blacklist_keeps_10x_when_it_downloads(tmp_path: Path):
    def download(url, dest: Path):
        assert "10xgenomics" in url
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(b"chr1\t0\t5\n")
        return dest

    path, source, digest = fetch_blacklist(tmp_path, download)
    assert source == "10x"
    assert digest == hashlib.sha256(b"chr1\t0\t5\n").hexdigest()
    assert path.read_bytes() == b"chr1\t0\t5\n"


def test_blacklist_falls_back_to_encode_v2(tmp_path: Path, monkeypatch):
    payload_bed = b"chr1\t0\t10\tLow Mappability\n"
    buf = gzip.compress(payload_bed)
    monkeypatch.setattr(
        "scripts.fetch_filter_beds.ENCODE_BLACKLIST_SHA256",
        hashlib.sha256(buf).hexdigest(),
    )
    seen: list[str] = []

    def download(url, dest: Path):
        seen.append(url)
        dest.parent.mkdir(parents=True, exist_ok=True)
        if "10xgenomics" in url:
            raise urllib.error.URLError("403")
        dest.write_bytes(buf)
        return dest

    path, source, digest = fetch_blacklist(tmp_path, download)
    assert source == "encode_v2"
    assert seen[0].startswith("https://cf.10xgenomics.com/")
    assert "Boyle-Lab/Blacklist" in seen[1]
    assert path.read_bytes() == payload_bed
    write_blacklist_source(tmp_path, source, seen[1], digest)
    note = (tmp_path / "blacklist.source.txt").read_text(encoding="utf-8")
    assert "source=encode_v2" in note
    assert "not the 10x SV blacklist" in note


def test_blacklist_rejects_a_bad_encode_hash(tmp_path: Path, monkeypatch):
    monkeypatch.setattr("scripts.fetch_filter_beds.ENCODE_BLACKLIST_SHA256", "0" * 64)

    def download(url, dest: Path):
        dest.parent.mkdir(parents=True, exist_ok=True)
        if "10xgenomics" in url:
            raise urllib.error.URLError("403")
        dest.write_bytes(gzip.compress(b"chr1\t0\t1\n"))
        return dest

    with pytest.raises(SystemExit, match="sha256"):
        fetch_blacklist(tmp_path, download)
