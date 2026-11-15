"""Haplotype simulator outputs and the simulated label-run check."""

from __future__ import annotations

from pathlib import Path

import pytest
from scripts.check_label_pilot import true_events_labeled
from scripts.simulate_haplotypes import simulate


def _tsv(path: Path, rows: list[dict[str, str]]) -> None:
    columns: list[str] = []
    for row in rows:
        for key in row:
            if key not in columns:
                columns.append(key)
    lines = ["\t".join(columns)]
    lines.extend("\t".join(row.get(col, "NA") for col in columns) for row in rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


def test_simulator_writes_both_haplotypes(tmp_path: Path):
    meta = simulate(tmp_path, length=20_000, n_events=2, coverage=2, seed=1)
    hap1 = (tmp_path / "hap1.fa").read_text(encoding="utf-8")
    hap2 = (tmp_path / "hap2.fa").read_text(encoding="utf-8")
    assert hap1.startswith(">chrSim\n") and hap2.startswith(">chrSim\n")
    assert len(hap1) > 1000 and len(hap2) > 1000
    assert (tmp_path / "reads.sam").is_file()
    if meta["bam"]:
        assert Path(str(meta["bam"]) + ".bai").is_file()


def test_true_event_must_be_label_one(tmp_path: Path):
    truth = tmp_path / "truth.tsv"
    labeled = tmp_path / "labeled.tsv"
    _tsv(
        truth,
        [{"chrom": "chrSim", "start": "1000", "end": "1800", "type": "DEL", "gt": "HET", "size": "800"}],
    )
    _tsv(
        labeled,
        [
            {
                "chrom": "chrSim",
                "minbp": "1000",
                "maxbp": "1800",
                "svtype": "DEL",
                "complextype": "NA",
                "label": "1",
            }
        ],
    )
    assert true_events_labeled(labeled, truth) == 1
    _tsv(
        labeled,
        [
            {
                "chrom": "chrSim",
                "minbp": "1000",
                "maxbp": "1800",
                "svtype": "DEL",
                "complextype": "NA",
                "label": "0",
            }
        ],
    )
    with pytest.raises(SystemExit, match="not 1"):
        true_events_labeled(labeled, truth)


def test_pilot_check_reports_misses_and_requires_three(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    from scripts.check_label_pilot import assess_pilot, main

    truth = tmp_path / "truth.tsv"
    labeled = tmp_path / "labeled.tsv"
    summary = tmp_path / "dataset.summary.json"
    summary.write_text("{}\n", encoding="utf-8")
    events = [
        {"chrom": "chrSim", "start": str(10_000 * (i + 1)), "end": str(10_000 * (i + 1) + 800), "type": kind, "gt": "HOM"}
        for i, kind in enumerate(("DEL", "INV", "DUP", "inv.del"))
    ]
    _tsv(truth, events)
    calls = [
        {
            "chrom": "chrSim",
            "minbp": event["start"],
            "maxbp": event["end"],
            "svtype": event["type"] if event["type"] != "inv.del" else "DEL",
            "complextype": event["type"] if event["type"] == "inv.del" else "NA",
            "label": "1",
        }
        for event in events
        if event["type"] != "DUP"
    ]
    _tsv(labeled, calls)
    result, mismatches = assess_pilot(labeled, truth)
    assert not mismatches
    assert result.matched == 3 and result.missed == ("DUP HOM at 30000-30800",)
    assert main([str(summary), str(labeled), str(truth)]) == 0
    assert "missed DUP HOM at 30000-30800" in capsys.readouterr().err
    _tsv(labeled, [call for call in calls if call["complextype"] != "inv.del" and call["svtype"] != "INV"])
    with pytest.raises(SystemExit, match="need 3"):
        main([str(summary), str(labeled), str(truth)])


def test_missing_unreliable_bed_is_an_error(tmp_path: Path):
    from typer.testing import CliRunner

    from arcsv2.cli import app

    result = CliRunner().invoke(
        app,
        [
            "ml-label",
            "--calls",
            str(tmp_path / "calls.tab"),
            "--reference",
            str(tmp_path / "ref.fa"),
            "--hap1",
            str(tmp_path / "hap1.fa"),
            "--hap2",
            str(tmp_path / "hap2.fa"),
            "--sample",
            "sim",
            "--unreliable",
            str(tmp_path / "missing.bed"),
            "-o",
            str(tmp_path / "out.tsv"),
        ],
    )
    assert result.exit_code != 0
    assert "unreliable" in (result.output or "").lower()


def test_fetch_data_rows_skips_ml_comment_preamble(tmp_path: Path):
    from scripts.fetch_data import _rows

    samples = tmp_path / "samples.ml.tsv"
    samples.write_text(
        "# ML-0 S2 pinned samples (HPRC Release 2 x NYGC 1KGP high-coverage).\n"
        "# Extra columns: superpopulation, family_id\n"
        "sample\tpopulation\tsex\tfamily_id\thap1_url\thap2_url\n"
        "HG01884\tACB\tfemale\tBB03\thttps://example/h1.gz\thttps://example/h2.gz\n",
        encoding="utf-8",
    )
    rows = _rows(samples)
    assert len(rows) == 1
    assert rows[0]["sample"] == "HG01884"
    assert rows[0]["family_id"] == "BB03"


def test_fetch_data_refuses_placeholder_urls(tmp_path: Path):
    from scripts.fetch_data import main

    samples = tmp_path / "samples.tsv"
    samples.write_text(
        "sample\tpopulation\tsex\tcram_url\tcrai_url\thap1_url\thap2_url\t"
        "cram_checksum\tcrai_checksum\thap1_checksum\thap2_checksum\n"
        "HG002\tAJ\tmale\tREPLACE_CRAM_URL\tREPLACE_CRAI_URL\tREPLACE_HAP1_URL\t"
        "REPLACE_HAP2_URL\tmd5:a\tmd5:b\tmd5:c\tmd5:d\n",
        encoding="utf-8",
    )
    with pytest.raises(SystemExit, match="placeholder"):
        main(["--samples", str(samples), "--dest", str(tmp_path / "data")])


def test_fetch_data_download_cache_resume_and_checksum(tmp_path: Path):
    import hashlib
    import threading
    from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer

    from scripts.fetch_data import _download, fetch_samples, resolve_url

    assert resolve_url("s3://1000genomes/foo/bar").endswith("/foo/bar")
    assert resolve_url("s3://human-pangenomics/foo/bar") == (
        "https://human-pangenomics.s3.amazonaws.com/foo/bar"
    )
    payload = b"abcdefghijklmnopqrstuvwxyz" * 40
    digest = hashlib.md5(payload).hexdigest()
    root = tmp_path / "www"
    root.mkdir()
    (root / "file.bin").write_bytes(payload)
    (root / "hap.fa.gz").write_bytes(payload)

    class Handler(SimpleHTTPRequestHandler):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, directory=str(root), **kwargs)

        def do_GET(self):  # noqa: N802
            range_header = self.headers.get("Range")
            path = root / self.path.lstrip("/")
            data = path.read_bytes()
            if range_header and range_header.startswith("bytes="):
                start = int(range_header.split("=", 1)[1].split("-", 1)[0] or "0")
                body = data[start:]
                self.send_response(206)
                self.send_header("Content-Range", f"bytes {start}-{len(data) - 1}/{len(data)}")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)
                return
            self.send_response(200)
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)

        def log_message(self, format, *args):  # noqa: A003
            return

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        url = f"http://127.0.0.1:{server.server_address[1]}/file.bin"
        dest = tmp_path / "out" / "file.bin"
        first = _download(url, dest, f"md5:{digest}")
        assert first["downloaded"] == "fresh" and dest.read_bytes() == payload
        cached = _download(url, dest, f"md5:{digest}")
        assert cached["downloaded"] == "cached"
        partial = dest.with_suffix(dest.suffix + ".partial")
        dest.rename(partial)
        partial.write_bytes(payload[:100])
        resumed = _download(url, dest, f"md5:{digest}")
        assert resumed["downloaded"] == "resumed" and dest.read_bytes() == payload
        with pytest.raises(SystemExit, match="checksum"):
            _download(url, tmp_path / "out" / "bad.bin", "md5:deadbeefdeadbeefdeadbeefdeadbeef")

        hap_url = f"http://127.0.0.1:{server.server_address[1]}/hap.fa.gz"
        rows = [
            {
                "sample": "HG002",
                "hap1_url": hap_url,
                "hap2_url": hap_url,
                "hap1_checksum": f"md5:{digest}",
                "hap2_checksum": f"md5:{digest}",
            }
        ]
        records = fetch_samples(rows, tmp_path / "data", fetch_cram=False)
        assert (tmp_path / "data" / "haps" / "HG002.hap1.fa.gz").is_file()
        assert (tmp_path / "data" / "haps" / "HG002.hap2.fa.gz").is_file()
        assert any(rec.get("workflow_path", "").endswith("HG002.hap1.fa.gz") for rec in records)

        none_dest = tmp_path / "out" / "index.crai"
        waived = _download(url, none_dest, "none")
        assert waived["checksum"] == "none" and none_dest.read_bytes() == payload
    finally:
        server.shutdown()


def test_catalog_intersection_and_round_robin_selection():
    from scripts.fetch_data import (
        cram_s3_url,
        intersect_candidates,
        parse_1kgp_index,
        parse_hprc_csv,
        select_samples,
    )

    fixtures = Path(__file__).resolve().parents[0] / "fixtures" / "catalog"
    unrelated = parse_1kgp_index(
        (fixtures / "1000G_2504_high_coverage.sequence.index").read_text(encoding="utf-8"),
        related=False,
    )
    related = parse_1kgp_index(
        (fixtures / "1000G_698_related_high_coverage.sequence.index").read_text(encoding="utf-8"),
        related=True,
    )
    kg = {**unrelated, **related}
    hprc = parse_hprc_csv(
        (fixtures / "assemblies_release2_v1.0.index.csv").read_text(encoding="utf-8")
    )
    assert "ONLY1KGP" in kg and "ONLY1KGP" not in hprc
    assert "ONLYHPRC" in hprc and "ONLYHPRC" not in kg
    assert "HG00408" in kg and kg["HG00408"]["related"] == "1"
    assert "HG00438" in kg and kg["HG00438"]["related"] == "0"
    candidates = intersect_candidates(kg, hprc)
    names = {row["sample"] for row in candidates}
    assert names == {"HG00408", "HG00438", "NA12718"}
    assert "additional_698_related" in cram_s3_url("HG00408", "ERR3988762", True)
    assert "/data/ERR3231111/HG00438.final.cram" in cram_s3_url("HG00438", "ERR3231111", False)
    hap1 = hprc["HG00408"]["hap1"]["assembly"]
    assert "_pat_" in hap1 or hap1.endswith("_pat.fa.gz")
    chosen = select_samples(candidates, names=None, n_samples=2, seed=0)
    assert len(chosen) == 2
    # Round-robin across CEU then CHS (sorted pops: CEU, CHS).
    assert chosen[0]["population"] == "CEU" and chosen[1]["population"] == "CHS"
    exact = select_samples(candidates, names=["HG00408"], n_samples=None, seed=0)
    assert [row["sample"] for row in exact] == ["HG00408"]


def test_write_label_auto_yaml_and_snakemake_dry(tmp_path: Path):
    from scripts.fetch_data import write_label_auto_yaml

    yaml_path = tmp_path / "label.auto.yaml"
    write_label_auto_yaml(
        yaml_path,
        samples=["HG00408"],
        chroms=["chr20"],
        cram={"HG00408": "https://1000genomes.s3.amazonaws.com/x.final.cram"},
        sequential=True,
        header_comment="# test",
    )
    text = yaml_path.read_text(encoding="utf-8")
    assert "data_root: /data" in text
    assert "sequential: true" in text
    assert "HG00408:" in text
    assert "simple_repeats: /data/ref/simple_repeats.bed" in text
    assert "segdups: /data/ref/segdups.bed" in text
    assert "blacklist: /data/ref/blacklist.bed" in text


def test_make_ref_beds_n_runs(tmp_path: Path):
    from scripts.make_ref_beds import main, n_runs

    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\nACGTNNNNAC\n>chr2\nNNNN\n", encoding="utf-8")
    assert n_runs(fasta) == [("chr1", 4, 8), ("chr2", 0, 4)]
    centro = tmp_path / "centro.txt"
    centro.write_text("#chrom start end\nchr1 10 20 centromere\n", encoding="utf-8")
    assert main(["--fasta", str(fasta), "--centromeres", str(centro), "--out-dir", str(tmp_path)]) == 0
    gaps = (tmp_path / "gaps.bed").read_text(encoding="utf-8").strip().splitlines()
    assert "chr1\t4\t8" in gaps and "chr2\t0\t4" in gaps
    assert "chr1\t10\t20" in (tmp_path / "unreliable.bed").read_text(encoding="utf-8")


def test_fetch_filter_beds_merge_intervals():
    from scripts.fetch_filter_beds import merge_intervals

    merged = merge_intervals(
        [("chr1", 0, 10), ("chr1", 8, 20), ("chr1", 30, 40), ("chr2", 0, 5)]
    )
    assert merged == [("chr1", 0, 20), ("chr1", 30, 40), ("chr2", 0, 5)]


def test_fetch_filter_beds_dry_run(tmp_path: Path):
    from scripts.fetch_filter_beds import main

    assert main(["--dest", str(tmp_path), "--dry-run"]) == 0
    assert not (tmp_path / "simple_repeats.bed").exists()


@pytest.mark.network
@pytest.mark.skipif(__import__("os").environ.get("ARCSV2_NETWORK") != "1", reason="set ARCSV2_NETWORK=1")
def test_head_one_public_cram_and_assembly():
    from scripts.fetch_data import head_size

    cram = (
        "https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage/"
        "additional_698_related/data/ERR3988762/HG00408.final.cram"
    )
    assembly = (
        "https://human-pangenomics.s3.amazonaws.com/working/HPRC/HG00408/"
        "assemblies/release2/HG00408_pat_hprc_r2_v1.0.1.fa.gz"
    )
    assert head_size(cram) > 1_000_000_000
    assert head_size(assembly) > 100_000_000
