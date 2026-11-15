"""Download pinned reference and sample files for labeling.

Manual mode (existing samples.tsv):

    python scripts/fetch_data.py --samples workflows/config/samples.tsv --dest /data

Auto mode (catalog intersection → config → download):

    python scripts/fetch_data.py --auto --n-samples 6 --chroms chr20,chr21,chr22
    python scripts/fetch_data.py --auto --samples HG00408 --chroms chr20 --dry-run

Writes ``<dest>/manifest.lock.json``. Resumes partial downloads with HTTP Range.
``s3://`` URLs rewrite to the public HTTPS endpoints. Haplotypes are also linked
under ``<dest>/haps/`` for ``label.smk`` when ``data_root`` is ``/data``.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import sys
import time
import urllib.error
import urllib.request
from collections import defaultdict
from pathlib import Path

REQUIRED = (
    "sample",
    "cram_url",
    "crai_url",
    "hap1_url",
    "hap2_url",
    "cram_checksum",
    "crai_checksum",
    "hap1_checksum",
    "hap2_checksum",
)

_CHECKSUM_ALIASES = {
    "cram_checksum": ("cram_sha256",),
    "crai_checksum": ("crai_sha256",),
    "hap1_checksum": ("hap1_sha256",),
    "hap2_checksum": ("hap2_sha256",),
}

S3_HTTPS = {
    "1000genomes": "https://1000genomes.s3.amazonaws.com",
    "human-pangenomics": "https://human-pangenomics.s3.amazonaws.com",
    "giab": "https://giab.s3.amazonaws.com",
}

# Pinned catalogs (round 13). HG002 is not in 1KGP; keep it out of --auto.
HPRC_COMMIT = "2a7e3a5b1ca10454762699176b1a0a1027ccbf21"
HPRC_CSV = (
    f"https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/"
    f"{HPRC_COMMIT}/data_tables/assemblies_release2_v1.0.index.csv"
)
INDEX_UNRELATED = (
    "s3://1000genomes/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index"
)
INDEX_RELATED = (
    "s3://1000genomes/1000G_2504_high_coverage/additional_698_related/"
    "1000G_698_related_high_coverage.sequence.index"
)
REF_FA = (
    "s3://1000genomes/technical/reference/GRCh38_reference_genome/"
    "GRCh38_full_analysis_set_plus_decoy_hla.fa"
)
REF_FAI = REF_FA + ".fai"
CENTROMERES = (
    "s3://1000genomes/technical/reference/GRCh38_reference_genome/"
    "20150713_location_of_centromeres_and_other_regions.txt"
)
SEGDUPS_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz"

# Paper merge-filter BEDs (Methods p. 27). Built by scripts/fetch_filter_beds.py.
FILTER_BED_URLS = {
    "simple_repeats": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz",
    "rmsk": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz",
    "microsat": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/microsat.txt.gz",
    "segdups": SEGDUPS_URL,
    "blacklist": "https://cf.10xgenomics.com/supp/genome/GRCh38/svblacklist.bed",
}

USER_AGENT = "arcsv2-fetch/1"


def _rows(path: Path) -> list[dict[str, str]]:
    """Parse a samples TSV, skipping `#` comment lines before/after the header.

    The header may be plain or `#`-prefixed. Extra columns (e.g. ``family_id``,
    ``superpopulation``) are kept; only REQUIRED keys are validated later.
    """
    lines = path.read_text(encoding="utf-8").splitlines()
    header: list[str] | None = None
    rows: list[dict[str, str]] = []
    for raw in lines:
        line = raw.lstrip("\ufeff")
        if not line.strip():
            continue
        if line.startswith("#"):
            # Allow a commented header row; ignore other comment lines (e.g. "pinned samples").
            maybe = [c.strip() for c in line.lstrip("#").strip().split("\t")]
            if header is None and maybe and maybe[0] == "sample":
                header = maybe
            continue
        cells = [c.strip() for c in line.split("\t")]
        if header is None:
            header = cells
            continue
        rows.append(dict(zip(header, cells)))
    return rows


def _checksum_value(row: dict[str, str], key: str) -> str:
    if row.get(key, "").strip():
        return row[key].strip()
    for alias in _CHECKSUM_ALIASES.get(key, ()):
        if row.get(alias, "").strip():
            return row[alias].strip()
    return ""


def _parse_checksum(expected: str) -> tuple[str, str] | None:
    """Return ``(algo, hex)`` or ``None`` when checksums are waived (``none``)."""
    text = (expected or "").strip()
    if not text or text.lower() == "none":
        return None
    if text.startswith("md5:"):
        return "md5", text[4:].lower()
    if text.startswith("sha256:"):
        return "sha256", text[7:].lower()
    if len(text) == 32 and all(c in "0123456789abcdefABCDEF" for c in text):
        return "md5", text.lower()
    return "sha256", text.lower()


def _digest_file(path: Path, algo: str) -> str:
    digest = hashlib.md5() if algo == "md5" else hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def resolve_url(url: str) -> str:
    """Rewrite ``s3://`` to the public HTTPS endpoint when known."""
    if not url.startswith("s3://"):
        return url
    rest = url[5:]
    bucket, _, key = rest.partition("/")
    base = S3_HTTPS.get(bucket)
    if base is None:
        raise SystemExit(
            f"unknown s3 bucket {bucket!r}. Known: {', '.join(sorted(S3_HTTPS))} "
            "or pass an https:// URL."
        )
    return f"{base}/{key}"


def _is_placeholder(url: str) -> bool:
    return (not url) or url.startswith("REPLACE_")


def _link_or_copy(src: Path, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() or dest.is_symlink():
        dest.unlink()
    try:
        os.symlink(src.resolve(), dest)
    except OSError:
        shutil.copyfile(src, dest)


def _download(url: str, dest: Path, expected: str) -> dict:
    dest.parent.mkdir(parents=True, exist_ok=True)
    started = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    parsed = _parse_checksum(expected)
    if parsed is None:
        algo, digest_expected = "", ""
    else:
        algo, digest_expected = parsed
    if dest.is_file():
        if digest_expected and _digest_file(dest, algo) == digest_expected:
            return {
                "url": url,
                "path": str(dest),
                "checksum": f"{algo}:{digest_expected}",
                "size": dest.stat().st_size,
                "downloaded": "cached",
                "time": started,
            }
        if not digest_expected:
            return {
                "url": url,
                "path": str(dest),
                "checksum": "none",
                "size": dest.stat().st_size,
                "downloaded": "cached",
                "time": started,
            }
    if _is_placeholder(url):
        raise SystemExit(
            f"refusing to download placeholder URL for {dest.name}. "
            "Fill workflows/config/samples.tsv with pinned URLs and checksums "
            "(see docs/PRODUCTION_PLAN.md §3.1)."
        )
    resolved = resolve_url(url)
    tmp = dest.with_suffix(dest.suffix + ".partial")
    existing = tmp.stat().st_size if tmp.is_file() else 0
    headers = {"User-Agent": USER_AGENT}
    if existing:
        headers["Range"] = f"bytes={existing}-"
    request = urllib.request.Request(resolved, headers=headers)
    try:
        response = urllib.request.urlopen(request, timeout=120)
    except Exception:
        if existing:
            tmp.unlink(missing_ok=True)
            existing = 0
            request = urllib.request.Request(resolved, headers={"User-Agent": USER_AGENT})
            response = urllib.request.urlopen(request, timeout=120)
        else:
            raise
    status = getattr(response, "status", 200)
    if existing and status == 200:
        tmp.unlink(missing_ok=True)
        existing = 0
    mode = "ab" if existing else "wb"
    with response, tmp.open(mode) as handle:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            handle.write(chunk)
    if digest_expected:
        digest = _digest_file(tmp, algo)
        if digest != digest_expected:
            tmp.unlink(missing_ok=True)
            raise SystemExit(
                f"checksum mismatch for {url}: got {algo}:{digest}, expected {algo}:{digest_expected}"
            )
        checksum = f"{algo}:{digest}"
    else:
        checksum = "none"
    tmp.replace(dest)
    return {
        "url": url,
        "resolved": resolved,
        "path": str(dest),
        "checksum": checksum,
        "size": dest.stat().st_size,
        "downloaded": "resumed" if existing else "fresh",
        "time": started,
    }


def head_size(url: str) -> int:
    """HTTPS HEAD; return Content-Length or 0."""
    resolved = resolve_url(url)
    request = urllib.request.Request(
        resolved, method="HEAD", headers={"User-Agent": USER_AGENT}
    )
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            return int(response.headers.get("Content-Length") or 0)
    except urllib.error.HTTPError as exc:
        raise SystemExit(f"HEAD failed for {resolved}: HTTP {exc.code}") from exc


def fetch_reference(dest: Path, ref_url: str, ref_sum: str, gaps_url: str, gaps_sum: str) -> list[dict]:
    records = []
    if ref_url:
        records.append(_download(ref_url, dest / "ref" / "GRCh38.fa", ref_sum))
    if gaps_url:
        records.append(_download(gaps_url, dest / "ref" / "gaps.bed", gaps_sum))
    return records


def _hap_name(sample: str, which: str, url: str) -> str:
    return f"{sample}.{which}.fa.gz" if url.rstrip("/").endswith(".gz") else f"{sample}.{which}.fa"


def fetch_samples(rows: list[dict[str, str]], dest: Path, fetch_cram: bool) -> list[dict]:
    records: list[dict] = []
    haps_root = dest / "haps"
    for row in rows:
        sample = row["sample"]
        base = dest / "samples" / sample
        hap1_name = _hap_name(sample, "hap1", row.get("hap1_url", ""))
        hap2_name = _hap_name(sample, "hap2", row.get("hap2_url", ""))
        pairs = [
            ("hap1_url", "hap1_checksum", base / "haps" / hap1_name, haps_root / hap1_name),
            ("hap2_url", "hap2_checksum", base / "haps" / hap2_name, haps_root / hap2_name),
        ]
        for url_key, sum_key, path, link in pairs:
            record = _download(row.get(url_key, ""), path, _checksum_value(row, sum_key))
            _link_or_copy(path, link)
            record["workflow_path"] = str(link)
            records.append(record)
        if fetch_cram:
            for url_key, sum_key, path in (
                ("cram_url", "cram_checksum", base / "reads" / f"{sample}.cram"),
                ("crai_url", "crai_checksum", base / "reads" / f"{sample}.cram.crai"),
            ):
                records.append(_download(row.get(url_key, ""), path, _checksum_value(row, sum_key)))
    return records


def parse_1kgp_index(text: str, *, related: bool) -> dict[str, dict[str, str]]:
    """SAMPLE_NAME → {run_id, population, md5, related}."""
    header: list[str] | None = None
    out: dict[str, dict[str, str]] = {}
    for line in text.splitlines():
        if not line or line.startswith("##"):
            continue
        if line.startswith("#ENA_FILE_PATH") or line.startswith("ENA_FILE_PATH"):
            header = line.lstrip("#").split("\t")
            continue
        parts = line.split("\t")
        if header is None:
            # Related index may omit the header; use the unrelated column order.
            header = [
                "ENA_FILE_PATH",
                "MD5SUM",
                "RUN_ID",
                "STUDY_ID",
                "STUDY_NAME",
                "CENTER_NAME",
                "SUBMISSION_ID",
                "SUBMISSION_DATE",
                "SAMPLE_ID",
                "SAMPLE_NAME",
                "POPULATION",
            ]
        row = {header[i]: parts[i] if i < len(parts) else "" for i in range(len(header))}
        sample = row.get("SAMPLE_NAME", "").strip()
        run_id = row.get("RUN_ID", "").strip()
        if not sample or not run_id:
            continue
        out[sample] = {
            "run_id": run_id,
            "population": row.get("POPULATION", "").strip() or "UNK",
            "md5": row.get("MD5SUM", "").strip().lower(),
            "related": "1" if related else "0",
        }
    return out


def cram_s3_url(sample: str, run_id: str, related: bool) -> str:
    if related:
        return (
            f"s3://1000genomes/1000G_2504_high_coverage/additional_698_related/data/"
            f"{run_id}/{sample}.final.cram"
        )
    return f"s3://1000genomes/1000G_2504_high_coverage/data/{run_id}/{sample}.final.cram"


def parse_hprc_csv(text: str) -> dict[str, dict[str, dict[str, str]]]:
    """sample_id → {hap1|hap2 → {assembly, assembly_md5, fai, gzi}}."""
    reader = csv.DictReader(text.splitlines())
    by_sample: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    for row in reader:
        sample = (row.get("sample_id") or "").strip()
        hap = (row.get("haplotype") or "").strip()
        if not sample or hap not in ("1", "2"):
            continue
        key = "hap1" if hap == "1" else "hap2"
        by_sample[sample][key] = {
            "assembly": (row.get("assembly") or "").strip(),
            "assembly_md5": (row.get("assembly_md5") or "").strip(),
            "fai": (row.get("assembly_fai") or "").strip(),
            "gzi": (row.get("assembly_gzi") or "").strip(),
        }
    return {sample: haps for sample, haps in by_sample.items() if "hap1" in haps and "hap2" in haps}


def intersect_candidates(
    kg: dict[str, dict[str, str]], hprc: dict[str, dict[str, dict[str, str]]]
) -> list[dict[str, str]]:
    """Samples in both catalogs with both haplotypes."""
    rows: list[dict[str, str]] = []
    for sample in sorted(set(kg) & set(hprc)):
        info = kg[sample]
        haps = hprc[sample]
        related = info["related"] == "1"
        cram = cram_s3_url(sample, info["run_id"], related)
        rows.append(
            {
                "sample": sample,
                "population": info["population"],
                "sex": "",
                "cram_url": resolve_url(cram),
                "crai_url": resolve_url(cram + ".crai"),
                "hap1_url": resolve_url(haps["hap1"]["assembly"]),
                "hap2_url": resolve_url(haps["hap2"]["assembly"]),
                "hap1_md5_url": resolve_url(haps["hap1"]["assembly_md5"]),
                "hap2_md5_url": resolve_url(haps["hap2"]["assembly_md5"]),
                "hap1_fai_url": resolve_url(haps["hap1"]["fai"]),
                "hap2_fai_url": resolve_url(haps["hap2"]["fai"]),
                "hap1_gzi_url": resolve_url(haps["hap1"]["gzi"]),
                "hap2_gzi_url": resolve_url(haps["hap2"]["gzi"]),
                "cram_checksum": f"md5:{info['md5']}" if info["md5"] else "none",
                "crai_checksum": "none",
                "hap1_checksum": "",
                "hap2_checksum": "",
                "related": info["related"],
                "run_id": info["run_id"],
            }
        )
    rows.sort(key=lambda r: (r["population"], r["sample"]))
    return rows


def select_samples(
    candidates: list[dict[str, str]],
    *,
    names: list[str] | None,
    n_samples: int | None,
    seed: int,
) -> list[dict[str, str]]:
    if names:
        by_name = {row["sample"]: row for row in candidates}
        missing = [name for name in names if name not in by_name]
        if missing:
            raise SystemExit(
                f"sample(s) not in 1KGP∩HPRC intersection: {', '.join(missing)}"
            )
        return [by_name[name] for name in names]
    if not n_samples or n_samples < 1:
        raise SystemExit("pass --n-samples N or --samples A,B,...")
    # Round-robin across populations for diversity; seed is a tie-breaker only.
    by_pop: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in candidates:
        by_pop[row["population"]].append(row)
    for pop in by_pop:
        by_pop[pop].sort(key=lambda r: (r["sample"], seed))
    pops = sorted(by_pop)
    chosen: list[dict[str, str]] = []
    indexes = {pop: 0 for pop in pops}
    while len(chosen) < n_samples and pops:
        progressed = False
        for pop in list(pops):
            idx = indexes[pop]
            if idx >= len(by_pop[pop]):
                pops.remove(pop)
                continue
            chosen.append(by_pop[pop][idx])
            indexes[pop] = idx + 1
            progressed = True
            if len(chosen) >= n_samples:
                break
        if not progressed:
            break
    if len(chosen) < n_samples:
        raise SystemExit(f"only {len(chosen)} candidates available (wanted {n_samples})")
    return chosen


def read_md5_sidecar(url: str) -> str:
    """Download a tiny ``.md5`` file and return ``md5:<hex>``."""
    resolved = resolve_url(url)
    request = urllib.request.Request(resolved, headers={"User-Agent": USER_AGENT})
    text = urllib.request.urlopen(request, timeout=60).read().decode("utf-8", "replace")
    for line in text.splitlines():
        parts = line.split()
        if parts and all(c in "0123456789abcdefABCDEF" for c in parts[0]) and len(parts[0]) == 32:
            return f"md5:{parts[0].lower()}"
    raise SystemExit(f"no MD5 found in {resolved}")


def write_samples_tsv(path: Path, rows: list[dict[str, str]], header_comment: str) -> None:
    cols = [
        "sample",
        "population",
        "sex",
        "cram_url",
        "crai_url",
        "hap1_url",
        "hap2_url",
        "cram_checksum",
        "crai_checksum",
        "hap1_checksum",
        "hap2_checksum",
    ]
    lines = [header_comment.rstrip(), "\t".join(cols)]
    for row in rows:
        lines.append("\t".join(row.get(col, "") for col in cols))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_label_auto_yaml(
    path: Path,
    *,
    samples: list[str],
    chroms: list[str],
    cram: dict[str, str],
    sequential: bool,
    header_comment: str,
) -> None:
    lines = [header_comment.rstrip(), "data_root: /data"]
    lines.append("reference: /data/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa")
    lines.append("gaps: /data/ref/gaps.bed")
    lines.append("unreliable: /data/ref/unreliable.bed")
    # Paper merge-filter BEDs (fetch with: python scripts/fetch_filter_beds.py --dest /data/ref).
    lines.append("simple_repeats: /data/ref/simple_repeats.bed")
    lines.append("segdups: /data/ref/segdups.bed")
    lines.append("blacklist: /data/ref/blacklist.bed")
    lines.append("chroms:")
    for chrom in chroms:
        lines.append(f"  - {chrom}")
    lines.append("samples:")
    for sample in samples:
        lines.append(f"  - {sample}")
    lines.append("cram:")
    for sample in samples:
        lines.append(f"  {sample}: {cram[sample]}")
    lines.append(f"sequential: {'true' if sequential else 'false'}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _cache_catalog(url: str, dest: Path, refresh: bool) -> tuple[Path, str]:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.is_file() and not refresh:
        digest = _digest_file(dest, "md5")
        return dest, digest
    record = _download(url, dest, "none")
    digest = _digest_file(dest, "md5")
    record["checksum"] = f"md5:{digest}"
    return dest, digest


def run_auto(args: argparse.Namespace) -> int:
    dest: Path = args.dest
    catalog_dir = dest / "catalog"
    kg_u_path, kg_u_md5 = _cache_catalog(
        INDEX_UNRELATED, catalog_dir / "1000G_2504_high_coverage.sequence.index", args.refresh_catalog
    )
    kg_r_path, kg_r_md5 = _cache_catalog(
        INDEX_RELATED,
        catalog_dir / "1000G_698_related_high_coverage.sequence.index",
        args.refresh_catalog,
    )
    hprc_path, hprc_md5 = _cache_catalog(
        HPRC_CSV, catalog_dir / f"assemblies_release2_v1.0.{HPRC_COMMIT[:12]}.csv", args.refresh_catalog
    )
    kg = parse_1kgp_index(kg_u_path.read_text(encoding="utf-8"), related=False)
    kg.update(parse_1kgp_index(kg_r_path.read_text(encoding="utf-8"), related=True))
    hprc = parse_hprc_csv(hprc_path.read_text(encoding="utf-8"))
    candidates = intersect_candidates(kg, hprc)
    names = [s.strip() for s in (args.sample_list or "").split(",") if s.strip()] or None
    chosen = select_samples(
        candidates, names=names, n_samples=args.n_samples, seed=args.seed
    )
    print(f"selected {len(chosen)} sample(s): " + ", ".join(r["sample"] for r in chosen))

    total_bytes = 0
    for row in chosen:
        for key in ("cram_url", "crai_url", "hap1_url", "hap2_url"):
            size = head_size(row[key])
            row[f"{key}_size"] = str(size)
            # Assemblies always download; CRAMs only with --fetch-cram.
            if key.startswith("hap") or args.fetch_cram:
                total_bytes += size
        if not args.dry_run:
            row["hap1_checksum"] = read_md5_sidecar(row["hap1_md5_url"])
            row["hap2_checksum"] = read_md5_sidecar(row["hap2_md5_url"])
        else:
            row["hap1_checksum"] = "md5:dryrun"
            row["hap2_checksum"] = "md5:dryrun"

    # Reference once.
    for url in (REF_FA, REF_FAI, CENTROMERES):
        total_bytes += head_size(url)
    if args.segdups or args.filter_beds:
        total_bytes += head_size(SEGDUPS_URL)
    if args.filter_beds:
        for key, url in FILTER_BED_URLS.items():
            if key == "segdups":
                continue  # already counted
            total_bytes += head_size(url)

    gb = total_bytes / (1024**3)
    print(f"plan: ~{gb:.1f} GB to download under {dest} (CRAMs streamed unless --fetch-cram)")

    chroms = [c.strip() for c in args.chroms.split(",") if c.strip()]
    if not chroms:
        raise SystemExit("--chroms must list at least one chromosome")
    comment = (
        f"# generated {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} by fetch_data.py --auto\n"
        f"# catalogs md5: unrelated={kg_u_md5} related={kg_r_md5} hprc={hprc_md5} "
        f"commit={HPRC_COMMIT}"
    )
    samples_tsv = Path(args.samples_out)
    label_yaml = Path(args.label_yaml)
    write_samples_tsv(samples_tsv, chosen, comment)
    write_label_auto_yaml(
        label_yaml,
        samples=[r["sample"] for r in chosen],
        chroms=chroms,
        cram={r["sample"]: r["cram_url"] for r in chosen},
        sequential=not args.no_sequential,
        header_comment=comment,
    )
    print(f"wrote {samples_tsv}")
    print(f"wrote {label_yaml}")

    if args.dry_run:
        print("dry-run: no downloads")
        return 0

    free = shutil.disk_usage(dest if dest.exists() else dest.parent).free
    if free < int(total_bytes * 1.2):
        raise SystemExit(
            f"free space under {dest} is {free / 1024**3:.1f} GB; "
            f"need ~{gb * 1.2:.1f} GB (1.2× plan)"
        )
    if not args.yes:
        answer = input(f"Download ~{gb:.1f} GB to {dest}? [y/N] ").strip().lower()
        if answer not in ("y", "yes"):
            print("aborted")
            return 1

    records: list[dict] = []
    ref_name = "GRCh38_full_analysis_set_plus_decoy_hla.fa"
    records.append(_download(REF_FA, dest / "ref" / ref_name, "none"))
    records.append(_download(REF_FAI, dest / "ref" / f"{ref_name}.fai", "none"))
    records.append(_download(CENTROMERES, dest / "ref" / "centromeres.txt", "none"))
    segdups_path = None
    if args.segdups or args.filter_beds:
        segdups_path = dest / "ref" / "genomicSuperDups.txt.gz"
        records.append(_download(SEGDUPS_URL, segdups_path, "none"))

    # Assembly indexes (fai/gzi) with waived checksums.
    for row in chosen:
        sample = row["sample"]
        base = dest / "samples" / sample
        for which, url_key, sum_key in (
            ("hap1", "hap1_url", "hap1_checksum"),
            ("hap2", "hap2_url", "hap2_checksum"),
        ):
            url = row[url_key]
            name = _hap_name(sample, which, url)
            path = base / "haps" / name
            record = _download(url, path, row[sum_key])
            _link_or_copy(path, dest / "haps" / name)
            record["workflow_path"] = str(dest / "haps" / name)
            records.append(record)
            records.append(_download(row[f"{which}_fai_url"], Path(str(path) + ".fai"), "none"))
            records.append(_download(row[f"{which}_gzi_url"], Path(str(path) + ".gzi"), "none"))
        if args.fetch_cram:
            records.append(
                _download(row["cram_url"], base / "reads" / f"{sample}.cram", row["cram_checksum"])
            )
            records.append(
                _download(row["crai_url"], base / "reads" / f"{sample}.cram.crai", "none")
            )

    import subprocess

    bed_cmd = [
        sys.executable,
        str(Path(__file__).resolve().parent / "make_ref_beds.py"),
        "--fasta",
        str(dest / "ref" / ref_name),
        "--centromeres",
        str(dest / "ref" / "centromeres.txt"),
        "--out-dir",
        str(dest / "ref"),
    ]
    if segdups_path is not None:
        bed_cmd.extend(["--segdups", str(segdups_path)])
    subprocess.check_call(bed_cmd)

    if args.filter_beds:
        filter_cmd = [
            sys.executable,
            str(Path(__file__).resolve().parent / "fetch_filter_beds.py"),
            "--dest",
            str(dest / "ref"),
        ]
        subprocess.check_call(filter_cmd)

    lock = {
        "created": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "mode": "auto",
        "hprc_commit": HPRC_COMMIT,
        "catalog_md5": {"unrelated": kg_u_md5, "related": kg_r_md5, "hprc": hprc_md5},
        "samples": samples_tsv.as_posix(),
        "label_yaml": label_yaml.as_posix(),
        "files": records,
    }
    dest.mkdir(parents=True, exist_ok=True)
    lock_path = dest / "manifest.lock.json"
    lock_path.write_text(json.dumps(lock, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {lock_path} ({len(records)} files)")
    print(
        "next: snakemake -s workflows/label.smk "
        f"--configfile {label_yaml.as_posix()} -j 2"
    )
    return 0


def run_manual(args: argparse.Namespace) -> int:
    if not args.samples.is_file():
        raise SystemExit(
            f"missing {args.samples}. Copy workflows/config/samples.tsv.example "
            "and fill pinned URLs and checksums, or use --auto."
        )
    rows = _rows(args.samples)
    if not rows:
        raise SystemExit(f"{args.samples} has no sample rows")
    for row in rows:
        missing = []
        for key in REQUIRED:
            if key.endswith("_checksum"):
                if not _checksum_value(row, key):
                    missing.append(key)
            elif key not in row or not str(row[key]).strip():
                missing.append(key)
        if missing:
            raise SystemExit(f"sample {row.get('sample', '?')} missing columns: {', '.join(missing)}")
    records = fetch_reference(
        args.dest, args.ref_url, args.ref_checksum, args.gaps_url, args.gaps_checksum
    )
    records.extend(fetch_samples(rows, args.dest, args.fetch_cram))
    lock = {
        "created": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "samples": args.samples.as_posix(),
        "files": records,
    }
    args.dest.mkdir(parents=True, exist_ok=True)
    lock_path = args.dest / "manifest.lock.json"
    lock_path.write_text(json.dumps(lock, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {lock_path} ({len(records)} files)")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--auto", action="store_true", help="Build config from public catalogs")
    parser.add_argument("--samples", type=Path, default=Path("workflows/config/samples.tsv"))
    parser.add_argument("--dest", type=Path, default=Path("/data"))
    parser.add_argument("--ref-url", default="", help="Pinned GRCh38 analysis-set FASTA URL")
    parser.add_argument("--ref-checksum", default="", help="md5:... or sha256:... or none")
    parser.add_argument("--gaps-url", default="", help="Pinned gaps BED URL")
    parser.add_argument("--gaps-checksum", default="")
    parser.add_argument(
        "--fetch-cram",
        action="store_true",
        help="Also download CRAMs. Pilots normally stream by region instead.",
    )
    parser.add_argument("--n-samples", type=int, default=None)
    parser.add_argument(
        "--sample-list",
        default="",
        dest="sample_list",
        help="Comma-separated sample IDs (alias of --samples for --auto)",
    )
    parser.add_argument("--chroms", default="chr20,chr21,chr22")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--yes", action="store_true", help="Skip the download confirmation")
    parser.add_argument("--refresh-catalog", action="store_true")
    parser.add_argument("--segdups", action="store_true", help="Also fetch UCSC segdups")
    parser.add_argument(
        "--filter-beds",
        action="store_true",
        help=(
            "Fetch paper merge-filter annotation BEDs into <dest>/ref "
            "(simple_repeats, segdups, blacklist; see scripts/fetch_filter_beds.py)"
        ),
    )
    parser.add_argument("--no-sequential", action="store_true")
    parser.add_argument(
        "--samples-out",
        type=Path,
        default=Path("workflows/config/samples.tsv"),
    )
    parser.add_argument(
        "--label-yaml",
        type=Path,
        default=Path("workflows/config/label.auto.yaml"),
    )
    args = parser.parse_args(argv)

    # Allow ``--samples HG00408,HG00438`` as a convenience alias of --sample-list.
    if args.auto and args.samples != Path("workflows/config/samples.tsv"):
        text = str(args.samples)
        if "," in text or (not Path(text).suffix and not Path(text).is_file()):
            args.sample_list = text
            args.samples = Path("workflows/config/samples.tsv")

    if args.auto:
        return run_auto(args)
    return run_manual(args)


if __name__ == "__main__":
    raise SystemExit(main())
