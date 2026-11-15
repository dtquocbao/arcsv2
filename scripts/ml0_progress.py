#!/usr/bin/env python3
"""ML-0 / label.smk progress snapshot from the data volume.

Reads sample/chrom lists from a label config YAML and reports filesystem
progress (BAMs, calls, merges, labels), call-table SV / genotype tallies,
optional live genotype lines from the host log, and an approximate DAG
step fraction.

Examples (inside ``arcsv2:ml0``, ``/data`` mounted)::

    python scripts/ml0_progress.py
    python scripts/ml0_progress.py --matrix
    python scripts/ml0_progress.py --by-sample HG01884

From Windows (repo path without apostrophe)::

    docker run --rm --entrypoint bash `
      -v arcsv2-data:/data -v ${PWD}:/work -w /work -e PYTHONPATH=/work/python `
      arcsv2:ml0 -lc "python scripts/ml0_progress.py --matrix"
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import sys
from collections import Counter, defaultdict
from pathlib import Path

CHROM_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX"]

LIVE_RE = re.compile(
    r"\[call:(?P<chrom>chr[\w]+)\s+(?P<sec>[\d.]+)s\]\s+"
    r"(?P<stage>load_bam|fit_libs|discover|blocks|build_graph|insertions|genotype)"
    r"(?::\s*(?P<detail>.*))?"
)
LOCUS_RE = re.compile(r"locus\s+(\d+)/(\d+)")
SITES_RE = re.compile(r"(\d+)/(\d+)\s+sites")


def _load_yaml(path: Path) -> dict:
    try:
        import yaml  # type: ignore
    except ImportError as exc:  # pragma: no cover
        raise SystemExit("PyYAML required (available in arcsv2:ml0)") from exc
    with path.open(encoding="utf-8") as fh:
        data = yaml.safe_load(fh) or {}
    if not isinstance(data, dict):
        raise SystemExit(f"expected mapping in {path}")
    return data


def _chrom_sort_key(chrom: str) -> tuple[int, str]:
    if chrom in CHROM_ORDER:
        return (CHROM_ORDER.index(chrom), chrom)
    return (len(CHROM_ORDER), chrom)


def _nonempty(path: Path) -> bool:
    try:
        return path.is_file() and path.stat().st_size > 0
    except OSError:
        return False


def _free_gb(root: Path) -> str:
    try:
        usage = shutil.disk_usage(root)
        return f"{usage.free / (1024**3):.0f}G free / {usage.total / (1024**3):.0f}G"
    except OSError:
        return "n/a"


def _primary_svtype(svtype: str) -> str:
    token = (svtype or "").split(",", 1)[0].strip().upper()
    if not token or token == "NA":
        return "OTHER"
    if token.startswith("DEL"):
        return "DEL"
    if token.startswith("DUP"):
        return "DUP"
    if token.startswith("INS"):
        return "INS"
    if token.startswith("INV"):
        return "INV"
    if token.startswith("BND"):
        return "BND"
    return token


def _tally_call_tab(path: Path) -> tuple[int, Counter[str], Counter[str], Counter[str]]:
    """Return (n_loci, svtype, gt, complextype) counters."""
    sv: Counter[str] = Counter()
    gt: Counter[str] = Counter()
    cx: Counter[str] = Counter()
    n = 0
    try:
        with path.open(encoding="utf-8", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                n += 1
                sv[_primary_svtype(row.get("svtype", ""))] += 1
                g = (row.get("gt") or "NA").strip() or "NA"
                gt[g] += 1
                c = (row.get("complextype") or "NA").strip() or "NA"
                if c and c != "NA":
                    cx[c] += 1
    except OSError:
        return 0, sv, gt, cx
    return n, sv, gt, cx


def _tally_labeled(path: Path) -> tuple[int, Counter[str], Counter[str]]:
    """Return (n_rows, label, sv_class-or-svtype)."""
    labels: Counter[str] = Counter()
    classes: Counter[str] = Counter()
    n = 0
    try:
        with path.open(encoding="utf-8", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                n += 1
                lab = (row.get("label") or "NA").strip() or "NA"
                labels[lab] += 1
                cls = (row.get("sv_class") or row.get("complextype") or row.get("svtype") or "NA")
                classes[_primary_svtype(str(cls).split(",", 1)[0])] += 1
    except OSError:
        return 0, labels, classes
    return n, labels, classes


def _live_jobs(host_log: Path, limit: int = 12) -> list[str]:
    if not host_log.is_file():
        return []
    # Read last ~256 KiB only.
    try:
        size = host_log.stat().st_size
        with host_log.open("rb") as fh:
            fh.seek(max(0, size - 262_144))
            text = fh.read().decode("utf-8", errors="replace")
    except OSError:
        return []
    latest: dict[str, tuple[float, str]] = {}
    for line in text.splitlines():
        m = LIVE_RE.search(line)
        if not m:
            continue
        chrom = m.group("chrom")
        sec = float(m.group("sec"))
        stage = m.group("stage")
        detail = (m.group("detail") or "").strip()
        extra = ""
        if stage == "genotype":
            lm = LOCUS_RE.search(detail)
            if lm:
                extra = f" locus {lm.group(1)}/{lm.group(2)}"
        elif stage == "insertions":
            sm = SITES_RE.search(detail)
            if sm:
                extra = f" sites {sm.group(1)}/{sm.group(2)}"
        elif stage == "load_bam" and "done" in detail:
            extra = " done"
        key = chrom  # last stage wins per chrom in the window
        latest[key] = (sec, f"{chrom}  {stage}{extra}  ({sec:.0f}s)")
    rows = sorted(latest.values(), key=lambda t: -t[0])[:limit]
    return [r[1] for r in rows]


def collect(cfg: dict, host_log: Path | None) -> dict:
    data_root = Path(str(cfg.get("data_root") or "/data"))
    results_dir = data_root / str(cfg.get("results_dir") or "results")
    labels_dir = data_root / str(cfg.get("labels_dir") or "labels")
    samples = list(cfg.get("samples") or [])
    chroms = list(cfg.get("chroms") or CHROM_ORDER)
    if not samples:
        raise SystemExit("config has no samples")

    n_s, n_c = len(samples), len(chroms)
    expected = {
        "fetch_region": n_s * n_c,
        "call": n_s * n_c,
        "index_haps": n_s * 2,
        "merge": n_s,
        "label": n_s,
        "dataset": 1,
    }
    expected["total"] = sum(expected.values())

    bam_ok: set[tuple[str, str]] = set()
    call_ok: set[tuple[str, str]] = set()
    merge_ok: set[str] = set()
    label_ok: set[str] = set()
    hap_ok = 0
    crams = 0

    sv_all: Counter[str] = Counter()
    gt_all: Counter[str] = Counter()
    cx_all: Counter[str] = Counter()
    loci_all = 0
    label_counts: Counter[str] = Counter()
    label_class: Counter[str] = Counter()
    label_rows = 0

    per_sample: dict[str, dict] = {}
    chrom_calls: Counter[str] = Counter()
    chrom_bams: Counter[str] = Counter()

    for sample in samples:
        s_bam = 0
        s_call = 0
        s_loci = 0
        s_sv: Counter[str] = Counter()
        for chrom in chroms:
            bam = data_root / "work" / sample / f"{chrom}.bam"
            bai = data_root / "work" / sample / f"{chrom}.bam.bai"
            if _nonempty(bam) and _nonempty(bai):
                bam_ok.add((sample, chrom))
                s_bam += 1
                chrom_bams[chrom] += 1
            tab = results_dir / sample / chrom / "arcsv_out.tab"
            if _nonempty(tab):
                call_ok.add((sample, chrom))
                s_call += 1
                chrom_calls[chrom] += 1
                n, sv, gt, cx = _tally_call_tab(tab)
                s_loci += n
                loci_all += n
                sv_all.update(sv)
                gt_all.update(gt)
                cx_all.update(cx)
                s_sv.update(sv)
        merged = results_dir / sample / "merged.tab"
        if _nonempty(merged):
            merge_ok.add(sample)
        lab = labels_dir / f"{sample}.labeled.tsv"
        if _nonempty(lab):
            label_ok.add(sample)
            nr, lc, lcl = _tally_labeled(lab)
            label_rows += nr
            label_counts.update(lc)
            label_class.update(lcl)
        for which in ("hap1", "hap2"):
            if _nonempty(data_root / "haps" / f"{sample}.{which}.mmi"):
                hap_ok += 1
        cram = data_root / "samples" / sample / "reads" / f"{sample}.cram"
        if _nonempty(cram):
            crams += 1
        per_sample[sample] = {
            "bam": s_bam,
            "call": s_call,
            "loci": s_loci,
            "sv": s_sv,
            "merged": sample in merge_ok,
            "labeled": sample in label_ok,
            "cram": _nonempty(cram),
        }

    dataset = _nonempty(labels_dir / "dataset.tsv")
    done = {
        "fetch_region": len(bam_ok),
        "call": len(call_ok),
        "index_haps": hap_ok,
        "merge": len(merge_ok),
        "label": len(label_ok),
        "dataset": 1 if dataset else 0,
    }
    done["total"] = sum(done.values())

    return {
        "data_root": data_root,
        "results_dir": results_dir,
        "labels_dir": labels_dir,
        "samples": samples,
        "chroms": chroms,
        "expected": expected,
        "done": done,
        "bam_ok": bam_ok,
        "call_ok": call_ok,
        "merge_ok": merge_ok,
        "label_ok": label_ok,
        "crams": crams,
        "disk": _free_gb(data_root),
        "sv_all": sv_all,
        "gt_all": gt_all,
        "cx_all": cx_all,
        "loci_all": loci_all,
        "label_rows": label_rows,
        "label_counts": label_counts,
        "label_class": label_class,
        "per_sample": per_sample,
        "chrom_calls": chrom_calls,
        "chrom_bams": chrom_bams,
        "live": _live_jobs(host_log) if host_log else [],
        "dataset": dataset,
    }


def _pct(done: int, total: int) -> str:
    if total <= 0:
        return "n/a"
    return f"{100.0 * done / total:5.1f}%"


def _fmt_counter(c: Counter[str], order: list[str] | None = None, top: int = 12) -> str:
    if not c:
        return "(none)"
    keys = order if order else [k for k, _ in c.most_common(top)]
    parts = []
    seen = set()
    for k in keys:
        if k in c:
            parts.append(f"{k}={c[k]}")
            seen.add(k)
    for k, v in c.most_common():
        if k not in seen and len(parts) < top:
            parts.append(f"{k}={v}")
            seen.add(k)
    return "  ".join(parts)


def print_report(st: dict, *, matrix: bool, by_sample: str | None) -> None:
    samples: list[str] = st["samples"]
    chroms: list[str] = st["chroms"]
    exp, done = st["expected"], st["done"]
    n_s, n_c = len(samples), len(chroms)

    print("=== ML-0 progress ===")
    print(f"data:    {st['data_root']}")
    print(f"results: {st['results_dir']}")
    print(f"labels:  {st['labels_dir']}")
    print(f"disk:    {st['disk']}")
    print(f"scope:   {n_s} samples × {n_c} chroms = {n_s * n_c} call jobs")
    print()
    print("Pipeline steps (filesystem vs expected DAG counts)")
    for key in ("fetch_region", "call", "index_haps", "merge", "label", "dataset"):
        print(f"  {key:14s}  {done[key]:5d} / {exp[key]:5d}  {_pct(done[key], exp[key])}")
    print(f"  {'TOTAL':14s}  {done['total']:5d} / {exp['total']:5d}  {_pct(done['total'], exp['total'])}")
    print(f"  CRAMs on disk: {st['crams']} / {n_s}")
    print(f"  dataset.tsv:   {'yes' if st['dataset'] else 'no'}")
    print()

    print("Calls (arcsv_out.tab)")
    print(f"  finished chrom-jobs: {done['call']} / {exp['call']}")
    print(f"  loci (rows):         {st['loci_all']:,}")
    print(f"  svtype:   {_fmt_counter(st['sv_all'], ['DEL', 'DUP', 'INS', 'INV', 'BND'])}")
    print(f"  genotype: {_fmt_counter(st['gt_all'], ['HOM', 'HET', 'REF', 'NA'])}")
    if st["cx_all"]:
        print(f"  complex:  {_fmt_counter(st['cx_all'], top=10)}")
    print()

    print("Labels (*.labeled.tsv)")
    print(f"  samples labeled: {done['label']} / {exp['label']}")
    print(f"  rows:            {st['label_rows']:,}")
    print(f"  label: {_fmt_counter(st['label_counts'], ['1', '0', 'NA'])}")
    if st["label_class"]:
        print(f"  class: {_fmt_counter(st['label_class'], ['DEL', 'DUP', 'INS', 'INV', 'BND'])}")
    print()

    print("Per chromosome (BAMs / calls across samples)")
    hdr = "  chrom   bam  call"
    print(hdr)
    for chrom in sorted(chroms, key=_chrom_sort_key):
        b = st["chrom_bams"].get(chrom, 0)
        c = st["chrom_calls"].get(chrom, 0)
        bar = "#" * (c * 20 // max(n_s, 1))
        print(f"  {chrom:5s}  {b:3d}/{n_s:<2d} {c:3d}/{n_s:<2d}  {bar}")
    print()

    if st["live"]:
        print("Live call stages (from recent host log)")
        for row in st["live"]:
            print(f"  {row}")
        print()

    if by_sample:
        name = by_sample
        if name not in st["per_sample"]:
            raise SystemExit(f"unknown sample {name!r}")
        info = st["per_sample"][name]
        print(f"Sample {name}")
        print(
            f"  bam={info['bam']}/{n_c}  call={info['call']}/{n_c}  "
            f"loci={info['loci']:,}  merged={info['merged']}  "
            f"labeled={info['labeled']}  cram={info['cram']}"
        )
        print(f"  svtype: {_fmt_counter(info['sv'], ['DEL', 'DUP', 'INS', 'INV', 'BND'])}")
        print("  chrom status: . none  b bam-only  C call done")
        line = []
        for chrom in sorted(chroms, key=_chrom_sort_key):
            if (name, chrom) in st["call_ok"]:
                mark = "C"
            elif (name, chrom) in st["bam_ok"]:
                mark = "b"
            else:
                mark = "."
            line.append(f"{chrom.replace('chr', '')}{mark}")
        print("  " + " ".join(line))
        print()

    if matrix:
        print("Matrix (rows=samples, cols=chroms; .=none b=bam C=call L=labeled-all-calls)")
        # header: abbreviated chrom numbers
        abbrev = [c.replace("chr", "") for c in sorted(chroms, key=_chrom_sort_key)]
        print("  sample       " + " ".join(f"{a:>2s}" for a in abbrev) + "  bam call lab")
        for sample in samples:
            cells = []
            for chrom in sorted(chroms, key=_chrom_sort_key):
                if (sample, chrom) in st["call_ok"]:
                    cells.append(" C")
                elif (sample, chrom) in st["bam_ok"]:
                    cells.append(" b")
                else:
                    cells.append(" .")
            info = st["per_sample"][sample]
            lab = "Y" if info["labeled"] else "n"
            print(
                f"  {sample:12s}"
                + "".join(cells)
                + f"  {info['bam']:2d}/{n_c:<2d} {info['call']:2d}/{n_c:<2d} {lab}"
            )


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--config",
        type=Path,
        default=Path("workflows/config/label.ml0_s4.yaml"),
        help="label.smk config (default: label.ml0_s4.yaml)",
    )
    p.add_argument(
        "--host-log",
        type=Path,
        default=Path("bench/logs/ml0/s4_host.log"),
        help="S4 host log for live genotype/stage lines (empty to skip)",
    )
    p.add_argument("--no-live", action="store_true", help="Skip host-log live section")
    p.add_argument("--matrix", action="store_true", help="Print sample×chrom status matrix")
    p.add_argument("--by-sample", default="", help="Detail one sample (e.g. HG01884)")
    args = p.parse_args(argv)

    if not args.config.is_file():
        raise SystemExit(f"missing config {args.config}")
    cfg = _load_yaml(args.config)
    host_log = None if args.no_live else args.host_log
    st = collect(cfg, host_log)
    print_report(st, matrix=args.matrix, by_sample=args.by_sample or None)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
