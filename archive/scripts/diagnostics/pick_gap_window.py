#!/usr/bin/env python3
"""Pick a chr20 window and compare original vs ARC-SV2 labeled calls.

Diagnostic only — does not modify the caller. Writes a TSV report under
``scratch/gap_analysis/``.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def load_labeled(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def as_int(value: str | None, default: int = 0) -> int:
    if value is None or value == "" or value.upper() == "NA":
        return default
    try:
        return int(float(value))
    except ValueError:
        return default


def as_float(value: str | None, default: float = float("nan")) -> float:
    if value is None or value == "" or value.upper() == "NA":
        return default
    try:
        return float(value)
    except ValueError:
        return default


def label_of(row: dict) -> str:
    raw = row.get("label", "")
    if raw is None or raw == "" or raw.lower() == "nan":
        return "NA"
    try:
        return str(int(float(raw)))
    except ValueError:
        return str(raw)


def overlaps(a_start: int, a_end: int, b_start: int, b_end: int, slop: int = 100) -> bool:
    return not (a_end + slop < b_start or b_end + slop < a_start)


def pick_window(orig: list[dict], v2: list[dict], chrom: str, width: int) -> tuple[int, int, dict]:
    """Slide windows; maximize orig-only true + v2-only false."""
    orig_true = [
        r
        for r in orig
        if r.get("chrom") == chrom and label_of(r) == "1"
    ]
    v2_false = [
        r
        for r in v2
        if r.get("chrom") == chrom and label_of(r) == "0"
    ]
    starts = sorted({as_int(r["minbp"]) // width * width for r in orig_true + v2_false})
    best = None
    for start in starts:
        end = start + width
        o_only = []
        for r in orig_true:
            lo, hi = as_int(r["minbp"]), as_int(r["maxbp"])
            if hi < start or lo > end:
                continue
            if not any(
                overlaps(lo, hi, as_int(v["minbp"]), as_int(v["maxbp"]))
                for v in v2
                if v.get("chrom") == chrom
            ):
                o_only.append(r)
        v_only_false = []
        for r in v2_false:
            lo, hi = as_int(r["minbp"]), as_int(r["maxbp"])
            if hi < start or lo > end:
                continue
            if not any(
                overlaps(lo, hi, as_int(o["minbp"]), as_int(o["maxbp"]))
                for o in orig
                if o.get("chrom") == chrom
            ):
                v_only_false.append(r)
        score = (len(o_only), len(v_only_false), -start)
        if best is None or score > best[0]:
            best = (score, start, end, o_only, v_only_false)
    assert best is not None
    _, start, end, o_only, v_only_false = best
    return start, end, {"orig_only_true": o_only, "v2_only_false": v_only_false}


def summarize(rows: list[dict]) -> dict:
    by_label: dict[str, int] = defaultdict(int)
    by_class: dict[str, int] = defaultdict(int)
    sr0 = 0
    for r in rows:
        by_label[label_of(r)] += 1
        by_class[r.get("sv_class") or r.get("svtype") or "?"] += 1
        if as_int(r.get("sr_support")) == 0:
            sr0 += 1
    return {"n": len(rows), "by_label": dict(by_label), "by_class": dict(by_class), "sr0": sr0}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True, help="arcsv2 repo root")
    parser.add_argument("--sample", default="HG01884")
    parser.add_argument("--chrom", default="chr20")
    parser.add_argument("--width", type=int, default=3_000_000)
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()
    root: Path = args.root
    out = args.out or (root / "scratch" / "gap_analysis")
    out.mkdir(parents=True, exist_ok=True)

    orig = load_labeled(root / "export" / "orig" / "labels" / f"{args.sample}.labeled.tsv")
    v2 = load_labeled(root / "export" / "labels" / f"{args.sample}.labeled.tsv")
    start, end, picks = pick_window(orig, v2, args.chrom, args.width)

    report = out / f"{args.sample}_{args.chrom}_{start}_{end}_window.txt"
    lines = [
        f"sample={args.sample} chrom={args.chrom} window={start}-{end} width={args.width}",
        f"orig labeled: {summarize(orig)}",
        f"v2 labeled: {summarize(v2)}",
        f"orig-only true in window: {len(picks['orig_only_true'])}",
        f"v2-only false in window: {len(picks['v2_only_false'])}",
        "",
        "=== original-only true (up to 10) ===",
    ]
    for r in picks["orig_only_true"][:10]:
        lines.append(
            f"  {r['svtype']}/{r.get('sv_class','')} {r['minbp']}-{r['maxbp']} "
            f"sr={r.get('sr_support')} pe={r.get('pe_support')} gt={r.get('gt')} "
            f"len={r.get('len_affected')} score={r.get('score_vs_ref')} reason={r.get('reason')}"
        )
    lines.append("")
    lines.append("=== ARC-SV2-only false (up to 10) ===")
    for r in picks["v2_only_false"][:10]:
        lines.append(
            f"  {r['svtype']}/{r.get('sv_class','')} {r['minbp']}-{r['maxbp']} "
            f"sr={r.get('sr_support')} pe={r.get('pe_support')} gt={r.get('gt')} "
            f"len={r.get('len_affected')} score={r.get('score_vs_ref')} reason={r.get('reason')}"
        )

    # Global mismatch lists for traces
    orig_true = [r for r in orig if label_of(r) == "1"]
    v2_all = v2
    orig_only_true = []
    for r in orig_true:
        lo, hi = as_int(r["minbp"]), as_int(r["maxbp"])
        if not any(
            overlaps(lo, hi, as_int(v["minbp"]), as_int(v["maxbp"])) for v in v2_all if v.get("chrom") == args.chrom
        ):
            orig_only_true.append(r)
    v2_false = [r for r in v2 if label_of(r) == "0"]
    v2_only_false = []
    for r in v2_false:
        lo, hi = as_int(r["minbp"]), as_int(r["maxbp"])
        if not any(
            overlaps(lo, hi, as_int(o["minbp"]), as_int(o["maxbp"])) for o in orig if o.get("chrom") == args.chrom
        ):
            v2_only_false.append(r)

    # Prefer high-SR orig-only true and high-PE v2-only false (matches gap signals)
    orig_only_true.sort(key=lambda r: (-as_int(r.get("sr_support")), as_int(r.get("pe_support"))))
    v2_only_false.sort(key=lambda r: (as_int(r.get("sr_support")), -as_int(r.get("pe_support"))))

    traces = out / f"{args.sample}_call_trace_candidates.tsv"
    with traces.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "set",
                "chrom",
                "minbp",
                "maxbp",
                "svtype",
                "sv_class",
                "gt",
                "sr_support",
                "pe_support",
                "len_affected",
                "score_vs_ref",
                "reason",
                "bp",
            ]
        )
        for r in orig_only_true[:5]:
            writer.writerow(
                [
                    "orig_only_true",
                    r.get("chrom"),
                    r.get("minbp"),
                    r.get("maxbp"),
                    r.get("svtype"),
                    r.get("sv_class"),
                    r.get("gt"),
                    r.get("sr_support"),
                    r.get("pe_support"),
                    r.get("len_affected"),
                    r.get("score_vs_ref"),
                    r.get("reason"),
                    r.get("bp"),
                ]
            )
        for r in v2_only_false[:5]:
            writer.writerow(
                [
                    "v2_only_false",
                    r.get("chrom"),
                    r.get("minbp"),
                    r.get("maxbp"),
                    r.get("svtype"),
                    r.get("sv_class"),
                    r.get("gt"),
                    r.get("sr_support"),
                    r.get("pe_support"),
                    r.get("len_affected"),
                    r.get("score_vs_ref"),
                    r.get("reason"),
                    r.get("bp"),
                ]
            )

    lines += [
        "",
        f"global orig-only true: {len(orig_only_true)} (median sr of top5: "
        f"{[as_int(r.get('sr_support')) for r in orig_only_true[:5]]})",
        f"global v2-only false: {len(v2_only_false)} (sr/pe of top5: "
        f"{[(as_int(r.get('sr_support')), as_int(r.get('pe_support'))) for r in v2_only_false[:5]]})",
        f"wrote {traces}",
    ]
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(report.read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
