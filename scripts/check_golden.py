#!/usr/bin/env python3
"""Compare caller outputs to the v2fix11 golden freeze (export/v2fix11).

Full chr20 check (end of each Part):
  python scripts/check_golden.py \\
    --results /data/results_v2fix11 \\
    --labels /data/labels_v2fix11

Self-test against the freeze itself:
  python scripts/check_golden.py \\
    --results export/v2fix11/results \\
    --labels export/v2fix11/labels

Fast locus check (minutes when paired with a window re-call; file compare is instant):
  python scripts/check_golden.py --fast --results /data/results_v2fix11
  python scripts/check_golden.py --fast --results export/v2fix11/results

Exit non-zero on any difference; print the first 20 differing rows.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

SAMPLES = ("HG01884", "NA19700")
MAX_DIFFS = 20

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_GOLDEN = REPO_ROOT / "export" / "v2fix11"
DEFAULT_FAST_DIR = REPO_ROOT / "tests" / "golden"


def _read_lines(path: Path) -> list[str]:
    if not path.is_file():
        raise FileNotFoundError(path)
    text = path.read_text(encoding="utf-8")
    return [ln for ln in text.splitlines() if ln != ""]


def _sorted_data_rows(path: Path) -> tuple[str | None, list[str]]:
    lines = _read_lines(path)
    if not lines:
        return None, []
    header, rows = lines[0], lines[1:]
    return header, sorted(rows)


def _summary_metrics(summary_path: Path) -> tuple[int, int, float]:
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    validated = 0
    labeled01 = 0
    for _cls, d in summary.get("classes", {}).items():
        validated += int(d.get("1", 0))
        labeled01 += int(d.get("0", 0)) + int(d.get("1", 0))
    prec = validated / labeled01 if labeled01 else 0.0
    return validated, labeled01, prec


def _print_diffs(label: str, only_a: list[str], only_b: list[str], header_mismatch: str | None) -> int:
    n = 0
    if header_mismatch:
        print(f"DIFF {label}: {header_mismatch}")
        n += 1
    for row in only_a[: max(0, MAX_DIFFS - n)]:
        print(f"DIFF {label} only_in_candidate: {row}")
        n += 1
    remain = MAX_DIFFS - n
    for row in only_b[: max(0, remain)]:
        print(f"DIFF {label} only_in_golden: {row}")
        n += 1
    extra = (len(only_a) + len(only_b) + (1 if header_mismatch else 0)) - n
    if extra > 0:
        print(f"... {extra} more differing rows for {label} (showing first {MAX_DIFFS})")
    return len(only_a) + len(only_b) + (1 if header_mismatch else 0)


def _diff_sorted_tables(label: str, cand: Path, gold: Path) -> int:
    try:
        ch, cr = _sorted_data_rows(cand)
        gh, gr = _sorted_data_rows(gold)
    except FileNotFoundError as exc:
        print(f"DIFF {label}: missing file {exc}")
        return 1
    header_mismatch = None
    if ch != gh:
        header_mismatch = f"header mismatch\n  candidate: {ch}\n  golden:    {gh}"
    only_c = sorted(set(cr) - set(gr))
    only_g = sorted(set(gr) - set(cr))
    # Also catch multiplicity differences via multiset compare
    from collections import Counter

    cc, gc = Counter(cr), Counter(gr)
    if cc != gc and not only_c and not only_g:
        for row, n in sorted((cc - gc).items()):
            only_c.extend([row] * n)
        for row, n in sorted((gc - cc).items()):
            only_g.extend([row] * n)
    if header_mismatch or only_c or only_g:
        return _print_diffs(label, only_c, only_g, header_mismatch)
    print(f"OK {label}: {len(cr)} data rows")
    return 0


def check_full(results: Path, labels: Path, golden: Path) -> int:
    diffs = 0
    gold_results = golden / "results"
    gold_labels = golden / "labels"
    for sample in SAMPLES:
        diffs += _diff_sorted_tables(
            f"merged.tab/{sample}",
            results / sample / "merged.tab",
            gold_results / sample / "merged.tab",
        )
    diffs += _diff_sorted_tables(
        "dataset.tsv",
        labels / "dataset.tsv",
        gold_labels / "dataset.tsv",
    )

    cand_summary = labels / "dataset.summary.json"
    gold_summary = gold_labels / "dataset.summary.json"
    try:
        cv, c01, cp = _summary_metrics(cand_summary)
        gv, g01, gp = _summary_metrics(gold_summary)
    except FileNotFoundError as exc:
        print(f"DIFF summary: missing file {exc}")
        return diffs + 1
    if (cv, c01) != (gv, g01) or abs(cp - gp) > 1e-12:
        print(
            "DIFF validated/precision summary:\n"
            f"  candidate: validated={cv} labeled_0_or_1={c01} precision={cp:.4f} ({cv}/{c01})\n"
            f"  golden:    validated={gv} labeled_0_or_1={g01} precision={gp:.4f} ({gv}/{g01})"
        )
        diffs += 1
    else:
        print(f"OK summary: validated={cv} precision={cp:.4f} ({cv}/{c01})")
    return diffs


def _near_rows(path: Path, center: int, slop: int) -> list[str]:
    lines = _read_lines(path)
    if not lines:
        return []
    header = lines[0]
    kept = [header]
    for line in lines[1:]:
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        start, end = int(parts[1]), int(parts[2])
        if abs(start - center) <= slop or abs(end - center) <= slop or start <= center <= end:
            kept.append(line)
    return kept


def check_fast(results: Path, golden_dir: Path) -> int:
    manifest_path = golden_dir / "manifest.json"
    if not manifest_path.is_file():
        print(f"DIFF fast: missing {manifest_path}")
        return 1
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    diffs = 0
    for entry in manifest:
        sample = entry["sample"]
        center = int(entry["center"])
        slop = int(entry.get("slop", 500))
        gold_path = golden_dir / entry["file"]
        cand_path = results / sample / "merged.tab"
        try:
            gold_lines = _read_lines(gold_path)
            cand_lines = _near_rows(cand_path, center, slop)
        except FileNotFoundError as exc:
            print(f"DIFF fast {entry['file']}: missing file {exc}")
            diffs += 1
            continue
        # Compare as sorted data rows (header must match)
        g_header, g_rows = (gold_lines[0], sorted(gold_lines[1:])) if gold_lines else (None, [])
        c_header, c_rows = (cand_lines[0], sorted(cand_lines[1:])) if cand_lines else (None, [])
        header_mismatch = None
        if g_header != c_header:
            header_mismatch = f"header mismatch\n  candidate: {c_header}\n  golden:    {g_header}"
        only_c = sorted(set(c_rows) - set(g_rows))
        only_g = sorted(set(g_rows) - set(c_rows))
        label = f"fast/{entry['file']}"
        if header_mismatch or only_c or only_g:
            diffs += _print_diffs(label, only_c, only_g, header_mismatch)
        else:
            print(f"OK {label}: {len(c_rows)} calls near {center} ({entry.get('note', '')})")
    return diffs


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        "--results",
        type=Path,
        default=None,
        help="Results root with {sample}/merged.tab (default: <golden>/results)",
    )
    p.add_argument(
        "--labels",
        type=Path,
        default=None,
        help="Labels dir with dataset.tsv + dataset.summary.json (default: <golden>/labels)",
    )
    p.add_argument(
        "--golden",
        type=Path,
        default=DEFAULT_GOLDEN,
        help=f"Golden export root (default: {DEFAULT_GOLDEN})",
    )
    p.add_argument(
        "--fast",
        action="store_true",
        help="Compare only tests/golden locus windows instead of full chr20",
    )
    p.add_argument(
        "--fast-dir",
        type=Path,
        default=DEFAULT_FAST_DIR,
        help=f"Fast golden dir (default: {DEFAULT_FAST_DIR})",
    )
    args = p.parse_args(argv)

    golden = args.golden.resolve()
    if not golden.is_dir():
        print(f"ERROR: golden export not found: {golden}", file=sys.stderr)
        print("Expected export/v2fix11 under the repo (/work) or pass --golden.", file=sys.stderr)
        return 2

    results = (args.results or (golden / "results")).resolve()
    labels = (args.labels or (golden / "labels")).resolve()

    if args.fast:
        diffs = check_fast(results, args.fast_dir.resolve())
        mode = "fast"
    else:
        diffs = check_full(results, labels, golden)
        mode = "full"

    if diffs:
        print(f"FAIL {mode} golden check: {diffs} difference group(s)")
        return 1
    print(f"PASS {mode} golden check")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
