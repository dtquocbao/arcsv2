"""Compare labeled call tables (orig vs ARC-SV2) for RESULTS.md rows."""
from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path


def load_rows(path: Path) -> list[dict[str, str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines:
        return []
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:] if line]


def sv_class(row: dict[str, str]) -> str:
    ct = row.get("complextype", "")
    if ct and ct not in ("", "NA"):
        return ct.split(",")[0]
    return (row.get("svtype") or "UNK").split(",")[0]


def key(row: dict[str, str]) -> tuple:
    return (row["chrom"], int(row["minbp"]), int(row["maxbp"]), sv_class(row))


def summarize(rows: list[dict[str, str]], name: str) -> dict:
    labeled = [r for r in rows if r.get("label") in ("0", "1")]
    by_class: dict[str, Counter] = defaultdict(Counter)
    for r in labeled:
        by_class[sv_class(r)][r["label"]] += 1
    precision = {}
    for kind, counts in sorted(by_class.items()):
        tot = counts["0"] + counts["1"]
        precision[kind] = {
            "n": tot,
            "tp": counts["1"],
            "fp": counts["0"],
            "precision": (counts["1"] / tot) if tot else None,
        }
    sr0 = sum(1 for r in rows if (r.get("sr_support") or "0").split(",")[0] in ("0", "0.0"))
    return {
        "name": name,
        "n_calls": len(rows),
        "n_labeled": len(labeled),
        "n_validated": sum(1 for r in labeled if r["label"] == "1"),
        "sr0": sr0,
        "sr0_frac": sr0 / len(rows) if rows else None,
        "precision_by_class": precision,
    }


def overlap(a: list[dict[str, str]], b: list[dict[str, str]], slop: int = 200) -> dict:
    """Approximate shared calls by chrom + class + interval within slop."""
    bins: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for r in b:
        bins[(r["chrom"], sv_class(r))].append(r)
    shared_a = []
    only_a = []
    used = set()
    for r in a:
        cands = bins.get((r["chrom"], sv_class(r)), [])
        start, end = int(r["minbp"]), int(r["maxbp"])
        hit = None
        for i, o in enumerate(cands):
            oid = id(o)
            if oid in used:
                continue
            os, oe = int(o["minbp"]), int(o["maxbp"])
            if abs(start - os) <= slop and abs(end - oe) <= slop:
                hit = oid
                break
            if start <= oe + slop and os <= end + slop:
                hit = oid
                break
        if hit is not None:
            used.add(hit)
            shared_a.append(r)
        else:
            only_a.append(r)
    only_b = [r for r in b if id(r) not in used]
    def prec(rows):
        lab = [r for r in rows if r.get("label") in ("0", "1")]
        if not lab:
            return None
        return sum(1 for r in lab if r["label"] == "1") / len(lab)
    return {
        "shared": len(shared_a),
        "shared_precision": prec(shared_a),
        "a_only": len(only_a),
        "a_only_precision": prec(only_a),
        "b_only": len(only_b),
        "b_only_precision": prec(only_b),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--orig", type=Path, required=True, help="orig labeled TSV or dataset")
    ap.add_argument("--v2", type=Path, required=True, help="v2 labeled TSV or dataset")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    orig = load_rows(args.orig)
    v2 = load_rows(args.v2)
    payload = {
        "orig": summarize(orig, "orig"),
        "v2": summarize(v2, "arcsv2"),
        "overlap": overlap(orig, v2),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
