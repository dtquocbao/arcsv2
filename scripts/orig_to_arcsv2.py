"""Add the ARC-SV2 locus columns to an original ARC-SV ``arcsv_out.tab``.

The original table has ``bp`` (block boundaries) and ``bp_uncertainty`` but no
``blocks`` / ``blocks_point`` / ``bp_point_source``, which ``arcsv2 ml-label``
needs to build allele sequences. Block i spans bp[i]-bp[i+1]. An internal
junction with uncertainty u is widened to a +/- u/2 window in ``blocks`` so the
labeler's breakpoint-shift search can use it; ``blocks_point`` keeps the
reported point. Chromosome names without ``chr`` get it when ``--add-chr``.

Usage: python scripts/orig_to_arcsv2.py IN.tab OUT.tab [--add-chr]
"""

from __future__ import annotations

import csv
import sys


def convert(row: dict[str, str], add_chr: bool) -> dict[str, str]:
    bps = [int(x) for x in row["bp"].split(",") if x.strip()]
    unc = [int(x) for x in (row.get("bp_uncertainty") or "").split(",") if x.strip()]
    unc += [0] * (len(bps) - len(unc))
    point = [f"{a}-{b}" for a, b in zip(bps, bps[1:])]
    raw = []
    for i, (a, b) in enumerate(zip(bps, bps[1:])):
        lo = a + (unc[i] - unc[i] // 2) if i > 0 else a
        hi = b - unc[i + 1] // 2 if i + 1 < len(bps) - 1 else b
        raw.append(f"{lo}-{max(lo, hi)}")
    out = dict(row)
    out["blocks"] = ",".join(raw)
    out["blocks_point"] = ",".join(point)
    out["bp_point_source"] = ",".join(["mid"] * max(0, len(bps) - 2)) or "NA"
    if add_chr and not out["chrom"].startswith("chr"):
        out["chrom"] = "chr" + out["chrom"]
    return out


def main() -> None:
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    add_chr = "--add-chr" in sys.argv
    if len(args) != 2:
        raise SystemExit(__doc__)
    src, dst = args
    with open(src, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = [convert(r, add_chr) for r in reader]
        fields = list(reader.fieldnames or [])
    for col in ("blocks", "blocks_point", "bp_point_source"):
        if col not in fields:
            fields.append(col)
    with open(dst, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"{len(rows)} calls -> {dst}")


if __name__ == "__main__":
    main()
