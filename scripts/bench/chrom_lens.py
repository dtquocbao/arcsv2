#!/usr/bin/env python3
"""Print GRCh38 autosome lengths from the analysis-set fai."""

from pathlib import Path

# GRCh38 analysis set, primary chromosomes. These are the lengths the fai
# printout reports for chr1 and chr20 (chr1 248956422, chr20 64444167).
GRCH38_BP = {
    "chr1": 248_956_422,
    "chr20": 64_444_167,
}


def main() -> None:
    fai = Path("/data/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai")
    autos = []
    for line in fai.read_text().splitlines():
        name, length, *_ = line.split("\t")
        if name.startswith("chr") and name[3:].isdigit() and 1 <= int(name[3:]) <= 22:
            autos.append((name, int(length)))
    for name, length in autos:
        print(f"{name}\t{length}")
    total = sum(length for _, length in autos)
    by_name = dict(autos)
    chr20 = by_name["chr20"]
    chr1 = by_name["chr1"]
    print(f"autosomes\t{len(autos)}\t{total}")
    print(f"total_over_chr20\t{total / chr20:.2f}")
    print(f"chr1_over_chr20\t{chr1 / chr20:.2f}")
    print(f"user_45x_chr20\t{45 * chr20}")


if __name__ == "__main__":
    main()
