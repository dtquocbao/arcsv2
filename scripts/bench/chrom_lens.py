#!/usr/bin/env python3
from pathlib import Path
fai = Path("/data/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai")
autos = []
for line in fai.read_text().splitlines():
    name, length, *_ = line.split("\t")
    if name.startswith("chr") and name[3:].isdigit() and 1 <= int(name[3:]) <= 22:
        autos.append((name, int(length)))
for name, length in autos:
    print(f"{name}\t{length}")
total = sum(L for _, L in autos)
chr20 = dict(autos)["chr20"]
chr1 = dict(autos)["chr1"]
print(f"autosomes\t{len(autos)}\t{total}")
print(f"total_over_chr20\t{total/chr20:.2f}")
print(f"chr1_over_chr20\t{chr1/chr20:.2f}")
print(f"user_45x_chr20\t{45*chr20}")
