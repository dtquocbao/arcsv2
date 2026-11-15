#!/usr/bin/env python3
"""Compare WG ML-0 labels on chr20 for golden samples vs gate numbers."""
from collections import Counter
from pathlib import Path

samples = {"HG01884", "NA19700"}
c_all = Counter()
c_gate = Counter()
c_all_class = {}
for p in Path("/data/labels_ml0_s3").glob("*.labeled.tsv"):
    lines = p.read_text(encoding="utf-8").splitlines()
    hdr = lines[0].split("\t")
    idx = {h: i for i, h in enumerate(hdr)}
    for line in lines[1:]:
        cols = line.split("\t")
        if cols[idx["chrom"]] != "chr20":
            continue
        lab = cols[idx["label"]]
        sv = cols[idx["sv_class"]]
        c_all[lab] += 1
        c_all_class.setdefault(sv, Counter())[lab] += 1
        if cols[idx["sample"]] in samples:
            c_gate[lab] += 1

def prec(c):
    n1, n0 = c.get("1", 0), c.get("0", 0)
    den = n1 + n0
    return n1 / den if den else None, n1, n0, c.get("NA", 0)

print("chr20 all40", prec(c_all), dict(c_all))
print("chr20 golden2", prec(c_gate), dict(c_gate))
for sv, c in sorted(c_all_class.items()):
    print("chr20", sv, prec(c))
