#!/usr/bin/env python3
from collections import Counter
from pathlib import Path

for sample in ("HG01884", "NA19700"):
    p = Path(f"/data/labels_ml0_s3/{sample}.labeled.tsv")
    lines = p.read_text(encoding="utf-8").splitlines()
    hdr = lines[0].split("\t")
    idx = {h: i for i, h in enumerate(hdr)}
    c = Counter()
    n20 = 0
    for line in lines[1:]:
        cols = line.split("\t")
        if cols[idx["chrom"]] != "chr20":
            continue
        n20 += 1
        c[cols[idx["label"]]] += 1
    print(sample, "chr20_rows", n20, dict(c))
    mt = Path(f"/data/results_ml0_s3/{sample}/merged.tab")
    mlines = mt.read_text(encoding="utf-8").splitlines()
    n_m = sum(1 for ln in mlines[1:] if ln.startswith("chr20\t"))
    print(sample, "merged_chr20", n_m, "merged_total", len(mlines) - 1)
