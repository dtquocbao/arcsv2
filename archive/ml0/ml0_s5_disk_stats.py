#!/usr/bin/env python3
from pathlib import Path

def dir_bytes(root: Path) -> int:
    total = 0
    for p in root.rglob("*"):
        if p.is_file():
            total += p.stat().st_size
    return total

rows = []
base = Path("/data/results_ml0_s3")
for d in sorted(base.iterdir()):
    if d.is_dir():
        rows.append((dir_bytes(d), d.name))
rows.sort()
print(f"results n={len(rows)}")
print(f"min_MB={rows[0][0]/1048576:.1f} {rows[0][1]}")
print(f"med_MB={rows[len(rows)//2][0]/1048576:.1f} {rows[len(rows)//2][1]}")
print(f"max_MB={rows[-1][0]/1048576:.1f} {rows[-1][1]}")
print(f"mean_MB={sum(t for t,_ in rows)/len(rows)/1048576:.1f}")
print(f"total_MB={sum(t for t,_ in rows)/1048576:.1f}")
lab = list(Path("/data/labels_ml0_s3").glob("*.labeled.tsv"))
print(f"labeled n={len(lab)} mean_MB={sum(p.stat().st_size for p in lab)/len(lab)/1048576:.2f}")
haps = Path("/data/haps")
if haps.is_dir():
    print(f"haps_total_GB={dir_bytes(haps)/1073741824:.1f}")
