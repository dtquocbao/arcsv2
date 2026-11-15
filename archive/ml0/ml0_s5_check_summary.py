#!/usr/bin/env python3
import json
from pathlib import Path

d = json.loads(Path("/data/ml/dataset_v0.summary.json").read_text())
print(json.dumps(d["precision_by_superpopulation"], indent=2))
print("n_rows", d["n_rows"])
print("label_totals", d["label_totals"])
