#!/usr/bin/env python3
from pathlib import Path
log = Path("bench/logs/ml0/s4_host.log")
first40 = None
first39 = None
last_progress = None
for i, line in enumerate(log.open(encoding="utf-8", errors="replace"), 1):
    if "S4_PROGRESS labeled=" in line:
        last_progress = (i, line.strip())
        if first40 is None and "labeled=40/40" in line:
            first40 = (i, line.strip())
        if first39 is None and "labeled=39/40" in line:
            first39 = (i, line.strip())
print("first39", first39)
print("first40", first40)
print("last", last_progress)
# nearby START timestamp: scan backward from first40
if first40:
    lines = log.read_text(encoding="utf-8", errors="replace").splitlines()
    for j in range(first40[0] - 1, max(0, first40[0] - 50), -1):
        if "START_S4" in lines[j] or "2026-" in lines[j][:40]:
            print("near", j + 1, lines[j][:200])
            if "START_S4" in lines[j] or "S4_PROGRESS" in lines[j]:
                break
