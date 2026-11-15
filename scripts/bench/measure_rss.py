#!/usr/bin/env python3
"""Run a command; print SAMPLE=... elapsed_sec=... max_rss_kb=... rc=...

Peak RSS = max VmRSS among the process and its descendants (1 Hz), matching
prior gate_lean style where fork workers share CoW pages with the parent.
"""
from __future__ import annotations

import os
import subprocess
import sys
import time
from pathlib import Path


def _proc_rss_and_ppid(pid: int) -> tuple[int, int | None]:
    try:
        st = Path(f"/proc/{pid}/status").read_text().splitlines()
    except (FileNotFoundError, PermissionError):
        return 0, None
    rss = 0
    ppid = None
    for line in st:
        if line.startswith("VmRSS:"):
            rss = int(line.split()[1])
        elif line.startswith("PPid:"):
            ppid = int(line.split()[1])
    return rss, ppid


def _tree_peak_rss_kb(root_pid: int) -> int:
    # Build pid -> (rss, ppid) for live procs; take max RSS in the subtree.
    table: dict[int, tuple[int, int | None]] = {}
    try:
        for entry in Path("/proc").iterdir():
            if not entry.name.isdigit():
                continue
            pid = int(entry.name)
            table[pid] = _proc_rss_and_ppid(pid)
    except FileNotFoundError:
        return _proc_rss_and_ppid(root_pid)[0]

    children: dict[int, list[int]] = {}
    for pid, (_rss, ppid) in table.items():
        if ppid is None:
            continue
        children.setdefault(ppid, []).append(pid)

    peak = 0
    stack = [root_pid]
    seen: set[int] = set()
    while stack:
        pid = stack.pop()
        if pid in seen:
            continue
        seen.add(pid)
        rss, _ = table.get(pid, (0, None))
        peak = max(peak, rss)
        stack.extend(children.get(pid, ()))
    return peak


def main() -> int:
    if len(sys.argv) < 3 or sys.argv[2] != "--":
        print("usage: measure_rss.py SAMPLE -- cmd...", file=sys.stderr)
        return 2
    sample = sys.argv[1]
    cmd = sys.argv[3:]
    t0 = time.time()
    proc = subprocess.Popen(cmd)
    peak = 0
    while proc.poll() is None:
        peak = max(peak, _tree_peak_rss_kb(proc.pid))
        time.sleep(1.0)
    peak = max(peak, _tree_peak_rss_kb(proc.pid))
    elapsed = time.time() - t0
    rc = proc.returncode if proc.returncode is not None else -1
    line = f"SAMPLE={sample} elapsed_sec={elapsed:.1f} max_rss_kb={peak} rc={rc}"
    print(line)
    out = os.environ.get("MEASURE_OUT")
    if out:
        Path(out).write_text(line + "\n")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
