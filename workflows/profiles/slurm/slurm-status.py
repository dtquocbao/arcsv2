#!/usr/bin/env python3
"""Map ``sacct`` state to Snakemake cluster status (success / failed / running)."""

from __future__ import annotations

import subprocess
import sys

SUCCESS = {"COMPLETED", "COMPLETING", "STAMPEDE"}
FAILED = {"FAILED", "CANCELLED", "TIMEOUT", "NODE_FAIL", "OUT_OF_MEMORY", "PREEMPTED"}
RUNNING = {"RUNNING", "PENDING", "CONFIGURING", "SUSPENDED", "REQUEUED"}


def main() -> int:
    if len(sys.argv) != 2:
        print("failed")
        return 1
    jobid = sys.argv[1].split(".")[0]
    try:
        proc = subprocess.run(
            ["sacct", "-j", jobid, "--format=State", "--noheader", "--parsable2"],
            capture_output=True,
            text=True,
            check=False,
            timeout=30,
        )
    except (OSError, subprocess.TimeoutExpired):
        print("running")
        return 0
    if proc.returncode != 0:
        print("running")
        return 0
    states = {line.strip().split("|")[0] for line in proc.stdout.splitlines() if line.strip()}
    if not states:
        print("running")
        return 0
    if states & FAILED:
        print("failed")
        return 0
    if states <= SUCCESS or states == {"COMPLETED"}:
        print("success")
        return 0
    if states & RUNNING or states - SUCCESS - FAILED:
        print("running")
        return 0
    print("success")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
