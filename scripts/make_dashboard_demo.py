"""Rebuild the dashboard's bundled demo results.

Simulates two runs of the same 2 Mb chromosome (30x and 15x), calls them
through the real SAM path, and writes one folder per run:

    python/arcsv2/dashboard/demo/<run>/arcsv_out.tab
    python/arcsv2/dashboard/demo/<run>/truth.tsv
    python/arcsv2/dashboard/demo/<run>/run.json

Reads and the reference stay in a temporary folder; only the small result
tables are kept.

    python scripts/make_dashboard_demo.py
"""

from __future__ import annotations

import json
import shutil
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "scripts"))

from simulate_haplotypes import simulate  # noqa: E402

from arcsv2 import __version__  # noqa: E402
from arcsv2.config import CallConfig  # noqa: E402
from arcsv2.pipeline import call_sam  # noqa: E402

DEMO = ROOT / "python" / "arcsv2" / "dashboard" / "demo"
RUNS = (
    ("simulated_2Mb_30x", 30),
    ("simulated_2Mb_15x", 15),
)


def build(length: int = 2_000_000, events: int = 40, seed: int = 7) -> None:
    DEMO.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        for name, coverage in RUNS:
            work = Path(tmp) / name
            meta = simulate(work, length=length, n_events=events, coverage=coverage, seed=seed)
            dest = DEMO / name
            dest.mkdir(parents=True, exist_ok=True)
            started = time.perf_counter()
            rows = call_sam(str(work / "reads.sam"), str(dest / "arcsv_out.tab"), CallConfig(), chrom="chrSim", chrom_len=length)
            seconds = time.perf_counter() - started
            shutil.copyfile(work / "truth.tsv", dest / "truth.tsv")
            run = {
                "name": name,
                "source": "simulation",
                "arcsv2_version": __version__,
                "calls": len(rows),
                "call_seconds": round(seconds, 1),
                **meta,
            }
            (dest / "run.json").write_text(json.dumps(run, indent=2) + "\n", encoding="utf-8", newline="\n")
            log = dest / "graph_log.txt"
            if log.exists():
                log.unlink()
            from arcsv2.dashboard import data as dash

            calls = dash.load_calls(dest / "arcsv_out.tab")
            truth = dash.load_truth(dest / "truth.tsv")
            events_scored, scored = dash.match_truth(calls, truth)
            summary = dash.summarize(events_scored, scored)
            print(
                f"{name}: {len(rows)} calls in {seconds:.1f} s; "
                f"{summary['detected']} detected / {summary['correct']} correct / "
                f"{summary['false_calls']} FP "
                f"({summary['recall']:.0%} detected, "
                f"{summary['genotype_concordance']:.0%} genotype concordance)"
            )


if __name__ == "__main__":
    build()
