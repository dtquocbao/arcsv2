"""Local dress rehearsal: 40 events on 2 Mb, then the label workflow.

Writes ``labels/dataset.tsv`` for the dashboard ML tab. This is not the CI
pilot and not an M2 timing run. Run it inside the Linux container.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from simulate_haplotypes import simulate

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    out = ROOT / "pilot-large"
    meta = simulate(out, length=2_000_000, n_events=40, coverage=20, seed=7)
    (out / "gaps.bed").write_text("", encoding="utf-8")
    (out / "unreliable.bed").write_text("", encoding="utf-8")
    if not meta.get("bam"):
        raise SystemExit("pilot-large needs samtools or pysam to write an indexed BAM")
    haps = ROOT / "haps"
    haps.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(out / "hap1.fa", haps / "large.hap1.fa")
    shutil.copyfile(out / "hap2.fa", haps / "large.hap2.fa")
    command = [
        "snakemake",
        "-s",
        "workflows/label.smk",
        "-j",
        "2",
        "--configfile",
        "workflows/config/label.large.yaml",
    ]
    completed = subprocess.run(command, cwd=ROOT, check=False)
    if completed.returncode != 0:
        return completed.returncode
    from check_label_pilot import main as check_pilot

    check_code = check_pilot(
        [
            str(ROOT / "labels" / "dataset.summary.json"),
            str(ROOT / "labels" / "large.labeled.tsv"),
            str(out / "truth.tsv"),
        ]
    )
    if check_code != 0:
        return check_code
    print("labels/dataset.tsv is ready. Open it with: arcsv2 dashboard --results .")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
