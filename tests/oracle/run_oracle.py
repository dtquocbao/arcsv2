"""Run the original ARC-SV on a BAM region and store ``arcsv_out.tab``.

Requires ``ARCSV_BIN`` (the oracle environment's ``arcsv``) plus ``-i -R -G -r -o``.
"""

from __future__ import annotations

import argparse
import os
import subprocess
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--bam", required=True)
    parser.add_argument("-R", "--reference", required=True)
    parser.add_argument("-G", "--gaps", required=True)
    parser.add_argument("-r", "--region", required=True)
    parser.add_argument("-o", "--out", required=True)
    args = parser.parse_args()
    binary = os.environ.get("ARCSV_BIN", "arcsv")
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [binary, "call", "-i", args.bam, "-r", args.region, "-R", args.reference, "-G", args.gaps, "-o", str(out)],
        check=True,
    )
    tab = out / "arcsv_out.tab"
    if not tab.exists():
        raise SystemExit(f"oracle did not write {tab}")
    print(tab)


if __name__ == "__main__":
    main()
