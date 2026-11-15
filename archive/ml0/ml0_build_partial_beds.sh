set -euo pipefail
export PYTHONPATH=
cd /work
python - <<'PY'
"""Build simple_repeats.bed and segdups.bed from already-downloaded raw files; skip blacklist."""
from pathlib import Path
import sys
sys.path.insert(0, "scripts")
from fetch_filter_beds import (
    _parse_ucsc_chrom_start_end,
    merge_intervals,
    write_bed,
    RMSK_REPEAT_CLASSES,
)

dest = Path("/data/ref")
raw = dest / "raw"
# simple repeats union
intervals = []
intervals.extend(_parse_ucsc_chrom_start_end(raw / "simpleRepeat.txt.gz", chrom_i=1, start_i=2, end_i=3))
intervals.extend(
    _parse_ucsc_chrom_start_end(
        raw / "rmsk.txt.gz",
        chrom_i=5,
        start_i=6,
        end_i=7,
        class_i=11,
        class_allow=RMSK_REPEAT_CLASSES,
    )
)
intervals.extend(_parse_ucsc_chrom_start_end(raw / "microsat.txt.gz", chrom_i=1, start_i=2, end_i=3))
merged = merge_intervals(intervals)
write_bed(dest / "simple_repeats.bed", merged)
print(f"wrote simple_repeats.bed n={len(merged)}")

# segdups: chrom/start/end at columns 1,2,3 in genomicSuperDups
intervals = _parse_ucsc_chrom_start_end(raw / "genomicSuperDups.txt.gz", chrom_i=1, start_i=2, end_i=3)
merged = merge_intervals(intervals)
write_bed(dest / "segdups.bed", merged)
print(f"wrote segdups.bed n={len(merged)}")
print("OK partial beds")
PY
ls -lh /data/ref/simple_repeats.bed /data/ref/segdups.bed