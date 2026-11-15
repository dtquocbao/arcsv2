"""Haplotype-aware paired-end simulator.

Each haplotype is a list of segments taken from the reference (forward or
reverse-complemented) or novel sequence. Read pairs are sampled from the
haplotypes and written as SAM against the reference, with soft clips at
junctions and unmapped mates inside novel sequence. The output exercises the
same SAM path as real data (``arcsv2.samio.load_region``).

Truth events use the types the caller reports: DEL, INV, DUP, INS, and the
complex ``inv.del`` (A B C D -> A C' D).

    python scripts/simulate_haplotypes.py --out demo_run --length 2000000 --events 40
"""

from __future__ import annotations

import argparse
import json
import random
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

EVENT_TYPES = ("DEL", "INV", "DUP", "INS", "inv.del")
_COMP = str.maketrans("ACGT", "TGCA")


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


@dataclass(frozen=True)
class Event:
    chrom: str
    start: int
    end: int
    kind: str
    hom: bool
    size: int

    def as_row(self) -> dict[str, str]:
        return {
            "chrom": self.chrom,
            "start": str(self.start),
            "end": str(self.end),
            "type": self.kind,
            "gt": "HOM" if self.hom else "HET",
            "size": str(self.size),
        }


# A segment is (ref_start, ref_end, strand) or (None, novel_sequence, 0).
Segment = tuple


def _segments(ref_len: int, events: list[Event], novel: dict[int, str]) -> list[Segment]:
    out: list[Segment] = []
    cursor = 0
    for ev in sorted(events, key=lambda e: e.start):
        out.append((cursor, ev.start, 1))
        if ev.kind == "DEL":
            cursor = ev.end
        elif ev.kind == "INV":
            out.append((ev.start, ev.end, -1))
            cursor = ev.end
        elif ev.kind == "DUP":
            out.append((ev.start, ev.end, 1))
            out.append((ev.start, ev.end, 1))
            cursor = ev.end
        elif ev.kind == "INS":
            out.append((None, novel[ev.start], 0))
            cursor = ev.start
        elif ev.kind == "inv.del":
            mid = (ev.start + ev.end) // 2
            out.append((mid, ev.end, -1))
            cursor = ev.end
        else:
            raise ValueError(ev.kind)
    out.append((cursor, ref_len, 1))
    return [seg for seg in out if seg[0] is None or seg[1] > seg[0]]


def _hap_sequence(ref: str, segs: list[Segment]) -> tuple[str, list[tuple[int, int, Segment]]]:
    parts: list[str] = []
    bounds: list[tuple[int, int, Segment]] = []
    offset = 0
    for seg in segs:
        if seg[0] is None:
            piece = seg[1]
        else:
            piece = ref[seg[0] : seg[1]]
            if seg[2] < 0:
                piece = revcomp(piece)
        parts.append(piece)
        bounds.append((offset, offset + len(piece), seg))
        offset += len(piece)
    return "".join(parts), bounds


def _align(bounds: list[tuple[int, int, Segment]], a0: int, a1: int, min_aligned: int = 30):
    """Reference position, CIGAR and strand of hap interval ``[a0, a1)``.

    The longest overlap with a reference segment is the alignment; the rest is
    soft-clipped. None when no reference segment covers ``min_aligned`` bases.
    """
    best = None
    for h0, h1, seg in bounds:
        if seg[0] is None:
            continue
        lo, hi = max(a0, h0), min(a1, h1)
        if hi - lo <= 0:
            continue
        if best is None or hi - lo > best[1] - best[0]:
            best = (lo, hi, h0, seg)
    if best is None or best[1] - best[0] < min_aligned:
        return None
    lo, hi, h0, seg = best
    left, right = lo - a0, a1 - hi
    if seg[2] > 0:
        pos = seg[0] + (lo - h0)
    else:
        pos = seg[1] - (hi - h0)
        left, right = right, left
    cigar = (f"{left}S" if left else "") + f"{hi - lo}M" + (f"{right}S" if right else "")
    return pos, cigar, seg[2]


def place_events(
    length: int,
    n_events: int,
    rng: random.Random,
    chrom: str,
    sizes: tuple[int, ...] = (300, 800, 1500, 3000),
) -> list[Event]:
    step = length // (n_events + 1)
    events: list[Event] = []
    for i in range(n_events):
        kind = EVENT_TYPES[i % len(EVENT_TYPES)]
        size = rng.choice(sizes)
        start = step * (i + 1)
        end = start if kind == "INS" else start + size
        events.append(Event(chrom, start, end, kind, rng.random() < 0.5, size))
    return events


def _write_indexed_bam(sam_path: Path, bam_path: Path) -> str | None:
    """Sort ``sam_path`` into ``bam_path`` and write ``bam_path.bai``.

    Uses pysam when it is installed, otherwise ``samtools`` on PATH.
    Returns None when neither is available.
    """
    try:
        import pysam
    except ImportError:
        pysam = None
    if pysam is not None:
        pysam.sort("-o", str(bam_path), str(sam_path))
        pysam.index(str(bam_path))
        return str(bam_path)
    if shutil.which("samtools"):
        subprocess.check_call(["samtools", "sort", "-o", str(bam_path), str(sam_path)])
        subprocess.check_call(["samtools", "index", str(bam_path)])
        return str(bam_path)
    return None


def simulate(
    out_dir: Path,
    *,
    length: int = 2_000_000,
    n_events: int = 40,
    coverage: int = 30,
    read_len: int = 100,
    insert_mean: int = 400,
    insert_sd: int = 40,
    seed: int = 7,
    chrom: str = "chrSim",
) -> dict:
    """Write reference, both haplotypes, truth, SAM, and an indexed BAM when possible."""
    rng = random.Random(seed)
    out_dir.mkdir(parents=True, exist_ok=True)
    ref = "".join(rng.choice("ACGT") for _ in range(length))
    events = place_events(length, n_events, rng, chrom)
    novel = {ev.start: "".join(rng.choice("ACGT") for _ in range(ev.size)) for ev in events if ev.kind == "INS"}
    haps = [
        _segments(length, events, novel),
        _segments(length, [ev for ev in events if ev.hom], novel),
    ]
    lines = ["@HD\tVN:1.6\tSO:unsorted", f"@SQ\tSN:{chrom}\tLN:{length}", "@RG\tID:rg0\tLB:lib0\tSM:sim"]
    n_pairs = 0
    hap_seqs: list[str] = []
    for segs in haps:
        seq, bounds = _hap_sequence(ref, segs)
        hap_seqs.append(seq)
        for _ in range(coverage * len(seq) // (2 * read_len) // 2):
            frag = max(2 * read_len, int(rng.gauss(insert_mean, insert_sd)))
            x = rng.randrange(0, len(seq) - frag)
            read1 = seq[x : x + read_len]
            read2 = revcomp(seq[x + frag - read_len : x + frag])
            a1 = _align(bounds, x, x + read_len)
            a2 = _align(bounds, x + frag - read_len, x + frag)
            if a1 is None and a2 is None:
                continue
            # Read 1 is forward on the haplotype, read 2 reverse. A reversed
            # segment flips the strand on the reference.
            rev1 = a1 is not None and a1[2] < 0
            rev2 = a2 is not None and a2[2] > 0
            recs = []
            for me, mate, read, rev, mate_rev, first in (
                (a1, a2, read1, rev1, rev2, True),
                (a2, a1, read2, rev2, rev1, False),
            ):
                flag = 0x1 | (0x40 if first else 0x80)
                if me is None:
                    flag |= 0x4
                    pos, cigar, mapq = mate[0], "*", 0
                else:
                    pos, cigar, mapq = me[0], me[1], 60
                    if rev:
                        flag |= 0x10
                if mate is None:
                    flag |= 0x8
                    mpos = pos
                else:
                    mpos = mate[0]
                    if mate_rev:
                        flag |= 0x20
                if me is not None and mate is not None:
                    span = max(pos, mpos) + read_len - min(pos, mpos)
                    tlen = span if pos <= mpos else -span
                else:
                    tlen = 0
                stored = revcomp(read) if (flag & 0x10) else read
                recs.append(
                    f"p{n_pairs}\t{flag}\t{chrom}\t{pos + 1}\t{mapq}\t{cigar}\t=\t{mpos + 1}\t{tlen}\t"
                    f"{stored}\t{'I' * read_len}\tRG:Z:rg0"
                )
            lines.extend(recs)
            n_pairs += 1
    (out_dir / "ref.fa").write_text(f">{chrom}\n{ref}\n", encoding="utf-8", newline="\n")
    for name, seq in zip(("hap1.fa", "hap2.fa"), hap_seqs, strict=True):
        (out_dir / name).write_text(f">{chrom}\n{seq}\n", encoding="utf-8", newline="\n")
    sam_path = out_dir / "reads.sam"
    sam_path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    bam = _write_indexed_bam(sam_path, out_dir / "reads.bam")
    header = ["chrom", "start", "end", "type", "gt", "size"]
    rows = ["\t".join(header)] + ["\t".join(ev.as_row()[k] for k in header) for ev in events]
    (out_dir / "truth.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8", newline="\n")
    meta = {
        "simulator": "scripts/simulate_haplotypes.py",
        "length": length,
        "events": n_events,
        "coverage": coverage,
        "read_len": read_len,
        "insert_mean": insert_mean,
        "insert_sd": insert_sd,
        "seed": seed,
        "pairs": n_pairs,
        "bam": bam,
    }
    (out_dir / "simulation.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8", newline="\n")
    return meta


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--length", type=int, default=2_000_000)
    parser.add_argument("--events", type=int, default=40)
    parser.add_argument("--coverage", type=int, default=30)
    parser.add_argument("--seed", type=int, default=7)
    args = parser.parse_args()
    meta = simulate(args.out, length=args.length, n_events=args.events, coverage=args.coverage, seed=args.seed)
    print(json.dumps(meta))


if __name__ == "__main__":
    main()
