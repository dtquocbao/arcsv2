"""Deterministic paired-end simulator for ARC-SV2 tests.

A fixed seed reproduces the SAM text byte for byte. ``--mosaic_fraction`` mixes
reference-haplotype reads into an otherwise alternate sample. ``--n_libraries``
emits distinct insert sizes and ``RG`` / ``LB`` tags.
"""

from __future__ import annotations

import argparse
import hashlib
import random
from pathlib import Path

BASES = "ACGT"

TEMPLATES = ("DEL", "DUP", "INV", "INS", "ddup", "ddup.del", "inv.del", "invddup", "dup.del")


def _dna(rng: random.Random, n: int) -> str:
    return "".join(rng.choice(BASES) for _ in range(n))


def _revcomp(seq: str) -> str:
    table = str.maketrans("ACGT", "TGCA")
    return seq.translate(table)[::-1]


def apply_template(ref: str, kind: str, start: int, end: int) -> str:
    """Apply one template to a haplotype sequence. Coordinates are half-open."""
    if kind == "DEL":
        return ref[:start] + ref[end:]
    if kind == "DUP":
        return ref[:end] + ref[start:end] + ref[end:]
    if kind == "INV":
        return ref[:start] + _revcomp(ref[start:end]) + ref[end:]
    if kind == "INS":
        return ref[:start] + "N" * (end - start) + ref[start:]
    if kind == "ddup":
        return ref[:end] + ref[start:end] + ref[start:end] + ref[end:]
    if kind == "ddup.del":
        return ref[:start] + ref[start:end] + ref[end:]
    if kind == "inv.del":
        mid = (start + end) // 2
        return ref[:start] + _revcomp(ref[start:mid]) + ref[end:]
    if kind == "invddup":
        return ref[:end] + _revcomp(ref[start:end]) + ref[end:]
    if kind == "dup.del":
        return ref[:end] + ref[start:end]
    raise ValueError(kind)


def simulate(
    out_dir: Path,
    *,
    seed: int = 3,
    length: int = 20_000,
    read_len: int = 100,
    coverage: int = 8,
    n_libraries: int = 1,
    mosaic_fraction: float = 1.0,
    template: str = "DEL",
    sv_start: int = 8_000,
    sv_end: int = 8_500,
    hom: bool = False,
) -> str:
    """Write ``ref.fa``, ``hap1.fa``, ``hap2.fa``, ``truth.tsv``, ``reads.sam``.

    Returns the SHA-256 of the SAM text.
    """
    rng = random.Random(seed)
    out_dir.mkdir(parents=True, exist_ok=True)
    ref = _dna(rng, length)
    alt = apply_template(ref, template, sv_start, sv_end)
    hap1 = alt
    hap2 = alt if hom else ref
    (out_dir / "ref.fa").write_bytes(f">chrSim\n{ref}\n".encode())
    (out_dir / "hap1.fa").write_bytes(f">chrSim\n{hap1}\n".encode())
    (out_dir / "hap2.fa").write_bytes(f">chrSim\n{hap2}\n".encode())
    (out_dir / "truth.tsv").write_bytes(
        (
            "chrom\tstart\tend\ttemplate\tgt\n"
            f"chrSim\t{sv_start}\t{sv_end}\t{template}\t{'HOM' if hom else 'HET'}\n"
        ).encode()
    )
    lines = ["@HD\tVN:1.6\tSO:unsorted", "@SQ\tSN:chrSim\tLN:" + str(length)]
    for lib in range(n_libraries):
        lines.append(f"@RG\tID:rg{lib}\tLB:lib{lib}\tSM:sim")
    insert_means = [400 + 50 * i for i in range(n_libraries)]
    n_pairs = max(1, coverage * length // (2 * read_len))
    for i in range(n_pairs):
        lib = i % n_libraries
        mean = insert_means[lib]
        use_alt = rng.random() < mosaic_fraction
        seq_src = hap1 if use_alt else hap2
        # Reads are placed on the reference coordinate for the concordant
        # background. A fraction of pairs are written as deletion-spanning
        # discordant pairs so the caller has something to cluster.
        if template == "DEL" and use_alt and i % 40 == 0:
            left = sv_start - read_len
            right = sv_end
            tlen = (right + read_len) - left
            lines.append(_sam(i, 1, left, read_len, 60, tlen, "lib", lib, seq_src, left, False))
            lines.append(_sam(i, 2, right, read_len, 60, -tlen, "lib", lib, seq_src, right, True))
            continue
        start = rng.randrange(0, max(1, length - mean - read_len))
        tlen = mean
        lines.append(_sam(i, 1, start, read_len, 60, tlen, "lib", lib, ref, start, False))
        lines.append(_sam(i, 2, start + tlen - read_len, read_len, 60, -tlen, "lib", lib, ref, start + tlen - read_len, True))
    text = "\n".join(lines) + "\n"
    data = text.encode("utf-8")
    (out_dir / "reads.sam").write_bytes(data)
    return hashlib.sha256(data).hexdigest()


def _sam(i, read, pos, rlen, mapq, tlen, _prefix, lib, seq, seq_pos, reverse) -> str:
    flag = 1 | 2 | 64 if read == 1 else 1 | 2 | 128
    if reverse:
        flag |= 16
    else:
        flag |= 32
    # Mate is reverse when this read is forward, and the reverse.
    chunk = seq[seq_pos : seq_pos + rlen]
    if len(chunk) < rlen:
        chunk = (chunk + "N" * rlen)[:rlen]
    if reverse:
        chunk = _revcomp(chunk)
    qual = "I" * rlen
    return (
        f"r{i}\t{flag}\tchrSim\t{pos + 1}\t{mapq}\t{rlen}M\t=\t{pos + 1}\t{tlen}\t{chunk}\t{qual}\t"
        f"RG:Z:rg{lib}\tLB:Z:lib{lib}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=3)
    parser.add_argument("--length", type=int, default=20_000)
    parser.add_argument("--mosaic_fraction", type=float, default=1.0)
    parser.add_argument("--n_libraries", type=int, default=1)
    parser.add_argument("--template", choices=TEMPLATES, default="DEL")
    args = parser.parse_args()
    digest = simulate(
        args.out,
        seed=args.seed,
        length=args.length,
        n_libraries=args.n_libraries,
        mosaic_fraction=args.mosaic_fraction,
        template=args.template,
    )
    print(digest)


if __name__ == "__main__":
    main()
