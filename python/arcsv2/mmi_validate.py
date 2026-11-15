"""Haplotype ``.mmi`` validation before labeling (U4 / ML-0 HG00140)."""

from __future__ import annotations

from pathlib import Path

from arcsv2.input_checks import InputCheckError

# Healthy HPRC indexes are ~5.5–7× the compressed FASTA size; corrupt HG00140 was ~4.6×.
MMI_BYTES_PER_FASTA_GZ = 5.2
MMI_BYTES_PER_FASTA = 0.85

SMOKE_CHROM = "chr20"
SMOKE_START = 30_000_000
SMOKE_END = 30_001_000

NO_HIT_MAX_FRACTION = 0.02


def min_expected_mmi_bytes(fasta: Path) -> int:
    size = fasta.stat().st_size
    if str(fasta).endswith(".gz"):
        return int(size * MMI_BYTES_PER_FASTA_GZ)
    return int(size * MMI_BYTES_PER_FASTA)


def validate_mmi_size(mmi: Path, fasta: Path, *, hap_name: str) -> None:
    if not mmi.is_file():
        raise InputCheckError(
            f"Haplotype index missing for {hap_name}: {mmi}",
            f"Rebuild with: minimap2 -d {mmi} {fasta}",
        )
    if not fasta.is_file():
        raise InputCheckError(f"Haplotype FASTA missing for {hap_name}: {fasta}", "Check --hap1/--hap2 paths.")
    expected = min_expected_mmi_bytes(fasta)
    actual = mmi.stat().st_size
    if actual < expected:
        raise InputCheckError(
            f"{hap_name} index looks too small ({actual / (1024**3):.2f} GiB; "
            f"expected at least {expected / (1024**3):.2f} GiB for {fasta.name}). "
            "Stale or corrupt indexes cause mass no_hit labels.",
            f"Rebuild with: minimap2 -d {mmi} {fasta}",
        )


def mmi_smoke_aligns(mmi: Path, reference: Path) -> bool:
    try:
        import mappy
        import pysam
    except ImportError:
        return True
    try:
        with pysam.FastaFile(str(reference)) as fa:
            if SMOKE_CHROM not in fa.references:
                return True
            probe = fa.fetch(SMOKE_CHROM, SMOKE_START, SMOKE_END)
    except (OSError, ValueError):
        return True
    if len(probe) < 100:
        return True
    aligner = mappy.Aligner(str(mmi), n_threads=1)
    if not aligner:
        return False
    for hit in aligner.map(probe):
        if hit.is_primary:
            return True
    return False


def validate_mmi_smoke(mmi: Path, reference: Path, *, hap_name: str) -> None:
    if not mmi_smoke_aligns(mmi, reference):
        raise InputCheckError(
            f"{hap_name} index failed the smoke alignment (no primary hit for GRCh38 {SMOKE_CHROM}:{SMOKE_START}-{SMOKE_END}).",
            f"Rebuild with: minimap2 -d {mmi} <haplotype.fa>",
        )


def _resolve_mmi(fasta: Path, index: Path | None) -> Path:
    if index is not None:
        return Path(index)
    if fasta.suffix == ".fa":
        return fasta.with_suffix(".mmi")
    return Path(f"{fasta}.mmi")


def validate_haplotype_indexes(
    hap1: Path,
    hap2: Path,
    hap1_index: Path | None,
    hap2_index: Path | None,
    reference: Path,
) -> None:
    for name, fasta, idx in (("hap1", hap1, hap1_index), ("hap2", hap2, hap2_index)):
        mmi_path = _resolve_mmi(Path(fasta), Path(idx) if idx else None)
        validate_mmi_size(mmi_path, Path(fasta), hap_name=name)
        validate_mmi_smoke(mmi_path, reference, hap_name=name)


def no_hit_fraction(labeled_rows: list[dict[str, str]]) -> tuple[float, int, int]:
    """Return (fraction, no_hit_count, eligible_count). Excludes ``no_sequence`` insertions."""
    eligible = [
        r
        for r in labeled_rows
        if str(r.get("label", "NA")) in ("0", "1")
        and str(r.get("reason", "")) != "no_sequence"
    ]
    if not eligible:
        return 0.0, 0, 0
    no_hit = sum(1 for r in eligible if str(r.get("reason")) == "no_hit")
    return no_hit / len(eligible), no_hit, len(eligible)


def check_no_hit_fraction(labeled_rows: list[dict[str, str]], *, max_fraction: float = NO_HIT_MAX_FRACTION) -> None:
    frac, no_hit, total = no_hit_fraction(labeled_rows)
    if total == 0:
        return
    if frac > max_fraction:
        pct = 100 * frac
        limit = 100 * max_fraction
        raise InputCheckError(
            f"Labeling quality check failed: no_hit is {no_hit}/{total} ({pct:.1f}%), above {limit:.0f}% limit. "
            "This usually means corrupt or stale .mmi indexes.",
            "Rebuild both haplotype indexes with minimap2 -d and re-run ml-label.",
        )
