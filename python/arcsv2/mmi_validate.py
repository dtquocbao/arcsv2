"""Haplotype ``.mmi`` validation before labeling (U4 / ML-0 HG00140)."""

from __future__ import annotations

from pathlib import Path

from arcsv2.input_checks import InputCheckError

# Healthy HPRC indexes are ~5.5–7× the compressed FASTA size; corrupt HG00140 was ~4.6×.
MMI_BYTES_PER_FASTA_GZ = 5.2
MMI_BYTES_PER_FASTA = 0.85

# Outside the chr20 centromere (GRCh38 chr20:26.44–30.04 Mb) and a second
# euchromatic window on chr1 (centromere ~122–125 Mb). Both are 1 kb.
SMOKE_PROBES = (
    ("chr20", 40_000_000, 40_001_000),
    ("chr1", 40_000_000, 40_001_000),
)

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
    """True only when every probe has a primary hit. Missing tools or an unreadable reference raise."""
    try:
        import mappy
        import pysam
    except ImportError as exc:
        raise InputCheckError(
            f"Smoke alignment cannot run without mappy and pysam ({exc}).",
            "Run arcsv2 mmi-validate inside the ARC-SV2 container.",
        ) from exc
    if not Path(mmi).is_file():
        raise InputCheckError(
            f"Haplotype index missing: {mmi}",
            f"Rebuild with: minimap2 -d {mmi} <haplotype.fa>",
        )
    if not Path(reference).is_file():
        raise InputCheckError(
            f"Reference FASTA missing for the smoke alignment: {reference}",
            "Pass the GRCh38 analysis-set FASTA as --reference.",
        )
    try:
        fasta = pysam.FastaFile(str(reference))
    except (OSError, ValueError) as exc:
        raise InputCheckError(
            f"Could not read reference FASTA {reference}: {exc}",
            "Pass a FASTA pysam can open, with an index if the file is bgzipped.",
        ) from exc
    try:
        references = set(fasta.references)
        missing = [chrom for chrom, _start, _end in SMOKE_PROBES if chrom not in references]
        if missing:
            raise InputCheckError(
                f"Reference {reference} has no {', '.join(missing)} for the smoke probe.",
                "Use the GRCh38 analysis set, which includes chr1 and chr20.",
            )
        sequences: list[tuple[str, int, int, str]] = []
        for chrom, start, end in SMOKE_PROBES:
            try:
                seq = fasta.fetch(chrom, start, end)
            except (OSError, ValueError) as exc:
                raise InputCheckError(
                    f"Could not fetch {chrom}:{start}-{end} from {reference}: {exc}",
                    "Use the GRCh38 analysis set and a readable FASTA index.",
                ) from exc
            if len(seq) < 100:
                raise InputCheckError(
                    f"Smoke probe {chrom}:{start}-{end} is too short ({len(seq)} bp).",
                    "The probes are chr20:40000000-40001000 and chr1:40000000-40001000.",
                )
            sequences.append((chrom, start, end, seq))
    finally:
        fasta.close()

    aligner = mappy.Aligner(str(mmi), n_threads=1)
    if not aligner:
        return False
    for _chrom, _start, _end, seq in sequences:
        if not any(hit.is_primary for hit in aligner.map(seq)):
            return False
    return True


def validate_mmi_smoke(mmi: Path, reference: Path, *, hap_name: str) -> None:
    if not mmi_smoke_aligns(mmi, reference):
        where = ", ".join(f"{chrom}:{start}-{end}" for chrom, start, end in SMOKE_PROBES)
        raise InputCheckError(
            f"{hap_name} index failed the smoke alignment (no primary hit for GRCh38 {where}).",
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
