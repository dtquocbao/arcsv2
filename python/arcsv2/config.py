"""Defaults from the ARC-SV2 roadmap §4.12. Paper mode is the default.

Every default is either the published value or the original code value the
roadmap names. Extensions in §8 stay off unless a flag turns them on.
"""

from __future__ import annotations

from pydantic import BaseModel, Field


class CallConfig(BaseModel):
    """Paper-mode caller settings (§4.12)."""

    mode: str = "paper"
    seed: int = 204178949

    # §4.1 read filters
    min_mapq: int = 20
    min_mapq_softclip: int = 10
    min_clipped_bases: int = 5
    min_clipped_qual: int = 15
    lowqual_trim_extra: int = 5
    min_softclip_support: int = 2
    split_read_leeway: int = 2

    # §4.2 library
    insert_sample_size: int = 1_000_000
    max_kde_samples: int = 1_000_000
    # C3: when True, fit the insert PMF on the first insert_sample_size inserts
    # per library (coordinate order). Default False — enable only after golden PASS.
    sample_insert_sizes: bool = False
    robustness_parameter: float = 1e-4
    # Paper: C = exp(3^2 / 2). The log-likelihood-ratio threshold is 4.5.
    discordance_log_c: float = 4.5
    insert_max_mu_multiple: float = 3.0

    # §4.3 breakpoints
    softclip_merge_distance: int = 5
    microhomology_distance: int = 25
    drp_insertion_slop: int = 20
    # orig:arcsv_call_options.max_ins_cluster_slop — post-clique Ins filter.
    drp_insertion_cluster_slop: int = 25
    null_reps: int = 4
    target_fdr: float = 0.1
    cluster_distance_sd: float = 3.0
    max_pecluster_size: int = 80
    bp_confidence: float = 0.95
    # Breakpoint merge always uses clip/split precedence over PE (original
    # breakpoint_merge.py). Paper Methods say "merge all overlapping"; see
    # docs/PARITY.md.
    # Novel skip edges: withhold PE support when the insert under the reference
    # path through middle blocks is still in [q_low, q_high]. Not in literal
    # orig:add_pair. Default False (literal add_pair): original kills the weak
    # ~80 bp pair-only DELs at genotype (ref preferred / score_vs_ref<=0), which
    # ARC-SV2 already has. Set True to restore the ARC-SV2-only ref-path veto.
    skip_edge_ref_gate: bool = False

    # §4.4 graph
    insert_qlow: float = 0.025
    insert_qhigh: float = 0.975
    min_edge_support: int = 2
    max_edge_span: int = 2_000_000
    # orig:arcsv_call_options.py. Pairs farther than this, or on another chromosome, are hanging.
    max_pair_distance: int = 2_000_000

    # §4.5 paths. Paper skips loci with more than 1,000 traversals.
    max_back_edges: int = 2
    max_paths: int = 1000
    # Original default is 0, so a locus over the cap is skipped. Values of 2, 4, …
    # raise min_edge_support by that much, two at a time, and retry.
    max_mes_extra: int = 0
    top_paths: int = 50
    # Loci are genotyped in separate processes when this is greater than 1.
    threads: int = 1

    # §4.7
    allele_fractions: list[float] = Field(default_factory=lambda: [0.5])

    # §4.9
    min_affected_length: int = 50

    # §4.10 labels
    validation_flank: int = 150
    validation_mapq: int = 60
    validation_identity: float = 0.90
    validation_max_gap: int = 48

    # Paper mode pools libraries into one insert model (§8.5).
    pool_libraries: bool = True

    # C5 optional smoove-style noise filters. Default OFF keeps golden identical.
    # When True, drop alignments with NM (edit distance) above max_nm.
    filter_high_nm: bool = False
    max_nm: int = 5
    # When True, drop pairs where both mates are soft-clipped.
    filter_both_mates_clipped: bool = False

    # M2: store pairs as numpy SoA instead of per-pair Python objects.
    # list[PairRecord] remains available when False.
    columnar_pairs: bool = True

    # R0+: python | rust evidence/kernels. Default rust after R1/R2 golden PASS
    # (BAM load + block prep/build_graph). engine=rust requires arcsv2_native.
    # --engine python keeps the pure-Python path.
    engine: str = "rust"


OUTPUT_COLUMNS: tuple[str, ...] = (
    "chrom",
    "minbp",
    "maxbp",
    "id",
    "svtype",
    "complextype",
    "num_sv",
    "bp",
    "bp_uncertainty",
    "reference",
    "rearrangement",
    "len_affected",
    "filter",
    "sv_bp",
    "sv_bp_uncertainty",
    "gt",
    "af",
    "inslen",
    "sr_support",
    "pe_support",
    "score_vs_ref",
    "score_vs_next",
    "rearrangement_next",
    "num_paths",
    # Reference intervals for the locus blocks, in name order. Not reconstructed from bp.
    "blocks",
    # Same blocks, closed at the clip position or the interval midpoint.
    "blocks_point",
    # Per internal junction: clip or mid (same order as ``blocks_point`` gaps).
    "bp_point_source",
)

FEATURE_SCHEMA = "arcsv2-f1"
