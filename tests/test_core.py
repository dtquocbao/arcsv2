"""Unit tests for the §4 kernels and a constructed deletion call."""

from __future__ import annotations

import random

from arcsv2.config import OUTPUT_COLUMNS, CallConfig
from arcsv2.engine.breakpoints import (
    Breakpoint,
    Drp,
    Softclip,
    clamp_breakpoint,
    classify_split_flags,
    cluster_drps,
    cluster_softclips,
    fdr_keep,
    pearson_q95,
    score_cluster,
    uniform_ci,
)
from arcsv2.engine.caller import (
    PairRecord,
    _drp_from_pair,
    _format_inslen,
    _insertion_position,
    _stamp_clip_positions,
    build_graph,
    call_evidence,
    insertion_length_is_bound,
    reference_block_lengths_midpoint,
)
from arcsv2.engine.classify import affected_length, complex_type
from arcsv2.engine.graph import (
    Block,
    Edge,
    EdgeRead,
    Graph,
    PathLayout,
    block_distance,
    create_blocks,
    enumerate_paths,
    placement_distances,
)
from arcsv2.engine.library import (
    double_int,
    eval_pmf,
    fit_library,
    pi_robust_one,
    pmf_kernel_smooth,
    sample_chunks,
)
from arcsv2.engine.likelihood import genotype_locus, reference_path, score_path, select_in_band
from arcsv2.engine.names import path_nodes_to_tokens


def test_pmf_sums_to_one_and_pi_is_small():
    inserts = [380, 390, 400, 400, 400, 410, 420] * 50
    pmf = pmf_kernel_smooth(inserts, 0, 1200)
    assert abs(sum(pmf) - 1.0) < 1e-12
    lib = fit_library(inserts)
    # A tight library sits near the median; the deletion cutoff is past it.
    assert lib.z_del_cut > lib.insert_median
    assert lib.z_ins_cut < lib.insert_median
    pi = pi_robust_one(pmf, 1e-4)
    assert 0 < pi <= max(pmf)


def test_distant_hanging_uses_the_distant_class():
    lib = fit_library([400] * 30, read_len_short=100, read_len_long=100)
    lib.class_probs = [0.5, 0.1, 0.05, 0.2, 0.15]
    graph = Graph(create_blocks([], 1_000))
    edge = graph.add_edge(
        0,
        1,
        support=0,
        read=EdgeRead(
            offset=10,
            lib_id=0,
            pmapped=0.5,
            hanging=True,
            mapq=60,
            kind="pe",
            mate_len=100,
        ),
    )
    graph.add_edge(
        0,
        1,
        support=0,
        read=EdgeRead(
            offset=10,
            lib_id=0,
            pmapped=0.5,
            hanging=True,
            mapq=60,
            kind="pe",
            distant=True,
            mate_len=100,
        ),
    )
    scored = score_path((0, 1), graph, [lib], [(edge, 0), (edge, 1)])
    assert abs(scored.read_p[0] - 0.1) < 1e-12
    assert abs(scored.read_p[1] - 0.2) < 1e-12


def test_sample_chunks_is_deterministic_and_stops():
    chunks = [[1, 2], [3], [4, 5, 6]]
    a = sample_chunks(chunks, 4, random.Random(1))
    b = sample_chunks(chunks, 4, random.Random(1))
    assert a == b and len(a) == 4
    try:
        sample_chunks([], 10, random.Random(0))
        raise AssertionError("expected an error")
    except RuntimeError as exc:
        assert "zero usable" in str(exc)


def test_double_int_matches_brute_force():
    pmf = [0.0, 0.2, 0.5, 0.3]
    from arcsv2.engine.library import cdf_and_cdf_sum

    _cdf, cdf_sum = cdf_and_cdf_sum(pmf)

    def brute(a1, b1, a2, b2):
        # Inclusive endpoints: fragment length is y - x + 1.
        # orig:double_int sums f(y - x + 1), not f(y - x).
        total = 0.0
        for x in range(a1, b1 + 1):
            for y in range(a2, b2 + 1):
                total += eval_pmf(pmf, y - x + 1)
        return total

    got = double_int(cdf_sum, 0, 5, 0, 5, xmin=0, xmax=5, ymin=0, ymax=5)
    assert abs(got - brute(0, 5, 0, 5)) < 1e-12


def test_softclip_order_invariant():
    clips = [
        Softclip(10, "R", 30, 8, 0, 10),
        Softclip(12, "R", 40, 8, 2, 12),
        Softclip(40, "L", 20, 6, 40, 50),
        Softclip(11, "R", 10, 9, 1, 11),
    ]
    shuffled = list(reversed(clips))
    a = [(c.side, c.pos, c.n_reads, c.sum_mapq) for c in cluster_softclips(clips, 5)]
    b = [(c.side, c.pos, c.n_reads, c.sum_mapq) for c in cluster_softclips(shuffled, 5)]
    assert a == b
    right = [c for c in cluster_softclips(clips, 5) if c.side == "R"]
    assert len(right) == 1
    # Highest total MAPQ: pos 10 has 30, pos 12 has 40, pos 11 has 10.
    assert right[0].pos == 12


def test_every_split_type():
    assert classify_split_flags(True, True, False, False) == "Del+"
    assert classify_split_flags(False, False, False, False) == "Del-"
    assert classify_split_flags(True, True, True, False) == "Dup+"
    assert classify_split_flags(False, False, True, False) == "Dup-"
    assert classify_split_flags(True, False, False, True) == "InvL+"
    assert classify_split_flags(True, False, False, False) == "InvL-"
    assert classify_split_flags(False, True, False, True) == "InvR+"
    assert classify_split_flags(False, True, False, False) == "InvR-"


def test_clique_and_fdr():
    pairs = [
        Drp("Del", 100, 500, 800, 0),
        Drp("Del", 110, 510, 790, 1),
        Drp("Del", 90, 490, 810, 2),
        Drp("Del", 5000, 5400, 800, 3),
    ]
    shuffled = list(reversed(pairs))
    a, _ = cluster_drps(pairs, max_distance=50)
    b, _ = cluster_drps(shuffled, max_distance=50)
    assert [tuple(p.index for p in cl) for cl in a] == [tuple(p.index for p in cl) for cl in b]
    assert any(len(cl) == 3 for cl in a)
    scores = [score_cluster(cl, 400, 30) for cl in a]
    kept = fdr_keep(scores, null_scores=[0.0, 0.0, 0.0, 0.0], null_reps=4, alpha=0.1)
    assert kept  # the real clusters outscore an all-zero null


def test_uniform_ci_known_quantile():
    assert pearson_q95(2) == 20
    lower, upper = uniform_ci([0, 10])
    assert lower[0] == int(__import__("math").floor((1 - 20) * 10))
    assert upper[1] == int(__import__("math").ceil(20 * 10))


def test_affected_length_matches_published_examples():
    # Lengths 100, 50, 150, 25, 100. orig:sv_affected_len.py::test_affected_len.
    lengths = {0: 100, 1: 50, 2: 150, 3: 25, 4: 100}
    ref = [(i, 1) for i in range(5)]
    assert affected_length(ref, ref, lengths)[0] == 0
    inverted = [(0, 1), (1, 1), (2, 1), (3, -1), (4, 1)]
    assert affected_length(ref, inverted, lengths)[0] == 25
    deletion = [(0, 1), (4, 1)]
    assert affected_length(ref, deletion, lengths)[0] == 225
    dup = [(0, 1), (1, 1), (2, 1), (3, 1), (2, 1), (3, 1), (4, 1)]
    assert affected_length(ref, dup, lengths)[0] == 175


def test_reference_block_lengths_include_median_gaps():
    """orig:sv_output blocks_midpoints — middle block 47 + gaps 7,7 → 54."""
    # Synthetic geometry mirroring HG01884 DEL 613786 (v2fix7 blocks).
    blocks = [
        Block(610416, 613783),
        Block(613790, 613837),  # uncovered length 47
        Block(613844, 614491),
        Block(614492, 628873),
    ]
    lengths = reference_block_lengths_midpoint(blocks, 0, 3)
    assert lengths[1] == 54
    assert len(blocks[1]) == 47
    ref = [(i, 1) for i in range(4)]
    alt = [(0, 1), (2, 1), (3, 1)]  # ACD
    assert affected_length(ref, alt, lengths)[0] == 54

    # Six HG01884 v2fix7 geometries from CURSOR_TASK_RECALL (affected block = 1).
    cases = [
        # name, block intervals, expected midpoint length of block B
        (
            [(610416, 613783), (613790, 613837), (613844, 614491), (614492, 628873)],
            54,
        ),
        (
            [(6598812, 6606006), (6606035, 6606062), (6606088, 6606325), (6606352, 6607258)],
            55,  # 27 + floor medians; orig had different BP topology → 54
        ),
        (
            [(9165835, 9168922), (9168927, 9168975), (9168980, 9172618)],
            53,
        ),
        (
            [(11004830, 11021356), (11021387, 11021422), (11021453, 11028694)],
            66,
        ),
        (
            [(14071279, 14072398), (14072438, 14072461), (14072501, 14080442)],
            63,
        ),
        (
            [(30599883, 30602113), (30602127, 30602176), (30602190, 30621650)],
            63,
        ),
    ]
    for intervals, expect in cases:
        blks = [Block(s, e) for s, e in intervals]
        mid = reference_block_lengths_midpoint(blks, 0, len(blks) - 1)
        assert mid[1] == expect, (intervals, mid[1], expect)
        assert mid[1] >= 50


def test_complex_type_table():
    assert complex_type((0, 1, 2, 3, 4, 5, 2, 3, 6, 7)) == "ddup"
    assert complex_type((0, 1, 2, 3, 4, 5, 2, 3, 8, 9)) == "ddup.del"
    assert complex_type((0, 1, 2, 3, 5, 4, 3, 2, 6, 7)) == "invddup.inv"


def test_simplify_blocks_diploid_collapses_unused_cut_to_abbc():
    """ABCD / ABCBCD with only reference B→C adjacencies → ABC / ABBC.

    Matches orig:sv_validate.simplify_blocks_diploid at HG01884 ~60792115.
    """
    from arcsv2.engine.classify import classify_simple, simplify_blocks_diploid

    blocks = [
        Block(0, 100),
        Block(110, 140),
        Block(140, 170),
        Block(180, 300),
    ]
    # Hap1 ABCBCD, hap2 reference ABCD — middle cut unused as a novel junction.
    path_alt = (0, 1, 2, 3, 4, 5, 2, 3, 4, 5, 6, 7)
    path_ref = (0, 1, 2, 3, 4, 5, 6, 7)
    nb, np1, np2, groups = simplify_blocks_diploid(blocks, path_alt, path_ref)
    assert len([b for b in nb if not b.is_insertion]) == 3
    assert groups == [[0], [1, 2], [3]]
    assert (nb[1].start, nb[1].end) == (110, 170)
    # Remapped: A B B C
    assert path_nodes_to_tokens(np1) == [(0, 1), (1, 1), (1, 1), (2, 1)]
    assert path_nodes_to_tokens(np2) == [(0, 1), (1, 1), (2, 1)]
    named_ref = [(0, 1), (1, 1), (2, 1)]
    named_alt = [(0, 1), (1, 1), (1, 1), (2, 1)]
    assert classify_simple(named_ref, named_alt) == "DUP"


def test_simplify_blocks_diploid_preserves_allele_base_coords():
    """ABCBCD → ABBC must keep the labeler's predicted allele bases identical.

    Gapped B|C cut (as at ~60792115): fusion renames paths but closed
    blocks_point spans along alt/ref must not change.
    """
    from arcsv2.engine.caller import _allele_base_coords
    from arcsv2.engine.classify import simplify_blocks_diploid

    blocks = [
        Block(0, 100),
        Block(110, 140),  # B
        Block(150, 170),  # C — gap 140-150 (extra cut)
        Block(180, 300),
    ]
    path_alt = (0, 1, 2, 3, 4, 5, 2, 3, 4, 5, 6, 7)  # ABCBCD
    path_ref = (0, 1, 2, 3, 4, 5, 6, 7)  # ABCD
    before_alt = _allele_base_coords(blocks, path_alt, 0, 3)
    before_ref = _allele_base_coords(blocks, path_ref, 0, 3)
    nb, np1, np2, groups = simplify_blocks_diploid(blocks, path_alt, path_ref)
    assert groups == [[0], [1, 2], [3]]
    after_alt = _allele_base_coords(nb, np1, 0, 2)
    after_ref = _allele_base_coords(nb, np2, 0, 2)
    assert after_alt == before_alt
    assert after_ref == before_ref
    # Sanity: duplicated middle is B∪C closed span, not a truncated half.
    assert len(before_alt) > len(before_ref)
    # Outer A/D flanks are unchanged, so 150 bp label queries match too.
    assert _allele_base_coords(blocks, path_alt, 0, 3, flank=150) == _allele_base_coords(
        nb, np1, 0, 2, flank=150
    )


def test_simplify_del_acd_changes_flank_trimmed_allele():
    """DEL ACD→AC keeps full bases but moves the right 150 bp flank cut.

    Fusion of kept C into outer D is why five window TPs lost label=1 after
    seq-preserving simplify; caller must detect this and skip simplify.
    """
    from arcsv2.engine.caller import _allele_base_coords
    from arcsv2.engine.classify import simplify_blocks_diploid

    blocks = [
        Block(0, 200),
        Block(210, 260),  # deleted B
        Block(270, 500),  # C (kept)
        Block(510, 800),  # D (outer)
    ]
    path_del = (0, 1, 4, 5, 6, 7)  # ACD
    path_ref = (0, 1, 2, 3, 4, 5, 6, 7)  # ABCD
    nb, np1, np2, groups = simplify_blocks_diploid(blocks, path_del, path_ref)
    assert groups == [[0], [1], [2, 3]]
    assert _allele_base_coords(blocks, path_del, 0, 3) == _allele_base_coords(
        nb, np1, 0, 2
    )
    assert _allele_base_coords(blocks, path_del, 0, 3, flank=150) != _allele_base_coords(
        nb, np1, 0, 2, flank=150
    )


def test_simplify_blocks_diploid_keeps_novel_cut():
    """A novel skip across the middle cut must not fuse the skipped block away."""
    from arcsv2.engine.classify import classify_simple, simplify_blocks_diploid

    blocks = [Block(0, 100), Block(110, 140), Block(140, 170), Block(180, 300)]
    # DEL of B: A → C → D. Ref still uses B; C→D is reference-only so C+D fuse.
    path_del = (0, 1, 4, 5, 6, 7)
    path_ref = (0, 1, 2, 3, 4, 5, 6, 7)
    nb, np1, np2, groups = simplify_blocks_diploid(blocks, path_del, path_ref)
    assert groups == [[0], [1], [2, 3]]
    assert path_nodes_to_tokens(np2) == [(0, 1), (1, 1), (2, 1)]
    assert path_nodes_to_tokens(np1) == [(0, 1), (2, 1)]
    assert classify_simple([(0, 1), (1, 1), (2, 1)], [(0, 1), (2, 1)]) == "DEL"


def test_breakpoint_interval_is_left_uncovered():
    blocks = create_blocks([(100, 110), (200, 205)], 300)
    assert [(block.start, block.end) for block in blocks] == [(0, 100), (110, 200), (205, 300)]


def test_find_block_binary_search():
    from arcsv2.engine.graph import find_block

    blocks = create_blocks([(100, 110), (200, 205)], 300)
    assert find_block(blocks, 0) == 0
    assert find_block(blocks, 99) == 0
    assert find_block(blocks, 100) == 1
    assert find_block(blocks, 150) == 1
    assert find_block(blocks, 200) == 2
    assert find_block(blocks, 299) == 2
    assert find_block(blocks, 300) is None
    # Trailing insertion blocks are ignored by the search prefix.
    blocks.append(Block(0, 0, is_insertion=True, insertion_len=50))
    assert find_block(blocks, 150) == 1
    assert find_block(blocks, 299) == 2


def test_overlapping_vertices_matches_linear_scan():
    from arcsv2.engine.graph import overlapping_vertices

    blocks = create_blocks([(100 * i, 100 * i + 10) for i in range(1, 200)], 25_000)
    blocks.append(Block(0, 0, is_insertion=True, insertion_len=40))

    def linear(start: int, end: int, reverse: bool) -> list[int]:
        hits = [
            i
            for i, block in enumerate(blocks)
            if not block.is_insertion and end > block.start and start < block.end
        ]
        if reverse:
            return [2 * block for block in reversed(hits)]
        return [2 * block + 1 for block in hits]

    cases = [
        (0, 50, False),
        (150, 160, False),
        (95, 215, True),
        (5000, 5010, False),
        (19_950, 20_100, True),
        (24_900, 25_000, False),
        (25_000, 25_100, False),
        (110, 110, False),
    ]
    for start, end, reverse in cases:
        assert overlapping_vertices(blocks, start, end, reverse) == linear(start, end, reverse)


def test_breakpoint_interval_is_clamped_to_the_chromosome():
    kept = clamp_breakpoint(Breakpoint(-1471, -1400, n_drp=4, drp_types=("Del",)), 20_000)
    assert kept is None
    clipped = clamp_breakpoint(Breakpoint(-10, 25, n_drp=2, drp_types=("Del",)), 20_000)
    assert clipped is not None
    assert (clipped.start, clipped.end) == (0, 25)


def test_blocks_point_meets_at_the_clip():
    from arcsv2.engine.breakpoints import Softclip
    from arcsv2.engine.caller import (
        _blocks_column,
        _blocks_point_column,
        _bp_point_source_column,
        _stamp_clip_positions_from_softclips,
    )

    blocks = [Block(0, 33_000), Block(36_000, 50_000), Block(80_000, 100_000)]
    blocks[0].clip_pos = 33_333
    blocks[1].clip_pos = 55_000
    assert _blocks_column(blocks, 0, 2) == "0-33000,36000-50000,80000-100000"
    assert _blocks_point_column(blocks, 0, 2) == "0-33333,33333-55000,55000-100000"
    assert _bp_point_source_column(blocks, 0, 2) == "clip,clip"
    gap_blocks = [Block(536_580, 536_794), Block(536_885, 537_000)]
    _stamp_clip_positions_from_softclips(
        gap_blocks,
        [Softclip(pos=536_880, side="R", mapq=60, clipped_len=10, aligned_start=536_700, aligned_end=536_880)],
    )
    assert gap_blocks[0].clip_pos == 536_880
    inserted = [Block(0, 1000), Block(0, 0, is_insertion=True), Block(1400, 2000)]
    inserted[0].clip_pos = 1200
    assert _blocks_point_column(inserted, 0, 2) == "0-1200,_,1200-2000"





def test_merge_precedence_keeps_split_ends_for_overlapping_drp_cis():
    """chr20:2883362 geometry: overlapping Del CIs + split ends → 2 BPs under precedence."""
    from arcsv2.engine.breakpoints import Breakpoint, merge_breakpoints

    left_split = Breakpoint(2_883_362, 2_883_376, n_split=20)
    right_split = Breakpoint(2_883_670, 2_883_684, n_split=20)
    left_drp = Breakpoint(2_883_365, 2_883_569, n_drp=9, drp_types=("Del",), cluster_ids=(0,))
    right_drp = Breakpoint(2_883_496, 2_883_682, n_drp=9, drp_types=("Del",), cluster_ids=(0,))
    pieces = [left_split, right_split, left_drp, right_drp]

    prec = merge_breakpoints(pieces, min_split_or_clip=2)
    assert len(prec) == 2
    assert (prec[0].start, prec[0].end) == (2_883_362, 2_883_376)
    assert (prec[1].start, prec[1].end) == (2_883_670, 2_883_684)
    assert prec[0].n_drp > 0 and prec[1].n_drp > 0
    assert prec[0].n_split >= 20 and prec[1].n_split >= 20


def test_score_cluster_del_uses_mean_insert():
    """Del size D matches pecluster.lr_del / Methods eq. 2 (mean z, not min z)."""
    from arcsv2.engine.breakpoints import Drp, score_cluster

    cluster = [
        Drp("Del", 100, 600, 500, 0),
        Drp("Del", 110, 610, 700, 1),
        Drp("Del", 90, 590, 900, 2),
    ]
    median, scale = 400.0, 30.0
    mean_z = (500 + 700 + 900) / 3
    max_size = 590 - 110
    d_mean = min(mean_z - median, float(max_size))
    d_min = min(500 - median, float(max_size))
    assert d_mean != d_min
    var = 2.0 * scale * scale
    expected = sum(((z - median) ** 2 - (z - median - d_mean) ** 2) / var for z in (500, 700, 900))
    assert abs(score_cluster(cluster, median, scale) - expected) < 1e-9


def test_reference_path_on_three_blocks():
    graph = Graph(create_blocks([(100, 101), (200, 201)], 300))
    graph.add_reference_edges(2)
    paths, status = enumerate_paths(graph, max_paths=10, min_support=2, end=2 * (graph.n_blocks - 1) + 1)
    assert status == "ok"
    assert reference_path(0, graph.n_blocks - 1) in paths


def test_deletion_is_called():
    inserts = [390, 395, 400, 400, 405, 410] * 40
    lib = fit_library(inserts, read_len_short=200, read_len_long=200)
    pairs = []
    # Concordant background, away from the deletion.
    for i in range(30):
        start = 200 + i * 30
        pairs.append(
            PairRecord(
                start1=start,
                end1=start + 200,
                start2=start + 200,
                end2=start + 400,
                reverse1=False,
                reverse2=True,
                mapq1=60,
                mapq2=60,
                insert=400,
            )
        )
    # Six deletion-spanning pairs. Implied insert lands on the library median.
    for i in range(6):
        pairs.append(
            PairRecord(
                start1=4800,
                end1=5000,
                start2=5500,
                end2=5700,
                reverse1=False,
                reverse2=True,
                mapq1=60,
                mapq2=60,
                insert=900,
                qname_hash=i,
            )
        )
    from arcsv2.engine.breakpoints import Breakpoint

    bps = [Breakpoint(5000, 5001, n_drp=6, drp_types=("Del",)), Breakpoint(5500, 5501, n_drp=6, drp_types=("Del",))]
    rows, _logs = call_evidence(
        pairs,
        chrom="chr20",
        chrom_len=20_000,
        cfg=CallConfig(),
        libs=[lib],
        breakpoints=bps,
    )
    assert rows, _logs
    row = rows[0].as_dict()
    assert list(row) == list(OUTPUT_COLUMNS)
    assert row["svtype"] == "DEL"
    assert row["reference"] == "ABC"
    assert row["rearrangement"] == "AC"
    from arcsv2.label import blocks_from_row

    assert blocks_from_row(row) == {
        "A": (0, 5000),
        "B": (5001, 5500),
        "C": (5501, 20_000),
    }
    # No clip on these intervals, so the point estimate is the midpoint and the
    # 1 bp gaps close. ``blocks`` still records the uncovered intervals.
    assert row["blocks_point"] == "0-5000,5000-5500,5500-20000"
    assert row["bp_point_source"] == "mid,mid"
    # Uncovered blocks column still records 5001-5500 (499 bp), but len_affected
    # uses midpoint-extended block lengths (Methods / orig:sv_output).
    assert int(row["len_affected"]) == 500
    assert row["gt"] == "HOM"
    assert float(row["score_vs_ref"]) > 0
    # Breakpoint intervals (5000, 5001) and (5500, 5501) are the internal gaps.
    assert row["bp_uncertainty"] == "0,1,1,0"
    assert row["sv_bp_uncertainty"] == "1,1"
    # One path is over a cap of 1. Raising support by 2 drops the reference
    # edges and the deletion path remains.
    skipped, skip_logs = call_evidence(
        pairs,
        chrom="chr20",
        chrom_len=20_000,
        cfg=CallConfig(max_paths=1, max_mes_extra=0),
        libs=[lib],
        breakpoints=bps,
    )
    assert skipped == []
    assert any("skipped_npaths" in line for line in skip_logs)
    recovered, _re_logs = call_evidence(
        pairs,
        chrom="chr20",
        chrom_len=20_000,
        cfg=CallConfig(max_paths=1, max_mes_extra=2),
        libs=[lib],
        breakpoints=bps,
    )
    assert recovered and recovered[0].as_dict()["svtype"] == "DEL"


def test_events_that_only_share_a_cut_point_stay_separate():
    """Three deletions 50 kb apart. Touching loci must not become one graph."""
    inserts = [390, 395, 400, 400, 405, 410] * 40
    lib = fit_library(inserts, read_len_short=200, read_len_long=200)
    events = [(50_000, 50_500), (100_000, 100_500), (150_000, 150_500)]
    pairs = []
    bps = []
    q = 0
    for left, right in events:
        for i in range(30):
            start = left - 3_000 + i * 30
            pairs.append(
                PairRecord(
                    start1=start,
                    end1=start + 200,
                    start2=start + 200,
                    end2=start + 400,
                    reverse1=False,
                    reverse2=True,
                    mapq1=60,
                    mapq2=60,
                    insert=400,
                    qname_hash=q,
                )
            )
            q += 1
        for i in range(6):
            pairs.append(
                PairRecord(
                    start1=left - 200,
                    end1=left,
                    start2=right,
                    end2=right + 200,
                    reverse1=False,
                    reverse2=True,
                    mapq1=60,
                    mapq2=60,
                    insert=900,
                    qname_hash=q,
                )
            )
            q += 1
        bps.append(Breakpoint(left, left + 1, n_drp=6, drp_types=("Del",)))
        bps.append(Breakpoint(right, right + 1, n_drp=6, drp_types=("Del",)))
    rows, logs = call_evidence(
        pairs,
        chrom="chrSim",
        chrom_len=200_000,
        cfg=CallConfig(),
        libs=[lib],
        breakpoints=bps,
    )
    assert len(rows) == 3, logs
    assert [row.as_dict()["svtype"] for row in rows] == ["DEL", "DEL", "DEL"]
    assert [row.as_dict()["gt"] for row in rows] == ["HOM", "HOM", "HOM"]
    parallel, parallel_logs = call_evidence(
        pairs,
        chrom="chrSim",
        chrom_len=200_000,
        cfg=CallConfig(threads=2),
        libs=[lib],
        breakpoints=bps,
    )
    assert [row.as_dict() for row in parallel] == [row.as_dict() for row in rows], parallel_logs


def _spike_library():
    return fit_library([400] * 400, read_len_short=200, read_len_long=200)


def _copies(pair: PairRecord, n: int) -> list[PairRecord]:
    return [
        PairRecord(
            pair.start1,
            pair.end1,
            pair.start2,
            pair.end2,
            pair.reverse1,
            pair.reverse2,
            pair.mapq1,
            pair.mapq2,
            pair.insert,
            qname_hash=i,
        )
        for i in range(n)
    ]


def test_inversion_drp_uses_the_original_ends():
    lib = _spike_library()
    both_rev = PairRecord(10000, 10200, 11000, 11200, True, True, 60, 60, 400)
    both_fwd = PairRecord(9800, 10000, 10800, 11000, False, False, 60, 60, 400)
    rev = _drp_from_pair(both_rev, lib, 0)
    fwd = _drp_from_pair(both_fwd, lib, 1)
    assert rev is not None and rev.type == "InvR"
    assert (rev.a, rev.b) == (10000, 11000)
    assert fwd is not None and fwd.type == "InvL"
    assert (fwd.a, fwd.b) == (10000, 11000)


def test_hanging_pair_is_not_an_inversion_drp():
    """orig:process_discordant_pair never sees hangings — only completed pairs."""
    lib = _spike_library()
    hang = PairRecord(10000, 10200, 0, 0, False, False, 60, 0, 0, hanging=True)
    assert _drp_from_pair(hang, lib, 0) is None
    # Same-orient completed pair still classifies.
    assert _drp_from_pair(
        PairRecord(9800, 10000, 10800, 11000, False, False, 60, 60, 400), lib, 1
    ).type == "InvL"


def test_dup_geometry_includes_nested_fr_pairs():
    """orig Dup: (fwd_start > rev_start) OR (fwd_end > rev_end); checked before Ins."""
    lib = _spike_library()
    # Nested FR with short signed insert (below z_ins) → Dup, not Ins.
    nested = PairRecord(1000, 1200, 1050, 1150, False, True, 60, 60, 150)
    assert nested.insert < lib.z_ins_cut
    drp = _drp_from_pair(nested, lib, 0)
    assert drp is not None and drp.type == "Dup"
    assert (drp.a, drp.b) == (1050, 1200)
    # Outward: signed insert is negative (rev upstream of fwd) → Dup.
    outward = PairRecord(2000, 2200, 1000, 1200, False, True, 60, 60, -800)
    out = _drp_from_pair(outward, lib, 1)
    assert out is not None and out.type == "Dup"


def test_concordant_opposite_insert_is_not_a_drp():
    """orig:pecluster early return — opposite + insert in [z_ins, z_del] → None."""
    lib = _spike_library()
    # Nested FR geometry but concordant insert → not Dup/Ins.
    concordant_nested = PairRecord(1000, 1200, 1050, 1150, False, True, 60, 60, 400)
    assert lib.z_ins_cut <= 400 <= lib.z_del_cut
    assert _drp_from_pair(concordant_nested, lib, 0) is None
    # Overlapping short-looking mates with concordant insert → not Ins.
    concordant_overlap = PairRecord(1000, 1180, 1100, 1280, False, True, 60, 60, 400)
    assert _drp_from_pair(concordant_overlap, lib, 1) is None


def test_ins_allows_overlapping_mates_and_filters_cluster_slop():
    """Ins needs no non-overlap; clusters with overlap > cluster_slop are dropped."""
    from arcsv2.engine.breakpoints import Drp, cluster_drps, ins_cluster_ok

    lib = _spike_library()
    # Short insert, overlapping mates (fwd_end > rev_start) — orig Ins, not dropped.
    overlap_pair = PairRecord(1000, 1180, 1100, 1280, False, True, 60, 60, 50)
    assert overlap_pair.insert < lib.z_ins_cut
    ins = _drp_from_pair(overlap_pair, lib, 0)
    assert ins is not None and ins.type == "Ins"
    assert ins.a == 1180 and ins.b == 1100  # pos1 > pos2 allowed

    # Pairwise overlap 35 > max_ins_pair_slop (20) → no clique (even if cluster_slop is large).
    too_far = [
        Drp("Ins", 100, 80, 50, 0),
        Drp("Ins", 105, 70, 50, 1),
    ]
    assert not ins_cluster_ok(too_far, cluster_slop=25)
    assert cluster_drps(too_far, max_distance=50, ins_slop=20, ins_cluster_slop=25)[0] == []
    assert cluster_drps(too_far, max_distance=50, ins_slop=20, ins_cluster_slop=100)[0] == []

    # Pairwise overlap 17 ≤ 20 → clique; cluster overlap also ≤ cluster_slop 25.
    ok_pairs = [
        Drp("Ins", 100, 90, 50, 0),
        Drp("Ins", 105, 88, 50, 1),
    ]
    assert ins_cluster_ok(ok_pairs, cluster_slop=25)
    kept, _ = cluster_drps(ok_pairs, max_distance=50, ins_slop=20, ins_cluster_slop=25)
    assert len(kept) == 1 and len(kept[0]) == 2


def test_concordant_and_out_of_range_pairs_are_stored_without_support():
    lib = _spike_library()
    blocks = [Block(0, 5000), Block(5000, 6000), Block(6000, 12_000)]
    inside = PairRecord(5200, 5400, 5400, 5600, False, True, 60, 60, 400)
    far = PairRecord(1000, 1200, 7000, 7200, False, True, 60, 60, 6000)
    graph = build_graph(blocks, _copies(inside, 3) + [far], lib, CallConfig())
    within = [edge for edge in graph.edges.values() if edge.u // 2 == edge.v // 2 == 1 and edge.reads]
    assert len(within) == 1 and within[0].support < 2 and len(within[0].reads) == 3
    distant = [edge for edge in graph.edges.values() if edge.reads and {edge.u // 2, edge.v // 2} == {0, 2}]
    assert len(distant) == 1 and distant[0].support < 2


def test_concordant_tail_is_not_a_duplication():
    lib = _spike_library()
    blocks = create_blocks([(10000, 10001), (11000, 11001)], 30_000)
    tail = PairRecord(10100, 10300, 10300, 10500, False, True, 60, 60, 200)
    # Outward geometry needs a non-concordant signed insert to keep Dup support.
    outward = PairRecord(10000, 10200, 10800, 11000, True, False, 60, 60, -600)
    graph = build_graph(blocks, _copies(tail, 4) + _copies(outward, 2), lib, CallConfig())
    dup_edges = [e for e in graph.edges.values() if e.u // 2 == e.v // 2 and not e.reference]
    assert len(dup_edges) == 1
    assert dup_edges[0].support == 2


def test_non_outward_pair_does_not_support_same_block_dup():
    """A concordant-orientation pair must not open a tandem edge.

    On a 200 bp block the no-dup offset is 180 (below q_low) and the
    duplication-path insert is 380 (inside the quantile). Literal add_pair
    would count it. The forward read is upstream of the reverse read, so
    the outward test leaves support at 0.
    """
    from arcsv2.engine.caller import _adjacency_insert, _dup_without_insert, _is_outward

    lib = fit_library(list(range(360, 441)) * 10, read_len_short=100, read_len_long=100)
    blocks = [Block(0, 2000), Block(2000, 2200), Block(2200, 5000)]
    pair = PairRecord(2010, 2110, 2090, 2190, False, True, 60, 60, 180)
    implied = _adjacency_insert(pair, blocks, 1, 1)
    alt = _dup_without_insert(pair, blocks, 1)
    assert not _is_outward(pair)
    assert lib.q_low <= implied <= lib.q_high
    assert alt < lib.q_low
    graph = build_graph(blocks, _copies(pair, 6), lib, CallConfig())
    dup = [e for e in graph.edges.values() if e.u // 2 == e.v // 2 == 1 and e.reads]
    assert len(dup) == 1
    assert dup[0].support == 0



def test_multiblock_pe_between_junction_credits_same_block_dup():
    """orig:add_pair flip_parity junction: multi-block mates credit B.out↔B.in.

    Mate1 spans C←B (rev), mate2 spans A→B (fwd). Stored PE is a skip;
    between support lands on the tandem edge when alt_dist < q_low.
    """
    from arcsv2.engine.caller import _between_junction, _between_junction_supported
    from arcsv2.engine.graph import overlapping_vertices

    lib = fit_library(list(range(360, 441)) * 10, read_len_short=100, read_len_long=100)
    blocks = [Block(0, 1000), Block(1000, 1050), Block(1050, 5000)]
    # Reverse covers B+C → [C.in, B.in]; forward covers A+B → [A.out, B.out].
    pair = PairRecord(1030, 1150, 800, 1020, True, False, 60, 60, 300)
    v1 = overlapping_vertices(blocks, pair.start1, pair.end1, pair.reverse1)
    v2 = overlapping_vertices(blocks, pair.start2, pair.end2, pair.reverse2)
    assert v1 == [4, 2] and v2 == [1, 3]
    assert _between_junction(v1, v2) == (2, 3)
    from arcsv2.engine.caller import _read_offset

    off1 = _read_offset(pair.start1, pair.end1, True, blocks[v1[0] // 2])
    off2 = _read_offset(pair.start2, pair.end2, False, blocks[v2[0] // 2])
    assert _between_junction_supported(pair, v1, v2, off1, off2, blocks, lib)
    graph = build_graph(blocks, _copies(pair, 2), lib, CallConfig(min_edge_support=2))
    tandem = [e for e in graph.edges.values() if e.key() == (2, 3)]
    assert tandem
    assert tandem[0].support >= 2
    assert any(r.kind == "pe" and r.supports_adjacency for r in tandem[0].reads)


def test_gapped_ci_blocks_same_block_support_from_split_and_pe():
    """NA19700 DUP 52764531 geometry with gapped half-open CI blocks.

    Half-open CI intervals leave 1-bp gaps around the tandem block. A Dup+
    split plus a multi-block PE (flip_parity between junction) must still
    reach same-block support >= 2 so decompose opens the locus.
    """
    from arcsv2.engine.breakpoints import Breakpoint
    from arcsv2.engine.caller import _split_path
    from arcsv2.engine.evidence import SplitHit
    from arcsv2.engine.graph import decompose, overlapping_vertices

    # Softclip CIs of width 1 at the two DUP ends (gapped create_blocks).
    cis = [
        Breakpoint(1000, 1001, n_softclip=2, clip_pos=1000, n_split=1),
        Breakpoint(1050, 1051, n_softclip=2, clip_pos=1050, n_split=1),
    ]
    gapped = create_blocks([(bp.start, bp.end) for bp in cis], 5000)
    assert [(b.start, b.end) for b in gapped] == [(0, 1000), (1001, 1050), (1051, 5000)]
    assert gapped[1].start - gapped[0].end == 1

    lib = fit_library(list(range(360, 441)) * 10, read_len_short=100, read_len_long=100)
    # Dup+ halves: first ends at tandem block end; second starts in the 1-bp gap.
    hit = SplitHit(
        start1=1002,
        end1=1050,
        reverse1=False,
        start2=1000,
        end2=1100,
        reverse2=False,
        mapq=60,
        split_type="Dup+",
        bp1=(1050, 1050),
        bp2=(1000, 1000),
        chrom="chrSim",
        primary_start=1002,
        primary_reverse=False,
        qname_hash=42,
    )
    parsed = _split_path(hit, gapped, CallConfig().split_read_leeway)
    assert parsed is not None
    verts, _off = parsed
    assert verts[0] // 2 == 1 and verts[1] // 2 == 1

    # Multi-block PE: rev covers B+C, fwd covers A+B → between (B.in, B.out).
    pair = PairRecord(1030, 1150, 800, 1020, True, False, 60, 60, 300)
    v1 = overlapping_vertices(gapped, pair.start1, pair.end1, pair.reverse1)
    v2 = overlapping_vertices(gapped, pair.start2, pair.end2, pair.reverse2)
    assert v1[0] // 2 == 2 and v1[-1] // 2 == 1
    assert v2[0] // 2 == 0 and v2[-1] // 2 == 1

    cfg = CallConfig(min_edge_support=2)
    graph = build_graph(gapped, _copies(pair, 1), lib, cfg, splits=[hit])
    tandem = [e for e in graph.edges.values() if e.u // 2 == e.v // 2 == 1 and abs(e.u - e.v) == 1]
    assert tandem
    assert tandem[0].support >= 2
    kinds = {r.kind for r in tandem[0].reads if r.supports_adjacency or r.kind == "split"}
    assert "split" in kinds and "pe" in {r.kind for r in tandem[0].reads if r.supports_adjacency}
    loci = decompose(graph, cfg.min_edge_support)
    assert any(lo <= 1 <= hi for lo, hi in loci)


def test_concordant_pair_adds_no_support_to_short_skip_or_same_block_dup():
    """With skip_edge_ref_gate on: ~80 bp middle + concordant FR → no skip support."""
    lib = fit_library(list(range(360, 441)) * 10, read_len_short=200, read_len_long=200)
    # q ≈ [357, 443]. With an 80 bp middle, both inserts are in-window only when
    # skip implied ∈ [357, 363] (then ref = skip+80 ∈ [437, 443]).
    blocks = [Block(0, 1000), Block(1000, 1080), Block(1080, 5000)]
    # offsets sum -640 → skip 360, ref 440 (both in window) → no skip support.
    concordant = PairRecord(800, 1000, 1080, 1240, False, True, 60, 60, 440)
    graph = build_graph(
        blocks,
        _copies(concordant, 4),
        lib,
        CallConfig(min_edge_support=2, skip_edge_ref_gate=True),
    )
    skip = [e for e in graph.edges.values() if e.key() == (1, 4)]  # A.out → C.in
    assert len(skip) == 1
    assert skip[0].support == 0
    assert all(not r.supports_adjacency for r in skip[0].reads if r.kind == "pe")

    # Same-block concordant FR on B: stored but no tandem-dup support.
    one = [Block(0, 5000)]
    inside = PairRecord(2000, 2200, 2200, 2400, False, True, 60, 60, 400)
    g2 = build_graph(one, _copies(inside, 4), lib, CallConfig())
    within = [e for e in g2.edges.values() if e.u // 2 == e.v // 2 == 0 and e.reads]
    assert within
    assert all(e.support == 0 for e in within)


def test_literal_add_pair_credits_concordant_short_skip():
    """Default gate off matches orig:add_pair: short skip gets PE support."""
    lib = fit_library(list(range(360, 441)) * 10, read_len_short=200, read_len_long=200)
    blocks = [Block(0, 1000), Block(1000, 1080), Block(1080, 5000)]
    concordant = PairRecord(800, 1000, 1080, 1240, False, True, 60, 60, 440)
    graph = build_graph(
        blocks,
        _copies(concordant, 4),
        lib,
        CallConfig(min_edge_support=2, skip_edge_ref_gate=False),
    )
    skip = [e for e in graph.edges.values() if e.key() == (1, 4)]
    assert len(skip) == 1
    assert skip[0].support >= 2


def test_multiblock_pe_mate_adds_within_read_novel_support():
    """orig:add_pair within-read walk: multi-block PE (no split) credits path edges.

    An insertion block is skipped by overlapping_vertices, so one mate that
    spans A|ins|C yields verts [A.out, C.out] and within-read support on the
    novel skip A.out→C.in — separate from the between-mate edge to D.
    """
    lib = fit_library(list(range(360, 441)) * 10, read_len_short=200, read_len_long=200)
    blocks = [
        Block(0, 1000),
        Block(1000, 1000, is_insertion=True, insertion_len=80),
        Block(1000, 3000),
        Block(3000, 5000),
    ]
    # Mate1 overlaps A and C (insertion skipped). Mate2 on D.
    # Between-mate anchor is A.out↔D.in; within-read novel is A.out→C.in.
    pair = PairRecord(800, 1200, 3200, 3400, False, True, 60, 60, 2200)
    graph = build_graph(blocks, _copies(pair, 2), lib, CallConfig(min_edge_support=2), splits=[])
    skip = [e for e in graph.edges.values() if e.key() == (1, 4)]  # A.out → C.in
    assert len(skip) == 1
    assert skip[0].support >= 2
    assert any(
        r.kind == "pe" and r.supports_adjacency and not r.in_likelihood for r in skip[0].reads
    )


def test_deletion_pair_still_supports_the_skip_edge():
    """True Del DRP: ref-path insert outside 95% → skip edge keeps PE support."""
    lib = fit_library(list(range(360, 441)) * 10, read_len_short=200, read_len_long=200)
    # 500 bp middle so ref = skip+500 is far outside q while skip stays in-window.
    blocks = [Block(0, 1000), Block(1000, 1500), Block(1500, 5000)]
    # offsets sum -600 → skip 400 (in), ref 900 (out) → keep PE support on skip.
    deletion = PairRecord(800, 1000, 1500, 1700, False, True, 60, 60, 900)
    graph = build_graph(blocks, _copies(deletion, 3), lib, CallConfig(min_edge_support=2))
    skip = [e for e in graph.edges.values() if e.key() == (1, 4)]
    assert len(skip) == 1
    assert skip[0].support >= 2
    assert sum(1 for r in skip[0].reads if r.supports_adjacency) >= 2


def test_pair_spanning_a_short_block_keeps_reference_probability():
    """A pair that jumps a short block is still possible on the reference path.

    The reference insert is the raw offset plus the distance through the short
    block, so f is small but not zero. Scoring only the direct edge would
    report probability 0 on the reference and treat the pair as pure deletion
    evidence. With the optional ref gate on, concordant pairs do not increment
    skip-edge support, but the stored reads remain for likelihood.
    """
    lib = fit_library(list(range(360, 441)) * 10, read_len_short=200, read_len_long=200)
    blocks = [Block(0, 1000), Block(1000, 1030), Block(1030, 5000)]
    pair = PairRecord(800, 1000, 1030, 1230, False, True, 60, 60, 430)
    graph = build_graph(
        blocks,
        _copies(pair, 2),
        lib,
        CallConfig(min_edge_support=2, skip_edge_ref_gate=True),
    )
    novel = [edge for edge in graph.edges.values() if edge.key() == (1, 4)]
    assert len(novel) == 1
    assert novel[0].support == 0  # gate on: concordant under ref → no PE support
    # block A start 0 minus read start 800, plus read end 1230 minus block C start 1030.
    assert novel[0].reads[0].offset == -600
    index = [(novel[0], i) for i in range(len(novel[0].reads))]
    ref_path = reference_path(0, 2)
    deletion = (0, 1, 4, 5)
    ref_scored = score_path(ref_path, graph, [lib], index)
    del_scored = score_path(deletion, graph, [lib], index)
    assert all(p > 0 for p in ref_scored.read_p)
    assert sum(del_scored.read_p) > sum(ref_scored.read_p)


def test_deletion_duplication_and_inversion_from_pairs():
    """No breakpoints are supplied. Discovery and the graph have to find the event."""
    # Wide enough that a breakpoint CI of a few bases still leaves the
    # implied insert inside the middle 95% of f.
    lib = fit_library(list(range(360, 441)) * 10, read_len_short=200, read_len_long=200)
    background = []
    for i in range(20):
        start = 200 + i * 40
        background.append(
            PairRecord(start, start + 200, start + 200, start + 400, False, True, 60, 60, 400)
        )
    deletion = _copies(PairRecord(4800, 5000, 5500, 5700, False, True, 60, 60, 900), 8)
    inversion = _copies(PairRecord(9800, 10000, 10800, 11000, False, False, 60, 60, 400), 8)
    inversion += _copies(PairRecord(10000, 10200, 11000, 11200, True, True, 60, 60, 400), 8)

    def one(pairs, chrom_len):
        rows, logs = call_evidence(
            background + pairs,
            chrom="chrSim",
            chrom_len=chrom_len,
            cfg=CallConfig(),
            libs=[lib],
        )
        assert rows, logs
        return rows[0].as_dict()

    deleted = one(deletion, 20_000)
    assert deleted["svtype"] == "DEL"
    assert deleted["rearrangement"] == "AC"
    assert deleted["gt"] == "HOM"

    def duplication_call(outward_n: int) -> dict[str, str]:
        # Concordant pairs in the flanks are the depth that keeps a third copy
        # from beating two copies. Fewer outward pairs are the heterozygous case.
        local = []
        for i in range(80):
            start = 200 + (i % 40) * 40
            local.append(
                PairRecord(start, start + 200, start + 200, start + 400, False, True, 60, 60, 400, qname_hash=10_000 + i)
            )
        # Outward Dup: rev upstream of fwd → negative signed insert (not concordant).
        local.extend(_copies(PairRecord(2500, 2700, 3300, 3500, True, False, 60, 60, -600), outward_n))
        rows, logs = call_evidence(
            local,
            chrom="chrSim",
            chrom_len=6_000,
            cfg=CallConfig(),
            libs=[lib],
        )
        assert rows, logs
        return rows[0].as_dict()

    homozygous = duplication_call(12)
    assert homozygous["svtype"] == "DUP"
    assert homozygous["rearrangement"] == "ABBC"
    assert homozygous["gt"] == "HOM"
    assert homozygous["reference"].count("B") == 1
    # The tied higher-copy genotype is not the runner-up.
    assert float(homozygous["score_vs_next"]) > 0

    heterozygous = duplication_call(4)
    assert heterozygous["svtype"] == "DUP"
    assert heterozygous["rearrangement"] == "ABBC"
    assert heterozygous["gt"] == "HET"

    inverted = one(inversion, 30_000)
    assert inverted["svtype"] == "INV"
    assert inverted["rearrangement"] == "AB'C"


def test_copy_number_band_does_not_flip_other_types():
    ref = (0, 1, 2, 3, 4, 5)
    deletion = (0, 1, 4, 5)
    het = (-10.0, deletion, ref, 0.5)
    hom = (-10.3, deletion, deletion, 0.5)
    chosen, _runner = select_in_band([het, hom], ref)
    assert chosen[1] != chosen[2]

    two = (0, 1, 2, 3, 2, 3, 4, 5)
    three = (0, 1, 2, 3, 2, 3, 2, 3, 4, 5)
    higher = (-10.0, three, ref, 0.5)
    lower = (-10.2, two, two, 0.5)
    reference = (-10.1, ref, ref, 0.5)
    chosen, _runner = select_in_band([higher, lower, reference], ref)
    assert chosen[1] == chosen[2] == two

    # More than LIKELIHOOD_BAND behind the peak, the higher-copy genotype wins.
    far_lower = (-13.5, two, two, 0.5)
    chosen, _runner = select_in_band([higher, far_lower], ref)
    assert chosen[1] == three


def test_split_only_tandem_dup_reports_pe_zero_and_two_copies():
    """60 bp tandem with splits + concordant inside → pe 0; prefer ≤2 copies."""
    from arcsv2.engine.caller import _support_for_junction, build_graph
    from arcsv2.engine.evidence import SplitHit
    from arcsv2.engine.likelihood import _max_block_copies

    lib = fit_library(list(range(360, 441)) * 10, read_len_short=100, read_len_long=100)
    blocks = [Block(0, 1000), Block(1000, 1060), Block(1060, 5000)]
    inside = [
        PairRecord(1010, 1030, 1030, 1050, False, True, 60, 60, 400, qname_hash=i)
        for i in range(8)
    ]
    flanks = []
    for i in range(40):
        s = 200 + (i % 20) * 30
        flanks.append(
            PairRecord(s, s + 100, s + 100, s + 200, False, True, 60, 60, 400, qname_hash=1000 + i)
        )
    hits = [
        SplitHit(
            start1=1030,
            end1=1060,
            reverse1=False,
            start2=1000,
            end2=1030,
            reverse2=False,
            mapq=60,
            split_type="Dup+",
            bp1=(1060, 1060),
            bp2=(1000, 1000),
            chrom="chrSim",
        )
        for _ in range(6)
    ]
    graph = build_graph(blocks, flanks + inside, lib, CallConfig(), hits)
    sr, pe = _support_for_junction(graph, 1, 1)
    assert pe == 0
    assert sr == 6
    ref = reference_path(0, 2)
    two = (0, 1, 2, 3, 2, 3, 4, 5)
    three = (0, 1, 2, 3, 2, 3, 2, 3, 4, 5)
    call = genotype_locus([ref, two, three], ref, graph, [lib], pi=1e-4)
    assert call is not None
    assert max(_max_block_copies(call.path1), _max_block_copies(call.path2)) <= 2


def test_insertion_clip_is_reported_and_a_censored_size_is_a_bound():
    blocks = create_blocks([(4800, 5001)], 8_000)
    _stamp_clip_positions(
        blocks,
        [Breakpoint(4800, 5001, n_softclip=4, n_clip_left=2, n_clip_right=2, clip_pos=5000)],
    )
    assert blocks[0].end == 4800
    assert _insertion_position(blocks, 0, 1) == 5000
    assert insertion_length_is_bound(316, [31, 100, 316], 500)
    assert not insertion_length_is_bound(100, [31, 100, 316], 500)
    assert insertion_length_is_bound(3162, [31, 100, 316, 1000, 3162], 400)
    flagged = Block(0, 0, is_insertion=True, insertion_len=316, length_is_lower_bound=True)
    assert _format_inslen([flagged], {0: 316}) == ">316"


def test_novel_insertion_is_searched_at_a_two_sided_clip():
    inserts = [390, 395, 400, 400, 405, 410] * 40
    lib = fit_library(inserts, read_len_short=200, read_len_long=200)
    pairs = [
        PairRecord(4800, 5000, 0, 0, False, False, 60, 0, 0, hanging=True, qname_hash=i)
        for i in range(4)
    ]
    rows, logs = call_evidence(
        pairs,
        chrom="chrIns",
        chrom_len=8_000,
        cfg=CallConfig(min_affected_length=20),
        libs=[lib],
        breakpoints=[Breakpoint(5000, 5001, n_clip_left=2, n_clip_right=2)],
    )
    assert rows, logs
    row = rows[0].as_dict()
    assert row["svtype"] == "INS"
    assert "_" in row["rearrangement"]
    assert row["inslen"].startswith(">")
    inslen = int(row["inslen"].lstrip(">").split(",")[0])
    assert inslen >= 31
    assert int(row["len_affected"]) == inslen


def test_only_a_discordant_pair_spans_the_insertion():
    from arcsv2.engine.caller import _pair_spans_boundary
    from arcsv2.engine.graph import create_blocks

    lib = fit_library([400] * 40, read_len_short=200, read_len_long=200)
    blocks = create_blocks([(5000, 5001)], 8_000)
    concordant = PairRecord(4800, 5000, 5001, 5201, False, True, 60, 60, 400)
    assert not _pair_spans_boundary(build_graph(blocks, [concordant], lib, CallConfig()), 0)
    discordant = PairRecord(4000, 4200, 5600, 5800, False, True, 60, 60, 2000)
    assert _pair_spans_boundary(build_graph(blocks, [discordant], lib, CallConfig()), 0)
    hanging = PairRecord(4800, 5000, 0, 0, False, False, 60, 0, 0, hanging=True)
    assert not _pair_spans_boundary(build_graph(blocks, [hanging], lib, CallConfig()), 0)


def test_hanging_ignores_off_path_insertions():
    """Hanging likelihood must not scan novel insertions outside the path."""
    from arcsv2.engine.likelihood import _hanging_probability, _hap_blocks, _path_has_insertion

    blocks = [
        Block(0, 100),
        Block(100, 200),
        Block(0, 0, is_insertion=True, insertion_len=50),
    ]
    graph = Graph(blocks)
    path = (0, 1, 2, 3)  # reference only; insertion block 2 is off-path
    hap = _hap_blocks(blocks)
    assert not _path_has_insertion(path, hap)
    read = EdgeRead(
        offset=10,
        lib_id=0,
        pmapped=1.0,
        hanging=True,
        mapq=60,
        kind="pe",
        distant=False,
        mate_len=100,
    )
    edge = Edge(1, 0, 1, reads=[read])
    lib = fit_library([400] * 40, read_len_short=100, read_len_long=100)
    # Old full-graph any(is_insertion) was True; path-local check is False -> class_probs[1].
    got = _hanging_probability(edge, read, path, hap, lib)
    assert abs(got - lib.class_probs[1]) < 1e-15
    scored = score_path(path, graph, [lib], [(edge, 0)], hap=hap)
    assert abs(scored.read_p[0] - lib.class_probs[1]) < 1e-15


def test_path_layout_matches_block_distance_and_placements():
    """C5 genotype index must stay bit-identical to the linear helpers."""
    blocks = [
        Block(0, 100),
        Block(100, 250),
        Block(260, 400),
        Block(400, 550),
    ]
    path = (0, 1, 2, 3, 2, 3, 4, 5, 6, 7)
    layout = PathLayout(path, blocks)
    for i in range(len(path)):
        for j in range(i, len(path)):
            assert layout.distance(i, j) == block_distance(path, blocks, i, j)
    for v1, v2 in ((1, 2), (1, 4), (0, 5), (3, 6)):
        assert layout.placement_distances(v1, v2) == placement_distances(
            path, blocks, v1, v2
        )
        adj = (2, 3, 2)
        assert layout.placement_distances(v1, v2, adj, None) == placement_distances(
            path, blocks, v1, v2, adj, None
        )


def test_c5_noise_filter_flags_default_off():
    cfg = CallConfig()
    assert cfg.filter_high_nm is False
    assert cfg.filter_both_mates_clipped is False
    assert cfg.max_nm == 5
