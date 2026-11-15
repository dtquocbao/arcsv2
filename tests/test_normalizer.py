"""Gap and insertion normalizer against orig:sv_inference_insertions.py tests."""

from arcsv2.engine.library import cdf_and_cdf_sum
from arcsv2.engine.normalizer import (
    HapBlock,
    compute_normalizing_constant,
    genome_blocks_gaps,
    get_gap_overlap_positions,
    get_insertion_overlap_positions,
)


def _opens(intervals) -> list[tuple[int, int]]:
    return [(iv.lo, iv.hi) for iv in intervals]


def _flat_cdf_sum(n_unit: int) -> list[float]:
    """``[0] + [1/n]*n``, the insert law in test_compute_normalizing_constant."""
    pmf = [0.0] + [1.0 / n_unit] * n_unit
    _cdf, cdf_sum = cdf_and_cdf_sum(pmf)
    return cdf_sum


def test_gap_overlap_matches_original_truth():
    read_len = 50
    blocks = [
        HapBlock(0, 100),
        HapBlock(100, 200),
        HapBlock(249, 300),
        HapBlock(350, 400),
        HapBlock(500, 600),
    ]
    paths = (
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
        [0, 1, 2, 3, 4, 5, 7, 6, 8, 9],
    )
    truth = ([(299, 301), (399, 451)], [(299, 326), (424, 451)])
    for path, expected in zip(paths, truth):
        assert _opens(get_gap_overlap_positions(path, blocks, read_len)) == expected

    blocks = [
        HapBlock(0, 100),
        HapBlock(200, 300),
        HapBlock(350, 400),
        HapBlock(0, 50, is_de_novo=True),
        HapBlock(0, 50, is_de_novo=True),
    ]
    path = [0, 1, 6, 7, 2, 3, 8, 9, 4, 5]
    assert _opens(get_gap_overlap_positions(path, blocks, read_len)) == [
        (99, 131),
        (169, 201),
        (349, 356),
        (394, 401),
    ]


def test_insertion_overlap_matches_original_truth():
    blocks = [
        HapBlock(0, 100),
        HapBlock(100, 200),
        HapBlock(210, 300),
        HapBlock(350, 360),
        HapBlock(370, 400),
        HapBlock(0, 100, is_de_novo=True),
        HapBlock(0, 10, is_de_novo=True),
    ]
    paths = (
        list(range(10)),
        [0, 1, 10, 11, 2, 3],
        [0, 1, 2, 3, 10, 11, 2, 3],
        [0, 1, 2, 3, 12, 13, 2, 3],
        [0, 1, 2, 3, 4, 5, 10, 11, 6, 7],
        [0, 1, 2, 3, 4, 5, 12, 13, 6, 7],
    )
    truth = [
        [],
        [(80, 170)],
        [(185, 275)],
        [],
        [(305, 395)],
        [],
    ]
    for path, expected in zip(paths, truth):
        found, _de_novo, _trans = get_insertion_overlap_positions(path, blocks, 50, 20)
        assert _opens(found) == expected

    blocks = [
        HapBlock(0, 100),
        HapBlock(200, 300),
        HapBlock(350, 400),
        HapBlock(0, 50, is_de_novo=True),
        HapBlock(0, 50, is_de_novo=True),
    ]
    path = [0, 1, 6, 7, 2, 3, 8, 9, 4, 5]
    found, _, _ = get_insertion_overlap_positions(path, blocks, 50, 20)
    assert _opens(found) == [(130, 170), (355, 395)]


def test_normalizer_gap_and_insertion_cases():
    cdf_sum = _flat_cdf_sum(200)
    classes = [1.0, 0.0, 0.0, 0.0]
    contiguous = [HapBlock(0, 1000), HapBlock(1000, 2000)]
    nc1 = compute_normalizing_constant(list(range(4)), contiguous, cdf_sum, classes, 100, 100)

    gapped = [HapBlock(0, 1000), HapBlock(1100, 2000)]
    nc3 = compute_normalizing_constant(list(range(4)), gapped, cdf_sum, classes, 100, 100)
    assert 0 < nc3 < nc1

    filled = [
        HapBlock(0, 1000),
        HapBlock(1000, 1940),
        HapBlock(0, 60, is_de_novo=True),
    ]
    nc_filled = compute_normalizing_constant([0, 1, 4, 5, 2, 3], filled, cdf_sum, classes, 100, 100)
    assert nc_filled == nc1

    short = [
        HapBlock(0, 1000),
        HapBlock(1000, 1938),
        HapBlock(0, 62, is_de_novo=True),
    ]
    nc_short = compute_normalizing_constant([0, 1, 4, 5, 2, 3], short, cdf_sum, classes, 100, 100)
    assert 0 < nc_short < nc1


def test_median_gap_is_inserted_between_blocks():
    blocks = [
        HapBlock(0, 100),
        HapBlock(105, 200),
        HapBlock(200, 300),
        HapBlock(305, 400),
        HapBlock(420, 500),
    ]
    pieces = genome_blocks_gaps(blocks, [0, 1, 4, 5, 8, 9])
    assert [len(block) for block in pieces] == [0, 100, 2, 100, 12, 80, 0]
