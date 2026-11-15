"""SAM column 9, mate-unmapped pairs, soft clips, and split-read edges."""

from __future__ import annotations

import pytest

from arcsv2.config import CallConfig
from arcsv2.engine.caller import PairRecord, build_graph
from arcsv2.engine.evidence import parse_split, split_breakpoints
from arcsv2.engine.graph import create_blocks
from arcsv2.engine.library import fit_library
from arcsv2.samio import load_region


def _write(tmp_path, body: str) -> str:
    text = "@SQ\tSN:chr1\tLN:20000\n" + body
    if not text.endswith("\n"):
        text += "\n"
    path = tmp_path / "reads.sam"
    path.write_bytes(text.encode("utf-8"))
    return str(path)


def test_signed_insert_not_tlen_and_unmapped_mates_hang(tmp_path):
    """Insert is outer signed span (orig:process_insert_len), not abs(TLEN)."""
    seq = "A" * 100
    qual = "I" * 100
    path = _write(
        tmp_path,
        "\n".join(
            [
                f"r1\t0\tchr1\t100\t60\t100M\t=\t300\t400\t{seq}\t{qual}\tTLEN:i:999",
                f"r1\t16\tchr1\t300\t60\t100M\t=\t100\t-400\t{seq}\t{qual}",
                f"hang\t8\tchr1\t500\t60\t100M\t*\t0\t0\t{seq}\t{qual}",
            ]
        ),
    )
    pairs, _extras, _splits, chrom, _end, _clips = load_region(path)
    assert chrom == "chr1"
    normal = [pair for pair in pairs if not pair.hanging]
    hanging = [pair for pair in pairs if pair.hanging]
    # 0-based: fwd [99,199), rev [299,399) → 399 − 99 = 300 (TLEN column is ignored).
    assert len(normal) == 1 and normal[0].insert == 300
    assert len(hanging) == 1 and hanging[0].start1 == 499 and hanging[0].insert == 0


def test_signed_insert_length_matches_process_insert_len():
    """Outer FR insert is last.ref_end − first.ref_start (+ soft clips); Dup can be negative."""
    from arcsv2.samio import SamRead, insert_length

    def read(start: int, cigar: str, reverse: bool, seq_len: int = 100) -> SamRead:
        return SamRead(
            qname="q",
            flag=16 if reverse else 0,
            rname="chr1",
            start=start,
            mapq=60,
            cigar=cigar,
            seq="A" * seq_len,
            qual="I" * seq_len,
            tlen=0,
            tags={},
        )

    # FR non-overlap: fwd 100–200, rev 300–400 → 400 − 100 = 300.
    assert insert_length(read(100, "100M", False), read(300, "100M", True)) == 300
    # Outward Dup: rev 100–200, fwd 300–400 → 200 − 300 = −100.
    assert insert_length(read(300, "100M", False), read(100, "100M", True)) == -100
    # Soft clips shrink aligned span (395 − 100 = 295) then + first 5' (10) + last 3' (5).
    assert insert_length(read(100, "10S90M", False, 100), read(300, "95M5S", True, 100)) == 310
    # Same orientation → None.
    assert insert_length(read(100, "100M", False), read(300, "100M", False)) is None


def test_region_window_keeps_pairs_that_overlap_it(tmp_path):
    seq = "A" * 100
    qual = "I" * 100
    path = _write(
        tmp_path,
        "\n".join(
            [
                f"mid\t0\tchr1\t1500\t60\t100M\t=\t1700\t400\t{seq}\t{qual}",
                f"mid\t16\tchr1\t1700\t60\t100M\t=\t1500\t-400\t{seq}\t{qual}",
                f"far\t0\tchr1\t5000\t60\t100M\t=\t5200\t400\t{seq}\t{qual}",
                f"far\t16\tchr1\t5200\t60\t100M\t=\t5000\t-400\t{seq}\t{qual}",
            ]
        ),
    )
    pairs, _extras, _splits, _chrom, _end, _clips = load_region(path, start=999, end=2000)
    assert len(pairs) == 1
    assert pairs[0].start1 == 1499


def test_softclips_become_breakpoints(tmp_path):
    seq = "A" * 20
    qual = "I" * 20
    path = _write(
        tmp_path,
        "\n".join(
            [
                f"s1\t0\tchr1\t2000\t40\t8S12M\t*\t0\t0\t{seq}\t{qual}",
                f"s2\t0\tchr1\t2000\t40\t8S12M\t*\t0\t0\t{seq}\t{qual}",
            ]
        ),
    )
    _pairs, extras, splits, _chrom, _end, _clips = load_region(path)
    assert splits == []
    assert len(extras) == 1
    assert extras[0].start == 1999
    assert extras[0].n_softclip == 2


def test_deletion_split_supports_the_junction():
    cfg = CallConfig()
    hit = parse_split(
        chrom="chr1",
        start=4800,
        mapq=60,
        cigar="100M100S",
        reverse=False,
        sa="chr1,5501,+,100S100M,60,0;",
        cfg=cfg,
    )
    assert hit is not None
    assert hit.split_type == "Del+"
    assert hit.bp1 == (4900, 4900)
    assert hit.bp2 == (5500, 5500)
    blocks = create_blocks([(bp.start, bp.end) for bp in split_breakpoints([hit])], 10_000)
    lib = fit_library([400] * 50, read_len_short=100, read_len_long=100)
    graph = build_graph(blocks, [], lib, cfg, [hit, hit])
    split_edges = [edge for edge in graph.edges.values() if any(read.kind == "split" for read in edge.reads)]
    assert len(split_edges) == 1
    assert split_edges[0].key() == (1, 4)
    assert split_edges[0].support == 2


def test_same_block_dup_split_makes_tandem_edge():
    """Dup+ both segments on one block → B.out↔B.in (orig:block_seq_to_path flip)."""
    from arcsv2.engine.evidence import SplitHit
    from arcsv2.engine.graph import Block, in_node, out_node

    cfg = CallConfig()
    lib = fit_library([400] * 50, read_len_short=100, read_len_long=100)
    # Block B = [1000, 1200]; first ends at B.end, second starts at B.start.
    blocks = [Block(0, 1000), Block(1000, 1200), Block(1200, 5_000)]
    b_idx = 1
    hit = SplitHit(
        start1=1100,
        end1=1200,
        reverse1=False,
        start2=1000,
        end2=1100,
        reverse2=False,
        mapq=60,
        split_type="Dup+",
        bp1=(1200, 1200),
        bp2=(1000, 1000),
        chrom="chr1",
    )
    graph = build_graph(blocks, [], lib, cfg, [hit, hit, hit])
    want = tuple(sorted((out_node(b_idx), in_node(b_idx))))
    split_edges = [
        edge
        for edge in graph.edges.values()
        if any(read.kind == "split" for read in edge.reads)
    ]
    assert len(split_edges) == 1
    assert split_edges[0].key() == want
    assert split_edges[0].support == 3


def test_hanging_read_stays_on_the_block_edge():
    lib = fit_library([400] * 50, read_len_short=100, read_len_long=100)
    blocks = create_blocks([], 5_000)
    hang = PairRecord(100, 200, 0, 0, False, False, 60, 0, 0, hanging=True)
    graph = build_graph(blocks, [hang], lib, CallConfig())
    stored = [read for edge in graph.edges.values() for read in edge.reads if read.hanging]
    assert len(stored) == 1
    assert graph.neighbors(0, 2) == []


def test_distant_and_unmapped_mates_are_hanging(tmp_path):
    seq = "A" * 100
    short = "A" * 80
    qual = "I" * 100
    short_qual = "I" * 80
    path = _write(
        tmp_path,
        "\n".join(
            [
                f"far\t0\tchr1\t100\t60\t100M\t=\t3000100\t0\t{seq}\t{qual}",
                f"far\t16\tchr1\t3000100\t60\t100M\t=\t100\t0\t{seq}\t{qual}",
                f"cross\t0\tchr1\t200\t60\t100M\tchr2\t100\t0\t{seq}\t{qual}",
                f"cross\t16\tchr2\t100\t60\t100M\tchr1\t200\t0\t{seq}\t{qual}",
                f"um\t0\tchr1\t400\t60\t100M\t*\t0\t0\t{seq}\t{qual}",
                f"um\t4\t*\t0\t0\t*\tchr1\t400\t0\t{short}\t{short_qual}",
            ]
        ),
    )
    pairs, _extras, _splits, chrom, _end, _clips = load_region(path, chrom="chr1")
    assert chrom == "chr1"
    assert not any(pair for pair in pairs if not pair.hanging)
    distant = [pair for pair in pairs if pair.distant]
    unmapped = [pair for pair in pairs if pair.hanging and not pair.distant]
    # Each end of the 3 Mb pair is its own distant hanging read. The chr2 mate is not.
    assert len(distant) == 3
    assert {pair.start1 for pair in distant} == {99, 199, 3_000_099}
    assert all(pair.mate_len == 100 for pair in distant)
    assert len(unmapped) == 1
    assert unmapped[0].start1 == 399
    assert unmapped[0].mate_len == 80


def _pair_key(pair):
    return (
        pair.start1,
        pair.end1,
        pair.start2,
        pair.end2,
        pair.reverse1,
        pair.reverse2,
        pair.mapq1,
        pair.mapq2,
        pair.insert,
        pair.hanging,
        pair.distant,
        pair.mate_len,
    )


@pytest.mark.io
def test_sam_and_bam_give_identical_calls(tmp_path):
    pysam = pytest.importorskip("pysam")
    from arcsv2.engine.caller import call_evidence

    seq = "A" * 100
    qual = "I" * 100
    body = "\n".join(
        [
            f"near\t99\tchr1\t1200\t60\t100M\t=\t1500\t400\t{seq}\t{qual}",
            f"near\t147\tchr1\t1500\t60\t100M\t=\t1200\t-400\t{seq}\t{qual}",
            f"far\t99\tchr1\t5000\t60\t100M\t=\t5300\t400\t{seq}\t{qual}",
            f"far\t147\tchr1\t5300\t60\t100M\t=\t5000\t-400\t{seq}\t{qual}",
        ]
    )
    sam_path = _write(tmp_path, body)
    bam_path = tmp_path / "reads.bam"
    with pysam.AlignmentFile(sam_path, "r") as sam, pysam.AlignmentFile(
        str(bam_path), "wb", header=sam.header
    ) as bam:
        for record in sam:
            bam.write(record)
    pysam.index(str(bam_path))

    def _loaded(path: str, start: int | None = None, end: int | None = None):
        pairs, breakpoints, splits, chrom, chrom_end, softclips = load_region(
            path, chrom="chr1", start=start, end=end
        )
        rows, _logs = call_evidence(pairs, chrom=chrom, chrom_len=max(chrom_end, 20_000))
        return pairs, breakpoints, splits, [row.as_dict() for row in rows]

    sam_all = _loaded(sam_path)
    bam_all = _loaded(str(bam_path))
    assert [_pair_key(pair) for pair in sam_all[0]] == [_pair_key(pair) for pair in bam_all[0]]
    assert len(sam_all[0]) == 2
    assert sam_all[3] == bam_all[3]
    sam_window = _loaded(sam_path, 999, 2000)
    bam_window = _loaded(str(bam_path), 999, 2000)
    assert len(sam_window[0]) == 1
    assert [_pair_key(pair) for pair in sam_window[0]] == [_pair_key(pair) for pair in bam_window[0]]
    assert sam_window[3] == bam_window[3]


@pytest.mark.io
def test_cram_without_a_reference_is_rejected(tmp_path):
    pytest.importorskip("pysam")
    path = tmp_path / "reads.cram"
    path.write_bytes(b"not a cram")
    with pytest.raises(ValueError, match="reference"):
        load_region(str(path))
