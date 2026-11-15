"""Synthetic assembly labels. Alignment cases skip when mappy is missing."""

from __future__ import annotations

import math
import random
import sys
from pathlib import Path

import pytest

from arcsv2.label import (
    LABEL_VERSION,
    MAX_LABEL_QUERY_LEN,
    Hit,
    ReferenceView,
    _blocks_point_from_junctions,
    _midpoint_candidates,
    _passes,
    blocks_from_row,
    dataset_tables,
    haplotype_sequences,
    label_call,
    label_from_passes,
    label_rows,
    queries_for,
)

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from simulate import apply_template  # noqa: E402


def test_blocks_from_row_keeps_uncovered_breakpoint_gaps():
    row = {
        "reference": "ABC",
        "blocks": "0-5000,5001-5500,5501-20000",
    }
    assert blocks_from_row(row) == {
        "A": (0, 5000),
        "B": (5001, 5500),
        "C": (5501, 20_000),
    }
    insertion = {"reference": "A_B", "blocks": "0-100,_,100-250"}
    assert blocks_from_row(insertion) == {"A": (0, 100), "B": (100, 250)}


def test_midpoint_candidates_cover_interval_edges():
    cands = _midpoint_candidates(536_794, 536_885)
    assert cands[0] == 536_794 and 536_884 in cands


def test_blocks_point_can_shift_a_midpoint_junction():
    row = {
        "blocks": "0-100,200-300",
        "blocks_point": "0-150,150-300",
        "bp_point_source": "mid",
    }
    assert _blocks_point_from_junctions(row, [160]) == "0-160,160-300"


def test_label_sequence_uses_the_point_and_150_bp_flanks():
    chrom = "ACGT" * 20_000
    row = {
        "chrom": "chr1",
        "reference": "ABC",
        "rearrangement": "AC",
        "blocks": "0-30000,33000-36000,40000-80000",
        "blocks_point": "0-33333,33333-36000,36000-80000",
    }
    alt, ref = haplotype_sequences(row, {"chr1": chrom}, flank=150)
    left = chrom[33_333 - 150 : 33_333]
    middle = chrom[33_333:36_000]
    right = chrom[36_000 : 36_000 + 150]
    assert ref == left + middle + right
    assert alt == left + right
    assert len(left) == 150 and len(right) == 150
    assert chrom[30_000:33_000] not in ref


def test_novel_sequence_and_unreliable_windows_are_missing():
    missing = label_call(
        {"rearrangement": "A_B", "reference": "AB", "blocks": "0-10,_,10-20", "chrom": "chrSim"},
        {"chrSim": "ACGT" * 10},
        None,
    )
    assert math.isnan(missing.label) and missing.reason == "no_sequence"
    unreliable = label_call(
        {
            "chrom": "chrSim",
            "reference": "ABC",
            "rearrangement": "AC",
            "blocks": "100-200,200-400,400-500",
            "svtype": "DEL",
            "gt": "HET",
        },
        {"chrSim": "A" * 800},
        None,
        {"chrSim": [(150, 160)]},
    )
    assert math.isnan(unreliable.label) and unreliable.reason == "unreliable_assembly"


def test_megabase_ref_query_is_skipped():
    """Multi-block loci can build multi-Mb ref alleles; skip before fetching FASTA."""
    flank = 150
    # Outer blocks trimmed to flank; one huge internal block makes ref too long.
    internal = MAX_LABEL_QUERY_LEN
    left_end = flank
    mid_end = left_end + internal
    right_end = mid_end + flank
    # Dict FASTA only needs to cover coords if we wrongly fetch; keep it tiny so
    # a regression that builds the sequence fails loudly.
    chrom = "A" * (2 * flank + 10)
    row = {
        "chrom": "chrSim",
        "reference": "ABC",
        "rearrangement": "AC",
        "blocks": f"0-{left_end},{left_end}-{mid_end},{mid_end}-{right_end}",
        "blocks_point": f"0-{left_end},{left_end}-{mid_end},{mid_end}-{right_end}",
    }
    early, queries = queries_for(row, {"chrSim": chrom}, flank=flank)
    assert early is not None and math.isnan(early.label)
    assert early.reason == "query_too_long"
    assert queries == {}
    # Just under the cap still builds queries (no aligner needed for this check).
    internal_ok = MAX_LABEL_QUERY_LEN - 2 * flank - 1
    mid_end_ok = left_end + internal_ok
    right_end_ok = mid_end_ok + flank
    chrom_ok = "A" * right_end_ok
    row_ok = {
        "chrom": "chrSim",
        "reference": "ABC",
        "rearrangement": "AC",
        "blocks": f"0-{left_end},{left_end}-{mid_end_ok},{mid_end_ok}-{right_end_ok}",
        "blocks_point": f"0-{left_end},{left_end}-{mid_end_ok},{mid_end_ok}-{right_end_ok}",
    }
    early_ok, queries_ok = queries_for(row_ok, {"chrSim": chrom_ok}, flank=flank)
    assert early_ok is None
    assert len(queries_ok["ref"]) <= MAX_LABEL_QUERY_LEN


def test_two_megabase_call_is_skipped_before_sequence_build():
    """Synthetic 2 Mb call must NA as query_too_long without needing a 2 Mb FASTA."""
    flank = 150
    span = 2_000_000
    row = {
        "chrom": "chrSim",
        "reference": "ABC",
        "rearrangement": "AC",
        "blocks": f"0-{flank},{flank}-{flank + span},{flank + span}-{2 * flank + span}",
        "blocks_point": f"0-{flank},{flank}-{flank + span},{flank + span}-{2 * flank + span}",
    }
    early, queries = queries_for(row, {"chrSim": "A"}, flank=flank)
    assert early is not None
    assert early.reason == "query_too_long"
    assert queries == {}


def test_either_strand_passes_and_a_heterozygote_is_not_ambiguous():
    hit = Hit("hap1", 0, 80, 80, 60, 80, 80, "80M", -1)
    ok, reason = _passes(hit)
    assert ok and reason == "validated"
    label, why, ref_ok = label_from_passes(
        {"hap1": True, "hap2": False},
        {"hap1": False, "hap2": True},
    )
    assert label == 1.0 and why == "validated" and ref_ok
    label, why, _ref_ok = label_from_passes(
        {"hap1": True, "hap2": True},
        {"hap1": False, "hap2": False},
    )
    assert label == 1.0 and why == "validated"
    label, why, _ref_ok = label_from_passes(
        {"hap1": True, "hap2": False},
        {"hap1": True, "hap2": False},
    )
    assert math.isnan(label) and why == "ambiguous"


def test_reference_window_is_opened_once(tmp_path, monkeypatch):
    path = tmp_path / "ref.fa"
    _fasta(path, "chrSim", "ACGT" * 50)
    opened = {"n": 0}
    original = ReferenceView.__init__

    def _counting(self, fasta):
        opened["n"] += 1
        original(self, fasta)

    monkeypatch.setattr(ReferenceView, "__init__", _counting)
    rows = [
        {"chrom": "chrSim", "rearrangement": "A_B", "reference": "AB", "blocks": "0-10,_,10-20"},
        {"chrom": "chrSim", "rearrangement": "A_B", "reference": "AB", "blocks": "0-10,_,10-20"},
    ]
    label_rows(rows, path, None, "sim")
    assert opened["n"] == 1


def test_dataset_drops_missing_labels_and_refuses_mixed_versions():
    kept, summary = dataset_tables(
        [
            [
                {
                    "chrom": "chr20",
                    "sample": "HG002",
                    "label": "1",
                    "reason": "validated",
                    "sv_class": "deletion",
                    "label_version": LABEL_VERSION,
                    "feature_schema": "arcsv2-f1",
                },
                {
                    "chrom": "chr20",
                    "sample": "HG002",
                    "label": "NA",
                    "reason": "ambiguous",
                    "sv_class": "deletion",
                    "label_version": LABEL_VERSION,
                    "feature_schema": "arcsv2-f1",
                },
            ]
        ]
    )
    assert len(kept) == 1 and kept[0]["group"] == "HG002"
    assert summary["classes"]["deletion"] == {"1": 1, "0": 0, "NA": 1, "reasons": {"validated": 1, "ambiguous": 1}}
    with pytest.raises(ValueError, match="mixed"):
        dataset_tables(
            [
                [{"label": "1", "label_version": LABEL_VERSION, "feature_schema": "arcsv2-f1"}],
                [{"label": "1", "label_version": "other", "feature_schema": "arcsv2-f1"}],
            ]
        )


def _fasta(path: Path, name: str, seq: str) -> None:
    path.write_text(f">{name}\n{seq}\n", encoding="utf-8")


def _call(blocks: str, rearrangement: str, svtype: str, gt: str) -> dict[str, str]:
    return {
        "chrom": "chrSim",
        "svtype": svtype,
        "gt": gt,
        "reference": "ABC",
        "rearrangement": rearrangement,
        "blocks": blocks,
    }


def _table(rows: list[dict[str, str]]) -> str:
    columns: list[str] = []
    for row in rows:
        for key in row:
            if key not in columns:
                columns.append(key)
    lines = ["\t".join(columns)]
    lines.extend("\t".join(row.get(col, "NA") for col in columns) for row in rows)
    return "\n".join(lines) + "\n"


def _revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


@pytest.mark.label
def test_synthetic_haplotypes_label_true_and_false_calls(tmp_path):
    pytest.importorskip("mappy")
    from arcsv2.label import HapAligner

    rng = random.Random(204178949)
    ref = "".join(rng.choice("ACGT") for _ in range(50_000))
    hap1 = apply_template(ref, "DEL", 10_000, 20_000)
    hap2 = apply_template(ref, "INV", 12_000, 18_000)
    ref_path = tmp_path / "ref.fa"
    hap1_path = tmp_path / "hap1.fa"
    hap2_path = tmp_path / "hap2.fa"
    _fasta(ref_path, "chrSim", ref)
    _fasta(hap1_path, "chrSim", hap1)
    _fasta(hap2_path, "chrSim", hap2)
    aligner = HapAligner(hap1_path, hap2_path)
    fasta = {"chrSim": ref}
    deletion = _call("9850-10000,10000-20000,20000-20150", "AC", "DEL", "HET")
    inversion = _call("11850-12000,12000-18000,18000-18150", "AB'C", "INV", "HET")
    shifted = _call("10350-10500,10500-20500,20500-20650", "AC", "DEL", "HET")
    novel = _call("0-100,_,100-200", "A_C", "INS", "HET")

    labeled = label_rows([deletion, inversion, shifted, novel], fasta, aligner, "sim")
    again = label_rows([deletion, inversion, shifted, novel], fasta, aligner, "sim")
    assert _table(labeled) == _table(again)

    assert labeled[0]["label"] == "1"
    assert labeled[0]["reason"] == "validated"
    assert labeled[0]["gt_concordant"] == "1"
    assert labeled[0]["ref_validated"] == "0"
    assert labeled[1]["label"] == "1"
    assert labeled[1]["gt_concordant"] == "1"
    assert labeled[2]["label"] == "0"
    assert labeled[2]["reason"] in ("large_gap", "clipped")
    assert labeled[3]["label"] == "NA" and labeled[3]["reason"] == "no_sequence"

    ref_hap = HapAligner(ref_path, ref_path)
    repeat = _call("1000-1500,1500-2500,2500-3000", "ABC", "DUP", "HOM")
    ambiguous = label_call(repeat, fasta, ref_hap)
    assert math.isnan(ambiguous.label) and ambiguous.reason == "ambiguous"
    assert ambiguous.ref_validated
    hit = ref_hap.best_hit(ref[1000:3000])
    assert hit is not None and hit.mapq >= 60 and hit.q_start == 0

    deleted = apply_template(ref, "DEL", 10_000, 20_000)
    het_fwd = HapAligner(hap1_path, ref_path)
    het_row = label_rows([deletion], fasta, het_fwd, "sim")[0]
    assert het_row["label"] == "1" and het_row["reason"] == "validated"
    assert het_row["gt_concordant"] == "1" and het_row["ref_validated"] == "1"

    rev_hap = tmp_path / "rev.fa"
    rev_ref = tmp_path / "revref.fa"
    _fasta(rev_hap, "chrSim", _revcomp(deleted))
    _fasta(rev_ref, "chrSim", _revcomp(ref))
    hom_rev = label_rows([_call(deletion["blocks"], "AC", "DEL", "HOM")], fasta, HapAligner(rev_hap, rev_hap), "sim")[0]
    assert hom_rev["label"] == "1" and hom_rev["gt_concordant"] == "1"
    assert hom_rev["strand"] == "-1"
    het_rev = label_rows([deletion], fasta, HapAligner(rev_hap, rev_ref), "sim")[0]
    assert het_rev["label"] == "1" and het_rev["reason"] == "validated"
    assert het_rev["gt_concordant"] == "1" and het_rev["strand"] == "-1"

    joint = label_rows([deletion, inversion], fasta, HapAligner(hap1_path, hap2_path), "sim")
    sequential = label_rows(
        [deletion, inversion],
        fasta,
        HapAligner(hap1_path, hap2_path, sequential=True),
        "sim",
    )
    assert [row["label"] for row in sequential] == [row["label"] for row in joint]
    assert [row["reason"] for row in sequential] == [row["reason"] for row in joint]


@pytest.mark.label
def test_sequential_loads_each_haplotype_once(tmp_path, monkeypatch):
    pytest.importorskip("mappy")
    import mappy

    from arcsv2.label import HapAligner

    rng = random.Random(204178949)
    ref = "".join(rng.choice("ACGT") for _ in range(8_000))
    hap = apply_template(ref, "DEL", 2_000, 4_000)
    ref_path = tmp_path / "ref.fa"
    hap_path = tmp_path / "hap.fa"
    _fasta(ref_path, "chrSim", ref)
    _fasta(hap_path, "chrSim", hap)
    loads = {"n": 0}
    real = mappy.Aligner

    def counting_aligner(*args, **kwargs):
        loads["n"] += 1
        return real(*args, **kwargs)

    monkeypatch.setattr(mappy, "Aligner", counting_aligner)
    rows = [
        _call("1850-2000,2000-4000,4000-4150", "AC", "DEL", "HET"),
        _call("1850-2000,2000-4000,4000-4150", "AC", "DEL", "HOM"),
        _call("0-100,_,100-200", "A_C", "INS", "HET"),
    ] * 3
    fasta = {"chrSim": ref}
    label_rows(rows, fasta, HapAligner(hap_path, hap_path, sequential=True), "sim")
    assert loads["n"] == 2
