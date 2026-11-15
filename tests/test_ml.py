"""Features, labels, merge-filter, and VCF."""

from pathlib import Path

import pytest

from arcsv2.beds import read_bed_by_chrom
from arcsv2.config import FEATURE_SCHEMA
from arcsv2.features import FEATURE_SCHEMA as SCHEMA
from arcsv2.features import featurize_row, fit_pca
from arcsv2.label import is_validated, largest_indel, sequence_from_call
from arcsv2.merge_filter import merge_calls, passes_filter
from arcsv2.predict import predict_rows
from arcsv2.train import choose_threshold
from arcsv2.vcf import format_record, vcf_header, write_vcf_gz


def _row(**kwargs):
    base = {
        "chrom": "chr1",
        "minbp": "1000",
        "maxbp": "1500",
        "svtype": "DEL",
        "complextype": "NA",
        "num_sv": "1",
        "bp": "0,1000,1500,4000",
        "bp_uncertainty": "0,0,0,0",
        "reference": "ABC",
        "rearrangement": "AC",
        "len_affected": "500",
        "gt": "HET",
        "af": "0.500",
        "inslen": "NA",
        "sr_support": "4",
        "pe_support": "8",
        "sv_bp": "1000,1500",
        "sv_bp_uncertainty": "0,0",
        "score_vs_ref": "20",
        "score_vs_next": "5",
        "rearrangement_next": "ABC/AC",
        "num_paths": "3",
        "filter": "PASS",
        "id": "chr1_0-4000",
    }
    base.update(kwargs)
    return base


def test_features_are_finite_for_three_classes():
    rows = [
        _row(),
        _row(svtype="BND,BND", complextype="ddup.del", reference="ABCDE", rearrangement="ADCDE", num_sv="2"),
        _row(svtype="BND", complextype="invddup.inv", reference="ABCD", rearrangement="AB'C'D", num_sv="2"),
    ]
    pca = fit_pca(rows)
    for row in rows:
        feats = featurize_row(row, pca)
        assert all(v == v and abs(v) != float("inf") for v in feats.values())
        assert SCHEMA == FEATURE_SCHEMA
        assert "len_affected" in feats
        assert "rearr_pc1" in feats


def test_label_reasons():
    ok, reason = is_validated(q_start=0, q_end=100, q_len=100, mapq=60, mlen=95, blen=100, cigar="100M")
    assert ok and reason == "validated"
    assert is_validated(q_start=1, q_end=100, q_len=100, mapq=60, mlen=95, blen=100, cigar="100M")[1] == "clipped"
    assert is_validated(q_start=0, q_end=100, q_len=100, mapq=10, mlen=95, blen=100, cigar="100M")[1] == "low_mapq"
    assert is_validated(q_start=0, q_end=100, q_len=100, mapq=60, mlen=80, blen=100, cigar="100M")[1] == "low_identity"
    assert is_validated(q_start=0, q_end=100, q_len=100, mapq=60, mlen=95, blen=100, cigar="50M60D50M")[1] == "large_gap"
    assert largest_indel("10M49D5M") == 49
    assert sequence_from_call("AB", "A_B", {"A": "ACGT", "B": "TTTT"}, "AA", "CC") is None


def test_merge_filter_and_allele_frequency():
    ok, reason = passes_filter(0, 1000, 200)
    assert ok and reason == "PASS"
    assert passes_filter(0, 100, 40)[1] == "short"
    repeats = [(0, 900)]
    assert passes_filter(0, 1000, 200, simple_repeats=repeats)[1] == "simple_repeat"
    assert passes_filter(0, 1000, 200, segdups=[(0, 600)])[1] == "segdup"
    merged = merge_calls(
        [
            {"chrom": "chr1", "start": 0, "end": 100, "sv_class": "DEL", "gt": "HET", "sample": "a"},
            {"chrom": "chr1", "start": 10, "end": 110, "sv_class": "DEL", "gt": "HOM", "sample": "b"},
            {"chrom": "chr1", "start": 5000, "end": 5100, "sv_class": "DEL", "gt": "HET", "sample": "a"},
        ],
        n_samples=2,
    )
    assert len(merged) == 2
    assert merged[0]["allele_count"] == 3
    assert merged[0]["allele_frequency"] == 3 / 4
    repeated = merge_calls(
        [
            {"chrom": "chr1", "start": 0, "end": 100, "sv_class": "DEL", "gt": "HET", "sample": "a"},
            {"chrom": "chr1", "start": 5, "end": 100, "sv_class": "DEL", "gt": "HET", "sample": "a"},
            {"chrom": "chr1", "start": 0, "end": 100, "sv_class": "DEL", "gt": "HOM", "sample": "a"},
        ],
        n_samples=1,
    )
    assert repeated[0]["allele_count"] == 2
    assert repeated[0]["n_samples_called"] == 1


def test_threshold_keeps_recall():
    y = [1, 1, 1, 1, 1, 0, 0, 0, 0, 0]
    scores = [0.95, 0.9, 0.85, 0.8, 0.7, 0.6, 0.4, 0.3, 0.2, 0.1]
    choice = choose_threshold(y, scores, min_recall=0.90)
    assert choice.recall >= 0.90
    assert choice.precision >= 0.8


def test_predict_refuses_bad_schema_and_keeps_insertions():
    manifest = {"feature_schema": FEATURE_SCHEMA, "thresholds": {"deletion": 0.5}, "models": {"deletion": "paper"}}
    rows = predict_rows([_row(), _row(svtype="INS", rearrangement="A_C")], manifest, pca=None)
    assert rows[0]["ml_model"] == "untrained"
    assert rows[0]["ml_prob"] == "NA"
    assert rows[1]["ml_prob"] == "NA"
    assert rows[1]["ml_model"] == "none"
    try:
        predict_rows([_row()], {"feature_schema": "other", "thresholds": {}}, None)
        raise AssertionError("schema mismatch must fail")
    except ValueError as exc:
        assert "FEATURE_SCHEMA" in str(exc)


def test_bed_intervals_stay_on_their_chromosome(tmp_path):
    bed = tmp_path / "tracks.bed"
    bed.write_text("chr1\t0\t100\nchr2\t0\t50\n", encoding="utf-8")
    grouped = read_bed_by_chrom(bed)
    assert grouped["chr1"] == [(0, 100)]
    assert grouped["chr2"] == [(0, 50)]


def test_breakend_brackets_follow_orientation():
    # Deletion AC: forward blocks, orig orients ('-', '+') → N[mate[ and ]mate]N.
    lines = format_record(_row(), "evt")
    assert "N[chr1:1501[" in lines[1]
    assert "]chr1:1000]N" in lines[2]
    # Inversion A to B': both breakends use ']' and the '-' coordinate.
    inv = format_record(
        _row(rearrangement="AB'C", reference="ABC", sv_bp="100,200;200,300", svtype="BND,BND"),
        "inv",
    )
    assert "N]chr1:200]" in inv[1]
    assert "N]chr1:100]" in inv[2]


def test_fitted_logistic_replaces_the_score_proxy():
    pytest.importorskip("sklearn")
    from arcsv2.train import fit_labeled_calls

    rows = []
    for i, label in enumerate([1, 1, 1, 1, 0, 0, 0, 0]):
        score = 40 if label else 1
        rows.append(_row(label=str(label), score_vs_ref=str(score), id=f"c{i}"))
    manifest = fit_labeled_calls(rows)
    assert manifest["fitted"] is True
    assert "deletion" in manifest["pca"]
    spec = manifest["classifiers"]["deletion"]
    assert spec["mean"] and spec["scale"] and spec["pca"]
    assert manifest["grouping"]["deletion"] == "stratified"
    assert manifest["chromosome_holdout"]["deletion"]["status"] == "skipped"
    assert manifest["thresholds"]["deletion"] > 0
    bound = featurize_row(_row(inslen=">316", rearrangement="A_B"), None)
    assert bound["inslen_lower_bound"] == 1.0
    assert bound["ins_total"] == 0.0
    assert passes_filter(0, 10, 10, length_is_lower_bound=True)[0]
    assert not passes_filter(0, 10, 10)[0]
    scored = predict_rows([_row(score_vs_ref="40"), _row(score_vs_ref="1")], manifest, None)
    assert scored[0]["ml_model"] == "logistic"
    assert float(scored[0]["ml_prob"]) > float(scored[1]["ml_prob"])


def test_vcf_contains_cpx_and_breakends():
    header = vcf_header("HG1", [("chr1", 100000)], model="model_v1")
    assert "##fileformat=VCFv4.2" in header
    assert "##arcsv2_ml_model=model_v1" in header
    assert "ID=CIPOS" in header
    assert "ID=SVSCORE" in header
    lines = format_record(_row(svtype="BND,BND", complextype="ddup.del"), "evt1")
    assert "<CPX>" in lines[0]
    assert "CPX_TYPE=ddup.del" in lines[0]
    assert "REARR=" in lines[0]
    assert "GT:SR:PE:SVSCORE:SVSCORE_NEXT" in lines[0]
    assert any("SVTYPE=BND" in line and "EVENT=evt1" in line for line in lines)


def test_vcf_simple_del_has_ci_and_format(tmp_path):
    row = _row(bp_point_source="mid", bp_uncertainty="10,0,0,12")
    lines = format_record(row, "del1")
    assert "SVTYPE=DEL" in lines[0]
    assert "CIPOS=" in lines[0]
    assert "IMPRECISE" in lines[0]
    assert "0/1:4:8:20:5" in lines[0]
    assert any("SVTYPE=BND" in line for line in lines[1:])
    gz = write_vcf_gz([row], tmp_path / "arcsv_out.vcf.gz", "HG1", [("chr1", 100000)])
    assert gz.is_file()
    tbi = Path(str(gz) + ".tbi")
    assert tbi.is_file()
    import pysam

    with pysam.VariantFile(str(gz)) as vf:
        records = list(vf.fetch())
    types = [r.info.get("SVTYPE") for r in records]
    assert "DEL" in types
    assert "BND" in types
