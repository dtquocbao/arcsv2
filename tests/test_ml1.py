"""ML-1 models stay off the default trainer and round-trip without pickle."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from arcsv2.features import assert_features_allowed
from arcsv2.ml1 import _outer_chrom_folds, evaluate, score_rows, training_rows
from arcsv2.predict import predict_rows
from arcsv2.train import class_of_row, fit_labeled_calls, load_manifest

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "train_perf"


def _rows() -> list[dict[str, str]]:
    with (FIXTURE / "rows.tsv").open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_labeler_columns_are_not_features():
    with pytest.raises(ValueError, match="labeler columns"):
        assert_features_allowed(["sr_support", "reason", "mapq_missing"])


def test_f2_manifest_refuses_f1_input_and_f1_still_scores():
    manifest = load_manifest(FIXTURE / "manifest.json")
    row = _rows()[0]
    scored = predict_rows([row], manifest, None)
    assert scored[0]["ml_model"] == "logistic"
    assert scored[0]["ml_prob"] != "NA"
    f2 = dict(manifest)
    f2["feature_schema"] = "arcsv2-f2"
    with pytest.raises(ValueError, match="arcsv2-f2"):
        predict_rows([row], f2, None)


def test_xgboost_round_trip_on_fixture():
    rows = [row for row in _rows() if class_of_row(row) == "inversion"]
    labels = {row["label"] for row in rows}
    assert labels == {"0", "1"}
    manifest = fit_labeled_calls(rows, features="f1", model="xgboost")
    assert manifest["feature_schema"] == "arcsv2-f1"
    assert manifest["classifiers"]["inversion"]["format"] == "xgboost-json"
    first = score_rows(rows, manifest)
    second = score_rows(rows, json.loads(json.dumps(manifest)))
    assert first == second
    assert any(value is not None for value in first)


def test_ml1_eval_is_deterministic_on_the_fixture():
    rows = _rows()
    first = evaluate(rows, n_boot=20)
    second = evaluate(rows, n_boot=20)
    assert json.dumps(first, sort_keys=True) == json.dumps(second, sort_keys=True)
    assert "sr>0" in first["classes"]["deletion"]["models"]
    assert "xgboost f1+context no-pca" in first["classes"]["complex"]["models"]
    assert first["protocol"]["outer_group"] == "sample"
    assert first["protocol"]["inner_group"] == "sample"


def test_chrom_folds_hold_out_whole_chromosomes():
    rows = training_rows(_rows())
    folds = _outer_chrom_folds(rows, 5)
    assert len(folds) == 5
    seen: list[str] = []
    by_chrom: dict[str, int] = {}
    for row in rows:
        by_chrom[row["chrom"]] = by_chrom.get(row["chrom"], 0) + 1
    sizes = []
    for train, test in folds:
        assert set(train).isdisjoint(test)
        seen.extend(test)
        sizes.append(sum(by_chrom[chrom] for chrom in test))
    assert sorted(seen) == sorted(by_chrom)
    assert max(sizes) <= 3 * min(sizes)
