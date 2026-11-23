"""ml-train speedup keeps the manifest the current trainer wrote."""

from __future__ import annotations

import csv
import importlib.util
from pathlib import Path

import numpy as np

from arcsv2.features import (
    feature_names,
    featurize_row,
    fit_pca,
    fit_pca_arrays,
    rearrangement_matrix,
    reference_indicator,
)
from arcsv2.train import (
    _class_cache,
    _design_indexed,
    _logistic_prob,
    _logistic_probs,
    fit_labeled_calls,
    write_manifest,
)

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "train_perf"


def _rows() -> list[dict[str, str]]:
    with (FIXTURE / "rows.tsv").open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _compare():
    path = ROOT / "scripts" / "compare_manifest.py"
    spec = importlib.util.spec_from_file_location("compare_manifest", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_fit_pca_arrays_matches_fit_pca():
    rows = _rows()
    stacked_r = np.stack(
        [
            rearrangement_matrix(r.get("reference", ""), r.get("rearrangement", "")).ravel()
            for r in rows
        ]
    )
    stacked_f = np.stack([reference_indicator(r.get("reference", "")) for r in rows])
    wrapped = fit_pca(rows).to_dict()
    direct = fit_pca_arrays(stacked_r, stacked_f).to_dict()
    assert wrapped.keys() == direct.keys()
    for key in wrapped:
        np.testing.assert_allclose(wrapped[key], direct[key], atol=1e-12, rtol=1e-12)


def test_cached_pca_columns_match_featurize_row():
    rows = _rows()[:40]
    pca = fit_pca(rows)
    cache = _class_cache(rows)
    matrix, names = _design_indexed(cache, np.arange(len(rows)), pca)
    assert names == feature_names()
    for i, row in enumerate(rows):
        feats = featurize_row(row, pca)
        for j, name in enumerate(names):
            got = matrix[i, j]
            expected = float(feats.get(name, 0.0))
            assert abs(got - expected) <= 1e-9


def test_vectorized_scores_match_logistic_prob_including_clamp():
    spec = {
        "coef": [1.0, -2.0, 0.5],
        "intercept": 0.0,
        "mean": [0.0, 1.0, 0.0],
        "scale": [2.0, 0.0, 1.0],
    }
    values = np.array(
        [
            [0.0, 1.0, 0.0],
            [4.0, 3.0, -2.0],
            [200.0, 1.0, 0.0],
            [0.0, 81.0, 0.0],
            [60.0, 1.0, 0.0],
            [120.0, 1.0, 0.0],
        ],
        dtype=float,
    )
    vector = _logistic_probs(spec, values)
    for i, row in enumerate(values):
        scalar = _logistic_prob(spec, row)
        assert abs(float(vector[i]) - scalar) <= 1e-12
    # z = 200/2 = 100 > 60 clamps to 1. z = -2*(80)/1 = -160 < -60 clamps to 0.
    # z = 60/2 = 30 stays on the sigmoid (strictly below the clamp, and not 1).
    # z = 120/2 = 60 is the clamp boundary and stays on the sigmoid too.
    assert float(vector[2]) == 1.0
    assert float(vector[3]) == 0.0
    assert float(vector[4]) != 1.0
    assert float(vector[4]) == _logistic_prob(spec, values[4])
    assert float(vector[5]) == _logistic_prob(spec, values[5])


def test_fit_labeled_calls_matches_saved_manifest(tmp_path):
    import json

    rows = _rows()
    samples = {row.get("sample") for row in rows}
    chroms = {row.get("chrom") for row in rows if row.get("chrom")}
    assert len(rows) >= 200
    assert len(samples) >= 3
    assert len(chroms) >= 3
    fresh = tmp_path / "manifest.json"
    write_manifest(fresh, fit_labeled_calls(rows))
    errors, stats = _compare().compare(
        json.loads((FIXTURE / "manifest.json").read_text(encoding="utf-8")),
        json.loads(fresh.read_text(encoding="utf-8")),
    )
    assert not errors, errors[:8]
    assert stats["max_abs"] <= 1e-9
    assert stats["max_rel"] <= 1e-9
