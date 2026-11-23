"""Per-class model training (§7.3).

Selection, on out-of-fold scores, picks the threshold that keeps recall ≥ 0.90
and then the candidate with the highest precision there. Tie-break is average
precision. Models are stored as XGBoost JSON or skops, never pickle.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from arcsv2.config import FEATURE_SCHEMA

CLASSES = ("deletion", "duplication", "inversion", "complex")

# Roadmap §4.11. min_samples_split is not an XGBoost parameter and is omitted.
PAPER_PRESETS: dict[str, dict] = {
    "complex": {"kind": "xgboost", "n_estimators": 160, "max_depth": 143, "gamma": 0, "seed": 559},
    "inversion": {
        "kind": "rf",
        "bootstrap": True,
        "max_depth": 215,
        "min_samples_leaf": 1,
        "min_samples_split": 2,
        "max_features": "sqrt",
    },
    "duplication": {"kind": "xgboost", "learning_rate": 0.3, "max_depth": 14, "n_estimators": 195},
    "deletion": {"kind": "xgboost", "learning_rate": 0.3, "max_depth": 84, "n_estimators": 360},
}


@dataclass
class ThresholdChoice:
    threshold: float
    precision: float
    recall: float
    average_precision: float


def _average_precision(y: list[int], scores: list[float]) -> float:
    order = sorted(range(len(scores)), key=lambda i: -scores[i])
    hit = 0
    total_pos = sum(y)
    if total_pos == 0:
        return 0.0
    acc = 0.0
    for rank, idx in enumerate(order, start=1):
        if y[idx] == 1:
            hit += 1
            acc += hit / rank
    return acc / total_pos


def choose_threshold(y: list[int], scores: list[float | None], min_recall: float = 0.90) -> ThresholdChoice:
    """Highest precision among thresholds that keep recall ≥ ``min_recall``.

    Scores of ``None`` were not produced (a one-class training fold) and are
    left out of the threshold. Candidates are the distinct scores, highest
    first, which is the same walk as sorting the unique values descending.
    """
    labels: list[int] = []
    values: list[float] = []
    for label, score in zip(y, scores):
        if score is None:
            continue
        labels.append(label)
        values.append(float(score))
    if not labels:
        return ThresholdChoice(1.0, 0.0, 0.0, 0.0)
    positives = sum(labels)
    average_precision = _average_precision(labels, values)
    # Stable sort so a tied score keeps the first original representative,
    # matching ``sorted(set(values), reverse=True)``.
    order = sorted(range(len(values)), key=lambda i: -values[i])
    best: ThresholdChoice | None = None
    tp = fp = 0
    cursor = 0
    n = len(order)
    while cursor < n:
        threshold = values[order[cursor]]
        while cursor < n and values[order[cursor]] == threshold:
            if labels[order[cursor]] == 1:
                tp += 1
            else:
                fp += 1
            cursor += 1
        recall = tp / positives if positives else 0.0
        precision = tp / (tp + fp) if (tp + fp) else 0.0
        if recall + 1e-12 < min_recall:
            continue
        if best is None or precision > best.precision + 1e-12:
            best = ThresholdChoice(threshold, precision, recall, average_precision)
        elif abs(precision - best.precision) <= 1e-12 and average_precision > best.average_precision:
            best = ThresholdChoice(threshold, precision, recall, average_precision)
    if best is None:
        # No threshold holds the recall. Report the most permissive one.
        return ThresholdChoice(
            min(values) if values else 0.0,
            0.0,
            1.0 if positives else 0.0,
            average_precision,
        )
    return best


def class_of_row(row: dict[str, str]) -> str | None:
    sv = (row.get("svtype") or "").split(",")[0]
    if sv == "DEL":
        return "deletion"
    if sv == "DUP":
        return "duplication"
    if sv == "INV":
        return "inversion"
    if "_" in (row.get("rearrangement") or "") or sv == "INS":
        return None  # insertions are not trained
    return "complex"


def write_manifest(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def load_manifest(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _design(rows: list[dict[str, str]], pca) -> tuple:
    """Feature matrix in ``feature_names()`` order."""
    from arcsv2.features import feature_names, featurize_row

    names = feature_names()
    matrix = []
    for row in rows:
        feats = featurize_row(row, pca)
        matrix.append([float(feats.get(name, 0.0)) for name in names])
    return np.asarray(matrix, dtype=float), names


def _scale_fit(x):
    """Column mean and standard deviation. A constant column is left unscaled."""
    from sklearn.preprocessing import StandardScaler

    scaler = StandardScaler()
    scaler.fit(x)
    scale = np.asarray(scaler.scale_, dtype=float).copy()
    scale[scale < 1e-12] = 1.0
    return np.asarray(scaler.mean_, dtype=float), scale


def _scale_apply(x, mean, scale):
    return (x - mean) / scale


def _logistic_spec(x, y, names: list[str]) -> dict:
    from sklearn.linear_model import LogisticRegression

    mean, scale = _scale_fit(x)
    model = LogisticRegression(max_iter=500)
    model.fit(_scale_apply(x, mean, scale), y)
    return {
        "coef": [float(v) for v in model.coef_.ravel()],
        "intercept": float(model.intercept_.ravel()[0]),
        "features": names,
        "mean": [float(v) for v in mean],
        "scale": [float(v) for v in scale],
    }


def _sample_groups(rows: list[dict[str, str]]) -> list[str] | None:
    samples = [row.get("sample") or "" for row in rows]
    if any(not sample for sample in samples):
        return None
    if len(set(samples)) < 2:
        return None
    return samples


def _class_cache(rows: list[dict[str, str]]):
    """Base features and PCA inputs, computed once per class.

    ``X_base`` holds every ``feature_names()`` column except the four PCA
    columns, in that same order. ``R`` and ``F`` are the stacked rearrangement
    and reference matrices those PCA columns are projected from.
    """
    from arcsv2.features import (
        PCA_FEATURES,
        feature_names,
        featurize_row,
        rearrangement_matrix,
        reference_indicator,
    )

    names = feature_names()
    pca_names = set(PCA_FEATURES)
    base_names = [name for name in names if name not in pca_names]
    base_cols = [names.index(name) for name in base_names]
    pca_cols = [names.index(name) for name in PCA_FEATURES]
    n = len(rows)
    x_base = np.empty((n, len(base_names)), dtype=float)
    rearr_rows = []
    ref_rows = []
    for i, row in enumerate(rows):
        feats = featurize_row(row, None)
        x_base[i] = [float(feats.get(name, 0.0)) for name in base_names]
        rearr_rows.append(
            rearrangement_matrix(row.get("reference", ""), row.get("rearrangement", "")).ravel()
        )
        ref_rows.append(reference_indicator(row.get("reference", "")))
    rearr = np.stack(rearr_rows) if rearr_rows else np.zeros((0, 0))
    refs = np.stack(ref_rows) if ref_rows else np.zeros((0, 0))
    return x_base, rearr, refs, names, base_cols, pca_cols


def _design_indexed(cache, idx, pca) -> tuple:
    """Design matrix for ``idx`` from a class cache and a fitted PCA."""
    x_base, rearr, refs, names, base_cols, pca_cols = cache
    chosen = np.asarray(idx, dtype=int)
    matrix = np.zeros((len(chosen), len(names)), dtype=float)
    matrix[:, base_cols] = x_base[chosen]
    # Same products as ``(R - mean) @ components.T``, reduced one row at a time
    # the way ``PcaModel.transform_*`` does (``components @ centered_row``).
    # A single gemm moves the centered PCA means by ~1e-16 and that is enough
    # to change a stored threshold.
    centered_rearr = rearr[chosen] - pca.rearr_mean
    centered_ref = refs[chosen] - pca.ref_mean
    rearr_scores = _project_rows(centered_rearr, pca.rearr_components)
    ref_scores = _project_rows(centered_ref, pca.ref_components)
    matrix[:, pca_cols[0]] = rearr_scores[:, 0]
    matrix[:, pca_cols[1]] = rearr_scores[:, 1]
    matrix[:, pca_cols[2]] = ref_scores[:, 0]
    matrix[:, pca_cols[3]] = ref_scores[:, 1]
    return matrix, names


def _project_rows(centered: np.ndarray, components: np.ndarray) -> np.ndarray:
    """``components @ centered_row`` for each row. Shape ``(n, n_components)``."""
    scores = np.empty((centered.shape[0], components.shape[0]), dtype=float)
    for i in range(centered.shape[0]):
        scores[i] = components @ centered[i]
    return scores


def _logistic_probs(spec: dict, values) -> np.ndarray:
    """Vector form of ``_logistic_prob``. The scalar function is unchanged."""
    mean = np.asarray(spec["mean"], dtype=float)
    scale = np.asarray(spec["scale"], dtype=float)
    scale_safe = scale.copy()
    scale_safe[scale_safe == 0.0] = 1.0
    coef = np.asarray(spec["coef"], dtype=float)
    scaled = (values - mean) / scale_safe
    # Multiply-then-add in feature order, matching ``_logistic_prob``.
    z = np.full(values.shape[0], float(spec["intercept"]), dtype=float)
    for i in range(coef.shape[0]):
        z += float(coef[i]) * scaled[:, i]
    out = np.empty(z.shape, dtype=float)
    high = z > 60.0
    low = z < -60.0
    mid = ~high & ~low
    out[high] = 1.0
    out[low] = 0.0
    out[mid] = 1.0 / (1.0 + np.exp(-z[mid]))
    return out


def _oof_probabilities(
    rows: list[dict[str, str]],
    y: list[int],
    folds: int = 4,
    cache=None,
) -> tuple[list[float | None], str]:
    """Out-of-fold probabilities. PCA and scaling are fit on each training fold.

    Folds are ``GroupKFold`` on sample id when every row has one and at least
    two samples are present. Otherwise the split is stratified by label.
    A training fold with one class leaves its test scores missing.
    """
    from sklearn.model_selection import GroupKFold, StratifiedKFold

    from arcsv2.features import fit_pca_arrays

    n = len(rows)
    scores: list[float | None] = [None] * n
    groups = _sample_groups(rows)
    n_pos = sum(y)
    n_neg = n - n_pos
    if groups is not None:
        n_splits = min(folds, len(set(groups)))
        grouping = "sample"
    else:
        n_splits = min(folds, n, n_pos, n_neg)
        grouping = "stratified"
    if cache is None:
        cache = _class_cache(rows)
    _x_base, rearr, refs, names, _base_cols, _pca_cols = cache
    if n_splits < 2 or len(set(y)) < 2:
        pca = fit_pca_arrays(rearr, refs)
        x, names = _design_indexed(cache, np.arange(n), pca)
        spec = _logistic_spec(x, y, names)
        return [float(v) for v in _logistic_probs(spec, x)], grouping
    indices = np.arange(n)
    if grouping == "sample":
        splitter = GroupKFold(n_splits=n_splits)
        splits = splitter.split(indices, y, groups=groups)
    else:
        splitter = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=204178949)
        splits = splitter.split(indices, y)
    for train_idx, test_idx in splits:
        train_y = [y[i] for i in train_idx]
        if len(set(train_y)) < 2:
            continue
        pca = fit_pca_arrays(rearr[train_idx], refs[train_idx])
        x_train, names = _design_indexed(cache, train_idx, pca)
        spec = _logistic_spec(x_train, train_y, names)
        x_test, _names = _design_indexed(cache, test_idx, pca)
        folded = _logistic_probs(spec, x_test)
        for local, idx in enumerate(test_idx):
            scores[int(idx)] = float(folded[local])
    return scores, grouping


def _logistic_prob(spec: dict, values) -> float:
    import math

    mean = spec.get("mean")
    scale = spec.get("scale")
    z = float(spec["intercept"])
    for i, (coef, value) in enumerate(zip(spec["coef"], values)):
        number = float(value)
        if mean is not None and scale is not None:
            denom = float(scale[i]) or 1.0
            number = (number - float(mean[i])) / denom
        z += float(coef) * number
    if z > 60.0:
        return 1.0
    if z < -60.0:
        return 0.0
    return 1.0 / (1.0 + math.exp(-z))


def _chromosome_holdout(
    rows: list[dict[str, str]],
    y: list[int],
    threshold: float,
    cache=None,
) -> dict:
    """Train on every chromosome but one and score the held-out chromosome.

    Rows with an empty chromosome stay in every training set. The chromosome
    index map is built once.
    """
    from arcsv2.features import fit_pca_arrays

    by_chrom: dict[str, list[int]] = {}
    for i, row in enumerate(rows):
        chrom = row.get("chrom") or ""
        if chrom:
            by_chrom.setdefault(chrom, []).append(i)
    names = sorted(by_chrom)
    if len(names) < 2:
        return {"status": "skipped", "reason": "fewer than 2 chromosomes"}
    if cache is None:
        cache = _class_cache(rows)
    _x_base, rearr, refs, _names, _base_cols, _pca_cols = cache
    n = len(rows)
    per: dict[str, dict] = {}
    for held in names:
        test_idx = by_chrom[held]
        train_mask = np.ones(n, dtype=bool)
        train_mask[test_idx] = False
        train_idx = np.flatnonzero(train_mask)
        train_y = [y[i] for i in train_idx]
        test_y = [y[i] for i in test_idx]
        if len(set(train_y)) < 2 or not test_idx:
            per[held] = {"status": "skipped"}
            continue
        pca = fit_pca_arrays(rearr[train_idx], refs[train_idx])
        x_train, names_ = _design_indexed(cache, train_idx, pca)
        spec = _logistic_spec(x_train, train_y, names_)
        x_test, _feature_names = _design_indexed(cache, test_idx, pca)
        scores = [float(v) for v in _logistic_probs(spec, x_test)]
        choice = choose_threshold(test_y, scores, min_recall=0.0)
        tp = fp = 0
        for label, score in zip(test_y, scores):
            if score >= threshold:
                if label == 1:
                    tp += 1
                else:
                    fp += 1
        positives = sum(test_y)
        per[held] = {
            "status": "ok",
            "n": len(test_idx),
            "recall_at_oof_threshold": (tp / positives) if positives else 0.0,
            "precision_at_oof_threshold": tp / (tp + fp) if (tp + fp) else 0.0,
            "average_precision": choice.average_precision,
        }
    return {"status": "ok", "chromosomes": per}


def fit_labeled_calls(rows: list[dict[str, str]]) -> dict:
    """Per-class logistic models, out-of-fold thresholds, and the training PCA.

    ``rows`` are call-table records with a ``label`` of 0 or 1. Insertions are
    left out. The coefficients are stored in the manifest; nothing is pickled.
    """
    from arcsv2.features import fit_pca_arrays

    grouped: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        kind = class_of_row(row)
        if kind is None:
            continue
        label = row.get("label", "")
        if label not in ("0", "1"):
            continue
        grouped.setdefault(kind, []).append(row)
    classifiers: dict[str, dict] = {}
    thresholds: dict[str, float] = {}
    grouping: dict[str, str] = {}
    holdout: dict[str, dict] = {}
    pca_by_class: dict[str, dict] = {}
    for kind, group in grouped.items():
        y = [int(row["label"]) for row in group]
        if len(set(y)) < 2:
            continue
        # One featurization. Every fold indexes into this cache.
        cache = _class_cache(group)
        _x_base, rearr, refs, _names, _base_cols, _pca_cols = cache
        class_pca = fit_pca_arrays(rearr, refs)
        oof, how = _oof_probabilities(group, y, cache=cache)
        thresholds[kind] = choose_threshold(y, oof).threshold
        x, names = _design_indexed(cache, np.arange(len(group)), class_pca)
        spec = _logistic_spec(x, y, names)
        spec["pca"] = class_pca.to_dict()
        classifiers[kind] = spec
        grouping[kind] = how
        holdout[kind] = _chromosome_holdout(group, y, thresholds[kind], cache=cache)
        pca_by_class[kind] = class_pca.to_dict()
    return {
        "feature_schema": FEATURE_SCHEMA,
        "fitted": bool(classifiers),
        "models": {kind: "logistic" for kind in classifiers},
        "pca": pca_by_class,
        "classifiers": classifiers,
        "thresholds": thresholds,
        "grouping": grouping,
        "chromosome_holdout": holdout,
    }
