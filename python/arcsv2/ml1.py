"""Per-class models beyond the default logistic f1 trainer.

The default ``fit_labeled_calls`` path is unchanged. This module is used only
when ``ml-train`` is given ``--features`` or ``--model``, and by ``ml-eval``.
Every fit (PCA, scaling, imputation, model, hyperparameters, threshold) uses
training rows only.
"""

from __future__ import annotations

import base64
import hashlib
import math
import sys
import tempfile
from pathlib import Path

import numpy as np

from arcsv2.config import FEATURE_SCHEMA
from arcsv2.features import (
    CONTEXT_FEATURES,
    FEATURE_SCHEMA_F2,
    LABELER_COLUMNS,
    assert_features_allowed,
    featurize_row,
)
from arcsv2.train import (
    CLASSES,
    PAPER_PRESETS,
    _class_cache,
    _design_indexed,
    _logistic_probs,
    _logistic_spec,
    choose_threshold,
    class_of_row,
)

# Same seed the default trainer uses for StratifiedKFold.
SEED = 204178949
N_JOBS = 1
INNER_SPLITS = 4
OUTER_SPLITS = 5
N_BOOT = 1000


def _seed() -> int:
    return SEED


def library_versions() -> dict[str, str]:
    import sklearn
    import skops
    import xgboost

    return {
        "scikit-learn": sklearn.__version__,
        "xgboost": xgboost.__version__,
        "skops": skops.__version__,
    }


def parse_context_value(raw: str | None, name: str) -> float:
    text = "" if raw is None else str(raw).strip()
    if text in ("", "NA", "na", "NaN", "nan", "None", "."):
        return float("nan")
    if name == "length_is_lower_bound":
        lowered = text.lower()
        if lowered in ("1", "true", "t", "yes"):
            return 1.0
        if lowered in ("0", "false", "f", "no"):
            return 0.0
    try:
        return float(text)
    except ValueError:
        return float("nan")


def context_matrix(rows: list[dict[str, str]]) -> np.ndarray:
    matrix = np.empty((len(rows), len(CONTEXT_FEATURES)), dtype=float)
    for i, row in enumerate(rows):
        for j, name in enumerate(CONTEXT_FEATURES):
            matrix[i, j] = parse_context_value(row.get(name), name)
    assert_features_allowed(list(CONTEXT_FEATURES))
    return matrix


def _medians(block: np.ndarray) -> np.ndarray:
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        med = np.nanmedian(block, axis=0)
    med = np.asarray(med, dtype=float)
    med[np.isnan(med)] = 0.0
    return med


def _impute(block: np.ndarray, medians: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    missing = np.isnan(block)
    filled = block.copy()
    for j in range(filled.shape[1]):
        filled[missing[:, j], j] = medians[j]
    return filled, missing.astype(float)


_PCA_NAMES = ("rearr_pc1", "rearr_pc2", "ref_pc1", "ref_pc2")


def design_matrix(
    cache,
    context: np.ndarray | None,
    idx,
    pca,
    feature_set: str,
    family: str,
    medians=None,
    drop_pca: bool = False,
):
    """Fold design. Imputation statistics come from ``idx`` unless ``medians`` is passed."""
    matrix, names = _design_indexed(cache, idx, pca)
    names = list(names)
    used = None
    if drop_pca:
        keep = [i for i, name in enumerate(names) if name not in _PCA_NAMES]
        matrix = matrix[:, keep]
        names = [names[i] for i in keep]
    if feature_set == "f1+context":
        if context is None:
            raise ValueError("f1+context requires the context columns")
        chosen = np.asarray(idx, dtype=int)
        block = context[chosen]
        if family == "xgboost":
            extra = list(CONTEXT_FEATURES)
            matrix = np.concatenate([matrix, block], axis=1)
        else:
            used = _medians(block) if medians is None else np.asarray(medians, dtype=float)
            filled, missing = _impute(block, used)
            extra = list(CONTEXT_FEATURES) + [f"{name}_missing" for name in CONTEXT_FEATURES]
            matrix = np.concatenate([matrix, filled, missing], axis=1)
        names.extend(extra)
    assert_features_allowed(names)
    if LABELER_COLUMNS.intersection(names):
        raise ValueError(f"labeler columns cannot be features: {sorted(LABELER_COLUMNS.intersection(names))}")
    return matrix, names, used


def sr_score(row: dict[str, str]) -> float:
    """1 when any split-read support count is > 0."""
    text = row.get("sr_support") or ""
    for part in str(text).split(","):
        piece = part.strip()
        if not piece or piece in ("NA", "na", "."):
            continue
        try:
            if float(piece) > 0.0:
                return 1.0
        except ValueError:
            continue
    return 0.0


def _xgb_candidates(kind: str) -> list[dict]:
    grid = []
    for depth in (3, 6, 10):
        for trees in (200, 500):
            for rate in (0.05, 0.1):
                for child in (1, 5):
                    grid.append(
                        {
                            "family": "xgboost",
                            "max_depth": depth,
                            "n_estimators": trees,
                            "learning_rate": rate,
                            "min_child_weight": child,
                            "subsample": 0.8,
                            "colsample_bytree": 0.8,
                            "gamma": 0.0,
                            "seed": _seed(),
                            "source": "grid",
                        }
                    )
    preset = PAPER_PRESETS.get(kind, {})
    if preset.get("kind") == "xgboost":
        grid.append(
            {
                "family": "xgboost",
                "max_depth": int(preset["max_depth"]),
                "n_estimators": int(preset["n_estimators"]),
                "learning_rate": float(preset.get("learning_rate", 0.3)),
                "min_child_weight": 1,
                "subsample": 0.8,
                "colsample_bytree": 0.8,
                "gamma": float(preset.get("gamma", 0.0)),
                "seed": int(preset.get("seed", _seed())),
                "source": "paper",
            }
        )
    return grid


def _rf_candidates(kind: str) -> list[dict]:
    grid = [
        {
            "family": "rf",
            "n_estimators": 300,
            "max_depth": None,
            "min_samples_leaf": 1,
            "min_samples_split": 2,
            "max_features": "sqrt",
            "bootstrap": True,
            "seed": _seed(),
            "source": "grid",
        }
    ]
    preset = PAPER_PRESETS.get(kind, {})
    if preset.get("kind") == "rf":
        grid.append(
            {
                "family": "rf",
                "n_estimators": int(preset.get("n_estimators", 100)),
                "max_depth": preset.get("max_depth"),
                "min_samples_leaf": int(preset.get("min_samples_leaf", 1)),
                "min_samples_split": int(preset.get("min_samples_split", 2)),
                "max_features": preset.get("max_features", "sqrt"),
                "bootstrap": bool(preset.get("bootstrap", True)),
                "seed": _seed(),
                "source": "paper",
            }
        )
    return grid


def _fit_xgb(x, y: np.ndarray, params: dict):
    import xgboost as xgb

    n_pos = int(np.sum(y == 1))
    n_neg = int(np.sum(y == 0))
    model = xgb.XGBClassifier(
        objective="binary:logistic",
        max_depth=int(params["max_depth"]),
        n_estimators=int(params["n_estimators"]),
        learning_rate=float(params["learning_rate"]),
        min_child_weight=float(params["min_child_weight"]),
        subsample=float(params.get("subsample", 0.8)),
        colsample_bytree=float(params.get("colsample_bytree", 0.8)),
        gamma=float(params.get("gamma", 0.0)),
        scale_pos_weight=(n_neg / n_pos) if n_pos else 1.0,
        tree_method="hist",
        n_jobs=N_JOBS,
        random_state=int(params.get("seed", _seed())),
        verbosity=0,
    )
    model.fit(np.asarray(x, dtype=float), np.asarray(y, dtype=int))
    return model


def _fit_rf(x, y: np.ndarray, params: dict):
    from sklearn.ensemble import RandomForestClassifier

    model = RandomForestClassifier(
        n_estimators=int(params.get("n_estimators", 100)),
        max_depth=params.get("max_depth"),
        min_samples_leaf=int(params.get("min_samples_leaf", 1)),
        min_samples_split=int(params.get("min_samples_split", 2)),
        max_features=params.get("max_features", "sqrt"),
        bootstrap=bool(params.get("bootstrap", True)),
        random_state=int(params.get("seed", _seed())),
        n_jobs=N_JOBS,
    )
    model.fit(np.asarray(x, dtype=float), np.asarray(y, dtype=int))
    return model


def _predict_fitted(kind: str, model, x) -> np.ndarray:
    if kind == "logistic":
        return np.asarray(_logistic_probs(model, x), dtype=float)
    return np.asarray(model.predict_proba(np.asarray(x, dtype=float))[:, 1], dtype=float)


def _pca_key(train_idx) -> bytes:
    return hashlib.sha1(np.asarray(train_idx, dtype=np.int32).tobytes()).digest()


def _fit_fold(
    family: str,
    params: dict,
    cache,
    context,
    train_idx,
    y_train,
    feature_set: str,
    drop_pca: bool = False,
    pca_cache: dict | None = None,
):
    from arcsv2.features import fit_pca_arrays

    _x_base, rearr, refs, _names, _base_cols, _pca_cols = cache
    key = _pca_key(train_idx) if pca_cache is not None else None
    pca = None if key is None else pca_cache.get(key)
    if pca is None:
        pca = fit_pca_arrays(rearr[train_idx], refs[train_idx])
        if key is not None:
            pca_cache[key] = pca
    matrix, names, medians = design_matrix(
        cache, context, train_idx, pca, feature_set, family, drop_pca=drop_pca
    )
    y = np.asarray(y_train, dtype=int)
    if family == "logistic":
        spec = _logistic_spec(matrix, y.tolist(), names)
        if medians is not None:
            spec["medians"] = [float(v) for v in medians]
        spec["pca"] = pca.to_dict()
        spec["context"] = "median" if feature_set == "f1+context" else "none"
        return "logistic", spec, pca, medians, names
    if family == "xgboost":
        model = _fit_xgb(matrix, y, params)
        return "xgboost", model, pca, medians, names
    model = _fit_rf(matrix, y, params)
    return "rf", model, pca, medians, names


def _score_fold(kind, model, cache, context, idx, pca, feature_set: str, medians, names, drop_pca: bool = False):
    matrix, _names, _med = design_matrix(
        cache,
        context,
        idx,
        pca,
        feature_set,
        kind,
        medians,
        drop_pca=drop_pca,
    )
    # Names are determined by the feature set; the stored list is the training order.
    del names
    return _predict_fitted(kind, model, matrix)


def _average_precision(y, scores) -> float | None:
    from sklearn.metrics import average_precision_score

    labels = np.asarray(y, dtype=int)
    values = np.asarray(scores, dtype=float)
    keep = np.isfinite(values)
    labels, values = labels[keep], values[keep]
    if len(labels) == 0 or len(set(labels.tolist())) < 2:
        return None
    return float(average_precision_score(labels, values))


def _roc_auc(y, scores) -> float | None:
    from sklearn.metrics import roc_auc_score

    labels = np.asarray(y, dtype=int)
    values = np.asarray(scores, dtype=float)
    keep = np.isfinite(values)
    labels, values = labels[keep], values[keep]
    if len(labels) == 0 or len(set(labels.tolist())) < 2:
        return None
    return float(roc_auc_score(labels, values))


def _at_threshold(y, scores, threshold: float) -> tuple[float, float]:
    tp = fp = positives = 0
    for label, score in zip(y, scores):
        if score is None or (isinstance(score, float) and not math.isfinite(score)):
            continue
        if int(label) == 1:
            positives += 1
        if float(score) >= threshold:
            if int(label) == 1:
                tp += 1
            else:
                fp += 1
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / positives if positives else 0.0
    return precision, recall


def choose_precision_threshold(y, scores, min_precision: float = 0.95):
    """Highest recall among thresholds with precision >= ``min_precision``.

    The scan is highest score first, so an equal recall keeps the higher threshold.
    Returns ``None`` when no threshold reaches the precision.
    """
    paired = [
        (int(label), float(score))
        for label, score in zip(y, scores)
        if score is not None and math.isfinite(float(score))
    ]
    if not paired:
        return None
    labels = [label for label, _score in paired]
    values = [score for _label, score in paired]
    positives = sum(labels)
    if positives == 0:
        return None
    order = sorted(range(len(values)), key=lambda i: -values[i])
    best = None
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
        precision = tp / (tp + fp) if (tp + fp) else 0.0
        recall = tp / positives
        if precision + 1e-12 < min_precision:
            continue
        if best is None or recall > best["recall"] + 1e-12:
            best = {"threshold": threshold, "precision": precision, "recall": recall}
    return best


def _inner_oof(
    rows,
    y,
    cache,
    context,
    train_idx,
    feature_set: str,
    family: str,
    params: dict,
    drop_pca: bool = False,
    pca_cache: dict | None = None,
    group_col: str = "sample",
):
    from sklearn.model_selection import GroupKFold

    groups = [rows[i].get(group_col) or "" for i in train_idx]
    n_groups = len(set(groups))
    n_splits = min(INNER_SPLITS, n_groups)
    scores: list[float | None] = [None] * len(train_idx)
    if n_splits < 2:
        return scores, "too few samples"
    y_train = [y[i] for i in train_idx]
    splitter = GroupKFold(n_splits=n_splits)
    dummy = np.zeros((len(train_idx), 1))
    for inner_train, inner_test in splitter.split(dummy, y_train, groups):
        global_train = [train_idx[i] for i in inner_train]
        global_test = [train_idx[i] for i in inner_test]
        fold_y = [y[i] for i in global_train]
        if len(set(fold_y)) < 2:
            continue
        try:
            kind, model, pca, medians, names = _fit_fold(
                family, params, cache, context, global_train, fold_y, feature_set, drop_pca, pca_cache
            )
            pred = _score_fold(
                kind, model, cache, context, global_test, pca, feature_set, medians, names, drop_pca
            )
        except Exception as exc:
            return None, f"{type(exc).__name__}: {exc}"
        for local, value in zip(inner_test, pred):
            scores[int(local)] = float(value)
    return scores, "ok"


def _select_params(
    rows,
    y,
    cache,
    context,
    train_idx,
    feature_set: str,
    family: str,
    candidates: list[dict],
    drop_pca: bool = False,
    pca_cache: dict | None = None,
    group_col: str = "sample",
):
    best = None
    errors = []
    y_local = [y[i] for i in train_idx]
    for params in candidates:
        if params.get("family") != "logistic":
            _log(
                "      candidate {source} depth={depth} trees={trees} rate={rate}".format(
                    source=params.get("source"),
                    depth=params.get("max_depth"),
                    trees=params.get("n_estimators"),
                    rate=params.get("learning_rate"),
                )
            )
        scores, status = _inner_oof(
            rows, y, cache, context, train_idx, feature_set, family, params, drop_pca, pca_cache, group_col
        )
        if scores is None:
            errors.append({"params": _public_params(params), "error": status})
            continue
        ap = _average_precision(y_local, [np.nan if s is None else s for s in scores])
        if ap is None:
            errors.append({"params": _public_params(params), "error": "undefined AP"})
            continue
        if best is None or ap > best["ap"] + 1e-12:
            best = {"ap": ap, "params": params, "scores": scores, "status": status}
    return best, errors


def _public_params(params: dict) -> dict:
    out = {}
    for key, value in params.items():
        if key == "family":
            continue
        if isinstance(value, (np.floating, np.integer)):
            value = value.item()
        out[key] = value
    return out


def _ci(values: list[float | None]) -> list[float | None]:
    clean = [float(v) for v in values if v is not None and math.isfinite(float(v))]
    if not clean:
        return [None, None]
    lo, hi = np.percentile(np.asarray(clean, dtype=float), [2.5, 97.5])
    return [float(lo), float(hi)]


def _at_rows(y, scores, thresholds) -> tuple[float | None, float | None]:
    """Precision and recall using each row's own training-fold threshold."""
    tp = fp = positives = 0
    used = False
    for label, score, threshold in zip(y, scores, thresholds):
        if threshold is None or score is None or not math.isfinite(float(score)):
            continue
        used = True
        if int(label) == 1:
            positives += 1
        if float(score) >= float(threshold):
            if int(label) == 1:
                tp += 1
            else:
                fp += 1
    if not used:
        return None, None
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / positives if positives else 0.0
    return precision, recall


def _bootstrap(y, scores, samples, thr_r, thr_p, n_boot: int, seed: int):
    by: dict[str, list[int]] = {}
    for i, sample in enumerate(samples):
        by.setdefault(sample, []).append(i)
    keys = sorted(by)
    rng = np.random.default_rng(seed)
    aps: list[float | None] = []
    precisions: list[float | None] = []
    recalls: list[float | None] = []
    prec95: list[float | None] = []
    rec95: list[float | None] = []
    for _ in range(n_boot):
        drawn = rng.choice(keys, size=len(keys), replace=True)
        idx: list[int] = []
        for sample in drawn:
            idx.extend(by[sample])
        y_b = [y[i] for i in idx]
        s_b = [scores[i] for i in idx]
        aps.append(_average_precision(y_b, s_b))
        precision, recall = _at_rows(y_b, s_b, [thr_r[i] for i in idx])
        precisions.append(precision)
        recalls.append(recall)
        precision, recall = _at_rows(y_b, s_b, [thr_p[i] for i in idx])
        prec95.append(precision)
        rec95.append(recall)
    return {
        "ap": _ci(aps),
        "precision_at_recall90": _ci(precisions),
        "recall_at_recall90": _ci(recalls),
        "precision_at_precision95": _ci(prec95),
        "recall_at_precision95": _ci(rec95),
        "n_boot": n_boot,
    }


def _gain(model, names: list[str]) -> dict[str, float]:
    booster = model.get_booster()
    raw = booster.get_score(importance_type="gain")
    mapped: dict[str, float] = {}
    for key, value in raw.items():
        if key.startswith("f") and key[1:].isdigit():
            index = int(key[1:])
            name = names[index] if index < len(names) else key
        else:
            name = key
        mapped[name] = float(value)
    return mapped


def _mean_gain(gains: list[dict[str, float]]) -> list[dict[str, float]]:
    totals: dict[str, float] = {}
    for gain in gains:
        for name, value in gain.items():
            totals[name] = totals.get(name, 0.0) + value
    if not gains:
        return []
    ranked = sorted(totals, key=lambda name: (-totals[name] / len(gains), name))
    return [{"feature": name, "gain": totals[name] / len(gains)} for name in ranked[:15]]


def _log(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


def _outer_sample_folds(rows: list[dict[str, str]], n_splits: int) -> list[tuple[list[str], list[str]]]:
    from sklearn.model_selection import StratifiedGroupKFold

    samples: list[str] = []
    seen: set[str] = set()
    superpop: dict[str, str] = {}
    for row in rows:
        sample = row.get("sample") or ""
        if sample not in seen:
            seen.add(sample)
            samples.append(sample)
        superpop[sample] = row.get("superpopulation") or "NA"
    y = [superpop[sample] for sample in samples]
    counts: dict[str, int] = {}
    for label in y:
        counts[label] = counts.get(label, 0) + 1
    usable = min(n_splits, min(counts.values()) if counts else 0)
    if usable < 2:
        raise ValueError(f"need at least 2 samples in every superpopulation, got {counts}")
    splitter = StratifiedGroupKFold(n_splits=usable, shuffle=True, random_state=_seed())
    dummy = np.zeros((len(samples), 1))
    folds = []
    for train_idx, test_idx in splitter.split(dummy, y, samples):
        folds.append(([samples[i] for i in train_idx], [samples[i] for i in test_idx]))
    return folds


def _outer_chrom_folds(rows: list[dict[str, str]], n_splits: int) -> list[tuple[list[str], list[str]]]:
    """GroupKFold over chromosomes. Folds are balanced by row count.

    sklearn assigns the heaviest remaining chromosome to the lightest fold,
    so a large chromosome is not paired with another large one.
    """
    from sklearn.model_selection import GroupKFold

    chroms = [row.get("chrom") or "" for row in rows]
    n_groups = len(set(chroms))
    usable = min(n_splits, n_groups)
    if usable < 2:
        raise ValueError(f"need at least 2 chromosomes, got {n_groups}")
    splitter = GroupKFold(n_splits=usable)
    dummy = np.zeros((len(rows), 1))
    folds = []
    for train_idx, test_idx in splitter.split(dummy, groups=chroms):
        train = sorted({chroms[i] for i in train_idx})
        test = sorted({chroms[i] for i in test_idx})
        folds.append((train, test))
    return folds


def _pack_model(kind: str, model, names: list[str], medians, pca, params: dict) -> dict:
    payload: dict = {
        "format": kind,
        "features": names,
        "params": _public_params(params),
    }
    if pca is not None:
        payload["pca"] = pca.to_dict()
    if kind == "logistic":
        payload.update(model)
        payload["format"] = "logistic"
        return payload
    if kind == "xgboost":
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "model.json"
            model.save_model(path)
            text = path.read_text(encoding="utf-8")
        payload["format"] = "xgboost-json"
        payload["model_json"] = text
        payload["sha256"] = hashlib.sha256(text.encode("utf-8")).hexdigest()
        return payload
    import skops.io as sio

    blob = sio.dumps(model)
    payload["format"] = "skops"
    payload["model_b64"] = base64.b64encode(blob).decode("ascii")
    payload["sha256"] = hashlib.sha256(blob).hexdigest()
    return payload


def _load_xgb(spec: dict):
    import xgboost as xgb

    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "model.json"
        path.write_text(spec["model_json"], encoding="utf-8")
        model = xgb.XGBClassifier()
        model.load_model(path)
    return model


def _load_rf(spec: dict):
    import skops.io as sio

    blob = base64.b64decode(spec["model_b64"])
    untrusted = sio.get_untrusted_types(data=blob)
    return sio.loads(blob, trusted=untrusted)


def load_table(path: Path) -> list[dict[str, str]]:
    """Read a TSV or TSV.GZ call table. The default ``ml-train`` reader is separate."""
    import csv
    import gzip

    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def subset_samples(
    rows: list[dict[str, str]], n: int, superpopulation: str | None = None
) -> tuple[list[dict[str, str]], str, list[str]]:
    """Keep ``n`` sorted sample ids from one superpopulation.

    The superpopulation with the most samples wins a tie by name. A 5-fold
    split needs every stratum to have at least 5 samples, so a mixed subset
    of five samples is not used.
    """
    by: dict[str, set[str]] = {}
    for row in rows:
        by.setdefault(row.get("superpopulation") or "NA", set()).add(row.get("sample") or "")
    if superpopulation is not None:
        chosen_sp = superpopulation
    else:
        chosen_sp = sorted(by, key=lambda name: (-len(by[name]), name))[0]
    chosen = sorted(by.get(chosen_sp, []))[:n]
    keep = set(chosen)
    return [row for row in rows if (row.get("sample") or "") in keep], chosen_sp, chosen


def training_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    kept = []
    for row in rows:
        if row.get("exclude") == "1":
            continue
        if row.get("label") not in ("0", "1"):
            continue
        if class_of_row(row) is None:
            continue
        kept.append(row)
    return kept


def _candidate_list(kind: str, family: str) -> list[dict]:
    if family == "logistic":
        return [{"family": "logistic", "source": "logistic", "seed": _seed()}]
    if family == "xgboost":
        return _xgb_candidates(kind)
    if family == "rf":
        return _rf_candidates(kind)
    raise ValueError(f"unknown model family {family}")


def apply_exclude_splits(
    rows: list[dict[str, str]], splits: list[str]
) -> tuple[list[dict[str, str]], dict]:
    """Drop rows whose ``split`` is in ``splits`` before any fold is built.

    An empty ``splits`` list leaves the table unchanged. Logistic regression
    later fills a missing ``coverage_ratio`` with the training-fold median and
    adds a missing indicator. XGBoost keeps that value as NaN.
    """
    names = [name for name in splits if name]
    if not names:
        return rows, {"splits": [], "rows": 0, "samples": 0, "sample_ids": []}
    if rows and "split" not in rows[0]:
        raise ValueError("dataset has no split column; --exclude-split cannot be applied")
    drop = set(names)
    excluded = [row for row in rows if row.get("split") in drop]
    kept = [row for row in rows if row.get("split") not in drop]
    sample_ids = sorted({row.get("sample") or "" for row in excluded})
    return kept, {
        "splits": list(names),
        "rows": len(excluded),
        "samples": len(sample_ids),
        "sample_ids": sample_ids,
    }


def evaluate(
    rows: list[dict[str, str]],
    *,
    n_splits: int = OUTER_SPLITS,
    n_boot: int = N_BOOT,
    include_rf: str = "inversion",
    outer_group: str = "sample",
    ablation: bool = True,
    exclude_splits: list[str] | None = None,
) -> dict:
    """Outer-fold report. ``include_rf`` is ``inversion`` or ``all`` or ``none``.

    ``outer_group`` is ``sample`` (stratified by superpopulation) or ``chrom``.
    Chromosome folds use ``GroupKFold`` on the outside and the inside. The
    bootstrap is still over samples. ``exclude_splits`` drops those rows before
    folds are built.

    Missing context values, including ``coverage_ratio``, are not imputed in
    the table. Logistic regression fills each one with the median of the
    training fold and appends a 0/1 missing indicator. XGBoost leaves NaN in
    the matrix so the tree can split on missingness.
    """
    if outer_group not in ("sample", "chrom"):
        raise ValueError(f"unknown outer group {outer_group}")
    rows, excluded = apply_exclude_splits(rows, list(exclude_splits or []))
    if excluded["rows"]:
        _log(
            f"excluded split={','.join(excluded['splits'])} "
            f"samples={excluded['samples']} rows={excluded['rows']}"
        )
    rows = training_rows(rows)
    if outer_group == "chrom":
        folds = _outer_chrom_folds(rows, n_splits)
        outer_name = "GroupKFold"
        stratify = None
        group_col = "chrom"
    else:
        folds = _outer_sample_folds(rows, n_splits)
        outer_name = "StratifiedGroupKFold"
        stratify = "superpopulation"
        group_col = "sample"
    _log(
        f"ml-eval rows={len(rows)} samples={len({r.get('sample') for r in rows})} "
        f"folds={len(folds)} outer_group={outer_group}"
    )
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        grouped.setdefault(class_of_row(row) or "", []).append(row)
    report_classes = {}
    for kind in CLASSES:
        group = grouped.get(kind, [])
        y = [int(row["label"]) for row in group]
        if len(group) < 2 or len(set(y)) < 2:
            report_classes[kind] = {"n": len(group), "status": "skipped"}
            continue
        _log(f"class {kind} n={len(group)} pos={sum(y)} neg={len(y) - sum(y)}")
        cache = _class_cache(group)
        context = context_matrix(group)
        by_group: dict[str, list[int]] = {}
        for i, row in enumerate(group):
            by_group.setdefault(row.get(group_col) or "", []).append(i)
        families = [
            ("logistic", "f1", False),
            ("logistic", "f1+context", False),
        ]
        if ablation:
            families.append(("logistic", "f1+context", True))
        families.extend([
            ("xgboost", "f1", False),
            ("xgboost", "f1+context", False),
        ])
        if ablation:
            families.append(("xgboost", "f1+context", True))
        if include_rf == "all" or (include_rf == "inversion" and kind == "inversion"):
            families.extend([("rf", "f1", False), ("rf", "f1+context", False)])
        models = {"sr>0": _rule_report(group, y, folds, by_group, n_boot)}
        importance = {}
        held: dict[str, dict] = {}
        pca_cache: dict = {}
        for family, feature_set, drop_pca in families:
            label = f"{family} {feature_set}" + (" no-pca" if drop_pca else "")
            _log(f"  {label}")
            pooled_y: list[int] = []
            pooled_s: list[float] = []
            pooled_sample: list[str] = []
            pooled_sp: list[str] = []
            chosen_params = []
            paper_errors = []
            gains = []
            pooled_thr_r: list[float | None] = []
            pooled_thr_p: list[float | None] = []
            pooled_i: list[int] = []
            candidates = _candidate_list(kind, family)
            for fold_i, (train_keys, test_keys) in enumerate(folds, start=1):
                train_idx = [i for key in train_keys for i in by_group.get(key, [])]
                test_idx = [i for key in test_keys for i in by_group.get(key, [])]
                if len(train_idx) < 2 or len({y[i] for i in train_idx}) < 2 or not test_idx:
                    continue
                best, errors = _select_params(
                    group,
                    y,
                    cache,
                    context,
                    train_idx,
                    feature_set,
                    family,
                    candidates,
                    drop_pca,
                    pca_cache,
                    group_col,
                )
                paper_errors.extend(errors)
                if best is None:
                    continue
                params = best["params"]
                chosen_params.append(_public_params(params))
                recall_choice = choose_threshold(
                    [y[train_idx[i]] for i in range(len(train_idx))],
                    best["scores"],
                    min_recall=0.90,
                )
                precision_choice = choose_precision_threshold(
                    [y[train_idx[i]] for i in range(len(train_idx))],
                    best["scores"],
                    min_precision=0.95,
                )
                thr_r = recall_choice.threshold
                thr_p = None if precision_choice is None else precision_choice["threshold"]
                train_y = [y[i] for i in train_idx]
                kind_m, model, pca, medians, names = _fit_fold(
                    family, params, cache, context, train_idx, train_y, feature_set, drop_pca, pca_cache
                )
                pred = _score_fold(
                    kind_m, model, cache, context, test_idx, pca, feature_set, medians, names, drop_pca
                )
                for local, value in enumerate(pred):
                    row_i = test_idx[local]
                    pooled_y.append(y[row_i])
                    pooled_s.append(float(value))
                    pooled_sample.append(group[row_i].get("sample") or "")
                    pooled_sp.append(group[row_i].get("superpopulation") or "NA")
                    pooled_thr_r.append(thr_r)
                    pooled_thr_p.append(thr_p)
                    pooled_i.append(row_i)
                if family == "xgboost":
                    gains.append(_gain(model, names))
                _log(f"    fold {fold_i}/{len(folds)} inner AP {best['ap']:.4f} source {params.get('source')}")
            models[label] = _finish_model(
                pooled_y,
                pooled_s,
                pooled_sample,
                pooled_sp,
                pooled_thr_r,
                pooled_thr_p,
                chosen_params,
                paper_errors,
                n_boot,
            )
            if family == "xgboost" and not drop_pca:
                importance[label] = _mean_gain(gains)
            held[label] = {
                "i": pooled_i,
                "y": pooled_y,
                "s": pooled_s,
                "sample": pooled_sample,
                "thr_r": pooled_thr_r,
                "thr_p": pooled_thr_p,
            }
        for base in ("logistic f1+context", "xgboost f1+context"):
            alt = base + " no-pca"
            if base in held and alt in models:
                models[alt]["delta_vs_with_pca"] = _paired_delta(held[base], held[alt], n_boot)
        report_classes[kind] = {
            "n": len(group),
            "n_pos": int(sum(y)),
            "n_neg": int(len(y) - sum(y)),
            "models": models,
            "importance": importance,
            "vs_sr": _vs_rule(models),
        }
    _apply_operating_points(report_classes)
    for kind, block in report_classes.items():
        if "models" in block:
            block["auto"] = _auto_decision(block["models"])
    return {
        "protocol": {
            "outer": outer_name,
            "outer_splits": len(folds),
            "outer_group": outer_group,
            "outer_stratify": stratify,
            "inner": "GroupKFold",
            "inner_splits": INNER_SPLITS,
            "inner_group": group_col,
            "bootstrap": "sample",
            "bootstrap_samples": n_boot,
            "seed": _seed(),
            "n_jobs": N_JOBS,
            "tree_method": "hist",
        },
        "notes": {
            "frac_blacklist": "constant 0.0 on every ML-0 row; not used as a feature",
            "coverage_ratio": "NA on every ML-0 row; kept, so logistic/RF see a missing indicator and XGBoost sees NaN",
            "paper_presets": (
                "Included as extra XGBoost/RF candidates. max_depth 143 (complex) and 84 (deletion) "
                "are far above the grid {3, 6, 10} and are not defaults."
            ),
            "labeler_columns_excluded": sorted(LABELER_COLUMNS),
        },
        "library_versions": library_versions(),
        "excluded": excluded,
        "folds": [{"train": train, "test": test} for train, test in folds],
        "classes": report_classes,
    }


def _align(left: dict, right: dict):
    common = [i for i in left["i"] if i in set(right["i"])]
    pos_r = {i: k for k, i in enumerate(right["i"])}
    pos_l = {i: k for k, i in enumerate(left["i"])}
    idx_l = [pos_l[i] for i in common]
    idx_r = [pos_r[i] for i in common]
    return idx_l, idx_r


def _paired_delta(with_pca: dict, no_pca: dict, n_boot: int) -> dict:
    """no-pca minus with-pca. A negative AP means dropping the four PCA columns hurt."""
    idx_l, idx_r = _align(with_pca, no_pca)
    y = [with_pca["y"][i] for i in idx_l]
    s_with = [with_pca["s"][i] for i in idx_l]
    s_no = [no_pca["s"][i] for i in idx_r]
    samples = [with_pca["sample"][i] for i in idx_l]
    thr_with_r = [with_pca["thr_r"][i] for i in idx_l]
    thr_no_r = [no_pca["thr_r"][i] for i in idx_r]
    thr_with_p = [with_pca["thr_p"][i] for i in idx_l]
    thr_no_p = [no_pca["thr_p"][i] for i in idx_r]

    def point(scores, thr_r, thr_p):
        ap = _average_precision(y, scores)
        p90, r90 = _at_rows(y, scores, thr_r)
        p95, r95 = _at_rows(y, scores, thr_p)
        return ap, p90, r90, r95

    ap_w, p90_w, r90_w, r95_w = point(s_with, thr_with_r, thr_with_p)
    ap_n, p90_n, r90_n, r95_n = point(s_no, thr_no_r, thr_no_p)

    def sub(a, b):
        if a is None or b is None:
            return None
        return float(a) - float(b)

    by: dict[str, list[int]] = {}
    for i, sample in enumerate(samples):
        by.setdefault(sample, []).append(i)
    keys = sorted(by)
    rng = np.random.default_rng(_seed())
    draws = {"ap": [], "precision_at_recall90": [], "recall_at_recall90": [], "recall_at_precision95": []}
    for _ in range(n_boot):
        drawn = rng.choice(keys, size=len(keys), replace=True)
        idx: list[int] = []
        for sample in drawn:
            idx.extend(by[sample])
        y_b = [y[i] for i in idx]
        ap_delta = sub(
            _average_precision(y_b, [s_no[i] for i in idx]),
            _average_precision(y_b, [s_with[i] for i in idx]),
        )
        p_no, r_no = _at_rows(y_b, [s_no[i] for i in idx], [thr_no_r[i] for i in idx])
        p_w, r_w = _at_rows(y_b, [s_with[i] for i in idx], [thr_with_r[i] for i in idx])
        _p95_no, r95_no = _at_rows(y_b, [s_no[i] for i in idx], [thr_no_p[i] for i in idx])
        _p95_w, r95_w = _at_rows(y_b, [s_with[i] for i in idx], [thr_with_p[i] for i in idx])
        draws["ap"].append(ap_delta)
        draws["precision_at_recall90"].append(sub(p_no, p_w))
        draws["recall_at_recall90"].append(sub(r_no, r_w))
        draws["recall_at_precision95"].append(sub(r95_no, r95_w))
    return {
        "definition": "no-pca minus with-pca",
        "average_precision": sub(ap_n, ap_w),
        "precision_at_recall90": sub(p90_n, p90_w),
        "recall_at_recall90": sub(r90_n, r90_w),
        "recall_at_precision95": sub(r95_n, r95_w),
        "ci": {key: _ci(values) for key, values in draws.items()},
    }


def _vs_rule(models: dict) -> list[dict]:
    rule = models.get("sr>0") or {}
    rule_ci = (rule.get("ci") or {}).get("ap") or [None, None]
    out = []
    for name, model in models.items():
        if name == "sr>0":
            continue
        ci = (model.get("ci") or {}).get("ap") or [None, None]
        beats = (
            ci[0] is not None
            and rule_ci[1] is not None
            and float(ci[0]) > float(rule_ci[1])
        )
        out.append(
            {
                "model": name,
                "beats_sr_ap_outside_ci": beats,
                "ap": model.get("average_precision"),
                "ap_ci": ci,
                "sr_ap_ci": rule_ci,
            }
        )
    return out


def _finish_model(y, scores, samples, superpops, thr_r, thr_p, params, errors, n_boot):
    """Pool outer test rows. Each row keeps the threshold chosen on its training fold."""
    precision90, recall90 = _at_rows(y, scores, thr_r)
    reachable = bool(thr_p) and all(threshold is not None for threshold in thr_p)
    if reachable:
        precision95, recall95 = _at_rows(y, scores, thr_p)
    else:
        precision95, recall95 = None, None
    by_sp = {}
    for sp in sorted(set(superpops)):
        mask = [i for i, value in enumerate(superpops) if value == sp]
        if not mask:
            continue
        prec, rec = _at_rows(
            [y[i] for i in mask],
            [scores[i] for i in mask],
            [thr_r[i] for i in mask],
        )
        by_sp[sp] = {"n": len(mask), "precision": prec, "recall": rec}
    paper_folds = sum(1 for item in params if item.get("source") == "paper")
    return {
        "n": len(y),
        "average_precision": _average_precision(y, scores),
        "roc_auc": _roc_auc(y, scores),
        "precision_at_recall90": precision90,
        "recall_at_recall90": recall90,
        "precision95_reachable": reachable,
        "precision_at_precision95": precision95,
        "recall_at_precision95": recall95 if reachable else None,
        "ci": _bootstrap(y, scores, samples, thr_r, thr_p, n_boot, _seed()),
        "by_superpopulation": by_sp,
        "chosen_params": params,
        "paper_preset_folds": paper_folds,
        "candidate_errors": errors[:8],
    }


def _apply_operating_points(report_classes: dict) -> None:
    return None


def _auto_decision(models: dict) -> dict:
    """Highest outer recall at precision >= 0.95; otherwise highest AP.

    A simpler model wins when its point estimate lies inside the leader's
    bootstrap CI. Simplicity order is logistic, then rf, then xgboost, and
    within a family ``f1`` before ``f1+context``. The ``sr>0`` rule is a
    baseline, not an auto candidate.
    """
    order = [
        "logistic f1",
        "logistic f1+context",
        "rf f1",
        "rf f1+context",
        "xgboost f1",
        "xgboost f1+context",
    ]
    present = [name for name in order if name in models and models[name].get("average_precision") is not None]
    reachable = [name for name in present if models[name].get("precision95_reachable")]
    metric = "recall_at_precision95" if reachable else "average_precision"
    pool = reachable or present
    if not pool:
        return {"choice": None, "why": "no fitted candidate"}

    def value(name: str) -> float:
        raw = models[name].get(metric)
        return float("-inf") if raw is None else float(raw)

    leader = max(pool, key=lambda name: (value(name), -order.index(name)))
    ci_key = "recall_at_precision95" if metric == "recall_at_precision95" else "ap"
    lo, hi = (models[leader].get("ci") or {}).get(ci_key) or (None, None)
    choice = leader
    if lo is not None and hi is not None:
        for name in order:
            if name not in pool:
                continue
            point = value(name)
            if lo - 1e-12 <= point <= hi + 1e-12:
                choice = name
                break
    if choice == leader:
        why = f"highest {metric} ({value(leader):.4f})"
    else:
        why = (
            f"{choice} ({value(choice):.4f}) is inside the {leader} {metric} "
            f"CI [{lo:.4f}, {hi:.4f}] and is simpler"
        )
    if not reachable:
        why = "no candidate reached precision 0.95; " + why
    return {"choice": choice, "metric": metric, "why": why, "leader": leader}


def _rule_report(group, y, folds, by_sample, n_boot: int) -> dict:
    pooled_y: list[int] = []
    pooled_s: list[float] = []
    pooled_sample: list[str] = []
    pooled_sp: list[str] = []
    thr_r: list[float | None] = []
    thr_p: list[float | None] = []
    for train_samples, test_samples in folds:
        train_idx = [i for sample in train_samples for i in by_sample.get(sample, [])]
        test_idx = [i for sample in test_samples for i in by_sample.get(sample, [])]
        if not test_idx or not train_idx:
            continue
        train_scores = [sr_score(group[i]) for i in train_idx]
        train_y = [y[i] for i in train_idx]
        recall_choice = choose_threshold(train_y, train_scores, min_recall=0.90)
        precision_choice = choose_precision_threshold(train_y, train_scores, min_precision=0.95)
        precision_threshold = None if precision_choice is None else precision_choice["threshold"]
        for i in test_idx:
            pooled_y.append(y[i])
            pooled_s.append(sr_score(group[i]))
            pooled_sample.append(group[i].get("sample") or "")
            pooled_sp.append(group[i].get("superpopulation") or "NA")
            thr_r.append(recall_choice.threshold)
            thr_p.append(precision_threshold)
    return _finish_model(pooled_y, pooled_s, pooled_sample, pooled_sp, thr_r, thr_p, [], [], n_boot)


def fit_selected(rows: list[dict[str, str]], *, features: str = "f1", model: str = "logistic") -> dict:
    """Train one model per class. Hyperparameters and the threshold use inner folds only."""
    if features not in ("f1", "f1+context"):
        raise ValueError(f"unknown features {features}")
    if model not in ("logistic", "xgboost", "rf", "auto"):
        raise ValueError(f"unknown model {model}")
    if features == "f1" and model == "logistic":
        raise ValueError("the default logistic f1 path stays in fit_labeled_calls")
    rows = [row for row in rows if row.get("exclude") != "1"]
    if model == "auto":
        report = evaluate(rows)
        choices = {
            kind: (block.get("auto") or {}).get("choice")
            for kind, block in report["classes"].items()
            if block.get("auto")
        }
    else:
        report = None
        choices = {kind: f"{model} {features}" for kind in CLASSES}
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        if row.get("label") not in ("0", "1"):
            continue
        kind = class_of_row(row)
        if kind is None:
            continue
        grouped.setdefault(kind, []).append(row)
    classifiers: dict[str, dict] = {}
    thresholds: dict[str, float] = {}
    pca_by_class: dict[str, dict] = {}
    decisions = {}
    for kind, group in grouped.items():
        choice = choices.get(kind)
        if not choice:
            continue
        family, feature_set = choice.split(" ", 1)
        y = [int(row["label"]) for row in group]
        if len(set(y)) < 2:
            continue
        cache = _class_cache(group)
        context = context_matrix(group) if feature_set == "f1+context" else None
        train_idx = list(range(len(group)))
        candidates = _candidate_list(kind, family)
        pca_cache: dict = {}
        best, _errors = _select_params(
            group, y, cache, context, train_idx, feature_set, family, candidates, False, pca_cache
        )
        if best is None:
            continue
        recall_choice = choose_threshold(y, best["scores"], min_recall=0.90)
        kind_m, fitted, pca, medians, names = _fit_fold(
            family, best["params"], cache, context, train_idx, y, feature_set, False, pca_cache
        )
        spec = _pack_model(kind_m, fitted, names, medians, pca, best["params"])
        spec["feature_set"] = feature_set
        classifiers[kind] = spec
        thresholds[kind] = recall_choice.threshold
        if pca is not None:
            pca_by_class[kind] = pca.to_dict()
        decisions[kind] = {"choice": choice, "params": _public_params(best["params"])}
    schema = FEATURE_SCHEMA_F2 if features == "f1+context" or any(
        (spec.get("feature_set") == "f1+context") for spec in classifiers.values()
    ) else FEATURE_SCHEMA
    # auto may mix feature sets. The manifest schema is f2 if any class uses context.
    if any(spec.get("feature_set") == "f1+context" for spec in classifiers.values()):
        schema = FEATURE_SCHEMA_F2
    return {
        "feature_schema": schema,
        "fitted": bool(classifiers),
        "models": {kind: classifiers[kind]["format"] for kind in classifiers},
        "features_flag": features,
        "model_flag": model,
        "library_versions": library_versions(),
        "pca": pca_by_class,
        "classifiers": classifiers,
        "thresholds": thresholds,
        "decisions": decisions,
        "evaluation": None if report is None else _public_report(report),
    }


def _public_report(report: dict) -> dict:
    """Evaluation without the sample fold lists, which are large but not secret."""
    return {key: value for key, value in report.items() if key != "folds"}


def score_rows(rows: list[dict[str, str]], manifest: dict) -> list[float | None]:
    """Probabilities aligned to ``rows``. Insertions and missing classes are None."""
    classifiers = manifest.get("classifiers") or {}
    out: list[float | None] = [None] * len(rows)
    by_kind: dict[str, list[int]] = {}
    for i, row in enumerate(rows):
        kind = class_of_row(row)
        if kind is None or kind not in classifiers:
            continue
        by_kind.setdefault(kind, []).append(i)
    for kind, indices in by_kind.items():
        spec = classifiers[kind]
        probs = _score_spec([rows[i] for i in indices], spec)
        for index, prob in zip(indices, probs):
            out[index] = float(prob)
    return out


def _score_spec(rows: list[dict[str, str]], spec: dict) -> np.ndarray:
    from arcsv2.features import PcaModel

    fmt = spec.get("format") or "logistic"
    pca = PcaModel.from_dict(spec["pca"]) if spec.get("pca") else None
    names = spec.get("features") or []
    if fmt == "logistic" and spec.get("context") != "median" and "model_json" not in spec:
        # Legacy f1 logistic: the caller's existing per-row path handles this.
        matrix = []
        for row in rows:
            feats = featurize_row(row, pca)
            matrix.append([float(feats.get(name, 0.0)) for name in names])
        return _predict_fitted("logistic", spec, np.asarray(matrix, dtype=float))
    kind = "logistic" if fmt == "logistic" else ("xgboost" if fmt == "xgboost-json" else "rf")
    if fmt == "xgboost-json":
        model = _load_xgb(spec)
    elif fmt == "skops":
        model = _load_rf(spec)
    else:
        model = spec
    medians = spec.get("medians")
    matrix = np.zeros((len(rows), len(names)), dtype=float)
    ctx_pos = {name: i for i, name in enumerate(CONTEXT_FEATURES)}
    for r, row in enumerate(rows):
        feats = featurize_row(row, pca)
        for c, name in enumerate(names):
            base = name[: -len("_missing")] if name.endswith("_missing") else name
            if name.endswith("_missing") and base in ctx_pos:
                raw = parse_context_value(row.get(base), base)
                matrix[r, c] = 1.0 if math.isnan(raw) else 0.0
            elif name in ctx_pos:
                raw = parse_context_value(row.get(name), name)
                if math.isnan(raw):
                    if kind == "xgboost":
                        matrix[r, c] = np.nan
                    else:
                        matrix[r, c] = float(medians[ctx_pos[name]]) if medians is not None else 0.0
                else:
                    matrix[r, c] = raw
            else:
                matrix[r, c] = float(feats.get(name, 0.0))
    if kind == "logistic":
        return _predict_fitted("logistic", spec, matrix)
    return _predict_fitted(kind, model, matrix)


def write_evaluation(
    labels: Path,
    out: Path,
    *,
    max_samples: int | None = None,
    n_boot: int = N_BOOT,
    outer_group: str = "sample",
    include_rf: str = "inversion",
    ablation: bool = True,
    exclude_splits: list[str] | None = None,
) -> dict:
    """Run ``evaluate`` and write ``report.json`` plus ``report.md``."""
    import json
    import time

    rows = load_table(labels)
    subset = None
    if max_samples is not None:
        rows, superpop, chosen = subset_samples(rows, max_samples)
        subset = {"superpopulation": superpop, "samples": chosen, "max_samples": max_samples}
        _log(f"subset {superpop} samples={','.join(chosen)} rows={len(rows)}")
    stop, peak = _rss_sampler()
    started = time.perf_counter()
    try:
        report = evaluate(
            rows,
            n_boot=n_boot,
            outer_group=outer_group,
            include_rf=include_rf,
            ablation=ablation,
            exclude_splits=exclude_splits,
        )
    finally:
        stop.set()
    report["wall_s"] = time.perf_counter() - started
    report["peak_rss_bytes"] = peak["rss"] or _current_rss()
    report["subset"] = subset
    report["labels"] = str(labels)
    out.mkdir(parents=True, exist_ok=True)
    (out / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (out / "report.md").write_text(render_markdown(report), encoding="utf-8")
    _log(f"wrote {out / 'report.json'} wall_s={report['wall_s']:.1f} peak_rss={report['peak_rss_bytes']}")
    return report


def _current_rss() -> int | None:
    try:
        import psutil
    except ImportError:
        return None
    return int(psutil.Process().memory_info().rss)


def _rss_sampler():
    import threading

    stop = threading.Event()
    peak = {"rss": 0}

    def loop() -> None:
        while not stop.is_set():
            rss = _current_rss()
            if rss and rss > peak["rss"]:
                peak["rss"] = rss
            stop.wait(0.5)

    threading.Thread(target=loop, daemon=True).start()
    return stop, peak


def _protocol_sentence(report: dict) -> str:
    protocol = report.get("protocol") or {}
    if protocol.get("outer_group") == "chrom":
        return (
            "Outer `GroupKFold` by chromosome, balanced by row count. "
            "Hyperparameters and thresholds use inner `GroupKFold` on the training chromosomes. "
            "Intervals are 95% bootstrap-over-samples percentiles."
        )
    return (
        "Outer `StratifiedGroupKFold` by sample, stratified by superpopulation. "
        "Hyperparameters and thresholds are chosen with inner `GroupKFold` on the "
        "training samples only. Intervals are 95% bootstrap-over-samples percentiles."
    )


def render_markdown(report: dict) -> str:
    lines = [
        "# ML-1 model comparison",
        "",
        _protocol_sentence(report),
        "",
        report["notes"]["frac_blacklist"] + ".",
        report["notes"]["coverage_ratio"] + ".",
        report["notes"]["paper_presets"],
        "",
    ]
    excluded = report.get("excluded") or {}
    if excluded.get("rows"):
        lines.append(
            f"Excluded split `{', '.join(excluded.get('splits') or [])}`: "
            f"{excluded['samples']} samples, {excluded['rows']} rows."
        )
        lines.append("")
    lines.extend([
        "## Per class",
        "",
    ])
    header = (
        "| Model | n | AP | AP CI | ROC AUC | Precision @ recall≥0.90 | Recall @ recall≥0.90 | "
        "Recall @ precision≥0.95 |"
    )
    sep = "|:--|---:|---:|:--|---:|---:|---:|---:|"
    for kind, block in report["classes"].items():
        lines.append(f"### {kind}")
        lines.append("")
        if "models" not in block:
            lines.append(f"Skipped ({block.get('status', 'no rows')}).")
            lines.append("")
            continue
        lines.append(
            f"n = {block['n']} ({block['n_pos']} positive / {block['n_neg']} negative)."
        )
        lines.append("")
        lines.append(header)
        lines.append(sep)
        for name, model in block["models"].items():
            lines.append(
                "| {name} | {n} | {ap} | {ci} | {auc} | {p90} | {r90} | {r95} |".format(
                    name=name,
                    n=model.get("n"),
                    ap=_fmt(model.get("average_precision")),
                    ci=_fmt_ci(model.get("ci", {}).get("ap")),
                    auc=_fmt(model.get("roc_auc")),
                    p90=_fmt(model.get("precision_at_recall90")),
                    r90=_fmt(model.get("recall_at_recall90")),
                    r95=_fmt(model.get("recall_at_precision95"))
                    if model.get("precision95_reachable")
                    else "not reachable",
                )
            )
        lines.append("")
        decision = block.get("auto") or {}
        lines.append(f"Auto choice: `{decision.get('choice')}` — {decision.get('why')}.")
        lines.append("")
        for name, model in block["models"].items():
            delta = model.get("delta_vs_with_pca")
            if not delta:
                continue
            lines.append(
                f"`{name}` minus the same model with PCA: AP {_fmt(delta.get('average_precision'))} "
                f"(CI {_fmt_ci(delta.get('ci', {}).get('ap'))}), "
                f"precision @ recall≥0.90 {_fmt(delta.get('precision_at_recall90'))} "
                f"(CI {_fmt_ci(delta.get('ci', {}).get('precision_at_recall90'))})."
            )
        if any(model.get("delta_vs_with_pca") for model in block["models"].values()):
            lines.append("")
        beats = [row for row in block.get("vs_sr") or [] if row.get("beats_sr_ap_outside_ci")]
        if kind == "inversion":
            if beats:
                names = ", ".join(row["model"] for row in beats)
                lines.append(f"Beats `sr>0` on AP outside the bootstrap CI: {names}.")
            else:
                lines.append(
                    "No model beats `sr>0` on AP outside the bootstrap CI. "
                    "Prefer the rule and more labeled inversions over a fragile model."
                )
            rule = block["models"].get("sr>0") or {}
            ci = (rule.get("ci") or {}).get("ap") or [None, None]
            if ci[0] is not None and ci[1] is not None:
                lines.append(f"`sr>0` AP CI width is {ci[1] - ci[0]:.3f} ([{ci[0]:.3f}, {ci[1]:.3f}]).")
            lines.append("")
        if block.get("importance"):
            lines.append("XGBoost gain (mean over outer training fits), top 15:")
            lines.append("")
            for label, rows in block["importance"].items():
                if not rows:
                    continue
                lines.append(f"`{label}`: " + ", ".join(f"{row['feature']} ({row['gain']:.3g})" for row in rows))
            lines.append("")
    lines.append("## Precision at the recall ≥ 0.90 threshold, by superpopulation")
    lines.append("")
    for kind, block in report["classes"].items():
        models = block.get("models") or {}
        if not models:
            continue
        lines.append(f"### {kind}")
        lines.append("")
        supers = sorted({sp for model in models.values() for sp in (model.get("by_superpopulation") or {})})
        lines.append("| Model | " + " | ".join(supers) + " |")
        lines.append("|:--|" + "|".join("---:" for _ in supers) + "|")
        for name, model in models.items():
            cells = []
            for sp in supers:
                cell = (model.get("by_superpopulation") or {}).get(sp)
                cells.append(_fmt(cell["precision"]) if cell else "")
            lines.append("| " + name + " | " + " | ".join(cells) + " |")
        lines.append("")
    return "\n".join(lines) + "\n"


def _fmt(value) -> str:
    if value is None:
        return ""
    return f"{float(value):.3f}"


def _fmt_ci(pair) -> str:
    if not pair or pair[0] is None:
        return ""
    return f"[{pair[0]:.3f}, {pair[1]:.3f}]"
