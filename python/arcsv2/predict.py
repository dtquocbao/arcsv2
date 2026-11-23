"""Apply a trained bundle (§7.4).

Refuses a mismatched ``FEATURE_SCHEMA``. Classes without a model (insertions)
keep ``ml_prob = NA`` and are kept by default. Passing rows get ``FILTER``
left as PASS; failures are marked so the VCF writer can set ``FILTER=MLFAIL``.
"""

from __future__ import annotations

import math

from arcsv2.config import FEATURE_SCHEMA
from arcsv2.features import CONTEXT_FEATURES, FEATURE_SCHEMA_F2, PcaModel, featurize_row
from arcsv2.train import class_of_row


def _logistic_prob(spec: dict, feats: dict[str, float]) -> float:
    mean = spec.get("mean")
    scale = spec.get("scale")
    z = float(spec["intercept"])
    for i, (name, coef) in enumerate(zip(spec["features"], spec["coef"])):
        value = float(feats.get(name, 0.0))
        if mean is not None and scale is not None:
            denom = float(scale[i]) or 1.0
            value = (value - float(mean[i])) / denom
        z += float(coef) * value
    if z > 60.0:
        return 1.0
    if z < -60.0:
        return 0.0
    return 1.0 / (1.0 + math.exp(-z))


def _has_context(rows: list[dict[str, str]]) -> bool:
    return bool(rows) and all(name in rows[0] for name in CONTEXT_FEATURES)


def _legacy_logistic(manifest: dict) -> bool:
    if manifest.get("feature_schema") != FEATURE_SCHEMA:
        return False
    for spec in (manifest.get("classifiers") or {}).values():
        if not isinstance(spec, dict):
            continue
        if spec.get("format") not in (None, "logistic"):
            return False
        if spec.get("context") == "median" or spec.get("model_json") or spec.get("model_b64"):
            return False
    return True


def predict_rows(rows: list[dict[str, str]], manifest: dict, pca=None) -> list[dict[str, str]]:
    schema = manifest.get("feature_schema")
    if schema not in (FEATURE_SCHEMA, FEATURE_SCHEMA_F2):
        raise ValueError(f"FEATURE_SCHEMA mismatch: bundle has {schema}, code has {FEATURE_SCHEMA}")
    if schema == FEATURE_SCHEMA_F2 and not _has_context(rows):
        raise ValueError(
            f"schema mismatch: manifest is {FEATURE_SCHEMA_F2} but input is {FEATURE_SCHEMA}"
        )
    if _legacy_logistic(manifest):
        return _predict_legacy(rows, manifest, pca)
    return _predict_fitted(rows, manifest, pca)


def _predict_legacy(rows: list[dict[str, str]], manifest: dict, pca=None) -> list[dict[str, str]]:
    raw_pca = manifest.get("pca")
    if pca is None and isinstance(raw_pca, dict) and "rearr_mean" in raw_pca:
        pca = PcaModel.from_dict(raw_pca)
    thresholds: dict[str, float] = manifest.get("thresholds", {})
    classifiers: dict = manifest.get("classifiers") or {}
    out: list[dict[str, str]] = []
    for row in rows:
        updated = dict(row)
        kind = class_of_row(row)
        updated["sv_class"] = kind or "insertion"
        if kind is None:
            updated["ml_prob"] = "NA"
            updated["ml_pass"] = "1"
            updated["ml_model"] = "none"
            out.append(updated)
            continue
        spec = classifiers.get(kind)
        if not spec:
            # An untrained bundle must not look like a probability from the caller score.
            updated["ml_prob"] = "NA"
            updated["ml_pass"] = "NA"
            updated["ml_model"] = "untrained"
            out.append(updated)
            continue
        class_pca = pca
        if spec.get("pca"):
            class_pca = PcaModel.from_dict(spec["pca"])
        feats = featurize_row(row, class_pca)
        prob = _logistic_prob(spec, feats)
        updated["ml_prob"] = f"{prob:.6f}"
        updated["ml_class"] = kind
        updated["ml_model"] = "logistic"
        passed = prob >= float(thresholds.get(kind, 1.0))
        updated["ml_pass"] = "1" if passed else "0"
        out.append(updated)
    return out


def _predict_fitted(rows: list[dict[str, str]], manifest: dict, pca=None) -> list[dict[str, str]]:
    from arcsv2.ml1 import score_rows

    del pca
    probabilities = score_rows(rows, manifest)
    thresholds: dict[str, float] = manifest.get("thresholds", {})
    out: list[dict[str, str]] = []
    for row, prob in zip(rows, probabilities):
        updated = dict(row)
        kind = class_of_row(row)
        updated["sv_class"] = kind or "insertion"
        if kind is None or prob is None:
            updated["ml_prob"] = "NA"
            updated["ml_pass"] = "1" if kind is None else "NA"
            updated["ml_model"] = "none" if kind is None else "untrained"
            out.append(updated)
            continue
        updated["ml_prob"] = f"{prob:.6f}"
        updated["ml_class"] = kind
        fmt = (manifest.get("classifiers") or {}).get(kind, {}).get("format", "logistic")
        updated["ml_model"] = fmt
        updated["ml_pass"] = "1" if prob >= float(thresholds.get(kind, 1.0)) else "0"
        out.append(updated)
    return out
