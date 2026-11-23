"""Call-table features (§7.2). Schema ``arcsv2-f1``.

The paper's Table S1 is not public. These features are the reconstruction in
the roadmap: Tier A from the call table, a fixed-layout block PCA, and optional
Tier B diagnostics. Nothing is padded to 48. Missing values are 0 with a
``*_missing`` flag. PCA is fit on training rows only.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass

import numpy as np

from arcsv2.config import FEATURE_SCHEMA

TOKEN_RE = re.compile(r"([A-Za-z]\d*)('?)|(_)|(=)")
_BLOCK_SLOTS = 12
# 2 orientations × 12 blocks, plus insertion and translocation.
_ALT_CHANNELS = 2 * _BLOCK_SLOTS + 2

TIER_A: tuple[str, ...] = (
    "num_sv",
    "n_blocks_ref",
    "n_tokens_alt",
    "n_deleted_blocks",
    "n_duplicated_blocks",
    "n_inverted_blocks",
    "n_novel_adjacencies",
    "len_affected",
    "log10_len_affected",
    "span",
    "log10_span",
    "ref_window",
    "inner_block_len_min",
    "inner_block_len_max",
    "inner_block_len_mean",
    "flank_len_min",
    "len_affected_over_span",
    "is_het",
    "af",
    "score_vs_ref",
    "score_vs_next",
    "score_vs_ref_per_bp",
    "score_vs_ref_per_kb",
    "score_ratio_next_ref",
    "num_paths",
    "log10_num_paths",
    "next_is_ref_het",
    "split_support_sum",
    "split_support_min",
    "split_support_max",
    "split_support_mean",
    "pair_support_sum",
    "pair_support_min",
    "pair_support_max",
    "pair_support_mean",
    "min_bp_total_support",
    "n_bp_no_split_support",
    "n_bp_no_pair_support",
    "bp_uncertainty_sum",
    "bp_uncertainty_min",
    "bp_uncertainty_max",
    "bp_uncertainty_mean",
    "sv_bp_uncertainty_max",
    "ins_total",
    "has_insertion",
    "inslen_lower_bound",
)

PCA_FEATURES = ("rearr_pc1", "rearr_pc2", "ref_pc1", "ref_pc2")


def _tokens(text: str) -> list[tuple[str, bool]]:
    out: list[tuple[str, bool]] = []
    for match in TOKEN_RE.finditer(text or ""):
        if match.group(1):
            out.append((match.group(1), bool(match.group(2))))
        elif match.group(3):
            out.append(("_", False))
        elif match.group(4):
            out.append(("=", False))
    return out


def _block_slot(name: str) -> int | None:
    if name in ("_", "="):
        return None
    # A..Z a..z then A1... Map the first 12 reference names only.
    order = "ABCDEFGHIJKLabcdefghijkl"
    if name in order:
        return order.index(name)
    return None


def _inslen(text: str) -> tuple[list[float], bool]:
    """Insertion lengths. A ``>`` prefix marks a lower bound, not an exact size."""
    bound = ">" in (text or "")
    values = _nums((text or "").replace(">", ""))
    return values, bound


def _nums(text: str) -> list[float]:
    if text is None or text in ("", "NA"):
        return []
    parts = re.split(r"[,;]", str(text))
    values = []
    for part in parts:
        part = part.strip()
        if not part or part == "NA":
            continue
        try:
            values.append(float(part))
        except ValueError:
            continue
    return values


def _stats(values: list[float]) -> tuple[float, float, float, float]:
    if not values:
        return 0.0, 0.0, 0.0, 0.0
    return (
        float(sum(values)),
        float(min(values)),
        float(max(values)),
        float(sum(values) / len(values)),
    )


def rearrangement_matrix(reference: str, rearrangement: str) -> np.ndarray:
    """12 positions × (2×12 blocks + '_' + '=')."""
    matrix = np.zeros((_BLOCK_SLOTS, _ALT_CHANNELS), dtype=float)
    for i, (name, inverted) in enumerate(_tokens(rearrangement)[:_BLOCK_SLOTS]):
        if name == "_":
            matrix[i, -2] = 1.0
        elif name == "=":
            matrix[i, -1] = 1.0
        else:
            slot = _block_slot(name)
            if slot is None or slot >= _BLOCK_SLOTS:
                continue
            col = slot + ( _BLOCK_SLOTS if inverted else 0)
            matrix[i, col] = 1.0
    return matrix


def reference_indicator(reference: str) -> np.ndarray:
    vec = np.zeros(_BLOCK_SLOTS, dtype=float)
    for i, (name, _inv) in enumerate(_tokens(reference)[:_BLOCK_SLOTS]):
        slot = _block_slot(name)
        if slot is not None and slot < _BLOCK_SLOTS:
            vec[i] = 1.0
    return vec


@dataclass
class PcaModel:
    """Centered SVD components fit on training calls only."""

    rearr_mean: np.ndarray
    rearr_components: np.ndarray  # (2, flattened)
    ref_mean: np.ndarray
    ref_components: np.ndarray

    def transform_rearr(self, matrix: np.ndarray) -> tuple[float, float]:
        flat = matrix.ravel() - self.rearr_mean
        scores = self.rearr_components @ flat
        return float(scores[0]), float(scores[1] if len(scores) > 1 else 0.0)

    def transform_ref(self, vec: np.ndarray) -> tuple[float, float]:
        centered = vec - self.ref_mean
        scores = self.ref_components @ centered
        return float(scores[0]), float(scores[1] if len(scores) > 1 else 0.0)

    def to_dict(self) -> dict[str, list]:
        return {
            "rearr_mean": self.rearr_mean.tolist(),
            "rearr_components": self.rearr_components.tolist(),
            "ref_mean": self.ref_mean.tolist(),
            "ref_components": self.ref_components.tolist(),
        }

    @classmethod
    def from_dict(cls, raw: dict) -> "PcaModel":
        return cls(
            rearr_mean=np.asarray(raw["rearr_mean"], dtype=float),
            rearr_components=np.asarray(raw["rearr_components"], dtype=float),
            ref_mean=np.asarray(raw["ref_mean"], dtype=float),
            ref_components=np.asarray(raw["ref_components"], dtype=float),
        )


def fit_pca_arrays(rearr: np.ndarray, refs: np.ndarray) -> PcaModel:
    """PCA on stacked rearrangement rows and reference-indicator rows.

    ``rearr`` is ``rearrangement_matrix(...).ravel()`` per call. ``refs`` is
    ``reference_indicator(...)`` per call. Each SVD runs once.
    """
    rearr_mean, rearr_components = _center_components(rearr)
    ref_mean, ref_components = _center_components(refs)
    return PcaModel(rearr_mean, rearr_components, ref_mean, ref_components)


def fit_pca(rows: list[dict[str, str]]) -> PcaModel:
    if rows:
        rearr = np.stack(
            [
                rearrangement_matrix(r.get("reference", ""), r.get("rearrangement", "")).ravel()
                for r in rows
            ]
        )
        refs = np.stack([reference_indicator(r.get("reference", "")) for r in rows])
    else:
        rearr = np.zeros((1, _BLOCK_SLOTS * _ALT_CHANNELS))
        refs = np.zeros((1, _BLOCK_SLOTS))
    return fit_pca_arrays(rearr, refs)


def _center_components(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mean = matrix.mean(axis=0)
    centered = matrix - mean
    # Economy SVD. Components are the first two right singular vectors.
    try:
        _u, _s, vt = np.linalg.svd(centered, full_matrices=False)
    except np.linalg.LinAlgError:
        vt = np.zeros((2, matrix.shape[1]))
    comps = np.zeros((2, matrix.shape[1]))
    n = min(2, vt.shape[0])
    comps[:n] = vt[:n]
    return mean, comps


def _finite(value: float) -> tuple[float, float]:
    if value is None or not math.isfinite(value):
        return 0.0, 1.0
    return float(value), 0.0


def featurize_row(row: dict[str, str], pca: PcaModel | None = None) -> dict[str, float]:
    """One call → a finite feature dict. PCA features are 0 until a model is supplied."""
    ref = _tokens(row.get("reference", ""))
    alt = _tokens(row.get("rearrangement", ""))
    ref_names = [n for n, _ in ref if n not in ("_", "=")]
    alt_names = [n for n, _inv in alt]
    n_deleted = len([n for n in ref_names if n not in set(alt_names)])
    n_dup = len([n for n in alt_names if alt_names.count(n) > ref_names.count(n) and n not in ("_", "=")])
    n_inv = len([1 for n, inv in alt if inv])
    n_novel = 0
    for (n1, s1), (n2, s2) in zip(alt, alt[1:]):
        slot1 = _block_slot(n1)
        slot2 = _block_slot(n2)
        consecutive = (
            s1 is False
            and s2 is False
            and slot1 is not None
            and slot2 is not None
            and slot2 == slot1 + 1
        )
        if not consecutive:
            n_novel += 1
    bp = _nums(row.get("bp", ""))
    span = (float(row["maxbp"]) - float(row["minbp"])) if row.get("maxbp") not in (None, "NA") else 0.0
    ref_window = (bp[-1] - bp[0]) if len(bp) >= 2 else 0.0
    # Inner blocks are the bp gaps excluding the two flanks.
    inner = [bp[i + 1] - bp[i] for i in range(1, len(bp) - 2)] if len(bp) >= 4 else []
    flanks = []
    if len(bp) >= 2:
        flanks.append(bp[1] - bp[0])
        flanks.append(bp[-1] - bp[-2])
    length = float(row.get("len_affected") or 0)
    num_sv = float(row.get("num_sv") or 1)
    score_ref = float(row.get("score_vs_ref") or 0)
    score_next = float(row.get("score_vs_next") or 0)
    num_paths = float(row.get("num_paths") or 0)
    af = float(row.get("af") or 0)
    sr = _nums(row.get("sr_support", ""))
    pe = _nums(row.get("pe_support", ""))
    unc = _nums(row.get("bp_uncertainty", ""))
    sv_unc = _nums(row.get("sv_bp_uncertainty", ""))
    ins, ins_bound = _inslen(row.get("inslen", ""))
    sr_s = _stats(sr)
    pe_s = _stats(pe)
    unc_s = _stats(unc)
    inner_s = _stats(inner)
    total_support = [a + b for a, b in zip(sr or [0], pe or [0])]
    features = {
        "num_sv": num_sv,
        "n_blocks_ref": float(len(ref_names)),
        "n_tokens_alt": float(len(alt)),
        "n_deleted_blocks": float(n_deleted),
        "n_duplicated_blocks": float(n_dup),
        "n_inverted_blocks": float(n_inv),
        "n_novel_adjacencies": float(n_novel),
        "len_affected": length,
        "log10_len_affected": math.log10(length) if length > 0 else 0.0,
        "span": span,
        "log10_span": math.log10(span) if span > 0 else 0.0,
        "ref_window": ref_window,
        "inner_block_len_min": inner_s[1],
        "inner_block_len_max": inner_s[2],
        "inner_block_len_mean": inner_s[3],
        "flank_len_min": float(min(flanks)) if flanks else 0.0,
        "len_affected_over_span": (length / span) if span else 0.0,
        "is_het": 1.0 if row.get("gt") == "HET" else 0.0,
        "af": af,
        "score_vs_ref": score_ref,
        "score_vs_next": score_next,
        "score_vs_ref_per_bp": score_ref / num_sv if num_sv else 0.0,
        "score_vs_ref_per_kb": score_ref / (length / 1000.0) if length else 0.0,
        "score_ratio_next_ref": score_next / score_ref if score_ref else 0.0,
        "num_paths": num_paths,
        "log10_num_paths": math.log10(num_paths) if num_paths > 0 else 0.0,
        "next_is_ref_het": 1.0 if "/" in (row.get("rearrangement_next") or "") and (row.get("reference") or "") in (row.get("rearrangement_next") or "").split("/") else 0.0,
        "split_support_sum": sr_s[0],
        "split_support_min": sr_s[1],
        "split_support_max": sr_s[2],
        "split_support_mean": sr_s[3],
        "pair_support_sum": pe_s[0],
        "pair_support_min": pe_s[1],
        "pair_support_max": pe_s[2],
        "pair_support_mean": pe_s[3],
        "min_bp_total_support": float(min(total_support)) if total_support else 0.0,
        "n_bp_no_split_support": float(sum(1 for v in sr if v == 0)),
        "n_bp_no_pair_support": float(sum(1 for v in pe if v == 0)),
        "bp_uncertainty_sum": unc_s[0],
        "bp_uncertainty_min": unc_s[1],
        "bp_uncertainty_max": unc_s[2],
        "bp_uncertainty_mean": unc_s[3],
        "sv_bp_uncertainty_max": float(max(sv_unc)) if sv_unc else 0.0,
        # A leading '>' means the number is a minimum. It is not used as a point length.
        "ins_total": 0.0 if ins_bound else (float(sum(ins)) if ins else 0.0),
        "has_insertion": 1.0 if "_" in (row.get("rearrangement") or "") else 0.0,
        "inslen_lower_bound": 1.0 if ins_bound else 0.0,
    }
    cleaned: dict[str, float] = {}
    for name, value in features.items():
        finite, missing = _finite(value)
        cleaned[name] = finite
        if missing:
            cleaned[f"{name}_missing"] = 1.0
    if pca is None:
        cleaned["rearr_pc1"] = 0.0
        cleaned["rearr_pc2"] = 0.0
        cleaned["ref_pc1"] = 0.0
        cleaned["ref_pc2"] = 0.0
    else:
        r1, r2 = pca.transform_rearr(rearrangement_matrix(row.get("reference", ""), row.get("rearrangement", "")))
        f1, f2 = pca.transform_ref(reference_indicator(row.get("reference", "")))
        cleaned["rearr_pc1"], cleaned["rearr_pc2"] = r1, r2
        cleaned["ref_pc1"], cleaned["ref_pc2"] = f1, f2
    cleaned["feature_schema_ok"] = 1.0
    return cleaned


def feature_names(sample: dict[str, float] | None = None) -> list[str]:
    names = list(TIER_A) + list(PCA_FEATURES)
    if sample:
        extra = [k for k in sample if k.endswith("_missing")]
        names.extend(sorted(extra))
    return names


SCHEMA = FEATURE_SCHEMA
