"""Library models (§4.2).

Insert pmf: Gaussian kernel on the integer grid, bandwidth
``max(1, 0.785 * IQR / n**(1/5))`` (orig:bamparser_streaming.py::pmf_kernel_smooth).
The density is the standard normal kernel convolved with the integer histogram,
then normalized to sum to 1. sklearn's KernelDensity is the oracle for 1e-12
checks; this is the same estimator without the ``rtol`` tree approximation.

Discordance cutoffs use the paper likelihood-ratio constant
``C = exp(9/2)`` (log C = 4.5), scanned the same way as
orig:bamparser_streaming.py::get_lr_cutoff.

Robustness π is orig:call_sv.py::compute_pi_robust (``p = 1e-4``), then the
median across libraries.
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass, field


def _percentile(values: list[float], pct: float) -> float:
    if not values:
        raise ValueError("empty sample")
    xs = sorted(values)
    if len(xs) == 1:
        return float(xs[0])
    rank = (len(xs) - 1) * (pct / 100.0)
    lo = int(math.floor(rank))
    hi = int(math.ceil(rank))
    if lo == hi:
        return float(xs[lo])
    w = rank - lo
    return float(xs[lo] * (1.0 - w) + xs[hi] * w)


def _median(values: list[float]) -> float:
    return _percentile(values, 50)


def pmf_kernel_smooth(
    inserts: list[int],
    xmin: int = 0,
    xmax: int | None = None,
    max_kde_samples: int = 1_000_000,
    rng: random.Random | None = None,
) -> list[float]:
    """Gaussian KDE on ``xmin .. xmax``, normalized to a pmf (orig:pmf_kernel_smooth)."""
    if not inserts:
        raise ValueError("no insert lengths")
    sample = list(inserts)
    if len(sample) > max_kde_samples:
        rng = rng or random.Random(0)
        sample = rng.sample(sample, max_kde_samples)
    if xmax is None:
        xmax = max(sample)
    truncated = [int(x) for x in sample if xmin <= int(x) <= xmax]
    if not truncated:
        raise ValueError("no insert lengths inside the KDE window")
    q1 = _percentile(truncated, 25)
    q3 = _percentile(truncated, 75)
    iqr = q3 - q1
    n = len(truncated)
    bw = max(1.0, 0.785 * iqr / (n ** 0.2))
    width = xmax - xmin + 1
    hist = [0] * width
    for x in truncated:
        hist[x - xmin] += 1
    # Kernel half-width. 8 bw covers > 1 - 1e-15 of a standard normal.
    radius = max(1, int(math.ceil(8.0 * bw)))
    inv = 1.0 / (bw * math.sqrt(2.0 * math.pi) * n)
    kernel = [math.exp(-0.5 * ((k / bw) ** 2)) for k in range(-radius, radius + 1)]
    pmf = [0.0] * width
    for i, count in enumerate(hist):
        if count == 0:
            continue
        weight = count * inv
        lo = max(0, i - radius)
        hi = min(width, i + radius + 1)
        k0 = radius + (lo - i)
        for j in range(lo, hi):
            pmf[j] += weight * kernel[k0 + (j - lo)]
    total = sum(pmf)
    if total <= 0:
        raise ValueError("KDE total mass is zero")
    return [p / total for p in pmf]


def location_scale(inserts: list[int]) -> tuple[float, float]:
    """m = median, s = IQR / 1.35 (§4.2)."""
    m = _median([float(x) for x in inserts])
    q1 = _percentile([float(x) for x in inserts], 25)
    q3 = _percentile([float(x) for x in inserts], 75)
    s = (q3 - q1) / 1.35
    return m, s


def cdf_and_cdf_sum(pmf: list[float]) -> tuple[list[float], list[float]]:
    """F and F̃. Past the last bin, F stays 1 and F̃ grows by 1 per step.

    The closed form used by double_int matches
    orig:call_sv.py::create_insert_cs: for ``t >= len(pmf)``,
    ``F̃(t) = F̃(len-1) + (t - len + 1)``.
    """
    cdf: list[float] = []
    acc = 0.0
    for p in pmf:
        acc += p
        cdf.append(acc)
    cdf_sum: list[float] = []
    acc = 0.0
    for f in cdf:
        acc += f
        cdf_sum.append(acc)
    return cdf, cdf_sum


def eval_cdf_sum(cdf_sum: list[float], t: int) -> float:
    """orig:call_sv.py::create_insert_cs."""
    if t < 0:
        return 0.0
    if t >= len(cdf_sum):
        return cdf_sum[-1] + (t - len(cdf_sum) + 1)
    return cdf_sum[t]


def eval_pmf(pmf: list[float], t: int) -> float:
    if t < 0 or t >= len(pmf):
        return 0.0
    return pmf[t]


def double_int(
    cdf_sum: list[float],
    a1: int,
    b1: int,
    a2: int,
    b2: int,
    xmin: int = 0,
    xmax: int | None = None,
    ymin: int = 0,
    ymax: int | None = None,
) -> float:
    """O(1) sum of f(y-x) over a rectangle. orig:sv_inference_insertions.py::double_int."""
    if xmax is None:
        xmax = 10**18
    if ymax is None:
        ymax = 10**18
    if a1 > xmax or a2 > ymax or b1 < xmin or b2 < ymin:
        return 0.0
    a1 = max(xmin, a1)
    a2 = max(ymin, a2)
    b1 = min(xmax, b1)
    b2 = min(ymax, b2)
    if a1 > b1 or a2 > b2:
        return 0.0
    return (
        eval_cdf_sum(cdf_sum, b2 - a1 + 1)
        + eval_cdf_sum(cdf_sum, a2 - b1 - 1)
        - eval_cdf_sum(cdf_sum, b2 - b1)
        - eval_cdf_sum(cdf_sum, a2 - a1)
    )


def discordance_cutoffs(pmf: list[float], log_c: float = 4.5) -> tuple[int, int]:
    """Return ``(z_ins, z_del)``.

    Deletion pairs have insert ``z > z_del``. Insertion pairs have ``z < z_ins``.
    ``z*`` is the first grid point at which the log density ratio from the mode
    crosses ``log C`` and stays there, scanned as in orig:get_lr_cutoff.
    For a normal library, ``log C = 4.5`` lands on ``m ± 3s``.
    """
    if not pmf:
        raise ValueError("empty pmf")
    mode = max(pmf)
    if mode <= 0:
        raise ValueError("pmf mode is zero")
    logmode = math.log(mode)
    z_del = len(pmf) - 1
    for i in range(len(pmf) - 2, -1, -1):
        if pmf[i] > 0 and logmode - math.log(pmf[i]) < log_c:
            z_del = i + 1
            break
    z_ins = 0
    for i in range(1, len(pmf)):
        if pmf[i] > 0 and logmode - math.log(pmf[i]) < log_c:
            z_ins = i - 1
            break
    return z_ins, z_del


def pi_robust_one(pmf: list[float], p: float = 1e-4) -> float:
    """orig:call_sv.py::compute_pi_robust."""
    ordered = sorted(pmf, reverse=True)
    acc = 0.0
    target = 1.0 - p
    for value in ordered:
        acc += value
        if acc >= target:
            return value
    return ordered[-1]


def pi_robust(pmfs: list[list[float]], p: float = 1e-4) -> float:
    values = [pi_robust_one(pmf, p) for pmf in pmfs]
    return _median(values)


def quantile_index(cdf: list[float], q: float) -> int:
    for i, value in enumerate(cdf):
        if value >= q:
            return i
    return len(cdf) - 1


@dataclass
class LibraryModel:
    """One sequencing library (§4.2)."""

    lib_id: int
    name: str
    read_len_short: int
    read_len_long: int
    insert_median: float
    insert_iqr: float
    insert_scale: float
    pmf: list[float]
    cdf: list[float]
    cdf_sum: list[float]
    z_ins_cut: int
    z_del_cut: int
    pi_robust: float
    q_low: int
    q_high: int
    # Mapping-outcome class probabilities, orig:conditional_mappable_model.py
    # order (1,1), (1,0), (0,1), (1,-1), (-1,1).
    class_probs: list[float] = field(default_factory=lambda: [1.0, 0.0, 0.0, 0.0, 0.0])

    def f(self, insert: int) -> float:
        return eval_pmf(self.pmf, int(insert))

    @property
    def pmapped(self) -> float:
        return self.class_probs[0]


def model_from_mapstats(counts: dict[tuple[int, int], int]) -> list[float]:
    """orig:conditional_mappable_model.py::model_from_mapstats.

    Adds the mirrored observation so the five pair classes stay symmetric.
    """
    classes = [(1, 1), (1, 0), (0, 1), (1, -1), (-1, 1)]
    mirrored = {c: counts.get(c, 0) + counts.get((c[1], c[0]), 0) for c in classes}
    total = sum(mirrored.values())
    if total == 0:
        return [1.0, 0.0, 0.0, 0.0, 0.0]
    return [mirrored[c] / total for c in classes]


def fit_library(
    inserts: list[int],
    *,
    name: str = "lib1",
    lib_id: int = 0,
    read_len_short: int = 150,
    read_len_long: int = 150,
    mapstats: dict[tuple[int, int], int] | None = None,
    log_c: float = 4.5,
    robustness: float = 1e-4,
    q_low: float = 0.025,
    q_high: float = 0.975,
    max_mu_multiple: float = 3.0,
) -> LibraryModel:
    if not inserts:
        raise ValueError("zero usable pairs")
    median, scale = location_scale(inserts)
    q1 = _percentile([float(x) for x in inserts], 25)
    q3 = _percentile([float(x) for x in inserts], 75)
    xmax = int(math.ceil(max_mu_multiple * median))
    # orig:bamparser_streaming.py caps the KDE at 3 × median. Inserts past
    # that window are dropped inside pmf_kernel_smooth.
    pmf = pmf_kernel_smooth(inserts, 0, max(xmax, 1))
    cdf, cdf_sum = cdf_and_cdf_sum(pmf)
    z_ins, z_del = discordance_cutoffs(pmf, log_c)
    pi = pi_robust_one(pmf, robustness)
    return LibraryModel(
        lib_id=lib_id,
        name=name,
        read_len_short=read_len_short,
        read_len_long=read_len_long,
        insert_median=median,
        insert_iqr=q3 - q1,
        insert_scale=scale,
        pmf=pmf,
        cdf=cdf,
        cdf_sum=cdf_sum,
        z_ins_cut=z_ins,
        z_del_cut=z_del,
        pi_robust=pi,
        q_low=quantile_index(cdf, q_low),
        q_high=quantile_index(cdf, q_high),
        class_probs=model_from_mapstats(mapstats or {}),
    )


def sample_chunks(
    chunks: list[list[int]],
    target: int,
    rng: random.Random,
) -> list[int]:
    """Shuffled grid of chunks without replacement (§C1.2).

    Stops at ``target`` or when the grid is exhausted. Errors if nothing usable
    was found (the original resamples forever).
    """
    order = list(range(len(chunks)))
    rng.shuffle(order)
    taken: list[int] = []
    for idx in order:
        taken.extend(chunks[idx])
        if len(taken) >= target:
            return taken[:target]
    if not taken:
        raise RuntimeError("zero usable pairs")
    return taken


def pool_inserts(groups: list[list[int]]) -> list[int]:
    pooled: list[int] = []
    for group in groups:
        pooled.extend(group)
    return pooled
