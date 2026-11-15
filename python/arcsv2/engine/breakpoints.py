"""Breakpoint discovery (§4.3).

Soft clips, split-read types (orig:constants.py::SPLIT_TYPES), discordant-pair
cliques, scores, permutation FDR, uniform confidence intervals
(orig:unif_ci.py), and the merge rule in §4.3d.
"""

from __future__ import annotations

import math
import random
from collections import defaultdict
from dataclasses import dataclass

# orig:constants.py
SPLIT_FIRST_PLUS = 0x1
SPLIT_SECOND_PLUS = 0x2
SPLIT_OVERLAP = 0x4
SPLIT_LEFT_FIRST = 0x8
SPLIT_TYPES = {
    SPLIT_FIRST_PLUS | SPLIT_SECOND_PLUS: "Del+",
    0x0: "Del-",
    SPLIT_FIRST_PLUS | SPLIT_SECOND_PLUS | SPLIT_OVERLAP: "Dup+",
    SPLIT_OVERLAP: "Dup-",
    SPLIT_FIRST_PLUS | SPLIT_LEFT_FIRST: "InvL+",
    SPLIT_FIRST_PLUS: "InvL-",
    SPLIT_SECOND_PLUS | SPLIT_LEFT_FIRST: "InvR+",
    SPLIT_SECOND_PLUS: "InvR-",
}

# 95% Pearson type VI quantiles, orig:unif_ci.py column 0 (n = 1 .. 100).
# Generated from that table; n > 100 reuses the n = 100 value, as the original
# indexes with min(n, 100).
_Q95: tuple[float, ...] = (
    1e100, 20.0, 4.47213595499958, 2.71441761659491, 2.11474252688113,
    1.82056420302608, 1.64754897244207, 1.53412740463439, 1.45421543344895,
    1.39495079396242, 1.34928284767356, 1.31303243448788, 1.28356884211252,
    1.25915510457657, 1.23859896844555, 1.22105530006757, 1.2059085510307,
    1.19269988180426, 1.18108035034134, 1.17077991372278, 1.16158634964154,
    1.15333058538934, 1.1458762736386, 1.13911224808052, 1.13294697233035,
    1.12730439408171, 1.1221208065875, 1.11734244323153, 1.11292361303261,
    1.10882524037588, 1.10501371035276, 1.10145994764715, 1.09813867568295,
    1.09502781619681, 1.09210799914855, 1.0893621600274, 1.08677520690405,
    1.08433374353668, 1.08202583782587, 1.07984082718679, 1.07776915415201,
    1.07580222686676, 1.07393230018905, 1.07215237393049, 1.07045610542357,
    1.06883773411671, 1.06729201631068, 1.06581416848037, 1.06439981789286,
    1.06304495944939, 1.06174591785498, 1.06049931436413, 1.05930203746972,
    1.05815121700027, 1.05704420117208, 1.05597853621039, 1.05495194821025,
    1.05396232695515, 1.05300771145129, 1.05208627696893, 1.05119632341098,
    1.05033626485263, 1.0495046201171, 1.04870000426947, 1.04792112092607,
    1.04716675528959, 1.04643576783136, 1.04572708855154, 1.04503971175671,
    1.04437269130105, 1.04372513624393, 1.04309620688201, 1.04248511111864,
    1.04189110113765, 1.04131347035208, 1.04075155060185, 1.04020470957686,
    1.03967234844462, 1.03915389966395, 1.03864882496761, 1.03815661349914,
    1.03767678009012, 1.03720886366573, 1.03675242576766, 1.03630704918429,
    1.03587233667926, 1.03544790981029, 1.03503340783072, 1.03462848666735,
    1.03423281796815, 1.03384608821464, 1.03346799789362, 1.03309826072387,
    1.03273660293345, 1.0323827625839, 1.03203648893777, 1.03169754186625,
    1.03136569129405, 1.03104071667873, 1.03072240652209,
)


def pearson_q95(n: int) -> float:
    """95% quantile used by orig:unif_ci.py::get_q. n is the cluster size."""
    if n < 1:
        raise ValueError("n must be positive")
    return _Q95[min(n, 100) - 1]


def uniform_ci(xs: list[int], slop: int = 5) -> tuple[tuple[int, int], tuple[int, int]]:
    """orig:unif_ci.py::uniform_ci at confidence 0.95.

    Returns ``(lower_ci, upper_ci)``, each a closed integer interval.
    """
    if not xs:
        raise ValueError("empty sample")
    m = min(xs)
    big = max(xs)
    q = pearson_q95(len(xs))
    lower = (int(math.floor((1 - q) * (big - m) + m)), m + slop)
    upper = (big - slop, int(math.ceil(q * (big - m) + m)))
    return lower, upper


def classify_split_flags(
    first_plus: bool,
    second_plus: bool,
    overlap: bool,
    left_first: bool,
) -> str:
    """orig:splitreads.py flag assembly into SPLIT_TYPES."""
    if first_plus == second_plus:
        flag = (
            int(first_plus) * SPLIT_FIRST_PLUS
            + int(second_plus) * SPLIT_SECOND_PLUS
            + int(overlap) * SPLIT_OVERLAP
        )
    else:
        flag = (
            int(first_plus) * SPLIT_FIRST_PLUS
            + int(second_plus) * SPLIT_SECOND_PLUS
            + int(left_first) * SPLIT_LEFT_FIRST
        )
    return SPLIT_TYPES[flag]


@dataclass(frozen=True)
class Softclip:
    pos: int
    side: str  # "L" or "R"
    mapq: int
    clipped_len: int
    aligned_start: int
    aligned_end: int

    def aligned_overlaps(self, other: "Softclip") -> bool:
        return self.aligned_start < other.aligned_end and other.aligned_start < self.aligned_end


@dataclass
class SoftclipCluster:
    side: str
    pos: int
    n_reads: int
    sum_mapq: int
    aligned_start: int
    aligned_end: int


def cluster_softclips(clips: list[Softclip], merge_distance: int = 5) -> list[SoftclipCluster]:
    """Merge same-side clips within ``merge_distance`` (§4.3a).

    The candidate position is the one whose reads have the highest total MAPQ.
    Ties break toward the smaller position. Input order does not matter.
    """
    by_side: dict[str, list[Softclip]] = {"L": [], "R": []}
    for clip in clips:
        by_side[clip.side].append(clip)
    clusters: list[SoftclipCluster] = []
    for side, group in by_side.items():
        group.sort(key=lambda c: (c.pos, c.mapq, c.aligned_start))
        bucket: list[Softclip] = []

        def flush() -> None:
            if not bucket:
                return
            pos_mapq: dict[int, int] = defaultdict(int)
            for clip in bucket:
                pos_mapq[clip.pos] += clip.mapq
            # Highest total MAPQ, then lowest position.
            pos = min(pos_mapq, key=lambda p: (-pos_mapq[p], p))
            clusters.append(
                SoftclipCluster(
                    side=side,
                    pos=pos,
                    n_reads=len(bucket),
                    sum_mapq=sum(c.mapq for c in bucket),
                    aligned_start=min(c.aligned_start for c in bucket),
                    aligned_end=max(c.aligned_end for c in bucket),
                )
            )

        for clip in group:
            if bucket and clip.pos - bucket[0].pos > merge_distance and clip.pos - bucket[-1].pos > merge_distance:
                # Chain-merge: a new clip joins if it is within merge_distance of any member.
                if min(abs(clip.pos - old.pos) for old in bucket) > merge_distance:
                    flush()
                    bucket = [clip]
                    continue
            if bucket and min(abs(clip.pos - old.pos) for old in bucket) > merge_distance:
                flush()
                bucket = [clip]
            else:
                bucket.append(clip)
        flush()
    clusters.sort(key=lambda c: (c.side, c.pos))
    return clusters


def merge_microhomology(
    clusters: list[SoftclipCluster],
    max_distance: int = 25,
) -> list[tuple[int, int, int, int]]:
    """Opposite-side clusters ≤ ``max_distance`` apart with overlapping aligned bases.

    Returns ``(start, end, n_left, n_right)`` intervals. A cluster is consumed
    by at most one pair (lowest right position wins).
    """
    left = [c for c in clusters if c.side == "L"]
    right = [c for c in clusters if c.side == "R"]
    used_l: set[int] = set()
    used_r: set[int] = set()
    intervals: list[tuple[int, int, int, int]] = []
    for i, lc in enumerate(left):
        for j, rc in enumerate(right):
            if i in used_l or j in used_r:
                continue
            if abs(lc.pos - rc.pos) > max_distance:
                continue
            if lc.aligned_start < rc.aligned_end and rc.aligned_start < lc.aligned_end:
                a, b = sorted((lc.pos, rc.pos))
                intervals.append((a, b if b > a else a + 1, lc.n_reads, rc.n_reads))
                used_l.add(i)
                used_r.add(j)
    # Unpaired clusters become 1-bp intervals.
    for i, lc in enumerate(left):
        if i not in used_l:
            intervals.append((lc.pos, lc.pos + 1, lc.n_reads, 0))
    for j, rc in enumerate(right):
        if j not in used_r:
            intervals.append((rc.pos, rc.pos + 1, 0, rc.n_reads))
    intervals.sort()
    return intervals


@dataclass(frozen=True)
class Drp:
    """One discordant pair. ``a`` is the left 3'-most position, ``b`` the right."""

    type: str  # Del, Ins, Dup, InvL, InvR
    a: int
    b: int
    insert: int
    index: int = 0


def compatible(p: Drp, q: Drp, max_distance: float, ins_slop: int = 20) -> bool:
    """§4.3c / orig compatibility. Insertions allow ``h = ins_slop`` of overlap."""
    if p.type != q.type:
        return False
    if abs(p.a - q.a) > max_distance or abs(p.b - q.b) > max_distance:
        return False
    if p.type == "Ins":
        # orig:is_ins_compatible — close ends and overlap ≤ max_ins_pair_slop.
        overlap = max(0, max(p.a, q.a) - min(p.b, q.b))
        return overlap <= ins_slop
    return max(p.a, q.a) < min(p.b, q.b)


def ins_cluster_ok(cluster: list[Drp], cluster_slop: int = 25) -> bool:
    """orig:is_ins_cluster_compatible — cluster-wide overlap ≤ max_ins_cluster_slop."""
    if not cluster or cluster[0].type != "Ins":
        return True
    overlap = max(0, max(p.a for p in cluster) - min(p.b for p in cluster))
    return overlap <= cluster_slop


def _bron_kerbosch(
    adj: list[set[int]],
    r: list[int],
    p: set[int],
    x: set[int],
    best_size: list[int],
    out: list[list[int]],
    node_cap: int,
    expansions: list[int],
    expansion_cap: int,
) -> bool:
    """Returns False if the expansion cap was hit."""
    if expansions[0] > expansion_cap:
        return False
    expansions[0] += 1
    if not p and not x:
        if len(r) > best_size[0]:
            best_size[0] = len(r)
            out.clear()
            out.append(sorted(r))
        elif len(r) == best_size[0] and len(r) >= 2:
            out.append(sorted(r))
        return True
    # Tomita pivoting: pick the vertex in P∪X with the most neighbors in P.
    universe = p | x
    if not universe:
        return True
    pivot = max(universe, key=lambda v: len(adj[v] & p))
    candidates = sorted(p - adj[pivot])
    for v in candidates:
        if expansions[0] > expansion_cap:
            return False
        if not _bron_kerbosch(
            adj,
            r + [v],
            p & adj[v],
            x & adj[v],
            best_size,
            out,
            node_cap,
            expansions,
            expansion_cap,
        ):
            return False
        p.remove(v)
        x.add(v)
    return True


def largest_cliques(
    pairs: list[Drp],
    max_distance: float,
    ins_slop: int = 20,
    node_cap: int = 40,
    expansion_cap: int = 100_000,
) -> tuple[list[list[Drp]], bool]:
    """Largest clique(s) of one connected component.

    Tie-break is the lexicographically smallest index set, and every maximum
    clique is returned (orig:pecluster.py::cluster_handle_component uses
    igraph.largest_cliques). ``hit_cap`` is true when the search was truncated.
    """
    n = len(pairs)
    if n > node_cap:
        # Huge component: keep it only when every pair is compatible (orig).
        ok = all(
            compatible(pairs[i], pairs[j], max_distance, ins_slop)
            for i in range(n)
            for j in range(i + 1, n)
        )
        return ([list(pairs)] if ok else []), True
    adj: list[set[int]] = [set() for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if compatible(pairs[i], pairs[j], max_distance, ins_slop):
                adj[i].add(j)
                adj[j].add(i)
    best_size = [0]
    out: list[list[int]] = []
    ok = _bron_kerbosch(
        adj, [], set(range(n)), set(), best_size, out, node_cap, [0], expansion_cap
    )
    cliques = [sorted((pairs[i] for i in clique), key=lambda d: d.index) for clique in out]
    cliques.sort(key=lambda cl: tuple(d.index for d in cl))
    # Singletons are ignored by the caller.
    cliques = [c for c in cliques if len(c) >= 2]
    return cliques, not ok


def cluster_drps(
    pairs: list[Drp],
    max_distance: float,
    ins_slop: int = 20,
    max_cluster_size: int = 80,
    ins_cluster_slop: int = 25,
) -> tuple[list[list[Drp]], list[str]]:
    """Connected components, then the largest clique in each (§4.3c).

    Components are grown in order of ``a`` so the partition does not depend on
    input order. Logs (as returned strings) when a size cap is hit. Insertion
    clusters are filtered by ``ins_cluster_ok`` (orig:is_ins_cluster_compatible).
    """
    ordered = sorted(pairs, key=lambda d: (d.a, d.b, d.insert, d.index))
    logs: list[str] = []
    n = len(ordered)
    seen = [False] * n
    clusters: list[list[Drp]] = []

    def neighbors(i: int) -> list[int]:
        # Pairs are sorted by ``a``. Compatibility requires |a_i - a_j| ≤ max_distance,
        # so only that window can match.
        out: list[int] = []
        ai = ordered[i].a
        j = i - 1
        while j >= 0 and not (ai - ordered[j].a > max_distance):
            if compatible(ordered[i], ordered[j], max_distance, ins_slop):
                out.append(j)
            j -= 1
        j = i + 1
        while j < n and not (ordered[j].a - ai > max_distance):
            if compatible(ordered[i], ordered[j], max_distance, ins_slop):
                out.append(j)
            j += 1
        return out

    for i in range(n):
        if seen[i]:
            continue
        # BFS component.
        stack = [i]
        seen[i] = True
        comp = [i]
        while stack:
            u = stack.pop()
            for v in neighbors(u):
                if not seen[v]:
                    seen[v] = True
                    stack.append(v)
                    comp.append(v)
        members = [ordered[k] for k in comp]
        if len(members) < 2:
            continue
        if len(members) > max_cluster_size:
            fully = all(
                compatible(members[a], members[b], max_distance, ins_slop)
                for a in range(len(members))
                for b in range(a + 1, len(members))
            )
            logs.append(f"component size {len(members)} exceeded cap {max_cluster_size}")
            if fully:
                clusters.append(sorted(members, key=lambda d: d.index))
            continue
        cliques, hit = largest_cliques(members, max_distance, ins_slop, node_cap=max_cluster_size)
        if hit:
            logs.append(f"clique search cap hit on component of {len(members)}")
        clusters.extend(cliques)
    if pairs and pairs[0].type == "Ins":
        clusters = [c for c in clusters if ins_cluster_ok(c, ins_cluster_slop)]
    clusters.sort(key=lambda cl: tuple(d.index for d in cl))
    return clusters, logs


def score_cluster(cluster: list[Drp], median: float, scale: float) -> float:
    """Log cluster score (§4.3c).

    Inversions and duplications score by size. Deletions use Methods eq. 2 /
    ``pecluster.lr_del``: ``D = min(mean(z) − μ, min b − max a)``.
    Insertions use the symmetric shift ``max(0, m - mean z)``.
    """
    if scale <= 0:
        scale = 1.0
    kind = cluster[0].type
    if kind in ("InvL", "InvR", "Dup"):
        return float(len(cluster))
    zs = [p.insert for p in cluster]
    var = 2.0 * scale * scale
    if kind == "Del":
        mean_z = sum(zs) / len(zs)
        max_size = min(p.b for p in cluster) - max(p.a for p in cluster)
        d_size = min(mean_z - median, float(max_size))
        d_size = max(0.0, d_size)
        return sum(((z - median) ** 2 - (z - median - d_size) ** 2) / var for z in zs)
    # Insertion.
    overlap = max(0, max(p.a for p in cluster) - min(p.b for p in cluster))
    ins = max(float(overlap), median - (sum(zs) / len(zs)))
    return sum(((z - median) ** 2 - (z - median + ins) ** 2) / var for z in zs)


def fdr_keep(
    scores: list[float],
    null_scores: list[float],
    null_reps: int,
    alpha: float = 0.1,
) -> list[int]:
    """Indices of clusters with score ≥ s* at target FDR (§4.3c).

    ``π0 = min(1, M0/M)`` with ``M0`` the mean null-cluster count per replicate.
    """
    m = len(scores)
    if m == 0:
        return []
    m0 = len(null_scores) / max(null_reps, 1)
    if m0 == 0:
        return list(range(m))
    pi0 = min(1.0, m0 / m)
    order = sorted(range(m), key=lambda i: (scores[i], i))
    null_sorted = sorted(null_scores)

    def p_null_exceeds(x: float) -> float:
        # Fraction of null scores >= x. bisect_left on the sorted nulls.
        lo, hi = 0, len(null_sorted)
        while lo < hi:
            mid = (lo + hi) // 2
            if null_sorted[mid] < x:
                lo = mid + 1
            else:
                hi = mid
        return 1.0 - lo / len(null_sorted)

    cutoff_rank = m  # reject nothing
    for rank, idx in enumerate(order):
        num_rej = m - rank
        est = m * pi0 * p_null_exceeds(scores[idx]) / num_rej
        if est <= alpha:
            cutoff_rank = rank
            break
    return [order[i] for i in range(cutoff_rank, m)]


def relocate_drps(
    pairs: list[Drp],
    chrom_len: int,
    rng: random.Random,
    max_insert: float | None = None,
) -> list[Drp]:
    """Uniform relocation along the non-gap chromosome, keeping insert and type."""
    moved: list[Drp] = []
    for i, pair in enumerate(pairs):
        if max_insert is not None and pair.insert > max_insert:
            continue
        span = pair.b - pair.a
        if abs(span) >= chrom_len:
            continue
        lo = max(0, -span)
        hi = chrom_len - max(0, span)
        if hi <= lo:
            continue
        # random.Random.randrange is deterministic for a fixed seed.
        new_a = rng.randrange(lo, hi)
        moved.append(Drp(pair.type, new_a, new_a + span, pair.insert, index=i))
    return moved


@dataclass
class Breakpoint:
    start: int
    end: int  # half-open
    n_softclip: int = 0
    n_split: int = 0
    n_drp: int = 0
    drp_types: tuple[str, ...] = ()
    n_clip_left: int = 0
    n_clip_right: int = 0
    # Reference base supported by soft clips, when this interval contains any.
    clip_pos: int | None = None
    # DRP cluster that produced this interval. Mate ends share one id.
    cluster_ids: tuple[int, ...] = ()

    @property
    def point(self) -> int:
        return self.start

    def overlaps(self, other: "Breakpoint") -> bool:
        return self.start < other.end and other.start < self.end

    def _kept_clip(self, other: "Breakpoint") -> int | None:
        """Clip coordinate that survives a merge. More clips win; a tie keeps the rightmost."""
        if self.clip_pos is None:
            return other.clip_pos
        if other.clip_pos is None:
            return self.clip_pos
        self_n = self.n_clip_left + self.n_clip_right
        other_n = other.n_clip_left + other.n_clip_right
        if other_n > self_n:
            return other.clip_pos
        if self_n > other_n:
            return self.clip_pos
        return max(self.clip_pos, other.clip_pos)

    def merged(self, other: "Breakpoint") -> "Breakpoint":
        return Breakpoint(
            start=min(self.start, other.start),
            end=max(self.end, other.end),
            n_softclip=self.n_softclip + other.n_softclip,
            n_split=self.n_split + other.n_split,
            n_drp=self.n_drp + other.n_drp,
            drp_types=tuple(sorted(set(self.drp_types + other.drp_types))),
            n_clip_left=self.n_clip_left + other.n_clip_left,
            n_clip_right=self.n_clip_right + other.n_clip_right,
            clip_pos=self._kept_clip(other),
            cluster_ids=tuple(sorted(set(self.cluster_ids + other.cluster_ids))),
        )


def cluster_to_intervals(
    cluster: list[Drp], cluster_id: int | None = None
) -> tuple[Breakpoint, Breakpoint]:
    """Two breakpoint intervals from a DRP cluster (orig:pecluster.py::cluster_to_ci)."""
    p1 = [p.a for p in cluster]
    p2 = [p.b for p in cluster]
    lower_1, upper_1 = uniform_ci(p1)
    lower_2, upper_2 = uniform_ci(p2)
    kind = cluster[0].type
    if kind in ("Del", "Ins"):
        left, right = upper_1, lower_2
    elif kind == "Dup":
        left, right = lower_1, upper_2
    elif kind == "InvL":
        left, right = upper_1, upper_2
    else:
        left, right = lower_1, lower_2
    ids = (cluster_id,) if cluster_id is not None else ()

    def as_bp(interval: tuple[int, int]) -> Breakpoint:
        a, b = interval
        if b < a:
            a, b = b, a
        # Half-open, at least 1 bp.
        return Breakpoint(
            a,
            b + 1 if b == a else b,
            n_drp=len(cluster),
            drp_types=(kind,),
            cluster_ids=ids,
        )

    return as_bp(left), as_bp(right)


def clamp_breakpoint(bp: Breakpoint, chrom_len: int) -> Breakpoint | None:
    """Drop a confidence interval that misses the chromosome, and clip one that crosses an end."""
    if chrom_len <= 0:
        return None
    start = min(max(bp.start, 0), chrom_len)
    end = min(max(bp.end, 0), chrom_len)
    if end <= start:
        return None
    return Breakpoint(
        start=start,
        end=end,
        n_softclip=bp.n_softclip,
        n_split=bp.n_split,
        n_drp=bp.n_drp,
        drp_types=bp.drp_types,
        n_clip_left=bp.n_clip_left,
        n_clip_right=bp.n_clip_right,
        clip_pos=bp.clip_pos if bp.clip_pos is not None and start <= bp.clip_pos <= end else None,
        cluster_ids=bp.cluster_ids,
    )


def _bp_precedence(bp: Breakpoint) -> int:
    """Original precedence: clip/split = 2, DRP-only = 1 (breakpoint_merge.Breakpoint.precedence)."""
    if bp.n_split > 0 or bp.n_softclip > 0 or bp.n_clip_left > 0 or bp.n_clip_right > 0:
        return 2
    if bp.clip_pos is not None:
        return 2
    if bp.n_drp > 0:
        return 1
    return 0


def _absorb_support(keeper: Breakpoint, other: Breakpoint) -> Breakpoint:
    """Combine support; keep ``keeper`` coordinates (orig ``Breakpoint.__add__`` when ranks differ)."""
    return Breakpoint(
        start=keeper.start,
        end=keeper.end,
        n_softclip=keeper.n_softclip + other.n_softclip,
        n_split=keeper.n_split + other.n_split,
        n_drp=keeper.n_drp + other.n_drp,
        drp_types=tuple(sorted(set(keeper.drp_types + other.drp_types))),
        n_clip_left=keeper.n_clip_left + other.n_clip_left,
        n_clip_right=keeper.n_clip_right + other.n_clip_right,
        clip_pos=keeper.clip_pos if keeper.clip_pos is not None else other.clip_pos,
        cluster_ids=tuple(sorted(set(keeper.cluster_ids + other.cluster_ids))),
    )


def _interval_dist(a: tuple[int, int], b: tuple[int, int]) -> int:
    """Gap between half-open intervals; 0 if they touch or overlap (orig closed-interval dist)."""
    return max(a[0] - b[1], b[0] - a[1], 0)


def _merge_nearby_intervals(
    objects: dict[tuple[int, int], list[Breakpoint]],
    mergefun,
    max_distance: int = 0,
) -> dict[tuple[int, int], Breakpoint]:
    """Port of ``helper.merge_nearby`` for interval keys."""
    if not objects:
        return {}
    locations = sorted(objects.keys())
    merged: dict[tuple[int, int], Breakpoint] = {}
    cur_locs: list[tuple[int, int]] = []
    cur_objects: list[Breakpoint] = []
    for loc in locations:
        if not cur_locs:
            cur_locs = [loc]
            cur_objects = list(objects[loc])
        elif any(_interval_dist(loc, existing) <= max_distance for existing in cur_locs):
            cur_locs.append(loc)
            cur_objects.extend(objects[loc])
        else:
            merged.update(mergefun(cur_locs, cur_objects))
            cur_locs = [loc]
            cur_objects = list(objects[loc])
    merged.update(mergefun(cur_locs, cur_objects))
    return merged


def _bp_mergefun_same_level(_locs: list[tuple[int, int]], bps: list[Breakpoint]):
    """Union merge for same-precedence BPs (orig ``bp_mergefun`` / ``sum``)."""
    combined = bps[0]
    for bp in bps[1:]:
        combined = combined.merged(bp)
    return {(combined.start, combined.end): combined}


def _bp_mergefun_precedence(
    _locs: list[tuple[int, int]],
    bps: list[Breakpoint],
    max_distance: int = 0,
) -> dict[tuple[int, int], Breakpoint]:
    """Port of ``breakpoint_merge.bp_mergefun_precedence``."""
    merged_above: list[Breakpoint] = []
    merged: dict[tuple[int, int], Breakpoint] = {}
    prec_list = sorted({_bp_precedence(bp) for bp in bps}, reverse=True)
    for i, prec in enumerate(prec_list):
        p_dict: dict[tuple[int, int], list[Breakpoint]] = defaultdict(list)
        for bp in bps:
            if _bp_precedence(bp) == prec and bp not in merged_above:
                p_dict[(bp.start, bp.end)].append(bp)
        if p_dict:
            p_mrg = _merge_nearby_intervals(p_dict, _bp_mergefun_same_level, max_distance)
            merged.update(p_mrg)
        if i + 1 < len(prec_list):
            lower_prec = prec_list[i + 1]
            lower_bps = [bp for bp in bps if _bp_precedence(bp) == lower_prec]
            for lower_bp in lower_bps:
                for loc, bp in list(merged.items()):
                    if _interval_dist((lower_bp.start, lower_bp.end), (bp.start, bp.end)) <= max_distance:
                        merged[loc] = _absorb_support(bp, lower_bp)
                        merged_above.append(lower_bp)
    return merged


def merge_breakpoints(
    pieces: list[Breakpoint],
    min_split_or_clip: int = 2,
) -> list[Breakpoint]:
    """Merge overlapping intervals with original clip/split precedence (§4.3d).

    Keeps clip/split coordinates when a PE interval touches them
    (``orig:breakpoint_merge``). Paper Methods say "merge all overlapping";
    see ``docs/PARITY.md``.

    Keep a merged interval if it has a DRP cluster or at least
    ``min_split_or_clip`` split/clipped reads.
    """
    if not pieces:
        return []
    # Collapse identical intervals, then connected-component precedence merge.
    by_iv: dict[tuple[int, int], Breakpoint] = {}
    for bp in pieces:
        key = (bp.start, bp.end)
        by_iv[key] = by_iv[key].merged(bp) if key in by_iv else bp
    objects = {key: [bp] for key, bp in by_iv.items()}
    merged_map = _merge_nearby_intervals(objects, _bp_mergefun_precedence, max_distance=0)
    merged_list = sorted(merged_map.values(), key=lambda b: (b.start, b.end))
    return [
        bp
        for bp in merged_list
        if bp.n_drp > 0 or (bp.n_softclip + bp.n_split) >= min_split_or_clip
    ]
