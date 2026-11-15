"""Real intervals with the openness rules used by orig:pyinter.

Two pieces become one when their union has no gap. ``[0, 100)`` and
``[100, 200)`` join into ``[0, 200)``. ``(5, 10)`` and ``(10, 15)`` stay
apart, because the shared endpoint is in neither piece. That is the
arithmetic behind ``get_gap_overlap_positions``.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Iv:
    lo: int
    hi: int
    lo_in: bool
    hi_in: bool

    def empty(self) -> bool:
        if self.lo > self.hi:
            return True
        if self.lo == self.hi:
            return not (self.lo_in and self.hi_in)
        return False

    def __iter_endpoints__(self) -> tuple[int, int, bool, bool]:
        return self.lo, self.hi, self.lo_in, self.hi_in


def open_iv(lo: int, hi: int) -> Iv:
    return Iv(lo, hi, False, False)


def closed_iv(lo: int, hi: int) -> Iv:
    return Iv(lo, hi, True, True)


def closedopen_iv(lo: int, hi: int) -> Iv:
    return Iv(lo, hi, True, False)


def _left_of(a: Iv, b: Iv) -> bool:
    if a.hi < b.lo:
        return True
    return a.hi == b.lo and not (a.hi_in and b.lo_in)


def overlaps(a: Iv, b: Iv) -> bool:
    return not _left_of(a, b) and not _left_of(b, a)


def _abut(a: Iv, b: Iv) -> bool:
    """True when the pieces touch and at least one of them owns the shared point."""
    if a.hi == b.lo and (a.hi_in or b.lo_in):
        return True
    return b.hi == a.lo and (b.hi_in or a.lo_in)


def joinable(a: Iv, b: Iv) -> bool:
    return overlaps(a, b) or _abut(a, b)


def _later_start(a: Iv, b: Iv) -> tuple[int, bool]:
    if a.lo > b.lo or (a.lo == b.lo and not a.lo_in and b.lo_in):
        return a.lo, a.lo_in
    if b.lo > a.lo or (a.lo == b.lo and not b.lo_in and a.lo_in):
        return b.lo, b.lo_in
    return a.lo, a.lo_in or b.lo_in


def _earlier_end(a: Iv, b: Iv) -> tuple[int, bool]:
    if a.hi < b.hi or (a.hi == b.hi and not a.hi_in and b.hi_in):
        return a.hi, a.hi_in
    if b.hi < a.hi or (a.hi == b.hi and not b.hi_in and a.hi_in):
        return b.hi, b.hi_in
    return a.hi, a.hi_in or b.hi_in


def intersection(a: Iv, b: Iv) -> Iv | None:
    if not overlaps(a, b):
        return None
    lo, lo_in = _later_start(a, b)
    hi, hi_in = _earlier_end(a, b)
    iv = Iv(lo, hi, lo_in, hi_in)
    return None if iv.empty() else iv


def union(a: Iv, b: Iv) -> Iv:
    if a.lo < b.lo or (a.lo == b.lo and a.lo_in and not b.lo_in):
        lo, lo_in = a.lo, a.lo_in
    elif b.lo < a.lo or (a.lo == b.lo and b.lo_in and not a.lo_in):
        lo, lo_in = b.lo, b.lo_in
    else:
        lo, lo_in = a.lo, a.lo_in or b.lo_in
    if a.hi > b.hi or (a.hi == b.hi and a.hi_in and not b.hi_in):
        hi, hi_in = a.hi, a.hi_in
    elif b.hi > a.hi or (a.hi == b.hi and b.hi_in and not a.hi_in):
        hi, hi_in = b.hi, b.hi_in
    else:
        hi, hi_in = a.hi, a.hi_in or b.hi_in
    return Iv(lo, hi, lo_in, hi_in)


def subtract(a: Iv, b: Iv) -> list[Iv]:
    """``a \\ b``. Abutting but non-overlapping ``b`` leaves ``a`` unchanged."""
    if a.empty() or not overlaps(a, b):
        return [] if a.empty() else [a]
    shared = intersection(a, b)
    if shared is None:
        return [a]
    pieces: list[Iv] = []
    if a.lo < shared.lo or (a.lo == shared.lo and a.lo_in and not shared.lo_in):
        left = Iv(a.lo, shared.lo, a.lo_in, not shared.lo_in)
        if not left.empty():
            pieces.append(left)
    if a.hi > shared.hi or (a.hi == shared.hi and a.hi_in and not shared.hi_in):
        right = Iv(shared.hi, a.hi, not shared.hi_in, a.hi_in)
        if not right.empty():
            pieces.append(right)
    return pieces


class IntervalSet:
    def __init__(self, items: list[Iv] | None = None) -> None:
        self.items: list[Iv] = []
        for item in items or []:
            self.add(item)

    def add(self, iv: Iv) -> None:
        if iv.empty():
            return
        pending = sorted([*self.items, iv], key=lambda item: (item.lo, not item.lo_in, item.hi))
        merged: list[Iv] = []
        for item in pending:
            if merged and joinable(merged[-1], item):
                merged[-1] = union(merged[-1], item)
            else:
                merged.append(item)
        self.items = merged

    def intersection(self, other: IntervalSet) -> IntervalSet:
        out = IntervalSet()
        for a in self.items:
            for b in other.items:
                shared = intersection(a, b)
                if shared is not None:
                    out.add(shared)
        return out

    def difference(self, other: IntervalSet) -> IntervalSet:
        out = IntervalSet()
        for a in self.items:
            pieces = [a]
            for b in other.items:
                nxt: list[Iv] = []
                for piece in pieces:
                    nxt.extend(subtract(piece, b))
                pieces = nxt
            for piece in pieces:
                out.add(piece)
        return out

    def __iter__(self):
        return iter(self.items)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, IntervalSet):
            return NotImplemented
        return self.items == other.items

    def __repr__(self) -> str:
        return f"IntervalSet({self.items!r})"
