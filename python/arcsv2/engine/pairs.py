"""Columnar pair storage (M2): numpy struct-of-arrays for model-eligible pairs.

Target layout keeps under ~100 bytes/pair. Thin PairView objects are created
only where call-site code still needs attribute access. The Python list of
PairRecord objects remains available when CallConfig.columnar_pairs is False.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterator, Sequence

import numpy as np

# Flag bit layout in flags uint8
_FLAG_REVERSE1 = 1 << 0
_FLAG_REVERSE2 = 1 << 1
_FLAG_HANGING = 1 << 2
_FLAG_DISTANT = 1 << 3


@dataclass(slots=True)
class PairView:
    """Thin attribute view over one row of a Pairs SoA. Matches PairRecord fields."""

    start1: int
    end1: int
    start2: int
    end2: int
    reverse1: bool
    reverse2: bool
    mapq1: int
    mapq2: int
    insert: int
    lib_id: int = 0
    hanging: bool = False
    qname_hash: int = 0
    distant: bool = False
    mate_len: int = 0


class Pairs:
    """Struct-of-arrays for pair evidence.

    Columns (bytes/row when packed):
      start1/end1/start2/end2: int64 x4 = 32
      insert: int32 = 4
      mate_len: int32 = 4
      qname_hash: int64 = 8
      mapq1/mapq2: uint8 x2 = 2
      lib_id: uint16 = 2   # samio hashes RG/LB into 0..999
      flags: uint8 = 1
      Total payload ~53 B + numpy overhead; well under 100 B/pair at scale.
    """

    __slots__ = (
        "start1",
        "end1",
        "start2",
        "end2",
        "insert",
        "mate_len",
        "qname_hash",
        "mapq1",
        "mapq2",
        "lib_id",
        "flags",
        "_n",
        "_cap",
    )

    def __init__(self, capacity: int = 0) -> None:
        cap = max(int(capacity), 0)
        self._n = 0
        self._cap = cap
        self.start1 = np.empty(cap, dtype=np.int64)
        self.end1 = np.empty(cap, dtype=np.int64)
        self.start2 = np.empty(cap, dtype=np.int64)
        self.end2 = np.empty(cap, dtype=np.int64)
        self.insert = np.empty(cap, dtype=np.int32)
        self.mate_len = np.empty(cap, dtype=np.int32)
        self.qname_hash = np.empty(cap, dtype=np.int64)
        self.mapq1 = np.empty(cap, dtype=np.uint8)
        self.mapq2 = np.empty(cap, dtype=np.uint8)
        self.lib_id = np.empty(cap, dtype=np.uint16)
        self.flags = np.empty(cap, dtype=np.uint8)

    def __len__(self) -> int:
        return self._n

    def release(self) -> None:
        """Drop column buffers (M4: after rust build_graph has copied them)."""
        empty_i64 = np.empty(0, dtype=np.int64)
        empty_i32 = np.empty(0, dtype=np.int32)
        empty_u8 = np.empty(0, dtype=np.uint8)
        empty_u16 = np.empty(0, dtype=np.uint16)
        self.start1 = empty_i64
        self.end1 = empty_i64
        self.start2 = empty_i64
        self.end2 = empty_i64
        self.insert = empty_i32
        self.mate_len = empty_i32
        self.qname_hash = empty_i64
        self.mapq1 = empty_u8
        self.mapq2 = empty_u8
        self.lib_id = empty_u16
        self.flags = empty_u8
        self._n = 0
        self._cap = 0

    def bytes_per_pair(self) -> float:
        if self._n == 0:
            return 0.0
        total = 0
        for name in self.__slots__:
            if name.startswith("_"):
                continue
            arr = getattr(self, name)
            total += arr[: self._n].nbytes
        return total / self._n

    def _grow(self, need: int) -> None:
        if need <= self._cap:
            return
        new_cap = max(need, max(1024, self._cap * 2))
        for name in (
            "start1",
            "end1",
            "start2",
            "end2",
            "insert",
            "mate_len",
            "qname_hash",
            "mapq1",
            "mapq2",
            "lib_id",
            "flags",
        ):
            old = getattr(self, name)
            fresh = np.empty(new_cap, dtype=old.dtype)
            if self._n:
                fresh[: self._n] = old[: self._n]
            setattr(self, name, fresh)
        self._cap = new_cap

    def append(
        self,
        *,
        start1: int,
        end1: int,
        start2: int,
        end2: int,
        reverse1: bool,
        reverse2: bool,
        mapq1: int,
        mapq2: int,
        insert: int,
        lib_id: int = 0,
        hanging: bool = False,
        qname_hash: int = 0,
        distant: bool = False,
        mate_len: int = 0,
    ) -> int:
        i = self._n
        self._grow(i + 1)
        self.start1[i] = start1
        self.end1[i] = end1
        self.start2[i] = start2
        self.end2[i] = end2
        self.insert[i] = insert
        self.mate_len[i] = mate_len
        self.qname_hash[i] = qname_hash
        self.mapq1[i] = mapq1 & 0xFF
        self.mapq2[i] = mapq2 & 0xFF
        self.lib_id[i] = lib_id & 0xFFFF
        flags = 0
        if reverse1:
            flags |= _FLAG_REVERSE1
        if reverse2:
            flags |= _FLAG_REVERSE2
        if hanging:
            flags |= _FLAG_HANGING
        if distant:
            flags |= _FLAG_DISTANT
        self.flags[i] = flags
        self._n = i + 1
        return i

    def trim(self) -> None:
        """Drop unused capacity after load."""
        n = self._n
        if n == self._cap:
            return
        for name in (
            "start1",
            "end1",
            "start2",
            "end2",
            "insert",
            "mate_len",
            "qname_hash",
            "mapq1",
            "mapq2",
            "lib_id",
            "flags",
        ):
            setattr(self, name, getattr(self, name)[:n].copy())
        self._cap = n

    @classmethod
    def from_arrays(
        cls,
        *,
        start1,
        end1,
        start2,
        end2,
        insert,
        mate_len,
        qname_hash,
        mapq1,
        mapq2,
        lib_id,
        flags,
    ) -> Pairs:
        """Wrap existing columnar arrays (Rust / zero-copy load path)."""
        p = cls(capacity=0)
        p.start1 = np.asarray(start1, dtype=np.int64)
        p.end1 = np.asarray(end1, dtype=np.int64)
        p.start2 = np.asarray(start2, dtype=np.int64)
        p.end2 = np.asarray(end2, dtype=np.int64)
        p.insert = np.asarray(insert, dtype=np.int32)
        p.mate_len = np.asarray(mate_len, dtype=np.int32)
        p.qname_hash = np.asarray(qname_hash, dtype=np.int64)
        p.mapq1 = np.asarray(mapq1, dtype=np.uint8)
        p.mapq2 = np.asarray(mapq2, dtype=np.uint8)
        p.lib_id = np.asarray(lib_id, dtype=np.uint16)
        p.flags = np.asarray(flags, dtype=np.uint8)
        p._n = int(p.start1.shape[0])
        p._cap = p._n
        return p

    def view_at(self, i: int) -> PairView:
        flags = int(self.flags[i])
        return PairView(
            start1=int(self.start1[i]),
            end1=int(self.end1[i]),
            start2=int(self.start2[i]),
            end2=int(self.end2[i]),
            reverse1=bool(flags & _FLAG_REVERSE1),
            reverse2=bool(flags & _FLAG_REVERSE2),
            mapq1=int(self.mapq1[i]),
            mapq2=int(self.mapq2[i]),
            insert=int(self.insert[i]),
            lib_id=int(self.lib_id[i]),
            hanging=bool(flags & _FLAG_HANGING),
            qname_hash=int(self.qname_hash[i]),
            distant=bool(flags & _FLAG_DISTANT),
            mate_len=int(self.mate_len[i]),
        )

    def __getitem__(self, i: int) -> PairView:
        if i < 0:
            i += self._n
        if i < 0 or i >= self._n:
            raise IndexError(i)
        return self.view_at(i)

    def __iter__(self) -> Iterator[PairView]:
        for i in range(self._n):
            yield self.view_at(i)

    def to_records(self) -> list:
        """Materialize PairRecord list (fallback / tests)."""
        from arcsv2.engine.caller import PairRecord

        out: list[PairRecord] = []
        for i in range(self._n):
            v = self.view_at(i)
            out.append(
                PairRecord(
                    start1=v.start1,
                    end1=v.end1,
                    start2=v.start2,
                    end2=v.end2,
                    reverse1=v.reverse1,
                    reverse2=v.reverse2,
                    mapq1=v.mapq1,
                    mapq2=v.mapq2,
                    insert=v.insert,
                    lib_id=v.lib_id,
                    hanging=v.hanging,
                    qname_hash=v.qname_hash,
                    distant=v.distant,
                    mate_len=v.mate_len,
                )
            )
        return out

    @classmethod
    def from_records(cls, records: Sequence) -> Pairs:
        p = cls(capacity=len(records))
        for r in records:
            p.append(
                start1=r.start1,
                end1=r.end1,
                start2=r.start2,
                end2=r.end2,
                reverse1=r.reverse1,
                reverse2=r.reverse2,
                mapq1=r.mapq1,
                mapq2=r.mapq2,
                insert=r.insert,
                lib_id=r.lib_id,
                hanging=r.hanging,
                qname_hash=r.qname_hash,
                distant=r.distant,
                mate_len=r.mate_len,
            )
        p.trim()
        return p
