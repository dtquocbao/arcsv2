#!/usr/bin/env python3
"""Compare two ml-train manifests.

Floats may differ by at most 1e-9 absolute and 1e-9 relative. Thresholds, keys,
grouping, and chromosome-holdout structure must match exactly.

Usage: python scripts/compare_manifest.py old.json new.json
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ABS_TOL = 1e-9
REL_TOL = 1e-9


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _is_threshold_key(key: str) -> bool:
    return key == "threshold" or key == "thresholds"


def _rel(a: float, b: float) -> float:
    scale = max(abs(a), abs(b))
    if scale == 0.0:
        return 0.0
    return abs(a - b) / scale


def _walk(old, new, path: str, errors: list[str], stats: dict) -> None:
    if type(old) is not type(new) and not (isinstance(old, (int, float)) and isinstance(new, (int, float)) and not isinstance(old, bool) and not isinstance(new, bool)):
        errors.append(f"{path}: type {type(old).__name__} vs {type(new).__name__}")
        return
    if isinstance(old, dict):
        old_keys = set(old)
        new_keys = set(new)
        if old_keys != new_keys:
            missing = sorted(old_keys - new_keys)
            extra = sorted(new_keys - old_keys)
            errors.append(f"{path}: keys differ missing={missing} extra={extra}")
        for key in sorted(old_keys & new_keys):
            child = f"{path}.{key}" if path else key
            if _is_threshold_key(key) and not isinstance(old[key], dict):
                if old[key] != new[key]:
                    errors.append(f"{child}: threshold {old[key]!r} != {new[key]!r}")
                continue
            _walk(old[key], new[key], child, errors, stats)
        return
    if isinstance(old, list):
        if len(old) != len(new):
            errors.append(f"{path}: length {len(old)} != {len(new)}")
            return
        for i, (left, right) in enumerate(zip(old, new)):
            _walk(left, right, f"{path}[{i}]", errors, stats)
        return
    if isinstance(old, bool) or isinstance(new, bool):
        if old != new:
            errors.append(f"{path}: {old!r} != {new!r}")
        return
    if isinstance(old, (int, float)) and isinstance(new, (int, float)):
        # Threshold objects are handled above. A bare number here is a model float
        # unless the path's last key is a threshold (lists of thresholds).
        leaf = path.rsplit(".", 1)[-1]
        leaf = leaf.split("[", 1)[0]
        if _is_threshold_key(leaf) or ".thresholds" in f".{path}" or path.startswith("thresholds"):
            if old != new:
                errors.append(f"{path}: threshold {old!r} != {new!r}")
            return
        a = float(old)
        b = float(new)
        if a != a and b != b:
            return
        diff = abs(a - b)
        rel = _rel(a, b)
        stats["max_abs"] = max(stats["max_abs"], diff)
        stats["max_rel"] = max(stats["max_rel"], rel)
        stats["n"] += 1
        if diff > ABS_TOL or rel > REL_TOL:
            errors.append(f"{path}: abs {diff:.3e} rel {rel:.3e} ({a!r} vs {b!r})")
        return
    if old != new:
        errors.append(f"{path}: {old!r} != {new!r}")


def compare(old: dict, new: dict) -> tuple[list[str], dict]:
    errors: list[str] = []
    stats = {"max_abs": 0.0, "max_rel": 0.0, "n": 0}
    _walk(old, new, "", errors, stats)
    return errors, stats


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print("usage: python scripts/compare_manifest.py old.json new.json", file=sys.stderr)
        return 2
    old = _load(Path(argv[1]))
    new = _load(Path(argv[2]))
    errors, stats = compare(old, new)
    print(f"max_abs {stats['max_abs']:.3e}")
    print(f"max_rel {stats['max_rel']:.3e}")
    print(f"floats {stats['n']}")
    if errors:
        print("FAIL")
        for line in errors[:40]:
            print(line)
        if len(errors) > 40:
            print(f"... {len(errors) - 40} more")
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
