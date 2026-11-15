"""Simulator checksum and oracle helpers."""

import hashlib
from pathlib import Path

from scripts.simulate import simulate


def test_simulator_is_deterministic(tmp_path: Path):
    a = tmp_path / "a"
    b = tmp_path / "b"
    da = simulate(a, seed=3, length=5_000, coverage=4, n_libraries=2, mosaic_fraction=0.5)
    db = simulate(b, seed=3, length=5_000, coverage=4, n_libraries=2, mosaic_fraction=0.5)
    assert da == db
    text_a = (a / "reads.sam").read_bytes()
    text_b = (b / "reads.sam").read_bytes()
    assert text_a == text_b
    assert hashlib.sha256(text_a).hexdigest() == da
    assert "LB:Z:lib1" in text_a.decode()
    assert (a / "truth.tsv").read_text(encoding="utf-8").count("DEL") == 1
    assert (a / "hap1.fa").read_text(encoding="utf-8").startswith(">chrSim")
