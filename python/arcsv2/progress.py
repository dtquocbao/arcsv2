"""Stderr progress lines for long-running call stages."""

from __future__ import annotations

import sys
import time
from dataclasses import dataclass, field


@dataclass
class Progress:
    """Timestamped stderr messages. Opt out with ``quiet=True``."""

    label: str = "arcsv2"
    quiet: bool = False
    t0: float = field(default_factory=time.monotonic)

    def log(self, message: str) -> None:
        if self.quiet:
            return
        elapsed = time.monotonic() - self.t0
        sys.stderr.write(f"[{self.label} {elapsed:7.1f}s] {message}\n")
        sys.stderr.flush()

    def stage(self, name: str, detail: str = "") -> None:
        msg = name if not detail else f"{name}: {detail}"
        self.log(msg)
