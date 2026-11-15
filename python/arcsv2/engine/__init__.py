"""Engine package. The Python modules are the reference implementation of §4.

The Rust crate ``arcsv-core`` implements the same kernels (KDE, clique search,
normalizer, path search, affected-length alignment). Paper mode is the default.
"""

from arcsv2.engine.caller import call_evidence

__all__ = ["call_evidence"]
