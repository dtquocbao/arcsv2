# arcsv-core audit (R0 / PERF2)

Existing pure-Rust crate under `core/arcsv-core` (library, drp, paths,
affected, names). Not wired to Python and not exercised against
caller-v1.4 golden outputs.

| Module | Role | vs Python caller-v1.4 |
|:--|:--|:--|
| `library.rs` | KDE / insert PMF | Unverified vs `engine/library.py`; keep for later kernel port |
| `drp.rs` | Discordant clustering sketches | Unverified vs `engine/breakpoints.py` |
| `paths.rs` | Path enumeration | Unverified vs `engine/graph.py` |
| `affected.rs` | Affected-length | Unverified vs `engine/classify.py` |
| `names.rs` | Block/path naming | Unverified vs `engine/names.py` |

Decision: **keep** the crate in the workspace for `cargo test` in CI. Do not
delete. New PyO3 surface is `native/arcsv2-native` (optional import). Reuse
arcsv-core pieces only after bit-identical tests against Python on golden
windows.
