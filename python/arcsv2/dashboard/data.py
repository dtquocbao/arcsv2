"""Load ARC-SV2 outputs for the dashboard and score them against a truth set.

Nothing here imports Streamlit, so the logic is unit-tested on its own.

A *results folder* holds one or more runs. A run is any folder that contains
``arcsv_out.tab`` (or ``merged.tab``). Optional files beside it:

- ``truth.tsv``: ``chrom start end type gt [size]`` for scoring a simulation
  or a benchmark region.
- ``run.json``: free-form run metadata (coverage, runtime, ...).
- ``labels/dataset.summary.json`` and ``model/manifest.json`` anywhere under
  the results folder: shown on the ML tab.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from arcsv2.config import OUTPUT_COLUMNS

CALL_FILES = ("arcsv_out.tab", "merged.tab")
NUMERIC = ("minbp", "maxbp", "len_affected", "num_sv", "score_vs_ref", "score_vs_next", "num_paths", "af")
OUTCOMES = ("correct", "wrong genotype", "wrong type", "missed")
SIZE_BINS = ((0, 500, "<500 bp"), (500, 1000, "500-999 bp"), (1000, 3000, "1-3 kb"), (3000, 10**12, "≥3 kb"))
TOKEN_RE = re.compile(r"([A-Za-z]\d*|_)('?)")


@dataclass
class Run:
    name: str
    calls_path: Path
    truth_path: Path | None = None
    meta: dict = field(default_factory=dict)


def discover_runs(root: str | Path) -> list[Run]:
    """Every folder under ``root`` with a call table, sorted by name."""
    base = Path(root)
    if not base.exists():
        return []
    found: dict[Path, Path] = {}
    for name in CALL_FILES:
        for path in sorted(base.rglob(name)):
            found.setdefault(path.parent, path)
    runs: list[Run] = []
    for folder, calls in sorted(found.items()):
        rel = folder.relative_to(base).as_posix()
        truth = folder / "truth.tsv"
        meta_path = folder / "run.json"
        meta = {}
        if meta_path.is_file():
            try:
                meta = json.loads(meta_path.read_text(encoding="utf-8"))
            except json.JSONDecodeError:
                meta = {"error": "run.json is not valid JSON"}
        runs.append(Run(rel if rel != "." else folder.name, calls, truth if truth.is_file() else None, meta))
    return runs


def call_type(row: dict | pd.Series) -> str:
    """``inv.del`` style name for complex calls, else the first svtype token."""
    complex_name = _text(row.get("complextype"))
    if complex_name not in ("", "NA"):
        return complex_name
    return _text(row.get("svtype")).split(",")[0] or "NA"


def _text(value) -> str:
    """Cell as text. Missing values (None, NaN) become an empty string."""
    if value is None or (isinstance(value, float) and value != value):
        return ""
    return str(value)


def load_calls(path: str | Path) -> pd.DataFrame:
    """Read a call table. Columns keep the caller's names; numbers are parsed."""
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    for col in OUTPUT_COLUMNS:
        if col not in frame.columns:
            frame[col] = "NA"
    for col in NUMERIC:
        frame[col] = pd.to_numeric(frame[col].replace("NA", None), errors="coerce")
    frame["call_type"] = [call_type(row) for row in frame.to_dict("records")]
    frame["call_id"] = [f"call{i + 1}" for i in range(len(frame))]
    return frame


def load_truth(path: str | Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    for col in ("start", "end"):
        frame[col] = pd.to_numeric(frame[col])
    if "size" in frame.columns:
        frame["size"] = pd.to_numeric(frame["size"], errors="coerce")
    else:
        frame["size"] = frame["end"] - frame["start"]
    frame["gt"] = frame["gt"].str.upper()
    frame["event_id"] = [f"event{i + 1}" for i in range(len(frame))]
    return frame


def size_bin(size: float) -> str:
    for low, high, label in SIZE_BINS:
        if low <= size < high:
            return label
    return SIZE_BINS[-1][2]


def match_truth(calls: pd.DataFrame, truth: pd.DataFrame, slop: int = 100) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Score calls against truth events.

    A call matches an event when it is on the same chromosome and its
    ``[minbp, maxbp]`` span overlaps the event widened by ``slop`` on each side
    (insertions are points). Each event takes its best matching call:

    - ``correct``: same type and genotype
    - ``wrong genotype``: same type, other genotype
    - ``wrong type``: overlapping call of another type
    - ``missed``: no overlapping call

    Calls that match no event are ``false call``. Returns ``(events, calls)``
    with an ``outcome`` column on each.
    """
    events = truth.copy()
    scored = calls.copy()
    scored["outcome"] = "false call"
    scored["event_id"] = ""
    outcomes: list[str] = []
    matched: list[str] = []
    rank = {name: i for i, name in enumerate(OUTCOMES)}
    for ev in events.to_dict("records"):
        lo, hi = ev["start"] - slop, max(ev["end"], ev["start"]) + slop
        mask = (scored["chrom"] == ev["chrom"]) & (scored["maxbp"] >= lo) & (scored["minbp"] <= hi)
        best_outcome, best_call = "missed", ""
        for idx, call in scored[mask].iterrows():
            if call["call_type"] == ev["type"]:
                outcome = "correct" if call["gt"] == ev["gt"] else "wrong genotype"
            else:
                outcome = "wrong type"
            if rank[outcome] < rank[best_outcome]:
                best_outcome, best_call = outcome, call["call_id"]
            scored.at[idx, "outcome"] = "matched"
            scored.at[idx, "event_id"] = ev["event_id"]
        outcomes.append(best_outcome)
        matched.append(best_call)
    events["outcome"] = outcomes
    events["call_id"] = matched
    events["size_bin"] = [size_bin(s) for s in events["size"]]
    return events, scored


def summarize(events: pd.DataFrame, calls: pd.DataFrame) -> dict:
    """Headline numbers for the overview tiles."""
    n = len(events)
    found = int(events["outcome"].isin(["correct", "wrong genotype"]).sum()) if n else 0
    correct = int((events["outcome"] == "correct").sum()) if n else 0
    return {
        "events": n,
        "calls": len(calls),
        "detected": found,
        "recall": found / n if n else float("nan"),
        "correct": correct,
        "genotype_concordance": correct / found if found else float("nan"),
        "false_calls": int((calls["outcome"] == "false call").sum()) if "outcome" in calls else 0,
    }


def tokens(text: str) -> list[tuple[str, bool]]:
    """``AB'C`` -> ``[('A', False), ('B', True), ('C', False)]``. ``_`` is a novel insertion."""
    return [(m.group(1), bool(m.group(2))) for m in TOKEN_RE.finditer(text or "")]


def block_lengths(row: dict | pd.Series) -> dict[str, int]:
    """Reference block name -> length, from the ``blocks`` column when present."""
    names = [name for name, _inv in tokens(str(row.get("reference") or "")) if name != "_"]
    spans: list[int] = []
    for part in str(row.get("blocks") or "").split(","):
        part = part.strip()
        if "-" in part:
            left, right = part.split("-", 1)
            try:
                spans.append(int(right) - int(left))
            except ValueError:
                pass
    if len(spans) != len(names):
        bounds = [int(x) for x in str(row.get("bp") or "").split(",") if x.strip().lstrip("-").isdigit()]
        spans = [b - a for a, b in zip(bounds, bounds[1:])]
    if len(spans) != len(names):
        spans = [1] * len(names)
    return dict(zip(names, spans))


def insertion_length(row: dict | pd.Series) -> int:
    text = str(row.get("inslen") or "NA").replace(">", "").split(",")[0]
    try:
        return int(float(text))
    except ValueError:
        return 0


def haplotype_layout(row: dict | pd.Series, flank_cap: int = 2000) -> pd.DataFrame:
    """Segments for drawing the reference and the called haplotypes.

    Outer flank blocks are capped at ``flank_cap`` bp so the event stays
    visible. HET calls draw one alternative and one reference haplotype.
    """
    lengths = block_lengths(row)
    names = list(lengths)
    shown = {
        name: (min(length, flank_cap) if i in (0, len(names) - 1) else length)
        for i, (name, length) in enumerate(lengths.items())
    }
    ins = max(insertion_length(row), 1)
    tracks = [("Reference", str(row.get("reference") or ""))]
    alt = str(row.get("rearrangement") or "")
    if str(row.get("gt") or "").upper() == "HOM":
        tracks += [("Haplotype 1", alt), ("Haplotype 2", alt)]
    else:
        tracks += [("Haplotype 1", alt), ("Haplotype 2", tracks[0][1])]
    rows: list[dict] = []
    for track, text in tracks:
        x = 0
        for name, inverted in tokens(text):
            width = ins if name == "_" else shown.get(name, 1)
            rows.append(
                {
                    "track": track,
                    "block": "novel" if name == "_" else name,
                    "label": ("ins" if name == "_" else name) + ("'" if inverted else ""),
                    "orientation": "inverted" if inverted else "forward",
                    "start": x,
                    "end": x + width,
                    "length_bp": ins if name == "_" else lengths.get(name, 0),
                    "capped": name in (names[0], names[-1]) if names else False,
                }
            )
            x += width
    return pd.DataFrame(rows)


def load_json(path: Path) -> dict | None:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def rearrangement_plain(row: dict | pd.Series) -> str:
    """Plain-language summary of a rearrangement string (U6)."""
    ref = _text(row.get("reference"))
    alt = _text(row.get("rearrangement"))
    sv = call_type(row)
    if alt in ("", "NA") and ref in ("", "NA"):
        return f"A {sv} variant was called at this locus."
    if alt == ref:
        return f"No change from the reference layout ({ref})."
    parts: list[str] = []
    for name, inverted in tokens(alt):
        if name == "_":
            parts.append("an inserted segment")
        elif name == "=":
            parts.append("an insertion of unknown sequence")
        else:
            seg = f"segment {name}"
            if inverted:
                seg += " (reverse orientation)"
            parts.append(seg)
    if len(parts) == 1 and sv in ("DEL", "Dup", "DUP", "INV"):
        return f"{sv}: {parts[0]} compared with reference {ref}."
    joined = ", ".join(parts[:-1]) + (" and " + parts[-1] if len(parts) > 1 else parts[0])
    return f"The alternate haplotype orders {joined} (reference layout {ref})."


def load_batch_status(path: Path) -> pd.DataFrame:
    """Read U4 ``batch_status.tsv``."""
    if not path.is_file():
        return pd.DataFrame(columns=["sample", "status", "message", "updated_utc"])
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def find_batch_status_files(root: str | Path) -> list[Path]:
    base = Path(root)
    if not base.is_dir():
        return []
    return sorted(base.rglob("batch_status.tsv"))


def find_ml_artifacts(root: str | Path) -> dict[str, Path]:
    """Label summaries, labeled tables and model manifests under ``root``."""
    base = Path(root)
    found: dict[str, Path] = {}
    if not base.exists():
        return found
    for path in sorted(base.rglob("dataset.summary.json")):
        found.setdefault("summary", path)
    for path in sorted(base.rglob("manifest.json")):
        found.setdefault("manifest", path)
    for path in sorted(base.rglob("dataset.tsv")):
        found.setdefault("dataset", path)
    return found


# --- Project status -------------------------------------------------------


def find_repo_root(start: str | Path | None = None) -> Path | None:
    """The arcsv2 source checkout, if the package runs from one."""
    candidates = [Path(start)] if start else []
    candidates.append(Path(__file__).resolve())
    candidates.append(Path.cwd())
    for origin in candidates:
        for folder in [origin, *origin.parents]:
            project = folder / "pyproject.toml"
            if project.is_file() and 'name = "arcsv2"' in project.read_text(encoding="utf-8", errors="ignore"):
                return folder
    return None


def _contains(path: Path, needle: str) -> bool:
    try:
        return needle in path.read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return False


def milestones(repo: Path | None, results: Path | None = None) -> pd.DataFrame:
    """Milestones from ``docs/LABELS_PLAN.md`` and ``docs/CROSS_PLATFORM.md``, detected from files."""
    rows: list[tuple[str, str, bool, str]] = []
    if repo is None:
        return pd.DataFrame(columns=["milestone", "item", "done", "evidence"])
    py = repo / "python" / "arcsv2"
    results_dir = results or repo
    ml = find_ml_artifacts(results_dir)
    manifest = load_json(ml["manifest"]) if "manifest" in ml else None
    if not isinstance(manifest, dict):
        manifest = None
    model_kinds = set((manifest or {}).get("models", {}).values())
    bench = repo / "bench" / "RESULTS.md"
    bench_rows = [
        line for line in (bench.read_text(encoding="utf-8").splitlines() if bench.is_file() else [])
        if line.startswith("| 20")
    ]
    checks = [
        ("M0", "BAM/CRAM input (pysam)", _contains(py / "samio.py", "def iter_bam"), "samio.py: iter_bam"),
        ("M0", "blocks column in call table", _contains(py / "engine" / "caller.py", '"blocks"'), "caller.py writes blocks"),
        ("M1", "ml-label with haplotype alignment", _contains(py / "label.py", "class HapAligner"), "label.py: HapAligner"),
        ("M1", "label tests", (repo / "tests" / "test_label.py").is_file(), "tests/test_label.py"),
        ("M1", "labels accept either strand", not _contains(py / "label.py", "hit.strand < 0"), "label.py strand rule"),
        ("Platform", "container environment", (repo / "environment.yml").is_file(), "environment.yml"),
        ("Platform", "dev container", (repo / ".devcontainer" / "devcontainer.json").is_file(), ".devcontainer/"),
        ("M2", "pilot timing recorded", bool(bench_rows), "bench/RESULTS.md rows"),
        ("M3", "grouped CV + scaling in training", _contains(py / "train.py", "GroupKFold"), "train.py: GroupKFold"),
        ("M4", "labeled dataset", "dataset" in ml or "summary" in ml, "dataset.tsv / summary under results"),
        ("M5", "model comparison", len(model_kinds - {"logistic"}) > 0, "manifest models beyond logistic"),
    ]
    for milestone, item, done, evidence in checks:
        rows.append((milestone, item, bool(done), evidence))
    return pd.DataFrame(rows, columns=["milestone", "item", "done", "evidence"])


def review_documents(repo: Path | None) -> list[tuple[str, Path]]:
    """``(title, path)`` for each review round, newest last."""
    if repo is None:
        return []
    folder = repo / "docs" / "reviews"
    if not folder.is_dir():
        return []
    docs: list[tuple[str, Path]] = []
    for path in sorted(folder.glob("*.md"), key=_review_order):
        title = path.stem
        for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
            if line.startswith("# "):
                title = line[2:].strip()
                break
        docs.append((title, path))
    return docs


def _review_order(path: Path) -> tuple[str, int]:
    match = re.search(r"round(\d+)", path.stem)
    return (path.stem.split("-")[0], int(match.group(1)) if match else 1)
