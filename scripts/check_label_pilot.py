"""Check a simulated label.smk run: true events that were called are label 1."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import NamedTuple


class PilotResult(NamedTuple):
    matched: int
    missed: tuple[str, ...]
    n_events: int


class PilotMismatch(NamedTuple):
    truth: str
    call_type: str
    minbp: str
    maxbp: str
    label: str
    reason: str
    shift: str


def _rows(path: Path) -> list[dict[str, str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:] if line]


def _call_type(row: dict[str, str]) -> str:
    complex_type = row.get("complextype") or ""
    if complex_type not in ("", "NA"):
        return complex_type
    return (row.get("svtype") or "").split(",")[0]


def _overlaps(call: dict[str, str], event: dict[str, str], slop: int = 500) -> bool:
    if call.get("chrom") != event.get("chrom"):
        return False
    left = int(call["minbp"])
    right = int(call["maxbp"])
    start = int(event["start"]) - slop
    end = int(event["end"]) + slop
    return left < end and start < right


def _event_name(event: dict[str, str]) -> str:
    gt = event.get("gt") or ""
    return f"{event['type']} {gt} at {event['start']}-{event['end']}"


def assess_pilot(labeled: Path, truth: Path) -> tuple[PilotResult, list[PilotMismatch]]:
    """Score non-insertion truth events against the labeled calls."""
    calls = _rows(labeled)
    events = [row for row in _rows(truth) if row.get("type") != "INS"]
    matched = 0
    missed: list[str] = []
    mismatches: list[PilotMismatch] = []
    for event in events:
        same = [
            call
            for call in calls
            if _call_type(call) == event["type"] and _overlaps(call, event)
        ]
        if not same:
            missed.append(_event_name(event))
            continue
        bad = [call for call in same if call.get("label") != "1"]
        if bad:
            for call in bad:
                mismatches.append(
                    PilotMismatch(
                        _event_name(event),
                        _call_type(call),
                        call.get("minbp", ""),
                        call.get("maxbp", ""),
                        call.get("label", ""),
                        call.get("reason", ""),
                        call.get("label_bp_shift", "0"),
                    )
                )
        else:
            matched += 1
    result = PilotResult(matched, tuple(missed), len(events))
    return result, mismatches


def true_events_labeled(labeled: Path, truth: Path) -> int:
    """How many non-insertion truth events overlap a label-1 call of the same type."""
    result, mismatches = assess_pilot(labeled, truth)
    if mismatches:
        row = mismatches[0]
        raise SystemExit(
            f"{row.truth} overlaps a call labeled {row.label}, not 1"
        )
    return result.matched


def _print_report(result: PilotResult, mismatches: list[PilotMismatch]) -> None:
    for name in result.missed:
        print(f"missed {name}", file=sys.stderr)
    if mismatches:
        print("truth\tcall_type\tminbp\tmaxbp\tlabel\treason\tlabel_bp_shift", file=sys.stderr)
        for row in mismatches:
            print(
                f"{row.truth}\t{row.call_type}\t{row.minbp}\t{row.maxbp}\t"
                f"{row.label}\t{row.reason}\t{row.shift}",
                file=sys.stderr,
            )


def main(argv: list[str] | None = None) -> int:
    args = argv if argv is not None else sys.argv[1:]
    summary = Path(args[0] if args else "labels/dataset.summary.json")
    labeled = Path(args[1] if len(args) > 1 else "labels/sim.labeled.tsv")
    truth = Path(args[2] if len(args) > 2 else "pilot/truth.tsv")
    if not summary.is_file():
        raise SystemExit(f"missing {summary}")
    json.loads(summary.read_text(encoding="utf-8"))
    result, mismatches = assess_pilot(labeled, truth)
    _print_report(result, mismatches)
    need = min(3, result.n_events)
    if mismatches:
        raise SystemExit(f"{len(mismatches)} labeled call(s) overlap truth but are not label 1")
    if result.matched < need:
        missed = ", ".join(result.missed) or "none"
        raise SystemExit(
            f"{result.matched} of {result.n_events} non-insertion events labeled 1 "
            f"(need {need}). Missed: {missed}"
        )
    print(f"{result.matched} of {result.n_events} non-insertion events labeled 1")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
