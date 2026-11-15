"""Dashboard data layer and a headless render of every tab."""

from pathlib import Path

import pytest

pd = pytest.importorskip("pandas")

from arcsv2.dashboard import data  # noqa: E402

DEMO = Path(__file__).resolve().parents[1] / "python" / "arcsv2" / "dashboard" / "demo"


def _calls(rows):
    frame = pd.DataFrame(rows)
    for col in ("complextype", "reference", "rearrangement", "blocks", "bp"):
        if col not in frame:
            frame[col] = "NA"
    frame["call_type"] = [data.call_type(r) for r in frame.to_dict("records")]
    frame["call_id"] = [f"call{i + 1}" for i in range(len(frame))]
    return frame


def test_match_truth_outcomes():
    truth = pd.DataFrame(
        [
            {"chrom": "c", "start": 1000, "end": 2000, "type": "DEL", "gt": "HOM", "size": 1000},
            {"chrom": "c", "start": 5000, "end": 6000, "type": "DUP", "gt": "HOM", "size": 1000},
            {"chrom": "c", "start": 9000, "end": 9500, "type": "INV", "gt": "HET", "size": 500},
            {"chrom": "c", "start": 20000, "end": 20000, "type": "INS", "gt": "HET", "size": 300},
            {"chrom": "c", "start": 30000, "end": 31000, "type": "DEL", "gt": "HET", "size": 1000},
        ]
    )
    truth["event_id"] = [f"event{i}" for i in range(len(truth))]
    calls = _calls(
        [
            {"chrom": "c", "minbp": 990, "maxbp": 2010, "svtype": "DEL", "gt": "HOM"},
            {"chrom": "c", "minbp": 5010, "maxbp": 5990, "svtype": "DUP", "gt": "HET"},
            {"chrom": "c", "minbp": 9020, "maxbp": 9480, "svtype": "BND,BND", "complextype": "other", "gt": "HET"},
            {"chrom": "c", "minbp": 20000, "maxbp": 20000, "svtype": "INS", "gt": "HET"},
            {"chrom": "c", "minbp": 50000, "maxbp": 50500, "svtype": "DEL", "gt": "HET"},
        ]
    )
    events, scored = data.match_truth(calls, truth)
    assert list(events["outcome"]) == ["correct", "wrong genotype", "wrong type", "correct", "missed"]
    assert list(scored["outcome"]) == ["matched"] * 4 + ["false call"]
    summary = data.summarize(events, scored)
    assert summary["detected"] == 3 and summary["correct"] == 2 and summary["false_calls"] == 1


def test_complex_type_uses_complextype():
    assert data.call_type({"svtype": "BND,BND", "complextype": "inv.del"}) == "inv.del"
    assert data.call_type({"svtype": "DEL", "complextype": "NA"}) == "DEL"


def test_haplotype_layout_het_inversion():
    row = {"reference": "ABC", "rearrangement": "AB'C", "gt": "HET", "blocks": "0-5000,5000-6000,6000-9000", "inslen": "NA"}
    layout = data.haplotype_layout(row, flank_cap=2000)
    hap1 = layout[layout["track"] == "Haplotype 1"]
    assert list(hap1["label"]) == ["A", "B'", "C"]
    assert list(hap1["orientation"]) == ["forward", "inverted", "forward"]
    assert list(hap1["length_bp"]) == [5000, 1000, 3000]
    assert hap1["end"].max() == 2000 + 1000 + 2000  # flanks capped
    hap2 = layout[layout["track"] == "Haplotype 2"]
    assert list(hap2["label"]) == ["A", "B", "C"]


def test_insertion_layout_uses_inslen():
    row = {"reference": "AB", "rearrangement": "A_B", "gt": "HOM", "blocks": "0-100,100-200", "inslen": ">316"}
    layout = data.haplotype_layout(row)
    novel = layout[(layout["track"] == "Haplotype 1") & (layout["block"] == "novel")]
    assert int(novel["length_bp"].iloc[0]) == 316


def test_demo_runs_are_bundled_and_score():
    runs = data.discover_runs(DEMO)
    assert [r.name for r in runs] == ["simulated_2Mb_15x", "simulated_2Mb_30x"]
    for run in runs:
        assert run.truth_path is not None and run.meta.get("coverage")
        calls = data.load_calls(run.calls_path)
        events, scored = data.match_truth(calls, data.load_truth(run.truth_path))
        assert len(events) == 40
        assert data.summarize(events, scored)["recall"] > 0.5


def test_milestones_detect_repo():
    repo = data.find_repo_root(Path(__file__))
    assert repo is not None
    table = data.milestones(repo)
    assert {"M0", "M1", "M3"} <= set(table["milestone"])


def test_app_renders_every_tab():
    pytest.importorskip("streamlit")
    try:
        import pyarrow  # noqa: F401
    except (ImportError, AttributeError) as exc:
        pytest.skip(f"pyarrow cannot be imported in this environment: {exc}")
    from streamlit.testing.v1 import AppTest

    app = Path(__file__).resolve().parents[1] / "python" / "arcsv2" / "dashboard" / "app.py"
    at = AppTest.from_file(str(app), default_timeout=60)
    at.run()
    assert not at.exception, at.exception
    assert len(at.tabs) == 6
    assert at.tabs[0].label == "Guide"
    assert any("True events" in m.label for m in at.metric)
