"""ARC-SV2 dashboard. Start it with ``arcsv2 dashboard``.

The results folder comes from ``ARCSV2_RESULTS`` (set by the CLI). Without it,
the bundled simulation demo is shown.
"""

from __future__ import annotations

import os
from pathlib import Path

import altair as alt
import pandas as pd
import streamlit as st

from arcsv2 import __version__
from arcsv2.dashboard import data

DEMO_DIR = Path(__file__).resolve().parent / "demo"

# Reference palette (dataviz skill): categorical slots in fixed order, and the
# reserved status palette. Status colors always ship with an icon + label.
BLOCK_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
NOVEL_COLOR = "#8a8985"
OUTCOME_STYLE = {
    "correct": ("#0ca30c", "✅"),
    "wrong genotype": ("#fab219", "⚠️"),
    "wrong type": ("#ec835a", "🔶"),
    "missed": ("#d03b3b", "❌"),
}
TYPE_ORDER = ["DEL", "DUP", "INV", "INS", "inv.del"]

st.set_page_config(page_title="ARC-SV2 dashboard", page_icon="🧬", layout="wide")


# --- data access ------------------------------------------------------------


@st.cache_data(show_spinner=False)
def _calls(path: str, mtime: float) -> pd.DataFrame:
    return data.load_calls(path)


@st.cache_data(show_spinner=False)
def _truth(path: str, mtime: float) -> pd.DataFrame:
    return data.load_truth(path)


def load_run(run: data.Run) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    calls = _calls(str(run.calls_path), run.calls_path.stat().st_mtime)
    if run.truth_path is None:
        return calls, None
    truth = _truth(str(run.truth_path), run.truth_path.stat().st_mtime)
    events, scored = data.match_truth(calls, truth)
    return scored, events


def results_folder() -> tuple[Path, bool]:
    configured = os.environ.get("ARCSV2_RESULTS", "").strip()
    if configured:
        return Path(configured), False
    return DEMO_DIR, True


# --- charts -----------------------------------------------------------------


def outcome_chart(events: pd.DataFrame, by: str, title: str, order: list[str]) -> alt.Chart:
    counts = events.groupby([by, "outcome"]).size().reset_index(name="events")
    return (
        alt.Chart(counts, title=title)
        .mark_bar(cornerRadiusEnd=4)
        .encode(
            y=alt.Y(f"{by}:N", sort=order, title=None, scale=alt.Scale(domain=order, paddingInner=0.35)),
            x=alt.X("events:Q", title="Events", axis=alt.Axis(tickMinStep=1)),
            color=alt.Color(
                "outcome:N",
                scale=alt.Scale(domain=list(OUTCOME_STYLE), range=[c for c, _ in OUTCOME_STYLE.values()]),
                legend=alt.Legend(title="Outcome", orient="bottom"),
            ),
            order=alt.Order("outcome_rank:Q"),
            tooltip=[alt.Tooltip(f"{by}:N", title=by.replace("_", " ")), "outcome:N", "events:Q"],
        )
        .transform_calculate(outcome_rank=f"indexof({list(OUTCOME_STYLE)}, datum.outcome)")
        .properties(height=alt.Step(34))
    )


def haplotype_chart(row: pd.Series) -> alt.Chart:
    layout = data.haplotype_layout(row)
    blocks = [b for b in dict.fromkeys(layout["block"]) if b != "novel"]
    domain = blocks + (["novel"] if "novel" in set(layout["block"]) else [])
    colors = [BLOCK_COLORS[i % len(BLOCK_COLORS)] for i in range(len(blocks))] + (
        [NOVEL_COLOR] if "novel" in domain else []
    )
    layout["mid"] = (layout["start"] + layout["end"]) / 2
    layout["glyph"] = layout["label"] + layout["orientation"].map({"forward": " ▶", "inverted": " ◀"})
    tracks = ["Reference", "Haplotype 1", "Haplotype 2"]
    base = alt.Chart(layout).encode(y=alt.Y("track:N", sort=tracks, title=None))
    bars = base.mark_bar(height=26, cornerRadius=3, stroke="white", strokeWidth=2).encode(
        x=alt.X("start:Q", title="Position along haplotype (bp, flanks capped at 2 kb)"),
        x2="end:Q",
        color=alt.Color("block:N", scale=alt.Scale(domain=domain, range=colors), legend=None),
        tooltip=[
            alt.Tooltip("track:N"),
            alt.Tooltip("label:N", title="block"),
            alt.Tooltip("orientation:N"),
            alt.Tooltip("length_bp:Q", title="length (bp)", format=","),
        ],
    )
    labels = base.mark_text(color="white", fontWeight="bold", fontSize=12).encode(x="mid:Q", text="glyph:N")
    return (bars + labels).properties(height=150)


# --- tabs -------------------------------------------------------------------


def tab_overview(runs: list[data.Run], run: data.Run, calls: pd.DataFrame, events: pd.DataFrame | None) -> None:
    if events is None:
        st.info(
            "No `truth.tsv` beside this call table — normal for real 1KGP/HPRC runs. "
            "Accuracy bars need a simulation or benchmark that ships truth. "
            "Call counts and the explorer still work; labels are on the ML tab."
        )
        c1, c2, c3 = st.columns(3)
        c1.metric("Calls", len(calls))
        c2.metric("Homozygous", int((calls["gt"] == "HOM").sum()))
        c3.metric("Median score vs reference", f"{calls['score_vs_ref'].median():.1f}" if len(calls) else "–")
        counts = calls.groupby("call_type").size().reset_index(name="calls")
        st.altair_chart(
            alt.Chart(counts, title="Calls by type")
            .mark_bar(cornerRadiusEnd=4, color=BLOCK_COLORS[0])
            .encode(y=alt.Y("call_type:N", title=None, sort="-x", scale=alt.Scale(paddingInner=0.35)), x=alt.X("calls:Q", title="Calls"), tooltip=["call_type", "calls"])
            .properties(height=alt.Step(34)),
            width="stretch",
        )
        return
    summary = data.summarize(events, calls)
    c1, c2, c3, c4, c5 = st.columns(5)
    c1.metric("True events", summary["events"])
    c2.metric("Detected (right type)", f"{summary['recall']:.0%}", help=f"{summary['detected']} of {summary['events']}")
    c3.metric("Genotype concordance", f"{summary['genotype_concordance']:.0%}", help="Among detected events")
    c4.metric("False calls", summary["false_calls"], help="Calls that overlap no true event")
    seconds = run.meta.get("call_seconds")
    c5.metric("Call time", f"{seconds:.0f} s" if isinstance(seconds, (int, float)) else "–", help=f"{run.meta.get('pairs', '?'):,} read pairs" if isinstance(run.meta.get("pairs"), int) else None)

    left, right = st.columns(2)
    types = [t for t in TYPE_ORDER if t in set(events["type"])] + sorted(set(events["type"]) - set(TYPE_ORDER))
    left.altair_chart(outcome_chart(events, "type", "Outcome by SV type", types), width="stretch")
    right.altair_chart(
        outcome_chart(events, "size_bin", "Outcome by event size", [label for _, _, label in data.SIZE_BINS]),
        width="stretch",
    )

    problems = events[events["outcome"] != "correct"]
    st.subheader(f"Events not called correctly ({len(problems)})")
    if problems.empty:
        st.success("✅ Every true event was called with the right type and genotype.")
    else:
        merged = problems.merge(
            calls[["call_id", "call_type", "gt", "rearrangement", "score_vs_ref", "score_vs_next"]].rename(columns={"gt": "called_gt"}),
            on="call_id",
            how="left",
        )
        merged["outcome"] = merged["outcome"].map(lambda o: f"{OUTCOME_STYLE[o][1]} {o}")
        for col in ("call_type", "called_gt", "rearrangement"):
            merged[col] = merged[col].fillna("–")
        st.dataframe(
            merged[["outcome", "type", "gt", "size", "start", "call_type", "called_gt", "rearrangement", "score_vs_ref", "score_vs_next"]],
            hide_index=True,
            width="stretch",
            column_config={
                "type": "true type",
                "gt": "true GT",
                "size": st.column_config.NumberColumn("size (bp)", format="%d"),
                "start": st.column_config.NumberColumn("position", format="%d"),
                "call_type": "called type",
                "called_gt": "called GT",
                "score_vs_ref": st.column_config.NumberColumn("score vs ref", format="%.1f"),
                "score_vs_next": st.column_config.NumberColumn("score vs next", format="%.1f"),
            },
        )

    if len(runs) > 1:
        st.subheader("All runs")
        table = []
        for other in runs:
            other_calls, other_events = load_run(other)
            row = {"run": other.name, "calls": len(other_calls), "coverage": other.meta.get("coverage"), "call seconds": other.meta.get("call_seconds")}
            if other_events is not None:
                s = data.summarize(other_events, other_calls)
                row.update({"detected": s["recall"], "genotype concordance": s["genotype_concordance"], "false calls": s["false_calls"]})
            table.append(row)
        st.dataframe(
            pd.DataFrame(table),
            hide_index=True,
            width="stretch",
            column_config={
                "detected": st.column_config.ProgressColumn("detected", format="percent", min_value=0, max_value=1),
                "genotype concordance": st.column_config.ProgressColumn("genotype concordance", format="percent", min_value=0, max_value=1),
            },
        )


PIPELINE_STEPS = [
    ("Reads", "SAM, BAM or CRAM for one region or chromosome.", "samio.py"),
    ("Library model", "Insert-size distribution (KDE) and the cut-offs that make a pair discordant.", "engine/library.py"),
    ("Breakpoints", "Soft clips, split reads and discordant-pair clusters, kept at 10% FDR.", "engine/breakpoints.py"),
    ("Blocks & graph", "Breakpoints cut the reference into blocks; read pairs support new joins.", "engine/graph.py"),
    ("Loci", "Cut points split the graph into independent regions.", "engine/graph.py"),
    ("Candidate haplotypes", "Every walk through a locus's blocks (up to 1,000).", "engine/graph.py"),
    ("Diploid genotyping", "Likelihood of all reads for each pair of haplotypes; best pair wins.", "engine/likelihood.py"),
    ("Call table", "Rearrangement, genotype, scores and support per call.", "arcsv_out.tab"),
    ("ML filter", "Per-class model keeps calls at precision chosen with recall ≥ 0.90.", "train.py / predict.py"),
    ("Filtered calls", "Final table and VCF with breakends.", "vcf.py"),
]


def narrative(row: pd.Series) -> list[str]:
    lengths = data.block_lengths(row)
    alt_tokens = data.tokens(str(row["rearrangement"]))
    ref_names = list(lengths)
    kept = {name for name, _ in alt_tokens}
    deleted = [n for n in ref_names if n not in kept]
    inverted = sorted({name for name, inv in alt_tokens if inv})
    counts = pd.Series([n for n, _ in alt_tokens if n != "_"]).value_counts()
    duplicated = [n for n, c in counts.items() if c > 1]
    changes = []
    if deleted:
        changes.append("deleted " + ", ".join(f"{n} ({lengths[n]:,} bp)" for n in deleted))
    if inverted:
        changes.append("inverted " + ", ".join(f"{n} ({lengths.get(n, 0):,} bp)" for n in inverted))
    if duplicated:
        changes.append("duplicated " + ", ".join(f"{n} ({lengths.get(n, 0):,} bp)" for n in duplicated))
    if "_" in kept:
        changes.append(f"inserted novel sequence ({row.get('inslen', 'NA')} bp)")
    sr = str(row.get("sr_support", "0"))
    pe = str(row.get("pe_support", "0"))
    return [
        f"**1. Breakpoints.** Evidence at `{row.get('sv_bp', 'NA')}` split the reference into blocks "
        f"{', '.join(f'{n} ({lengths[n]:,} bp)' for n in ref_names)}.",
        f"**2. Graph.** Read pairs and split reads that join blocks in a new way become adjacencies. "
        f"Per new junction: split-read support {sr}, read-pair support {pe}.",
        f"**3. Candidate haplotypes.** The path search found {int(row['num_paths']) if pd.notna(row['num_paths']) else '?'} "
        f"ways to walk the blocks in this locus.",
        f"**4. Genotype.** The best diploid genotype is **{row['gt']}** for `{row['reference']} → {row['rearrangement']}`"
        f" ({'; '.join(changes) or 'no change'}). It beats the reference by {row['score_vs_ref']:.1f} log-likelihood units and the next "
        f"alternative (`{row.get('rearrangement_next', 'NA')}`) by {row['score_vs_next']:.1f}.",
        f"**5. Type.** Reported as `{row['call_type']}` with {row.get('len_affected', 'NA'):,.0f} bp affected."
        if pd.notna(row.get("len_affected"))
        else f"**5. Type.** Reported as `{row['call_type']}`.",
    ]


def call_label(row: pd.Series) -> str:
    return f"{row['call_id']} · {row['call_type']} {row['gt']} · {row['chrom']}:{int(row['minbp']):,} · {row['reference']} → {row['rearrangement']}"


def call_detail(row: pd.Series) -> None:
    st.altair_chart(haplotype_chart(row), width="stretch")
    st.caption("Each color is one reference block. ◀ marks an inverted block; grey is novel inserted sequence. For a HET call, haplotype 2 keeps the reference.")
    for line in narrative(row):
        st.markdown(line)
    with st.expander("All columns for this call"):
        st.dataframe(row.drop(labels=["call_id"]).astype(str).rename("value"), width="stretch")


def tab_how(calls: pd.DataFrame) -> None:
    st.markdown("#### Pipeline")
    st.markdown(
        "ARC-SV2 splits the genome at candidate breakpoints, links the pieces with read evidence, "
        "tries every way the pieces can be walked, and keeps the pair of haplotypes with the highest "
        "likelihood. An ML filter then scores each call."
    )
    for row_start in (0, 5):
        cols = st.columns(5)
        for offset, col in enumerate(cols):
            number = row_start + offset
            title, detail, where = PIPELINE_STEPS[number]
            with col.container(border=True):
                st.markdown(f"**{number + 1} · {title}**")
                st.caption(detail)
                st.markdown(f"<span style='font-size:0.75rem;opacity:0.7'>`{where}`</span>", unsafe_allow_html=True)
    st.markdown(
        "Training labels come from a side branch: calls are aligned to phased long-read assemblies "
        "(HPRC, HG002 Q100) with `arcsv2 ml-label`, and the labeled calls train step 9."
    )
    st.markdown("#### Worked example")
    if calls.empty:
        st.info("No calls in this run.")
        return
    # Complex events first: they show the method best.
    ordered = calls.assign(_complex=~calls["call_type"].isin(TYPE_ORDER[:4])).sort_values(
        ["_complex", "call_type", "score_vs_ref"], ascending=[False, True, False]
    )
    labels = {row["call_id"]: call_label(row) for _, row in ordered.iterrows()}
    choice = st.selectbox("Pick a call", list(labels), format_func=labels.get)
    call_detail(calls[calls["call_id"] == choice].iloc[0])


def tab_explorer(calls: pd.DataFrame, events: pd.DataFrame | None) -> None:
    f1, f2, f3, f4 = st.columns([2, 1, 2, 2])
    types = sorted(calls["call_type"].unique())
    chosen_types = f1.multiselect("SV type", types, default=types)
    gt = f2.selectbox("Genotype", ["All", "HOM", "HET"])
    max_score = float(calls["score_vs_ref"].max()) if len(calls) else 0.0
    min_score = f3.slider("Min score vs reference", 0.0, max(1.0, float(round(max_score))), 0.0)
    outcomes = ["All"] + (["matched", "false call"] if events is not None else [])
    outcome = f4.selectbox("Truth match", outcomes, disabled=events is None)
    view = calls[calls["call_type"].isin(chosen_types) & (calls["score_vs_ref"].fillna(0) >= min_score)]
    if gt != "All":
        view = view[view["gt"] == gt]
    if outcome != "All" and "outcome" in view:
        view = view[view["outcome"] == outcome]
    columns = ["call_id", "chrom", "minbp", "maxbp", "call_type", "gt", "reference", "rearrangement", "len_affected", "inslen", "sr_support", "pe_support", "score_vs_ref", "score_vs_next", "num_paths"]
    if "outcome" in view:
        columns.append("outcome")
    st.caption(f"{len(view)} of {len(calls)} calls. Select a row to see how it was called.")
    selection = st.dataframe(
        view[columns],
        hide_index=True,
        width="stretch",
        on_select="rerun",
        selection_mode="single-row",
        column_config={
            "minbp": st.column_config.NumberColumn("start", format="%d"),
            "maxbp": st.column_config.NumberColumn("end", format="%d"),
            "call_type": "type",
            "len_affected": st.column_config.NumberColumn("affected bp", format="%d"),
            "score_vs_ref": st.column_config.NumberColumn("score vs ref", format="%.1f"),
            "score_vs_next": st.column_config.NumberColumn("score vs next", format="%.1f"),
        },
    )
    st.download_button("Download filtered calls (TSV)", view.to_csv(sep="\t", index=False), file_name="arcsv2_calls.tsv")
    rows = selection.selection.rows if selection is not None else []
    if rows:
        st.divider()
        call_detail(view.iloc[rows[0]])


def tab_ml(results: Path) -> None:
    found = data.find_ml_artifacts(results)
    if not found:
        st.info(
            "No labels or trained model under this results folder yet.\n\n"
            "Next steps: `./scripts/arcsv2.ps1 fetch --n-samples 2 --chroms chr20 --yes`, "
            "then Snakemake via the wrapper with `workflows/config/label.auto.yaml` "
            "(see the **Guide** tab). Point the dashboard at the folder that holds "
            "`labels/dataset.tsv`."
        )
        return
    if "summary" in found:
        st.markdown(f"#### Label summary · `{found['summary'].name}`")
        summary = data.load_json(found["summary"]) or {}
        st.json(summary, expanded=False)
    if "dataset" in found:
        labeled = pd.read_csv(found["dataset"], sep="\t", dtype=str, keep_default_na=False)
        if {"sv_class", "label"} <= set(labeled.columns):
            counts = labeled.groupby(["sv_class", "label"]).size().reset_index(name="calls")
            counts["label"] = counts["label"].map({"1": "real (1)", "0": "false (0)"}).fillna("NA")
            st.altair_chart(
                alt.Chart(counts, title="Labeled calls by class")
                .mark_bar(cornerRadiusEnd=4, size=24)
                .encode(
                    y=alt.Y("sv_class:N", title=None),
                    x=alt.X("calls:Q", title="Calls"),
                    color=alt.Color("label:N", scale=alt.Scale(domain=["real (1)", "false (0)", "NA"], range=[BLOCK_COLORS[0], BLOCK_COLORS[1], NOVEL_COLOR])),
                    tooltip=["sv_class", "label", "calls"],
                ),
                width="stretch",
            )
    if "manifest" in found:
        manifest = data.load_json(found["manifest"]) or {}
        st.markdown(f"#### Model · `{found['manifest']}`")
        fitted = manifest.get("fitted", False)
        st.markdown(("✅ **Trained**" if fitted else "⏳ **Untrained presets**") + f" · feature schema `{manifest.get('feature_schema', '?')}`")
        rows = [
            {"class": kind, "model": manifest.get("models", {}).get(kind, "–"), "threshold": manifest.get("thresholds", {}).get(kind), "CV grouping": manifest.get("grouping", {}).get(kind, "–")}
            for kind in sorted(set(manifest.get("models", {})) | set(manifest.get("thresholds", {})))
        ]
        if rows:
            st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")
        if manifest.get("chromosome_holdout"):
            with st.expander("Chromosome hold-out check"):
                st.json(manifest["chromosome_holdout"], expanded=False)


def tab_guide() -> None:
    """Step-by-step runbook for Windows/macOS/Linux via Compose."""
    st.markdown("#### Run ARC-SV2 end to end")
    st.caption(
        "Heavy tools (Snakemake, samtools, mappy) run inside the Linux container. "
        "Use `scripts/arcsv2.ps1` on Windows or `scripts/arcsv2.sh` on macOS and Linux."
    )
    platform = st.radio(
        "Your host OS",
        ["Windows (PowerShell)", "macOS / Linux (bash)"],
        horizontal=True,
        key="guide_platform",
    )
    windows = platform.startswith("Windows")
    wrap = "./scripts/arcsv2.ps1" if windows else "./scripts/arcsv2.sh"
    lang = "powershell" if windows else "bash"
    cont = "`" if windows else "\\"

    def code(*lines: str) -> str:
        return f"```{lang}\n" + "\n".join(lines) + "\n```"

    snake = code(
        f"{wrap} --entrypoint snakemake {cont}",
        f"  -s workflows/label.smk {cont}",
        "  --configfile workflows/config/label.auto.yaml -j 2",
    )
    steps = [
        (
            "1 · Build the image once",
            "Creates the Compose project `arcsv2` and the `arcsv2:dev` image. "
            "Needs Docker Desktop (Windows/macOS) or Docker Engine (Linux).",
            code(f"{wrap} build"),
        ),
        (
            "2 · Plan the download (optional)",
            "Intersects HPRC Release 2 with 1KGP 30× CRAMs, prints sample IDs and GB. "
            "Does not download assemblies.",
            code(f"{wrap} fetch --dry-run --n-samples 2 --chroms chr20"),
        ),
        (
            "3 · Fetch reference + haplotype assemblies",
            "Writes gitignored `workflows/config/samples.tsv` and `label.auto.yaml`, "
            "downloads ~6 GB into the `arcsv2-data` volume at `/data`. "
            "CRAMs are streamed later by chromosome (not downloaded whole).",
            code(f"{wrap} fetch --n-samples 2 --chroms chr20 --yes"),
        ),
        (
            "4 · Call and label (Snakemake)",
            "Always use the wrapper — do not call `snakemake` on the host PATH. "
            "Omit `--use-conda` inside the image.",
            snake,
        ),
        (
            "5 · Open this dashboard on your results",
            "Fetch/Snakemake write under `/data` (Compose volume). "
            "The wrapper defaults to that path. Point at the repo for the simulated pilot.",
            code(
                f"{wrap} dashboard                 # defaults to /data",
                f"{wrap} dashboard --results /data",
                f"{wrap} dashboard --results /work   # simulated pilot under the repo",
                '# or natively, after: pip install -e ".[dashboard]"',
                "arcsv2 dashboard --results .",
            ),
        ),
        (
            "6 · Record timing",
            "Fill wall time, peak RSS, call count, and git revision in `bench/RESULTS.md`.",
            "See `docs/PRODUCTION_PLAN.md` §7 and `docs/CROSS_PLATFORM.md`.",
        ),
    ]
    for title, detail, body in steps:
        with st.container(border=True):
            st.markdown(f"**{title}**")
            st.markdown(detail)
            st.markdown(body)

    st.markdown("#### Quick checks")
    st.markdown(
        f"- **PermissionError on `/data`:** wrappers run `init-data` first; "
        "or `docker compose --profile init run --rm init-data`.\n"
        f"- **`snakemake` not found on the host:** use `{wrap} --entrypoint snakemake …`.\n"
        "- **Generated configs:** `samples.tsv` and `label.auto.yaml` are gitignored; "
        "only `samples.tsv.example` is committed.\n"
        "- **Simulated dress rehearsal (no download):** inside the container, "
        "`python scripts/pilot_large.py` or the CI 200 kb pilot in `docs/CROSS_PLATFORM.md`.\n"
        "- **Apple Silicon:** if a package lacks arm64, build with "
        "`docker compose build --platform linux/amd64 cli` (slower under emulation)."
    )


def tab_status(results: Path) -> None:
    repo = data.find_repo_root(os.environ.get("ARCSV2_REPO") or None)
    if repo is None:
        st.info("Project status needs the source checkout. Run the dashboard from the repo, or set `ARCSV2_REPO`.")
        return
    table = data.milestones(repo, results)
    done = int(table["done"].sum())
    st.progress(done / len(table) if len(table) else 0.0, text=f"{done} of {len(table)} tracked items done (auto-detected from the repository)")
    table["status"] = table["done"].map({True: "✅ done", False: "⏳ open"})
    st.dataframe(table[["milestone", "item", "status", "evidence"]], hide_index=True, width="stretch")
    reviews = data.review_documents(repo)
    if reviews:
        st.markdown("#### Review history")
        titles = [title for title, _ in reviews]
        pick = st.selectbox("Review", titles, index=len(titles) - 1)
        path = dict(reviews)[pick]
        text = path.read_text(encoding="utf-8")
        # Step headings down two levels so the review sits under this section.
        text = "\n".join("##" + line if line.startswith("#") else line for line in text.splitlines())
        with st.container(border=True):
            st.markdown(text)


# --- page -------------------------------------------------------------------


def main() -> None:
    results, is_demo = results_folder()
    st.title("ARC-SV2")
    st.caption(f"Complex structural-variant caller with an ML filter · version {__version__}")

    with st.sidebar:
        st.header("Data")
        st.markdown(f"**Results folder**  \n`{results}`")
        if is_demo:
            st.info("Showing the bundled simulation demo. Use `arcsv2 dashboard --results <folder>` for your own runs.")
        runs = data.discover_runs(results)
        if not runs:
            st.error("No `arcsv_out.tab` or `merged.tab` found under this folder.")
            st.stop()
        names = [r.name for r in runs]
        default = max(range(len(runs)), key=lambda i: (runs[i].meta.get("coverage") or 0, i))
        run = runs[names.index(st.selectbox("Run", names, index=default))]
        if run.meta:
            st.markdown("**Run details**")
            details = {k: v for k, v in run.meta.items() if isinstance(v, (int, float, str)) and k != "name"}
            st.dataframe(
                pd.DataFrame({"field": list(details), "value": [str(v) for v in details.values()]}),
                hide_index=True,
                width="stretch",
            )

    calls, events = load_run(run)
    tabs = st.tabs(
        ["Guide", "Overview", "How it works", "Call explorer", "ML & labels", "Project status"]
    )
    with tabs[0]:
        tab_guide()
    with tabs[1]:
        tab_overview(runs, run, calls, events)
    with tabs[2]:
        tab_how(calls)
    with tabs[3]:
        tab_explorer(calls, events)
    with tabs[4]:
        tab_ml(results)
    with tabs[5]:
        tab_status(results)


main()
