"""Streamlit page bodies for ARC-SV Studio."""

from __future__ import annotations

import json
import os
from pathlib import Path

import pandas as pd
import streamlit as st

from arcsv2 import __version__
from arcsv2.dashboard import data
from arcsv2.dashboard.studio import checks, example, jobs, progress


def results_path() -> Path:
    raw = os.environ.get("ARCSV2_RESULTS", "").strip()
    if raw:
        return Path(raw).expanduser().resolve()
    return Path(__file__).resolve().parent / "demo"


def page_setup() -> None:
    st.header("Setup check")
    st.caption("Green means you are ready. Red rows include a fix hint.")
    root = results_path()
    for item in checks.run_all_checks(root):
        icon = "✅" if item.ok else "❌"
        with st.container(border=True):
            st.markdown(f"{icon} **{item.name}** — {item.detail}")
            if not item.ok and item.fix:
                st.info(f"**Fix:** {item.fix}")

    st.subheader("Download reference (GRCh38)")
    st.caption("Needed once for real CRAM/BAM runs. About 3 GB. The example run uses a tiny bundled reference and does not need this.")
    prog = st.progress(0, text="Idle")
    status = st.empty()
    if st.button("Download reference (GRCh38)", type="primary"):
        try:
            cache = Path(os.environ.get("ARCSV2_CACHE", "") or (root / "ref_cache"))

            def on_progress(done: int, total: int | None) -> None:
                if total and total > 0:
                    pct = min(100, int(100 * done / total))
                    prog.progress(pct / 100.0, text=f"Downloading… {pct}% ({done // (1024**2)} / {total // (1024**2)} MiB)")
                else:
                    prog.progress(min(0.99, done / (3.2 * 1024**3)), text=f"Downloading… {done // (1024**2)} MiB")

            fasta, gaps = checks.download_grch38(cache_dir=cache, on_progress=on_progress)
            prog.progress(1.0, text="Done")
            status.success(f"Reference ready: `{fasta}` (gaps: `{gaps.name}`)")
            st.rerun()
        except Exception as exc:
            status.error(str(exc))


def page_example() -> None:
    st.header("Try the example")
    st.markdown(
        "Runs **`arcsv2 run`** on a small bundled BAM (about 100 kb of simulated reads). "
        "Expect **2–5 minutes**. When it finishes, open **Results** to explore the calls."
    )
    root = results_path()
    ready = example.example_results_ready(root)
    if ready:
        st.success(f"Example results are ready under `{ready}`.")
        if st.button("Open in Results", type="primary"):
            st.session_state["studio_results_focus"] = str(ready)
            st.info("Go to **Results** in the sidebar — the example run is selected when available.")
    running = any(j.get("kind") == "example" and not j.get("cancelled") for j in jobs.load_jobs(root))
    if st.button("Run example (2–5 min)", type="primary", disabled=bool(ready) and not st.checkbox("Re-run even if results exist", value=False)):
        with st.spinner("Starting example… preparing BAM if needed"):
            job = example.start_example_run(root)
        st.success(f"Example started (`{job['id']}`). Watch progress under **My analyses**.")
        st.caption(f"Log: `{job['log']}`")
        st.rerun()
    if running and not ready:
        st.info("Example is running. Refresh this page or open **My analyses** for progress.")


def page_new_analysis() -> None:
    st.header("New analysis")
    root = results_path()
    files = jobs.list_read_files(root)
    if not files:
        st.warning(
            f"No BAM/CRAM files found in your data folder yet. "
            f"Copy read files into **{root}** (the folder you chose when starting Studio), then click Refresh."
        )
        if st.button("Refresh file list"):
            st.rerun()
        return
    labels: list[str] = []
    for p in files:
        try:
            labels.append(p.relative_to(root).as_posix())
        except ValueError:
            labels.append(p.name)
    picked = st.multiselect("Select input files", labels, default=labels[:1])
    mode = st.radio("Mode", ["Quick test (chr20)", "Whole genome (chr1–22 + X)"], horizontal=True)
    chroms = "chr20" if mode.startswith("Quick") else ""
    total_bytes = sum(files[labels.index(x)].stat().st_size for x in picked) if picked else 0
    time_est, disk_est = jobs.estimate_run(chroms, max(1, len(picked)), total_bytes // max(1, len(picked) or 1))
    st.markdown(f"**Estimated time:** {time_est}  \n**Estimated disk:** {disk_est}")
    if st.button("Start analysis", type="primary", disabled=not picked):
        batch_root = root / "analyses" / jobs._utc().replace(":", "").replace("-", "")[:15]
        manifest = batch_root / "manifest.tsv"
        rows = []
        for label in picked:
            path = files[labels.index(label)]
            sample = path.stem.replace(".cram", "").replace(".bam", "")
            rows.append((sample, path))
        jobs.write_batch_manifest(rows, manifest)
        job = jobs.start_batch_job(manifest=manifest, out_dir=batch_root, chroms=chroms)
        st.success(f"Started analysis `{job['id']}`. Track progress under **My analyses**.")
        st.rerun()


def page_analyses() -> None:
    st.header("My analyses")
    root = results_path()
    rows = progress.list_analysis_rows(root)
    if not rows:
        st.info("No analyses yet. Start one from **Try example** or **New analysis**.")
        return
    frame = pd.DataFrame(rows)
    show = frame[["sample", "stage", "percent", "eta_utc", "status", "message", "batch_dir", "job_id"]].copy()
    show = show.rename(
        columns={
            "stage": "Stage",
            "percent": "%",
            "eta_utc": "ETA (UTC)",
            "status": "Status",
            "message": "Note",
            "sample": "Sample",
        }
    )
    st.dataframe(show.drop(columns=["batch_dir", "job_id"], errors="ignore"), hide_index=True, width="stretch")

    for row in rows:
        sample = row["sample"]
        batch_dir = Path(row["batch_dir"])
        failed = batch_dir / sample / "FAILED.txt"
        if row.get("stage") == "failed" and failed.is_file():
            st.error(f"**{sample}:** {failed.read_text(encoding='utf-8').strip()}")

    job_ids = sorted({r["job_id"] for r in rows if r.get("job_id")})
    col_a, col_b = st.columns(2)
    with col_a:
        pick_cancel = st.selectbox("Cancel job", [""] + job_ids, key="cancel_job")
        if st.button("Cancel", disabled=not pick_cancel) and pick_cancel:
            msg = progress.cancel_job(pick_cancel, root)
            st.warning(msg)
            st.rerun()
    with col_b:
        batch_dirs = sorted({r["batch_dir"] for r in rows})
        pick_resume = st.selectbox("Resume batch folder", [""] + batch_dirs, key="resume_batch")
        if st.button("Resume", disabled=not pick_resume) and pick_resume:
            job = progress.resume_batch(Path(pick_resume))
            st.success(f"Resumed as `{job['id']}`")
            st.rerun()

    if st.button("Refresh progress"):
        st.rerun()

    with st.expander("Details (technical)", expanded=False):
        st.json(jobs.load_jobs(root))


def _size_bin_counts(calls: pd.DataFrame) -> dict[str, int]:
    bins = {"50–499 bp": 0, "500–4,999 bp": 0, "5,000–49,999 bp": 0, "≥50,000 bp": 0}
    for length in calls["len_affected"].fillna(0):
        try:
            n = int(float(length))
        except (TypeError, ValueError):
            continue
        if n < 500:
            bins["50–499 bp"] += 1
        elif n < 5000:
            bins["500–4,999 bp"] += 1
        elif n < 50000:
            bins["5,000–49,999 bp"] += 1
        else:
            bins["≥50,000 bp"] += 1
    return bins


def _render_qc_table(payload: dict) -> None:
    flags = payload.get("flags") or []
    rows = [
        {"Item": "Sample", "Value": payload.get("sample", "—")},
        {"Item": "Runtime (s)", "Value": payload.get("runtime_sec", "—")},
        {"Item": "Median read length (bp)", "Value": payload.get("median_read_length", "—")},
        {"Item": "Median insert size (bp)", "Value": payload.get("median_insert_size", "—")},
        {"Item": "Coverage (×)", "Value": payload.get("coverage_x", "—")},
        {"Item": "Merged calls", "Value": payload.get("merged_calls", "—")},
        {"Item": "Flags", "Value": "; ".join(flags) if flags else "None"},
    ]
    st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")
    by_type = payload.get("calls_by_type") or {}
    if by_type:
        st.markdown("**Calls by type (QC)**")
        st.dataframe(
            pd.DataFrame({"Type": list(by_type), "Count": list(by_type.values())}),
            hide_index=True,
        )


def page_results() -> None:
    from arcsv2.dashboard.app import haplotype_chart, load_run

    st.header("Results")
    root = results_path()
    runs = data.discover_runs(root)
    if not runs:
        st.info("No results yet. Try the **example** or run an analysis.")
        return
    names = [r.name for r in runs]
    focus = st.session_state.get("studio_results_focus", "")
    default_idx = 0
    if focus:
        for i, r in enumerate(runs):
            if focus in str(r.calls_path) or focus in r.name:
                default_idx = i
                break
    run = runs[names.index(st.selectbox("Run", names, index=default_idx))]
    calls, _events = load_run(run)

    st.subheader("Summary")
    c1, c2, c3 = st.columns(3)
    c1.metric("Total calls", int(len(calls)))
    by_type = calls.groupby("call_type").size().to_dict()
    type_bits = ", ".join(f"{k}: {v}" for k, v in sorted(by_type.items()))
    c2.metric("By type", len(by_type), help=type_bits)
    size_bins = _size_bin_counts(calls)
    c3.metric("Size bins", sum(1 for v in size_bins.values() if v), help=", ".join(f"{k}: {v}" for k, v in size_bins.items()))
    st.dataframe(
        pd.DataFrame({"Type": list(by_type), "Count": list(by_type.values())}),
        hide_index=True,
    )
    st.dataframe(
        pd.DataFrame({"Size bin": list(size_bins), "Count": list(size_bins.values())}),
        hide_index=True,
    )

    st.subheader("Calls")
    all_types = sorted(calls["call_type"].dropna().unique().tolist())
    filt_type = st.multiselect("Type", all_types, default=all_types)
    max_len = int(calls["len_affected"].max() or 5000)
    min_len = st.slider("Minimum affected length (bp)", 0, max(max_len, 50), 50)
    view = calls[(calls["call_type"].isin(filt_type)) & (calls["len_affected"] >= min_len)]
    st.dataframe(
        view[
            [
                "call_id",
                "chrom",
                "minbp",
                "maxbp",
                "call_type",
                "gt",
                "len_affected",
                "score_vs_ref",
            ]
        ],
        hide_index=True,
        height=320,
    )

    qc = run.calls_path.parent / "qc" / "report.json"
    if not qc.is_file():
        # sample-level qc next to merged.tab parent
        alt = run.calls_path.parent.parent / "qc" / "report.json"
        if alt.is_file():
            qc = alt
    if qc.is_file():
        st.subheader("QC report")
        try:
            _render_qc_table(json.loads(qc.read_text(encoding="utf-8")))
        except json.JSONDecodeError:
            st.warning("QC report could not be read.")

    st.subheader("Single variant")
    if view.empty:
        return
    cid = st.selectbox("Call", view["call_id"].tolist())
    row = view[view["call_id"] == cid].iloc[0]
    st.markdown(f"**In plain language:** {data.rearrangement_plain(row)}")
    st.altair_chart(haplotype_chart(row), use_container_width=True)
    merged = run.calls_path
    vcf = merged.parent / "arcsv_out.vcf.gz"
    if not vcf.is_file() and merged.name == "merged.tab":
        # per-chrom vcf under first chrom folder
        for child in merged.parent.iterdir():
            cand = child / "arcsv_out.vcf.gz"
            if cand.is_file():
                vcf = cand
                break
    if vcf.is_file():
        st.download_button("Download VCF", vcf.read_bytes(), file_name=vcf.name)
    st.download_button("Download .tab", merged.read_bytes(), file_name=merged.name)
    st.download_button("Download CSV", view.to_csv(index=False).encode("utf-8"), file_name="calls.csv")


def page_help() -> None:
    st.header("Help")
    st.markdown(
        f"""
### What is ARC-SV Studio?
A guided app for calling **structural variants** (large rearrangements) from CRAM/BAM files.
Version **{__version__}**.

### Quick glossary
| Term | Meaning |
|:--|:--|
| **SV** | DNA segment that is deleted, duplicated, inverted, or moved |
| **CRAM/BAM** | Aligned sequencing reads |
| **Genotype (GT)** | HET = one copy changed; HOM = both copies changed |
| **Score** | Higher = more confident call |

### FAQ
**Does closing the browser stop my run?** No — analyses started here run in the background.

**Where are files stored?** In the data folder you chose at launch (default **Documents/ARC-SV**).

**Command line instead?** See `arcsv2 run` and `docs/CLUSTER.md` for SLURM.

**Something failed?** Open **My analyses** for the plain-language error, or see `docs/TROUBLESHOOTING.md`.
"""
    )
