"""Command line. ``--mode paper`` is the default and matches §4.12."""

from __future__ import annotations

from pathlib import Path

import typer

from arcsv2 import __version__
from arcsv2.beds import read_bed, read_bed_by_chrom
from arcsv2.config import CallConfig
from arcsv2.pipeline import call_genome, call_sam, fasta_length, parse_region

app = typer.Typer(add_completion=False, no_args_is_help=True)


def version_callback(value: bool) -> None:
    if value:
        typer.echo(__version__)
        raise typer.Exit()


@app.callback()
def main(version: bool = typer.Option(False, "--version", callback=version_callback, is_eager=True)) -> None:
    """ARC-SV2 structural-variant caller."""


@app.command("call")
def call(
    inputs: Path = typer.Option(..., "--input", "-i", help="SAM, BAM, or CRAM"),
    reference: Path = typer.Option(..., "--reference", "-R"),
    gaps: Path = typer.Option(..., "--gaps", "-G"),
    out: Path = typer.Option(..., "--out", "-o"),
    region: str = typer.Option("", "--region"),
    threads: int = typer.Option(1, "--threads"),
    seed: int = typer.Option(204178949, "--seed"),
    mode: str = typer.Option("paper", "--mode"),
    max_paths: int = typer.Option(1000, "--max-paths"),
    max_mes_extra: int = typer.Option(0, "--max-mes-extra"),
    filter_high_nm: bool = typer.Option(
        False, "--filter-high-nm", help="C5: drop alignments with NM > --max-nm"
    ),
    max_nm: int = typer.Option(5, "--max-nm", help="Edit-distance cap when --filter-high-nm"),
    filter_both_mates_clipped: bool = typer.Option(
        False,
        "--filter-both-mates-clipped",
        help="C5: drop pairs where both mates are soft-clipped",
    ),
    columnar_pairs: bool = typer.Option(
        True,
        "--columnar-pairs/--no-columnar-pairs",
        help="M2: store pairs as numpy SoA (default on; --no-columnar-pairs for list fallback)",
    ),
    engine: str = typer.Option(
        "rust",
        "--engine",
        help="BAM/graph engine: rust (requires arcsv2_native) or python",
    ),
    debug: bool = typer.Option(False, "--debug", help="Show Python traceback on failure"),
) -> None:
    """Call SVs in one region. Paper defaults are on unless a flag changes them."""
    if threads < 1:
        raise typer.BadParameter("--threads must be at least 1")
    if engine not in ("python", "rust"):
        raise typer.BadParameter("--engine must be python or rust")
    from arcsv2.input_checks import InputCheckError, validate_call_inputs

    try:
        validate_call_inputs(inputs, reference, gaps, region=region)
    except InputCheckError as exc:
        if debug:
            raise
        raise typer.BadParameter(str(exc)) from None
    cfg = CallConfig(
        seed=seed,
        mode=mode,
        max_paths=max_paths,
        max_mes_extra=max_mes_extra,
        threads=threads,
        filter_high_nm=filter_high_nm,
        max_nm=max_nm,
        filter_both_mates_clipped=filter_both_mates_clipped,
        columnar_pairs=columnar_pairs,
        engine=engine,
    )
    chrom, start, end = parse_region(region)
    dest = out / "arcsv_out.tab" if out.suffix == "" else out
    if dest.suffix == "":
        dest.mkdir(parents=True, exist_ok=True)
        dest = dest / "arcsv_out.tab"
    chrom_len = fasta_length(str(reference), chrom) if reference.is_file() else None
    gap_list = read_bed(gaps, chrom) if gaps.is_file() else []
    rows = call_sam(
        str(inputs),
        str(dest),
        cfg,
        chrom=chrom,
        chrom_len=chrom_len,
        gaps=gap_list,
        start=start,
        end=end,
        reference=str(reference) if reference.is_file() else None,
    )
    vcf_gz = dest.with_name("arcsv_out.vcf.gz")
    extra = f" and {vcf_gz}" if vcf_gz.is_file() else ""
    typer.echo(f"wrote {len(rows)} calls to {dest}{extra}")


@app.command("call-genome")
def call_genome_cmd(
    manifest: Path = typer.Option(..., "--manifest", help="TSV chrom<path"),
    out: Path = typer.Option(..., "--out", "-o"),
    jobs: int = typer.Option(1, "--jobs"),
    threads_per_job: int = typer.Option(1, "--threads-per-job"),
) -> None:
    """Chromosome jobs, largest first. An existing tab is left in place."""
    del threads_per_job
    paths: dict[str, str] = {}
    for line in manifest.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        chrom, path = line.split("\t")[:2]
        paths[chrom] = path
    written = call_genome(paths, str(out), jobs=jobs)
    typer.echo(f"finished {len(written)} chromosomes")


@app.command("run")
def run_cmd(
    inputs: Path = typer.Option(..., "--input", "-i", help="CRAM, BAM, or SAM"),
    reference: str = typer.Option(
        "GRCh38",
        "--reference",
        "-R",
        help="GRCh38 (download/cache) or path to reference FASTA",
    ),
    out: Path = typer.Option(..., "--out", "-o", help="Output directory (sample subfolder created)"),
    chroms: str = typer.Option("", "--chroms", help="Comma/space list; default = primary autosomes + chrX"),
    threads: str = typer.Option("auto", "--threads", help="Thread count or auto"),
    seed: int = typer.Option(204178949, "--seed"),
    engine: str = typer.Option("rust", "--engine", help="rust or python"),
    debug: bool = typer.Option(False, "--debug", help="Show Python traceback on failure"),
) -> None:
    """Whole-sample SV calling: extract chromosomes, call, merge-filter, summarize."""
    if engine not in ("python", "rust"):
        raise typer.BadParameter("--engine must be python or rust")
    from arcsv2.input_checks import InputCheckError
    from arcsv2.run import print_run_summary, run_sample

    try:
        meta = run_sample(
            reads=inputs,
            reference=reference,
            out_dir=out,
            chroms=chroms or None,
            threads=threads,
            seed=seed,
            engine=engine,
        )
    except InputCheckError as exc:
        if debug:
            raise
        raise typer.BadParameter(str(exc)) from None
    except Exception as exc:
        if debug:
            raise
        raise typer.BadParameter(str(exc)) from None
    print_run_summary(meta)


@app.command("batch")
def batch_cmd(
    manifest: Path = typer.Option(..., "--manifest", help="TSV: sample, input [, reference]"),
    out: Path = typer.Option(..., "--out", "-o", help="Batch output directory"),
    reference: str = typer.Option("GRCh38", "--reference", "-R", help="Default reference when manifest omits it"),
    chroms: str = typer.Option("", "--chroms", help="Chromosome list (default primary + chrX)"),
    threads: str = typer.Option("auto", "--threads"),
    seed: int = typer.Option(204178949, "--seed"),
    engine: str = typer.Option("rust", "--engine"),
    retries: int = typer.Option(3, "--retries", help="Retries per sample for transient errors"),
    debug: bool = typer.Option(False, "--debug"),
) -> None:
    """Run ``arcsv2 run`` per manifest row; failures are isolated with FAILED.txt + status table."""
    from arcsv2.batch import run_batch
    from arcsv2.input_checks import InputCheckError

    try:
        summary = run_batch(
            manifest,
            out,
            reference_default=reference,
            chroms=chroms or None,
            threads=threads,
            seed=seed,
            engine=engine,
            max_retries=retries,
        )
    except InputCheckError as exc:
        if debug:
            raise
        raise typer.BadParameter(str(exc)) from None
    except Exception as exc:
        if debug:
            raise
        raise typer.BadParameter(str(exc)) from None
    n_ok = len(summary.get("done", []))
    n_bad = len(summary.get("failed", []))
    typer.echo(f"batch finished: {n_ok} done, {n_bad} failed (see {out / 'batch_status.tsv'})")
    if n_bad:
        raise typer.Exit(code=1)


@app.command("merge-filter")
def merge_filter_cmd(
    tables: list[Path] = typer.Argument(..., help="arcsv_out.tab files"),
    out: Path = typer.Option(..., "--out", "-o"),
    min_size: int = typer.Option(50, "--min-size"),
    simple_repeats: Path | None = typer.Option(None, "--simple-repeats"),
    segdups: Path | None = typer.Option(None, "--segdups"),
    blacklist: Path | None = typer.Option(None, "--blacklist"),
    gaps: Path | None = typer.Option(None, "--gaps"),
) -> None:
    """Apply the paper merge-filter (§4.9) to one or more call tables."""
    from arcsv2.config import OUTPUT_COLUMNS
    from arcsv2.merge_filter import inslen_is_lower_bound, merge_calls, passes_filter

    repeat_iv = read_bed_by_chrom(simple_repeats) if simple_repeats and simple_repeats.is_file() else {}
    segdup_iv = read_bed_by_chrom(segdups) if segdups and segdups.is_file() else {}
    black_iv = read_bed_by_chrom(blacklist) if blacklist and blacklist.is_file() else {}
    gap_iv = read_bed_by_chrom(gaps) if gaps and gaps.is_file() else {}

    records = []
    for table in tables:
        lines = table.read_text(encoding="utf-8").splitlines()
        if not lines:
            continue
        header = lines[0].split("\t")
        for line in lines[1:]:
            row = dict(zip(header, line.split("\t")))
            start = int(row["minbp"])
            end = int(row["maxbp"])
            affected = int(float(row["len_affected"]))
            chrom = row["chrom"]
            ok, _reason = passes_filter(
                start,
                end,
                affected,
                simple_repeats=repeat_iv.get(chrom, []),
                segdups=segdup_iv.get(chrom, []),
                blacklist=black_iv.get(chrom, []),
                gaps=gap_iv.get(chrom, []),
                min_affected=min_size,
                length_is_lower_bound=inslen_is_lower_bound(row.get("inslen")),
            )
            if not ok:
                continue
            records.append(
                {
                    "chrom": row["chrom"],
                    "start": start,
                    "end": end,
                    "sv_class": row.get("complextype") if row.get("complextype") not in ("", "NA") else row.get("svtype"),
                    "gt": row.get("gt", "HET"),
                    "sample": table.parent.name,
                    "row": row,
                }
            )
    merged = merge_calls(records, n_samples=max(1, len(tables)))
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(OUTPUT_COLUMNS) + "\n")
        for rec in merged:
            row = dict(rec["row"])
            row["af"] = f"{rec['allele_frequency']:.6f}"
            handle.write("\t".join(row.get(col, "NA") for col in OUTPUT_COLUMNS) + "\n")
    typer.echo(f"wrote {len(merged)} merged calls to {out}")


@app.command("ml-label-hits")
def ml_label_hits_cmd(
    hits: Path = typer.Option(..., "--hits", help="TSV q_start q_end q_len mapq mlen blen cigar"),
    out: Path = typer.Option(..., "--out", "-o"),
) -> None:
    """Apply the paper alignment criteria (§4.10) to a table of minimap2 hits."""
    from arcsv2.label import is_validated

    lines = hits.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join([*header, "label", "reason"]) + "\n")
        for line in lines[1:]:
            if not line:
                continue
            row = dict(zip(header, line.split("\t")))
            ok, reason = is_validated(
                q_start=int(row["q_start"]),
                q_end=int(row["q_end"]),
                q_len=int(row["q_len"]),
                mapq=int(row["mapq"]),
                mlen=int(row["mlen"]),
                blen=int(row["blen"]),
                cigar=row.get("cigar", ""),
            )
            label = "1" if ok else ("nan" if reason == "novel_insertion" else "0")
            handle.write(line + f"\t{label}\t{reason}\n")
    typer.echo(f"wrote labels to {out}")


@app.command("ml-label")
def ml_label_cmd(
    calls: Path = typer.Option(..., "--calls"),
    reference: Path = typer.Option(..., "--reference", "-R"),
    hap1: Path = typer.Option(..., "--hap1"),
    hap2: Path = typer.Option(..., "--hap2"),
    out: Path = typer.Option(..., "--out", "-o"),
    sample: str = typer.Option(..., "--sample"),
    hap1_index: Path | None = typer.Option(None, "--hap1-index"),
    hap2_index: Path | None = typer.Option(None, "--hap2-index"),
    unreliable: Path | None = typer.Option(None, "--unreliable"),
    flank: int = typer.Option(150, "--flank"),
    threads: int = typer.Option(1, "--threads"),
    sequential: bool = typer.Option(
        False,
        "--sequential",
        help="Load one haplotype index at a time (about half the RAM).",
    ),
) -> None:
    """Label one sample. Two whole-genome indexes need well over 16 GB of RAM.

    Hits on either strand count. A heterozygote is label 1 when the alt passes
    on one haplotype and the reference passes on the other. ``--sequential``
    maps against hap1 then hap2 so only one index sits in memory.
    """
    from arcsv2.atomic_io import atomic_write_text
    from arcsv2.beds import read_bed_by_chrom
    from arcsv2.input_checks import InputCheckError
    from arcsv2.label import LABEL_COLUMNS, HapAligner, label_rows
    from arcsv2.mmi_validate import check_no_hit_fraction, validate_haplotype_indexes

    if unreliable is not None and not unreliable.is_file():
        raise typer.BadParameter(f"--unreliable was given but not found: {unreliable}")
    try:
        validate_haplotype_indexes(hap1, hap2, hap1_index, hap2_index, reference)
    except InputCheckError as exc:
        raise typer.BadParameter(str(exc)) from None
    lines = calls.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    rows = [dict(zip(header, line.split("\t"))) for line in lines[1:] if line]
    aligner = HapAligner(
        hap1, hap2, hap1_index, hap2_index, threads=threads, sequential=sequential
    )
    beds = read_bed_by_chrom(unreliable) if unreliable is not None else None
    labeled = label_rows(rows, reference, aligner, sample, beds, flank)
    try:
        check_no_hit_fraction(labeled)
    except InputCheckError as exc:
        failed = out.parent / "FAILED.txt"
        from arcsv2.batch import write_failed_txt

        write_failed_txt(out.parent, str(exc))
        raise typer.BadParameter(f"{exc} (also wrote {failed})") from None
    columns = list(header)
    for name in LABEL_COLUMNS:
        if name not in columns:
            columns.append(name)
    lines_out = ["\t".join(columns)]
    for row in labeled:
        lines_out.append("\t".join(row.get(col, "NA") for col in columns))
    atomic_write_text(out, "\n".join(lines_out) + "\n")
    typer.echo(f"wrote {len(labeled)} labels to {out}")


@app.command("mmi-validate")
def mmi_validate_cmd(
    mmi: Path = typer.Argument(..., help="minimap2 haplotype index (.mmi)"),
    fasta: Path | None = typer.Option(
        None,
        "--fasta",
        help="FASTA this index was built from. Enables the size check (HG00140).",
    ),
    reference: Path | None = typer.Option(
        None,
        "--reference",
        "-R",
        help="Reference FASTA for a one-locus smoke alignment.",
    ),
) -> None:
    """Fail if a haplotype index is missing, too small, or fails a smoke alignment.

    ``ml-label`` runs the same checks on both haplotypes before mapping.
    """
    from arcsv2.input_checks import InputCheckError
    from arcsv2.mmi_validate import validate_mmi_size, validate_mmi_smoke

    try:
        if fasta is not None:
            validate_mmi_size(mmi, fasta, hap_name=mmi.name)
        if reference is not None:
            validate_mmi_smoke(mmi, reference, hap_name=mmi.name)
        if fasta is None and reference is None and not mmi.is_file():
            raise InputCheckError(f"Haplotype index missing: {mmi}", "Rebuild with minimap2 -d.")
    except InputCheckError as exc:
        raise typer.BadParameter(str(exc)) from None
    typer.echo(f"ok {mmi}")


@app.command("ml-dataset")
def ml_dataset_cmd(
    tables: list[Path] = typer.Argument(...),
    out: Path = typer.Option(..., "--out", "-o"),
) -> None:
    """Concatenate labeled tables, drop NA labels, and write a summary JSON."""
    from arcsv2.label import dataset_tables, write_dataset

    parsed: list[list[dict[str, str]]] = []
    for table in tables:
        lines = table.read_text(encoding="utf-8").splitlines()
        if not lines:
            continue
        header = lines[0].split("\t")
        parsed.append([dict(zip(header, line.split("\t"))) for line in lines[1:] if line])
    kept, summary = dataset_tables(parsed)
    write_dataset(kept, summary, out)
    typer.echo(f"wrote {len(kept)} labeled calls to {out}")


@app.command("ml-annotate")
def ml_annotate_cmd(
    calls: Path = typer.Option(..., "--calls"),
    out: Path = typer.Option(..., "--out", "-o"),
    segdups: Path | None = typer.Option(None, "--segdups"),
    simple_repeats: Path | None = typer.Option(None, "--simple-repeats"),
    low_map: Path | None = typer.Option(None, "--mappability"),
) -> None:
    """Add Tier C overlap fractions (§7.2) to a call table."""
    from arcsv2.annotate import annotate_span

    lines = calls.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    extra = [
        "frac_segdup",
        "frac_simple_repeat",
        "frac_low_mappability",
        "frac_SINE",
        "frac_LINE",
        "frac_LTR",
        "frac_satellite",
        "bp_overlap_segdup",
    ]
    seg = read_bed_by_chrom(segdups) if segdups and segdups.is_file() else {}
    repeats = read_bed_by_chrom(simple_repeats) if simple_repeats and simple_repeats.is_file() else {}
    low = read_bed_by_chrom(low_map) if low_map and low_map.is_file() else {}
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join([*header, *extra]) + "\n")
        for line in lines[1:]:
            if not line:
                continue
            row = dict(zip(header, line.split("\t")))
            chrom = row.get("chrom", "")
            feats = annotate_span(
                int(row["minbp"]),
                int(row["maxbp"]),
                segdups=seg.get(chrom, []),
                simple_repeats=repeats.get(chrom, []),
                low_mappability=low.get(chrom, []),
            )
            handle.write(line + "\t" + "\t".join(f"{feats[name]:.6f}" for name in extra) + "\n")
    typer.echo(f"wrote annotations to {out}")


@app.command("ml-train")
def ml_train_cmd(
    out: Path = typer.Option(..., "--out", "-o"),
    labels: Path | None = typer.Option(None, "--labels", help="TSV class label score"),
) -> None:
    """Write a manifest. A call table with ``label`` fits per-class logistics.

    A three-column ``class label score`` file only sets thresholds on the paper
    presets. Fitted coefficients and the training PCA are stored in the manifest.
    """
    from arcsv2.config import FEATURE_SCHEMA
    from arcsv2.train import PAPER_PRESETS, choose_threshold, fit_labeled_calls, write_manifest

    out.mkdir(parents=True, exist_ok=True)
    if labels and labels.is_file():
        text = labels.read_text(encoding="utf-8").splitlines()
        header = text[0].split("\t") if text else []
        if "label" in header and "svtype" in header:
            rows = [dict(zip(header, line.split("\t"))) for line in text[1:] if line]
            payload = fit_labeled_calls(rows)
            write_manifest(out / "manifest.json", payload)
            typer.echo(f"wrote {out / 'manifest.json'}")
            return
        grouped: dict[str, tuple[list[int], list[float]]] = {}
        for line in text[1:]:
            kind, label_s, score_s = line.split("\t")[:3]
            ys, scores = grouped.setdefault(kind, ([], []))
            ys.append(int(label_s))
            scores.append(float(score_s))
        thresholds: dict[str, float] = {}
        for kind, (ys, scores) in grouped.items():
            thresholds[kind] = choose_threshold(ys, scores).threshold
    else:
        thresholds = {}
    payload = {
        "feature_schema": FEATURE_SCHEMA,
        "fitted": False,
        "models": {kind: preset["kind"] for kind, preset in PAPER_PRESETS.items()},
        "presets": PAPER_PRESETS,
        "thresholds": thresholds,
    }
    write_manifest(out / "manifest.json", payload)
    typer.echo(f"wrote {out / 'manifest.json'}")


@app.command("ml-predict")
def ml_predict_cmd(
    model: Path = typer.Option(..., "--model"),
    calls: Path = typer.Option(..., "--calls"),
    out: Path = typer.Option(..., "--out", "-o"),
) -> None:
    """Score calls with a manifest from ``ml-train`` (§7.4)."""
    from arcsv2.predict import predict_rows
    from arcsv2.train import load_manifest

    manifest_path = model / "manifest.json" if model.is_dir() else model
    manifest = load_manifest(manifest_path)
    lines = calls.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    rows = [dict(zip(header, line.split("\t"))) for line in lines[1:] if line]
    scored = predict_rows(rows, manifest, None)
    extra = ["sv_class", "ml_prob", "ml_pass", "ml_model"]
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join([*header, *extra]) + "\n")
        for row in scored:
            handle.write("\t".join(row.get(col, "NA") for col in [*header, *extra]) + "\n")
    typer.echo(f"wrote {len(scored)} scored calls to {out}")


@app.command("dashboard")
def dashboard_cmd(
    results: Path | None = typer.Option(None, "--results", "-r", help="Folder with arcsv_out.tab runs. Default: bundled demo."),
    port: int = typer.Option(8501, "--port"),
    host: str = typer.Option("localhost", "--host", help="Use 0.0.0.0 inside a container."),
    legacy: bool = typer.Option(False, "--legacy", help="Open the developer dashboard instead of ARC-SV Studio."),
) -> None:
    """Open ARC-SV Studio (default) or the legacy developer dashboard. Needs arcsv2[dashboard]."""
    import os
    import subprocess
    import sys

    try:
        import streamlit  # noqa: F401
    except ImportError as exc:
        raise typer.BadParameter("The dashboard needs Streamlit. Install arcsv2[dashboard].") from exc
    dash = Path(__file__).resolve().parent / "dashboard"
    app_path = dash / ("app.py" if legacy else "studio_app.py")
    env = dict(os.environ)
    if results is not None:
        if not results.exists():
            raise typer.BadParameter(f"results folder not found: {results}")
        env["ARCSV2_RESULTS"] = str(results.resolve())
    command = [
        sys.executable, "-m", "streamlit", "run", str(app_path),
        "--server.port", str(port),
        "--server.address", host,
        "--server.headless", "true" if host != "localhost" else "false",
        "--browser.gatherUsageStats", "false",
        "--theme.primaryColor", "#2a78d6",
    ]
    raise typer.Exit(subprocess.call(command, env=env))


if __name__ == "__main__":
    app()
