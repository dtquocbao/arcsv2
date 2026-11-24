# Label calls from haplotype assemblies.
# Run this inside the Linux container (docs/CROSS_PLATFORM.md). Paths below
# are container paths. Fill workflows/config/label.yaml before a pilot.
# Inside the container the tools are already on PATH, so omit --use-conda.
# Outside it, `snakemake --use-conda` builds environment.yml per rule.
# fetch_region sorts and indexes the BAM; pysam fetch requires the .bai.
#
# Set data_root: /data in label.yaml (or label.auto.yaml) so haps/work/results/
# labels live on the Docker volume. Pilot configs keep data_root: "." .
#
# Paper merge-filter BEDs (Methods p. 27): set simple_repeats / segdups /
# blacklist / gaps in the config (fetch with scripts/fetch_filter_beds.py).
# Leave those keys empty to keep min-size-only filtering (with-vs-without).
#
# C4: call is one Snakemake job per (sample, chrom); merge joins chroms for
# each sample. chr20-only configs (one entry in chroms) keep working. For
# whole-genome call-only (no labels), see workflows/Snakefile.

configfile: "workflows/config/label.yaml"

import math
import os
import shlex

# Global container: every rule runs inside the SIF when --sdm apptainer is set.
# Empty config and unset ARCSV2_SIF leave CI/laptop runs on the host (Snakemake
# ignores container: unless that deployment method is passed). Paths in the
# config must be the same on the host and in the container (Grace binds /scratch).
_SIF = str(config.get("container") or os.environ.get("ARCSV2_SIF") or "").strip()
if _SIF:
    container: _SIF

# Rules declare conda: for host runs. Inside the SIF the tools are already on
# PATH and the image has micromamba, not conda; Snakemake 8 still runs
# `conda info` in the container when a rule has conda:, so drop it there.
_CONDA_ENV = None if _SIF else "../environment.yml"

DATA = config.get("data_root", ".").rstrip("/") or "."
RESULTS = config.get("results_dir", "results").rstrip("/") or "results"
LABELS = config.get("labels_dir", "labels").rstrip("/") or "labels"

# M4: length-scaled mem_mb = ceil(chrom_Mb * MB_PER_MB * 1.2). See Snakefile.
REF_FAI = str(config.get("reference", "")) + ".fai"
CHROM_LEN_BP: dict[str, int] = {}
if os.path.isfile(REF_FAI):
    with open(REF_FAI, encoding="utf-8") as _fh:
        for _line in _fh:
            _name, _length, *_rest = _line.split("\t")
            CHROM_LEN_BP[_name] = int(_length)
MEM_MB_PER_MB = float(config.get("mem_mb_per_mb", 30.68))
MEM_MB_FLOOR = int(config.get("mem_mb_floor", 4000))

# Call walltime is minutes: max(floor, ceil(chrom_Mb * min_per_mb)).
# chr1 (248.96 Mb) at 0.5 min/Mb is ceil(124.5) = 125 min, not the old
# max(180, Mb*12) ≈ 50 h. Retries multiply runtime and mem_mb by attempt
# (attempt 1 is the base). The SLURM plugin wants runtime in minutes.
#
# Grace walltime caps (minutes). short is 2 h, so a 120 min merge cannot
# retry: attempt 2 would ask for 240 min and sbatch rejects it. _cap keeps
# every request inside the partition the slurm-grace profile assigns
# (merge and dataset on short; the other rules on medium).
PARTITION_LIMIT_MIN = {"short": 120, "medium": 1440, "long": 10080}
FETCH_RUNTIME_MIN = 90
MERGE_RUNTIME_MIN = 30
DATASET_RUNTIME_MIN = 30
COVERAGE_RUNTIME_MIN = 30


def _cap(minutes, partition) -> int:
    """Minutes that fit the named Grace partition (short 2 h, medium 1 day, long 7 days)."""
    limit = PARTITION_LIMIT_MIN[str(partition)]
    return min(int(minutes), int(limit))


def call_runtime_min(wildcards) -> int:
    """Base call runtime in minutes. Unknown contigs use the floor."""
    floor = int(config.get("call_runtime_floor_min", 60))
    length = CHROM_LEN_BP.get(wildcards.chrom)
    if length is None:
        return floor
    per_mb = float(config.get("call_min_per_mb", 0.5))
    mb = length / 1_000_000.0
    return max(floor, int(math.ceil(mb * per_mb)))


def index_runtime_min(_wildcards) -> int:
    return int(config.get("index_runtime_min", 360))


def label_runtime_min(_wildcards) -> int:
    return int(config.get("label_runtime_min", 480))


def call_mem_mb(wildcards) -> int:
    if str(config.get("mem_mb_mode", "scaled")).lower() == "flat":
        return int(config.get("mem_mb", 8000))
    length = CHROM_LEN_BP.get(wildcards.chrom)
    if length is None:
        return int(config.get("mem_mb", 8000))
    return max(
        MEM_MB_FLOOR,
        int(math.ceil((length / 1_000_000.0) * MEM_MB_PER_MB * 1.2)),
    )


def data(path: str) -> str:
    return path if DATA in (".", "") else f"{DATA}/{path}"


def results(path: str) -> str:
    return data(f"{RESULTS}/{path}")


def labels(path: str) -> str:
    return data(f"{LABELS}/{path}")


def hap_fasta(sample, which):
    """Prefer uncompressed FASTA; accept .fa.gz from scripts/fetch_data.py."""
    for ext in (".fa", ".fa.gz"):
        path = data(f"haps/{sample}.{which}{ext}")
        if os.path.exists(path):
            return path
    return data(f"haps/{sample}.{which}.fa")


def coverage_flags() -> str:
    """Reference and gaps for ``arcsv2 ml-coverage``. Empty paths are omitted."""
    parts: list[str] = []
    reference = str(config.get("reference") or "").strip()
    if reference:
        parts.extend(["--reference", shlex.quote(reference)])
    gaps = str(config.get("gaps") or "").strip()
    if gaps:
        parts.extend(["--gaps", shlex.quote(gaps)])
    return " ".join(parts)


def merge_filter_bed_args():
    """Optional paper merge-filter BED flags from config (empty = skip).

    Annotation keys empty => min-size only (matches prior v2fix* tables).
    When any of simple_repeats/segdups/blacklist is set, also pass ``gaps``
    so the full Methods p. 27 filter is applied together.
    """
    parts = []
    annotation = (
        ("--simple-repeats", "simple_repeats"),
        ("--segdups", "segdups"),
        ("--blacklist", "blacklist"),
    )
    any_annotation = False
    for flag, key in annotation:
        path = str(config.get(key) or "").strip()
        if path:
            any_annotation = True
            parts.extend([flag, shlex.quote(path)])
    if any_annotation:
        gaps = str(config.get("gaps") or "").strip()
        if gaps:
            parts.extend(["--gaps", shlex.quote(gaps)])
    return " ".join(parts)


rule all:
    input:
        labels("dataset.tsv"),
        labels("dataset.summary.json"),

rule fetch_region:
    output:
        # temp() so the per-chromosome BAM and index are removed once call and
        # coverage have both finished. coverage lists them as inputs, so they
        # stay until the depth table is written.
        bam=temp(data("work/{sample}/{chrom}.bam")),
        bai=temp(data("work/{sample}/{chrom}.bam.bai")),
    resources:
        mem_mb=lambda wc, attempt: 4000 * attempt,
        runtime=lambda wc, attempt: _cap(FETCH_RUNTIME_MIN * attempt, "medium"),
    params:
        reference=config["reference"],
        cram=lambda wildcards: config["cram"][wildcards.sample],
    conda:
        _CONDA_ENV
    shell:
        "samtools view -u -T {params.reference:q} {params.cram:q} {wildcards.chrom:q} "
        "| samtools sort -o {output.bam:q}\n"
        "samtools index {output.bam:q} {output.bai:q}"

rule call:
    input:
        bam=data("work/{sample}/{chrom}.bam"),
        bai=data("work/{sample}/{chrom}.bam.bai"),
    output:
        results("{sample}/{chrom}/arcsv_out.tab"),
    threads: int(config.get("threads", 4))
    resources:
        mem_mb=lambda wc, attempt: call_mem_mb(wc) * attempt,
        runtime=lambda wc, attempt: _cap(call_runtime_min(wc) * attempt, "medium"),
    params:
        reference=config["reference"],
        gaps=config["gaps"],
        engine=config.get("engine", "rust"),
        out=lambda wildcards: results(f"{wildcards.sample}/{wildcards.chrom}"),
    conda:
        _CONDA_ENV
    shell:
        "arcsv2 call -i {input.bam:q} -R {params.reference:q} -G {params.gaps:q} "
        "-o {params.out:q} --region {wildcards.chrom:q} --threads {threads} "
        "--engine {params.engine}"

rule coverage:
    """Depth columns for this chromosome. Does not rewrite arcsv_out.tab.

    PYTHONPATH=python loads this repo's arcsv2 so the frozen 1.6.1 image,
    which has no ml-coverage command, still runs the working tree. The call
    rule does not set PYTHONPATH, so the caller stays the image build.
    """
    input:
        calls=results("{sample}/{chrom}/arcsv_out.tab"),
        bam=data("work/{sample}/{chrom}.bam"),
        bai=data("work/{sample}/{chrom}.bam.bai"),
    output:
        results("{sample}/{chrom}/coverage.tsv"),
    resources:
        mem_mb=lambda wc, attempt: 4000 * attempt,
        runtime=lambda wc, attempt: _cap(COVERAGE_RUNTIME_MIN * attempt, "short"),
    params:
        extra=coverage_flags(),
    conda:
        _CONDA_ENV
    shell:
        "PYTHONPATH=python arcsv2 ml-coverage --calls {input.calls:q} "
        "--reads {input.bam:q} {params.extra} -o {output:q}"


rule merge:
    input:
        tabs=lambda wildcards: expand(
            results("{sample}/{chrom}/arcsv_out.tab"),
            sample=wildcards.sample,
            chrom=config["chroms"],
        ),
        # Coverage must be in this DAG so the temp BAM is not deleted first.
        # merge-filter still sees only the call tables, so merged.tab is unchanged.
        cov=lambda wildcards: expand(
            results("{sample}/{chrom}/coverage.tsv"),
            sample=wildcards.sample,
            chrom=config["chroms"],
        ),
    output:
        results("{sample}/merged.tab"),
    resources:
        mem_mb=lambda wc, attempt: 8000 * attempt,
        runtime=lambda wc, attempt: _cap(MERGE_RUNTIME_MIN * attempt, "short"),
    params:
        beds=merge_filter_bed_args(),
    conda:
        _CONDA_ENV
    shell:
        "arcsv2 merge-filter {input.tabs:q} {params.beds} -o {output:q}"

rule index_haps:
    """Build one haplotype index, then reject a too-small or unloadable .mmi.

    The size check is the HG00140 failure (stale index, mass no_hit labels).
    ``ml-label`` repeats it, with a smoke alignment, before mapping.
    """
    input:
        lambda wildcards: hap_fasta(wildcards.sample, wildcards.which),
    output:
        # Removed once label, the only consumer, has finished.
        temp(data("haps/{sample}.{which}.mmi")),
    resources:
        mem_mb=lambda wc, attempt: 16000 * attempt,
        runtime=lambda wc, attempt: _cap(index_runtime_min(wc) * attempt, "medium"),
    params:
        reference=config["reference"],
    conda:
        _CONDA_ENV
    shell:
        "minimap2 -d {output:q} {input:q} && "
        "arcsv2 mmi-validate {output:q} --fasta {input:q} --reference {params.reference:q}"

def label_mem_mb(_wildcards) -> int:
    return int(config.get("label_mem_mb", 32000))


def label_threads(_wildcards) -> int:
    return int(config.get("label_threads", config.get("threads", 4)))


rule label:
    input:
        calls=results("{sample}/merged.tab"),
        hap1=lambda wildcards: hap_fasta(wildcards.sample, "hap1"),
        hap2=lambda wildcards: hap_fasta(wildcards.sample, "hap2"),
        idx1=data("haps/{sample}.hap1.mmi"),
        idx2=data("haps/{sample}.hap2.mmi"),
    output:
        labels("{sample}.labeled.tsv"),
    threads: label_threads
    resources:
        mem_mb=lambda wc, attempt: label_mem_mb(wc) * attempt,
        runtime=lambda wc, attempt: _cap(label_runtime_min(wc) * attempt, "medium"),
    params:
        reference=config["reference"],
        unreliable=config.get("unreliable", ""),
        sequential="--sequential" if config.get("sequential") else "",
    conda:
        _CONDA_ENV
    shell:
        "arcsv2 ml-label --calls {input.calls:q} --reference {params.reference:q} "
        "--hap1 {input.hap1:q} --hap2 {input.hap2:q} "
        "--hap1-index {input.idx1:q} --hap2-index {input.idx2:q} "
        "--unreliable {params.unreliable:q} --threads {threads} {params.sequential} "
        "--sample {wildcards.sample:q} -o {output:q}"

rule dataset:
    input:
        expand(labels("{sample}.labeled.tsv"), sample=config["samples"]),
    output:
        tsv=labels("dataset.tsv"),
        summary=labels("dataset.summary.json"),
    resources:
        mem_mb=lambda wc, attempt: 4000 * attempt,
        runtime=lambda wc, attempt: _cap(DATASET_RUNTIME_MIN * attempt, "short"),
    conda:
        _CONDA_ENV
    shell:
        "arcsv2 ml-dataset {input:q} -o {output.tsv:q}"
