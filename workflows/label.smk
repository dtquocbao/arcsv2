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


def _runtime_hms(minutes: int) -> str:
    minutes = max(1, int(minutes))
    hours, mins = divmod(minutes, 60)
    return f"{hours}:{mins:02d}:00"


def call_runtime_slurm(wildcards) -> str:
    """SLURM ``--time`` from chromosome length (generous for whole-genome call)."""
    length = CHROM_LEN_BP.get(wildcards.chrom)
    if length is None:
        return _runtime_hms(240)
    mb = length / 1_000_000.0
    return _runtime_hms(max(180, int(mb * 12)))


def fetch_runtime_slurm(wildcards) -> str:
    del wildcards
    return _runtime_hms(90)


def merge_runtime_slurm(_wildcards) -> str:
    return _runtime_hms(120)


def index_runtime_slurm(_wildcards) -> str:
    return _runtime_hms(360)


def label_runtime_slurm(_wildcards) -> str:
    return _runtime_hms(480)


def dataset_runtime_slurm(_wildcards) -> str:
    return _runtime_hms(60)


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
        bam=data("work/{sample}/{chrom}.bam"),
        bai=data("work/{sample}/{chrom}.bam.bai"),
    resources:
        mem_mb=4000,
        runtime=fetch_runtime_slurm,
    params:
        reference=config["reference"],
        cram=lambda wildcards: config["cram"][wildcards.sample],
    conda:
        "../environment.yml"
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
        mem_mb=call_mem_mb,
        runtime=call_runtime_slurm,
    params:
        reference=config["reference"],
        gaps=config["gaps"],
        engine=config.get("engine", "rust"),
        out=lambda wildcards: results(f"{wildcards.sample}/{wildcards.chrom}"),
    conda:
        "../environment.yml"
    shell:
        "arcsv2 call -i {input.bam:q} -R {params.reference:q} -G {params.gaps:q} "
        "-o {params.out:q} --region {wildcards.chrom:q} --threads {threads} "
        "--engine {params.engine}"

rule merge:
    input:
        lambda wildcards: expand(
            results("{sample}/{chrom}/arcsv_out.tab"),
            sample=wildcards.sample,
            chrom=config["chroms"],
        ),
    output:
        results("{sample}/merged.tab"),
    resources:
        mem_mb=8000,
        runtime=merge_runtime_slurm,
    params:
        beds=merge_filter_bed_args(),
    conda:
        "../environment.yml"
    shell:
        "arcsv2 merge-filter {input:q} {params.beds} -o {output:q}"

rule index_haps:
    input:
        lambda wildcards: hap_fasta(wildcards.sample, wildcards.which),
    output:
        data("haps/{sample}.{which}.mmi"),
    resources:
        mem_mb=16000,
        runtime=index_runtime_slurm,
    conda:
        "../environment.yml"
    shell:
        "minimap2 -d {output:q} {input:q}"

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
        mem_mb=label_mem_mb,
        runtime=label_runtime_slurm,
    params:
        reference=config["reference"],
        unreliable=config.get("unreliable", ""),
        sequential="--sequential" if config.get("sequential") else "",
    conda:
        "../environment.yml"
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
        mem_mb=4000,
        runtime=dataset_runtime_slurm,
    conda:
        "../environment.yml"
    shell:
        "arcsv2 ml-dataset {input:q} -o {output.tsv:q}"
