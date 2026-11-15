#!/usr/bin/env bash
# S2 dry-run: list HPRC x 1KGP candidates, select 40, pin hap MD5s, write samples.ml.tsv.
# Does NOT download CRAMs or assemblies.
set -euo pipefail
LOG=/work/bench/logs/ml0/s2_sample_select.log
mkdir -p /work/bench/logs/ml0 /tmp/ml0_s2
exec > >(tee -a "$LOG") 2>&1

echo "=== S2 sample dry-run $(date -Is) ==="
python /work/scripts/fetch_data.py --auto --dest /data --n-samples 200 --chroms chr20 --dry-run \
  --samples-out /tmp/ml0_s2/samples.all.tsv \
  --label-yaml /tmp/ml0_s2/label.all.yaml

python - <<'PY'
from __future__ import annotations

import csv
import importlib.util
import json
import sys
import urllib.request
from collections import Counter
from pathlib import Path

spec = importlib.util.spec_from_file_location("fetch_data", "/work/scripts/fetch_data.py")
fd = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules["fetch_data"] = fd
spec.loader.exec_module(fd)

# --- constants ---
USER_AGENT = "arcsv2-ml0-s2"
PEDIGREE_URL = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/"
    "20130606_sample_info/20130606_sample_info.txt"
)
# fallback known sexes for pilot samples if pedigree fetch fails
KNOWN_SEX = {"HG01884": "female", "NA19700": "female"}
POP_TO_SUPER = {
    "ACB": "AFR", "ASW": "AFR", "ESN": "AFR", "GWD": "AFR", "LWK": "AFR", "MSL": "AFR", "YRI": "AFR",
    "CLM": "AMR", "MXL": "AMR", "PEL": "AMR", "PUR": "AMR",
    "CDX": "EAS", "CHB": "EAS", "CHS": "EAS", "JPT": "EAS", "KHV": "EAS",
    "CEU": "EUR", "FIN": "EUR", "GBR": "EUR", "IBS": "EUR", "TSI": "EUR",
    "BEB": "SAS", "GIH": "SAS", "ITU": "SAS", "PJL": "SAS", "STU": "SAS",
}
MUST_INCLUDE = ["HG01884", "NA19700"]
EXCLUDE = {"HG002"}
N_TARGET = 40
SUPERS = ["AFR", "AMR", "EAS", "EUR", "SAS"]

# --- load candidates via same catalogs as fetch --auto ---
dest = Path("/data")
catalog = dest / "catalog"
kg_u, _ = fd._cache_catalog(fd.INDEX_UNRELATED, catalog / "1000G_2504_high_coverage.sequence.index", False)
kg_r, _ = fd._cache_catalog(
    fd.INDEX_RELATED, catalog / "1000G_698_related_high_coverage.sequence.index", False
)
hprc_p, _ = fd._cache_catalog(
    fd.HPRC_CSV, catalog / f"assemblies_release2_v1.0.{fd.HPRC_COMMIT[:12]}.csv", False
)
kg = fd.parse_1kgp_index(kg_u.read_text(encoding="utf-8"), related=False)
kg.update(fd.parse_1kgp_index(kg_r.read_text(encoding="utf-8"), related=True))
hprc = fd.parse_hprc_csv(hprc_p.read_text(encoding="utf-8"))
candidates = fd.intersect_candidates(kg, hprc)
print(f"intersection_n={len(candidates)}")

# --- sex from 1KGP sample info ---
sex_map: dict[str, str] = dict(KNOWN_SEX)
try:
    req = urllib.request.Request(PEDIGREE_URL, headers={"User-Agent": USER_AGENT})
    text = urllib.request.urlopen(req, timeout=120).read().decode("utf-8", "replace")
    # file is tab; Sample / Gender columns
    reader = csv.DictReader(text.splitlines(), delimiter="\t")
    for row in reader:
        sid = (row.get("Sample") or row.get("sample") or "").strip()
        g = (row.get("Gender") or row.get("gender") or row.get("Sex") or "").strip().lower()
        if not sid or not g:
            continue
        if g.startswith("f"):
            sex_map[sid] = "female"
        elif g.startswith("m"):
            sex_map[sid] = "male"
    print(f"pedigree_sex_n={len(sex_map)}")
except Exception as exc:  # noqa: BLE001
    print(f"WARN pedigree fetch failed: {exc}; using KNOWN_SEX only")

# annotate
for row in candidates:
    pop = row["population"]
    row["superpopulation"] = POP_TO_SUPER.get(pop, "UNK")
    row["sex"] = sex_map.get(row["sample"], "")

# pool stats
by_super = Counter(r["superpopulation"] for r in candidates)
by_sex = Counter(r["sex"] or "NA" for r in candidates)
print("pool_by_super", dict(by_super))
print("pool_by_sex", dict(by_sex))
assert "HG002" not in {r["sample"] for r in candidates}

by_name = {r["sample"]: r for r in candidates}
for must in MUST_INCLUDE:
    assert must in by_name, must

# --- select 40: force includes, then balance super x sex ---
chosen: list[dict[str, str]] = []
chosen_ids: set[str] = set()

def add(row: dict[str, str]) -> None:
    if row["sample"] in chosen_ids or row["sample"] in EXCLUDE:
        return
    chosen.append(row)
    chosen_ids.add(row["sample"])

for must in MUST_INCLUDE:
    add(by_name[must])

# target ~8 per superpop
per_super = N_TARGET // len(SUPERS)  # 8
remaining_slots = {s: per_super for s in SUPERS}
for row in chosen:
    s = row["superpopulation"]
    if s in remaining_slots:
        remaining_slots[s] = max(0, remaining_slots[s] - 1)

# prefer filling sex balance within each super
pool = [r for r in candidates if r["sample"] not in chosen_ids and r["sample"] not in EXCLUDE]
# sort: prefer known sex, then sample id for determinism
pool.sort(key=lambda r: (r["superpopulation"], r["sex"] or "z", r["population"], r["sample"]))

for superpop in SUPERS:
    need = remaining_slots[superpop]
    if need <= 0:
        continue
    sub = [r for r in pool if r["superpopulation"] == superpop]
    # aim for half female / half male among remaining
    females = [r for r in sub if r["sex"] == "female"]
    males = [r for r in sub if r["sex"] == "male"]
    unknown = [r for r in sub if r["sex"] not in ("female", "male")]
    # already chosen sex counts in this super
    already = [r for r in chosen if r["superpopulation"] == superpop]
    n_f = sum(1 for r in already if r["sex"] == "female")
    n_m = sum(1 for r in already if r["sex"] == "male")
    picked_here: list[dict[str, str]] = []
    while len(picked_here) < need and (females or males or unknown):
        want_f = n_f <= n_m
        if want_f and females:
            row = females.pop(0)
        elif (not want_f) and males:
            row = males.pop(0)
        elif females:
            row = females.pop(0)
        elif males:
            row = males.pop(0)
        else:
            row = unknown.pop(0)
        picked_here.append(row)
        if row["sex"] == "female":
            n_f += 1
        elif row["sex"] == "male":
            n_m += 1
    for row in picked_here:
        add(row)

# fill any shortfall from remaining pool round-robin by super
if len(chosen) < N_TARGET:
    leftover = [r for r in pool if r["sample"] not in chosen_ids]
    leftover.sort(key=lambda r: (r["superpopulation"], r["sex"] or "z", r["sample"]))
    for row in leftover:
        if len(chosen) >= N_TARGET:
            break
        add(row)

assert len(chosen) == N_TARGET, len(chosen)
chosen.sort(key=lambda r: (r["superpopulation"], r["population"], r["sample"]))

# pin hap MD5 sidecars (tiny downloads only)
print("fetching hap MD5 sidecars for", len(chosen), "samples")
for row in chosen:
    row["hap1_checksum"] = fd.read_md5_sidecar(row["hap1_md5_url"])
    row["hap2_checksum"] = fd.read_md5_sidecar(row["hap2_md5_url"])

# write samples.ml.tsv (pinned; commit this)
out_tsv = Path("/work/workflows/config/samples.ml.tsv")
comment = (
    "# ML-0 S2 pinned samples (HPRC Release 2 x NYGC 1KGP high-coverage).\n"
    f"# catalogs commit={fd.HPRC_COMMIT}; AF.6f caller-v1.6.1; DO NOT download until confirmed.\n"
    "# columns match samples.tsv.example; superpopulation documented in docs/ML_SAMPLES.md"
)
fd.write_samples_tsv(out_tsv, chosen, comment)
print(f"wrote {out_tsv}")

# summary + markdown
counts_super = Counter(r["superpopulation"] for r in chosen)
counts_sex = Counter(r["sex"] or "NA" for r in chosen)
counts_super_sex = Counter((r["superpopulation"], r["sex"] or "NA") for r in chosen)
print("selected_by_super", dict(counts_super))
print("selected_by_sex", dict(counts_sex))
print("selected_super_sex", {f"{a}/{b}": n for (a, b), n in sorted(counts_super_sex.items())})

md = Path("/work/docs/ML_SAMPLES.md")
lines = [
    "# ML-0 sample set (S2)",
    "",
    "Pinned list: `workflows/config/samples.ml.tsv` (URLs + checksums).",
    "Caller freeze: `caller-v1.6.1` / image `arcsv2:ml0`.",
    "",
    "## Selection rule",
    "",
    "- Universe: HPRC Release 2 assemblies intersect NYGC 1KGP high-coverage CRAMs",
    "  (`scripts/fetch_data.py --auto`).",
    "- Size: **40** samples.",
    "- Balance: target 8 per superpopulation (AFR, AMR, EAS, EUR, SAS).",
    "- Sex: both sexes within each superpopulation when the intersection allows.",
    "- Always include **HG01884** (ACB/AFR) and **NA19700** (ASW/AFR).",
    "- Exclude **HG002** (not in 1KGP; out of `--auto`).",
    "- Sex from 1000 Genomes `20130606_sample_info.txt` Gender column.",
    "- Hap checksums pinned from HPRC `.md5` sidecars; CRAM MD5 from the 1KGP index.",
    "- **No download until user confirms this list.**",
    "",
    "## Counts",
    "",
    "| superpopulation | n | female | male | other/NA |",
    "|:--|--:|--:|--:|--:|",
]
for s in SUPERS:
    n = counts_super.get(s, 0)
    nf = counts_super_sex.get((s, "female"), 0)
    nm = counts_super_sex.get((s, "male"), 0)
    no = n - nf - nm
    lines.append(f"| {s} | {n} | {nf} | {nm} | {no} |")
lines += [
    f"| **total** | **{len(chosen)}** | **{counts_sex.get('female', 0)}** | **{counts_sex.get('male', 0)}** | **{sum(1 for r in chosen if r['sex'] not in ('female','male'))}** |",
    "",
    "## Samples",
    "",
    "| sample | population | superpopulation | sex |",
    "|:--|:--|:--|:--|",
]
for r in chosen:
    lines.append(
        f"| {r['sample']} | {r['population']} | {r['superpopulation']} | {r['sex'] or 'NA'} |"
    )
lines.append("")
md.write_text("\n".join(lines) + "\n", encoding="utf-8")
print(f"wrote {md}")

# machine-readable summary for the agent report
summary = {
    "n": len(chosen),
    "by_super": dict(counts_super),
    "by_sex": dict(counts_sex),
    "samples": [
        {
            "sample": r["sample"],
            "population": r["population"],
            "superpopulation": r["superpopulation"],
            "sex": r["sex"],
        }
        for r in chosen
    ],
}
Path("/work/bench/logs/ml0/s2_selection.json").write_text(json.dumps(summary, indent=2) + "\n")
print("=== S2 selection ready; awaiting user confirm before download ===")
PY
