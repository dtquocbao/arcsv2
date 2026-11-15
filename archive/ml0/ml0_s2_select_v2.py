#!/usr/bin/env python3
"""Revise ML-0 S2 sample selection: diversity + relatedness.

Does NOT download CRAMs or haplotype assemblies — only tiny catalogs/ped/md5.
"""
from __future__ import annotations

import csv
import importlib.util
import json
import sys
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

WORK = Path("/work")
DEST = Path("/data")
LOG_DIR = WORK / "bench" / "logs" / "ml0"
USER_AGENT = "arcsv2-ml0-s2-v2"

# Sex / pedigree sources (small text)
SAMPLE_INFO_URL = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/"
    "20130606_sample_info/20130606_sample_info.txt"
)
PED_URLS = [
    # Classic 1KGP pedigree (2504 + related; 3502 individuals). Covers the
    # NYGC 3202 high-coverage panel samples used here. A dedicated 3202 ped
    # was not found at known public FTP paths (404 as of 2026-10-05).
    (
        "20130606_g1k.ped",
        "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/"
        "20130606_sample_info/20130606_g1k.ped",
    ),
]

POP_TO_SUPER = {
    "ACB": "AFR",
    "ASW": "AFR",
    "ESN": "AFR",
    "GWD": "AFR",
    "LWK": "AFR",
    "MSL": "AFR",
    "YRI": "AFR",
    "CLM": "AMR",
    "MXL": "AMR",
    "PEL": "AMR",
    "PUR": "AMR",
    "CDX": "EAS",
    "CHB": "EAS",
    "CHS": "EAS",
    "JPT": "EAS",
    "KHV": "EAS",
    "CEU": "EUR",
    "FIN": "EUR",
    "GBR": "EUR",
    "IBS": "EUR",
    "TSI": "EUR",
    "BEB": "SAS",
    "GIH": "SAS",
    "ITU": "SAS",
    "PJL": "SAS",
    "STU": "SAS",
}
MUST_INCLUDE = ["HG01884", "NA19700"]
EXCLUDE = {"HG002"}
N_TARGET = 40
SUPERS = ["AFR", "AMR", "EAS", "EUR", "SAS"]
MAX_PER_POP = 3
PER_SUPER = 8


def load_fetch_data():
    spec = importlib.util.spec_from_file_location("fetch_data", WORK / "scripts" / "fetch_data.py")
    fd = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules["fetch_data"] = fd
    spec.loader.exec_module(fd)
    return fd


def fetch_text(url: str, dest: Path, timeout: int = 120) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.is_file() and dest.stat().st_size > 0:
        print(f"cached {dest.name} ({dest.stat().st_size} B)")
        return dest
    print(f"fetching {url}")
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    text = urllib.request.urlopen(req, timeout=timeout).read()
    dest.write_bytes(text)
    print(f"wrote {dest} ({len(text)} B)")
    return dest


def parse_sex_map(path: Path) -> dict[str, str]:
    sex_map: dict[str, str] = {}
    text = path.read_text(encoding="utf-8", errors="replace")
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
    return sex_map


def parse_ped(path: Path) -> tuple[dict[str, str], set[frozenset[str]]]:
    """Return family_id per sample and undirected first-degree pairs.

    First-degree = parent-child or full siblings (share both parents when known).
    PED columns: family_id individual_id paternal_id maternal_id sex phenotype
    (optionally with a leading # header).
    """
    family_of: dict[str, str] = {}
    father: dict[str, str] = {}
    mother: dict[str, str] = {}
    children_of: dict[tuple[str, str], list[str]] = defaultdict(list)

    text = path.read_text(encoding="utf-8", errors="replace")
    for line in text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 4:
            continue
        # Some files have FID IID PAT MAT SEX PHENO; skip non-id header-ish rows
        fid, iid, pat, mat = parts[0], parts[1], parts[2], parts[3]
        if iid.lower() in {"individual_id", "iid", "sample"}:
            continue
        family_of[iid] = fid
        if pat not in {"0", "NA", ""}:
            father[iid] = pat
        if mat not in {"0", "NA", ""}:
            mother[iid] = mat
        if pat not in {"0", "NA", ""} and mat not in {"0", "NA", ""}:
            children_of[(pat, mat)].append(iid)

    pairs: set[frozenset[str]] = set()
    # parent-child
    for child, dad in father.items():
        pairs.add(frozenset((child, dad)))
    for child, mom in mother.items():
        pairs.add(frozenset((child, mom)))
    # full siblings
    for sibs in children_of.values():
        for i, a in enumerate(sibs):
            for b in sibs[i + 1 :]:
                pairs.add(frozenset((a, b)))
    return family_of, pairs


def first_degree_in_set(ids: set[str], pairs: set[frozenset[str]]) -> list[tuple[str, str]]:
    hits = []
    for pair in pairs:
        if len(pair) != 2:
            continue
        a, b = tuple(pair)
        if a in ids and b in ids:
            hits.append(tuple(sorted((a, b))))
    return sorted(hits)


def annotate(candidates: list[dict], sex_map: dict[str, str], family_of: dict[str, str]) -> None:
    for row in candidates:
        pop = row["population"]
        row["superpopulation"] = POP_TO_SUPER.get(pop, "UNK")
        row["sex"] = sex_map.get(row["sample"], "")
        row["family_id"] = family_of.get(row["sample"], row["sample"])


def pool_tables(candidates: list[dict]) -> tuple[list[str], list[str]]:
    """Return markdown lines for pool by superpop and by population."""
    by_super: dict[str, Counter] = defaultdict(Counter)
    by_pop: dict[str, Counter] = defaultdict(Counter)
    pop_super: dict[str, str] = {}
    for r in candidates:
        s = r["superpopulation"]
        p = r["population"]
        sex = r["sex"] or "NA"
        by_super[s][sex] += 1
        by_super[s]["n"] += 1
        by_pop[p][sex] += 1
        by_pop[p]["n"] += 1
        pop_super[p] = s

    lines_super = [
        "| superpopulation | n | female | male | other/NA |",
        "|:--|--:|--:|--:|--:|",
    ]
    for s in SUPERS + [x for x in sorted(by_super) if x not in SUPERS]:
        c = by_super[s]
        n = c["n"]
        nf = c.get("female", 0)
        nm = c.get("male", 0)
        no = n - nf - nm
        lines_super.append(f"| {s} | {n} | {nf} | {nm} | {no} |")
    total_n = len(candidates)
    total_f = sum(1 for r in candidates if r["sex"] == "female")
    total_m = sum(1 for r in candidates if r["sex"] == "male")
    lines_super.append(
        f"| **total** | **{total_n}** | **{total_f}** | **{total_m}** | **{total_n - total_f - total_m}** |"
    )

    lines_pop = [
        "| population | superpopulation | n | female | male | other/NA |",
        "|:--|:--|--:|--:|--:|--:|",
    ]
    for p in sorted(by_pop, key=lambda x: (POP_TO_SUPER.get(x, "ZZZ"), x)):
        c = by_pop[p]
        n = c["n"]
        nf = c.get("female", 0)
        nm = c.get("male", 0)
        no = n - nf - nm
        lines_pop.append(f"| {p} | {pop_super[p]} | {n} | {nf} | {nm} | {no} |")
    return lines_super, lines_pop


def conflicts_with(
    row: dict,
    chosen: list[dict],
    pairs: set[frozenset[str]],
    pop_counts: Counter,
) -> bool:
    if pop_counts[row["population"]] >= MAX_PER_POP:
        return True
    sid = row["sample"]
    chosen_ids = {r["sample"] for r in chosen}
    for other in chosen_ids:
        if frozenset((sid, other)) in pairs:
            return True
    return False


def select_forty(
    candidates: list[dict],
    pairs: set[frozenset[str]],
) -> list[dict]:
    by_name = {r["sample"]: r for r in candidates}
    for must in MUST_INCLUDE:
        assert must in by_name, must

    chosen: list[dict] = []
    chosen_ids: set[str] = set()
    pop_counts: Counter = Counter()

    def add(row: dict) -> bool:
        if row["sample"] in chosen_ids or row["sample"] in EXCLUDE:
            return False
        if conflicts_with(row, chosen, pairs, pop_counts):
            return False
        chosen.append(row)
        chosen_ids.add(row["sample"])
        pop_counts[row["population"]] += 1
        return True

    # Force anchors first (even if they are first-degree to each other —
    # relatedness will surface and family_id will be set).
    for must in MUST_INCLUDE:
        row = by_name[must]
        # Bypass relatedness for forced anchors only against each other
        if row["sample"] in chosen_ids or row["sample"] in EXCLUDE:
            continue
        if pop_counts[row["population"]] >= MAX_PER_POP and row["sample"] not in MUST_INCLUDE:
            raise SystemExit(f"cannot add forced {row['sample']}: pop cap")
        # Allow forced include even if related to another forced include
        related_block = False
        for other in chosen_ids:
            if other in MUST_INCLUDE:
                continue
            if frozenset((row["sample"], other)) in pairs:
                related_block = True
                break
        if related_block:
            raise SystemExit(f"forced {row['sample']} conflicts with non-forced chosen")
        if pop_counts[row["population"]] >= MAX_PER_POP:
            # still allow forced
            pass
        chosen.append(row)
        chosen_ids.add(row["sample"])
        pop_counts[row["population"]] += 1

    # Group pool by superpopulation
    pool = [r for r in candidates if r["sample"] not in chosen_ids and r["sample"] not in EXCLUDE]

    for superpop in SUPERS:
        already = [r for r in chosen if r["superpopulation"] == superpop]
        need = PER_SUPER - len(already)
        if need <= 0:
            continue
        sub = [r for r in pool if r["superpopulation"] == superpop]
        # Cover as many populations as possible: round-robin over pops present
        by_pop: dict[str, list[dict]] = defaultdict(list)
        for r in sub:
            by_pop[r["population"]].append(r)
        for pop in by_pop:
            # prefer known sex, then sample id
            by_pop[pop].sort(key=lambda r: (0 if r["sex"] in ("female", "male") else 1, r["sex"], r["sample"]))

        # Prefer untouched populations first so forced anchors do not starve coverage.
        def pop_order() -> list[str]:
            return sorted(by_pop, key=lambda p: (pop_counts[p], p))

        # Sex targets within this super: aim 4F/4M among the final 8
        n_f = sum(1 for r in already if r["sex"] == "female")
        n_m = sum(1 for r in already if r["sex"] == "male")
        target_f = 4
        target_m = 4

        picked = 0
        # Passes: cover zero-count pops, then fill with sex balance under pop cap.
        pass_idx = 0
        while picked < need and by_pop:
            progressed = False
            for pop in pop_order():
                if picked >= need:
                    break
                if pop_counts[pop] >= MAX_PER_POP:
                    continue
                # On early passes, skip pops that already have ≥1 until all
                # available pops in this super have been touched.
                if pass_idx == 0 and pop_counts[pop] > 0:
                    untouched = [
                        p
                        for p in by_pop
                        if pop_counts[p] == 0
                        and any(
                            r["sample"] not in chosen_ids
                            and not conflicts_with(r, chosen, pairs, pop_counts)
                            for r in by_pop[p]
                        )
                    ]
                    if untouched:
                        continue
                cands = [
                    r
                    for r in by_pop[pop]
                    if r["sample"] not in chosen_ids and not conflicts_with(r, chosen, pairs, pop_counts)
                ]
                if not cands:
                    continue
                want_f = n_f < target_f and (n_f <= n_m or n_m >= target_m)
                want_m = n_m < target_m and (n_m < n_f or n_f >= target_f)
                pick = None
                if want_f:
                    pick = next((r for r in cands if r["sex"] == "female"), None)
                if pick is None and want_m:
                    pick = next((r for r in cands if r["sex"] == "male"), None)
                if pick is None:
                    pick = cands[0]
                if add(pick):
                    if pick["sex"] == "female":
                        n_f += 1
                    elif pick["sex"] == "male":
                        n_m += 1
                    picked += 1
                    progressed = True
                    # one add per pop per pass helps spread coverage
                    if pass_idx == 0:
                        continue
            if not progressed:
                # relax sex preference
                for pop in pop_order():
                    if picked >= need:
                        break
                    if pop_counts[pop] >= MAX_PER_POP:
                        continue
                    cands = [
                        r
                        for r in by_pop[pop]
                        if r["sample"] not in chosen_ids
                        and not conflicts_with(r, chosen, pairs, pop_counts)
                    ]
                    if not cands:
                        continue
                    if add(cands[0]):
                        row = chosen[-1]
                        if row["sex"] == "female":
                            n_f += 1
                        elif row["sex"] == "male":
                            n_m += 1
                        picked += 1
                        progressed = True
            if not progressed:
                break
            pass_idx += 1
            if pass_idx > 20:
                break

        if picked < need:
            print(f"WARN {superpop}: only filled {len(already)+picked}/{PER_SUPER}")

    if len(chosen) < N_TARGET:
        # fill shortfall from any remaining non-conflicting
        leftover = [r for r in pool if r["sample"] not in chosen_ids]
        leftover.sort(key=lambda r: (r["superpopulation"], r["population"], r["sample"]))
        for row in leftover:
            if len(chosen) >= N_TARGET:
                break
            # prefer supers still under PER_SUPER
            super_counts = Counter(r["superpopulation"] for r in chosen)
            if super_counts[row["superpopulation"]] >= PER_SUPER:
                continue
            add(row)

    if len(chosen) != N_TARGET:
        raise SystemExit(f"selected {len(chosen)} != {N_TARGET}")
    chosen.sort(key=lambda r: (r["superpopulation"], r["population"], r["sample"]))
    return chosen


def write_samples_ml_tsv(path: Path, rows: list[dict], fd, comment: str) -> None:
    cols = [
        "sample",
        "population",
        "superpopulation",
        "sex",
        "family_id",
        "cram_url",
        "crai_url",
        "hap1_url",
        "hap2_url",
        "cram_checksum",
        "crai_checksum",
        "hap1_checksum",
        "hap2_checksum",
    ]
    lines = [comment.rstrip(), "\t".join(cols)]
    for row in rows:
        lines.append("\t".join(row.get(col, "") for col in cols))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    # also keep compatibility note: fetch_data REQUIRED ignores extra cols


def main() -> int:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    catalog_cache = DEST / "catalog"
    ped_cache = catalog_cache / "pedigree"

    fd = load_fetch_data()
    kg_u, _ = fd._cache_catalog(
        fd.INDEX_UNRELATED, catalog_cache / "1000G_2504_high_coverage.sequence.index", False
    )
    kg_r, _ = fd._cache_catalog(
        fd.INDEX_RELATED,
        catalog_cache / "1000G_698_related_high_coverage.sequence.index",
        False,
    )
    hprc_p, _ = fd._cache_catalog(
        fd.HPRC_CSV, catalog_cache / f"assemblies_release2_v1.0.{fd.HPRC_COMMIT[:12]}.csv", False
    )
    kg = fd.parse_1kgp_index(kg_u.read_text(encoding="utf-8"), related=False)
    kg.update(fd.parse_1kgp_index(kg_r.read_text(encoding="utf-8"), related=True))
    hprc = fd.parse_hprc_csv(hprc_p.read_text(encoding="utf-8"))
    candidates = fd.intersect_candidates(kg, hprc)
    print(f"intersection_n={len(candidates)}")
    assert "HG002" not in {r["sample"] for r in candidates}

    sample_info = fetch_text(SAMPLE_INFO_URL, ped_cache / "20130606_sample_info.txt")
    sex_map = parse_sex_map(sample_info)
    print(f"sex_map_n={len(sex_map)}")

    family_of: dict[str, str] = {}
    all_pairs: set[frozenset[str]] = set()
    ped_used = []
    for name, url in PED_URLS:
        try:
            path = fetch_text(url, ped_cache / name)
        except Exception as exc:  # noqa: BLE001
            print(f"WARN ped fetch failed {name}: {exc}")
            continue
        fam, pairs = parse_ped(path)
        family_of.update(fam)
        all_pairs |= pairs
        ped_used.append(name)
        print(f"ped {name}: families={len(set(fam.values()))} samples={len(fam)} first_degree_pairs={len(pairs)}")

    if not ped_used:
        raise SystemExit("no pedigree files available")

    annotate(candidates, sex_map, family_of)

    lines_super, lines_pop = pool_tables(candidates)
    print("pool_by_super:")
    print("\n".join(lines_super))

    chosen = select_forty(candidates, all_pairs)
    chosen_ids = {r["sample"] for r in chosen}
    related_hits = first_degree_in_set(chosen_ids, all_pairs)
    print(f"first_degree_pairs_in_40={related_hits}")

    # If any first-degree remain (forced anchors), ensure shared family_id
    for a, b in related_hits:
        fa = family_of.get(a, a)
        fb = family_of.get(b, b)
        fid = fa if fa == fb else f"{fa}|{fb}"
        for r in chosen:
            if r["sample"] in (a, b):
                r["family_id"] = fid

    # Pin hap MD5 sidecars only (tiny)
    print(f"fetching hap MD5 sidecars for {len(chosen)} samples")
    for row in chosen:
        row["hap1_checksum"] = fd.read_md5_sidecar(row["hap1_md5_url"])
        row["hap2_checksum"] = fd.read_md5_sidecar(row["hap2_md5_url"])
        row.setdefault("crai_checksum", "none")

    out_tsv = WORK / "workflows" / "config" / "samples.ml.tsv"
    comment = (
        "# ML-0 S2 pinned samples (HPRC Release 2 x NYGC 1KGP high-coverage).\n"
        f"# catalogs commit={fd.HPRC_COMMIT}; AF.6f caller-v1.6.1; DO NOT download until confirmed.\n"
        "# Extra columns: superpopulation, family_id (for grouped train/test if related)."
    )
    write_samples_ml_tsv(out_tsv, chosen, fd, comment)
    print(f"wrote {out_tsv}")

    counts_super = Counter(r["superpopulation"] for r in chosen)
    counts_sex = Counter(r["sex"] or "NA" for r in chosen)
    counts_super_sex = Counter((r["superpopulation"], r["sex"] or "NA") for r in chosen)
    counts_pop = Counter(r["population"] for r in chosen)

    # family_id column: only keep non-trivial when related or multi-member families in set
    fam_members = defaultdict(list)
    for r in chosen:
        fam_members[r["family_id"]].append(r["sample"])
    multi_fams = {fid: mems for fid, mems in fam_members.items() if len(mems) > 1}

    md = WORK / "docs" / "ML_SAMPLES.md"
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
        "- Balance: **8 per superpopulation** (AFR, AMR, EAS, EUR, SAS).",
        "- Population diversity: **at most 3 per population**; cover as many",
        "  populations per superpop as the pool allows.",
        "- Sex: aim **4 female / 4 male** per superpopulation when possible.",
        "- Always include **HG01884** (ACB/AFR) and **NA19700** (ASW/AFR).",
        "- Exclude **HG002** (not in 1KGP; out of `--auto`).",
        "- Relatedness: avoid first-degree pairs (parent–child, full siblings)",
        "  using 1KGP ped files; if a forced pair remains, share `family_id`",
        "  so train/test can keep families together.",
        "- Sex from 1000 Genomes `20130606_sample_info.txt` Gender column.",
        "- Hap checksums pinned from HPRC `.md5` sidecars; CRAM MD5 from the 1KGP index.",
        "- **No download until user confirms this list.**",
        "",
        "## Full pool (HPRC R2 ∩ NYGC 1KGP HC)",
        "",
        f"Intersection size: **{len(candidates)}** samples.",
        "",
        "### By superpopulation",
        "",
        *lines_super,
        "",
        "### By population",
        "",
        *lines_pop,
        "",
        "## Selected 40 — counts",
        "",
        "| superpopulation | n | female | male | other/NA | populations |",
        "|:--|--:|--:|--:|--:|:--|",
    ]
    for s in SUPERS:
        n = counts_super.get(s, 0)
        nf = counts_super_sex.get((s, "female"), 0)
        nm = counts_super_sex.get((s, "male"), 0)
        no = n - nf - nm
        pops = sorted({r["population"] for r in chosen if r["superpopulation"] == s})
        lines.append(f"| {s} | {n} | {nf} | {nm} | {no} | {', '.join(pops)} |")
    lines += [
        f"| **total** | **{len(chosen)}** | **{counts_sex.get('female', 0)}** | **{counts_sex.get('male', 0)}** | **{sum(1 for r in chosen if r['sex'] not in ('female','male'))}** | |",
        "",
        "Population counts in the 40 (max 3): "
        + ", ".join(f"{p}={n}" for p, n in sorted(counts_pop.items())),
        "",
        "## Relatedness",
        "",
        f"- Pedigree sources used: {', '.join(ped_used) if ped_used else 'none'}.",
        f"- First-degree pairs within the 40: "
        + (
            "**none**"
            if not related_hits
            else ", ".join(f"{a}–{b}" for a, b in related_hits)
        )
        + ".",
    ]
    if multi_fams:
        lines.append(
            "- Shared `family_id` for multi-member families in the set: "
            + "; ".join(f"{fid} ({', '.join(mems)})" for fid, mems in sorted(multi_fams.items()))
            + "."
        )
    else:
        lines.append(
            "- Each sample has its own `family_id` (equal to sample id or ped FID);"
            " no multi-member families in the 40."
        )
    lines += [
        "",
        "## Samples",
        "",
        "| sample | population | superpopulation | sex | family_id |",
        "|:--|:--|:--|:--|:--|",
    ]
    for r in chosen:
        lines.append(
            f"| {r['sample']} | {r['population']} | {r['superpopulation']} | "
            f"{r['sex'] or 'NA'} | {r['family_id']} |"
        )
    lines.append("")
    md.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {md}")

    summary = {
        "pool_n": len(candidates),
        "selected_n": len(chosen),
        "by_super": dict(counts_super),
        "by_sex": dict(counts_sex),
        "by_pop": dict(counts_pop),
        "first_degree_pairs": related_hits,
        "multi_family_ids": multi_fams,
        "ped_used": ped_used,
        "samples": [
            {
                "sample": r["sample"],
                "population": r["population"],
                "superpopulation": r["superpopulation"],
                "sex": r["sex"],
                "family_id": r["family_id"],
            }
            for r in chosen
        ],
    }
    (LOG_DIR / "s2_selection_v2.json").write_text(json.dumps(summary, indent=2) + "\n")
    print("=== S2 v2 selection ready; awaiting user confirm before download ===")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
