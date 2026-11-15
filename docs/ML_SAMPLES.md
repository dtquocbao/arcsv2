# ML-0 sample set (S2)

Pinned list: `workflows/config/samples.ml.tsv` (URLs + checksums).
Caller freeze: `caller-v1.6.1` / image `arcsv2:ml0`.

## Selection rule

- Universe: HPRC Release 2 assemblies intersect NYGC 1KGP high-coverage CRAMs
  (`scripts/fetch_data.py --auto`).
- Size: **40** samples.
- Balance: **8 per superpopulation** (AFR, AMR, EAS, EUR, SAS).
- Population diversity: **at most 3 per population**; cover as many
  populations per superpop as the pool allows.
- Sex: aim **4 female / 4 male** per superpopulation when possible.
- Always include **HG01884** (ACB/AFR) and **NA19700** (ASW/AFR).
- Exclude **HG002** (not in 1KGP; out of `--auto`).
- Relatedness: avoid first-degree pairs (parent–child, full siblings)
  using 1KGP ped files; if a forced pair remains, share `family_id`
  so train/test can keep families together.
- Sex from 1000 Genomes `20130606_sample_info.txt` Gender column.
- Hap checksums pinned from HPRC `.md5` sidecars; CRAM MD5 from the 1KGP index.
- **No download until user confirms this list.**

## Full pool (HPRC R2 ∩ NYGC 1KGP HC)

Intersection size: **223** samples.

### By superpopulation

| superpopulation | n | female | male | other/NA |
|:--|--:|--:|--:|--:|
| AFR | 64 | 41 | 23 | 0 |
| AMR | 43 | 22 | 21 | 0 |
| EAS | 50 | 15 | 35 | 0 |
| EUR | 30 | 17 | 13 | 0 |
| SAS | 36 | 15 | 21 | 0 |
| **total** | **223** | **110** | **113** | **0** |

Note: EUR pool has FIN/GBR/IBS/TSI only (no CEU in this intersection). AFR
covers all 7 1KGP AFR populations present in the pool.

### By population

| population | superpopulation | n | female | male | other/NA |
|:--|:--|--:|--:|--:|--:|
| ACB | AFR | 9 | 6 | 3 | 0 |
| ASW | AFR | 6 | 4 | 2 | 0 |
| ESN | AFR | 10 | 6 | 4 | 0 |
| GWD | AFR | 17 | 11 | 6 | 0 |
| LWK | AFR | 6 | 4 | 2 | 0 |
| MSL | AFR | 8 | 4 | 4 | 0 |
| YRI | AFR | 8 | 6 | 2 | 0 |
| CLM | AMR | 10 | 4 | 6 | 0 |
| MXL | AMR | 2 | 1 | 1 | 0 |
| PEL | AMR | 15 | 10 | 5 | 0 |
| PUR | AMR | 16 | 7 | 9 | 0 |
| CDX | EAS | 5 | 3 | 2 | 0 |
| CHB | EAS | 5 | 2 | 3 | 0 |
| CHS | EAS | 11 | 6 | 5 | 0 |
| JPT | EAS | 16 | 2 | 14 | 0 |
| KHV | EAS | 13 | 2 | 11 | 0 |
| FIN | EUR | 9 | 5 | 4 | 0 |
| GBR | EUR | 10 | 8 | 2 | 0 |
| IBS | EUR | 3 | 2 | 1 | 0 |
| TSI | EUR | 8 | 2 | 6 | 0 |
| BEB | SAS | 10 | 6 | 4 | 0 |
| GIH | SAS | 8 | 4 | 4 | 0 |
| ITU | SAS | 3 | 2 | 1 | 0 |
| PJL | SAS | 10 | 3 | 7 | 0 |
| STU | SAS | 5 | 0 | 5 | 0 |

## Selected 40 — counts

| superpopulation | n | female | male | other/NA | populations |
|:--|--:|--:|--:|--:|:--|
| AFR | 8 | 4 | 4 | 0 | ACB, ASW, ESN, GWD, LWK, MSL, YRI |
| AMR | 8 | 4 | 4 | 0 | CLM, MXL, PEL, PUR |
| EAS | 8 | 4 | 4 | 0 | CDX, CHB, CHS, JPT, KHV |
| EUR | 8 | 4 | 4 | 0 | FIN, GBR, IBS, TSI |
| SAS | 8 | 4 | 4 | 0 | BEB, GIH, ITU, PJL, STU |
| **total** | **40** | **20** | **20** | **0** | |

Population counts in the 40 (max 3): ACB=2, ASW=1, BEB=2, CDX=2, CHB=2, CHS=2, CLM=2, ESN=1, FIN=2, GBR=2, GIH=2, GWD=1, IBS=2, ITU=2, JPT=1, KHV=1, LWK=1, MSL=1, MXL=2, PEL=2, PJL=1, PUR=2, STU=1, TSI=2, YRI=1

## Relatedness

- Pedigree source: `20130606_g1k.ped` (FTP
  `technical/working/20130606_sample_info/`; 1997 families / 3502 individuals /
  1410 first-degree pairs). This covers the 2504 unrelated panel plus related
  samples overlapping the NYGC 3202 high-coverage set. A dedicated public
  `3202.ped` was not found at known IGSR FTP paths (404 as of 2026-10-05).
- First-degree pairs within the 40 (parent–child or full siblings): **none**.
- `family_id` is the ped FID (or sample id when absent). No multi-member
  families appear in this 40, so train/test grouping by `family_id` is
  equivalent to grouping by sample for now.

## Samples

| sample | population | superpopulation | sex | family_id |
|:--|:--|:--|:--|:--|
| HG01884 | ACB | AFR | female | BB03 |
| HG02055 | ACB | AFR | male | BB15 |
| NA19700 | ASW | AFR | male | 2367 |
| HG02922 | ESN | AFR | female | NG06 |
| HG02572 | GWD | AFR | male | GB21 |
| NA19036 | LWK | AFR | female | NA19036 |
| HG03098 | MSL | AFR | male | SL16 |
| NA18505 | YRI | AFR | female | Y005 |
| HG01150 | CLM | AMR | female | CLM12 |
| HG01346 | CLM | AMR | female | CLM27 |
| NA19682 | MXL | AMR | male | m010 |
| NA19776 | MXL | AMR | female | m034 |
| HG01928 | PEL | AMR | male | PEL011 |
| HG01940 | PEL | AMR | female | PEL015 |
| HG00642 | PUR | AMR | male | PR04 |
| HG00738 | PUR | AMR | male | PR07 |
| HG02155 | CDX | EAS | female | HG02155 |
| HG02391 | CDX | EAS | male | HG02391 |
| NA18565 | CHB | EAS | female | NA18565 |
| NA18608 | CHB | EAS | male | NA18608 |
| HG00408 | CHS | EAS | female | SH002 |
| HG00609 | CHS | EAS | male | SH062 |
| NA18940 | JPT | EAS | male | NA18940 |
| HG02080 | KHV | EAS | female | VN068 |
| HG00272 | FIN | EUR | female | HG00272 |
| HG00320 | FIN | EUR | female | HG00320 |
| HG00126 | GBR | EUR | male | HG00126 |
| HG00140 | GBR | EUR | male | HG00140 |
| HG01784 | IBS | EUR | female | IBS053 |
| HG01786 | IBS | EUR | female | IBS054 |
| NA20752 | TSI | EUR | male | NA20752 |
| NA20762 | TSI | EUR | male | NA20762 |
| HG03804 | BEB | SAS | female | BD18 |
| HG03816 | BEB | SAS | female | BD22 |
| NA20850 | GIH | SAS | male | NA20850 |
| NA21102 | GIH | SAS | female | NA21102 |
| HG03742 | ITU | SAS | male | IT015 |
| HG03784 | ITU | SAS | female | IT023 |
| HG02602 | PJL | SAS | male | PK09 |
| HG03688 | STU | SAS | male | ST016 |

