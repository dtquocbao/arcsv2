# DRP classification parity (orig vs ARC-SV2)

Rule table from `CURSOR_TASK_DRP_PARITY.md` tasks 1–3. Original:
`src/arcsv/arcsv/pecluster.py` `process_discordant_pair` / clustering.
ARC-SV2: `engine/caller.py` `_drp_from_pair`, `engine/breakpoints.py` `cluster_drps`.

## Shared filters

| Rule | Original | ARC-SV2 (after fix) |
|:--|:--|:--|
| Secondary / duplicate / supplementary | Drop | Drop (`samio.py`) |
| MAPQ ≥ 20 both mates | Always for DRPs | Always for DRPs (hangings skipped) |
| Hanging / unmapped mate / other chrom | Never enters DRP dicts | **Not classified as DRP** (graph only) |
| Concordant opposite + insert in [min,max] | Return None first | **Early return** before Del/Dup/Ins (`z_ins ≤ insert ≤ z_del`) |
| Insert length | Signed `process_insert_len` (outer + soft clips) | Same (`samio.insert_length`); KDE keeps `ilen ≥ 0` |

## Per-type typing

| Type | Original rule | ARC-SV2 (after fix) |
|:--|:--|:--|
| **Del** | Opposite; `ilen > max`; `fwd_end < rev_start` | Same |
| **Dup** | Opposite; `(fwd_start > rev_start) OR (fwd_end > rev_end)`; **before** Ins | Same order + geometry |
| **Ins** | Opposite; `ilen < min`; **no** non-overlap requirement; after Dup | Same; overlapping mates keep `(fwd_end, rev_start)` even if `a > b` |
| **InvL** | Both-forward completed pairs; 3′ ends | Same; **not** hangings |
| **InvR** | Both-reverse completed pairs; 5′ starts | Same; **not** hangings |

## Clustering (Ins)

| Knob | Original | ARC-SV2 (after fix) |
|:--|:--|:--|
| Pair compatibility | Close + overlap ≤ `max_ins_pair_slop` (20) | Same (`compatible`) |
| Post-clique filter | `is_ins_cluster_compatible` overlap ≤ 25 | `ins_cluster_ok` / `drp_insertion_cluster_slop=25` |
| Min size / max cliques | ≥2; all largest cliques | Same |

## Why counts diverged (pre-fix)

1. **InvL/InvR ~11×:** `_drp_from_pair` treated `hanging or same-orient` as Inv; distant/unmapped-mate hangings flooded Inv. Original never does.
2. **Dup under-count:** only outward `rev_start < fwd_start`; nested FR Dup clause missing; Ins checked before Dup.
3. **Ins clusters ~25× (v2fix):** non-overlap gate cut DRPs, remaining sparse size-2 cliques all passed FDR; no cluster-slop filter.
4. **Ins DRPs ~27× (v2fix2):** after removing the Ins non-overlap gate, `abs(TLEN)` + missing concordant early return flooded Ins with overlapping/read-through pairs. Fixed by signed `process_insert_len` + early return (see `CURSOR_TASK_INS_DUP.md`).

## Window 5.9–6.1 Mb (expected after fix)

Re-dump with `scripts/diagnostics/dump_drp_window.py` after the chr20 re-run: Inv should collapse toward original’s ~0 in that Del-rich window; nested Dups should reappear; Ins cluster count should fall sharply.
