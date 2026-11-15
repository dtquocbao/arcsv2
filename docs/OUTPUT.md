# Call outputs: `.tab` and VCF

`arcsv2 call` and `arcsv2 run` write both by default (U1):

| File | Role |
|:--|:--|
| `arcsv_out.tab` | Canonical ARC-SV call table (golden / parity artifact) |
| `arcsv_out.vcf.gz` + `.tbi` | VCF 4.2 view of the **same** calls |

The `.tab` is byte-stable for golden checks. The VCF is derived from it and may
gain INFO/FORMAT fields over time without changing the table.

## `.tab` columns

See `OUTPUT_COLUMNS` in `python/arcsv2/config.py`. Key fields: `chrom`,
`minbp`/`maxbp` (0-based half-open span), `id`, `svtype` / `complextype`,
`rearrangement`, `gt`, `af`, `len_affected`, `inslen`, `sr_support`,
`pe_support`, `score_vs_ref`, `score_vs_next`, `bp_uncertainty`,
`bp_point_source`, `blocks` / `blocks_point`.

## VCF ↔ tab mapping

| VCF | Source (`.tab`) | Notes |
|:--|:--|:--|
| CHROM | `chrom` | |
| POS | `minbp + 1` | VCF 1-based |
| ID | `id` | |
| REF / ALT | `N` / `<SVTYPE>` | Symbolic alts |
| QUAL | `score_vs_ref` | |
| FILTER | `PASS` or `MLFAIL` | ML filter when present |
| INFO/END | `maxbp` | |
| INFO/SVTYPE | `svtype` (first token) or `CPX` | `DEL\|DUP\|INV\|INS\|CPX` |
| INFO/SVLEN | from `len_affected` / `inslen` | Signed negative for DEL |
| INFO/CIPOS, CIEND | `bp_uncertainty` | Symmetric CI from recorded widths |
| INFO/IMPRECISE | `bp_point_source` all `mid` | DRP/midpoint breakpoints |
| INFO/CPX_TYPE | `complextype` | Complex only |
| INFO/REARR | `rearrangement` | |
| INFO/BLOCKS | `blocks` (else `blocks_point`) | |
| INFO/AFFECTED_LEN | `len_affected` | |
| FORMAT/GT | `gt` | HOM→`1/1`, else `0/1` |
| FORMAT/SR | sum of `sr_support` | |
| FORMAT/PE | sum of `pe_support` | |
| FORMAT/SVSCORE | `score_vs_ref` | |
| FORMAT/SVSCORE_NEXT | `score_vs_next` | |

Complex events also emit **BND mate pairs** per novel adjacency (`sv_bp` +
rearrangement orients), linked with `EVENT` / `MATEID`, so breakend-only tools
can consume the call.

## Size scope

Merge-filter and paper evaluation use **≥ 50 bp** affected length (`min_size`).
Insertion lengths written as `>N` are treated as lower bounds and are not
dropped for being short.

## Validation

```bash
# inside arcsv2:dev
python -c "import pysam; print(list(pysam.VariantFile('arcsv_out.vcf.gz')))"
tabix -p vcf arcsv_out.vcf.gz   # if re-indexing
```

Prefer `pysam.VariantFile` in CI (image may not ship `bcftools`). When
`bcftools` is available: `bcftools view -H arcsv_out.vcf.gz | head`.
