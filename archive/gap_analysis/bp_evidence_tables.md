# BP + evidence tables (HG01884 ABCBC sites)

## Site ~60792115
### Calls
| tool | span | type | rearr | ref | len | bp |
|---|---|---|---|---|--:|---|
| ORIGINAL | 60792108-60792176 | DUP | ABBC | ABC | 68 | `60787005,60792108,60792176,60795177` |
| ARC-SV2 (pre-simplify) | 60792115-60792170 | BND | ABCBCD | ABCD | 68 | `60787006,60792102,60792144,60792170,60795177` |

### ORIGINAL breakpoints (merged+filtered, ±500 bp of center)
| interval | evidence | notes |
|---|---|---|
| [60792102,60792114] | Dup-split, left-clip×2, split×4 | blockinterval [60787005, 60792102) |
| [60792144] | left-clip×2 | blockinterval [60792114, 60792144) |
| [60792170,60792182] | Dup-split, right-clip×2, split×4 | blockinterval [60792144, 60792170) |

### Side-by-side cut inventory
| pos / interval | ORIGINAL | ARC-SV2 current |
|---|---|---|
| 60787005/6 (locus left) | point cut (flank) | point/CI left flank |
| [60792102,60792114] | **split×4 Dup-** + left-clip×2 → CI; median→60792108 | CI [60792102,60792115) clip |
| **60792144** | **point; left-clip×2 only, Spl=0** (kept in BLOCKS) | **extra cut** → block split; path ABCBCD |
| [60792170,60792182] | **split×4 Dup-** + right-clip×2 → CI; median→60792176 | CI [60792170,60792183) clip |
| 60795177 (locus right) | point cut (flank) | point/CI right flank |

**Why ORIGINAL reports DUP ABBC:** BLOCKS still include the 60792144 cut (4 middle blocks), genotype path is `ABCBCD`, but `simplify_blocks_diploid` (`orig:sv_validate.py`) fuses B+C because B_out↔C_in is only used as reference adjacency (degree 1) in both haplotypes → output `ABC`/`ABBC`. ARC-SV2 lacked this step → emitted BND ABCBCD.

## Site ~64177566
### Calls
| tool | span | type | rearr | ref | len | bp |
|---|---|---|---|---|--:|---|
| ORIGINAL | 64173541-64173754 | BND,BND,BND,BND | ABCDCBCDCBCDE | ABCDE | 213 | `64171612,64173541,64173620,64173699,64173754,64176227` |
| ORIGINAL | 64176227-64176509 | DUP | ABBBC | ABC | 282 | `64173832,64176227,64176509,64177536` |
| ORIGINAL | 64177553-64177671 | DUP | ABBC | ABC | 118 | `64176606,64177553,64177671,64178724` |
| ARC-SV2 (pre-simplify) | 64173542-64173832 | BND,BND,BND | ABCDEFGHDBCDEFGHCDEFGHI | ABCDEFGHI | 291 | `64172765,64173541,64173595,64173624,64173699,64173706,64173717,64173735,64173832,64176227` |
| ARC-SV2 (pre-simplify) | 64177566-64177660 | BND | ABCDCDE | ABCDE | 118 | `64176281,64177536,64177542,64177654,64177660,64178724` |

### ORIGINAL breakpoints (merged+filtered, ±500 bp of center)
| interval | evidence | notes |
|---|---|---|
| [64177536] | left-clip×3 | blockinterval [64176606, 64177536) |
| [64177542,64177565] | Dup-split, left-clip×2, split×7 | blockinterval [64177536, 64177542) |
| [64177654] | right-clip×3 | blockinterval [64177565, 64177654) |
| [64177660,64177683] | Dup-split, left-clip×2, right-clip×22, split×7 | blockinterval [64177654, 64177660) |

### Side-by-side cut inventory
| pos / interval | ORIGINAL | ARC-SV2 current |
|---|---|---|
| 64176606 (locus left) | point/CI flank | flank |
| 64177536 | point; left-clip×3, Spl=0 | cut kept |
| [64177542,64177565] | CI (split Dup support on outer) | CI cut |
| 64177654 | point; right-clip×3, Spl=0 | cut kept |
| [64177660,64177683] | CI | CI cut |
| 64178724 (locus right) | flank | flank |

**Why ORIGINAL reports DUP ABBC:** same `simplify_blocks_diploid` fusion of unused internal cuts on path `ABCDCDE` → `ABC`/`ABBC`. ARC-SV2 kept all half-open CI cuts → BND `ABCDCDE`.
