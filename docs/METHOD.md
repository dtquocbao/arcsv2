# Method notes

Positions are 0-based and half-open. Block names follow `orig:helper.py::block_idx_to_name`.

## Library (§4.2)

The insert pmf is a Gaussian kernel on the integer grid. Bandwidth is `max(1, 0.785 * IQR / n^(1/5))`, the same rule as `orig:bamparser_streaming.py::pmf_kernel_smooth`. Samples are histogrammed and convolved, which is the KDE for integer inserts. The pmf is normalized to sum to 1. `F` and `F̃` follow `orig:call_sv.py::create_insert_cs`, including the linear extrapolation past the last bin.

The discordance threshold is the paper constant `C = exp(4.5)`. The scan direction matches `orig:get_lr_cutoff`. For a normal library this sits at about `m ± 3s`.

`π` is `orig:call_sv.py::compute_pi_robust` with `p = 1e-4`, then the median across libraries. Paper mode pools libraries before the fit (§8.5).

Chunk sampling shuffles the grid and draws without replacement. An empty draw raises an error.

## Breakpoints (§4.3)

Soft clips on one side within 5 bp become one cluster. The reported position is the one with the highest total MAPQ. Opposite-side clusters at most 25 bp apart with overlapping aligned bases are one breakpoint interval.

Split classes are `orig:constants.py::SPLIT_TYPES`.

Discordant pairs form a compatibility graph. Each component contributes its largest clique or cliques. The search is Bron–Kerbosch; a component above the size cap is kept only when it is already a clique, and the cap hit is logged. Scores for deletions and insertions are the Gaussian likelihood ratio in the roadmap; inversions and duplications score by size. The null relocates pairs with a seeded `random.Random`. FDR uses `π0 = min(1, M0/M)` and `α = 0.1`.

Breakpoint intervals use `orig:unif_ci.py` at 95%, including that file's Pearson type VI table. Overlapping intervals merge. A merged interval is kept when it has a discordant cluster or at least two split or clipped reads.

## Graph and paths (§4.4–§4.5)

Nodes are `2·block` (in) and `2·block+1` (out). Reference adjacencies are added at `min_edge_support` and are omitted across assembly gaps (`orig:GenomeGraph.add_ref_path_support`). An inter-block edge whose span is above 2 Mb is removed.

A block is spanned when any edge jumps over it. Cut points and minimal subgraphs follow `orig:sv_inference.py::decompose_graph`. The right-hand index is inclusive (`end_out = 2·end+1`). Loci of one or two reference blocks are discarded. Paths are the iterative DFS in `get_paths_iterative`. Neighbor lists come from an adjacency index and are sorted by node id. A novel insertion is tested where a breakpoint has clips on both sides or an insertion-type discordant pair. The sizes are `10^1.5, 10^2, …` up to the search width. The insertion block is appended and linked between the two reference blocks, and path search may step through it. Expanded loci merge when one starts strictly inside the previous locus. A shared cut-point block is not an overlap (`orig` uses `pyinter.open`). More than `max_paths` traversals (paper default 1,000) skips the locus unless `max_mes_extra` allows a higher edge-support retry.

Every read pair is stored on its adjacency. Support is 1 only when the implied insert falls in the middle 95% of `f`, and a same-block tandem duplication also requires the no-duplication insert to fall outside that window. Path search ignores edges below `min_edge_support`. Within-block reads on every block of the locus stay in the likelihood. The running original passes `list(range(start, end+1))` into `get_edges_in_range` (`sv_inference.py` around the call that builds the read list), so flank within-block edges are included. `duplicated_blocks` is defined and never called. Dropping those flank reads changes the duplication copy number.

## Likelihood and genotyping (§4.6–§4.7)

The normalizer is `orig:compute_normalizing_constant`: `double_int` on `F̃`, plus the gap and novel-insertion corrections from `get_gap_overlap_positions` and `get_insertion_overlap_positions`. Contiguous half-open pieces join into one interval, matching the original gap-overlap cases. Each read contributes the sum of `f(raw offset + block distance) · pmapped` over every placement of its two vertices (`orig:get_block_distances_between_nodes`). `EdgeRead.offset` is that raw offset. `n` is every read in the locus, including the flank within-block reads. A hanging read with no insertion overlap contributes `class_probs[1]` when the mate is unmapped and `class_probs[3]` when the mate is on another chromosome or farther than 2 Mb (`orig:max_pair_distance`). The hanging length is the mate query length when that line is in the file. When the path contains a novel insertion, the overlap term from `orig:compute_hanging_edge_likelihood` is added for the first placement of the anchored node.

At allele fraction 1/2 the diploid likelihood is the paper form `Π [2π + (1-π)(p1+p2)] / (G1+G2)^n`. Other fractions use `orig:diploid_likelihood_frac`. The search keeps the top 50 haplotypes plus the reference and maximizes over pairs.

A homozygous k-copy haplotype and a heterozygous (k+1)-copy haplotype can sit within a fraction of a nat of each other once the uncovered breakpoint interval shortens the extra copy. Genotypes within 0.5 nats of the best log-likelihood are one band. The lower-copy, then homozygous, preference applies only when those paths collapse to the same block sequence. A deletion or inversion keeps whichever zygosity has the higher likelihood. `score_vs_next` is the log ratio against the best genotype outside that band. `score_vs_ref` is the log ratio against the homozygous reference.

Short homozygous duplications (300–800 bp) are still sometimes called heterozygous when that likelihood leads by more than 0.5 nats. A 100 bp homozygous insertion can do the same. Those genotypes are left as the likelihood ranks them until they can be compared with the original caller on the same reads.

An insertion junction is reported at the soft-clip base when that base sits inside the breakpoint interval. A tested insertion length is written `>N` when it is the largest size on the grid, longer than the upper insert quantile, or when no discordant pair has ends on both sides of the insertion. A concordant pair that only steps across the reference breakpoint does not contain the novel sequence. Hanging and clipped reads set a lower bound, which is the homozygous insertion longer than the fragment. Features do not treat that value as an exact length, and the short-call filter does not drop it.

When a locus exceeds `max_paths`, edge support is raised by 2 and the search is repeated, up to `max_mes_extra` (`orig:sv_inference.py`). The default is 0, so the paper path still skips the locus.

`bp_uncertainty` is 0 at the two locus ends. Each internal boundary uses `orig:block_gap`: `next.start - this.end`, the uncovered breakpoint interval. `sv_bp_uncertainty` is that gap at each junction endpoint, written `u1,u2` and joined with `;` when a call has several junctions. The original sv string subtracts 2 from each closed interval (`hi - lo - 2`). This column reports the half-open gap width.

## Classification (§4.8)

Complex names are the table in `orig:sv_classify.py::classify_complex`, also matched after `path_reversed`. Affected length is Needleman–Wunsch on `(block, orientation)` tokens with match 1000, mismatch −1, gap −1. Every co-optimal traceback is kept, the strings are swapped, and a block is affected if any of those alignments leaves it unmatched.

Output columns keep the original `arcsv_out.tab` header. `blocks` is appended: each locus block as `start-end`, or `_` for an insertion block, in name order, with breakpoint intervals left uncovered. `blocks_point` is the same list closed at the soft-clip base when the interval has one, otherwise at the interval midpoint. A lone soft clip inside an existing interval can set that point even when discovery needed two reads to open the interval; clip positions are sorted and looked up with `bisect`. `bp_point_source` marks each internal junction as `clip` or `mid`. Breakpoint intervals closer than `min_affected_length` merge before blocks are built. A discordant-only interval within 100 bp of a clip- or split-supported breakpoint is absorbed into that breakpoint when it is not the mate end of the same DRP cluster and contains no soft clip.

## Labels and models (§4.10, §7)

Validation uses `blocks_point` and 150 bp of each outer block next to the event. A shorter outer block is padded from the adjacent reference until that side is 150 bp. When a junction is `mid`, labeling tries positions inside that uncovered interval (edges, midpoint, and a 10 bp grid, at most three such junctions per call) and records `label_bp_shift` when a shift validates the call. A hit must cover the whole query, have MAPQ ≥ 60, BLAST identity above 0.90, and no CIGAR indel of 49 bp or longer, on either strand. Assembly contigs are not oriented to the reference. `ml-label` aligns the predicted haplotype and the reference window to both assemblies with minimap2 through `mappy`, default parameters and no preset. `--sequential` maps base alt/ref against hap1 then hap2 (two loads). Mid-shift
queries run only when the base label is 0 (joint: immediately; sequential: a
second pair of loads). The reference FASTA is opened once and each call fetches its window. The label is 1 when the alt passes on at least one haplotype and the alt and reference do not both pass on the same haplotype, so a heterozygote whose other haplotype carries the reference allele is 1. It is 0 when the alt passes nowhere, and NA when both alleles pass on one haplotype, the call has `_` or `=`, or the window overlaps an unreliable-region BED. `gt_concordant` is recorded separately: homozygous needs the alt on both haplotypes, heterozygous on exactly one. Two whole-genome indexes plus the reference need well over 16 GB; label one sample at a time. `ml-dataset` concatenates labeled tables, drops NA labels, and writes counts of 1, 0, and NA. Feature schema `arcsv2-f1` is the §7.2 reconstruction. A labeled call table fits one logistic regression per class on standardized features. PCA is fit per class, both in each training fold and in the saved model. Folds are grouped by sample id when that column is present, and a chromosome hold-out is stored beside the threshold. A one-class fold does not contribute a score. Thresholds use the remaining out-of-fold probabilities, keep recall ≥ 0.90, and then maximize precision. Without labels, `ml-train` writes the paper presets and `ml-predict` reports `ml_model=untrained` with `ml_prob=NA`. BAM and CRAM are read with `pysam` over the file index; a CRAM needs the reference FASTA. The chromosome edge index is built once, and each worker receives the edges of its locus. Graph construction stays serial until a chromosome-scale timing shows it is worth splitting.
