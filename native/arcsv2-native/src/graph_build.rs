//! R2: block prep + build_graph kernel (bit-identical to Python caller/graph).
//!
//! Pair-end resolution uses rayon indexed par_iter (row order preserved).
//! used_splits claims stay serial in pair order. Support counts and edge
//! insertion order match Python dict / list append order.

use indexmap::IndexMap;
use rayon::prelude::*;

pub const FLAG_REVERSE1: u8 = 1 << 0;
pub const FLAG_REVERSE2: u8 = 1 << 1;
pub const FLAG_HANGING: u8 = 1 << 2;
pub const FLAG_DISTANT: u8 = 1 << 3;

/// EdgeRead flag bits in the columnar export.
pub const ER_HANGING: u8 = 1 << 0;
pub const ER_DISCORDANT: u8 = 1 << 1;
pub const ER_DISTANT: u8 = 1 << 2;
pub const ER_SUPPORTS_ADJ: u8 = 1 << 3;
pub const ER_IN_LIKELIHOOD: u8 = 1 << 4;

pub const KIND_PE: u8 = 0;
pub const KIND_SPLIT: u8 = 1;

#[derive(Clone, Debug)]
pub struct BlockIn {
    pub start: i64,
    pub end: i64,
    pub is_insertion: bool,
    pub insertion_len: i32,
    pub gap_left: bool,
    pub clip_pos: i64, // -1 = None
    pub clip_reads: i32,
}

impl BlockIn {
    fn len(&self) -> i32 {
        if self.is_insertion {
            self.insertion_len
        } else {
            (self.end - self.start).max(0) as i32
        }
    }
}

#[derive(Clone, Debug)]
pub struct SplitIn {
    pub start1: i64,
    pub end1: i64,
    pub reverse1: bool,
    pub start2: i64,
    pub end2: i64,
    pub reverse2: bool,
    pub mapq: u8,
    pub primary_start: i64,
    pub primary_reverse: bool,
    pub qname_hash: i64,
}

#[derive(Clone, Debug)]
pub struct GraphConfig {
    pub min_mapq: u8,
    pub min_edge_support: i32,
    pub max_edge_span: i64,
    pub split_read_leeway: i32,
    pub pool_libraries: bool,
    pub skip_edge_ref_gate: bool,
    pub pmapped: f64,
    pub q_low: i32,
    pub q_high: i32,
}

#[derive(Clone, Debug)]
pub struct PairsIn<'a> {
    pub start1: &'a [i64],
    pub end1: &'a [i64],
    pub start2: &'a [i64],
    pub end2: &'a [i64],
    pub insert: &'a [i32],
    pub mate_len: &'a [i32],
    pub qname_hash: &'a [i64],
    pub mapq1: &'a [u8],
    pub mapq2: &'a [u8],
    pub lib_id: &'a [u16],
    pub flags: &'a [u8],
}

impl PairsIn<'_> {
    pub fn len(&self) -> usize {
        self.start1.len()
    }
}

#[derive(Clone, Debug)]
struct EdgeReadRec {
    offset: i32,
    lib_id: u16,
    pmapped: f64,
    mapq: u8,
    kind: u8,
    er_flags: u8,
    mate_len: i32,
    adj1: Option<Vec<i32>>,
    adj2: Option<Vec<i32>>,
}

#[derive(Clone, Debug)]
struct EdgeRec {
    u: i32,
    v: i32,
    support: i32,
    reference: bool,
    reads: Vec<EdgeReadRec>,
}

#[derive(Clone, Debug)]
pub struct GraphColumns {
    pub edge_u: Vec<i32>,
    pub edge_v: Vec<i32>,
    pub edge_support: Vec<i32>,
    pub edge_reference: Vec<u8>,
    pub edge_read_ptr: Vec<i32>, // len = n_edges + 1 CSR
    pub read_offset: Vec<i32>,
    pub read_lib_id: Vec<u16>,
    pub read_pmapped: Vec<f64>,
    pub read_mapq: Vec<u8>,
    pub read_kind: Vec<u8>,
    pub read_flags: Vec<u8>,
    pub read_mate_len: Vec<i32>,
    pub read_adj1_ptr: Vec<i32>, // per-read; -1 = None; else index into adj_flat
    pub read_adj1_len: Vec<i32>,
    pub read_adj2_ptr: Vec<i32>,
    pub read_adj2_len: Vec<i32>,
    pub adj_flat: Vec<i32>,
    /// Updated clip_pos / clip_reads after block prep (same length as input blocks).
    pub block_clip_pos: Vec<i64>,
    pub block_clip_reads: Vec<i32>,
}

fn in_node(block: i32) -> i32 {
    2 * block
}
fn out_node(block: i32) -> i32 {
    2 * block + 1
}
fn read_vertex(block_idx: i32, is_reverse: bool) -> i32 {
    if is_reverse {
        in_node(block_idx)
    } else {
        out_node(block_idx)
    }
}

fn find_block(blocks: &[BlockIn], pos: i64) -> Option<usize> {
    let mut n = blocks.len();
    while n > 0 && blocks[n - 1].is_insertion {
        n -= 1;
    }
    let mut lo = 0usize;
    let mut hi = n;
    while lo < hi {
        let mid = (lo + hi) / 2;
        if pos < blocks[mid].end {
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }
    if lo < n {
        Some(lo)
    } else {
        None
    }
}

fn in_one_bp_gap(blocks: &[BlockIn], pos: i64) -> Option<(usize, usize)> {
    let idx = find_block(blocks, pos)?;
    if idx == 0 {
        return None;
    }
    let prev = &blocks[idx - 1];
    let nxt = &blocks[idx];
    if prev.is_insertion || nxt.is_insertion {
        return None;
    }
    if nxt.start - prev.end == 1 && prev.end <= pos && pos < nxt.start {
        Some((idx - 1, idx))
    } else {
        None
    }
}

fn overlapping_vertices(blocks: &[BlockIn], start: i64, end: i64, reverse: bool) -> Vec<i32> {
    if end <= start {
        return Vec::new();
    }
    let Some(first) = find_block(blocks, start) else {
        return Vec::new();
    };
    let mut hits: Vec<usize> = Vec::new();
    for i in first..blocks.len() {
        let block = &blocks[i];
        if block.is_insertion {
            continue;
        }
        if block.start >= end {
            break;
        }
        if start < block.end && end > block.start {
            hits.push(i);
        }
    }
    if hits.is_empty() {
        if let Some(gap) = in_one_bp_gap(blocks, start) {
            hits.push(gap.1);
        }
    }
    if reverse {
        hits.into_iter()
            .rev()
            .map(|b| 2 * b as i32)
            .collect()
    } else {
        hits.into_iter()
            .map(|b| 2 * b as i32 + 1)
            .collect()
    }
}

fn read_offset(start: i64, end: i64, reverse: bool, block: &BlockIn) -> i32 {
    if reverse {
        (end - block.start) as i32
    } else {
        (block.start - start) as i32
    }
}

fn node_gap(blocks: &[BlockIn], node: i32) -> i32 {
    let idx = (node / 2) as usize;
    let block = &blocks[idx];
    if block.is_insertion {
        return 0;
    }
    if node % 2 == 1 {
        if idx + 1 < blocks.len() && !blocks[idx + 1].is_insertion {
            return (blocks[idx + 1].start - blocks[idx].end) as i32;
        }
        return 0;
    }
    if idx > 0 && !blocks[idx - 1].is_insertion {
        return (blocks[idx].start - blocks[idx - 1].end) as i32;
    }
    0
}

fn block_distance(path: &[i32], blocks: &[BlockIn], i: usize, j: usize) -> i32 {
    if i / 2 == j / 2 {
        return 0;
    }
    let mut dist: i32 = 0;
    if (path[i] % 2) as usize == i % 2 {
        dist += blocks[(path[i] / 2) as usize].len();
    }
    if (path[j] % 2) as usize != j % 2 {
        dist += blocks[(path[j] / 2) as usize].len();
    }
    let start = if i % 2 == 0 { i + 2 } else { i + 1 };
    let end = if j % 2 == 1 { j - 1 } else { j };
    let mut k = start;
    while k < end {
        dist += blocks[(path[k] / 2) as usize].len();
        k += 2;
    }
    let mut gap: i32 = 0;
    // range(start - 1, end + 1) in Python — careful with underflow
    let gap_lo = start as i64 - 1;
    let gap_hi = end as i64 + 1;
    let mut gk = gap_lo;
    while gk < gap_hi {
        if gk >= 0 && (gk as usize) < path.len() {
            gap += node_gap(blocks, path[gk as usize]);
        }
        gk += 1;
    }
    dist += gap / 2;
    if i % 2 == 0 {
        dist = -dist;
    }
    dist
}

fn block_seq_to_path(seq: &[i32]) -> Vec<i32> {
    if seq.is_empty() {
        return Vec::new();
    }
    let mut path = vec![seq[0]];
    if seq.len() == 1 {
        return path;
    }
    for &node in &seq[1..seq.len() - 1] {
        let flipped = if node % 2 == 0 { node + 1 } else { node - 1 };
        path.push(flipped);
        path.push(node);
    }
    let last = seq[seq.len() - 1];
    let flipped = if last % 2 == 0 { last + 1 } else { last - 1 };
    path.push(flipped);
    path
}

fn adj_of(verts: &[i32]) -> Option<Vec<i32>> {
    if verts.len() <= 1 {
        None
    } else {
        Some(block_seq_to_path(verts))
    }
}

fn end_vertices(
    blocks: &[BlockIn],
    start: i64,
    end: i64,
    reverse: bool,
) -> Option<(Vec<i32>, i32)> {
    let verts = overlapping_vertices(blocks, start, end, reverse);
    if verts.is_empty() {
        let block_idx = find_block(blocks, start)?;
        let v = read_vertex(block_idx as i32, reverse);
        let off = read_offset(start, end, reverse, &blocks[block_idx]);
        return Some((vec![v], off));
    }
    let off = read_offset(start, end, reverse, &blocks[(verts[0] / 2) as usize]);
    Some((verts, off))
}

fn segment_agrees(
    verts: &[i32],
    is_first: bool,
    reverse: bool,
    start: i64,
    end: i64,
    blocks: &[BlockIn],
    leeway: i32,
) -> bool {
    let vertex = if is_first {
        *verts.last().unwrap()
    } else {
        verts[0]
    };
    let block = &blocks[(vertex / 2) as usize];
    if is_first != reverse {
        (block.end - end) as i32 <= leeway
    } else {
        (start - block.start) as i32 <= leeway
    }
}

fn split_path(
    hit: &SplitIn,
    blocks: &[BlockIn],
    leeway: i32,
) -> Option<(Vec<i32>, i32)> {
    let first = overlapping_vertices(blocks, hit.start1, hit.end1, hit.reverse1);
    let second = overlapping_vertices(blocks, hit.start2, hit.end2, hit.reverse2);
    if first.is_empty() || second.is_empty() {
        return None;
    }
    if !segment_agrees(
        &first,
        true,
        hit.reverse1,
        hit.start1,
        hit.end1,
        blocks,
        leeway,
    ) {
        return None;
    }
    if !segment_agrees(
        &second,
        false,
        hit.reverse2,
        hit.start2,
        hit.end2,
        blocks,
        leeway,
    ) {
        return None;
    }
    let offset = read_offset(
        hit.start1,
        hit.end1,
        hit.reverse1,
        &blocks[(first[0] / 2) as usize],
    );
    let mut verts = first;
    verts.extend(second);
    Some((verts, offset))
}

fn extended_badj(verts: &[i32]) -> Vec<i32> {
    let mut badj = block_seq_to_path(verts);
    if badj.len() > 1 {
        let last = *badj.last().unwrap();
        badj.push(last ^ 1);
    }
    badj
}

fn between_junction(verts1: &[i32], verts2: &[i32]) -> Option<(i32, i32)> {
    let badj1 = extended_badj(verts1);
    let badj2 = extended_badj(verts2);
    if badj1.is_empty() || badj2.is_empty() {
        return None;
    }
    Some((*badj1.last().unwrap(), *badj2.last().unwrap()))
}

fn within_read_path_edges(verts: &[i32]) -> Vec<(i32, i32)> {
    if verts.len() <= 1 {
        return Vec::new();
    }
    let path = block_seq_to_path(verts);
    let mut edges = Vec::new();
    let mut i = 0;
    while i + 1 < path.len() {
        edges.push((path[i], path[i + 1]));
        i += 2;
    }
    edges
}

fn pair_path(v1: i32, v2: i32) -> Vec<i32> {
    vec![v1 ^ 1, v1, v2, v2 ^ 1]
}

fn adjacency_insert(
    start1: i64,
    end1: i64,
    rev1: bool,
    start2: i64,
    end2: i64,
    rev2: bool,
    blocks: &[BlockIn],
    b1: i32,
    b2: i32,
) -> i32 {
    let v1 = read_vertex(b1, rev1);
    let v2 = read_vertex(b2, rev2);
    let mut offset = read_offset(start1, end1, rev1, &blocks[b1 as usize]);
    offset += read_offset(start2, end2, rev2, &blocks[b2 as usize]);
    let path = pair_path(v1, v2);
    block_distance(&path, blocks, 1, path.len() - 2) + offset
}

fn is_outward(start1: i64, rev1: bool, start2: i64, rev2: bool) -> bool {
    if rev1 == rev2 {
        return false;
    }
    let rev_start = if rev1 { start1 } else { start2 };
    let fwd_start = if rev1 { start2 } else { start1 };
    rev_start < fwd_start
}

fn dup_without_insert(
    start1: i64,
    end1: i64,
    rev1: bool,
    start2: i64,
    end2: i64,
    rev2: bool,
    blocks: &[BlockIn],
    block: i32,
) -> i32 {
    let v1 = read_vertex(block, rev1);
    let v2 = read_vertex(block, rev2);
    let mut offset = read_offset(start1, end1, rev1, &blocks[block as usize]);
    offset += read_offset(start2, end2, rev2, &blocks[block as usize]);
    let full = pair_path(v1, v2);
    let alt = &full[2..];
    if alt.len() < 2 {
        return offset;
    }
    block_distance(alt, blocks, 1, alt.len() - 2) + offset
}

fn ref_path_insert(
    start1: i64,
    end1: i64,
    rev1: bool,
    start2: i64,
    end2: i64,
    rev2: bool,
    blocks: &[BlockIn],
    b1: i32,
    b2: i32,
) -> i32 {
    let left = b1.min(b2);
    let right = b1.max(b2);
    let v1 = read_vertex(b1, rev1);
    let v2 = read_vertex(b2, rev2);
    let mut offset = read_offset(start1, end1, rev1, &blocks[b1 as usize]);
    offset += read_offset(start2, end2, rev2, &blocks[b2 as usize]);
    if v1 == out_node(left) && v2 == in_node(right) {
        let mut path = vec![v1 ^ 1, v1];
        for b in (left + 1)..right {
            path.push(in_node(b));
            path.push(out_node(b));
        }
        path.push(v2);
        path.push(v2 ^ 1);
        return block_distance(&path, blocks, 1, path.len() - 2) + offset;
    }
    if v1 == out_node(right) && v2 == in_node(left) {
        let mut path = vec![v1 ^ 1, v1];
        let mut b = right - 1;
        while b > left {
            path.push(out_node(b));
            path.push(in_node(b));
            b -= 1;
        }
        path.push(v2);
        path.push(v2 ^ 1);
        return block_distance(&path, blocks, 1, path.len() - 2) + offset;
    }
    adjacency_insert(start1, end1, rev1, start2, end2, rev2, blocks, b1, b2)
}

fn pair_supported(
    start1: i64,
    end1: i64,
    rev1: bool,
    start2: i64,
    end2: i64,
    rev2: bool,
    blocks: &[BlockIn],
    v1: i32,
    v2: i32,
    cfg: &GraphConfig,
) -> bool {
    let bi = v1 / 2;
    let bk = v2 / 2;
    let same_block_dup = bi == bk && (v1 - v2).abs() == 1;
    let implied = adjacency_insert(start1, end1, rev1, start2, end2, rev2, blocks, bi, bk);
    if !(cfg.q_low <= implied && implied <= cfg.q_high) {
        return false;
    }
    if same_block_dup {
        if !is_outward(start1, rev1, start2, rev2)
            || dup_without_insert(start1, end1, rev1, start2, end2, rev2, blocks, bi) >= cfg.q_low
        {
            return false;
        }
        return true;
    }
    if cfg.skip_edge_ref_gate && (bi - bk).abs() > 1 {
        let ref_implied =
            ref_path_insert(start1, end1, rev1, start2, end2, rev2, blocks, bi, bk);
        if cfg.q_low <= ref_implied && ref_implied <= cfg.q_high {
            return false;
        }
    }
    true
}

fn between_junction_supported(
    start1: i64,
    _end1: i64,
    rev1: bool,
    start2: i64,
    _end2: i64,
    rev2: bool,
    verts1: &[i32],
    verts2: &[i32],
    off1: i32,
    off2: i32,
    blocks: &[BlockIn],
    cfg: &GraphConfig,
) -> bool {
    let badj1 = extended_badj(verts1);
    let badj2 = extended_badj(verts2);
    if badj1.is_empty() || badj2.is_empty() {
        return false;
    }
    let mut combined = badj1.clone();
    combined.extend(badj2.iter().rev().copied());
    if combined.len() < 2 {
        return false;
    }
    let mut full = vec![combined[0] ^ 1];
    full.extend_from_slice(&combined);
    full.push(*combined.last().unwrap() ^ 1);
    let dist = block_distance(&full, blocks, 1, full.len() - 2) + off1 + off2;
    if !(cfg.q_low <= dist && dist <= cfg.q_high) {
        return false;
    }
    let j1 = *badj1.last().unwrap();
    let j2 = *badj2.last().unwrap();
    let same_block = j1 / 2 == j2 / 2 && (j1 - j2).abs() == 1;
    if same_block {
        let m = badj1.len() - 1;
        let mut alt_full = Vec::new();
        alt_full.extend_from_slice(&full[..m]);
        alt_full.extend_from_slice(&full[(m + 2)..]);
        if alt_full.len() < 2 {
            return false;
        }
        let alt_dist = block_distance(&alt_full, blocks, 1, alt_full.len() - 2) + off1 + off2;
        if alt_dist >= cfg.q_low {
            return false;
        }
        if verts1.len() == 1 && verts2.len() == 1 && !is_outward(start1, rev1, start2, rev2) {
            return false;
        }
    }
    true
}

fn add_edge(
    edges: &mut IndexMap<(i32, i32), EdgeRec>,
    u: i32,
    v: i32,
    support: i32,
    read: Option<EdgeReadRec>,
    reference: bool,
) {
    let key = (u.min(v), u.max(v));
    let edge = edges.entry(key).or_insert_with(|| EdgeRec {
        u: key.0,
        v: key.1,
        support: 0,
        reference: false,
        reads: Vec::new(),
    });
    if reference {
        edge.reference = true;
        edge.support = edge.support.max(support);
    } else {
        edge.support += support;
    }
    if let Some(r) = read {
        edge.reads.push(r);
    }
}

fn er_flags(
    hanging: bool,
    discordant: bool,
    distant: bool,
    supports_adj: bool,
    in_likelihood: bool,
) -> u8 {
    let mut f = 0u8;
    if hanging {
        f |= ER_HANGING;
    }
    if discordant {
        f |= ER_DISCORDANT;
    }
    if distant {
        f |= ER_DISTANT;
    }
    if supports_adj {
        f |= ER_SUPPORTS_ADJ;
    }
    if in_likelihood {
        f |= ER_IN_LIKELIHOOD;
    }
    f
}

fn add_reference_edges(edges: &mut IndexMap<(i32, i32), EdgeRec>, blocks: &[BlockIn], min_support: i32) {
    for i in 0..blocks.len().saturating_sub(1) {
        if blocks[i + 1].gap_left {
            continue;
        }
        if blocks[i].is_insertion || blocks[i + 1].is_insertion {
            continue;
        }
        add_edge(
            edges,
            out_node(i as i32),
            in_node((i + 1) as i32),
            min_support,
            None,
            true,
        );
    }
}

fn add_hanging(
    edges: &mut IndexMap<(i32, i32), EdgeRec>,
    pairs: &PairsIn,
    i: usize,
    blocks: &[BlockIn],
    cfg: &GraphConfig,
) {
    if pairs.mapq1[i] < cfg.min_mapq {
        return;
    }
    let Some(block_idx) = find_block(blocks, pairs.start1[i]) else {
        return;
    };
    let rev1 = pairs.flags[i] & FLAG_REVERSE1 != 0;
    let distant = pairs.flags[i] & FLAG_DISTANT != 0;
    let anchor = in_node(block_idx as i32);
    let lib_id = if cfg.pool_libraries {
        0
    } else {
        pairs.lib_id[i]
    };
    add_edge(
        edges,
        anchor,
        anchor + 1,
        0,
        Some(EdgeReadRec {
            offset: read_offset(
                pairs.start1[i],
                pairs.end1[i],
                rev1,
                &blocks[block_idx],
            ),
            lib_id,
            pmapped: cfg.pmapped,
            mapq: pairs.mapq1[i],
            kind: KIND_PE,
            er_flags: er_flags(true, false, distant, false, true),
            mate_len: pairs.mate_len[i],
            adj1: None,
            adj2: None,
        }),
        false,
    );
}

fn split_agrees(
    is_first: bool,
    reverse: bool,
    start: i64,
    end: i64,
    block: &BlockIn,
    leeway: i32,
) -> bool {
    if is_first != reverse {
        (block.end - end) as i32 <= leeway
    } else {
        (start - block.start) as i32 <= leeway
    }
}

fn split_support_edges(hit: &SplitIn, blocks: &[BlockIn], cfg: &GraphConfig) -> Vec<(i32, i32)> {
    if let Some((verts, _)) = split_path(hit, blocks, cfg.split_read_leeway) {
        let path = block_seq_to_path(&verts);
        let mut edges = Vec::new();
        let mut i = 0;
        while i + 1 < path.len() {
            edges.push((path[i], path[i + 1]));
            i += 2;
        }
        if !edges.is_empty() {
            return edges;
        }
    }
    let Some(b1) = find_block(blocks, hit.start1) else {
        return Vec::new();
    };
    let Some(b2) = find_block(blocks, hit.start2) else {
        return Vec::new();
    };
    if !split_agrees(
        true,
        hit.reverse1,
        hit.start1,
        hit.end1,
        &blocks[b1],
        cfg.split_read_leeway,
    ) {
        return Vec::new();
    }
    if !split_agrees(
        false,
        hit.reverse2,
        hit.start2,
        hit.end2,
        &blocks[b2],
        cfg.split_read_leeway,
    ) {
        return Vec::new();
    }
    let v1 = read_vertex(b1 as i32, hit.reverse1);
    let v2 = read_vertex(b2 as i32, hit.reverse2) ^ 1;
    if v1 == v2 {
        return Vec::new();
    }
    vec![(v1, v2)]
}

fn add_split(
    edges: &mut IndexMap<(i32, i32), EdgeRec>,
    hit: &SplitIn,
    blocks: &[BlockIn],
    cfg: &GraphConfig,
    in_likelihood: bool,
) {
    let edge_list = split_support_edges(hit, blocks, cfg);
    if edge_list.is_empty() {
        return;
    }
    let Some(b1) = find_block(blocks, hit.start1) else {
        return;
    };
    let Some(b2) = find_block(blocks, hit.start2) else {
        return;
    };
    let mut offset = read_offset(hit.start1, hit.end1, hit.reverse1, &blocks[b1]);
    offset += read_offset(hit.start2, hit.end2, hit.reverse2, &blocks[b2]);
    let mut scored = false;
    for (v1, v2) in edge_list {
        if v1 == v2 {
            continue;
        }
        let left_b = (v1 / 2).min(v2 / 2);
        let right_b = (v1 / 2).max(v2 / 2);
        if left_b != right_b
            && blocks[right_b as usize].start - blocks[left_b as usize].end > cfg.max_edge_span
        {
            continue;
        }
        let il = in_likelihood && !scored;
        add_edge(
            edges,
            v1,
            v2,
            1,
            Some(EdgeReadRec {
                offset,
                lib_id: 0,
                pmapped: cfg.pmapped,
                mapq: hit.mapq,
                kind: KIND_SPLIT,
                er_flags: er_flags(false, false, false, true, il),
                mate_len: 0,
                adj1: None,
                adj2: None,
            }),
            false,
        );
        scored = true;
    }
}

fn add_within_read_pe_support(
    edges: &mut IndexMap<(i32, i32), EdgeRec>,
    verts: &[i32],
    offset: i32,
    blocks: &[BlockIn],
    cfg: &GraphConfig,
    mapq: u8,
) {
    for (u, v) in within_read_path_edges(verts) {
        if u == v {
            continue;
        }
        let left_b = (u / 2).min(v / 2);
        let right_b = (u / 2).max(v / 2);
        if left_b != right_b
            && blocks[right_b as usize].start - blocks[left_b as usize].end > cfg.max_edge_span
        {
            continue;
        }
        add_edge(
            edges,
            u,
            v,
            1,
            Some(EdgeReadRec {
                offset,
                lib_id: 0,
                pmapped: cfg.pmapped,
                mapq,
                kind: KIND_PE,
                er_flags: er_flags(false, false, false, true, false),
                mate_len: 0,
                adj1: None,
                adj2: None,
            }),
            false,
        );
    }
}

fn prune(edges: &mut IndexMap<(i32, i32), EdgeRec>, blocks: &[BlockIn], min_support: i32, max_span: i64) {
    let keys: Vec<(i32, i32)> = edges.keys().copied().collect();
    let mut drop: Vec<(i32, i32)> = Vec::new();
    for key in keys {
        let edge = edges.get(&key).unwrap();
        if edge.support < min_support {
            if edge.reads.is_empty() {
                drop.push(key);
                continue;
            }
            let b1 = edge.u / 2;
            let b2 = edge.v / 2;
            if b1 != b2 {
                let left = b1.min(b2);
                let right = b1.max(b2);
                let span = blocks[right as usize].start - blocks[left as usize].end;
                if span > max_span {
                    drop.push(key);
                }
            }
            continue;
        }
        let b1 = edge.u / 2;
        let b2 = edge.v / 2;
        if b1 == b2 {
            continue;
        }
        let left = b1.min(b2);
        let right = b1.max(b2);
        if blocks[left as usize].is_insertion || blocks[right as usize].is_insertion {
            continue;
        }
        let span = blocks[right as usize].start - blocks[left as usize].end;
        if span > max_span {
            drop.push(key);
        }
    }
    for key in drop {
        edges.shift_remove(&key);
    }
}

/// Resolved end after used_splits serial claim (verts, offset, adj, from_split).
type ResolvedEnd = Option<(Vec<i32>, i32, Option<Vec<i32>>, bool)>;

/// Parallel end resolution given a pre-chosen split index (-1 = no split).
fn resolve_end_with_split(
    start: i64,
    end: i64,
    reverse: bool,
    split_idx: i32,
    splits: &[SplitIn],
    blocks: &[BlockIn],
    cfg: &GraphConfig,
) -> ResolvedEnd {
    if split_idx >= 0 {
        let hit = &splits[split_idx as usize];
        if let Some((verts, offset)) = split_path(hit, blocks, cfg.split_read_leeway) {
            let adj = adj_of(&verts);
            return Some((verts, offset, adj, true));
        }
    }
    let (verts, offset) = end_vertices(blocks, start, end, reverse)?;
    let adj = adj_of(&verts);
    Some((verts, offset, adj, false))
}

/// Serial pass: claim used_splits in pair order; return per-pair split idxs for each end.
fn claim_used_splits(
    pairs: &PairsIn,
    splits: &[SplitIn],
    splits_by_q: &IndexMap<i64, Vec<usize>>,
    blocks: &[BlockIn],
    cfg: &GraphConfig,
) -> (Vec<(i32, i32)>, Vec<bool>) {
    let n = pairs.len();
    let mut claims = vec![(-1i32, -1i32); n];
    let mut used = vec![false; splits.len()];
    for i in 0..n {
        let flags = pairs.flags[i];
        if flags & FLAG_HANGING != 0 {
            continue;
        }
        if pairs.mapq1[i] < cfg.min_mapq || pairs.mapq2[i] < cfg.min_mapq {
            continue;
        }
        let qh = pairs.qname_hash[i];
        let hits = match splits_by_q.get(&qh) {
            Some(h) => h.as_slice(),
            None => &[],
        };
        let rev1 = flags & FLAG_REVERSE1 != 0;
        let rev2 = flags & FLAG_REVERSE2 != 0;
        // left end
        let mut left_claim = -1i32;
        {
            let matched: Vec<usize> = hits
                .iter()
                .copied()
                .filter(|&si| {
                    splits[si].primary_start == pairs.start1[i]
                        && splits[si].primary_reverse == rev1
                })
                .collect();
            if matched.len() == 1 {
                let si = matched[0];
                if split_path(&splits[si], blocks, cfg.split_read_leeway).is_some() {
                    used[si] = true;
                    left_claim = si as i32;
                }
            }
        }
        // right end
        let mut right_claim = -1i32;
        {
            let matched: Vec<usize> = hits
                .iter()
                .copied()
                .filter(|&si| {
                    splits[si].primary_start == pairs.start2[i]
                        && splits[si].primary_reverse == rev2
                })
                .collect();
            if matched.len() == 1 {
                let si = matched[0];
                if split_path(&splits[si], blocks, cfg.split_read_leeway).is_some() {
                    used[si] = true;
                    right_claim = si as i32;
                }
            }
        }
        claims[i] = (left_claim, right_claim);
    }
    (claims, used)
}

fn process_pair_edges(
    pairs: &PairsIn,
    i: usize,
    left: &ResolvedEnd,
    right: &ResolvedEnd,
    blocks: &[BlockIn],
    cfg: &GraphConfig,
    edges: &mut IndexMap<(i32, i32), EdgeRec>,
) {
    let flags = pairs.flags[i];
    if flags & FLAG_HANGING != 0 {
        add_hanging(edges, pairs, i, blocks, cfg);
        return;
    }
    if pairs.mapq1[i] < cfg.min_mapq || pairs.mapq2[i] < cfg.min_mapq {
        return;
    }
    let (Some(l), Some(r)) = (left, right) else {
        return;
    };
    let (verts1, off1, adj1, from_split1) = l;
    let (verts2, off2, adj2, from_split2) = r;
    let v1 = verts1[0];
    let v2 = verts2[0];
    if v1 == v2 {
        return;
    }
    let b1 = v1 / 2;
    let b2 = v2 / 2;
    let span_lo = b1.min(b2);
    let span_hi = b1.max(b2);
    if span_lo != span_hi {
        let span = blocks[span_hi as usize].start - blocks[span_lo as usize].end;
        if span > cfg.max_edge_span {
            return;
        }
    }
    let rev1 = flags & FLAG_REVERSE1 != 0;
    let rev2 = flags & FLAG_REVERSE2 != 0;
    let orient_ok = v1 == read_vertex(b1, rev1) && v2 == read_vertex(b2, rev2);
    let junction = between_junction(verts1, verts2);
    let store_key = (v1.min(v2), v1.max(v2));
    let junction_key = junction.map(|(a, b)| (a.min(b), a.max(b)));
    let supported = if junction_key == Some(store_key) {
        if orient_ok {
            pair_supported(
                pairs.start1[i],
                pairs.end1[i],
                rev1,
                pairs.start2[i],
                pairs.end2[i],
                rev2,
                blocks,
                v1,
                v2,
                cfg,
            )
        } else {
            false
        }
    } else {
        false
    };
    let implied = if orient_ok {
        adjacency_insert(
            pairs.start1[i],
            pairs.end1[i],
            rev1,
            pairs.start2[i],
            pairs.end2[i],
            rev2,
            blocks,
            b1,
            b2,
        )
    } else {
        pairs.insert[i]
    };
    let (adj_lo, adj_hi) = if v1 < v2 {
        (adj1.clone(), adj2.clone())
    } else {
        (adj2.clone(), adj1.clone())
    };
    let lib_id = if cfg.pool_libraries {
        0
    } else {
        pairs.lib_id[i]
    };
    let discordant = !(cfg.q_low <= implied && implied <= cfg.q_high);
    add_edge(
        edges,
        v1,
        v2,
        if supported { 1 } else { 0 },
        Some(EdgeReadRec {
            offset: off1 + off2,
            lib_id,
            pmapped: cfg.pmapped,
            mapq: pairs.mapq1[i].min(pairs.mapq2[i]),
            kind: KIND_PE,
            er_flags: er_flags(false, discordant, false, supported, true),
            mate_len: 0,
            adj1: adj_lo,
            adj2: adj_hi,
        }),
        false,
    );
    if let Some((ju, jv)) = junction {
        let jk = (ju.min(jv), ju.max(jv));
        if jk != store_key
            && ju != jv
            && between_junction_supported(
                pairs.start1[i],
                pairs.end1[i],
                rev1,
                pairs.start2[i],
                pairs.end2[i],
                rev2,
                verts1,
                verts2,
                *off1,
                *off2,
                blocks,
                cfg,
            )
        {
            let left_b = (ju / 2).min(jv / 2);
            let right_b = (ju / 2).max(jv / 2);
            if left_b == right_b
                || blocks[right_b as usize].start - blocks[left_b as usize].end
                    <= cfg.max_edge_span
            {
                add_edge(
                    edges,
                    ju,
                    jv,
                    1,
                    Some(EdgeReadRec {
                        offset: off1 + off2,
                        lib_id,
                        pmapped: cfg.pmapped,
                        mapq: pairs.mapq1[i].min(pairs.mapq2[i]),
                        kind: KIND_PE,
                        er_flags: er_flags(false, false, false, true, false),
                        mate_len: 0,
                        adj1: None,
                        adj2: None,
                    }),
                    false,
                );
            }
        }
    }
    if !from_split1 && verts1.len() > 1 {
        add_within_read_pe_support(edges, verts1, *off1, blocks, cfg, pairs.mapq1[i]);
    }
    if !from_split2 && verts2.len() > 1 {
        add_within_read_pe_support(edges, verts2, *off2, blocks, cfg, pairs.mapq2[i]);
    }
}

fn edges_to_columns(
    edges: IndexMap<(i32, i32), EdgeRec>,
    block_clip_pos: Vec<i64>,
    block_clip_reads: Vec<i32>,
) -> GraphColumns {
    let n_edges = edges.len();
    let mut edge_u = Vec::with_capacity(n_edges);
    let mut edge_v = Vec::with_capacity(n_edges);
    let mut edge_support = Vec::with_capacity(n_edges);
    let mut edge_reference = Vec::with_capacity(n_edges);
    let mut edge_read_ptr = Vec::with_capacity(n_edges + 1);
    edge_read_ptr.push(0);
    let mut read_offset = Vec::new();
    let mut read_lib_id = Vec::new();
    let mut read_pmapped = Vec::new();
    let mut read_mapq = Vec::new();
    let mut read_kind = Vec::new();
    let mut read_flags = Vec::new();
    let mut read_mate_len = Vec::new();
    let mut read_adj1_ptr = Vec::new();
    let mut read_adj1_len = Vec::new();
    let mut read_adj2_ptr = Vec::new();
    let mut read_adj2_len = Vec::new();
    let mut adj_flat = Vec::new();

    for (_, edge) in edges {
        edge_u.push(edge.u);
        edge_v.push(edge.v);
        edge_support.push(edge.support);
        edge_reference.push(if edge.reference { 1 } else { 0 });
        for r in edge.reads {
            read_offset.push(r.offset);
            read_lib_id.push(r.lib_id);
            read_pmapped.push(r.pmapped);
            read_mapq.push(r.mapq);
            read_kind.push(r.kind);
            read_flags.push(r.er_flags);
            read_mate_len.push(r.mate_len);
            match r.adj1 {
                Some(a) if !a.is_empty() => {
                    read_adj1_ptr.push(adj_flat.len() as i32);
                    read_adj1_len.push(a.len() as i32);
                    adj_flat.extend(a);
                }
                _ => {
                    read_adj1_ptr.push(-1);
                    read_adj1_len.push(0);
                }
            }
            match r.adj2 {
                Some(a) if !a.is_empty() => {
                    read_adj2_ptr.push(adj_flat.len() as i32);
                    read_adj2_len.push(a.len() as i32);
                    adj_flat.extend(a);
                }
                _ => {
                    read_adj2_ptr.push(-1);
                    read_adj2_len.push(0);
                }
            }
        }
        edge_read_ptr.push(read_offset.len() as i32);
    }

    GraphColumns {
        edge_u,
        edge_v,
        edge_support,
        edge_reference,
        edge_read_ptr,
        read_offset,
        read_lib_id,
        read_pmapped,
        read_mapq,
        read_kind,
        read_flags,
        read_mate_len,
        read_adj1_ptr,
        read_adj1_len,
        read_adj2_ptr,
        read_adj2_len,
        adj_flat,
        block_clip_pos,
        block_clip_reads,
    }
}

/// Stamp clip positions from breakpoint clip_pos arrays onto blocks (mutates copies).
pub fn stamp_clip_positions(
    blocks: &mut [BlockIn],
    bp_clip_pos: &[i64],
    bp_support: &[i32],
) {
    for (&clip_pos, &support) in bp_clip_pos.iter().zip(bp_support.iter()) {
        if clip_pos < 0 {
            continue;
        }
        for i in 0..blocks.len().saturating_sub(1) {
            if blocks[i].is_insertion || blocks[i + 1].is_insertion {
                continue;
            }
            if blocks[i].end <= clip_pos
                && clip_pos <= blocks[i + 1].start
                && support >= blocks[i].clip_reads
            {
                blocks[i].clip_pos = clip_pos;
                blocks[i].clip_reads = support;
            }
        }
    }
}

/// Stamp from softclip positions (sorted). Mode of positions in each gap.
pub fn stamp_clip_positions_from_softclips(blocks: &mut [BlockIn], positions: &[i64]) {
    if positions.is_empty() {
        return;
    }
    // positions must be sorted
    for i in 0..blocks.len().saturating_sub(1) {
        if blocks[i].is_insertion || blocks[i + 1].is_insertion {
            continue;
        }
        let lo = blocks[i].end;
        let hi = blocks[i + 1].start;
        if hi <= lo {
            continue;
        }
        let left = positions.partition_point(|&p| p < lo);
        let right = positions.partition_point(|&p| p <= hi);
        if left >= right {
            continue;
        }
        // mode (most common); ties -> first encountered in Counter.most_common order
        // Python Counter.most_common: highest count, ties broken by first-seen order.
        let slice = &positions[left..right];
        let mut best_pos = slice[0];
        let mut best_count = 0i32;
        let mut j = 0;
        while j < slice.len() {
            let pos = slice[j];
            let mut k = j + 1;
            while k < slice.len() && slice[k] == pos {
                k += 1;
            }
            let count = (k - j) as i32;
            // first max wins (sorted positions: equal values contiguous; among
            // different values with same count, Python Counter keeps insertion
            // order of first occurrence which for sorted input is ascending).
            // Counter.most_common uses heap by (-count, first_index).
            // For sorted input, first occurrence order = ascending pos order.
            // So among ties, smaller pos that appeared first wins — which is
            // the first we see when scanning left-to-right with > not >=.
            if count > best_count {
                best_count = count;
                best_pos = pos;
            }
            j = k;
        }
        if best_count >= blocks[i].clip_reads {
            blocks[i].clip_pos = best_pos;
            blocks[i].clip_reads = best_count;
        }
    }
}

pub fn build_graph_columns(
    mut blocks: Vec<BlockIn>,
    pairs: PairsIn<'_>,
    splits: &[SplitIn],
    cfg: &GraphConfig,
    bp_clip_pos: &[i64],
    bp_support: &[i32],
    softclip_positions: &[i64],
) -> GraphColumns {
    stamp_clip_positions(&mut blocks, bp_clip_pos, bp_support);
    let mut sc = softclip_positions.to_vec();
    sc.sort_unstable();
    stamp_clip_positions_from_softclips(&mut blocks, &sc);

    let block_clip_pos: Vec<i64> = blocks.iter().map(|b| b.clip_pos).collect();
    let block_clip_reads: Vec<i32> = blocks.iter().map(|b| b.clip_reads).collect();

    let mut edges: IndexMap<(i32, i32), EdgeRec> = IndexMap::new();
    add_reference_edges(&mut edges, &blocks, cfg.min_edge_support);

    let mut splits_by_q: IndexMap<i64, Vec<usize>> = IndexMap::new();
    for (si, hit) in splits.iter().enumerate() {
        if hit.qname_hash < 0 {
            continue;
        }
        splits_by_q
            .entry(hit.qname_hash)
            .or_default()
            .push(si);
    }

    // Serial used_splits claims (pair order).
    let (claims, used) = claim_used_splits(&pairs, splits, &splits_by_q, &blocks, cfg);

    // Parallel resolve ends in chunks, then serial merge (deterministic).
    // Collecting all 29M ResolvedEnd rows at once spiked RSS by several GB on chr1.
    const RESOLVE_CHUNK: usize = 250_000;
    let n = pairs.len();
    let mut chunk_start = 0usize;
    while chunk_start < n {
        let chunk_end = (chunk_start + RESOLVE_CHUNK).min(n);
        let resolved: Vec<(ResolvedEnd, ResolvedEnd)> = (chunk_start..chunk_end)
            .into_par_iter()
            .map(|i| {
                let flags = pairs.flags[i];
                if flags & FLAG_HANGING != 0
                    || pairs.mapq1[i] < cfg.min_mapq
                    || pairs.mapq2[i] < cfg.min_mapq
                {
                    return (None, None);
                }
                let (lc, rc) = claims[i];
                let rev1 = flags & FLAG_REVERSE1 != 0;
                let rev2 = flags & FLAG_REVERSE2 != 0;
                let left = resolve_end_with_split(
                    pairs.start1[i],
                    pairs.end1[i],
                    rev1,
                    lc,
                    splits,
                    &blocks,
                    cfg,
                );
                let right = resolve_end_with_split(
                    pairs.start2[i],
                    pairs.end2[i],
                    rev2,
                    rc,
                    splits,
                    &blocks,
                    cfg,
                );
                (left, right)
            })
            .collect();
        for (offset, (left, right)) in resolved.into_iter().enumerate() {
            let i = chunk_start + offset;
            process_pair_edges(&pairs, i, &left, &right, &blocks, cfg, &mut edges);
        }
        chunk_start = chunk_end;
    }

    for (si, hit) in splits.iter().enumerate() {
        add_split(&mut edges, hit, &blocks, cfg, !used[si]);
    }

    prune(
        &mut edges,
        &blocks,
        cfg.min_edge_support,
        cfg.max_edge_span,
    );

    edges_to_columns(edges, block_clip_pos, block_clip_reads)
}
