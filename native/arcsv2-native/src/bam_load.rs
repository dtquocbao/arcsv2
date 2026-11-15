//! BAM/CRAM → columnar Pairs SoA + take-records for softclip/split discovery.
//!
//! Uses rust-htslib (system/bundled htslib) for CRAM + multithreaded BGZF.
//! Mate pairing is always a single ordered pass over the BAM/htslib fetch
//! stream (same encounter order as Python ``samio.load_region``).
//! ``LoadConfig.threads`` only sets htslib BGZF decompression threads — it
//! must not shard the chrom into rayon region fetches. Region-split loading
//! drops unmapped mates that chrom fetch returns, so hanging-pair fields
//! diverge (same n_pairs, different multiset / order).

use crate::cigar::{
    aligned_length, cigar_has_softclip, query_len_from_cigar, softclip_counts,
};
use crate::hashutil::{lib_id, qhash};
use indexmap::IndexMap;
use rust_htslib::bam::{self, record::Cigar, Read};
use std::path::Path;

pub const FLAG_UNMAPPED: u16 = 0x4;
pub const FLAG_MATE_UNMAPPED: u16 = 0x8;
pub const FLAG_REVERSE: u16 = 0x10;
pub const FLAG_SECONDARY: u16 = 0x100;
pub const FLAG_QCFAIL: u16 = 0x200;
pub const FLAG_DUP: u16 = 0x400;
pub const FLAG_SUPPLEMENTARY: u16 = 0x800;
pub const SKIP_FLAGS: u16 = FLAG_SECONDARY | FLAG_DUP | FLAG_QCFAIL | FLAG_SUPPLEMENTARY;

pub const FLAG_REVERSE1: u8 = 1 << 0;
pub const FLAG_REVERSE2: u8 = 1 << 1;
pub const FLAG_HANGING: u8 = 1 << 2;
pub const FLAG_DISTANT: u8 = 1 << 3;

#[derive(Clone, Debug)]
pub struct CompactRead {
    pub qname: String,
    pub flag: u16,
    pub rname: String,
    pub start: i64,
    pub mapq: u8,
    pub cigar: String,
    pub qual: String,
    pub sa: String,
    pub rg_key: String,
    pub nm: Option<i64>,
}

#[derive(Clone, Debug, Default)]
pub struct PairRow {
    pub start1: i64,
    pub end1: i64,
    pub start2: i64,
    pub end2: i64,
    pub insert: i32,
    pub mate_len: i32,
    pub qname_hash: i64,
    pub mapq1: u8,
    pub mapq2: u8,
    pub lib_id: u16,
    pub flags: u8,
}

#[derive(Clone, Debug)]
pub struct TakeRecord {
    pub qname: String,
    pub rname: String,
    pub start: i64,
    pub mapq: u8,
    pub cigar: String,
    pub reverse: bool,
    pub sa: String,
    pub qual: String,
}

#[derive(Clone, Debug)]
pub struct LoadConfig {
    pub max_pair_distance: i64,
    pub threads: usize,
    pub filter_high_nm: bool,
    pub max_nm: i64,
    pub filter_both_mates_clipped: bool,
}

impl Default for LoadConfig {
    fn default() -> Self {
        Self {
            max_pair_distance: 2_000_000,
            threads: 1,
            filter_high_nm: false,
            max_nm: 5,
            filter_both_mates_clipped: false,
        }
    }
}

#[derive(Debug)]
pub struct LoadResult {
    pub pairs: Vec<PairRow>,
    pub takes: Vec<TakeRecord>,
    pub chrom_name: String,
    pub chrom_end: i64,
    pub n_aln: u64,
    pub n_insert_stat: u64,
}

fn cigar_string(rec: &bam::Record) -> String {
    let cig = rec.cigar();
    if cig.is_empty() {
        return "*".to_string();
    }
    let mut s = String::new();
    for c in cig.iter() {
        match c {
            Cigar::Match(n)
            | Cigar::Ins(n)
            | Cigar::Del(n)
            | Cigar::RefSkip(n)
            | Cigar::SoftClip(n)
            | Cigar::HardClip(n)
            | Cigar::Pad(n)
            | Cigar::Equal(n)
            | Cigar::Diff(n) => {
                let op = match c {
                    Cigar::Match(_) => 'M',
                    Cigar::Ins(_) => 'I',
                    Cigar::Del(_) => 'D',
                    Cigar::RefSkip(_) => 'N',
                    Cigar::SoftClip(_) => 'S',
                    Cigar::HardClip(_) => 'H',
                    Cigar::Pad(_) => 'P',
                    Cigar::Equal(_) => '=',
                    Cigar::Diff(_) => 'X',
                };
                s.push_str(&n.to_string());
                s.push(op);
            }
        }
    }
    if s.is_empty() {
        "*".to_string()
    } else {
        s
    }
}

fn tag_str(rec: &bam::Record, tag: &[u8; 2]) -> Option<String> {
    use rust_htslib::bam::record::Aux;
    match rec.aux(tag) {
        Ok(Aux::String(s)) => Some(s.to_string()),
        _ => None,
    }
}

fn tag_nm(rec: &bam::Record) -> Option<i64> {
    use rust_htslib::bam::record::Aux;
    match rec.aux(b"NM") {
        Ok(Aux::U8(v)) => Some(v as i64),
        Ok(Aux::U16(v)) => Some(v as i64),
        Ok(Aux::U32(v)) => Some(v as i64),
        Ok(Aux::I8(v)) => Some(v as i64),
        Ok(Aux::I16(v)) => Some(v as i64),
        Ok(Aux::I32(v)) => Some(v as i64),
        _ => None,
    }
}

fn qual_string(rec: &bam::Record) -> String {
    let q = rec.qual();
    if q.is_empty() {
        return String::new();
    }
    let mut s = String::with_capacity(q.len());
    for &phred in q {
        s.push(char::from(phred.saturating_add(33)));
    }
    s
}

fn aligned_end(read: &CompactRead) -> i64 {
    read.start + aligned_length(&read.cigar)
}

fn in_window(read: &CompactRead, start: Option<i64>, end: Option<i64>) -> bool {
    if start.is_none() && end.is_none() {
        return true;
    }
    let read_end = aligned_end(read);
    if let Some(s) = start {
        if read_end <= s {
            return false;
        }
    }
    if let Some(e) = end {
        if read.start >= e {
            return false;
        }
    }
    true
}

fn on_target(
    read: &CompactRead,
    chrom: Option<&str>,
    start: Option<i64>,
    end: Option<i64>,
) -> bool {
    if read.flag & FLAG_UNMAPPED != 0 {
        return false;
    }
    if let Some(c) = chrom {
        if read.rname != c {
            return false;
        }
    }
    in_window(read, start, end)
}

fn insert_length(a: &CompactRead, b: &CompactRead) -> Option<i64> {
    let a_rev = a.flag & FLAG_REVERSE != 0;
    let b_rev = b.flag & FLAG_REVERSE != 0;
    if a_rev == b_rev {
        return None;
    }
    let which_minus: usize = if a_rev { 0 } else { 1 };
    let which_first = 1 - which_minus;
    let which_last = 1 - which_first;
    let pair = [a, b];
    let first = pair[which_first];
    let last = pair[which_last];
    let _first_end = first.start + aligned_length(&first.cigar);
    let last_end = last.start + aligned_length(&last.cigar);
    let mut ilen = last_end - first.start;
    let (first_left, _) = softclip_counts(&first.cigar);
    let (_, last_right) = softclip_counts(&last.cigar);
    ilen += first_left + last_right;
    Some(ilen)
}

fn pair_flags(reverse1: bool, reverse2: bool, hanging: bool, distant: bool) -> u8 {
    let mut f = 0u8;
    if reverse1 {
        f |= FLAG_REVERSE1;
    }
    if reverse2 {
        f |= FLAG_REVERSE2;
    }
    if hanging {
        f |= FLAG_HANGING;
    }
    if distant {
        f |= FLAG_DISTANT;
    }
    f
}

fn push_pair_fields(
    pairs: &mut Vec<PairRow>,
    start1: i64,
    end1: i64,
    start2: i64,
    end2: i64,
    reverse1: bool,
    reverse2: bool,
    mapq1: u8,
    mapq2: u8,
    insert: i32,
    lib: u16,
    hanging: bool,
    qh: i64,
    distant: bool,
    mate_len: i32,
) {
    pairs.push(PairRow {
        start1,
        end1,
        start2,
        end2,
        insert,
        mate_len,
        qname_hash: qh,
        mapq1,
        mapq2,
        lib_id: lib,
        flags: pair_flags(reverse1, reverse2, hanging, distant),
    });
}

fn hanging_push(
    pairs: &mut Vec<PairRow>,
    anchor: &CompactRead,
    mate: &CompactRead,
    distant: bool,
) {
    let mate_unmapped = mate.flag & FLAG_UNMAPPED != 0;
    let insert = if mate_unmapped {
        0
    } else {
        match insert_length(anchor, mate) {
            Some(v) => v as i32,
            None => 0,
        }
    };
    push_pair_fields(
        pairs,
        anchor.start,
        anchor.start + aligned_length(&anchor.cigar),
        if mate_unmapped { 0 } else { mate.start },
        if mate_unmapped {
            0
        } else {
            mate.start + aligned_length(&mate.cigar)
        },
        anchor.flag & FLAG_REVERSE != 0,
        mate.flag & FLAG_REVERSE != 0,
        anchor.mapq,
        if mate_unmapped { 0 } else { mate.mapq },
        insert,
        lib_id(&anchor.rg_key),
        true,
        qhash(&anchor.qname),
        distant,
        query_len_from_cigar(&mate.cigar) as i32,
    );
}

fn maybe_take(takes: &mut Vec<TakeRecord>, read: &CompactRead) {
    let has_sa = !read.sa.is_empty();
    let has_clip = cigar_has_softclip(&read.cigar);
    if !has_sa && !has_clip {
        return;
    }
    takes.push(TakeRecord {
        qname: read.qname.clone(),
        rname: read.rname.clone(),
        start: read.start,
        mapq: read.mapq,
        cigar: read.cigar.clone(),
        reverse: read.flag & FLAG_REVERSE != 0,
        sa: read.sa.clone(),
        qual: read.qual.clone(),
    });
}

struct Emitter<'a> {
    cfg: &'a LoadConfig,
    chrom: Option<&'a str>,
    start: Option<i64>,
    end: Option<i64>,
    pairs: Vec<PairRow>,
    takes: Vec<TakeRecord>,
    chrom_name: String,
    chrom_end: i64,
    n_insert_stat: u64,
}

impl<'a> Emitter<'a> {
    fn new(
        cfg: &'a LoadConfig,
        chrom: Option<&'a str>,
        start: Option<i64>,
        end: Option<i64>,
    ) -> Self {
        Self {
            cfg,
            chrom,
            start,
            end,
            pairs: Vec::new(),
            takes: Vec::new(),
            chrom_name: chrom.unwrap_or("").to_string(),
            chrom_end: 0,
            n_insert_stat: 0,
        }
    }

    fn note_read(&mut self, read: &CompactRead) {
        if self.chrom_name.is_empty() {
            self.chrom_name = read.rname.clone();
        }
        self.chrom_end = self.chrom_end.max(aligned_end(read));
    }

    fn emit_pair(&mut self, mut a: CompactRead, mut b: CompactRead) {
        if a.flag & FLAG_REVERSE != 0 && b.flag & FLAG_REVERSE == 0 {
            std::mem::swap(&mut a, &mut b);
        }
        let ilen = insert_length(&a, &b);
        if let Some(v) = ilen {
            if v >= 0 {
                self.n_insert_stat += 1;
            }
        }
        if self.cfg.filter_both_mates_clipped
            && cigar_has_softclip(&a.cigar)
            && cigar_has_softclip(&b.cigar)
        {
            return;
        }
        push_pair_fields(
            &mut self.pairs,
            a.start,
            a.start + aligned_length(&a.cigar),
            b.start,
            b.start + aligned_length(&b.cigar),
            a.flag & FLAG_REVERSE != 0,
            b.flag & FLAG_REVERSE != 0,
            a.mapq,
            b.mapq,
            ilen.unwrap_or(0) as i32,
            lib_id(&a.rg_key),
            false,
            qhash(&a.qname),
            false,
            0,
        );
        maybe_take(&mut self.takes, &a);
        maybe_take(&mut self.takes, &b);
    }

    fn emit_joined(&mut self, a: CompactRead, b: CompactRead) {
        let a_ok = on_target(&a, self.chrom, self.start, self.end);
        let b_ok = on_target(&b, self.chrom, self.start, self.end);
        if !a_ok && !b_ok {
            return;
        }
        let a_unmapped = a.flag & FLAG_UNMAPPED != 0;
        let b_unmapped = b.flag & FLAG_UNMAPPED != 0;
        if a_unmapped && b_unmapped {
            return;
        }
        if a_unmapped || b_unmapped {
            let (anchor, mate) = if a_unmapped { (&b, &a) } else { (&a, &b) };
            if on_target(anchor, self.chrom, self.start, self.end) {
                hanging_push(&mut self.pairs, anchor, mate, false);
                maybe_take(&mut self.takes, anchor);
            }
            return;
        }
        let distant = a.rname != b.rname
            || (a.start - b.start).abs() > self.cfg.max_pair_distance;
        if distant {
            if a_ok {
                hanging_push(&mut self.pairs, &a, &b, true);
                maybe_take(&mut self.takes, &a);
            }
            if b_ok {
                hanging_push(&mut self.pairs, &b, &a, true);
                maybe_take(&mut self.takes, &b);
            }
            return;
        }
        self.emit_pair(a, b);
    }

    fn finish_pending(&mut self, pending: IndexMap<String, CompactRead>) {
        // IndexMap preserves insertion order like Python 3.7+ dict.
        for (_, read) in pending {
            if !on_target(&read, self.chrom, self.start, self.end) {
                continue;
            }
            maybe_take(&mut self.takes, &read);
            if read.flag & FLAG_MATE_UNMAPPED == 0 {
                continue;
            }
            push_pair_fields(
                &mut self.pairs,
                read.start,
                read.start + aligned_length(&read.cigar),
                0,
                0,
                read.flag & FLAG_REVERSE != 0,
                false,
                read.mapq,
                0,
                0,
                lib_id(&read.rg_key),
                true,
                qhash(&read.qname),
                false,
                0,
            );
        }
    }
}

fn record_to_compact(
    rec: &bam::Record,
    header: &bam::HeaderView,
    max_nm: Option<i64>,
) -> Option<CompactRead> {
    let flag = rec.flags();
    if flag & SKIP_FLAGS != 0 {
        return None;
    }
    let nm = tag_nm(rec);
    if let Some(cap) = max_nm {
        if let Some(v) = nm {
            if v > cap {
                return None;
            }
        }
    }
    let unmapped = flag & FLAG_UNMAPPED != 0;
    let tid = rec.tid();
    let rname = if unmapped || tid < 0 {
        "*".to_string()
    } else {
        String::from_utf8_lossy(header.tid2name(tid as u32)).into_owned()
    };
    let start = if unmapped || rec.pos() < 0 {
        0
    } else {
        rec.pos()
    };
    let cigar = cigar_string(rec);
    let want_qual = cigar_has_softclip(&cigar);
    let qual = if want_qual {
        qual_string(rec)
    } else {
        String::new()
    };
    let sa = tag_str(rec, b"SA").unwrap_or_default();
    let rg_key = tag_str(rec, b"LB")
        .or_else(|| tag_str(rec, b"RG"))
        .unwrap_or_else(|| "lib1".to_string());
    let qname = String::from_utf8_lossy(rec.qname()).into_owned();
    Some(CompactRead {
        qname,
        flag,
        rname,
        start,
        mapq: rec.mapq(),
        cigar,
        qual,
        sa,
        rg_key,
        nm,
    })
}

fn open_reader(
    path: &Path,
    reference: Option<&str>,
    threads: usize,
) -> Result<bam::IndexedReader, String> {
    let mut reader = if let Some(fasta) = reference {
        // CRAM path: rust-htslib IndexedReader with reference.
        let mut r = bam::IndexedReader::from_path(path).map_err(|e| e.to_string())?;
        r.set_reference(Path::new(fasta)).map_err(|e| e.to_string())?;
        r
    } else {
        bam::IndexedReader::from_path(path).map_err(|e| e.to_string())?
    };
    if threads > 1 {
        reader.set_threads(threads).map_err(|e| e.to_string())?;
    }
    Ok(reader)
}

/// Load a BAM/CRAM region into columnar pairs.
///
/// Single ordered htslib fetch + mate-pair pass (Python stream order).
/// ``cfg.threads`` is BGZF decompression only (like pysam ``threads=``).
/// Region-sharding was removed: it cannot keep unmapped mates in BAM order.
///
/// M4: stream records into the mate map — do not buffer the whole chromosome
/// as ``Vec<CompactRead>`` (that alone was ~20 GB RSS on chr1).
pub fn load_bam_region(
    path: &str,
    chrom: Option<&str>,
    start: Option<i64>,
    end: Option<i64>,
    reference: Option<&str>,
    cfg: &LoadConfig,
) -> Result<LoadResult, String> {
    let path = Path::new(path);
    let max_nm = if cfg.filter_high_nm {
        Some(cfg.max_nm)
    } else {
        None
    };
    let bam_threads = cfg.threads.max(1);
    let mut reader = open_reader(path, reference, bam_threads)?;
    let header = reader.header().clone();
    if let Some(c) = chrom {
        let s = start.unwrap_or(0);
        if let Some(e) = end {
            let region = format!("{}:{}-{}", c, s + 1, e);
            reader.fetch(&region).map_err(|e| e.to_string())?;
        } else {
            reader.fetch(c).map_err(|e| e.to_string())?;
        }
    } else {
        return Err("rust loader requires an explicit chromosome".into());
    }

    let mut emitter = Emitter::new(cfg, chrom, start, end);
    let mut pending: IndexMap<String, CompactRead> = IndexMap::new();
    let mut n_aln = 0u64;
    for result in reader.records() {
        let rec = result.map_err(|e| e.to_string())?;
        let Some(compact) = record_to_compact(&rec, &header, max_nm) else {
            continue;
        };
        emitter.note_read(&compact);
        n_aln += 1;
        let qn = compact.qname.clone();
        if let Some(earlier) = pending.shift_remove(&qn) {
            emitter.emit_joined(earlier, compact);
        } else {
            pending.insert(qn, compact);
        }
    }
    emitter.finish_pending(pending);
    let mut result = LoadResult {
        pairs: emitter.pairs,
        takes: emitter.takes,
        chrom_name: emitter.chrom_name,
        chrom_end: emitter.chrom_end,
        n_aln,
        n_insert_stat: emitter.n_insert_stat,
    };
    if result.chrom_name.is_empty() {
        result.chrom_name = chrom.unwrap_or("").to_string();
    }
    Ok(result)
}
