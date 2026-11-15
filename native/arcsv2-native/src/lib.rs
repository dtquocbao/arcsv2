//! Optional Rust kernels for ARC-SV2.
//!
//! R1: BAM/CRAM -> columnar Pairs SoA (bit-identical to samio.load_region).
//! R2: block prep + build_graph -> CSR edges (bit-identical to caller.build_graph).
//! Float reductions stay ordered (no rayon unordered sum). No fast-math.

mod bam_load;
mod cigar;
mod graph_build;
mod hashutil;

use bam_load::{load_bam_region, LoadConfig, LoadResult};
use graph_build::{
    build_graph_columns, BlockIn, GraphConfig, PairsIn, SplitIn, ER_DISCORDANT, ER_DISTANT,
    ER_HANGING, ER_IN_LIKELIHOOD, ER_SUPPORTS_ADJ, KIND_PE, KIND_SPLIT,
};
use numpy::{PyArray1, PyReadonlyArray1};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

/// Package version string (matches Cargo.toml).
#[pyfunction]
fn native_version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

/// True when this extension was built with the rust engine path available.
#[pyfunction]
fn rust_engine_available() -> bool {
    true
}

fn load_result_to_py(py: Python<'_>, mut result: LoadResult) -> PyResult<PyObject> {
    let out = PyDict::new_bound(py);

    // Convert pairs first so peak RSS is not (pair rows + take dicts).
    let pairs = std::mem::take(&mut result.pairs);
    let n = pairs.len();
    let mut start1 = Vec::with_capacity(n);
    let mut end1 = Vec::with_capacity(n);
    let mut start2 = Vec::with_capacity(n);
    let mut end2 = Vec::with_capacity(n);
    let mut insert = Vec::with_capacity(n);
    let mut mate_len = Vec::with_capacity(n);
    let mut qname_hash = Vec::with_capacity(n);
    let mut mapq1 = Vec::with_capacity(n);
    let mut mapq2 = Vec::with_capacity(n);
    let mut lib_id = Vec::with_capacity(n);
    let mut flags = Vec::with_capacity(n);
    for p in pairs {
        start1.push(p.start1);
        end1.push(p.end1);
        start2.push(p.start2);
        end2.push(p.end2);
        insert.push(p.insert);
        mate_len.push(p.mate_len);
        qname_hash.push(p.qname_hash);
        mapq1.push(p.mapq1);
        mapq2.push(p.mapq2);
        lib_id.push(p.lib_id);
        flags.push(p.flags);
    }
    out.set_item("start1", PyArray1::from_vec_bound(py, start1))?;
    out.set_item("end1", PyArray1::from_vec_bound(py, end1))?;
    out.set_item("start2", PyArray1::from_vec_bound(py, start2))?;
    out.set_item("end2", PyArray1::from_vec_bound(py, end2))?;
    out.set_item("insert", PyArray1::from_vec_bound(py, insert))?;
    out.set_item("mate_len", PyArray1::from_vec_bound(py, mate_len))?;
    out.set_item("qname_hash", PyArray1::from_vec_bound(py, qname_hash))?;
    out.set_item("mapq1", PyArray1::from_vec_bound(py, mapq1))?;
    out.set_item("mapq2", PyArray1::from_vec_bound(py, mapq2))?;
    out.set_item("lib_id", PyArray1::from_vec_bound(py, lib_id))?;
    out.set_item("flags", PyArray1::from_vec_bound(py, flags))?;
    out.set_item("chrom_name", result.chrom_name)?;
    out.set_item("chrom_end", result.chrom_end)?;
    out.set_item("n_aln", result.n_aln)?;
    out.set_item("n_insert_stat", result.n_insert_stat)?;

    // M4: columnar takes (no per-take PyDict). Parallel string lists + arrays.
    let takes_rust = std::mem::take(&mut result.takes);
    let nt = takes_rust.len();
    let mut t_start = Vec::with_capacity(nt);
    let mut t_mapq = Vec::with_capacity(nt);
    let mut t_reverse = Vec::with_capacity(nt);
    let t_qname = PyList::empty_bound(py);
    let t_rname = PyList::empty_bound(py);
    let t_cigar = PyList::empty_bound(py);
    let t_sa = PyList::empty_bound(py);
    let t_qual = PyList::empty_bound(py);
    for t in takes_rust {
        t_start.push(t.start);
        t_mapq.push(t.mapq);
        t_reverse.push(u8::from(t.reverse));
        t_qname.append(t.qname)?;
        t_rname.append(t.rname)?;
        t_cigar.append(t.cigar)?;
        t_sa.append(t.sa)?;
        t_qual.append(t.qual)?;
    }
    out.set_item("take_start", PyArray1::from_vec_bound(py, t_start))?;
    out.set_item("take_mapq", PyArray1::from_vec_bound(py, t_mapq))?;
    out.set_item("take_reverse", PyArray1::from_vec_bound(py, t_reverse))?;
    out.set_item("take_qname", t_qname)?;
    out.set_item("take_rname", t_rname)?;
    out.set_item("take_cigar", t_cigar)?;
    out.set_item("take_sa", t_sa)?;
    out.set_item("take_qual", t_qual)?;
    out.set_item("n_takes", nt)?;
    Ok(out.into())
}

/// Load BAM/CRAM into M2 Pairs SoA arrays plus softclip/split take-records.
#[pyfunction]
#[pyo3(signature = (
    path,
    chrom=None,
    start=None,
    end=None,
    reference=None,
    max_pair_distance=2_000_000,
    threads=1,
    filter_high_nm=false,
    max_nm=5,
    filter_both_mates_clipped=false,
))]
fn load_bam_pairs(
    py: Python<'_>,
    path: &str,
    chrom: Option<&str>,
    start: Option<i64>,
    end: Option<i64>,
    reference: Option<&str>,
    max_pair_distance: i64,
    threads: usize,
    filter_high_nm: bool,
    max_nm: i64,
    filter_both_mates_clipped: bool,
) -> PyResult<PyObject> {
    let cfg = LoadConfig {
        max_pair_distance,
        threads: threads.max(1),
        filter_high_nm,
        max_nm,
        filter_both_mates_clipped,
    };
    let result = py
        .allow_threads(|| load_bam_region(path, chrom, start, end, reference, &cfg))
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e))?;
    load_result_to_py(py, result)
}

fn arr_i64(a: PyReadonlyArray1<'_, i64>) -> PyResult<Vec<i64>> {
    Ok(a.as_slice()?.to_vec())
}
fn arr_i32(a: PyReadonlyArray1<'_, i32>) -> PyResult<Vec<i32>> {
    Ok(a.as_slice()?.to_vec())
}
fn arr_u8(a: PyReadonlyArray1<'_, u8>) -> PyResult<Vec<u8>> {
    Ok(a.as_slice()?.to_vec())
}
fn arr_u16(a: PyReadonlyArray1<'_, u16>) -> PyResult<Vec<u16>> {
    Ok(a.as_slice()?.to_vec())
}

/// R2: block prep + build_graph. Returns CSR-style columnar edges + flat EdgeReads.
#[pyfunction]
#[pyo3(signature = (
    block_start,
    block_end,
    block_is_insertion,
    block_insertion_len,
    block_gap_left,
    block_clip_pos,
    block_clip_reads,
    start1,
    end1,
    start2,
    end2,
    insert,
    mate_len,
    qname_hash,
    mapq1,
    mapq2,
    lib_id,
    flags,
    split_start1,
    split_end1,
    split_reverse1,
    split_start2,
    split_end2,
    split_reverse2,
    split_mapq,
    split_primary_start,
    split_primary_reverse,
    split_qname_hash,
    bp_clip_pos,
    bp_clip_support,
    softclip_pos,
    pmapped,
    q_low,
    q_high,
    min_mapq,
    min_edge_support,
    max_edge_span,
    split_read_leeway,
    pool_libraries,
    skip_edge_ref_gate,
))]
#[allow(clippy::too_many_arguments)]
fn build_graph(
    py: Python<'_>,
    block_start: PyReadonlyArray1<'_, i64>,
    block_end: PyReadonlyArray1<'_, i64>,
    block_is_insertion: PyReadonlyArray1<'_, u8>,
    block_insertion_len: PyReadonlyArray1<'_, i32>,
    block_gap_left: PyReadonlyArray1<'_, u8>,
    block_clip_pos: PyReadonlyArray1<'_, i64>,
    block_clip_reads: PyReadonlyArray1<'_, i32>,
    start1: PyReadonlyArray1<'_, i64>,
    end1: PyReadonlyArray1<'_, i64>,
    start2: PyReadonlyArray1<'_, i64>,
    end2: PyReadonlyArray1<'_, i64>,
    insert: PyReadonlyArray1<'_, i32>,
    mate_len: PyReadonlyArray1<'_, i32>,
    qname_hash: PyReadonlyArray1<'_, i64>,
    mapq1: PyReadonlyArray1<'_, u8>,
    mapq2: PyReadonlyArray1<'_, u8>,
    lib_id: PyReadonlyArray1<'_, u16>,
    flags: PyReadonlyArray1<'_, u8>,
    split_start1: PyReadonlyArray1<'_, i64>,
    split_end1: PyReadonlyArray1<'_, i64>,
    split_reverse1: PyReadonlyArray1<'_, u8>,
    split_start2: PyReadonlyArray1<'_, i64>,
    split_end2: PyReadonlyArray1<'_, i64>,
    split_reverse2: PyReadonlyArray1<'_, u8>,
    split_mapq: PyReadonlyArray1<'_, u8>,
    split_primary_start: PyReadonlyArray1<'_, i64>,
    split_primary_reverse: PyReadonlyArray1<'_, u8>,
    split_qname_hash: PyReadonlyArray1<'_, i64>,
    bp_clip_pos: PyReadonlyArray1<'_, i64>,
    bp_clip_support: PyReadonlyArray1<'_, i32>,
    softclip_pos: PyReadonlyArray1<'_, i64>,
    pmapped: f64,
    q_low: i32,
    q_high: i32,
    min_mapq: u8,
    min_edge_support: i32,
    max_edge_span: i64,
    split_read_leeway: i32,
    pool_libraries: bool,
    skip_edge_ref_gate: bool,
) -> PyResult<PyObject> {
    let bs = arr_i64(block_start)?;
    let be = arr_i64(block_end)?;
    let bi = arr_u8(block_is_insertion)?;
    let bil = arr_i32(block_insertion_len)?;
    let bg = arr_u8(block_gap_left)?;
    let bcp = arr_i64(block_clip_pos)?;
    let bcr = arr_i32(block_clip_reads)?;
    if bs.len() != be.len()
        || bs.len() != bi.len()
        || bs.len() != bil.len()
        || bs.len() != bg.len()
        || bs.len() != bcp.len()
        || bs.len() != bcr.len()
    {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "block array length mismatch",
        ));
    }
    let blocks: Vec<BlockIn> = (0..bs.len())
        .map(|i| BlockIn {
            start: bs[i],
            end: be[i],
            is_insertion: bi[i] != 0,
            insertion_len: bil[i],
            gap_left: bg[i] != 0,
            clip_pos: bcp[i],
            clip_reads: bcr[i],
        })
        .collect();

    let s1 = start1.as_slice()?;
    let e1 = end1.as_slice()?;
    let s2 = start2.as_slice()?;
    let e2 = end2.as_slice()?;
    let ins = insert.as_slice()?;
    let ml = mate_len.as_slice()?;
    let qh = qname_hash.as_slice()?;
    let mq1 = mapq1.as_slice()?;
    let mq2 = mapq2.as_slice()?;
    let lid = lib_id.as_slice()?;
    let fl = flags.as_slice()?;
    let n = s1.len();
    if e1.len() != n
        || s2.len() != n
        || e2.len() != n
        || ins.len() != n
        || ml.len() != n
        || qh.len() != n
        || mq1.len() != n
        || mq2.len() != n
        || lid.len() != n
        || fl.len() != n
    {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "pairs array length mismatch",
        ));
    }

    let ss1 = split_start1.as_slice()?;
    let se1 = split_end1.as_slice()?;
    let sr1 = split_reverse1.as_slice()?;
    let ss2 = split_start2.as_slice()?;
    let se2 = split_end2.as_slice()?;
    let sr2 = split_reverse2.as_slice()?;
    let smq = split_mapq.as_slice()?;
    let sps = split_primary_start.as_slice()?;
    let spr = split_primary_reverse.as_slice()?;
    let sqh = split_qname_hash.as_slice()?;
    let ns = ss1.len();
    if se1.len() != ns
        || sr1.len() != ns
        || ss2.len() != ns
        || se2.len() != ns
        || sr2.len() != ns
        || smq.len() != ns
        || sps.len() != ns
        || spr.len() != ns
        || sqh.len() != ns
    {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "split array length mismatch",
        ));
    }
    let splits: Vec<SplitIn> = (0..ns)
        .map(|i| SplitIn {
            start1: ss1[i],
            end1: se1[i],
            reverse1: sr1[i] != 0,
            start2: ss2[i],
            end2: se2[i],
            reverse2: sr2[i] != 0,
            mapq: smq[i],
            primary_start: sps[i],
            primary_reverse: spr[i] != 0,
            qname_hash: sqh[i],
        })
        .collect();

    let bp_cp = bp_clip_pos.as_slice()?;
    let bp_sup = bp_clip_support.as_slice()?;
    let sc_pos = softclip_pos.as_slice()?;

    let cfg = GraphConfig {
        min_mapq,
        min_edge_support,
        max_edge_span,
        split_read_leeway,
        pool_libraries,
        skip_edge_ref_gate,
        pmapped,
        q_low,
        q_high,
    };

    let cols = py.allow_threads(|| {
        let pairs = PairsIn {
            start1: s1,
            end1: e1,
            start2: s2,
            end2: e2,
            insert: ins,
            mate_len: ml,
            qname_hash: qh,
            mapq1: mq1,
            mapq2: mq2,
            lib_id: lid,
            flags: fl,
        };
        build_graph_columns(blocks, pairs, &splits, &cfg, bp_cp, bp_sup, sc_pos)
    });

    let out = PyDict::new_bound(py);
    out.set_item("edge_u", PyArray1::from_vec_bound(py, cols.edge_u))?;
    out.set_item("edge_v", PyArray1::from_vec_bound(py, cols.edge_v))?;
    out.set_item("edge_support", PyArray1::from_vec_bound(py, cols.edge_support))?;
    out.set_item(
        "edge_reference",
        PyArray1::from_vec_bound(py, cols.edge_reference),
    )?;
    out.set_item("edge_read_ptr", PyArray1::from_vec_bound(py, cols.edge_read_ptr))?;
    out.set_item("read_offset", PyArray1::from_vec_bound(py, cols.read_offset))?;
    out.set_item("read_lib_id", PyArray1::from_vec_bound(py, cols.read_lib_id))?;
    out.set_item("read_pmapped", PyArray1::from_vec_bound(py, cols.read_pmapped))?;
    out.set_item("read_mapq", PyArray1::from_vec_bound(py, cols.read_mapq))?;
    out.set_item("read_kind", PyArray1::from_vec_bound(py, cols.read_kind))?;
    out.set_item("read_flags", PyArray1::from_vec_bound(py, cols.read_flags))?;
    out.set_item("read_mate_len", PyArray1::from_vec_bound(py, cols.read_mate_len))?;
    out.set_item("read_adj1_ptr", PyArray1::from_vec_bound(py, cols.read_adj1_ptr))?;
    out.set_item("read_adj1_len", PyArray1::from_vec_bound(py, cols.read_adj1_len))?;
    out.set_item("read_adj2_ptr", PyArray1::from_vec_bound(py, cols.read_adj2_ptr))?;
    out.set_item("read_adj2_len", PyArray1::from_vec_bound(py, cols.read_adj2_len))?;
    out.set_item("adj_flat", PyArray1::from_vec_bound(py, cols.adj_flat))?;
    out.set_item(
        "block_clip_pos",
        PyArray1::from_vec_bound(py, cols.block_clip_pos),
    )?;
    out.set_item(
        "block_clip_reads",
        PyArray1::from_vec_bound(py, cols.block_clip_reads),
    )?;
    // Export flag constants for Python rebuild.
    out.set_item("KIND_PE", KIND_PE)?;
    out.set_item("KIND_SPLIT", KIND_SPLIT)?;
    out.set_item("ER_HANGING", ER_HANGING)?;
    out.set_item("ER_DISCORDANT", ER_DISCORDANT)?;
    out.set_item("ER_DISTANT", ER_DISTANT)?;
    out.set_item("ER_SUPPORTS_ADJ", ER_SUPPORTS_ADJ)?;
    out.set_item("ER_IN_LIKELIHOOD", ER_IN_LIKELIHOOD)?;
    Ok(out.into())
}

#[pymodule]
fn arcsv2_native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(native_version, m)?)?;
    m.add_function(wrap_pyfunction!(rust_engine_available, m)?)?;
    m.add_function(wrap_pyfunction!(load_bam_pairs, m)?)?;
    m.add_function(wrap_pyfunction!(build_graph, m)?)?;
    Ok(())
}
