//! ARC-SV2 numerical kernels.
//!
//! These match the Python reference engine in `python/arcsv2/engine`. Paper
//! defaults live in the Python config; this crate does not hard-code sample
//! paths and does not use `unsafe`.

pub mod affected;
pub mod drp;
pub mod library;
pub mod names;
pub mod paths;

pub const VERSION: &str = "0.1.0";
