//! Secondary structure representation for tRNA
//!
//! This module implements secondary structure functions from konings.c:
//! - Trace2KHS: Convert traceback to Konings/Hogeweg structure string
//! - KHS2ct: Convert structure string to connect table
//! - IsRNAComplement: Check if bases can pair
//! - Align2kh: Convert alignment to structure string

mod konings;

pub use konings::{
    align2kh,
    is_rna_complement,
    khs2ct,
    trace2khs,
};
