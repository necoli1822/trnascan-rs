//! Viterbi algorithm module for tRNAscan-SE
//!
//! This module implements the core Viterbi dynamic programming algorithm
//! for aligning sequences to covariance models. It is a direct port of
//! the original C code from viterbi.c and model.c.
//!
//! # Main Components
//!
//! - `model.rs`: RearrangeCM - converts CM to integer log-odds states
//! - `viterbi.rs`: ViterbiAlign, matrix fill, and traceback
//!
//! # Algorithm Overview
//!
//! The Viterbi algorithm uses two matrices:
//! - `amx[j][diff][y]`: Main score matrix
//! - `bmx[y][j][diff]`: BEGIN state score matrix (for bifurcation handling)
//!
//! Where:
//! - j = sequence position (0..N)
//! - diff = j - i + 1 (diagonal difference)
//! - y = state index (0..statenum-1)

mod dbscan;
mod model;
mod small;
mod viterbi;

pub use dbscan::*;
pub use model::*;
pub use small::*;
pub use viterbi::*;
