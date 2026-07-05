//! Maximum Likelihood Covariance Model Construction
//!
//! This module implements the MaxModelMaker algorithm from maxmodelmaker.c
//! for constructing covariance models from multiple sequence alignments.
//!
//! # Overview
//!
//! MaxModelMaker uses a dynamic programming algorithm to find the optimal
//! tree structure that generates a given alignment with maximum likelihood.
//! The algorithm:
//!
//! 1. Transposes the alignment for efficient column access
//! 2. Pre-calculates singlet emission scores for each column
//! 3. Fills a scoring matrix via DP recursion
//! 4. Traces back to construct the consensus tree structure
//! 5. Builds a covariance model from the trace
//!
//! # Module Structure
//!
//! - `types` - Core types (MaxMxCell, node constants)
//! - `matrix` - Scoring matrix allocation and initialization
//! - `prior` - Prior probability distributions for regularization
//! - `emission` - Singlet and pair emission score calculations
//! - `transition` - State transition table construction
//! - `recursion` - Main DP recursion
//! - `traceback` - Traceback to construct consensus tree
//!
//! # Example
//!
//! ```ignore
//! use trnascan_rs::maxmodel::{maxmodelmaker, Prior};
//!
//! let sequences = vec![b"ACGU".as_slice(), b"ACGU".as_slice()];
//! let weights = vec![1.0, 1.0];
//! let prior = Prior::new();
//!
//! let result = maxmodelmaker(&sequences, &weights, 0.5, &prior);
//! ```

pub mod emission;
pub mod matrix;
pub mod prior;
pub mod recursion;
pub mod traceback;
pub mod transition;
pub mod types;

// Re-export main types and functions
pub use emission::{is_gap, pair_emissioncost, singlet_emissions, symbol_index, transpose_alignment};
pub use matrix::MaxMx;
pub use prior::Prior;
pub use recursion::recurse_maxmx;
pub use traceback::{trace_maxmx, MasterTrace};
pub use types::{MaxMxCell, MAXINSERT, MAXMX_NODETYPES, maxmx_node};

use crate::types::constants::*;

/// Error type for MaxModelMaker operations
#[derive(Debug, Clone)]
pub enum MaxModelError {
    /// Empty alignment provided
    EmptyAlignment,
    /// Alignment sequences have different lengths
    LengthMismatch,
    /// Number of weights doesn't match number of sequences
    WeightMismatch,
    /// Model construction failed
    ConstructionFailed(String),
}

impl std::fmt::Display for MaxModelError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MaxModelError::EmptyAlignment => write!(f, "Empty alignment provided"),
            MaxModelError::LengthMismatch => {
                write!(f, "Alignment sequences have different lengths")
            }
            MaxModelError::WeightMismatch => {
                write!(f, "Number of weights doesn't match number of sequences")
            }
            MaxModelError::ConstructionFailed(msg) => {
                write!(f, "Model construction failed: {}", msg)
            }
        }
    }
}

impl std::error::Error for MaxModelError {}

/// Result of MaxModelMaker
#[derive(Debug)]
pub struct MaxModelResult {
    /// The consensus master trace tree
    pub trace: MasterTrace,

    /// Information content of the alignment (bits)
    pub info_content: f64,

    /// Number of nodes in the consensus tree
    pub node_count: usize,
}

/// Create a maximally likely model structure from a multiple sequence alignment
///
/// This is the main entry point for the MaxModelMaker algorithm.
///
/// From maxmodelmaker.c Maxmodelmaker (lines 142-256)
///
/// # Arguments
///
/// * `aseqs` - Flushed sequence alignment; each sequence is 0..alen-1
/// * `weights` - Weights assigned to each sequence (usually 1.0)
/// * `gapthresh` - Heuristic: fractional occupancy threshold for MAT assignment
/// * `prior` - Prior probability distributions
///
/// # Returns
///
/// `Ok(MaxModelResult)` containing the consensus trace and info content,
/// or `Err(MaxModelError)` on failure.
///
/// # Example
///
/// ```ignore
/// let seqs = vec![b"ACGU".as_slice(), b"ACGU".as_slice()];
/// let weights = vec![1.0, 1.0];
/// let prior = Prior::new();
/// let result = maxmodelmaker(&seqs, &weights, 0.5, &prior)?;
/// println!("Info content: {} bits", result.info_content);
/// ```
pub fn maxmodelmaker(
    aseqs: &[&[u8]],
    weights: &[f32],
    gapthresh: f64,
    prior: &Prior,
) -> Result<MaxModelResult, MaxModelError> {
    // Validate inputs
    if aseqs.is_empty() {
        return Err(MaxModelError::EmptyAlignment);
    }

    let nseq = aseqs.len();
    let alen = aseqs[0].len();

    if alen == 0 {
        return Err(MaxModelError::EmptyAlignment);
    }

    // Check all sequences have same length
    for seq in aseqs.iter() {
        if seq.len() != alen {
            return Err(MaxModelError::LengthMismatch);
        }
    }

    // Check weights match sequences
    if weights.len() != nseq {
        return Err(MaxModelError::WeightMismatch);
    }

    // Step 1: Transpose alignment to [position][seq] indexing
    let aseqs_t = transpose_alignment(aseqs, alen, is_gap, symbol_index);

    // Step 2: Pre-calculate singlet emission scores
    let (mscore, gapcount) = singlet_emissions(&aseqs_t, weights, prior);

    // Step 3: Allocate and initialize scoring matrix
    let mut mmx = MaxMx::new(alen);
    mmx.initialize(nseq, prior, &mscore, &gapcount);

    // Step 4: Fill matrix via recursion
    recurse_maxmx(
        &aseqs_t,
        weights,
        prior,
        &mscore,
        &gapcount,
        gapthresh,
        &mut mmx,
    );

    // Step 5: Calculate information content
    // Score is stored in mmx[alen][0].sc[MATP] (using MATP slot for ROOT)
    let info_content = mmx.get(alen, 0).sc[maxmx_node::MATP] as f64 / INTPRECISION;

    // Step 6: Traceback to construct consensus tree
    let mut trace = trace_maxmx(&mmx);

    // Step 7: Number the nodes in the tree
    let node_count = trace.count_nodes();
    trace.number_nodes(0);

    Ok(MaxModelResult {
        trace,
        info_content,
        node_count,
    })
}

/// Convenience function to create weights array of 1.0
pub fn uniform_weights(nseq: usize) -> Vec<f32> {
    vec![1.0; nseq]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_maxmodelmaker_simple() {
        let seq1: &[u8] = b"ACGU";
        let seq2: &[u8] = b"ACGU";
        let aseqs: Vec<&[u8]> = vec![seq1, seq2];
        let weights = uniform_weights(2);
        let prior = Prior::new();

        let result = maxmodelmaker(&aseqs, &weights, 0.5, &prior);
        assert!(result.is_ok());

        let result = result.unwrap();
        assert!(result.node_count > 0);
        assert_eq!(result.trace.node_type, ROOT_NODE);
    }

    #[test]
    fn test_maxmodelmaker_with_gaps() {
        let seq1: &[u8] = b"AC-U";
        let seq2: &[u8] = b"ACGU";
        let aseqs: Vec<&[u8]> = vec![seq1, seq2];
        let weights = uniform_weights(2);
        let prior = Prior::new();

        let result = maxmodelmaker(&aseqs, &weights, 0.5, &prior);
        assert!(result.is_ok());
    }

    #[test]
    fn test_maxmodelmaker_trna_prior() {
        let seq1: &[u8] = b"ACGUACGU";
        let seq2: &[u8] = b"ACGUACGU";
        let aseqs: Vec<&[u8]> = vec![seq1, seq2];
        let weights = uniform_weights(2);
        let prior = Prior::new_trna();

        let result = maxmodelmaker(&aseqs, &weights, 0.5, &prior);
        assert!(result.is_ok());
    }

    #[test]
    fn test_maxmodelmaker_empty_alignment() {
        let aseqs: Vec<&[u8]> = vec![];
        let weights: Vec<f32> = vec![];
        let prior = Prior::new();

        let result = maxmodelmaker(&aseqs, &weights, 0.5, &prior);
        assert!(matches!(result, Err(MaxModelError::EmptyAlignment)));
    }

    #[test]
    fn test_maxmodelmaker_length_mismatch() {
        let seq1: &[u8] = b"ACGU";
        let seq2: &[u8] = b"ACG"; // Different length
        let aseqs: Vec<&[u8]> = vec![seq1, seq2];
        let weights = uniform_weights(2);
        let prior = Prior::new();

        let result = maxmodelmaker(&aseqs, &weights, 0.5, &prior);
        assert!(matches!(result, Err(MaxModelError::LengthMismatch)));
    }

    #[test]
    fn test_maxmodelmaker_weight_mismatch() {
        let seq1: &[u8] = b"ACGU";
        let seq2: &[u8] = b"ACGU";
        let aseqs: Vec<&[u8]> = vec![seq1, seq2];
        let weights = vec![1.0f32]; // Only one weight
        let prior = Prior::new();

        let result = maxmodelmaker(&aseqs, &weights, 0.5, &prior);
        assert!(matches!(result, Err(MaxModelError::WeightMismatch)));
    }

    #[test]
    fn test_uniform_weights() {
        let weights = uniform_weights(5);
        assert_eq!(weights.len(), 5);
        for w in weights {
            assert_eq!(w, 1.0);
        }
    }

    #[test]
    fn test_maxmodelmaker_longer_alignment() {
        // Test with a more realistic tRNA-like alignment
        let seq1: &[u8] = b"GCGGAUUUAGCUCAGDDGGGAGAGCGCCAGACUGAAYAUCUGGAGGUCCUGUGTPCGAUCCACAGAAUUCGCACCA";
        let seq2: &[u8] = b"GCGGAUUUAGCUCAGDDGGGAGAGCGCCAGACUGAAYAUCUGGAGGUCCUGUGTPCGAUCCACAGAAUUCGCACCA";
        let aseqs: Vec<&[u8]> = vec![seq1, seq2];
        let weights = uniform_weights(2);
        let prior = Prior::new_trna();

        let result = maxmodelmaker(&aseqs, &weights, 0.5, &prior);
        assert!(result.is_ok());

        let result = result.unwrap();
        assert!(result.node_count > 1);
    }
}
