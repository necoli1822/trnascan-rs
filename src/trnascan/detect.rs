//! Signal detection and scoring functions for tRNAscan.
//!
//! This module implements the core signal detection algorithms
//! from the Fichant-Burks tRNA scanning method.

use crate::trnascan::matrix::{ConsensusMatrix, MAX_SIGNAL_LEN};
use crate::trnascan::seq_utils::base_to_index;

/// Result of signal detection.
#[derive(Debug, Clone)]
pub struct SignalResult {
    /// Whether a signal was found (enough invariant bases matched)
    pub found: bool,
    /// Weight/frequency values at each position
    pub weight: Vec<f32>,
    /// Number of invariant bases found in the sequence
    pub ninv: i32,
}

/// Check for the presence of a signal at a given position.
///
/// This implements the readsignal function from trnascan.c (lines 1441-1557).
/// It checks if the sequence at the given position matches the consensus
/// matrix well enough to be considered a potential signal.
///
/// # Arguments
/// * `seq` - The DNA sequence (lowercase)
/// * `pos` - Position in sequence to start checking
/// * `matrix` - The consensus matrix to match against
/// * `threshold_inv` - Number of invariant bases allowed NOT to match
///
/// # Returns
/// SignalResult containing whether found, weight values, and invariant count
pub fn read_signal(
    seq: &[u8],
    pos: usize,
    matrix: &ConsensusMatrix,
    threshold_inv: i32,
) -> SignalResult {
    let lsig = matrix.lsig as usize;
    let ktot = matrix.ktot;
    let mut weight = vec![0.0f32; MAX_SIGNAL_LEN];
    let mut ninv = 0i32;

    // Check if we have enough sequence
    if pos + lsig > seq.len() {
        return SignalResult {
            found: false,
            weight,
            ninv: 0,
        };
    }

    // Count invariant bases that match
    // Note: table_inv is 1-indexed in original C code
    if ktot > 0 {
        for k in 1..=(ktot as usize) {
            let inv_pos = matrix.table_inv[k][0] as usize;
            let inv_base = matrix.table_inv[k][1];

            if pos + inv_pos < seq.len() {
                let j = base_to_index(seq[pos + inv_pos]);

                // Match if base matches invariant, or if ambiguous base (j == -1)
                // (original allows ambiguous bases to match anything)
                if j == inv_base || j == -1 {
                    ninv += 1;
                }
            }
        }
    }

    // If not enough invariant bases match, signal not found
    if ninv < ktot - threshold_inv {
        return SignalResult {
            found: false,
            weight,
            ninv,
        };
    }

    // Build weight table
    let mut match1 = true;
    let mut i = 0;

    while match1 && i < lsig {
        if pos + i >= seq.len() {
            match1 = false;
            break;
        }

        let tab_i = base_to_index(seq[pos + i]);

        // Handle ambiguous bases
        if tab_i == -1 {
            // Ambiguous base: assume match (weight = 1)
            weight[i] = 1.0;
        } else {
            weight[i] = matrix.table_cons[i][tab_i as usize];
        }

        // If no invariant bases and weight is 0, discard
        if ktot == 0 && weight[i] == 0.0 {
            match1 = false;
        } else {
            i += 1;
        }
    }

    SignalResult {
        found: match1,
        weight,
        ninv,
    }
}

/// Calculate similarity score for a potential signal.
///
/// This implements the scoring function from trnascan.c (lines 1565-1601).
/// The score is computed as the normalized sum of weights, adjusted for
/// invariant bases.
///
/// # Arguments
/// * `weight` - Weight values from read_signal
/// * `lsig` - Signal length
/// * `max` - Sum of maximum frequencies (maxtot from matrix)
/// * `ktot` - Number of invariant bases in matrix
/// * `threshold` - Minimum score threshold
/// * `ninv` - Number of invariant bases found
///
/// # Returns
/// (passes_threshold, score)
pub fn scoring(
    weight: &[f32],
    lsig: i32,
    max: f32,
    ktot: i32,
    threshold: f32,
    ninv: i32,
) -> (bool, f32) {
    // Sum up weights
    let mut tot: f32 = 0.0;
    for i in 0..(lsig as usize) {
        tot += weight[i];
    }

    // Subtract invariant contributions
    tot -= ninv as f32;
    let adjusted_max = max - ktot as f32;

    // Compute normalized score
    let score = if adjusted_max > 0.0 {
        tot / adjusted_max
    } else {
        0.0
    };

    (score >= threshold, score)
}

/// Count base pairings in a stem structure.
///
/// This implements the basepairing function from trnascan.c (lines 1608-1654).
/// It counts Watson-Crick base pairs (A-T, G-C) and wobble pairs (G-T).
///
/// # Arguments
/// * `seq` - The DNA sequence (starting at first base of stem)
/// * `npair` - Number of base pairs to check
/// * `lpair` - Distance between first position of 5' strand and last position of 3' strand
///
/// # Returns
/// Number of valid base pairs found
pub fn basepairing(seq: &[u8], npair: i32, lpair: i32) -> i32 {
    let mut ncomp = 0;

    for n in 0..(npair as usize) {
        let pos1 = n;
        let pos2 = lpair as usize - n;

        if pos2 >= seq.len() {
            continue;
        }

        let b1 = seq.get(pos1).copied().unwrap_or(b'n');
        let b2 = seq.get(pos2).copied().unwrap_or(b'n');

        // Handle N's - if one is N and other is not N, count as match
        // (conservative calling of N's as base pairing matches)
        if b2 == b'n' && b1 != b'n' {
            ncomp += 1;
            continue;
        }

        // Standard base pairing rules
        match b1 {
            b'a' => {
                if b2 == b't' {
                    ncomp += 1;
                }
            }
            b'c' => {
                if b2 == b'g' {
                    ncomp += 1;
                }
            }
            b'g' | b'r' => {
                // G pairs with C or T (wobble)
                if b2 == b'c' || b2 == b't' {
                    ncomp += 1;
                }
            }
            b't' | b'y' => {
                // T pairs with A or G (wobble)
                if b2 == b'a' || b2 == b'g' {
                    ncomp += 1;
                }
            }
            b'n' => {
                // N matches anything except N
                if b2 != b'n' {
                    ncomp += 1;
                }
            }
            _ => {}
        }
    }

    ncomp
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::trnascan::matrix::ConsensusMatrix;

    #[test]
    fn test_read_signal_basic() {
        let matrix = ConsensusMatrix::d_signal();
        // D signal expects: T at pos 0, G at pos 2, A at pos 6
        let seq = b"tggcccag"; // Should match D signal invariants

        let result = read_signal(seq, 0, &matrix, 1);
        assert!(result.found);
        assert_eq!(result.ninv, 3); // All 3 invariants should match
    }

    #[test]
    fn test_read_signal_partial_match() {
        let matrix = ConsensusMatrix::d_signal();
        // Missing one invariant
        let seq = b"aggcccag"; // A instead of T at position 0

        let result = read_signal(seq, 0, &matrix, 1);
        assert!(result.found); // Should still pass with threshold_inv = 1
        assert_eq!(result.ninv, 2);
    }

    #[test]
    fn test_read_signal_too_short() {
        let matrix = ConsensusMatrix::d_signal();
        let seq = b"tggc"; // Too short

        let result = read_signal(seq, 0, &matrix, 1);
        assert!(!result.found);
    }

    #[test]
    fn test_scoring() {
        let weight = [0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 1.0, 0.5];
        let (passes, score) = scoring(&weight, 8, 8.0, 2, 0.4, 2);

        // Total = 5.0, ninv = 2, so adjusted = 3.0
        // max = 8.0, ktot = 2, so adjusted_max = 6.0
        // score = 3.0 / 6.0 = 0.5
        assert!(passes);
        assert!((score - 0.5).abs() < 0.01);
    }

    #[test]
    fn test_basepairing_perfect() {
        // Perfect stem: GCGCG...CGCGC
        let seq = b"gcgcgnnnnncgcgc";
        // npair=5, lpair=14 (positions 0-14)
        let ncomp = basepairing(seq, 5, 14);
        assert_eq!(ncomp, 5);
    }

    #[test]
    fn test_basepairing_wobble() {
        // G-T wobble pair
        let seq = b"gnnnnnnnnnnnnnt";
        let ncomp = basepairing(seq, 1, 14);
        assert_eq!(ncomp, 1);
    }

    #[test]
    fn test_basepairing_mismatch() {
        // A-A mismatch
        let seq = b"annnnnnnnnnnnna";
        let ncomp = basepairing(seq, 1, 14);
        assert_eq!(ncomp, 0);
    }

    #[test]
    fn test_basepairing_with_n() {
        // N matches anything
        let seq = b"nnnnnnnnnnnnnna";
        let ncomp = basepairing(seq, 1, 14);
        assert_eq!(ncomp, 1); // N-A should count as match
    }
}
