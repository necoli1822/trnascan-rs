//! Consensus matrix handling for tRNAscan signal detection.
//!
//! This module implements the lectval function from the C code,
//! which reads and processes consensus matrices for signal detection.

use crate::trnascan::signals::{D_SIGNAL, TPC_SIGNAL};

/// Maximum signal length supported (positions in consensus matrix).
pub const MAX_SIGNAL_LEN: usize = 30;

/// Consensus matrix for signal detection.
///
/// This structure holds the frequency matrix and invariant base
/// information needed for signal detection in tRNA scanning.
#[derive(Debug, Clone)]
pub struct ConsensusMatrix {
    /// Frequency of each base (A=0, C=1, G=2, T=3) at each position.
    /// table_cons[position][base]
    pub table_cons: [[f32; 4]; MAX_SIGNAL_LEN],
    /// Position and nature of invariant bases.
    /// table_inv[k][0] = position, table_inv[k][1] = base code
    /// Base codes: A=0, C=1, G=2, T=3
    pub table_inv: [[i32; 2]; MAX_SIGNAL_LEN],
    /// Signal length (number of positions)
    pub lsig: i32,
    /// Number of invariant bases (positions with frequency 1.0)
    pub ktot: i32,
    /// Sum of maximum frequencies at each position
    pub maxtot: f32,
}

impl Default for ConsensusMatrix {
    fn default() -> Self {
        Self {
            table_cons: [[0.0; 4]; MAX_SIGNAL_LEN],
            table_inv: [[0; 2]; MAX_SIGNAL_LEN],
            lsig: 0,
            ktot: 0,
            maxtot: 0.0,
        }
    }
}

impl ConsensusMatrix {
    /// Create a new empty consensus matrix.
    pub fn new() -> Self {
        Self::default()
    }

    /// Load the T-Psi-C signal consensus matrix.
    pub fn tpc_signal() -> Self {
        load_consensus_matrix(TPC_SIGNAL)
    }

    /// Load the D signal consensus matrix.
    pub fn d_signal() -> Self {
        load_consensus_matrix(D_SIGNAL)
    }
}

/// Load a consensus matrix from string content.
///
/// This implements the lectval function from trnascan.c (lines 1138-1195).
/// The format is 4 floats per line (A, C, G, T frequencies).
///
/// # Arguments
/// * `content` - String containing the consensus matrix data
///
/// # Returns
/// A populated ConsensusMatrix structure
pub fn load_consensus_matrix(content: &str) -> ConsensusMatrix {
    let mut matrix = ConsensusMatrix::new();
    let mut i = 0;

    // Parse each line of frequencies
    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        let values: Vec<f32> = line
            .split_whitespace()
            .filter_map(|s| s.parse::<f32>().ok())
            .collect();

        if values.len() == 4 && i < MAX_SIGNAL_LEN {
            for (j, &val) in values.iter().enumerate() {
                matrix.table_cons[i][j] = val;
            }
            i += 1;
        }
    }

    matrix.lsig = i as i32;

    // Find invariant bases (positions where one base has frequency 1.0)
    // Note: ktot is 1-indexed in original C code
    let mut k = 0;
    for pos in 0..matrix.lsig as usize {
        for base in 0..4 {
            if (matrix.table_cons[pos][base] - 1.0).abs() < 0.0001 {
                k += 1;
                if k < MAX_SIGNAL_LEN as i32 {
                    matrix.table_inv[k as usize][0] = pos as i32;
                    matrix.table_inv[k as usize][1] = base as i32;
                }
            }
        }
    }
    matrix.ktot = k;

    // Calculate sum of maximum frequencies
    matrix.maxtot = 0.0;
    for pos in 0..matrix.lsig as usize {
        let max_freq = matrix.table_cons[pos]
            .iter()
            .cloned()
            .fold(f32::NEG_INFINITY, f32::max);
        matrix.maxtot += max_freq;
    }

    matrix
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_load_tpc_signal() {
        let matrix = ConsensusMatrix::tpc_signal();
        assert_eq!(matrix.lsig, 15);
        // TPC signal has 4 invariant positions: G at pos 5, T at pos 7, C at pos 8, C at pos 13
        assert_eq!(matrix.ktot, 4);
        // Check first position frequencies
        assert!((matrix.table_cons[0][0] - 0.0124).abs() < 0.001);
        assert!((matrix.table_cons[0][1] - 0.8100).abs() < 0.001);
    }

    #[test]
    fn test_load_d_signal() {
        let matrix = ConsensusMatrix::d_signal();
        assert_eq!(matrix.lsig, 8);
        // D signal has 3 invariant positions: T at pos 0, G at pos 2, A at pos 6
        assert_eq!(matrix.ktot, 3);
    }

    #[test]
    fn test_invariant_positions_tpc() {
        let matrix = ConsensusMatrix::tpc_signal();
        // Position 5 (0-indexed) should be invariant G
        assert!((matrix.table_cons[5][2] - 1.0).abs() < 0.001);
        // Position 7 should be invariant T
        assert!((matrix.table_cons[7][3] - 1.0).abs() < 0.001);
        // Position 8 should be invariant C
        assert!((matrix.table_cons[8][1] - 1.0).abs() < 0.001);
        // Position 13 should be invariant C
        assert!((matrix.table_cons[13][1] - 1.0).abs() < 0.001);
    }

    #[test]
    fn test_invariant_positions_d() {
        let matrix = ConsensusMatrix::d_signal();
        // Position 0 should be invariant T
        assert!((matrix.table_cons[0][3] - 1.0).abs() < 0.001);
        // Position 2 should be invariant G
        assert!((matrix.table_cons[2][2] - 1.0).abs() < 0.001);
        // Position 6 should be invariant A
        assert!((matrix.table_cons[6][0] - 1.0).abs() < 0.001);
    }

    #[test]
    fn test_maxtot_calculation() {
        let matrix = ConsensusMatrix::tpc_signal();
        // maxtot should be close to 15 since we have 15 positions
        // and invariant positions contribute 1.0 each
        assert!(matrix.maxtot > 10.0);
        assert!(matrix.maxtot <= 15.0);
    }
}
