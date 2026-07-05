//! Embedded consensus signal data for tRNAscan.
//!
//! This module contains the T-Psi-C and D signal consensus matrices
//! embedded as constants for use without external file dependencies.

/// Embedded TPCsignal data (from original/lib/models/TPCsignal).
///
/// Each line contains 4 floats representing the frequency of A, C, G, T
/// at each position of the T-Psi-C signal (15 positions total).
/// Invariant bases have frequency 1.0 for one nucleotide.
pub const TPC_SIGNAL: &str = r#"0.0124 0.8100 0.0083 0.1694
0.1618 0.3320 0.4772 0.0290
0.0992 0.4008 0.2934 0.2066
0.1818 0.1818 0.5083 0.1281
0.1612 0.0289 0.7768 0.0331
0.0000 0.0000 1.0000 0.0000
0.0581 0.0000 0.0000 0.9419
0.0000 0.0000 0.0000 1.0000
0.0000 1.0000 0.0000 0.0000
0.1818 0.0000 0.8182 0.0000
0.9876 0.0000 0.0083 0.0041
0.4959 0.0620 0.1694 0.2727
0.0620 0.1859 0.0083 0.7438
0.0000 1.0000 0.0000 0.0000
0.0289 0.7686 0.0330 0.1694
"#;

/// Embedded Dsignal data (from original/lib/models/Dsignal).
///
/// Each line contains 4 floats representing the frequency of A, C, G, T
/// at each position of the D signal (8 positions total).
/// Invariant bases have frequency 1.0 for one nucleotide.
pub const D_SIGNAL: &str = r#"0.0000 0.0000 0.0000 1.0000
0.5100 0.0400 0.4200 0.0200
0.0000 0.0000 1.0000 0.0000
0.0400 0.6400 0.0500 0.2800
0.0500 0.2500 0.3000 0.4000
0.0900 0.4900 0.1300 0.2900
1.0000 0.0000 0.0000 0.0000
0.1800 0.0100 0.7500 0.0600
"#;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tpc_signal_lines() {
        let lines: Vec<&str> = TPC_SIGNAL.lines().filter(|l| !l.is_empty()).collect();
        assert_eq!(lines.len(), 15);
    }

    #[test]
    fn test_d_signal_lines() {
        let lines: Vec<&str> = D_SIGNAL.lines().filter(|l| !l.is_empty()).collect();
        assert_eq!(lines.len(), 8);
    }

    #[test]
    fn test_tpc_signal_format() {
        for line in TPC_SIGNAL.lines().filter(|l| !l.is_empty()) {
            let values: Vec<f32> = line
                .split_whitespace()
                .map(|s| s.parse::<f32>().unwrap())
                .collect();
            assert_eq!(values.len(), 4);
            // Sum should be approximately 1.0
            let sum: f32 = values.iter().sum();
            assert!((sum - 1.0).abs() < 0.01, "Sum was {}", sum);
        }
    }

    #[test]
    fn test_d_signal_format() {
        for line in D_SIGNAL.lines().filter(|l| !l.is_empty()) {
            let values: Vec<f32> = line
                .split_whitespace()
                .map(|s| s.parse::<f32>().unwrap())
                .collect();
            assert_eq!(values.len(), 4);
            let sum: f32 = values.iter().sum();
            assert!((sum - 1.0).abs() < 0.01, "Sum was {}", sum);
        }
    }
}
