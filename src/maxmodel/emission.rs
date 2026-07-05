//! Emission score calculations for MaxModelMaker
//!
//! This module implements singlet_emissions and pair_emissioncost functions
//! from maxmodelmaker.c lines 989-1154.

use crate::maxmodel::prior::{Prior, U_MATL_ST, U_MATR_ST};
use crate::maxmodel::types::*;
use crate::types::constants::*;

/// Calculate dot product of count and probability vectors
///
/// Returns the log-likelihood score as an integer.
///
/// From maxmodelmaker.c dot_score (lines 1733-1744)
pub fn dot_score(cvec: &[f64], pvec: &[f64]) -> i32 {
    let mut score = 0.0;
    for (c, p) in cvec.iter().zip(pvec.iter()) {
        if *p > 0.0 {
            score += c * p.ln();
        }
    }
    (INTPRECISION * score) as i32
}

/// Calculate dot product for 2D matrices (flattened)
pub fn dot_score_2d(cmat: &[[f64; STATETYPES]; STATETYPES], pmat: &[[f64; STATETYPES]; STATETYPES]) -> i32 {
    let mut score = 0.0;
    for i in 0..STATETYPES {
        for j in 0..STATETYPES {
            if pmat[i][j] > 0.0 {
                score += cmat[i][j] * pmat[i][j].ln();
            }
        }
    }
    (INTPRECISION * score) as i32
}

/// Singlet emission scores for each alignment column
///
/// Pre-calculate expected singlet emission scores for each column and
/// count gaps in each column.
///
/// From maxmodelmaker.c singlet_emissions (lines 1008-1071)
///
/// # Arguments
/// * `aseqs_t` - Transposed alignment [1..alen][0..nseq-1], -1 for gaps, 0-3 for bases
/// * `weights` - Weights on sequences (usually 1.0 for each)
/// * `prior` - Prior probability distributions
///
/// # Returns
/// * `mscore` - Array of singlet emission scores for each column [0..alen]
/// * `gapcount` - Weighted counts of gaps in each column [0..alen+1]
pub fn singlet_emissions(
    aseqs_t: &[Vec<i8>],
    weights: &[f32],
    prior: &Prior,
) -> (Vec<i32>, Vec<f64>) {
    let alen = aseqs_t.len() - 2; // aseqs_t has guard columns at 0 and alen+1
    let nseq = weights.len();

    // Allocate output arrays
    let mut mscore = vec![0i32; alen + 2];
    let mut gapcount = vec![0.0f64; alen + 2];

    // Count symbol occurrences in each column
    let mut emcounts = vec![[0.0f64; ALPHASIZE]; alen + 1];

    for i in 1..=alen {
        gapcount[i] = 0.0;
        for sym in 0..ALPHASIZE {
            emcounts[i][sym] = 0.0;
        }

        for (idx, weight) in weights.iter().enumerate().take(nseq) {
            let sym = aseqs_t[i][idx];
            if sym >= 0 {
                emcounts[i][sym as usize] += *weight as f64;
            } else {
                gapcount[i] += *weight as f64;
            }
        }
    }
    gapcount[0] = 0.0;
    gapcount[alen + 1] = 0.0;

    // For each column, create emission probability vector and calculate score
    for i in 1..=alen {
        let mut emvec = [0.0f64; ALPHASIZE];
        for sym in 0..ALPHASIZE {
            emvec[sym] = emcounts[i][sym];
        }

        // Apply prior regularization
        prior.probify_singlet_emission(&mut emvec, U_MATL_ST);

        // Score = dot product of counts and log probabilities
        mscore[i] = dot_score(&emcounts[i], &emvec);

        // Add log(ALPHASIZE) for each non-gap symbol (log-odds conversion)
        let non_gap_count = nseq as f64 - gapcount[i];
        mscore[i] += (non_gap_count * INTPRECISION * (ALPHASIZE as f64).ln()) as i32;
    }

    (mscore, gapcount)
}

/// Calculate pair emission cost for a MATP node assignment
///
/// Count emission statistics for columns i and j, and assign a score.
///
/// From maxmodelmaker.c pair_emissioncost (lines 1089-1154)
///
/// # Arguments
/// * `coli` - Column i from transposed alignment [0..nseq-1]
/// * `colj` - Column j from transposed alignment [0..nseq-1]
/// * `weights` - Weights on sequences
/// * `prior` - Prior probability distributions
///
/// # Returns
/// Integer emission cost for the (i,j) column pair
pub fn pair_emissioncost(
    coli: &[i8],
    colj: &[i8],
    weights: &[f32],
    prior: &Prior,
) -> i32 {
    let nseq = weights.len();

    // Count arrays
    let mut matp_count = [[0.0f64; ALPHASIZE]; ALPHASIZE];
    let mut matl_count = [0.0f64; ALPHASIZE];
    let mut matr_count = [0.0f64; ALPHASIZE];

    // Count pairs and singlets
    for idx in 0..nseq {
        let symi = coli[idx];
        let symj = colj[idx];

        if symi == -1 {
            // Gap at i
            if symj != -1 {
                // Symbol at j only -> MATL (right emitter relative to gap)
                matl_count[symj as usize] += weights[idx] as f64;
            }
            // Both gaps -> no emission
        } else if symj == -1 {
            // Gap at j, symbol at i -> MATR
            matr_count[symi as usize] += weights[idx] as f64;
        } else {
            // Symbols at both -> MATP
            matp_count[symi as usize][symj as usize] += weights[idx] as f64;
        }
    }

    // Create probability matrices
    let mut matl_emit = [0.0f64; ALPHASIZE];
    let mut matr_emit = [0.0f64; ALPHASIZE];
    let mut matp_emit = [[0.0f64; ALPHASIZE]; ALPHASIZE];

    copy_singlet(&mut matl_emit, &matl_count);
    copy_singlet(&mut matr_emit, &matr_count);
    copy_pairwise(&mut matp_emit, &matp_count);

    prior.probify_singlet_emission(&mut matl_emit, U_MATL_ST);
    prior.probify_singlet_emission(&mut matr_emit, U_MATR_ST);
    prior.probify_pair_emission(&mut matp_emit);

    // Convert probabilities to log-odds
    for symi in 0..ALPHASIZE {
        matl_emit[symi] *= ALPHASIZE as f64;
        matr_emit[symi] *= ALPHASIZE as f64;
        for symj in 0..ALPHASIZE {
            matp_emit[symi][symj] *= (ALPHASIZE * ALPHASIZE) as f64;
        }
    }

    // Score is sum of dot products
    let mut sc = dot_score(&matl_count, &matl_emit);
    sc += dot_score(&matr_count, &matr_emit);

    // Flatten pair matrices for dot product
    let matp_count_flat: Vec<f64> = matp_count.iter().flatten().copied().collect();
    let matp_emit_flat: Vec<f64> = matp_emit.iter().flatten().copied().collect();
    sc += dot_score(&matp_count_flat, &matp_emit_flat);

    sc
}

/// Transpose alignment from [seq][pos] to [pos][seq] indexing
///
/// Also shifts to 1-based column indexing and converts symbols to indices.
/// Gaps are represented as -1.
///
/// From maxmodelmaker.c transpose_alignment (lines 953-983)
///
/// # Arguments
/// * `aseqs` - Original alignment [0..nseq-1][0..alen-1]
/// * `is_gap` - Function to check if a character is a gap
/// * `symbol_index` - Function to convert character to index (0-3)
///
/// # Returns
/// Transposed alignment [0..alen+1][0..nseq-1] with guard columns
pub fn transpose_alignment<F, G>(
    aseqs: &[&[u8]],
    alen: usize,
    is_gap: F,
    symbol_index: G,
) -> Vec<Vec<i8>>
where
    F: Fn(u8) -> bool,
    G: Fn(u8) -> i8,
{
    let nseq = aseqs.len();

    // Allocate [0..alen+1][0..nseq-1]
    let mut aseqs_t = Vec::with_capacity(alen + 2);
    for _ in 0..=(alen + 1) {
        aseqs_t.push(vec![-1i8; nseq]);
    }

    // Guard columns 0 and alen+1 are already -1

    // Fill in the actual alignment
    for (seqidx, seq) in aseqs.iter().enumerate() {
        for acol in 0..alen {
            let ch = seq[acol];
            aseqs_t[acol + 1][seqidx] = if is_gap(ch) { -1 } else { symbol_index(ch) };
        }
    }

    aseqs_t
}

/// Default gap check function
pub fn is_gap(c: u8) -> bool {
    matches!(c, b'-' | b'.' | b'_' | b' ')
}

/// Default symbol to index function (A=0, C=1, G=2, T/U=3)
pub fn symbol_index(c: u8) -> i8 {
    match c.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' | b'U' => 3,
        _ => -1, // Unknown -> treat as gap
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dot_score() {
        let counts = [1.0, 2.0, 3.0, 4.0];
        let probs = [0.25, 0.25, 0.25, 0.25];

        let score = dot_score(&counts, &probs);
        // All same probability -> score depends on counts * log(0.25)
        let expected = (10.0 * 0.25f64.ln() * INTPRECISION) as i32;
        assert_eq!(score, expected);
    }

    #[test]
    fn test_is_gap() {
        assert!(is_gap(b'-'));
        assert!(is_gap(b'.'));
        assert!(!is_gap(b'A'));
        assert!(!is_gap(b'G'));
    }

    #[test]
    fn test_symbol_index() {
        assert_eq!(symbol_index(b'A'), 0);
        assert_eq!(symbol_index(b'a'), 0);
        assert_eq!(symbol_index(b'C'), 1);
        assert_eq!(symbol_index(b'G'), 2);
        assert_eq!(symbol_index(b'T'), 3);
        assert_eq!(symbol_index(b'U'), 3);
        assert_eq!(symbol_index(b'N'), -1);
    }

    #[test]
    fn test_transpose_alignment() {
        let seq1: &[u8] = b"ACGT";
        let seq2: &[u8] = b"A-GT";
        let aseqs: Vec<&[u8]> = vec![seq1, seq2];

        let aseqs_t = transpose_alignment(&aseqs, 4, is_gap, symbol_index);

        // Should have 6 columns (0..5 for alen=4)
        assert_eq!(aseqs_t.len(), 6);

        // Check guard columns
        assert_eq!(aseqs_t[0][0], -1);
        assert_eq!(aseqs_t[0][1], -1);
        assert_eq!(aseqs_t[5][0], -1);
        assert_eq!(aseqs_t[5][1], -1);

        // Check column 1 (A's)
        assert_eq!(aseqs_t[1][0], 0); // A
        assert_eq!(aseqs_t[1][1], 0); // A

        // Check column 2 (C and gap)
        assert_eq!(aseqs_t[2][0], 1);  // C
        assert_eq!(aseqs_t[2][1], -1); // gap

        // Check column 3 (G's)
        assert_eq!(aseqs_t[3][0], 2); // G
        assert_eq!(aseqs_t[3][1], 2); // G

        // Check column 4 (T's)
        assert_eq!(aseqs_t[4][0], 3); // T
        assert_eq!(aseqs_t[4][1], 3); // T
    }

    #[test]
    fn test_singlet_emissions() {
        let prior = Prior::new();

        // Create simple alignment
        let seq1: &[u8] = b"AAAA";
        let seq2: &[u8] = b"AAAA";
        let aseqs: Vec<&[u8]> = vec![seq1, seq2];
        let weights = vec![1.0f32, 1.0f32];

        let aseqs_t = transpose_alignment(&aseqs, 4, is_gap, symbol_index);
        let (mscore, gapcount) = singlet_emissions(&aseqs_t, &weights, &prior);

        // No gaps expected
        for i in 1..=4 {
            assert_eq!(gapcount[i], 0.0);
        }

        // All A's -> consistent scores
        assert_eq!(mscore[1], mscore[2]);
        assert_eq!(mscore[2], mscore[3]);
    }

    #[test]
    fn test_pair_emissioncost() {
        let prior = Prior::new();

        // Column with A's
        let coli = vec![0i8, 0]; // A, A
        // Column with U's
        let colj = vec![3i8, 3]; // U, U
        let weights = vec![1.0f32, 1.0f32];

        let score = pair_emissioncost(&coli, &colj, &weights, &prior);

        // Should get a valid score (positive or negative)
        // With Watson-Crick pairs and uniform prior, score should be reasonable
        assert!(score != NEGINFINITY);
    }
}
