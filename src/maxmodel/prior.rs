//! Prior probability distributions for regularization
//!
//! This module implements the prior_s structure and probify functions from
//! structs.h lines 139-151 and probify.c.
//!
//! The prior distributions are used for Bayesian regularization when
//! converting counts to probabilities. This helps deal with small sample
//! statistics problems in model training.

use crate::types::constants::*;

/// Unique state type flags (used in probify functions)
pub const U_MATL_ST: usize = 0;
pub const U_MATR_ST: usize = 1;
pub const U_INSL_ST: usize = 2;
pub const U_INSR_ST: usize = 3;

/// Prior probability distributions for covariance model regularization
///
/// From structs.h struct prior_s (lines 139-151)
#[derive(Clone, Debug)]
pub struct Prior {
    /// State transition priors [from_node_type][to_node_type][from_state][to_state]
    /// Dimensions: [7][4][6][6] (NODETYPES x 4 x STATETYPES x STATETYPES)
    pub tprior: Vec<Vec<[[f64; STATETYPES]; STATETYPES]>>,

    /// MATP_ST pairwise emission prior [left_base][right_base]
    pub matp_prior: [[f64; ALPHASIZE]; ALPHASIZE],

    /// MATL_ST singlet emission prior
    pub matl_prior: [f64; ALPHASIZE],

    /// MATR_ST singlet emission prior
    pub matr_prior: [f64; ALPHASIZE],

    /// INSL_ST singlet emission prior
    pub insl_prior: [f64; ALPHASIZE],

    /// INSR_ST singlet emission prior
    pub insr_prior: [f64; ALPHASIZE],

    /// Alpha weights for state transitions (pseudocount weights)
    pub talpha: [f64; STATETYPES],

    /// Alpha weights for symbol emissions
    pub emalpha: [f64; STATETYPES],

    /// Background symbol frequencies for random model
    pub rfreq: [f64; ALPHASIZE],
}

impl Prior {
    /// Create a new Prior with default uniform distributions
    pub fn new() -> Self {
        // Initialize with uniform priors
        let uniform_singlet = [0.25; ALPHASIZE];
        let mut uniform_pair = [[0.0; ALPHASIZE]; ALPHASIZE];
        for i in 0..ALPHASIZE {
            for j in 0..ALPHASIZE {
                uniform_pair[i][j] = 1.0 / (ALPHASIZE * ALPHASIZE) as f64;
            }
        }

        // Initialize transition priors (7 from-node types, 4 to-node types)
        let mut tprior = Vec::with_capacity(NODETYPES);
        for _ in 0..NODETYPES {
            let mut to_vec = Vec::with_capacity(4);
            for _ in 0..4 {
                // Default: uniform transitions
                let mut trans = [[0.0; STATETYPES]; STATETYPES];
                for from in 0..STATETYPES {
                    for to in 0..STATETYPES {
                        trans[from][to] = 1.0 / STATETYPES as f64;
                    }
                }
                to_vec.push(trans);
            }
            tprior.push(to_vec);
        }

        Self {
            tprior,
            matp_prior: uniform_pair,
            matl_prior: uniform_singlet,
            matr_prior: uniform_singlet,
            insl_prior: uniform_singlet,
            insr_prior: uniform_singlet,
            talpha: [1.0; STATETYPES],
            emalpha: [1.0; STATETYPES],
            rfreq: [0.25; ALPHASIZE],
        }
    }

    /// Create a Prior with biologically-informed priors for tRNA
    ///
    /// This sets priors that reflect expectations for tRNA structure:
    /// - Watson-Crick base pairs are favored
    /// - Insertions/deletions are relatively rare
    pub fn new_trna() -> Self {
        let mut prior = Self::new();

        // Bias toward Watson-Crick pairs (A-U, U-A, G-C, C-G)
        // and wobble pairs (G-U, U-G)
        let wc_prob = 0.20;  // Watson-Crick pairs
        let wobble_prob = 0.05;  // Wobble pairs
        let other_prob = 0.01;  // Other pairs

        for i in 0..ALPHASIZE {
            for j in 0..ALPHASIZE {
                prior.matp_prior[i][j] = other_prob;
            }
        }
        // A-U, U-A (indices: A=0, C=1, G=2, U/T=3)
        prior.matp_prior[0][3] = wc_prob;
        prior.matp_prior[3][0] = wc_prob;
        // G-C, C-G
        prior.matp_prior[2][1] = wc_prob;
        prior.matp_prior[1][2] = wc_prob;
        // G-U, U-G wobble
        prior.matp_prior[2][3] = wobble_prob;
        prior.matp_prior[3][2] = wobble_prob;

        // Normalize pair prior
        let mut sum = 0.0;
        for i in 0..ALPHASIZE {
            for j in 0..ALPHASIZE {
                sum += prior.matp_prior[i][j];
            }
        }
        for i in 0..ALPHASIZE {
            for j in 0..ALPHASIZE {
                prior.matp_prior[i][j] /= sum;
            }
        }

        // Higher alpha for inserts (favor prior distribution)
        prior.emalpha[INSL_ST] = 10.0;
        prior.emalpha[INSR_ST] = 10.0;

        prior
    }

    /// Convert transition matrix from counts to probabilities
    ///
    /// From probify.c ProbifyTransitionMatrix (lines 129-153)
    pub fn probify_transition_matrix(
        &self,
        tmx: &mut [[f64; STATETYPES]; STATETYPES],
        from_node: usize,
        to_node: usize,
    ) {
        let to_idx = to_node.min(3); // Clamp to valid index

        for i in 0..STATETYPES {
            // If no transitions to DEL in prior, this must be an unused vector
            if self.tprior[from_node][to_idx][i][0] > 0.0 {
                let mut denom = 0.0;
                for j in 0..STATETYPES {
                    tmx[i][j] = tmx[i][j]
                        + self.talpha[i] * self.tprior[from_node][to_idx][i][j];
                    denom += tmx[i][j];
                }
                if denom > 0.0 {
                    for j in 0..STATETYPES {
                        tmx[i][j] /= denom;
                    }
                }
            }
        }
    }

    /// Convert singlet emission vector from counts to probabilities
    ///
    /// From probify.c ProbifySingletEmission (lines 167-195)
    pub fn probify_singlet_emission(&self, emvec: &mut [f64; ALPHASIZE], statetype: usize) {
        let em_prior = match statetype {
            U_MATL_ST => &self.matl_prior,
            U_MATR_ST => &self.matr_prior,
            U_INSL_ST => &self.insl_prior,
            U_INSR_ST => &self.insr_prior,
            _ => &self.matl_prior, // Default fallback
        };

        let alpha_idx = match statetype {
            U_MATL_ST => MATL_ST,
            U_MATR_ST => MATR_ST,
            U_INSL_ST => INSL_ST,
            U_INSR_ST => INSR_ST,
            _ => MATL_ST,
        };

        let mut denom = 0.0;
        for x in 0..ALPHASIZE {
            emvec[x] = emvec[x] + self.emalpha[alpha_idx] * em_prior[x];
            denom += emvec[x];
        }
        if denom > 0.0 {
            for x in 0..ALPHASIZE {
                emvec[x] /= denom;
            }
        }
    }

    /// Convert pairwise emission matrix from counts to probabilities
    ///
    /// From probify.c ProbifyPairEmission (lines 207-225)
    pub fn probify_pair_emission(&self, emx: &mut [[f64; ALPHASIZE]; ALPHASIZE]) {
        let mut denom = 0.0;
        for x in 0..ALPHASIZE {
            for y in 0..ALPHASIZE {
                emx[x][y] = emx[x][y] + self.emalpha[MATP_ST] * self.matp_prior[x][y];
                denom += emx[x][y];
            }
        }
        if denom > 0.0 {
            for x in 0..ALPHASIZE {
                for y in 0..ALPHASIZE {
                    emx[x][y] /= denom;
                }
            }
        }
    }
}

impl Default for Prior {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_prior_new() {
        let prior = Prior::new();

        // Check uniform singlet priors
        for x in 0..ALPHASIZE {
            assert!((prior.matl_prior[x] - 0.25).abs() < 1e-10);
            assert!((prior.matr_prior[x] - 0.25).abs() < 1e-10);
            assert!((prior.insl_prior[x] - 0.25).abs() < 1e-10);
            assert!((prior.insr_prior[x] - 0.25).abs() < 1e-10);
            assert!((prior.rfreq[x] - 0.25).abs() < 1e-10);
        }

        // Check alpha values
        for st in 0..STATETYPES {
            assert_eq!(prior.talpha[st], 1.0);
            assert_eq!(prior.emalpha[st], 1.0);
        }
    }

    #[test]
    fn test_prior_trna() {
        let prior = Prior::new_trna();

        // Watson-Crick pairs should have higher probability
        // A-U (0-3) should be higher than A-A (0-0)
        assert!(prior.matp_prior[0][3] > prior.matp_prior[0][0]);

        // Insert alphas should be higher
        assert!(prior.emalpha[INSL_ST] > prior.emalpha[MATL_ST]);
    }

    #[test]
    fn test_probify_singlet() {
        let prior = Prior::new();
        let mut emvec = [10.0, 5.0, 3.0, 2.0]; // Raw counts

        prior.probify_singlet_emission(&mut emvec, U_MATL_ST);

        // Should sum to 1.0
        let sum: f64 = emvec.iter().sum();
        assert!((sum - 1.0).abs() < 1e-10);

        // Higher counts should still have higher probability
        assert!(emvec[0] > emvec[1]);
        assert!(emvec[1] > emvec[2]);
        assert!(emvec[2] > emvec[3]);
    }

    #[test]
    fn test_probify_pair() {
        let prior = Prior::new();
        let mut emx = [[0.0; ALPHASIZE]; ALPHASIZE];

        // Add some counts
        emx[0][3] = 10.0; // A-U
        emx[2][1] = 8.0;  // G-C

        prior.probify_pair_emission(&mut emx);

        // Should sum to 1.0
        let mut sum = 0.0;
        for x in 0..ALPHASIZE {
            for y in 0..ALPHASIZE {
                sum += emx[x][y];
            }
        }
        assert!((sum - 1.0).abs() < 1e-10);

        // Higher counts should have higher probability
        assert!(emx[0][3] > emx[0][0]);
    }

    #[test]
    fn test_probify_transition() {
        let prior = Prior::new();
        let mut tmx = [[0.0; STATETYPES]; STATETYPES];

        // Add some transition counts
        tmx[MATL_ST][END_ST] = 10.0;
        tmx[DEL_ST][END_ST] = 2.0;

        prior.probify_transition_matrix(&mut tmx, MATL_NODE, END_NODE);

        // Row sums should be 1.0 or 0.0 (for unused rows)
        for i in 0..STATETYPES {
            let row_sum: f64 = tmx[i].iter().sum();
            if row_sum > 0.0 {
                assert!((row_sum - 1.0).abs() < 1e-10);
            }
        }
    }
}
