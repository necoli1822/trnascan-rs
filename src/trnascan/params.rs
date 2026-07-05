//! Search parameters for tRNAscan first-pass scanning.
//!
//! This module implements the Fichant-Burks algorithm parameter sets
//! as described in J. Mol. Biol. (1991) 220:659-671.

/// Search parameters for tRNAscan first-pass.
///
/// These parameters control the sensitivity and specificity of the
/// tRNA gene detection algorithm.
#[derive(Debug, Clone, Copy)]
pub struct SearchParams {
    /// T-Psi-C signal sequence matrix score threshold (default 0.40)
    pub tpc_sig_thresh: f32,
    /// D signal sequence matrix score threshold (default 0.40 strict, 0.30 relaxed)
    pub d_sig_thresh: f32,
    /// General score (SG) cutoff for final tRNA prediction
    pub sg_cutoff: i32,
    /// Number of TPC matrix invariant bases allowed NOT to match
    pub tpc_inv: i32,
    /// Number of base pairs in TPC stem required to increment SG
    pub tpc_incsg: i32,
    /// Number of base pairs in TPC stem required to keep candidate
    pub tpc_keep: i32,
    /// Number of D matrix invariant bases allowed NOT to match
    pub d_inv: i32,
    /// Minimum SG required to begin looking for anticodon loop
    pub look_for_acloop_sg: i32,
    /// Minimum base pairs required in anticodon loop
    pub acloop_min: i32,
    /// Number of base pairs in amino acyl stem needed to increment SG
    pub aa_incsg: i32,
    /// Number of base pairs in amino acyl stem needed to keep candidate
    pub aa_keep: i32,
}

impl Default for SearchParams {
    fn default() -> Self {
        Self::strict()
    }
}

impl SearchParams {
    /// Create strict (original tRNAscan 1.3) parameters.
    ///
    /// These are the default parameters from the original tRNAscan publication.
    pub fn strict() -> Self {
        Self {
            tpc_sig_thresh: 0.40,
            d_sig_thresh: 0.40,
            sg_cutoff: 5,
            tpc_inv: 2,
            tpc_incsg: 5,
            tpc_keep: 4,
            d_inv: 1,
            look_for_acloop_sg: 4,
            acloop_min: 4,
            aa_incsg: 7,
            aa_keep: 6,
        }
    }

    /// Create relaxed parameters for use with tRNAscan-SE.
    ///
    /// These parameters make tRNAscan into a rough pre-filter for
    /// the covariance model-based tRNA prediction.
    pub fn relaxed() -> Self {
        Self {
            tpc_sig_thresh: 0.40,
            d_sig_thresh: 0.30,
            sg_cutoff: 5,
            tpc_inv: 2,
            tpc_incsg: 4,
            tpc_keep: 2,
            d_inv: 2,
            look_for_acloop_sg: 3,
            acloop_min: 3,
            aa_incsg: 5,
            aa_keep: 4,
        }
    }

    /// Create alternate parameters for experimentation.
    pub fn alternate() -> Self {
        Self {
            tpc_sig_thresh: 0.40,
            d_sig_thresh: 0.30,
            sg_cutoff: 4,
            tpc_inv: 2,
            tpc_incsg: 4,
            tpc_keep: 2,
            d_inv: 2,
            look_for_acloop_sg: 3,
            acloop_min: 3,
            aa_incsg: 6,
            aa_keep: 4,
        }
    }
}

/// Constants for tRNA structure constraints
pub const MIN_VAR_LOOP: i32 = 28;
pub const MAX_INTRON_LEN: i32 = 60;
pub const MIN_SEQ_LEN: usize = 70;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_strict_params() {
        let params = SearchParams::strict();
        assert_eq!(params.sg_cutoff, 5);
        assert!((params.tpc_sig_thresh - 0.40).abs() < 0.001);
        assert!((params.d_sig_thresh - 0.40).abs() < 0.001);
    }

    #[test]
    fn test_relaxed_params() {
        let params = SearchParams::relaxed();
        assert_eq!(params.sg_cutoff, 5);
        assert!((params.d_sig_thresh - 0.30).abs() < 0.001);
        assert_eq!(params.tpc_keep, 2);
    }

    #[test]
    fn test_alternate_params() {
        let params = SearchParams::alternate();
        assert_eq!(params.sg_cutoff, 4);
    }
}
