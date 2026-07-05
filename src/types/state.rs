//! State structures for CM alignment algorithms
//!
//! This module implements IState and PState structures from structs.h lines 102-128.
//! These are optimized representations of CM states used in alignment algorithms.

use crate::types::constants::*;

/// Integer state structure from structs.h lines 111-118
///
/// In the alignment algorithms, a CM is converted to an array of states.
/// Each state contains probability info as integers instead of floating point.
///
/// CRITICAL DIFFERENCES FROM pstate_s:
/// - NO bifr field (only in pstate_s)
/// - Integer log-odds probabilities (vs floats in pstate_s)
/// - Rearranged transition order: INSL, INSR, DEL, MATP, MATL, MATR
///
/// The transition vector order is different than in the CM:
/// INSL and INSR are first: INSL, INSR, DEL, MATP, MATL, MATR
#[derive(Clone, Debug)]
pub struct IState {
    /// Index of node this state belongs to
    pub nodeidx: i32,

    /// Unique ID for type of this state (U_MATP_ST, U_MATL_ST, etc.)
    /// These are bit flags from constants.rs (U_DEL_ST = 1, U_MATP_ST = 2, etc.)
    pub statetype: u32,

    /// Offset in state array to first INS state
    pub offset: i32,

    /// Number of elements in tmx transition vector
    pub connectnum: i32,

    /// Rearranged transition vector (integer log-odds)
    /// Order: INSL, INSR, DEL, MATP, MATL, MATR
    pub tmx: [i32; STATETYPES],

    /// Integer log-odds emission vector
    /// Size: 4 for single emissions (MATL, MATR, INSL, INSR)
    ///       16 for pair emissions (MATP)
    pub emit: [i32; ALPHASIZE * ALPHASIZE],
}

impl IState {
    /// Create a new IState with default values
    pub fn new() -> Self {
        Self {
            nodeidx: 0,
            statetype: 0,
            offset: 0,
            connectnum: 0,
            tmx: [NEGINFINITY; STATETYPES],
            emit: [NEGINFINITY; ALPHASIZE * ALPHASIZE],
        }
    }

    /// Check if this state is of a specific type
    pub fn is_type(&self, state_flag: u32) -> bool {
        self.statetype & state_flag != 0
    }

    /// Get emission probability for a single base
    /// Used for MATL, MATR, INSL, INSR states
    pub fn single_emit(&self, base: usize) -> i32 {
        if base < ALPHASIZE {
            self.emit[base]
        } else {
            NEGINFINITY
        }
    }

    /// Get emission probability for a base pair
    /// Used for MATP states
    pub fn pair_emit(&self, left_base: usize, right_base: usize) -> i32 {
        if left_base < ALPHASIZE && right_base < ALPHASIZE {
            self.emit[left_base * ALPHASIZE + right_base]
        } else {
            NEGINFINITY
        }
    }
}

impl Default for IState {
    fn default() -> Self {
        Self::new()
    }
}

/// Probability state structure from structs.h lines 120-128
///
/// Similar to IState but with floating-point probabilities instead of integers.
///
/// CRITICAL DIFFERENCE FROM istate_s:
/// - HAS bifr field (NOT in istate_s!)
/// - Floating-point probabilities (vs integers in istate_s)
#[derive(Clone, Debug)]
pub struct PState {
    /// Index of node this state belongs to
    pub nodeidx: i32,

    /// Unique ID for type of this state (U_MATP_ST, U_MATL_ST, etc.)
    pub statetype: u32,

    /// Offset in state array to first INS state
    pub offset: i32,

    /// Number of elements in tmx transition vector
    pub connectnum: i32,

    /// (U_BIFURC_ST only) Index of right connection
    /// CRITICAL: This field exists ONLY in pstate_s, NOT in istate_s!
    pub bifr: i32,

    /// Rearranged transition vector (floating-point probabilities)
    /// Order: INSL, INSR, DEL, MATP, MATL, MATR
    pub tmx: [f64; STATETYPES],

    /// Emission vector (floating-point probabilities)
    /// Size: 4 for single emissions, 16 for pair emissions
    pub emit: [f64; ALPHASIZE * ALPHASIZE],
}

impl PState {
    /// Create a new PState with default values
    pub fn new() -> Self {
        Self {
            nodeidx: 0,
            statetype: 0,
            offset: 0,
            connectnum: 0,
            bifr: 0,
            tmx: [0.0; STATETYPES],
            emit: [0.0; ALPHASIZE * ALPHASIZE],
        }
    }

    /// Check if this state is of a specific type
    pub fn is_type(&self, state_flag: u32) -> bool {
        self.statetype & state_flag != 0
    }

    /// Get emission probability for a single base
    pub fn single_emit(&self, base: usize) -> f64 {
        if base < ALPHASIZE {
            self.emit[base]
        } else {
            0.0
        }
    }

    /// Get emission probability for a base pair
    pub fn pair_emit(&self, left_base: usize, right_base: usize) -> f64 {
        if left_base < ALPHASIZE && right_base < ALPHASIZE {
            self.emit[left_base * ALPHASIZE + right_base]
        } else {
            0.0
        }
    }
}

impl Default for PState {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_istate_creation() {
        let istate = IState::new();
        assert_eq!(istate.nodeidx, 0);
        assert_eq!(istate.statetype, 0);
        assert_eq!(istate.offset, 0);
        assert_eq!(istate.connectnum, 0);

        // Verify all transitions initialized to NEGINFINITY
        for i in 0..STATETYPES {
            assert_eq!(istate.tmx[i], NEGINFINITY);
        }

        // Verify all emissions initialized to NEGINFINITY
        for i in 0..(ALPHASIZE * ALPHASIZE) {
            assert_eq!(istate.emit[i], NEGINFINITY);
        }
    }

    #[test]
    fn test_pstate_creation() {
        let pstate = PState::new();
        assert_eq!(pstate.nodeidx, 0);
        assert_eq!(pstate.statetype, 0);
        assert_eq!(pstate.offset, 0);
        assert_eq!(pstate.connectnum, 0);
        assert_eq!(pstate.bifr, 0); // ONLY in PState!

        // Verify all probabilities initialized to 0.0
        for i in 0..STATETYPES {
            assert_eq!(pstate.tmx[i], 0.0);
        }

        for i in 0..(ALPHASIZE * ALPHASIZE) {
            assert_eq!(pstate.emit[i], 0.0);
        }
    }

    #[test]
    fn test_istate_vs_pstate_bifr() {
        // Verify IState does NOT have bifr field
        let _istate = IState::new();
        // Should not compile: istate.bifr (field doesn't exist)

        // Verify PState DOES have bifr field
        let mut pstate = PState::new();
        pstate.bifr = 42;
        assert_eq!(pstate.bifr, 42);
    }

    #[test]
    fn test_state_type_checking() {
        let mut istate = IState::new();
        istate.statetype = U_MATP_ST;
        assert!(istate.is_type(U_MATP_ST));
        assert!(!istate.is_type(U_MATL_ST));

        // Test flag combinations
        istate.statetype = U_MATP_ST | U_MATL_ST;
        assert!(istate.is_type(U_MATP_ST));
        assert!(istate.is_type(U_MATL_ST));
        assert!(!istate.is_type(U_INSL_ST));
    }

    #[test]
    fn test_istate_emissions() {
        let mut istate = IState::new();

        // Test single emission
        istate.emit[0] = 100;
        istate.emit[1] = 200;
        assert_eq!(istate.single_emit(0), 100);
        assert_eq!(istate.single_emit(1), 200);
        assert_eq!(istate.single_emit(99), NEGINFINITY);

        // Test pair emission
        istate.emit[0 * ALPHASIZE + 0] = 1000; // A-A
        istate.emit[0 * ALPHASIZE + 3] = 2000; // A-T
        assert_eq!(istate.pair_emit(0, 0), 1000);
        assert_eq!(istate.pair_emit(0, 3), 2000);
        assert_eq!(istate.pair_emit(99, 0), NEGINFINITY);
    }

    #[test]
    fn test_pstate_emissions() {
        let mut pstate = PState::new();

        // Test single emission
        pstate.emit[0] = 0.25;
        pstate.emit[1] = 0.5;
        assert_eq!(pstate.single_emit(0), 0.25);
        assert_eq!(pstate.single_emit(1), 0.5);
        assert_eq!(pstate.single_emit(99), 0.0);

        // Test pair emission
        pstate.emit[0 * ALPHASIZE + 0] = 0.1;
        pstate.emit[0 * ALPHASIZE + 3] = 0.2;
        assert_eq!(pstate.pair_emit(0, 0), 0.1);
        assert_eq!(pstate.pair_emit(0, 3), 0.2);
        assert_eq!(pstate.pair_emit(99, 0), 0.0);
    }

    #[test]
    fn test_state_array_sizes() {
        let istate = IState::default();
        assert_eq!(istate.tmx.len(), STATETYPES);
        assert_eq!(istate.emit.len(), ALPHASIZE * ALPHASIZE);

        let pstate = PState::default();
        assert_eq!(pstate.tmx.len(), STATETYPES);
        assert_eq!(pstate.emit.len(), ALPHASIZE * ALPHASIZE);
    }
}
