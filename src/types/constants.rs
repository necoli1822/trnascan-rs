//! Constants from tRNAscan-SE C codebase (structs.h, squid.h, save.c)
//!
//! This module provides exact constant values matching the original C implementation.
//! All values are cross-verified against structs.h, squid.h, and other header files.

// ============================================================================
// Core Algorithm Constants (structs.h)
// ============================================================================

/// Alphabet size (4 for DNA: A, C, G, T/U)
pub const ALPHASIZE: usize = 4;

/// Number of state types in the covariance model
pub const STATETYPES: usize = 6;

/// Number of node types in the CM tree structure
pub const NODETYPES: usize = 7;

/// Integer precision scaling factor for log probabilities
/// Used to convert floating-point log values to integer representation
pub const INTPRECISION: f64 = 1000.0;

// ============================================================================
// Integer Infinity Values (structs.h lines 238-239)
// ============================================================================

/// Negative infinity for integer log probabilities
/// CRITICAL: Must be exactly -999999 to match C code
pub const NEGINFINITY: i32 = -999999;

/// Positive infinity for integer log probabilities
/// CRITICAL: Must be exactly 999999 to match C code
pub const POSINFINITY: i32 = 999999;

// ============================================================================
// Memory Allocation Constants (structs.h lines 185, 198)
// ============================================================================

/// Block size for trace memory allocation
pub const TMEM_BLOCK: usize = 256;

/// Block size for trace stack allocation
pub const TSTACK_BLOCK: usize = 64;

// ============================================================================
// State Type Indices (structs.h lines 46-52)
// ============================================================================
// These are array indices (0-5) used to index into state arrays

/// Delete state index
pub const DEL_ST: usize = 0;

/// Match-pair state index (emits paired bases)
pub const MATP_ST: usize = 1;

/// Match-left state index (emits left base only)
pub const MATL_ST: usize = 2;

/// Match-right state index (emits right base only)
pub const MATR_ST: usize = 3;

/// Insert-left state index (left insertion)
pub const INSL_ST: usize = 4;

/// Insert-right state index (right insertion)
pub const INSR_ST: usize = 5;

// ============================================================================
// State Type Aliases (structs.h lines 54-56)
// ============================================================================

/// Begin state (alias for DEL_ST)
pub const BEGIN_ST: usize = DEL_ST;

/// Bifurcation state (alias for DEL_ST)
pub const BIFURC_ST: usize = DEL_ST;

/// End state (alias for DEL_ST)
pub const END_ST: usize = DEL_ST;

// ============================================================================
// Unique State Type Flags - BIT MASKS (structs.h lines 61-69)
// ============================================================================
// These are bit flags used for state type identification and filtering

/// Unique delete state flag (bit 0)
pub const U_DEL_ST: u32 = 1 << 0;      // = 1

/// Unique match-pair state flag (bit 1)
pub const U_MATP_ST: u32 = 1 << 1;     // = 2

/// Unique match-left state flag (bit 2)
pub const U_MATL_ST: u32 = 1 << 2;     // = 4

/// Unique match-right state flag (bit 3)
pub const U_MATR_ST: u32 = 1 << 3;     // = 8

/// Unique insert-left state flag (bit 4)
pub const U_INSL_ST: u32 = 1 << 4;     // = 16

/// Unique insert-right state flag (bit 5)
pub const U_INSR_ST: u32 = 1 << 5;     // = 32

/// Unique begin state flag (bit 6)
pub const U_BEGIN_ST: u32 = 1 << 6;    // = 64

/// Unique end state flag (bit 7)
pub const U_END_ST: u32 = 1 << 7;      // = 128

/// Unique bifurcation state flag (bit 8)
pub const U_BIFURC_ST: u32 = 1 << 8;   // = 256

// ============================================================================
// Node Types (structs.h lines 29-36)
// ============================================================================

/// Bifurcation node (splits tree into two subtrees)
pub const BIFURC_NODE: usize = 0;

/// Match-pair node (models Watson-Crick base pair)
pub const MATP_NODE: usize = 1;

/// Match-left node (models single left base)
pub const MATL_NODE: usize = 2;

/// Match-right node (models single right base)
pub const MATR_NODE: usize = 3;

/// Begin-left node (left subtree start)
pub const BEGINL_NODE: usize = 4;

/// Begin-right node (right subtree start)
pub const BEGINR_NODE: usize = 5;

/// Root node (tree root)
pub const ROOT_NODE: usize = 6;

/// End node (alias for BIFURC_NODE)
pub const END_NODE: usize = BIFURC_NODE;

// ============================================================================
// Sequence Format Constants (squid.h lines 96-117)
// ============================================================================

/// Unknown sequence format
pub const K_UNKNOWN: i32 = 0;

/// IntelliGenetics format
pub const K_IG: i32 = 1;

/// GenBank flat file format
pub const K_GENBANK: i32 = 2;

/// EMBL format
pub const K_EMBL: i32 = 3;

/// GCG format
pub const K_GCG: i32 = 4;

/// Strider format
pub const K_STRIDER: i32 = 5;

/// Fitch format
pub const K_FITCH: i32 = 6;

/// Phylip 3.2 format
pub const K_PHYLIP: i32 = 7;

/// Phylip interleaved format
pub const K_PHYLIPI: i32 = 8;

/// ASN.1 format
pub const K_ASN1: i32 = 9;

/// NBRF/PIR format
pub const K_PIR: i32 = 10;

/// Zuker MFOLD format
pub const K_ZUKER: i32 = 11;

/// Multiple sequence alignment format
pub const K_MSF: i32 = 12;

/// Clustal format
pub const K_CLUSTAL: i32 = 13;

/// FASTA format
pub const K_FASTA: i32 = 14;

/// Stockholm format
pub const K_STOCKHOLM: i32 = 15;

/// SELEX format
pub const K_SELEX: i32 = 16;

/// EPS (Encapsulated PostScript) format
pub const K_EPS: i32 = 17;

// ============================================================================
// Binary File Magic Number (save.c line 20)
// ============================================================================

/// Magic number for version 2.0 binary CM files
/// Used to identify and validate binary file format
pub const V20_MAGIC: u32 = 0xe3edb2b0;

// ============================================================================
// ILOG2 Function (structs.h line 240)
// ============================================================================

/// Natural logarithm of 2 (used in ILOG2 calculation)
const LN2: f64 = 0.69314718;

/// Integer log base 2 function - EXACT formula from C code
///
/// Converts a probability to an integer-scaled log2 value.
/// Formula: (log(a) / log(2)) * INTPRECISION
///
/// # Arguments
/// * `a` - Input probability value
///
/// # Returns
/// * Integer-scaled log2 value, or NEGINFINITY if a <= 0
///
/// # Example
/// ```
/// use trnascan_rs::types::constants::ilog2;
/// assert_eq!(ilog2(1.0), 0);       // log2(1) = 0
/// assert_eq!(ilog2(2.0), 1000);    // log2(2) = 1 * 1000
/// assert_eq!(ilog2(0.5), -1000);   // log2(0.5) = -1 * 1000
/// ```
#[inline]
pub fn ilog2(a: f64) -> i32 {
    if a > 0.0 {
        ((a.ln() / LN2) * INTPRECISION) as i32
    } else {
        NEGINFINITY
    }
}

// ============================================================================
// Unit Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ilog2_basic() {
        // Test exact values from golden file
        assert_eq!(ilog2(0.0), NEGINFINITY);
        assert_eq!(ilog2(-1.0), NEGINFINITY);
        assert_eq!(ilog2(0.5), -1000);
        assert_eq!(ilog2(1.0), 0);
        assert_eq!(ilog2(2.0), 1000);
        assert_eq!(ilog2(4.0), 2000);
        assert_eq!(ilog2(8.0), 3000);
    }

    #[test]
    fn test_state_indices() {
        // Verify state type indices are in range [0, STATETYPES)
        assert!(DEL_ST < STATETYPES);
        assert!(MATP_ST < STATETYPES);
        assert!(MATL_ST < STATETYPES);
        assert!(MATR_ST < STATETYPES);
        assert!(INSL_ST < STATETYPES);
        assert!(INSR_ST < STATETYPES);
    }

    #[test]
    fn test_state_aliases() {
        // Verify state aliases
        assert_eq!(BEGIN_ST, DEL_ST);
        assert_eq!(BIFURC_ST, DEL_ST);
        assert_eq!(END_ST, DEL_ST);
    }

    #[test]
    fn test_unique_state_flags() {
        // Verify flags are unique powers of 2
        assert_eq!(U_DEL_ST, 1);
        assert_eq!(U_MATP_ST, 2);
        assert_eq!(U_MATL_ST, 4);
        assert_eq!(U_MATR_ST, 8);
        assert_eq!(U_INSL_ST, 16);
        assert_eq!(U_INSR_ST, 32);
        assert_eq!(U_BEGIN_ST, 64);
        assert_eq!(U_END_ST, 128);
        assert_eq!(U_BIFURC_ST, 256);

        // Verify flags can be combined without collision
        let combined = U_MATP_ST | U_MATL_ST | U_MATR_ST;
        assert_eq!(combined, 14);
    }

    #[test]
    fn test_node_types() {
        // Verify node type indices are in range [0, NODETYPES)
        assert!(BIFURC_NODE < NODETYPES);
        assert!(MATP_NODE < NODETYPES);
        assert!(MATL_NODE < NODETYPES);
        assert!(MATR_NODE < NODETYPES);
        assert!(BEGINL_NODE < NODETYPES);
        assert!(BEGINR_NODE < NODETYPES);
        assert!(ROOT_NODE < NODETYPES);
    }

    #[test]
    fn test_node_alias() {
        assert_eq!(END_NODE, BIFURC_NODE);
    }
}
