//! Type definitions for MaxModelMaker algorithm
//!
//! This module implements the MaxMx scoring matrix cell and related types
//! from maxmodelmaker.c lines 40-55.

use crate::types::constants::*;

/// Maximum INSL or INSR path length between MATP nodes
pub const MAXINSERT: usize = 6;

/// Number of scoreable node types in MaxMx (excludes ROOT_NODE)
/// MATP=0, MATL=1, MATR=2, BIFURC=3, BEGINL=4, BEGINR=5
pub const MAXMX_NODETYPES: usize = 6;

/// Node type indices for MaxMx scoring (different from CM node types!)
/// These are used as indices into the sc[] array in MaxMxCell
pub mod maxmx_node {
    /// Match-pair node index in MaxMx
    pub const MATP: usize = 0;
    /// Match-left node index in MaxMx
    pub const MATL: usize = 1;
    /// Match-right node index in MaxMx
    pub const MATR: usize = 2;
    /// Bifurcation node index in MaxMx (also used as END)
    pub const BIFURC: usize = 3;
    /// Begin-left node index in MaxMx
    pub const BEGINL: usize = 4;
    /// Begin-right node index in MaxMx
    pub const BEGINR: usize = 5;
}

/// Scoring matrix cell for MaxModelMaker algorithm
///
/// One per cell of the 2D diagonal matrix of alignment against itself.
/// Each cell stores scores for assigning each possible node type,
/// plus traceback pointers to reconstruct the optimal tree.
///
/// From maxmodelmaker.c struct maxmx_s (lines 40-55)
#[derive(Clone, Debug)]
pub struct MaxMxCell {
    /// Scores of assigning each possible node type
    /// Indices: MATP=0, MATL=1, MATR=2, BIFURC=3, BEGINL=4, BEGINR=5
    pub sc: [i32; MAXMX_NODETYPES],

    // === Traceback info ===

    /// MATP assignment connects to node ftype at (i2, j2)
    pub matp_i2: i16,
    /// MATP j2 coordinate
    pub matp_j2: i16,
    /// MATP following node type
    pub matp_ftype: u8,

    /// MATL assignment connects to node ftype at (i2, j)
    pub matl_i2: i16,
    /// MATL following node type
    pub matl_ftype: u8,

    /// MATR assignment connects to node ftype at (i, j2)
    pub matr_j2: i16,
    /// MATR following node type
    pub matr_ftype: u8,

    /// BEGINL assignment connects to node ftype at (i, j)
    pub begl_ftype: u8,

    /// BEGINR assignment connects to node ftype at (i2, j)
    pub begr_i2: i16,
    /// BEGINR following node type
    pub begr_ftype: u8,

    /// Best bifurcation midpoint: splits into (i, mid) and (mid+1, j)
    pub bifurc_mid: i16,
}

impl MaxMxCell {
    /// Create a new MaxMxCell with default initialization
    ///
    /// All scores are set to NEGINFINITY and traceback pointers
    /// point at the cell itself.
    pub fn new(i: i16, j: i16) -> Self {
        Self {
            sc: [NEGINFINITY; MAXMX_NODETYPES],
            matp_i2: i,
            matp_j2: j,
            matp_ftype: maxmx_node::MATP as u8,
            matl_i2: i,
            matl_ftype: maxmx_node::MATL as u8,
            matr_j2: j,
            matr_ftype: maxmx_node::MATR as u8,
            begl_ftype: maxmx_node::BEGINL as u8,
            begr_i2: i,
            begr_ftype: maxmx_node::BEGINR as u8,
            bifurc_mid: i,
        }
    }
}

impl Default for MaxMxCell {
    fn default() -> Self {
        Self::new(0, 0)
    }
}

/// Transition table type (6x6 matrix of state transition counts/probabilities)
pub type TransTable = [[f64; STATETYPES]; STATETYPES];

/// Create a zeroed transition table
#[inline]
pub fn zero_transtable() -> TransTable {
    [[0.0; STATETYPES]; STATETYPES]
}

/// Copy one transition table to another
#[inline]
pub fn copy_transtable(dest: &mut TransTable, src: &TransTable) {
    for fy in 0..STATETYPES {
        for ty in 0..STATETYPES {
            dest[fy][ty] = src[fy][ty];
        }
    }
}

/// Singlet emission vector type
pub type SingletVec = [f64; ALPHASIZE];

/// Pairwise emission matrix type
pub type PairMatrix = [[f64; ALPHASIZE]; ALPHASIZE];

/// Zero a singlet emission vector
#[inline]
pub fn zero_singlet(vec: &mut SingletVec) {
    for i in 0..ALPHASIZE {
        vec[i] = 0.0;
    }
}

/// Zero a pairwise emission matrix
#[inline]
pub fn zero_pairwise(mx: &mut PairMatrix) {
    for i in 0..ALPHASIZE {
        for j in 0..ALPHASIZE {
            mx[i][j] = 0.0;
        }
    }
}

/// Copy a singlet emission vector
#[inline]
pub fn copy_singlet(dest: &mut SingletVec, src: &SingletVec) {
    for i in 0..ALPHASIZE {
        dest[i] = src[i];
    }
}

/// Copy a pairwise emission matrix
#[inline]
pub fn copy_pairwise(dest: &mut PairMatrix, src: &PairMatrix) {
    for i in 0..ALPHASIZE {
        for j in 0..ALPHASIZE {
            dest[i][j] = src[i][j];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_maxmxcell_new() {
        let cell = MaxMxCell::new(5, 10);

        // Verify scores are NEGINFINITY
        for y in 0..MAXMX_NODETYPES {
            assert_eq!(cell.sc[y], NEGINFINITY);
        }

        // Verify traceback pointers point to self
        assert_eq!(cell.matp_i2, 5);
        assert_eq!(cell.matp_j2, 10);
        assert_eq!(cell.matl_i2, 5);
        assert_eq!(cell.matr_j2, 10);
        assert_eq!(cell.begr_i2, 5);
        assert_eq!(cell.bifurc_mid, 5);
    }

    #[test]
    fn test_transtable_operations() {
        let mut table = zero_transtable();

        // Verify all zeros
        for fy in 0..STATETYPES {
            for ty in 0..STATETYPES {
                assert_eq!(table[fy][ty], 0.0);
            }
        }

        // Modify and copy
        table[0][1] = 1.5;
        table[2][3] = 2.5;

        let mut copy = zero_transtable();
        copy_transtable(&mut copy, &table);

        assert_eq!(copy[0][1], 1.5);
        assert_eq!(copy[2][3], 2.5);
    }

    #[test]
    fn test_singlet_pairwise_operations() {
        let mut singlet: SingletVec = [1.0, 2.0, 3.0, 4.0];
        let mut dest: SingletVec = [0.0; ALPHASIZE];

        copy_singlet(&mut dest, &singlet);
        assert_eq!(dest, [1.0, 2.0, 3.0, 4.0]);

        zero_singlet(&mut singlet);
        assert_eq!(singlet, [0.0; ALPHASIZE]);

        let mut pair: PairMatrix = [[0.0; ALPHASIZE]; ALPHASIZE];
        pair[1][2] = 5.0;

        let mut pair_copy: PairMatrix = [[0.0; ALPHASIZE]; ALPHASIZE];
        copy_pairwise(&mut pair_copy, &pair);
        assert_eq!(pair_copy[1][2], 5.0);

        zero_pairwise(&mut pair);
        assert_eq!(pair[1][2], 0.0);
    }

    #[test]
    fn test_maxmx_node_indices() {
        // Verify node indices are in expected range
        assert!(maxmx_node::MATP < MAXMX_NODETYPES);
        assert!(maxmx_node::MATL < MAXMX_NODETYPES);
        assert!(maxmx_node::MATR < MAXMX_NODETYPES);
        assert!(maxmx_node::BIFURC < MAXMX_NODETYPES);
        assert!(maxmx_node::BEGINL < MAXMX_NODETYPES);
        assert!(maxmx_node::BEGINR < MAXMX_NODETYPES);

        // Verify indices are unique and sequential
        assert_eq!(maxmx_node::MATP, 0);
        assert_eq!(maxmx_node::MATL, 1);
        assert_eq!(maxmx_node::MATR, 2);
        assert_eq!(maxmx_node::BIFURC, 3);
        assert_eq!(maxmx_node::BEGINL, 4);
        assert_eq!(maxmx_node::BEGINR, 5);
    }
}
