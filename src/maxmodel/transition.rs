//! Transition table construction functions for MaxModelMaker
//!
//! This module implements the from*_transtable and to*_transtable functions
//! from maxmodelmaker.c lines 1183-1717.
//!
//! These functions calculate state transition count tables for different
//! node type assignments and convert master tables to specific node type tables.

use crate::maxmodel::types::*;
use crate::types::constants::*;

/// Assign state type based on symbols at positions i and j
///
/// Given a cell (i,j) assigned to MATP for the whole alignment,
/// returns the actual state assignment for a sequence based on gaps.
///
/// From maxmodelmaker.c assign_cell (lines 1587-1598)
///
/// # Arguments
/// * `i` - Column i coordinate
/// * `j` - Column j coordinate
/// * `symi` - Symbol at position i (-1 for gap, 0-3 for base)
/// * `symj` - Symbol at position j (-1 for gap, 0-3 for base)
///
/// # Returns
/// State type: DEL_ST, MATP_ST, MATR_ST, MATL_ST, or END_ST
#[inline]
pub fn assign_cell(i: usize, j: usize, symi: i8, symj: i8) -> usize {
    if i > j {
        END_ST
    } else if symi >= 0 {
        if symj >= 0 {
            MATP_ST
        } else {
            MATL_ST
        }
    } else if symj >= 0 {
        MATR_ST
    } else {
        DEL_ST
    }
}

/// Calculate transition table from MATP node at (i,j) to MATP at (i2,j2)
///
/// From maxmodelmaker.c frommatp_transtable (lines 1183-1240)
///
/// # Arguments
/// * `aseqs_t` - Transposed alignment [1..alen][0..nseq-1]
/// * `weights` - Sequence weights
/// * `i`, `j` - Starting cell coordinates
/// * `i2`, `j2` - Ending cell coordinates
/// * `accum_insl` - Accumulated INSL counts from i to i2 for each sequence
/// * `accum_insr` - Accumulated INSR counts from j2 to j for each sequence
/// * `trans` - Output transition table (filled in)
pub fn frommatp_transtable(
    aseqs_t: &[Vec<i8>],
    weights: &[f32],
    i: usize,
    j: usize,
    i2: usize,
    j2: usize,
    accum_insl: &[i32],
    accum_insr: &[i32],
    trans: &mut TransTable,
) {
    let nseq = weights.len();

    // Zero the counter array
    *trans = zero_transtable();

    // For each sequence, assign states and accumulate transitions
    for idx in 0..nseq {
        let fy = assign_cell(i, j, aseqs_t[i][idx], aseqs_t[j][idx]);
        let ty = assign_cell(i2, j2, aseqs_t[i2][idx], aseqs_t[j2][idx]);

        let insl = accum_insl[idx];
        let insr = accum_insr[idx];
        let w = weights[idx] as f64;

        if insl == 0 && insr == 0 {
            trans[fy][ty] += w;
        } else if insl > 0 {
            trans[fy][INSL_ST] += w;
            trans[INSL_ST][INSL_ST] += (insl - 1) as f64 * w;
            if insr > 0 {
                trans[INSL_ST][INSR_ST] += w;
                trans[INSR_ST][INSR_ST] += (insr - 1) as f64 * w;
                trans[INSR_ST][ty] += w;
            } else {
                trans[INSL_ST][ty] += w;
            }
        } else if insr > 0 {
            trans[fy][INSR_ST] += w;
            trans[INSR_ST][INSR_ST] += (insr - 1) as f64 * w;
            trans[INSR_ST][ty] += w;
        }
    }
}

/// Calculate transition table from MATL node at (i,j) to MATP at (i2,j)
///
/// From maxmodelmaker.c frommatl_transtable (lines 1266-1307)
pub fn frommatl_transtable(
    aseqs_t: &[Vec<i8>],
    weights: &[f32],
    i: usize,
    j: usize,
    i2: usize,
    accum_insl: &[i32],
    trans: &mut TransTable,
) {
    let nseq = weights.len();

    // Zero the counter array
    *trans = zero_transtable();

    for idx in 0..nseq {
        let fy = if aseqs_t[i][idx] == -1 { DEL_ST } else { MATL_ST };
        let ty = assign_cell(i2, j, aseqs_t[i2][idx], aseqs_t[j][idx]);

        let insl = accum_insl[idx];
        let w = weights[idx] as f64;

        if insl == 0 {
            trans[fy][ty] += w;
        } else {
            trans[fy][INSL_ST] += w;
            trans[INSL_ST][INSL_ST] += (insl - 1) as f64 * w;
            trans[INSL_ST][ty] += w;
        }
    }
}

/// Calculate transition table from MATR node at (i,j) to MATP at (i,j2)
///
/// From maxmodelmaker.c frommatr_transtable (lines 1333-1373)
pub fn frommatr_transtable(
    aseqs_t: &[Vec<i8>],
    weights: &[f32],
    i: usize,
    j: usize,
    j2: usize,
    accum_insr: &[i32],
    trans: &mut TransTable,
) {
    let nseq = weights.len();

    // Zero the counter array
    *trans = zero_transtable();

    for idx in 0..nseq {
        let fy = if aseqs_t[j][idx] == -1 { DEL_ST } else { MATR_ST };
        let ty = assign_cell(i, j2, aseqs_t[i][idx], aseqs_t[j2][idx]);

        let insr = accum_insr[idx];
        let w = weights[idx] as f64;

        if insr == 0 {
            trans[fy][ty] += w;
        } else {
            trans[fy][INSR_ST] += w;
            trans[INSR_ST][INSR_ST] += (insr - 1) as f64 * w;
            trans[INSR_ST][ty] += w;
        }
    }
}

/// Calculate transition table from BEGINR node at (i,j) to MATP at (i2,j)
///
/// From maxmodelmaker.c frombeginr_transtable (lines 1399-1437)
pub fn frombeginr_transtable(
    aseqs_t: &[Vec<i8>],
    weights: &[f32],
    j: usize,
    i2: usize,
    accum_insl: &[i32],
    trans: &mut TransTable,
) {
    let nseq = weights.len();

    // Zero the counter array
    *trans = zero_transtable();

    let fy = BEGIN_ST;

    for idx in 0..nseq {
        let ty = assign_cell(i2, j, aseqs_t[i2][idx], aseqs_t[j][idx]);

        let insl = accum_insl[idx];
        let w = weights[idx] as f64;

        if insl == 0 {
            trans[fy][ty] += w;
        } else {
            trans[fy][INSL_ST] += w;
            trans[INSL_ST][INSL_ST] += (insl - 1) as f64 * w;
            trans[INSL_ST][ty] += w;
        }
    }
}

/// Calculate transition table from BEGINL node at (i,j) to same cell
///
/// From maxmodelmaker.c frombeginl_transtable (lines 1459-1486)
pub fn frombeginl_transtable(
    aseqs_t: &[Vec<i8>],
    weights: &[f32],
    i: usize,
    j: usize,
    trans: &mut TransTable,
) {
    let nseq = weights.len();

    // Zero the counter array
    *trans = zero_transtable();

    let fy = BEGIN_ST;

    for idx in 0..nseq {
        let ty = assign_cell(i, j, aseqs_t[i][idx], aseqs_t[j][idx]);
        trans[fy][ty] += weights[idx] as f64;
    }
}

/// Calculate transition table from ROOT to (i2,j2)
///
/// From maxmodelmaker.c fromroot_transtable (lines 1511-1564)
pub fn fromroot_transtable(
    aseqs_t: &[Vec<i8>],
    weights: &[f32],
    i2: usize,
    j2: usize,
    accum_insl: &[i32],
    accum_insr: &[i32],
    trans: &mut TransTable,
) {
    let nseq = weights.len();

    // Zero the counter array
    *trans = zero_transtable();

    for idx in 0..nseq {
        let ty = assign_cell(i2, j2, aseqs_t[i2][idx], aseqs_t[j2][idx]);

        let insl = accum_insl[idx];
        let insr = accum_insr[idx];
        let w = weights[idx] as f64;

        if insl == 0 && insr == 0 {
            trans[BEGIN_ST][ty] += w;
        } else if insl > 0 {
            trans[BEGIN_ST][INSL_ST] += w;
            trans[INSL_ST][INSL_ST] += (insl - 1) as f64 * w;
            if insr > 0 {
                trans[INSL_ST][INSR_ST] += w;
                trans[INSR_ST][INSR_ST] += (insr - 1) as f64 * w;
                trans[INSR_ST][ty] += w;
            } else {
                trans[INSL_ST][ty] += w;
            }
        } else if insr > 0 {
            trans[BEGIN_ST][INSR_ST] += w;
            trans[INSR_ST][INSR_ST] += (insr - 1) as f64 * w;
            trans[INSR_ST][ty] += w;
        }
    }
}

/// Convert master table to MATP target node table (identity copy)
///
/// From maxmodelmaker.c to_matp_transtable (lines 1615-1626)
#[inline]
pub fn to_matp_transtable(master_table: &TransTable, trans: &mut TransTable) {
    copy_transtable(trans, master_table);
}

/// Convert master table to MATR target node table
///
/// Combines certain state transitions for MATR assignment.
///
/// From maxmodelmaker.c to_matr_transtable (lines 1641-1656)
pub fn to_matr_transtable(master_table: &TransTable, trans: &mut TransTable) {
    for fy in 0..STATETYPES {
        // DEL absorbs MATL (both have gap on right)
        trans[fy][DEL_ST] = master_table[fy][DEL_ST] + master_table[fy][MATL_ST];
        trans[fy][MATP_ST] = 0.0;
        trans[fy][MATL_ST] = 0.0;
        // MATR absorbs MATP (both have emission on right)
        trans[fy][MATR_ST] = master_table[fy][MATR_ST] + master_table[fy][MATP_ST];
        trans[fy][INSL_ST] = master_table[fy][INSL_ST];
        trans[fy][INSR_ST] = master_table[fy][INSR_ST];
    }
}

/// Convert master table to MATL target node table
///
/// Combines certain state transitions for MATL assignment.
///
/// From maxmodelmaker.c to_matl_transtable (lines 1671-1686)
pub fn to_matl_transtable(master_table: &TransTable, trans: &mut TransTable) {
    for fy in 0..STATETYPES {
        // DEL absorbs MATR (both have gap on left)
        trans[fy][DEL_ST] = master_table[fy][DEL_ST] + master_table[fy][MATR_ST];
        trans[fy][MATP_ST] = 0.0;
        // MATL absorbs MATP (both have emission on left)
        trans[fy][MATL_ST] = master_table[fy][MATL_ST] + master_table[fy][MATP_ST];
        trans[fy][MATR_ST] = 0.0;
        trans[fy][INSL_ST] = master_table[fy][INSL_ST];
        trans[fy][INSR_ST] = master_table[fy][INSR_ST];
    }
}

/// Convert master table to BIFURC target node table
///
/// Combines all match states into BIFURC.
///
/// From maxmodelmaker.c to_bifurc_transtable (lines 1701-1717)
pub fn to_bifurc_transtable(master_table: &TransTable, trans: &mut TransTable) {
    for fy in 0..STATETYPES {
        // BIFURC absorbs all match/delete states
        trans[fy][BIFURC_ST] = master_table[fy][DEL_ST]
            + master_table[fy][MATR_ST]
            + master_table[fy][MATL_ST]
            + master_table[fy][MATP_ST];
        trans[fy][MATP_ST] = 0.0;
        trans[fy][MATL_ST] = 0.0;
        trans[fy][MATR_ST] = 0.0;
        trans[fy][INSL_ST] = master_table[fy][INSL_ST];
        trans[fy][INSR_ST] = master_table[fy][INSR_ST];
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_assign_cell() {
        // Both symbols present -> MATP
        assert_eq!(assign_cell(1, 5, 0, 2), MATP_ST);

        // Gap at j -> MATL
        assert_eq!(assign_cell(1, 5, 0, -1), MATL_ST);

        // Gap at i -> MATR
        assert_eq!(assign_cell(1, 5, -1, 2), MATR_ST);

        // Both gaps -> DEL
        assert_eq!(assign_cell(1, 5, -1, -1), DEL_ST);

        // i > j -> END
        assert_eq!(assign_cell(5, 3, 0, 0), END_ST);
    }

    #[test]
    fn test_to_matr_transtable() {
        let mut master = zero_transtable();
        master[MATP_ST][MATP_ST] = 5.0;
        master[MATP_ST][MATL_ST] = 3.0;
        master[MATP_ST][MATR_ST] = 2.0;
        master[MATP_ST][DEL_ST] = 1.0;

        let mut trans = zero_transtable();
        to_matr_transtable(&master, &mut trans);

        // MATP absorbs into MATR
        assert_eq!(trans[MATP_ST][MATR_ST], 7.0); // 5 + 2
        // MATL absorbs into DEL
        assert_eq!(trans[MATP_ST][DEL_ST], 4.0); // 1 + 3
        // MATP and MATL zeroed
        assert_eq!(trans[MATP_ST][MATP_ST], 0.0);
        assert_eq!(trans[MATP_ST][MATL_ST], 0.0);
    }

    #[test]
    fn test_to_matl_transtable() {
        let mut master = zero_transtable();
        master[MATP_ST][MATP_ST] = 5.0;
        master[MATP_ST][MATL_ST] = 3.0;
        master[MATP_ST][MATR_ST] = 2.0;
        master[MATP_ST][DEL_ST] = 1.0;

        let mut trans = zero_transtable();
        to_matl_transtable(&master, &mut trans);

        // MATP absorbs into MATL
        assert_eq!(trans[MATP_ST][MATL_ST], 8.0); // 5 + 3
        // MATR absorbs into DEL
        assert_eq!(trans[MATP_ST][DEL_ST], 3.0); // 1 + 2
        // MATP and MATR zeroed
        assert_eq!(trans[MATP_ST][MATP_ST], 0.0);
        assert_eq!(trans[MATP_ST][MATR_ST], 0.0);
    }

    #[test]
    fn test_to_bifurc_transtable() {
        let mut master = zero_transtable();
        master[MATP_ST][MATP_ST] = 5.0;
        master[MATP_ST][MATL_ST] = 3.0;
        master[MATP_ST][MATR_ST] = 2.0;
        master[MATP_ST][DEL_ST] = 1.0;

        let mut trans = zero_transtable();
        to_bifurc_transtable(&master, &mut trans);

        // All absorb into BIFURC
        assert_eq!(trans[MATP_ST][BIFURC_ST], 11.0); // 5 + 3 + 2 + 1
        // All others zeroed
        assert_eq!(trans[MATP_ST][MATP_ST], 0.0);
        assert_eq!(trans[MATP_ST][MATL_ST], 0.0);
        assert_eq!(trans[MATP_ST][MATR_ST], 0.0);
    }

    #[test]
    fn test_frombeginl_transtable() {
        // Simple 2-sequence alignment
        let aseqs_t = vec![
            vec![-1i8, -1],  // Guard column 0
            vec![0, 0],      // Column 1: A, A
            vec![2, 2],      // Column 2: G, G
            vec![-1, -1],    // Guard column 3
        ];
        let weights = vec![1.0f32, 1.0];

        let mut trans = zero_transtable();
        frombeginl_transtable(&aseqs_t, &weights, 1, 2, &mut trans);

        // Both sequences have symbols at both positions -> MATP
        assert_eq!(trans[BEGIN_ST][MATP_ST], 2.0);
    }
}
