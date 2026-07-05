//! MaxMx scoring matrix allocation and initialization
//!
//! This module implements the alloc_maxmx and init_maxmx functions from
//! maxmodelmaker.c lines 270-417.

use crate::maxmodel::prior::Prior;
use crate::maxmodel::types::*;
use crate::types::constants::*;

/// MaxMx scoring matrix
///
/// Lower diagonal matrix with inverted indexing: mmx[j][i] where i <= j+1.
/// Row j has j+2 columns (indices 0 to j+1).
///
/// The extra off-diagonal column (j+1, j) stores END node boundary conditions.
pub struct MaxMx {
    /// Matrix storage: mmx[j] is row j with j+2 cells
    pub data: Vec<Vec<MaxMxCell>>,
    /// Alignment length
    pub alen: usize,
}

impl MaxMx {
    /// Allocate and initialize a new MaxMx matrix
    ///
    /// Creates a lower diagonal matrix for alignment of length `alen`.
    /// Matrix is indexed as mmx[j][i] where:
    /// - j ranges from 0 to alen (inclusive)
    /// - i ranges from 0 to j+1 (inclusive)
    ///
    /// From maxmodelmaker.c alloc_maxmx (lines 270-304)
    pub fn new(alen: usize) -> Self {
        let mut data = Vec::with_capacity(alen + 1);

        for j in 0..=alen {
            let mut row = Vec::with_capacity(j + 2);
            for i in 0..=(j + 1) {
                // Initialize cell with traceback pointers pointing to itself
                row.push(MaxMxCell::new(i as i16, j as i16));
            }
            data.push(row);
        }

        Self { data, alen }
    }

    /// Get cell reference at (j, i)
    #[inline]
    pub fn get(&self, j: usize, i: usize) -> &MaxMxCell {
        &self.data[j][i]
    }

    /// Get mutable cell reference at (j, i)
    #[inline]
    pub fn get_mut(&mut self, j: usize, i: usize) -> &mut MaxMxCell {
        &mut self.data[j][i]
    }

    /// Initialize the off-diagonal and diagonal boundary conditions
    ///
    /// Off-diagonal (j+1, j): END score = 0
    /// Diagonal (j, j): MATL calculated, BEGINL/BEGINR derived
    ///
    /// From maxmodelmaker.c init_maxmx (lines 331-417)
    pub fn initialize(
        &mut self,
        nseq: usize,
        prior: &Prior,
        mscore: &[i32],
        gapcount: &[f64],
    ) {
        // Do the offdiagonal (j+1, j)
        // Set BIFURC/END alignment costs to zero
        // Everything else is left at NEGINFINITY
        for j in 0..=self.alen {
            self.data[j][j + 1].sc[maxmx_node::BIFURC] = 0;
        }

        // Do the diagonal (j, j)
        // MATL is calculated; then BEGINL, BEGINR are calculated
        for j in 1..=self.alen {
            let nseq_f = nseq as f64;

            // Make transition matrix for MATL -> END
            let mut trans = zero_transtable();
            trans[MATL_ST][END_ST] = nseq_f - gapcount[j];
            trans[DEL_ST][END_ST] = gapcount[j];

            // Apply prior regularization
            prior.probify_transition_matrix(&mut trans, MATL_NODE, END_NODE);

            // Score = sum P(j | MATL) + sum T(END | j,j,(DEL|MATL))
            let matl_trans_score = ((trans[MATL_ST][END_ST].ln() * (nseq_f - gapcount[j]))
                + (trans[DEL_ST][END_ST].ln() * gapcount[j]))
                * INTPRECISION;

            self.data[j][j].sc[maxmx_node::MATL] = mscore[j] + matl_trans_score as i32;
            self.data[j][j].matl_i2 = (j + 1) as i16;
            self.data[j][j].matl_ftype = END_NODE as u8;

            // MATR_NODE scores are exactly the same as MATL_NODE on diagonal
            self.data[j][j].sc[maxmx_node::MATR] = self.data[j][j].sc[maxmx_node::MATL];
            self.data[j][j].matr_j2 = (j as i16) - 1;
            self.data[j][j].matr_ftype = END_NODE as u8;

            // Calculate BEGINL -> MATL
            let mut trans = zero_transtable();
            trans[BEGIN_ST][MATL_ST] = nseq_f - gapcount[j];
            trans[BEGIN_ST][DEL_ST] = gapcount[j];
            prior.probify_transition_matrix(&mut trans, BEGINL_NODE, MATL_NODE);

            let beginl_trans_score = ((trans[BEGIN_ST][MATL_ST].ln() * (nseq_f - gapcount[j]))
                + (trans[BEGIN_ST][DEL_ST].ln() * gapcount[j]))
                * INTPRECISION;

            self.data[j][j].sc[maxmx_node::BEGINL] =
                self.data[j][j].sc[maxmx_node::MATL] + beginl_trans_score as i32;
            self.data[j][j].begl_ftype = maxmx_node::MATL as u8;

            // Calculate BEGINR -> MATL
            let mut trans = zero_transtable();
            trans[BEGIN_ST][DEL_ST] = gapcount[j];
            trans[BEGIN_ST][MATL_ST] = nseq_f - gapcount[j];
            prior.probify_transition_matrix(&mut trans, BEGINR_NODE, MATL_NODE);

            let beginr_trans_score = ((trans[BEGIN_ST][MATL_ST].ln() * (nseq_f - gapcount[j]))
                + (trans[BEGIN_ST][DEL_ST].ln() * gapcount[j]))
                * INTPRECISION;

            self.data[j][j].sc[maxmx_node::BEGINR] =
                self.data[j][j].sc[maxmx_node::MATL] + beginr_trans_score as i32;
            self.data[j][j].begr_i2 = j as i16;
            self.data[j][j].begr_ftype = maxmx_node::MATL as u8;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_maxmx_allocation() {
        let mmx = MaxMx::new(10);

        assert_eq!(mmx.alen, 10);
        assert_eq!(mmx.data.len(), 11); // 0..10 inclusive

        // Verify row sizes
        for j in 0..=10 {
            assert_eq!(mmx.data[j].len(), j + 2);
        }
    }

    #[test]
    fn test_maxmx_cell_initialization() {
        let mmx = MaxMx::new(5);

        // Check that cells are initialized with NEGINFINITY scores
        for j in 0..=5 {
            for i in 0..=(j + 1) {
                let cell = mmx.get(j, i);
                for y in 0..MAXMX_NODETYPES {
                    assert_eq!(cell.sc[y], NEGINFINITY);
                }
            }
        }
    }

    #[test]
    fn test_maxmx_traceback_pointers() {
        let mmx = MaxMx::new(3);

        // Check that traceback pointers point to self
        let cell = mmx.get(2, 1);
        assert_eq!(cell.matp_i2, 1);
        assert_eq!(cell.matp_j2, 2);
        assert_eq!(cell.matl_i2, 1);
        assert_eq!(cell.matr_j2, 2);
        assert_eq!(cell.begr_i2, 1);
        assert_eq!(cell.bifurc_mid, 1);
    }

    #[test]
    fn test_maxmx_get_mut() {
        let mut mmx = MaxMx::new(5);

        // Modify a cell
        mmx.get_mut(3, 2).sc[maxmx_node::MATP] = 100;

        assert_eq!(mmx.get(3, 2).sc[maxmx_node::MATP], 100);
    }
}
