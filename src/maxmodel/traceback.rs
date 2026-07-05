//! Traceback for MaxModelMaker algorithm
//!
//! This module implements the trace_maxmx function from maxmodelmaker.c
//! (lines 838-930) which constructs a consensus tree from the filled matrix.

use crate::maxmodel::matrix::MaxMx;
use crate::maxmodel::types::maxmx_node;
use crate::types::constants::*;

/// Master trace tree node for consensus structure
///
/// This represents a node in the consensus tree. The fields are used
/// slightly differently than in normal alignment tracebacks:
/// - emitl, emitr: 0..alen-1 column coordinates
/// - nodeidx: unused during construction, filled in later by NumberMasterTrace
/// - node_type: holds a node type (MATP_NODE, MATL_NODE, etc.)
#[derive(Debug, Clone)]
pub struct MasterTrace {
    /// Left position in alignment (0..alen-1) or -1 if nothing emitted
    pub emitl: i32,

    /// Right position in alignment (0..alen-1) or -1 if nothing emitted
    pub emitr: i32,

    /// Index of node (filled in by NumberMasterTrace)
    pub nodeidx: i32,

    /// Node type (MATP_NODE, MATL_NODE, MATR_NODE, BIFURC_NODE, etc.)
    pub node_type: usize,

    /// Left (or only) child
    pub nxtl: Option<Box<MasterTrace>>,

    /// Right child (BIFURC only)
    pub nxtr: Option<Box<MasterTrace>>,
}

impl MasterTrace {
    /// Create a new trace node
    pub fn new(emitl: i32, emitr: i32, node_type: usize) -> Self {
        Self {
            emitl,
            emitr,
            nodeidx: 0,
            node_type,
            nxtl: None,
            nxtr: None,
        }
    }

    /// Count total nodes in tree (iterative to avoid stack overflow)
    pub fn count_nodes(&self) -> usize {
        let mut count = 0;
        let mut stack: Vec<&MasterTrace> = vec![self];

        while let Some(node) = stack.pop() {
            count += 1;
            if let Some(ref left) = node.nxtl {
                stack.push(left);
            }
            if let Some(ref right) = node.nxtr {
                stack.push(right);
            }
        }

        count
    }

    /// Number nodes in pre-order traversal (iterative)
    /// Returns the total number of nodes
    pub fn number_nodes(&mut self, start_idx: i32) -> i32 {
        // We need to traverse and modify, so we use a work-stealing approach
        // First collect all node pointers, then assign indices
        let mut idx = start_idx;

        // Since we can't easily do in-place iteration with mutable borrows,
        // we use a simple recursive helper (tree depth is bounded by alignment length)
        fn number_helper(node: &mut MasterTrace, idx: &mut i32) {
            node.nodeidx = *idx;
            *idx += 1;

            if let Some(ref mut left) = node.nxtl {
                number_helper(left, idx);
            }
            if let Some(ref mut right) = node.nxtr {
                number_helper(right, idx);
            }
        }

        number_helper(self, &mut idx);
        idx
    }
}

impl Default for MasterTrace {
    fn default() -> Self {
        Self::new(0, 0, END_NODE)
    }
}

/// Traceback from filled MaxMx matrix to construct consensus tree
///
/// From maxmodelmaker.c trace_maxmx (lines 838-930)
///
/// Uses iterative approach to avoid stack overflow on deep trees.
///
/// The mmx scoring matrix traceback pointers are 1..alen coordinates.
/// They are converted to 0..alen-1 for the traceback tree.
///
/// ROOT alignment info is stored in mmx[alen][0] using MATP_NODE slots.
///
/// # Arguments
/// * `mmx` - Filled scoring matrix
///
/// # Returns
/// Root of the consensus traceback tree
pub fn trace_maxmx(mmx: &MaxMx) -> MasterTrace {
    let alen = mmx.alen;

    // Create root node spanning entire alignment
    let mut root = MasterTrace::new(0, (alen - 1) as i32, ROOT_NODE);

    // Get first segment from ROOT (stored in mmx[alen][0])
    let first_i = mmx.get(alen, 0).matp_i2 as i32 - 1;
    let first_j = mmx.get(alen, 0).matp_j2 as i32 - 1;
    let first_type = mmx.get(alen, 0).matp_ftype as usize;

    // Create placeholder for first child
    root.nxtl = Some(Box::new(MasterTrace::new(first_i, first_j, first_type)));

    // Use a worklist to iteratively build the tree
    // Each item contains: (path to parent, emitl, emitr, node_type, is_right_child)
    let mut worklist: Vec<(Vec<u8>, i32, i32, usize, bool)> = Vec::new();

    // Add first child to worklist
    worklist.push((vec![0], first_i, first_j, first_type, false));

    while let Some((path, emitl, emitr, node_type, _is_right)) = worklist.pop() {
        // Skip if off-diagonal or END
        if emitl > emitr || node_type == END_NODE {
            continue;
        }

        let i = (emitl + 1) as usize;
        let j = (emitr + 1) as usize;

        // Determine children based on node type
        let children = get_children(mmx, i, j, node_type);

        // Navigate to current node and set children
        let current = navigate_to_node(&mut root, &path);

        match children {
            TraceChildren::None => {
                // No children - leaf node
            }
            TraceChildren::Left(nxti, nxtj, nxt_type) => {
                current.nxtl = Some(Box::new(MasterTrace::new(nxti, nxtj, nxt_type)));
                let mut new_path = path.clone();
                new_path.push(0);
                worklist.push((new_path, nxti, nxtj, nxt_type, false));
            }
            TraceChildren::Both(li, lj, ltype, ri, rj, rtype) => {
                // BIFURC: add both children
                current.nxtl = Some(Box::new(MasterTrace::new(li, lj, ltype)));
                current.nxtr = Some(Box::new(MasterTrace::new(ri, rj, rtype)));

                // Add right first (so left is processed first due to stack order)
                let mut right_path = path.clone();
                right_path.push(1);
                worklist.push((right_path, ri, rj, rtype, true));

                let mut left_path = path.clone();
                left_path.push(0);
                worklist.push((left_path, li, lj, ltype, false));
            }
        }
    }

    root
}

/// Children information from traceback
enum TraceChildren {
    None,
    Left(i32, i32, usize),
    Both(i32, i32, usize, i32, i32, usize),
}

/// Get children for a node based on traceback pointers
fn get_children(mmx: &MaxMx, i: usize, j: usize, node_type: usize) -> TraceChildren {
    match node_type {
        MATP_NODE => {
            let nxti = mmx.get(j, i).matp_i2 as i32 - 1;
            let nxtj = mmx.get(j, i).matp_j2 as i32 - 1;
            let nxt_type = mmx.get(j, i).matp_ftype as usize;
            if nxti <= nxtj {
                TraceChildren::Left(nxti, nxtj, nxt_type)
            } else {
                TraceChildren::None
            }
        }

        MATL_NODE => {
            let nxti = mmx.get(j, i).matl_i2 as i32 - 1;
            let nxtj = (j as i32) - 1;
            let nxt_type = mmx.get(j, i).matl_ftype as usize;
            if nxti <= nxtj {
                TraceChildren::Left(nxti, nxtj, nxt_type)
            } else {
                TraceChildren::None
            }
        }

        MATR_NODE => {
            let nxti = (i as i32) - 1;
            let nxtj = mmx.get(j, i).matr_j2 as i32 - 1;
            let nxt_type = mmx.get(j, i).matr_ftype as usize;
            if nxti <= nxtj {
                TraceChildren::Left(nxti, nxtj, nxt_type)
            } else {
                TraceChildren::None
            }
        }

        BIFURC_NODE => {
            let mid = mmx.get(j, i).bifurc_mid as i32;
            let li = (i as i32) - 1;
            let lj = mid - 1;
            let ri = mid;
            let rj = (j as i32) - 1;
            TraceChildren::Both(li, lj, BEGINL_NODE, ri, rj, BEGINR_NODE)
        }

        BEGINL_NODE => {
            let nxti = (i as i32) - 1;
            let nxtj = (j as i32) - 1;
            let nxt_type = mmx.get(j, i).begl_ftype as usize;
            if nxti <= nxtj {
                TraceChildren::Left(nxti, nxtj, nxt_type)
            } else {
                TraceChildren::None
            }
        }

        BEGINR_NODE => {
            let nxti = mmx.get(j, i).begr_i2 as i32 - 1;
            let nxtj = (j as i32) - 1;
            let nxt_type = mmx.get(j, i).begr_ftype as usize;
            if nxti <= nxtj {
                TraceChildren::Left(nxti, nxtj, nxt_type)
            } else {
                TraceChildren::None
            }
        }

        _ => TraceChildren::None,
    }
}

/// Navigate to a node given path from root
fn navigate_to_node<'a>(root: &'a mut MasterTrace, path: &[u8]) -> &'a mut MasterTrace {
    let mut current = root;
    for &direction in path.iter() {
        current = if direction == 0 {
            current.nxtl.as_mut().unwrap().as_mut()
        } else {
            current.nxtr.as_mut().unwrap().as_mut()
        };
    }
    current
}

/// Convert maxmx_node index to CM node type constant
#[inline]
pub fn maxmx_to_cm_node(maxmx_type: u8) -> usize {
    match maxmx_type as usize {
        maxmx_node::MATP => MATP_NODE,
        maxmx_node::MATL => MATL_NODE,
        maxmx_node::MATR => MATR_NODE,
        maxmx_node::BIFURC => BIFURC_NODE,
        maxmx_node::BEGINL => BEGINL_NODE,
        maxmx_node::BEGINR => BEGINR_NODE,
        _ => END_NODE,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::maxmodel::emission::{is_gap, singlet_emissions, symbol_index, transpose_alignment};
    use crate::maxmodel::matrix::MaxMx;
    use crate::maxmodel::prior::Prior;
    use crate::maxmodel::recursion::recurse_maxmx;

    #[test]
    fn test_master_trace_new() {
        let trace = MasterTrace::new(5, 10, MATP_NODE);
        assert_eq!(trace.emitl, 5);
        assert_eq!(trace.emitr, 10);
        assert_eq!(trace.node_type, MATP_NODE);
        assert!(trace.nxtl.is_none());
        assert!(trace.nxtr.is_none());
    }

    #[test]
    fn test_master_trace_count_nodes() {
        // Single node
        let single = MasterTrace::new(0, 0, MATL_NODE);
        assert_eq!(single.count_nodes(), 1);

        // Tree with children
        let mut root = MasterTrace::new(0, 5, ROOT_NODE);
        root.nxtl = Some(Box::new(MasterTrace::new(1, 4, MATP_NODE)));
        assert_eq!(root.count_nodes(), 2);

        // Add bifurcation
        let mut bifurc = MasterTrace::new(1, 4, BIFURC_NODE);
        bifurc.nxtl = Some(Box::new(MasterTrace::new(1, 2, BEGINL_NODE)));
        bifurc.nxtr = Some(Box::new(MasterTrace::new(3, 4, BEGINR_NODE)));
        root.nxtl = Some(Box::new(bifurc));
        assert_eq!(root.count_nodes(), 4);
    }

    #[test]
    fn test_master_trace_number_nodes() {
        let mut root = MasterTrace::new(0, 5, ROOT_NODE);
        let mut child1 = MasterTrace::new(1, 4, MATP_NODE);
        child1.nxtl = Some(Box::new(MasterTrace::new(2, 3, MATL_NODE)));
        root.nxtl = Some(Box::new(child1));

        let next_idx = root.number_nodes(0);
        assert_eq!(next_idx, 3);
        assert_eq!(root.nodeidx, 0);
        assert_eq!(root.nxtl.as_ref().unwrap().nodeidx, 1);
        assert_eq!(
            root.nxtl
                .as_ref()
                .unwrap()
                .nxtl
                .as_ref()
                .unwrap()
                .nodeidx,
            2
        );
    }

    #[test]
    fn test_trace_maxmx() {
        // Build a simple alignment and run the full pipeline
        let seq1: &[u8] = b"ACGT";
        let seq2: &[u8] = b"ACGT";
        let aseqs: Vec<&[u8]> = vec![seq1, seq2];
        let weights = vec![1.0f32, 1.0];
        let prior = Prior::new();

        let aseqs_t = transpose_alignment(&aseqs, 4, is_gap, symbol_index);
        let (mscore, gapcount) = singlet_emissions(&aseqs_t, &weights, &prior);

        let mut mmx = MaxMx::new(4);
        mmx.initialize(2, &prior, &mscore, &gapcount);
        recurse_maxmx(&aseqs_t, &weights, &prior, &mscore, &gapcount, 0.5, &mut mmx);

        let trace = trace_maxmx(&mmx);

        // Root should exist and span the alignment
        assert_eq!(trace.node_type, ROOT_NODE);
        assert_eq!(trace.emitl, 0);
        assert_eq!(trace.emitr, 3); // alen - 1
        assert!(trace.nxtl.is_some());
    }

    #[test]
    fn test_maxmx_to_cm_node() {
        assert_eq!(maxmx_to_cm_node(maxmx_node::MATP as u8), MATP_NODE);
        assert_eq!(maxmx_to_cm_node(maxmx_node::MATL as u8), MATL_NODE);
        assert_eq!(maxmx_to_cm_node(maxmx_node::MATR as u8), MATR_NODE);
        assert_eq!(maxmx_to_cm_node(maxmx_node::BIFURC as u8), BIFURC_NODE);
        assert_eq!(maxmx_to_cm_node(maxmx_node::BEGINL as u8), BEGINL_NODE);
        assert_eq!(maxmx_to_cm_node(maxmx_node::BEGINR as u8), BEGINR_NODE);
    }
}
