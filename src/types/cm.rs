//! Covariance Model (CM) structures for tRNAscan-SE
//!
//! This module implements the core CM data structures from structs.h lines 76-99.
//! These represent the probabilistic model used for tRNA detection.

use crate::types::constants::*;

/// Node structure from structs.h lines 76-88
///
/// Each node in the CM tree contains:
/// - State type information
/// - Transition probabilities between states
/// - Emission probabilities for different state types
/// - Tree structure pointers (nxt, nxt2)
#[derive(Clone, Debug)]
pub struct Node {
    /// Node type (BIFURC_NODE, MATP_NODE, MATL_NODE, MATR_NODE, etc.)
    pub node_type: i32,

    /// Connection to left child node
    pub nxt: i32,

    /// Connection to right child node (BIFURC_NODE only)
    pub nxt2: i32,

    /// Transition probability matrix [from_state][to_state]
    /// Up to 49 transition probabilities (6x6 matrix, but sparsely populated)
    pub tmx: [[f64; STATETYPES]; STATETYPES],

    /// MATP (match-pair) emission probabilities [left_base][right_base]
    /// 4x4 matrix for paired base emissions
    pub mp_emit: [[f64; ALPHASIZE]; ALPHASIZE],

    /// INSL (insert-left) emission probabilities
    /// 4-element vector for single base emissions
    pub il_emit: [f64; ALPHASIZE],

    /// INSR (insert-right) emission probabilities
    /// 4-element vector for single base emissions
    pub ir_emit: [f64; ALPHASIZE],

    /// MATL (match-left) emission probabilities
    /// 4-element vector for single base emissions
    pub ml_emit: [f64; ALPHASIZE],

    /// MATR (match-right) emission probabilities
    /// 4-element vector for single base emissions
    pub mr_emit: [f64; ALPHASIZE],
}

impl Node {
    /// Create a new node with default (zero) probabilities
    pub fn new(node_type: i32) -> Self {
        Self {
            node_type,
            nxt: 0,
            nxt2: 0,
            tmx: [[0.0; STATETYPES]; STATETYPES],
            mp_emit: [[0.0; ALPHASIZE]; ALPHASIZE],
            il_emit: [0.0; ALPHASIZE],
            ir_emit: [0.0; ALPHASIZE],
            ml_emit: [0.0; ALPHASIZE],
            mr_emit: [0.0; ALPHASIZE],
        }
    }
}

impl Default for Node {
    fn default() -> Self {
        Self::new(0)
    }
}

/// Covariance Model (CM) structure from structs.h lines 92-99
///
/// A CM is a tree-based probabilistic model for RNA structure.
/// It consists of a collection of nodes arranged in a binary tree.
#[derive(Clone, Debug)]
pub struct CM {
    /// Number of nodes in the model
    pub nodes: usize,

    /// Array of nodes [0..nodes-1]
    pub nd: Vec<Node>,
}

impl CM {
    /// Create a new CM with the specified number of nodes
    pub fn new(nodes: usize) -> Self {
        Self {
            nodes,
            nd: vec![Node::default(); nodes],
        }
    }

    /// Get a reference to a node by index
    pub fn node(&self, idx: usize) -> Option<&Node> {
        self.nd.get(idx)
    }

    /// Get a mutable reference to a node by index
    pub fn node_mut(&mut self, idx: usize) -> Option<&mut Node> {
        self.nd.get_mut(idx)
    }
}

impl Default for CM {
    fn default() -> Self {
        Self::new(0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_node_creation() {
        let node = Node::new(MATP_NODE as i32);
        assert_eq!(node.node_type, MATP_NODE as i32);
        assert_eq!(node.nxt, 0);
        assert_eq!(node.nxt2, 0);

        // Verify all probabilities are initialized to zero
        for i in 0..STATETYPES {
            for j in 0..STATETYPES {
                assert_eq!(node.tmx[i][j], 0.0);
            }
        }

        for i in 0..ALPHASIZE {
            assert_eq!(node.il_emit[i], 0.0);
            assert_eq!(node.ir_emit[i], 0.0);
            assert_eq!(node.ml_emit[i], 0.0);
            assert_eq!(node.mr_emit[i], 0.0);

            for j in 0..ALPHASIZE {
                assert_eq!(node.mp_emit[i][j], 0.0);
            }
        }
    }

    #[test]
    fn test_cm_creation() {
        let cm = CM::new(10);
        assert_eq!(cm.nodes, 10);
        assert_eq!(cm.nd.len(), 10);

        // Verify all nodes are default-initialized
        for i in 0..10 {
            assert_eq!(cm.nd[i].node_type, 0);
        }
    }

    #[test]
    fn test_cm_node_access() {
        let mut cm = CM::new(5);

        // Test immutable access
        assert!(cm.node(0).is_some());
        assert!(cm.node(4).is_some());
        assert!(cm.node(5).is_none());

        // Test mutable access
        if let Some(node) = cm.node_mut(2) {
            node.node_type = MATP_NODE as i32;
        }
        assert_eq!(cm.nd[2].node_type, MATP_NODE as i32);
    }

    #[test]
    fn test_node_sizes() {
        // Verify array dimensions match constants
        let node = Node::default();
        assert_eq!(node.tmx.len(), STATETYPES);
        assert_eq!(node.tmx[0].len(), STATETYPES);
        assert_eq!(node.mp_emit.len(), ALPHASIZE);
        assert_eq!(node.mp_emit[0].len(), ALPHASIZE);
        assert_eq!(node.il_emit.len(), ALPHASIZE);
        assert_eq!(node.ir_emit.len(), ALPHASIZE);
        assert_eq!(node.ml_emit.len(), ALPHASIZE);
        assert_eq!(node.mr_emit.len(), ALPHASIZE);
    }
}
