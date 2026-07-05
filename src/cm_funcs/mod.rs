//! Covariance Model (CM) Functions
//!
//! This module implements functions for manipulating covariance models,
//! ported from the original C codebase (cm.c, maxmodelmaker.c).
//!
//! # Functions
//!
//! - `zero_cm` - Zero out all counts/probabilities
//! - `verify_cm` - Validate CM structure
//! - `normalize_cm` - Normalize probabilities to sum to 1.0
//! - `copy_cm` - Deep copy a CM
//! - `logoddsify_cm` - Convert probabilities to log-odds scores
//! - `probify_cm` - Convert counts to probabilities using priors
//! - `probify_emissions` - Convert emission counts to probabilities
//! - `probify_transitions` - Convert transition counts to probabilities
//! - `topofy_new_cm` - Set up CM topology from a master trace
//! - `topofy_from_guide` - Build CM topology from a guide tree
//! - `transmogrify` - Create individual trace from master trace
//! - `number_master_trace` - Assign sequential node numbers to master trace
//! - `trace_count` - Add weighted counts from a trace to CM
//! - `structure_stacks` - Parse WUSS notation to identify base-paired stems
//! - `cm_from_gparse` - Build CM from grammar parse

use crate::maxmodel::prior::Prior;
use crate::maxmodel::traceback::MasterTrace;
use crate::types::cm::{Node, CM};
use crate::types::constants::*;
use crate::types::trace::Trace;

// ============================================================================
// ZeroCM - Zero out all counts/probabilities
// ============================================================================

/// Zero out all probability parameters and counts in a CM
///
/// From cm.c CMZero() (lines 303-333)
///
/// Sets all emission counts, transition counts, and scores to zero.
/// This is typically called after allocating a new CM or before
/// re-training a model.
///
/// # Arguments
///
/// * `cm` - The CM to zero out
///
/// # Example
///
/// ```ignore
/// let mut cm = CM::new(10);
/// zero_cm(&mut cm);
/// ```
pub fn zero_cm(cm: &mut CM) {
    for node in cm.nd.iter_mut() {
        // Zero transition matrix
        for i in 0..STATETYPES {
            for j in 0..STATETYPES {
                node.tmx[i][j] = 0.0;
            }
        }

        // Zero emission probabilities
        for i in 0..ALPHASIZE {
            node.il_emit[i] = 0.0;
            node.ir_emit[i] = 0.0;
            node.ml_emit[i] = 0.0;
            node.mr_emit[i] = 0.0;
            for j in 0..ALPHASIZE {
                node.mp_emit[i][j] = 0.0;
            }
        }
    }
}

// ============================================================================
// VerifyCM - Validate CM structure
// ============================================================================

/// Validation error types for CM verification
#[derive(Debug, Clone, PartialEq)]
pub enum CmVerifyError {
    /// CM has no nodes
    EmptyModel,
    /// Invalid node connectivity (nxt out of range)
    InvalidConnectivity { node: usize, nxt: i32 },
    /// Emission probabilities don't sum to 1.0 (within tolerance)
    EmissionSumError {
        node: usize,
        state: &'static str,
        sum: f64,
    },
    /// Transition probabilities don't sum to 1.0 (within tolerance)
    TransitionSumError { node: usize, from_state: usize, sum: f64 },
    /// Invalid node type
    InvalidNodeType { node: usize, node_type: i32 },
    /// Missing root node
    MissingRoot,
}

impl std::fmt::Display for CmVerifyError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CmVerifyError::EmptyModel => write!(f, "CM has no nodes"),
            CmVerifyError::InvalidConnectivity { node, nxt } => {
                write!(f, "Node {} has invalid connectivity nxt={}", node, nxt)
            }
            CmVerifyError::EmissionSumError { node, state, sum } => {
                write!(
                    f,
                    "Node {} {} emissions sum to {} (expected 1.0)",
                    node, state, sum
                )
            }
            CmVerifyError::TransitionSumError {
                node,
                from_state,
                sum,
            } => {
                write!(
                    f,
                    "Node {} transitions from state {} sum to {} (expected 1.0)",
                    node, from_state, sum
                )
            }
            CmVerifyError::InvalidNodeType { node, node_type } => {
                write!(f, "Node {} has invalid type {}", node, node_type)
            }
            CmVerifyError::MissingRoot => write!(f, "CM has no root node"),
        }
    }
}

impl std::error::Error for CmVerifyError {}

/// Tolerance for probability sum checks
const PROB_TOLERANCE: f64 = 1e-6;

/// Verify that a CM is structurally valid
///
/// Checks:
/// - CM has at least one node
/// - Node connectivity is valid (nxt indices in range)
/// - Emission probabilities sum to approximately 1.0 for emitting nodes
/// - Transition probabilities sum to approximately 1.0 for non-leaf nodes
/// - Node types are valid
///
/// # Arguments
///
/// * `cm` - The CM to verify
///
/// # Returns
///
/// `Ok(())` if valid, `Err(CmVerifyError)` describing the first problem found
///
/// # Example
///
/// ```ignore
/// let cm = CM::new(10);
/// match verify_cm(&cm) {
///     Ok(()) => println!("CM is valid"),
///     Err(e) => println!("Invalid CM: {}", e),
/// }
/// ```
pub fn verify_cm(cm: &CM) -> Result<(), CmVerifyError> {
    // Check for empty model
    if cm.nodes == 0 || cm.nd.is_empty() {
        return Err(CmVerifyError::EmptyModel);
    }

    // Check that first node is a valid root type (ROOT_NODE = 6)
    let root_type = cm.nd[0].node_type;
    if root_type != ROOT_NODE as i32 && root_type != BEGINL_NODE as i32 {
        // Allow BEGINL_NODE for sub-CMs
        return Err(CmVerifyError::MissingRoot);
    }

    for (idx, node) in cm.nd.iter().enumerate() {
        // Check node type validity
        if node.node_type < 0 || node.node_type >= NODETYPES as i32 {
            return Err(CmVerifyError::InvalidNodeType {
                node: idx,
                node_type: node.node_type,
            });
        }

        // Check connectivity
        if node.nxt < -1 || node.nxt >= cm.nodes as i32 {
            return Err(CmVerifyError::InvalidConnectivity {
                node: idx,
                nxt: node.nxt,
            });
        }

        // For BIFURC_NODE, check nxt2 as well
        if node.node_type == BIFURC_NODE as i32 {
            if node.nxt2 < -1 || node.nxt2 >= cm.nodes as i32 {
                return Err(CmVerifyError::InvalidConnectivity {
                    node: idx,
                    nxt: node.nxt2,
                });
            }
        }

        // Check emission probability sums for emitting node types
        match node.node_type as usize {
            MATP_NODE => {
                // Check pair emissions
                let mut sum = 0.0;
                for i in 0..ALPHASIZE {
                    for j in 0..ALPHASIZE {
                        sum += node.mp_emit[i][j];
                    }
                }
                if sum > PROB_TOLERANCE && (sum - 1.0).abs() > PROB_TOLERANCE {
                    return Err(CmVerifyError::EmissionSumError {
                        node: idx,
                        state: "MATP",
                        sum,
                    });
                }
            }
            MATL_NODE => {
                // Check left match emissions
                let sum: f64 = node.ml_emit.iter().sum();
                if sum > PROB_TOLERANCE && (sum - 1.0).abs() > PROB_TOLERANCE {
                    return Err(CmVerifyError::EmissionSumError {
                        node: idx,
                        state: "MATL",
                        sum,
                    });
                }
            }
            MATR_NODE => {
                // Check right match emissions
                let sum: f64 = node.mr_emit.iter().sum();
                if sum > PROB_TOLERANCE && (sum - 1.0).abs() > PROB_TOLERANCE {
                    return Err(CmVerifyError::EmissionSumError {
                        node: idx,
                        state: "MATR",
                        sum,
                    });
                }
            }
            _ => {
                // BIFURC_NODE, BEGINL_NODE, BEGINR_NODE, ROOT_NODE don't emit
            }
        }

        // Check transition probability sums
        // Only check for nodes that have outgoing transitions
        if node.node_type != END_NODE as i32 {
            for from_state in 0..STATETYPES {
                let sum: f64 = node.tmx[from_state].iter().sum();
                // Only check if there are any transitions from this state
                if sum > PROB_TOLERANCE && (sum - 1.0).abs() > PROB_TOLERANCE {
                    return Err(CmVerifyError::TransitionSumError {
                        node: idx,
                        from_state,
                        sum,
                    });
                }
            }
        }
    }

    Ok(())
}

/// Simple boolean verification (compatible with original C interface)
///
/// Returns `true` if the CM is valid, `false` otherwise.
pub fn verify_cm_bool(cm: &CM) -> bool {
    verify_cm(cm).is_ok()
}

// ============================================================================
// NormalizeCM - Normalize probabilities
// ============================================================================

/// Normalize all probability distributions in a CM to sum to 1.0
///
/// From cm.c CMRenormalize() (lines 335-361)
///
/// Normalizes:
/// - Emission probabilities for each emitting state type
/// - Transition probabilities from each state
///
/// # Arguments
///
/// * `cm` - The CM to normalize
///
/// # Example
///
/// ```ignore
/// let mut cm = create_cm_with_counts();
/// normalize_cm(&mut cm);
/// ```
pub fn normalize_cm(cm: &mut CM) {
    for node in cm.nd.iter_mut() {
        // Normalize based on node type
        match node.node_type as usize {
            MATP_NODE => {
                // Normalize pair emissions
                let mut sum = 0.0;
                for i in 0..ALPHASIZE {
                    for j in 0..ALPHASIZE {
                        sum += node.mp_emit[i][j];
                    }
                }
                if sum > 0.0 {
                    for i in 0..ALPHASIZE {
                        for j in 0..ALPHASIZE {
                            node.mp_emit[i][j] /= sum;
                        }
                    }
                }

                // Also normalize insertion emissions for MATP nodes
                normalize_singlet(&mut node.il_emit);
                normalize_singlet(&mut node.ir_emit);
            }
            MATL_NODE => {
                normalize_singlet(&mut node.ml_emit);
                normalize_singlet(&mut node.il_emit);
            }
            MATR_NODE => {
                normalize_singlet(&mut node.mr_emit);
                normalize_singlet(&mut node.ir_emit);
            }
            _ => {
                // BIFURC, BEGINL, BEGINR, ROOT don't emit
            }
        }

        // Normalize transitions for all non-end nodes
        if node.node_type != END_NODE as i32 {
            for from_state in 0..STATETYPES {
                let sum: f64 = node.tmx[from_state].iter().sum();
                if sum > 0.0 {
                    for to_state in 0..STATETYPES {
                        node.tmx[from_state][to_state] /= sum;
                    }
                }
            }
        }
    }
}

/// Helper: normalize a singlet emission vector
fn normalize_singlet(emit: &mut [f64; ALPHASIZE]) {
    let sum: f64 = emit.iter().sum();
    if sum > 0.0 {
        for e in emit.iter_mut() {
            *e /= sum;
        }
    }
}

// ============================================================================
// CopyCM - Deep copy
// ============================================================================

/// Create a deep copy of a CM
///
/// Copies all fields including nodes, emission/transition probabilities.
///
/// # Arguments
///
/// * `src` - Source CM to copy
///
/// # Returns
///
/// A new CM that is an exact copy of the source
///
/// # Example
///
/// ```ignore
/// let original = CM::new(10);
/// let copy = copy_cm(&original);
/// assert_eq!(copy.nodes, original.nodes);
/// ```
pub fn copy_cm(src: &CM) -> CM {
    CM {
        nodes: src.nodes,
        nd: src.nd.clone(), // Vec<Node> implements Clone
    }
}

// ============================================================================
// LogoddsifyCM - Convert probabilities to log-odds
// ============================================================================

/// Convert CM probabilities to log-odds scores
///
/// From cm.c CMLogoddsify() (lines 774-950)
///
/// Converts:
/// - Emission probabilities to log-odds vs null model
/// - Transition probabilities to log probabilities
///
/// Log-odds score = log2(prob / null_prob)
///
/// # Arguments
///
/// * `cm` - The CM to convert (modified in place)
/// * `null` - Null model probabilities (uniform if None)
///
/// # Example
///
/// ```ignore
/// let mut cm = create_normalized_cm();
/// let null = [0.25, 0.25, 0.25, 0.25]; // Uniform null model
/// logoddsify_cm(&mut cm, Some(&null));
/// ```
pub fn logoddsify_cm(cm: &mut CM, null: Option<&[f64; ALPHASIZE]>) {
    // Default null model: uniform probabilities
    let default_null = [0.25f64; ALPHASIZE];
    let null = null.unwrap_or(&default_null);

    for node in cm.nd.iter_mut() {
        match node.node_type as usize {
            MATP_NODE => {
                // Convert pair emissions to log-odds
                for i in 0..ALPHASIZE {
                    for j in 0..ALPHASIZE {
                        let prob = node.mp_emit[i][j];
                        let null_prob = null[i] * null[j];
                        node.mp_emit[i][j] = prob_to_log_odds(prob, null_prob);
                    }
                }

                // Convert insertion emissions
                logoddsify_singlet(&mut node.il_emit, null);
                logoddsify_singlet(&mut node.ir_emit, null);
            }
            MATL_NODE => {
                logoddsify_singlet(&mut node.ml_emit, null);
                logoddsify_singlet(&mut node.il_emit, null);
            }
            MATR_NODE => {
                logoddsify_singlet(&mut node.mr_emit, null);
                logoddsify_singlet(&mut node.ir_emit, null);
            }
            _ => {}
        }

        // Convert transitions to log probabilities
        if node.node_type != END_NODE as i32 {
            for from_state in 0..STATETYPES {
                for to_state in 0..STATETYPES {
                    let prob = node.tmx[from_state][to_state];
                    node.tmx[from_state][to_state] = prob_to_log(prob);
                }
            }
        }
    }
}

/// Convert probability to log-odds score (log2(prob/null))
#[inline]
fn prob_to_log_odds(prob: f64, null_prob: f64) -> f64 {
    if prob > 0.0 && null_prob > 0.0 {
        (prob / null_prob).log2()
    } else if prob <= 0.0 {
        NEGINFINITY as f64 / INTPRECISION
    } else {
        0.0
    }
}

/// Convert probability to log probability (log2(prob))
#[inline]
fn prob_to_log(prob: f64) -> f64 {
    if prob > 0.0 {
        prob.log2()
    } else {
        NEGINFINITY as f64 / INTPRECISION
    }
}

/// Helper: convert singlet emission vector to log-odds
fn logoddsify_singlet(emit: &mut [f64; ALPHASIZE], null: &[f64; ALPHASIZE]) {
    for i in 0..ALPHASIZE {
        emit[i] = prob_to_log_odds(emit[i], null[i]);
    }
}

// ============================================================================
// ProbifyCM - Convert counts to probabilities
// ============================================================================

/// Convert counts in a CM to probabilities using Dirichlet priors
///
/// From probify.c Probify() and related functions
///
/// Applies Bayesian regularization using the provided prior distributions
/// to convert raw counts to smoothed probability estimates.
///
/// # Arguments
///
/// * `cm` - The CM with counts (modified in place to probabilities)
/// * `prior` - Prior probability distributions for regularization
///
/// # Example
///
/// ```ignore
/// let mut cm = create_cm_with_counts();
/// let prior = Prior::new();
/// probify_cm(&mut cm, &prior);
/// ```
pub fn probify_cm(cm: &mut CM, prior: &Prior) {
    probify_emissions(cm, prior);
    probify_transitions(cm, prior);
}

// ============================================================================
// ProbifyEmissions - Emission probability conversion
// ============================================================================

/// Convert emission counts to probabilities with Dirichlet priors
///
/// From probify.c ProbifyEmissions
///
/// # Arguments
///
/// * `cm` - The CM to convert (modified in place)
/// * `prior` - Prior probability distributions
pub fn probify_emissions(cm: &mut CM, prior: &Prior) {
    for (idx, node) in cm.nd.iter_mut().enumerate() {
        match node.node_type as usize {
            MATP_NODE => {
                // Apply pair emission prior
                prior.probify_pair_emission(&mut node.mp_emit);

                // Also probify insertion emissions
                probify_singlet_emission(&mut node.il_emit, prior, INSL_ST);
                probify_singlet_emission(&mut node.ir_emit, prior, INSR_ST);
            }
            MATL_NODE => {
                probify_singlet_emission(&mut node.ml_emit, prior, MATL_ST);
                probify_singlet_emission(&mut node.il_emit, prior, INSL_ST);
            }
            MATR_NODE => {
                probify_singlet_emission(&mut node.mr_emit, prior, MATR_ST);
                probify_singlet_emission(&mut node.ir_emit, prior, INSR_ST);
            }
            _ => {
                // Non-emitting nodes (BIFURC, BEGINL, BEGINR, ROOT) need no emission probify
                let _ = idx; // Silence unused variable warning
            }
        }
    }
}

/// Helper: apply prior to singlet emission vector
fn probify_singlet_emission(emit: &mut [f64; ALPHASIZE], prior: &Prior, state_type: usize) {
    // Get the appropriate prior and alpha for this state type
    let (em_prior, alpha) = match state_type {
        MATL_ST => (&prior.matl_prior, prior.emalpha[MATL_ST]),
        MATR_ST => (&prior.matr_prior, prior.emalpha[MATR_ST]),
        INSL_ST => (&prior.insl_prior, prior.emalpha[INSL_ST]),
        INSR_ST => (&prior.insr_prior, prior.emalpha[INSR_ST]),
        _ => (&prior.matl_prior, prior.emalpha[MATL_ST]),
    };

    // Add prior pseudocounts and normalize
    let mut denom = 0.0;
    for x in 0..ALPHASIZE {
        emit[x] = emit[x] + alpha * em_prior[x];
        denom += emit[x];
    }
    if denom > 0.0 {
        for x in 0..ALPHASIZE {
            emit[x] /= denom;
        }
    }
}

// ============================================================================
// ProbifyTransitions - Transition probability conversion
// ============================================================================

/// Convert transition counts to probabilities with Dirichlet priors
///
/// From probify.c ProbifyTransitions
///
/// # Arguments
///
/// * `cm` - The CM to convert (modified in place)
/// * `prior` - Prior probability distributions
pub fn probify_transitions(cm: &mut CM, prior: &Prior) {
    // First pass: collect node types and nxt indices
    let node_info: Vec<(i32, i32)> = cm.nd.iter()
        .map(|n| (n.node_type, n.nxt))
        .collect();

    for idx in 0..cm.nodes {
        let (node_type, nxt) = node_info[idx];

        // Skip END nodes (no outgoing transitions)
        if node_type == END_NODE as i32 {
            continue;
        }

        // Determine the next node type for transition prior lookup
        let from_node_type = node_type as usize;
        let to_node_type = if nxt >= 0 && (nxt as usize) < cm.nodes {
            node_info[nxt as usize].0 as usize
        } else {
            END_NODE
        };

        // Apply prior to transition matrix
        prior.probify_transition_matrix(&mut cm.nd[idx].tmx, from_node_type, to_node_type);
    }
}

// ============================================================================
// TopofyNewCM - Set node topology from trace
// ============================================================================

/// Set up CM node connectivity from a master trace tree
///
/// From maxmodelmaker.c TopofyNewCM
///
/// Assigns node types and sets up the nxt/nxt2 connectivity pointers
/// based on the structure of the master trace.
///
/// # Arguments
///
/// * `cm` - The CM to set up (should have nodes allocated)
/// * `trace` - The master trace tree defining the structure
///
/// # Example
///
/// ```ignore
/// let trace = create_master_trace();
/// let mut cm = CM::new(trace.count_nodes());
/// topofy_new_cm(&mut cm, &trace);
/// ```
pub fn topofy_new_cm(cm: &mut CM, trace: &MasterTrace) {
    // First pass: assign node types from trace
    let mut stack: Vec<(&MasterTrace, Option<usize>)> = vec![(trace, None)];

    while let Some((node, parent_idx)) = stack.pop() {
        let idx = node.nodeidx as usize;
        if idx >= cm.nodes {
            continue;
        }

        // Set node type from trace
        cm.nd[idx].node_type = node.node_type as i32;

        // Set parent's nxt pointer to this node
        if let Some(pidx) = parent_idx {
            if cm.nd[pidx].nxt < 0 {
                cm.nd[pidx].nxt = idx as i32;
            } else if cm.nd[pidx].node_type == BIFURC_NODE as i32 {
                cm.nd[pidx].nxt2 = idx as i32;
            }
        }

        // Add children to stack (right first so left is processed first)
        if let Some(ref right) = node.nxtr {
            stack.push((right, Some(idx)));
        }
        if let Some(ref left) = node.nxtl {
            stack.push((left, Some(idx)));
        }
    }
}

// ============================================================================
// TopofyFromGuide - Set topology from guide tree
// ============================================================================

/// Guide tree node for CM construction
#[derive(Debug, Clone)]
pub struct GuideNode {
    /// Node type (MATP_NODE, MATL_NODE, etc.)
    pub node_type: usize,
    /// Left child index (-1 if none)
    pub left: i32,
    /// Right child index (BIFURC only, -1 if none)
    pub right: i32,
    /// Emit positions in alignment (left, right)
    pub emit_pos: (i32, i32),
}

/// Guide tree for CM construction
#[derive(Debug)]
pub struct GuideTree {
    /// Nodes in the guide tree
    pub nodes: Vec<GuideNode>,
}

impl GuideTree {
    /// Create a new empty guide tree
    pub fn new() -> Self {
        GuideTree { nodes: Vec::new() }
    }

    /// Add a node to the guide tree
    pub fn add_node(&mut self, node: GuideNode) -> usize {
        let idx = self.nodes.len();
        self.nodes.push(node);
        idx
    }
}

impl Default for GuideTree {
    fn default() -> Self {
        Self::new()
    }
}

/// Build CM topology from a guide tree
///
/// Sets up node types and connectivity based on the guide tree structure.
///
/// # Arguments
///
/// * `cm` - The CM to set up
/// * `guide` - The guide tree defining the structure
pub fn topofy_from_guide(cm: &mut CM, guide: &GuideTree) {
    // Ensure CM has enough nodes
    if cm.nodes < guide.nodes.len() {
        // Resize CM
        cm.nd.resize(guide.nodes.len(), Node::default());
        cm.nodes = guide.nodes.len();
    }

    for (idx, gnode) in guide.nodes.iter().enumerate() {
        cm.nd[idx].node_type = gnode.node_type as i32;
        cm.nd[idx].nxt = gnode.left;
        cm.nd[idx].nxt2 = gnode.right;
    }
}

// ============================================================================
// Transmogrify - Create individual trace from master
// ============================================================================

/// Trace pool for efficient trace allocation
#[derive(Debug, Default)]
pub struct TracePool {
    /// Pre-allocated trace nodes
    traces: Vec<Trace>,
}

impl TracePool {
    /// Create a new trace pool
    pub fn new() -> Self {
        Self::default()
    }

    /// Get next available trace (or create new one)
    pub fn alloc(&mut self) -> Trace {
        self.traces.pop().unwrap_or_else(|| Trace::new(0, 0, 0, 0))
    }

    /// Return a trace to the pool
    pub fn free(&mut self, trace: Trace) {
        self.traces.push(trace);
    }
}

/// Create a sequence-specific trace from a master trace
///
/// From maxmodelmaker.c Transmogrify (lines 400-500)
///
/// "Transmogrify" converts the consensus master trace into a trace
/// specific to a particular sequence, by:
/// - Mapping alignment columns to sequence positions
/// - Handling insertions and deletions
///
/// # Arguments
///
/// * `mtr` - The master trace (consensus structure)
/// * `seq` - The sequence (may include gaps in alignment form)
///
/// # Returns
///
/// A tuple of (sequence-specific Trace, TracePool)
pub fn transmogrify(mtr: &MasterTrace, seq: &[u8]) -> (Trace, TracePool) {
    let pool = TracePool::new();

    // Build a mapping from alignment position to sequence position
    let mut aseq_to_seq: Vec<i32> = Vec::with_capacity(seq.len() + 1);
    let mut seqpos: i32 = 0;

    for (i, &c) in seq.iter().enumerate() {
        if !is_gap_char(c) {
            seqpos += 1;
        }
        aseq_to_seq.push(seqpos);
        let _ = i; // Mark as used
    }

    // Recursively transmogrify the trace tree
    let trace = transmogrify_node(mtr, &aseq_to_seq, seq);

    (trace, pool)
}

/// Check if a character is a gap
#[inline]
fn is_gap_char(c: u8) -> bool {
    matches!(c, b'-' | b'.' | b'_' | b'~')
}

/// Recursively convert a master trace node to a sequence-specific trace
fn transmogrify_node(mtr: &MasterTrace, aseq_to_seq: &[i32], seq: &[u8]) -> Trace {
    // Map alignment positions to sequence positions
    let emitl = if mtr.emitl >= 0 && (mtr.emitl as usize) < seq.len() {
        if !is_gap_char(seq[mtr.emitl as usize]) {
            aseq_to_seq[mtr.emitl as usize]
        } else {
            0 // Gap at this position
        }
    } else {
        0
    };

    let emitr = if mtr.emitr >= 0 && (mtr.emitr as usize) < seq.len() {
        if !is_gap_char(seq[mtr.emitr as usize]) {
            aseq_to_seq[mtr.emitr as usize]
        } else {
            0
        }
    } else {
        0
    };

    // Convert node type to trace type (unique state flag)
    let trace_type = node_type_to_trace_type(mtr.node_type);

    let mut trace = Trace::new(emitl, emitr, mtr.nodeidx, trace_type);

    // Recursively process children
    if let Some(ref left) = mtr.nxtl {
        trace.nxtl = Some(Box::new(transmogrify_node(left, aseq_to_seq, seq)));
    }

    if let Some(ref right) = mtr.nxtr {
        trace.nxtr = Some(Box::new(transmogrify_node(right, aseq_to_seq, seq)));
    }

    trace
}

/// Convert node type to unique state flag for trace
fn node_type_to_trace_type(node_type: usize) -> u32 {
    match node_type {
        MATP_NODE => U_MATP_ST,
        MATL_NODE => U_MATL_ST,
        MATR_NODE => U_MATR_ST,
        BIFURC_NODE => U_BIFURC_ST,
        BEGINL_NODE | BEGINR_NODE => U_BEGIN_ST,
        ROOT_NODE => U_BEGIN_ST,
        _ => U_END_ST,
    }
}

// ============================================================================
// NumberMasterTrace - Assign node numbers
// ============================================================================

/// Assign sequential node indices to a master trace tree
///
/// From maxmodelmaker.c NumberMasterTrace
///
/// Traverses the tree in pre-order and assigns sequential indices
/// starting from 0.
///
/// # Arguments
///
/// * `mtr` - The master trace to number (modified in place)
///
/// # Returns
///
/// Total number of nodes in the tree
///
/// # Example
///
/// ```ignore
/// let mut mtr = create_master_trace();
/// let count = number_master_trace(&mut mtr);
/// println!("Tree has {} nodes", count);
/// ```
pub fn number_master_trace(mtr: &mut MasterTrace) -> usize {
    // Use the existing method on MasterTrace
    let next_idx = mtr.number_nodes(0);
    next_idx as usize
}

// ============================================================================
// TraceCount - Count emissions/transitions
// ============================================================================

/// Add weighted counts from a trace to a CM
///
/// From maxmodelmaker.c TraceCount (lines 550-700)
///
/// Walks through a trace, adding emission and transition counts
/// to the CM. The counts are weighted by the provided weight.
///
/// # Arguments
///
/// * `cm` - The CM to add counts to
/// * `seq` - The sequence (used for emission counting)
/// * `weight` - Weight for this sequence (usually 1.0)
/// * `tr` - The trace to count from
///
/// # Returns
///
/// `true` if counting succeeded, `false` on error
pub fn trace_count(cm: &mut CM, seq: &[u8], weight: f32, tr: &Trace) -> bool {
    trace_count_recursive(cm, seq, weight, tr)
}

/// Recursive helper for trace_count
fn trace_count_recursive(cm: &mut CM, seq: &[u8], weight: f32, tr: &Trace) -> bool {
    let nodeidx = tr.nodeidx as usize;
    if nodeidx >= cm.nodes {
        return false;
    }

    let node = &mut cm.nd[nodeidx];
    let weight_f64 = weight as f64;

    // Count emissions based on state type
    if tr.is_type(U_MATP_ST) {
        // Pair emission
        if tr.emitl > 0 && tr.emitr > 0 {
            let left_pos = (tr.emitl - 1) as usize;
            let right_pos = (tr.emitr - 1) as usize;
            if left_pos < seq.len() && right_pos < seq.len() {
                if let (Some(li), Some(ri)) = (base_to_index(seq[left_pos]), base_to_index(seq[right_pos])) {
                    node.mp_emit[li][ri] += weight_f64;
                }
            }
        }
    } else if tr.is_type(U_MATL_ST) {
        // Left singlet emission
        if tr.emitl > 0 {
            let pos = (tr.emitl - 1) as usize;
            if pos < seq.len() {
                if let Some(bi) = base_to_index(seq[pos]) {
                    node.ml_emit[bi] += weight_f64;
                }
            }
        }
    } else if tr.is_type(U_MATR_ST) {
        // Right singlet emission
        if tr.emitr > 0 {
            let pos = (tr.emitr - 1) as usize;
            if pos < seq.len() {
                if let Some(bi) = base_to_index(seq[pos]) {
                    node.mr_emit[bi] += weight_f64;
                }
            }
        }
    } else if tr.is_type(U_INSL_ST) {
        // Left insertion emission
        if tr.emitl > 0 {
            let pos = (tr.emitl - 1) as usize;
            if pos < seq.len() {
                if let Some(bi) = base_to_index(seq[pos]) {
                    node.il_emit[bi] += weight_f64;
                }
            }
        }
    } else if tr.is_type(U_INSR_ST) {
        // Right insertion emission
        if tr.emitr > 0 {
            let pos = (tr.emitr - 1) as usize;
            if pos < seq.len() {
                if let Some(bi) = base_to_index(seq[pos]) {
                    node.ir_emit[bi] += weight_f64;
                }
            }
        }
    }

    // Count transitions
    // Transition from this node to child node
    if let Some(ref left) = tr.nxtl {
        let child_idx = left.nodeidx as usize;
        if child_idx < cm.nodes {
            // Determine state indices for transition
            let from_state = trace_type_to_state_index(tr.trace_type);
            let to_state = trace_type_to_state_index(left.trace_type);
            if from_state < STATETYPES && to_state < STATETYPES {
                node.tmx[from_state][to_state] += weight_f64;
            }
        }

        // Recurse to left child
        if !trace_count_recursive(cm, seq, weight, left) {
            return false;
        }
    }

    // For bifurcation, also process right child
    if let Some(ref right) = tr.nxtr {
        if !trace_count_recursive(cm, seq, weight, right) {
            return false;
        }
    }

    true
}

/// Convert base character to index (A=0, C=1, G=2, T/U=3)
#[inline]
fn base_to_index(base: u8) -> Option<usize> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' | b'U' => Some(3),
        _ => None,
    }
}

/// Convert trace type flag to state index
fn trace_type_to_state_index(trace_type: u32) -> usize {
    if trace_type & U_DEL_ST != 0 { DEL_ST }
    else if trace_type & U_MATP_ST != 0 { MATP_ST }
    else if trace_type & U_MATL_ST != 0 { MATL_ST }
    else if trace_type & U_MATR_ST != 0 { MATR_ST }
    else if trace_type & U_INSL_ST != 0 { INSL_ST }
    else if trace_type & U_INSR_ST != 0 { INSR_ST }
    else { DEL_ST }
}

// ============================================================================
// StructureStacks - Identify stems in structure
// ============================================================================

/// Represents a base-paired stem (stack) in secondary structure
#[derive(Debug, Clone, PartialEq)]
pub struct Stack {
    /// Left position of the stem (5' end)
    pub left_start: usize,
    /// Right position of the stem (3' end)
    pub right_start: usize,
    /// Number of base pairs in the stem
    pub length: usize,
}

/// Parse WUSS notation to identify base-paired stems
///
/// WUSS (Washington University Secondary Structure) notation uses:
/// - `(` and `)` for base pairs
/// - `<` and `>` for pseudoknot base pairs
/// - `.` for unpaired bases
/// - `:` for gaps
///
/// This function identifies contiguous stretches of base pairs (stems).
///
/// # Arguments
///
/// * `wuss` - The WUSS secondary structure string
///
/// # Returns
///
/// Vector of Stack structures representing the stems
///
/// # Example
///
/// ```ignore
/// let wuss = "(((....)))";
/// let stacks = structure_stacks(wuss);
/// assert_eq!(stacks.len(), 1);
/// assert_eq!(stacks[0].length, 3);
/// ```
pub fn structure_stacks(wuss: &str) -> Vec<Stack> {
    let mut stacks = Vec::new();
    let chars: Vec<char> = wuss.chars().collect();
    let n = chars.len();

    // First, compute the pairing partners
    let mut partner: Vec<Option<usize>> = vec![None; n];
    let mut stack: Vec<usize> = Vec::new();

    for (i, &c) in chars.iter().enumerate() {
        match c {
            '(' | '<' | '[' | '{' => {
                stack.push(i);
            }
            ')' | '>' | ']' | '}' => {
                if let Some(j) = stack.pop() {
                    partner[i] = Some(j);
                    partner[j] = Some(i);
                }
            }
            _ => {}
        }
    }

    // Now identify stems (contiguous base pairs)
    let mut used = vec![false; n];

    for i in 0..n {
        if used[i] || partner[i].is_none() {
            continue;
        }

        let j = partner[i].unwrap();
        if j <= i {
            continue; // Only process left partners
        }

        // Found a base pair, extend to find the full stem
        let mut left = i;
        let mut right = j;
        let mut length = 1;

        // Extend inward while maintaining consecutive base pairs
        while left + 1 < right - 1 {
            let next_left = left + 1;
            let next_right = right - 1;

            if partner[next_left] == Some(next_right) {
                left = next_left;
                right = next_right;
                length += 1;
                used[next_left] = true;
                used[next_right] = true;
            } else {
                break;
            }
        }

        used[i] = true;
        used[j] = true;

        stacks.push(Stack {
            left_start: i,
            right_start: j,
            length,
        });
    }

    stacks
}

// ============================================================================
// CMFromGParse - Build CM from grammar parse
// ============================================================================

/// Grammar parse node for CM construction
#[derive(Debug, Clone)]
pub struct GrammarParseNode {
    /// Node type in the grammar
    pub node_type: usize,
    /// Left position in sequence
    pub left_pos: i32,
    /// Right position in sequence
    pub right_pos: i32,
    /// Children indices
    pub children: Vec<usize>,
}

/// Grammar parse tree for CM construction
#[derive(Debug)]
pub struct GrammarParse {
    /// Nodes in the parse tree
    pub nodes: Vec<GrammarParseNode>,
    /// Root node index
    pub root: usize,
}

impl GrammarParse {
    /// Create a new empty grammar parse
    pub fn new() -> Self {
        GrammarParse {
            nodes: Vec::new(),
            root: 0,
        }
    }

    /// Add a node to the parse
    pub fn add_node(&mut self, node: GrammarParseNode) -> usize {
        let idx = self.nodes.len();
        self.nodes.push(node);
        idx
    }
}

impl Default for GrammarParse {
    fn default() -> Self {
        Self::new()
    }
}

/// Build a CM from a grammar parse tree
///
/// Constructs a covariance model based on the structure defined
/// in the grammar parse.
///
/// # Arguments
///
/// * `gparse` - The grammar parse tree
///
/// # Returns
///
/// A new CM with structure matching the parse
pub fn cm_from_gparse(gparse: &GrammarParse) -> CM {
    let num_nodes = gparse.nodes.len();
    let mut cm = CM::new(num_nodes);

    // Build CM structure from parse tree
    for (idx, pnode) in gparse.nodes.iter().enumerate() {
        cm.nd[idx].node_type = pnode.node_type as i32;

        // Set up connectivity
        if !pnode.children.is_empty() {
            cm.nd[idx].nxt = pnode.children[0] as i32;
        }

        if pnode.children.len() > 1 && pnode.node_type == BIFURC_NODE {
            cm.nd[idx].nxt2 = pnode.children[1] as i32;
        }
    }

    cm
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_zero_cm() {
        let mut cm = CM::new(3);
        // Set some non-zero values
        cm.nd[0].tmx[0][0] = 1.0;
        cm.nd[1].mp_emit[0][1] = 0.5;
        cm.nd[2].ml_emit[2] = 0.25;

        zero_cm(&mut cm);

        // Verify all are zero
        for node in &cm.nd {
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
            }
        }
    }

    #[test]
    fn test_verify_cm_empty() {
        let cm = CM::new(0);
        assert!(matches!(verify_cm(&cm), Err(CmVerifyError::EmptyModel)));
    }

    #[test]
    fn test_verify_cm_valid_root() {
        let mut cm = CM::new(2);
        cm.nd[0].node_type = ROOT_NODE as i32;
        cm.nd[0].nxt = 1;
        cm.nd[1].node_type = END_NODE as i32;
        cm.nd[1].nxt = -1;

        // Should pass basic checks (empty probabilities are OK for a new CM)
        // Actually, empty probs won't pass sum checks, so let's just test structure
        assert!(verify_cm(&cm).is_ok());
    }

    #[test]
    fn test_normalize_singlet() {
        let mut emit = [1.0, 2.0, 3.0, 4.0];
        normalize_singlet(&mut emit);

        let sum: f64 = emit.iter().sum();
        assert!((sum - 1.0).abs() < 1e-10);
        assert!((emit[0] - 0.1).abs() < 1e-10);
        assert!((emit[3] - 0.4).abs() < 1e-10);
    }

    #[test]
    fn test_copy_cm() {
        let mut original = CM::new(2);
        original.nd[0].node_type = ROOT_NODE as i32;
        original.nd[0].nxt = 1;
        original.nd[1].mp_emit[0][0] = 0.5;

        let copy = copy_cm(&original);

        assert_eq!(copy.nodes, original.nodes);
        assert_eq!(copy.nd[0].node_type, original.nd[0].node_type);
        assert_eq!(copy.nd[0].nxt, original.nd[0].nxt);
        assert_eq!(copy.nd[1].mp_emit[0][0], original.nd[1].mp_emit[0][0]);
    }

    #[test]
    fn test_prob_to_log_odds() {
        // Equal prob and null -> log odds = 0
        assert!((prob_to_log_odds(0.25, 0.25) - 0.0).abs() < 1e-10);

        // Double null -> log odds = 1
        assert!((prob_to_log_odds(0.5, 0.25) - 1.0).abs() < 1e-10);

        // Half null -> log odds = -1
        assert!((prob_to_log_odds(0.125, 0.25) - (-1.0)).abs() < 1e-10);

        // Zero prob -> large negative
        let lo = prob_to_log_odds(0.0, 0.25);
        assert!(lo < -900.0);
    }

    #[test]
    fn test_structure_stacks_simple() {
        let wuss = "(((....)))";
        let stacks = structure_stacks(wuss);

        assert_eq!(stacks.len(), 1);
        assert_eq!(stacks[0].left_start, 0);
        assert_eq!(stacks[0].right_start, 9);
        assert_eq!(stacks[0].length, 3);
    }

    #[test]
    fn test_structure_stacks_multiple() {
        let wuss = "((..))((...))";
        let stacks = structure_stacks(wuss);

        assert_eq!(stacks.len(), 2);
    }

    #[test]
    fn test_structure_stacks_nested() {
        let wuss = "((((....))))";
        let stacks = structure_stacks(wuss);

        assert_eq!(stacks.len(), 1);
        assert_eq!(stacks[0].length, 4);
    }

    #[test]
    fn test_base_to_index() {
        assert_eq!(base_to_index(b'A'), Some(0));
        assert_eq!(base_to_index(b'a'), Some(0));
        assert_eq!(base_to_index(b'C'), Some(1));
        assert_eq!(base_to_index(b'G'), Some(2));
        assert_eq!(base_to_index(b'T'), Some(3));
        assert_eq!(base_to_index(b'U'), Some(3));
        assert_eq!(base_to_index(b'N'), None);
        assert_eq!(base_to_index(b'-'), None);
    }

    #[test]
    fn test_is_gap_char() {
        assert!(is_gap_char(b'-'));
        assert!(is_gap_char(b'.'));
        assert!(is_gap_char(b'_'));
        assert!(is_gap_char(b'~'));
        assert!(!is_gap_char(b'A'));
        assert!(!is_gap_char(b'N'));
    }

    #[test]
    fn test_guide_tree() {
        let mut guide = GuideTree::new();

        let root_idx = guide.add_node(GuideNode {
            node_type: ROOT_NODE,
            left: 1,
            right: -1,
            emit_pos: (0, 10),
        });

        let child_idx = guide.add_node(GuideNode {
            node_type: MATP_NODE,
            left: -1,
            right: -1,
            emit_pos: (1, 9),
        });

        assert_eq!(root_idx, 0);
        assert_eq!(child_idx, 1);
        assert_eq!(guide.nodes.len(), 2);
    }

    #[test]
    fn test_cm_from_gparse() {
        let mut gparse = GrammarParse::new();

        let root = gparse.add_node(GrammarParseNode {
            node_type: ROOT_NODE,
            left_pos: 0,
            right_pos: 10,
            children: vec![1],
        });

        gparse.add_node(GrammarParseNode {
            node_type: END_NODE,
            left_pos: 0,
            right_pos: 0,
            children: vec![],
        });

        gparse.root = root;

        let cm = cm_from_gparse(&gparse);

        assert_eq!(cm.nodes, 2);
        assert_eq!(cm.nd[0].node_type, ROOT_NODE as i32);
        assert_eq!(cm.nd[0].nxt, 1);
        assert_eq!(cm.nd[1].node_type, END_NODE as i32);
    }

    #[test]
    fn test_number_master_trace() {
        let mut trace = MasterTrace::new(0, 5, ROOT_NODE);
        trace.nxtl = Some(Box::new(MasterTrace::new(1, 4, MATP_NODE)));

        let count = number_master_trace(&mut trace);

        assert_eq!(count, 2);
        assert_eq!(trace.nodeidx, 0);
        assert_eq!(trace.nxtl.as_ref().unwrap().nodeidx, 1);
    }

    #[test]
    fn test_trace_pool() {
        let mut pool = TracePool::new();

        let t1 = pool.alloc();
        assert_eq!(t1.nodeidx, 0);

        pool.free(t1);

        let t2 = pool.alloc();
        assert_eq!(t2.nodeidx, 0); // Should reuse the freed trace
    }
}
