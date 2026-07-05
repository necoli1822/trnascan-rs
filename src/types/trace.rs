//! Trace tree structures for alignment traceback
//!
//! This module implements the trace_s structure from structs.h lines 155-170.
//! Trace trees are binary trees used to store alignment tracebacks.

use crate::types::constants::*;

/// Trace tree structure from structs.h lines 160-170
///
/// Binary tree structure for storing a traceback of an alignment.
/// Also used for tracebacks of model constructions.
///
/// The tree represents the parse of a sequence through the CM:
/// - Leaf nodes represent emission of bases
/// - Internal nodes represent state transitions
/// - Bifurcation nodes have two children (nxtl and nxtr)
#[derive(Debug, Clone)]
pub struct Trace {
    /// Left position in sequence (1..N) or 0 if nothing emitted
    pub emitl: i32,

    /// Right position in sequence (1..N) or 0 if nothing emitted
    pub emitr: i32,

    /// Index of node responsible for this alignment
    pub nodeidx: i32,

    /// Type of substate used (U_MATP_ST, U_MATL_ST, etc.)
    /// These are unique state type flags from constants.rs
    pub trace_type: u32,

    /// Pointer to left (or only) branch, or None for leaf nodes
    pub nxtl: Option<Box<Trace>>,

    /// Pointer to right branch (BIFURC nodes only), else None
    pub nxtr: Option<Box<Trace>>,
}

impl Trace {
    /// Create a new trace node
    pub fn new(emitl: i32, emitr: i32, nodeidx: i32, trace_type: u32) -> Self {
        Self {
            emitl,
            emitr,
            nodeidx,
            trace_type,
            nxtl: None,
            nxtr: None,
        }
    }

    /// Create a leaf trace node (no children)
    pub fn leaf(emitl: i32, emitr: i32, nodeidx: i32, trace_type: u32) -> Self {
        Self::new(emitl, emitr, nodeidx, trace_type)
    }

    /// Create a trace node with one child (left)
    pub fn with_left(
        emitl: i32,
        emitr: i32,
        nodeidx: i32,
        trace_type: u32,
        left: Trace,
    ) -> Self {
        Self {
            emitl,
            emitr,
            nodeidx,
            trace_type,
            nxtl: Some(Box::new(left)),
            nxtr: None,
        }
    }

    /// Create a bifurcation trace node with two children
    pub fn bifurc(
        emitl: i32,
        emitr: i32,
        nodeidx: i32,
        left: Trace,
        right: Trace,
    ) -> Self {
        Self {
            emitl,
            emitr,
            nodeidx,
            trace_type: U_BIFURC_ST,
            nxtl: Some(Box::new(left)),
            nxtr: Some(Box::new(right)),
        }
    }

    /// Check if this is a leaf node
    pub fn is_leaf(&self) -> bool {
        self.nxtl.is_none() && self.nxtr.is_none()
    }

    /// Check if this is a bifurcation node
    pub fn is_bifurc(&self) -> bool {
        self.nxtl.is_some() && self.nxtr.is_some()
    }

    /// Check if this trace node is of a specific type
    pub fn is_type(&self, state_flag: u32) -> bool {
        self.trace_type & state_flag != 0
    }

    /// Get the span of this trace node
    pub fn span(&self) -> i32 {
        if self.emitr > 0 && self.emitl > 0 {
            self.emitr - self.emitl + 1
        } else {
            0
        }
    }

    /// Count total nodes in trace tree (self + descendants)
    pub fn count_nodes(&self) -> usize {
        let mut count = 1;

        if let Some(ref left) = self.nxtl {
            count += left.count_nodes();
        }

        if let Some(ref right) = self.nxtr {
            count += right.count_nodes();
        }

        count
    }

    /// Get depth of trace tree
    pub fn depth(&self) -> usize {
        let left_depth = self.nxtl.as_ref().map_or(0, |t| t.depth());
        let right_depth = self.nxtr.as_ref().map_or(0, |t| t.depth());
        1 + left_depth.max(right_depth)
    }

    /// Traverse trace tree in pre-order, calling visitor on each node
    pub fn traverse<F>(&self, visitor: &mut F)
    where
        F: FnMut(&Trace),
    {
        visitor(self);

        if let Some(ref left) = self.nxtl {
            left.traverse(visitor);
        }

        if let Some(ref right) = self.nxtr {
            right.traverse(visitor);
        }
    }
}

/// TraceStack for traversing trace trees
///
/// Implementation of the pushdown stack from trace.c lines 245-283.
/// Used for non-recursive traversal of trace trees.
#[derive(Debug)]
pub struct TraceStack {
    stack: Vec<Trace>,
}

impl TraceStack {
    /// Initialize a new trace stack
    /// Corresponds to InitTracestack() in trace.c lines 245-255
    pub fn new() -> Self {
        Self {
            stack: Vec::with_capacity(crate::types::constants::TSTACK_BLOCK),
        }
    }

    /// Push a trace onto the stack
    /// Corresponds to PushTracestack() in trace.c lines 257-267
    pub fn push(&mut self, tr: Trace) {
        self.stack.push(tr);
    }

    /// Pop a trace from the stack
    /// Corresponds to PopTracestack() in trace.c lines 268-277
    pub fn pop(&mut self) -> Option<Trace> {
        self.stack.pop()
    }

    /// Check if stack is empty
    pub fn is_empty(&self) -> bool {
        self.stack.is_empty()
    }

    /// Get the number of items on the stack
    pub fn len(&self) -> usize {
        self.stack.len()
    }
}

impl Default for TraceStack {
    fn default() -> Self {
        Self::new()
    }
}

/// Initialize a new trace tree with BEGIN and END nodes
/// Corresponds to InitTrace() in trace.c lines 44-78
pub fn init_trace() -> Trace {
    // Create END node (leaf)
    let end = Trace::new(-1, -1, -1, U_END_ST);

    // Create BEGIN node with END as child
    Trace {
        emitl: -1,
        emitr: -1,
        nodeidx: 0,
        trace_type: U_BEGIN_ST,
        nxtl: Some(Box::new(end)),
        nxtr: None,
    }
}

/// Attach a child trace to a parent node
/// Corresponds to AttachTrace() in trace.c lines 94-145
///
/// This creates a new trace node with the given parameters and attaches it
/// to the parent's left branch. For BIFURC nodes, the mechanics ensure
/// right branch is attached first, then left.
///
/// Returns a raw pointer to the newly attached node for further manipulation.
/// The pointer is valid as long as the parent trace exists.
pub fn attach_trace(
    parent: &mut Trace,
    emitl: i32,
    emitr: i32,
    nodeidx: i32,
    nodetype: u32,
) -> *mut Trace {
    // Create new trace node
    let new_trace = Trace::new(emitl, emitr, nodeidx, nodetype);

    // If parent already has a non-empty left branch (nxtl with children),
    // swap it to right and create new END on left
    if let Some(ref nxtl) = parent.nxtl {
        if nxtl.nxtl.is_some() {
            // Swap left to right
            parent.nxtr = parent.nxtl.take();

            // Create new END node for left
            let end = Trace::new(-1, -1, -1, U_END_ST);
            parent.nxtl = Some(Box::new(end));
        }
    }

    // Insert new trace between parent and its left child
    if let Some(old_left) = parent.nxtl.take() {
        let mut new_boxed = Box::new(new_trace);
        new_boxed.nxtl = Some(old_left);
        parent.nxtl = Some(new_boxed);
    } else {
        // No left child, just attach
        parent.nxtl = Some(Box::new(new_trace));
    }

    // Return pointer to the attached node
    parent
        .nxtl
        .as_mut()
        .map(|b| b.as_mut() as *mut Trace)
        .unwrap_or(std::ptr::null_mut())
}

/// Print trace tree for debugging
/// Corresponds to PrintTrace functionality
pub fn print_trace(tr: &Trace) {
    print_trace_recursive(tr, 0);
}

fn print_trace_recursive(tr: &Trace, indent: usize) {
    let indent_str = "  ".repeat(indent);
    let type_name = match tr.trace_type {
        U_DEL_ST => "DEL",
        U_MATP_ST => "MATP",
        U_MATL_ST => "MATL",
        U_MATR_ST => "MATR",
        U_INSL_ST => "INSL",
        U_INSR_ST => "INSR",
        U_BEGIN_ST => "BEGIN",
        U_END_ST => "END",
        U_BIFURC_ST => "BIFURC",
        _ => "UNKNOWN",
    };

    println!(
        "{}Node[{}]: type={}, emitl={}, emitr={}",
        indent_str, tr.nodeidx, type_name, tr.emitl, tr.emitr
    );

    if let Some(ref left) = tr.nxtl {
        println!("{}  Left:", indent_str);
        print_trace_recursive(left, indent + 2);
    }

    if let Some(ref right) = tr.nxtr {
        println!("{}  Right:", indent_str);
        print_trace_recursive(right, indent + 2);
    }
}

/// Collect all emitted positions from a trace tree
/// Returns (left_positions, right_positions)
pub fn collect_emissions(tr: &Trace) -> (Vec<i32>, Vec<i32>) {
    let mut left_pos = Vec::new();
    let mut right_pos = Vec::new();

    tr.traverse(&mut |node| {
        if node.emitl >= 0 {
            left_pos.push(node.emitl);
        }
        if node.emitr >= 0 {
            right_pos.push(node.emitr);
        }
    });

    (left_pos, right_pos)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_trace_stack_basic() {
        let mut stack = TraceStack::new();
        assert!(stack.is_empty());

        let trace1 = Trace::leaf(1, 2, 0, U_MATP_ST);
        let trace2 = Trace::leaf(3, 4, 1, U_MATL_ST);

        stack.push(trace1);
        stack.push(trace2);

        assert_eq!(stack.len(), 2);
        assert!(!stack.is_empty());

        let popped = stack.pop().unwrap();
        assert_eq!(popped.nodeidx, 1);

        let popped = stack.pop().unwrap();
        assert_eq!(popped.nodeidx, 0);

        assert!(stack.is_empty());
        assert!(stack.pop().is_none());
    }

    #[test]
    fn test_init_trace() {
        let trace = init_trace();
        assert_eq!(trace.trace_type, U_BEGIN_ST);
        assert_eq!(trace.nodeidx, 0);
        assert_eq!(trace.emitl, -1);
        assert!(trace.nxtl.is_some());

        let end = trace.nxtl.as_ref().unwrap();
        assert_eq!(end.trace_type, U_END_ST);
        assert_eq!(end.nodeidx, -1);
    }

    #[test]
    fn test_collect_emissions() {
        let left = Trace::leaf(1, 10, 1, U_MATP_ST);
        let right = Trace::leaf(5, 6, 2, U_MATP_ST);
        let root = Trace::bifurc(1, 10, 0, left, right);

        let (left_pos, right_pos) = collect_emissions(&root);

        // Root: 1, 10; Left: 1, 10; Right: 5, 6
        assert!(left_pos.contains(&1));
        assert!(left_pos.contains(&5));
        assert!(right_pos.contains(&10));
        assert!(right_pos.contains(&6));
    }

    #[test]
    fn test_trace_leaf() {
        let trace = Trace::leaf(5, 10, 3, U_MATP_ST);
        assert_eq!(trace.emitl, 5);
        assert_eq!(trace.emitr, 10);
        assert_eq!(trace.nodeidx, 3);
        assert_eq!(trace.trace_type, U_MATP_ST);
        assert!(trace.is_leaf());
        assert!(!trace.is_bifurc());
    }

    #[test]
    fn test_trace_with_left_child() {
        let child = Trace::leaf(1, 2, 1, U_MATL_ST);
        let parent = Trace::with_left(1, 5, 0, U_BEGIN_ST, child);

        assert_eq!(parent.emitl, 1);
        assert_eq!(parent.emitr, 5);
        assert!(parent.nxtl.is_some());
        assert!(parent.nxtr.is_none());
        assert!(!parent.is_leaf());
        assert!(!parent.is_bifurc());
    }

    #[test]
    fn test_trace_bifurc() {
        let left = Trace::leaf(1, 3, 1, U_MATL_ST);
        let right = Trace::leaf(4, 6, 2, U_MATR_ST);
        let bifurc = Trace::bifurc(1, 6, 0, left, right);

        assert_eq!(bifurc.trace_type, U_BIFURC_ST);
        assert!(bifurc.nxtl.is_some());
        assert!(bifurc.nxtr.is_some());
        assert!(!bifurc.is_leaf());
        assert!(bifurc.is_bifurc());
    }

    #[test]
    fn test_trace_type_checking() {
        let trace = Trace::leaf(1, 2, 0, U_MATP_ST);
        assert!(trace.is_type(U_MATP_ST));
        assert!(!trace.is_type(U_MATL_ST));

        // Test combined flags
        let multi_trace = Trace::leaf(1, 2, 0, U_MATP_ST | U_MATL_ST);
        assert!(multi_trace.is_type(U_MATP_ST));
        assert!(multi_trace.is_type(U_MATL_ST));
    }

    #[test]
    fn test_trace_span() {
        let trace1 = Trace::leaf(5, 10, 0, U_MATP_ST);
        assert_eq!(trace1.span(), 6); // 10 - 5 + 1

        let trace2 = Trace::leaf(0, 0, 0, U_DEL_ST);
        assert_eq!(trace2.span(), 0); // No emission

        let trace3 = Trace::leaf(7, 7, 0, U_MATL_ST);
        assert_eq!(trace3.span(), 1); // Single base
    }

    #[test]
    fn test_trace_count_nodes() {
        // Single leaf
        let leaf = Trace::leaf(1, 1, 0, U_MATL_ST);
        assert_eq!(leaf.count_nodes(), 1);

        // Parent with one child
        let child = Trace::leaf(2, 2, 1, U_MATL_ST);
        let parent = Trace::with_left(1, 3, 0, U_BEGIN_ST, child);
        assert_eq!(parent.count_nodes(), 2);

        // Bifurcation with two children
        let left = Trace::leaf(1, 1, 1, U_MATL_ST);
        let right = Trace::leaf(2, 2, 2, U_MATR_ST);
        let bifurc = Trace::bifurc(1, 2, 0, left, right);
        assert_eq!(bifurc.count_nodes(), 3);
    }

    #[test]
    fn test_trace_depth() {
        // Single leaf
        let leaf = Trace::leaf(1, 1, 0, U_MATL_ST);
        assert_eq!(leaf.depth(), 1);

        // Chain of depth 2
        let child = Trace::leaf(2, 2, 1, U_MATL_ST);
        let parent = Trace::with_left(1, 3, 0, U_BEGIN_ST, child);
        assert_eq!(parent.depth(), 2);

        // Bifurcation with unbalanced children
        let left = Trace::leaf(1, 1, 1, U_MATL_ST);
        let right_child = Trace::leaf(3, 3, 3, U_MATL_ST);
        let right = Trace::with_left(2, 4, 2, U_MATR_ST, right_child);
        let bifurc = Trace::bifurc(1, 4, 0, left, right);
        assert_eq!(bifurc.depth(), 3); // 1 (root) + 2 (right subtree)
    }

    #[test]
    fn test_trace_traverse() {
        // Build a small tree:
        //       root
        //      /    \
        //    left   right
        let left = Trace::leaf(1, 2, 1, U_MATL_ST);
        let right = Trace::leaf(3, 4, 2, U_MATR_ST);
        let root = Trace::bifurc(1, 4, 0, left, right);

        let mut visited = Vec::new();
        root.traverse(&mut |t| visited.push(t.nodeidx));

        // Pre-order traversal: root, left, right
        assert_eq!(visited, vec![0, 1, 2]);
    }

    #[test]
    fn test_trace_tree_structure() {
        // Test complex tree structure
        let leaf1 = Trace::leaf(1, 1, 3, U_MATL_ST);
        let _leaf2 = Trace::leaf(2, 2, 4, U_MATL_ST);
        let subtree = Trace::with_left(1, 2, 2, U_MATP_ST, leaf1);

        let left_branch = Trace::with_left(1, 3, 1, U_BEGIN_ST, subtree);
        let right_branch = Trace::leaf(4, 5, 5, U_MATR_ST);
        let root = Trace::bifurc(1, 5, 0, left_branch, right_branch);

        assert_eq!(root.count_nodes(), 5);
        assert_eq!(root.depth(), 4);
        assert!(root.is_bifurc());
    }
}
