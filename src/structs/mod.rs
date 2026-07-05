//! Data structure utilities
//!
//! This module implements utility functions for data structures from structs.c.
//! In Rust, most of these become trivial thanks to Vec and RAII, but we provide
//! the API for completeness and compatibility.

/// Allocate a 2D array
///
/// In Rust, this is trivial using Vec<Vec<T>>. This function is provided
/// for API compatibility with the C code.
///
/// # Example
/// ```
/// use trnascan_rs::structs::alloc_2d_array;
/// let arr: Vec<Vec<f64>> = alloc_2d_array(3, 4);
/// assert_eq!(arr.len(), 3);
/// assert_eq!(arr[0].len(), 4);
/// ```
pub fn alloc_2d_array<T: Default + Clone>(rows: usize, cols: usize) -> Vec<Vec<T>> {
    vec![vec![T::default(); cols]; rows]
}

/// Allocate a 3D array
///
/// # Example
/// ```
/// use trnascan_rs::structs::alloc_3d_array;
/// let arr: Vec<Vec<Vec<i32>>> = alloc_3d_array(2, 3, 4);
/// assert_eq!(arr.len(), 2);
/// assert_eq!(arr[0].len(), 3);
/// assert_eq!(arr[0][0].len(), 4);
/// ```
pub fn alloc_3d_array<T: Default + Clone>(d1: usize, d2: usize, d3: usize) -> Vec<Vec<Vec<T>>> {
    vec![vec![vec![T::default(); d3]; d2]; d1]
}

/// Copy a 2D array (deep copy)
///
/// # Example
/// ```
/// use trnascan_rs::structs::{alloc_2d_array, copy_2d_array};
/// let mut arr: Vec<Vec<i32>> = alloc_2d_array(2, 3);
/// arr[0][0] = 42;
/// let copy = copy_2d_array(&arr);
/// assert_eq!(copy[0][0], 42);
/// ```
pub fn copy_2d_array<T: Clone>(src: &[Vec<T>]) -> Vec<Vec<T>> {
    src.iter().map(|row| row.clone()).collect()
}

/// Copy a 3D array (deep copy)
pub fn copy_3d_array<T: Clone>(src: &[Vec<Vec<T>>]) -> Vec<Vec<Vec<T>>> {
    src.iter()
        .map(|plane| plane.iter().map(|row| row.clone()).collect())
        .collect()
}

/// Fill a 2D array with a value
pub fn fill_2d_array<T: Clone>(arr: &mut [Vec<T>], value: T) {
    for row in arr.iter_mut() {
        for cell in row.iter_mut() {
            *cell = value.clone();
        }
    }
}

/// Fill a 3D array with a value
pub fn fill_3d_array<T: Clone>(arr: &mut [Vec<Vec<T>>], value: T) {
    for plane in arr.iter_mut() {
        for row in plane.iter_mut() {
            for cell in row.iter_mut() {
                *cell = value.clone();
            }
        }
    }
}

/// Convert a state type to array index
///
/// Corresponds to StatetypeIndex() in structs.c lines 41-59.
/// Converts unique state type identifiers to valid array indices.
pub fn statetype_index(state_type: u32) -> usize {
    use crate::types::constants::*;

    match state_type {
        U_BEGIN_ST => BEGIN_ST,
        U_BIFURC_ST => BIFURC_ST,
        U_DEL_ST => DEL_ST,
        U_END_ST => END_ST,
        U_MATP_ST => MATP_ST,
        U_MATL_ST => MATL_ST,
        U_MATR_ST => MATR_ST,
        U_INSR_ST => INSR_ST,
        U_INSL_ST => INSL_ST,
        _ => panic!("Unknown state type: {}", state_type),
    }
}

/// Convert array index to unique state type
///
/// Corresponds to UniqueStatetype() in structs.c lines 67-88.
/// Converts an array index statetype into a unique statetype,
/// using additional information about the node type.
pub fn unique_statetype(node_type: i32, state_idx: usize) -> u32 {
    use crate::types::constants::*;

    match state_idx {
        DEL_ST => match node_type {
            -1 => U_END_ST,
            nt if nt == BIFURC_NODE as i32 => U_BIFURC_ST,
            nt if nt == BEGINL_NODE as i32 || nt == BEGINR_NODE as i32 => U_BEGIN_ST,
            _ => U_DEL_ST,
        },
        MATP_ST => U_MATP_ST,
        MATL_ST => U_MATL_ST,
        MATR_ST => U_MATR_ST,
        INSR_ST => U_INSR_ST,
        INSL_ST => U_INSL_ST,
        _ => panic!("Unknown state index: {}", state_idx),
    }
}

/// Alignment structure for storing sequence-model alignments
///
/// Corresponds to struct align_s in structs.c.
/// Used as a linked list for alignment of a model to a sequence.
#[derive(Debug, Clone)]
pub struct Align {
    /// Position in sequence
    pub pos: i32,
    /// ACGU base character
    pub sym: char,
    /// <.> secondary structure character
    pub ss: char,
    /// Node index in model
    pub nodeidx: i32,
    /// State type used
    pub state_type: i32,
}

impl Align {
    pub fn new(pos: i32, sym: char, ss: char, nodeidx: i32, state_type: i32) -> Self {
        Self {
            pos,
            sym,
            ss,
            nodeidx,
            state_type,
        }
    }
}

/// Alignment list builder
///
/// Rust-idiomatic replacement for the linked list alignment in C.
#[derive(Debug, Default)]
pub struct AlignList {
    items: Vec<Align>,
}

impl AlignList {
    pub fn new() -> Self {
        Self { items: Vec::new() }
    }

    /// Insert an alignment item after position `after_idx`
    pub fn insert_after(&mut self, after_idx: usize, item: Align) {
        if after_idx >= self.items.len() {
            self.items.push(item);
        } else {
            self.items.insert(after_idx + 1, item);
        }
    }

    /// Append an alignment item
    pub fn push(&mut self, item: Align) {
        self.items.push(item);
    }

    /// Remove item at index
    pub fn remove(&mut self, idx: usize) -> Option<Align> {
        if idx < self.items.len() {
            Some(self.items.remove(idx))
        } else {
            None
        }
    }

    /// Get all items
    pub fn items(&self) -> &[Align] {
        &self.items
    }

    /// Get mutable reference to all items
    pub fn items_mut(&mut self) -> &mut [Align] {
        &mut self.items
    }

    /// Number of items
    pub fn len(&self) -> usize {
        self.items.len()
    }

    /// Check if empty
    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }

    /// Clear all items
    pub fn clear(&mut self) {
        self.items.clear();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_alloc_2d_array() {
        let arr: Vec<Vec<i32>> = alloc_2d_array(3, 4);
        assert_eq!(arr.len(), 3);
        assert_eq!(arr[0].len(), 4);
        assert_eq!(arr[0][0], 0); // Default for i32
    }

    #[test]
    fn test_alloc_3d_array() {
        let arr: Vec<Vec<Vec<f64>>> = alloc_3d_array(2, 3, 4);
        assert_eq!(arr.len(), 2);
        assert_eq!(arr[0].len(), 3);
        assert_eq!(arr[0][0].len(), 4);
        assert_eq!(arr[0][0][0], 0.0);
    }

    #[test]
    fn test_copy_2d_array() {
        let mut arr: Vec<Vec<i32>> = alloc_2d_array(2, 3);
        arr[0][0] = 42;
        arr[1][2] = 99;

        let copy = copy_2d_array(&arr);
        assert_eq!(copy[0][0], 42);
        assert_eq!(copy[1][2], 99);

        // Verify deep copy
        assert!(std::ptr::eq(&arr[0], &arr[0])); // Same
        assert!(!std::ptr::eq(&arr[0], &copy[0])); // Different
    }

    #[test]
    fn test_fill_2d_array() {
        let mut arr: Vec<Vec<i32>> = alloc_2d_array(2, 3);
        fill_2d_array(&mut arr, 7);

        for row in &arr {
            for &cell in row {
                assert_eq!(cell, 7);
            }
        }
    }

    #[test]
    fn test_statetype_index() {
        use crate::types::constants::*;

        assert_eq!(statetype_index(U_DEL_ST), DEL_ST);
        assert_eq!(statetype_index(U_MATP_ST), MATP_ST);
        assert_eq!(statetype_index(U_MATL_ST), MATL_ST);
        assert_eq!(statetype_index(U_MATR_ST), MATR_ST);
        assert_eq!(statetype_index(U_INSL_ST), INSL_ST);
        assert_eq!(statetype_index(U_INSR_ST), INSR_ST);
    }

    #[test]
    fn test_unique_statetype() {
        use crate::types::constants::*;

        assert_eq!(unique_statetype(-1, DEL_ST), U_END_ST);
        assert_eq!(unique_statetype(BIFURC_NODE as i32, DEL_ST), U_BIFURC_ST);
        assert_eq!(unique_statetype(MATP_NODE as i32, DEL_ST), U_DEL_ST);
        assert_eq!(unique_statetype(0, MATP_ST), U_MATP_ST);
    }

    #[test]
    fn test_align_list() {
        let mut list = AlignList::new();
        assert!(list.is_empty());

        list.push(Align::new(0, 'A', '<', 0, 1));
        list.push(Align::new(1, 'C', '.', 1, 2));

        assert_eq!(list.len(), 2);
        assert_eq!(list.items()[0].sym, 'A');
        assert_eq!(list.items()[1].sym, 'C');

        // Insert between
        list.insert_after(0, Align::new(2, 'G', '>', 2, 3));
        assert_eq!(list.len(), 3);
        assert_eq!(list.items()[1].sym, 'G');

        // Remove
        let removed = list.remove(1);
        assert!(removed.is_some());
        assert_eq!(removed.unwrap().sym, 'G');
        assert_eq!(list.len(), 2);
    }
}
