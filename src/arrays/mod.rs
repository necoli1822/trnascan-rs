//! Array containers for tRNAscan-SE
//!
//! This module provides collection types for managing multiple tRNA objects
//! and CM scan results with array operations.
//!
//! Ported from:
//! - tRNAscanSE::ArraytRNA.pm (672 lines)
//! - tRNAscanSE::ArrayCMscanResults.pm (352 lines)

use crate::trna::TRna;
use crate::cm_scan::CMSearchHit;

// ============================================================================
// ArrayTRna - Collection of tRNA objects
// ============================================================================

/// Collection of tRNA objects with array operations
///
/// Provides storage and manipulation of multiple tRNA predictions,
/// including sorting, filtering, searching, and deduplication.
///
/// Ported from tRNAscanSE::ArraytRNA.pm
#[derive(Debug, Clone, Default)]
pub struct ArrayTRna {
    /// Internal storage of tRNA objects
    items: Vec<TRna>,
}

impl ArrayTRna {
    /// Create a new empty array
    pub fn new() -> Self {
        ArrayTRna { items: Vec::new() }
    }

    /// Create with preallocated capacity
    pub fn with_capacity(capacity: usize) -> Self {
        ArrayTRna {
            items: Vec::with_capacity(capacity),
        }
    }

    /// Add a tRNA to the collection
    pub fn push(&mut self, trna: TRna) {
        self.items.push(trna);
    }

    /// Get a reference to a tRNA at index
    pub fn get(&self, index: usize) -> Option<&TRna> {
        self.items.get(index)
    }

    /// Get a mutable reference to a tRNA at index
    pub fn get_mut(&mut self, index: usize) -> Option<&mut TRna> {
        self.items.get_mut(index)
    }

    /// Get number of tRNAs in collection
    pub fn len(&self) -> usize {
        self.items.len()
    }

    /// Check if collection is empty
    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }

    /// Iterate over tRNAs
    pub fn iter(&self) -> impl Iterator<Item = &TRna> {
        self.items.iter()
    }

    /// Iterate over tRNAs mutably
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut TRna> {
        self.items.iter_mut()
    }

    /// Clear all tRNAs
    pub fn clear(&mut self) {
        self.items.clear();
    }

    /// Remove tRNA at index and return it
    ///
    /// # Panics
    /// Panics if index is out of bounds
    pub fn remove(&mut self, index: usize) -> TRna {
        self.items.remove(index)
    }

    /// Insert tRNA at specific index
    pub fn insert(&mut self, index: usize, trna: TRna) {
        self.items.insert(index, trna);
    }

    // ========================================================================
    // Sorting Operations
    // ========================================================================

    /// Sort by genomic position (seqname, then start)
    ///
    /// Corresponds to Perl's sort_by_seqname_start and sort_by_coord
    pub fn sort_by_position(&mut self) {
        self.items.sort_by(|a, b| {
            a.ordered_seqname
                .cmp(&b.ordered_seqname)
                .then_with(|| a.start.cmp(&b.start))
        });
    }

    /// Sort by score (highest first)
    ///
    /// Corresponds to Perl's sort_by_score_coord
    pub fn sort_by_score(&mut self) {
        self.items.sort_by(|a, b| {
            a.is_mito()
                .cmp(&b.is_mito())
                .then_with(|| {
                    b.score
                        .partial_cmp(&a.score)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .then_with(|| a.ordered_seqname.cmp(&b.ordered_seqname))
                .then_with(|| a.start.cmp(&b.start))
        });
    }

    /// Sort by isotype alphabetically, then anticodon, then score
    ///
    /// Corresponds to Perl's sort_by_isotype and sort_by_ac_matscore
    pub fn sort_by_isotype(&mut self) {
        self.items.sort_by(|a, b| {
            a.is_mito()
                .cmp(&b.is_mito())
                .then_with(|| a.is_numt().cmp(&b.is_numt()))
                .then_with(|| a.isotype.cmp(&b.isotype))
                .then_with(|| a.anticodon.cmp(&b.anticodon))
                .then_with(|| {
                    b.mat_score
                        .partial_cmp(&a.mat_score)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .then_with(|| a.seqname.cmp(&b.seqname))
                .then_with(|| a.start.cmp(&b.start))
        });
    }

    /// Sort by strand, then position
    ///
    /// Corresponds to Perl's sort_by_strand_coord
    pub fn sort_by_strand(&mut self) {
        self.items.sort_by(|a, b| {
            a.strand
                .to_char()
                .cmp(&b.strand.to_char())
                .then_with(|| a.ordered_seqname.cmp(&b.ordered_seqname))
                .then_with(|| a.start.cmp(&b.start))
        });
    }

    /// Sort for tRNAscan-SE output format
    ///
    /// Special sorting: forward strand by start, reverse strand by end (descending)
    /// Corresponds to Perl's sort_by_tRNAscanSE_output
    pub fn sort_for_output(&mut self) {
        self.items.sort_by(|a, b| {
            use crate::trna::Strand;

            // First by sequence
            let seq_cmp = a.ordered_seqname.cmp(&b.ordered_seqname);
            if seq_cmp != std::cmp::Ordering::Equal {
                return seq_cmp;
            }

            // Then by strand
            match (&a.strand, &b.strand) {
                (Strand::Plus, Strand::Plus) => a.start.cmp(&b.start),
                (Strand::Minus, Strand::Minus) => b.end.cmp(&a.end),
                (Strand::Plus, Strand::Minus) => std::cmp::Ordering::Less,
                (Strand::Minus, Strand::Plus) => std::cmp::Ordering::Greater,
                _ => std::cmp::Ordering::Equal,
            }
        });
    }

    /// Sort by tRNAscan ID
    ///
    /// Corresponds to Perl's sort_by_tRNAscan_id
    pub fn sort_by_id(&mut self) {
        self.items.sort_by(|a, b| a.trnascan_id.cmp(&b.trnascan_id));
    }

    // ========================================================================
    // Filtering Operations
    // ========================================================================

    /// Filter by minimum score, returning a new array
    pub fn filter_by_score(&self, min_score: f64) -> Self {
        ArrayTRna {
            items: self
                .items
                .iter()
                .filter(|t| t.score >= min_score)
                .cloned()
                .collect(),
        }
    }

    /// Filter by isotype, returning a new array
    pub fn filter_by_isotype(&self, isotype: &str) -> Self {
        ArrayTRna {
            items: self
                .items
                .iter()
                .filter(|t| t.isotype == isotype)
                .cloned()
                .collect(),
        }
    }

    /// Filter by anticodon, returning a new array
    pub fn filter_by_anticodon(&self, anticodon: &str) -> Self {
        ArrayTRna {
            items: self
                .items
                .iter()
                .filter(|t| t.anticodon == anticodon)
                .cloned()
                .collect(),
        }
    }

    /// Filter by strand
    pub fn filter_by_strand(&self, strand: char) -> Self {
        ArrayTRna {
            items: self
                .items
                .iter()
                .filter(|t| t.strand.to_char() == strand)
                .cloned()
                .collect(),
        }
    }

    /// Filter by sequence name
    pub fn filter_by_seqname(&self, seqname: &str) -> Self {
        ArrayTRna {
            items: self
                .items
                .iter()
                .filter(|t| t.seqname == seqname)
                .cloned()
                .collect(),
        }
    }

    /// Filter pseudogenes
    pub fn filter_pseudogenes(&self) -> Self {
        ArrayTRna {
            items: self.items.iter().filter(|t| t.is_pseudo).cloned().collect(),
        }
    }

    /// Filter non-pseudogenes
    pub fn filter_non_pseudogenes(&self) -> Self {
        ArrayTRna {
            items: self
                .items
                .iter()
                .filter(|t| !t.is_pseudo)
                .cloned()
                .collect(),
        }
    }

    // ========================================================================
    // Merging and Deduplication
    // ========================================================================

    /// Merge another ArrayTRna into this one
    pub fn merge(&mut self, other: Self) {
        self.items.extend(other.items);
    }

    /// Deduplicate tRNAs based on identical positions
    ///
    /// Keeps the first occurrence, removes subsequent duplicates.
    /// Array should be sorted first for best results.
    pub fn deduplicate(&mut self) {
        let mut seen = std::collections::HashSet::new();
        self.items.retain(|t| {
            let key = (t.seqname.clone(), t.start, t.end, t.strand.to_char());
            seen.insert(key)
        });
    }

    /// Deduplicate by keeping highest-scoring hit at each position
    pub fn deduplicate_by_score(&mut self) {
        use std::collections::HashMap;

        // Group by position
        let mut position_map: HashMap<(String, i64, i64, char), TRna> = HashMap::new();

        for trna in self.items.drain(..) {
            let key = (
                trna.seqname.clone(),
                trna.start,
                trna.end,
                trna.strand.to_char(),
            );

            position_map
                .entry(key)
                .and_modify(|existing| {
                    if trna.score > existing.score {
                        *existing = trna.clone();
                    }
                })
                .or_insert(trna);
        }

        // Rebuild array from unique entries
        self.items = position_map.into_values().collect();
    }

    // ========================================================================
    // Searching Operations
    // ========================================================================

    /// Binary search for tRNA by ID (assumes sorted by tRNAscan_id)
    ///
    /// Corresponds to Perl's bsearch_id
    pub fn binary_search_by_id(&self, id: &str) -> Result<usize, usize> {
        self.items
            .binary_search_by(|t| t.trnascan_id.as_str().cmp(id))
    }

    /// Linear search for tRNA by ID
    pub fn find_by_id(&self, id: &str) -> Option<usize> {
        self.items.iter().position(|t| t.trnascan_id == id)
    }

    /// Linear search for tRNA by position
    pub fn find_by_position(&self, start: i64, end: i64) -> Option<usize> {
        self.items
            .iter()
            .position(|t| t.start == start && t.end == end)
    }

    // ========================================================================
    // ID Management
    // ========================================================================

    /// Reorder all tRNA IDs sequentially by position
    ///
    /// Corresponds to Perl's reorder_all_tRNA_id
    /// Sorts by output format, then renumbers IDs sequentially per sequence
    pub fn reorder_ids(&mut self) {
        self.sort_for_output();

        let mut ct = 0;
        let mut prev_seqname = String::new();

        for trna in &mut self.items {
            if trna.seqname != prev_seqname {
                ct = 0;
                prev_seqname = trna.seqname.clone();
            }
            ct += 1;
            trna.id = ct;
            trna.trnascan_id = format!(
                "{}.tRNA{}-{}{}",
                trna.seqname, ct, trna.isotype, trna.anticodon
            );
        }
    }

    /// Reorder IDs for a specific sequence
    ///
    /// Corresponds to Perl's reorder_tRNA_id
    pub fn reorder_ids_for_seq(&mut self, seqname: &str) {
        self.sort_for_output();

        let mut ct = 0;
        for trna in &mut self.items {
            if trna.seqname == seqname {
                ct += 1;
                trna.id = ct;
                trna.trnascan_id = format!(
                    "{}.tRNA{}-{}{}",
                    trna.seqname, ct, trna.isotype, trna.anticodon
                );
            }
        }
    }

    // ========================================================================
    // Utility Methods
    // ========================================================================

    /// Get count by isotype
    pub fn count_by_isotype(&self, isotype: &str) -> usize {
        self.items.iter().filter(|t| t.isotype == isotype).count()
    }

    /// Get count by anticodon
    pub fn count_by_anticodon(&self, anticodon: &str) -> usize {
        self.items
            .iter()
            .filter(|t| t.anticodon == anticodon)
            .count()
    }

    /// Get total count of pseudogenes
    pub fn count_pseudogenes(&self) -> usize {
        self.items.iter().filter(|t| t.is_pseudo).count()
    }

    /// Convert to Vec for consumption
    pub fn into_vec(self) -> Vec<TRna> {
        self.items
    }

    /// Create from Vec
    pub fn from_vec(items: Vec<TRna>) -> Self {
        ArrayTRna { items }
    }
}

// ============================================================================
// ArrayCMScanResults - Collection of CMscan results
// ============================================================================

/// Collection of CM scan results
///
/// Provides storage and manipulation of multiple CMsearch hits,
/// including filtering, sorting, and merging operations.
///
/// Ported from tRNAscanSE::ArrayCMscanResults.pm
#[derive(Debug, Clone, Default)]
pub struct ArrayCMScanResults {
    /// Internal storage of CM search hits
    items: Vec<CMSearchHit>,
}

impl ArrayCMScanResults {
    /// Create a new empty array
    pub fn new() -> Self {
        ArrayCMScanResults { items: Vec::new() }
    }

    /// Create with preallocated capacity
    pub fn with_capacity(capacity: usize) -> Self {
        ArrayCMScanResults {
            items: Vec::with_capacity(capacity),
        }
    }

    /// Add a CM search hit to the collection
    pub fn push(&mut self, hit: CMSearchHit) {
        self.items.push(hit);
    }

    /// Get a reference to a hit at index
    pub fn get(&self, index: usize) -> Option<&CMSearchHit> {
        self.items.get(index)
    }

    /// Get a mutable reference to a hit at index
    pub fn get_mut(&mut self, index: usize) -> Option<&mut CMSearchHit> {
        self.items.get_mut(index)
    }

    /// Get number of hits in collection
    pub fn len(&self) -> usize {
        self.items.len()
    }

    /// Check if collection is empty
    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }

    /// Iterate over hits
    pub fn iter(&self) -> impl Iterator<Item = &CMSearchHit> {
        self.items.iter()
    }

    /// Iterate over hits mutably
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut CMSearchHit> {
        self.items.iter_mut()
    }

    /// Clear all hits
    pub fn clear(&mut self) {
        self.items.clear();
    }

    /// Remove hit at index and return it
    pub fn remove(&mut self, index: usize) -> CMSearchHit {
        self.items.remove(index)
    }

    // ========================================================================
    // Filtering Operations
    // ========================================================================

    /// Filter by maximum E-value
    pub fn filter_by_evalue(&self, max_evalue: f64) -> Self {
        ArrayCMScanResults {
            items: self
                .items
                .iter()
                .filter(|h| h.evalue <= max_evalue)
                .cloned()
                .collect(),
        }
    }

    /// Filter by minimum score
    pub fn filter_by_score(&self, min_score: f64) -> Self {
        ArrayCMScanResults {
            items: self
                .items
                .iter()
                .filter(|h| h.score >= min_score)
                .cloned()
                .collect(),
        }
    }

    /// Filter by target name (sequence name)
    pub fn filter_by_target(&self, target_name: &str) -> Self {
        ArrayCMScanResults {
            items: self
                .items
                .iter()
                .filter(|h| h.target_name == target_name)
                .cloned()
                .collect(),
        }
    }

    /// Filter by query model name
    pub fn filter_by_model(&self, model_name: &str) -> Self {
        ArrayCMScanResults {
            items: self
                .items
                .iter()
                .filter(|h| h.query_name == model_name)
                .cloned()
                .collect(),
        }
    }

    /// Filter by strand
    pub fn filter_by_strand(&self, strand: char) -> Self {
        ArrayCMScanResults {
            items: self
                .items
                .iter()
                .filter(|h| h.strand == strand)
                .cloned()
                .collect(),
        }
    }

    // ========================================================================
    // Sorting Operations
    // ========================================================================

    /// Sort by genomic position (target, then seq_from)
    ///
    /// Corresponds to Perl's sort_by_tRNAscanSE_output
    pub fn sort_by_position(&mut self) {
        self.items.sort_by(|a, b| {
            a.target_name
                .cmp(&b.target_name)
                .then_with(|| {
                    if a.strand == '+' && b.strand == '+' {
                        a.seq_from.cmp(&b.seq_from)
                    } else if a.strand == '-' && b.strand == '-' {
                        b.seq_to.cmp(&a.seq_to)
                    } else {
                        a.strand.cmp(&b.strand).then_with(|| a.seq_from.cmp(&b.seq_from))
                    }
                })
                .then_with(|| {
                    b.score
                        .partial_cmp(&a.score)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
        });
    }

    /// Sort by score (highest first)
    pub fn sort_by_score(&mut self) {
        self.items.sort_by(|a, b| {
            b.score
                .partial_cmp(&a.score)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }

    /// Sort by E-value (lowest first)
    pub fn sort_by_evalue(&mut self) {
        self.items.sort_by(|a, b| {
            a.evalue
                .partial_cmp(&b.evalue)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }

    // ========================================================================
    // Selection Operations
    // ========================================================================

    /// Get the highest-scoring hit
    pub fn best_hit(&self) -> Option<&CMSearchHit> {
        self.items
            .iter()
            .max_by(|a, b| {
                a.score
                    .partial_cmp(&b.score)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    }

    /// Get the hit with lowest E-value
    pub fn best_evalue_hit(&self) -> Option<&CMSearchHit> {
        self.items
            .iter()
            .min_by(|a, b| {
                a.evalue
                    .partial_cmp(&b.evalue)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    }

    // ========================================================================
    // Merging and Deduplication
    // ========================================================================

    /// Merge another ArrayCMScanResults into this one
    pub fn merge(&mut self, other: Self) {
        self.items.extend(other.items);
    }

    /// Remove overlapping hits, keeping higher-scoring ones
    ///
    /// Corresponds to Perl's merge_indexes logic
    /// Assumes hits are sorted by position
    pub fn remove_overlaps(&mut self, overlap_range: usize) {
        if self.items.len() < 2 {
            return;
        }

        let mut i = 0;
        while i < self.items.len() - 1 {
            let mut j = i + 1;
            while j < self.items.len() {
                let overlap = self.check_overlap(i, j, overlap_range);

                if overlap {
                    // Remove the lower-scoring hit
                    if self.items[i].score >= self.items[j].score {
                        self.items.remove(j);
                        // Don't increment j, check same position again
                    } else {
                        self.items.remove(i);
                        // Don't increment i, it's now pointing to next element
                        break;
                    }
                } else {
                    j += 1;
                }
            }
            i += 1;
        }
    }

    /// Check if two hits overlap
    fn check_overlap(&self, i: usize, j: usize, overlap_range: usize) -> bool {
        let a = &self.items[i];
        let b = &self.items[j];

        // Must be same target and strand
        if a.target_name != b.target_name || a.strand != b.strand {
            return false;
        }

        // Calculate overlap
        let a_start = a.seq_from.min(a.seq_to);
        let a_end = a.seq_from.max(a.seq_to);
        let b_start = b.seq_from.min(b.seq_to);
        let b_end = b.seq_from.max(b.seq_to);

        // Check for intersection
        if a_start > b_end || b_start > a_end {
            return false;
        }

        // Calculate overlap size
        let overlap_start = a_start.max(b_start);
        let overlap_end = a_end.min(b_end);
        let overlap_size = (overlap_end - overlap_start + 1) as usize;

        overlap_size >= overlap_range
    }

    // ========================================================================
    // Utility Methods
    // ========================================================================

    /// Convert to Vec for consumption
    pub fn into_vec(self) -> Vec<CMSearchHit> {
        self.items
    }

    /// Create from Vec
    pub fn from_vec(items: Vec<CMSearchHit>) -> Self {
        ArrayCMScanResults { items }
    }

    /// Get count by model name
    pub fn count_by_model(&self, model_name: &str) -> usize {
        self.items
            .iter()
            .filter(|h| h.query_name == model_name)
            .count()
    }

    /// Get count by target
    pub fn count_by_target(&self, target_name: &str) -> usize {
        self.items
            .iter()
            .filter(|h| h.target_name == target_name)
            .count()
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_array_trna_new() {
        let arr = ArrayTRna::new();
        assert!(arr.is_empty());
        assert_eq!(arr.len(), 0);
    }

    #[test]
    fn test_array_trna_push_get() {
        let mut arr = ArrayTRna::new();
        let mut trna = TRna::new();
        trna.seqname = "chr1".to_string();
        trna.start = 100;
        trna.end = 175;
        arr.push(trna.clone());

        assert_eq!(arr.len(), 1);
        assert_eq!(arr.get(0).unwrap().seqname, "chr1");
    }

    #[test]
    fn test_array_trna_filter_by_isotype() {
        let mut arr = ArrayTRna::new();

        let mut t1 = TRna::new();
        t1.isotype = "Ala".to_string();
        arr.push(t1);

        let mut t2 = TRna::new();
        t2.isotype = "Gly".to_string();
        arr.push(t2);

        let mut t3 = TRna::new();
        t3.isotype = "Ala".to_string();
        arr.push(t3);

        let filtered = arr.filter_by_isotype("Ala");
        assert_eq!(filtered.len(), 2);
    }

    #[test]
    fn test_array_trna_sort_by_position() {
        let mut arr = ArrayTRna::new();

        let mut t1 = TRna::new();
        t1.start = 300;
        t1.ordered_seqname = 1;
        arr.push(t1);

        let mut t2 = TRna::new();
        t2.start = 100;
        t2.ordered_seqname = 1;
        arr.push(t2);

        let mut t3 = TRna::new();
        t3.start = 200;
        t3.ordered_seqname = 1;
        arr.push(t3);

        arr.sort_by_position();

        assert_eq!(arr.get(0).unwrap().start, 100);
        assert_eq!(arr.get(1).unwrap().start, 200);
        assert_eq!(arr.get(2).unwrap().start, 300);
    }

    #[test]
    fn test_array_trna_deduplicate() {
        let mut arr = ArrayTRna::new();

        let mut t1 = TRna::new();
        t1.seqname = "chr1".to_string();
        t1.start = 100;
        t1.end = 175;
        arr.push(t1.clone());
        arr.push(t1.clone()); // Duplicate

        let mut t2 = TRna::new();
        t2.seqname = "chr1".to_string();
        t2.start = 200;
        t2.end = 275;
        arr.push(t2);

        arr.deduplicate();
        assert_eq!(arr.len(), 2);
    }

    #[test]
    fn test_array_cmscan_new() {
        let arr = ArrayCMScanResults::new();
        assert!(arr.is_empty());
        assert_eq!(arr.len(), 0);
    }

    #[test]
    fn test_array_cmscan_filter_by_score() {
        let mut arr = ArrayCMScanResults::new();

        let mut h1 = CMSearchHit::default();
        h1.score = 80.0;
        arr.push(h1);

        let mut h2 = CMSearchHit::default();
        h2.score = 50.0;
        arr.push(h2);

        let mut h3 = CMSearchHit::default();
        h3.score = 60.0;
        arr.push(h3);

        let filtered = arr.filter_by_score(55.0);
        assert_eq!(filtered.len(), 2);
    }

    #[test]
    fn test_array_cmscan_best_hit() {
        let mut arr = ArrayCMScanResults::new();

        let mut h1 = CMSearchHit::default();
        h1.score = 50.0;
        arr.push(h1);

        let mut h2 = CMSearchHit::default();
        h2.score = 80.0;
        arr.push(h2);

        let mut h3 = CMSearchHit::default();
        h3.score = 60.0;
        arr.push(h3);

        let best = arr.best_hit().unwrap();
        assert_eq!(best.score, 80.0);
    }

    #[test]
    fn test_array_cmscan_merge() {
        let mut arr1 = ArrayCMScanResults::new();
        arr1.push(CMSearchHit::default());

        let mut arr2 = ArrayCMScanResults::new();
        arr2.push(CMSearchHit::default());
        arr2.push(CMSearchHit::default());

        arr1.merge(arr2);
        assert_eq!(arr1.len(), 3);
    }
}
