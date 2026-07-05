//! Scan results module for tRNAscan-SE.
//!
//! This module provides the `ScanResult` struct for collecting and managing
//! tRNA scan results. It corresponds to the Perl `tRNAscanSE::ScanResult.pm` module.
//!
//! # Overview
//!
//! The `ScanResult` struct aggregates all tRNA hits found during a scan,
//! providing methods for:
//! - Adding and organizing tRNA hits
//! - Filtering by isotype, score, or strand
//! - Sorting by various criteria
//! - Computing statistics
//! - Merging results from multiple scans
//! - Detecting and resolving overlapping hits
//! - Writing output in various formats
//!
//! # Example
//!
//! ```rust,ignore
//! use trnascan_rs::scan_result::ScanResult;
//! use trnascan_rs::pipeline::scanner::TrnaHit;
//!
//! let mut result = ScanResult::new();
//! result.set_source("genome.fasta".to_string());
//! result.set_seq_name("chr1".to_string());
//!
//! // Add tRNA hits
//! let trna = TrnaHit::new("chr1", 1, 1000, 1075, "Leu", "CAA");
//! result.add_trna(trna);
//!
//! // Calculate statistics
//! result.calculate_stats();
//!
//! // Get distribution
//! let isotype_dist = result.get_isotype_distribution();
//! ```

use std::collections::HashMap;
use std::io::{self, Write};

use crate::pipeline::scanner::TrnaHit;

/// A container for tRNA scan results from a sequence or file.
///
/// This struct collects all tRNA hits found during scanning, tracks
/// pseudogenes separately, and provides methods for filtering, sorting,
/// statistics, and output generation.
#[derive(Debug, Clone, Default)]
pub struct ScanResult {
    // Sequence information
    /// Name of the sequence being scanned
    pub seq_name: String,
    /// Length of the sequence
    pub seq_len: usize,

    // Results collections
    /// All tRNA hits found (non-pseudogene)
    pub trnas: Vec<TrnaHit>,
    /// Pseudogene hits (kept separate for reporting)
    pub pseudogenes: Vec<TrnaHit>,

    // Statistics (computed by calculate_stats)
    /// Total tRNAs found (non-pseudogene)
    pub total_found: usize,
    /// Total pseudogenes found
    pub total_pseudo: usize,
    /// Count per isotype
    pub isotype_counts: HashMap<String, usize>,
    /// Count per anticodon
    pub anticodon_counts: HashMap<String, usize>,

    // Metadata
    /// Processing time in seconds
    pub processing_time: f64,
    /// Source file or sequence name
    pub source: String,

    // Search options used
    /// Search mode used (e.g., "Eukaryotic", "Bacterial")
    pub mode: String,
    /// Organellar mode (e.g., "general", "mito", "chloroplast")
    pub org_mode: String,

    // Internal tracking
    /// Next tRNA ID to assign
    next_id: i32,
    /// Whether stats have been calculated
    stats_valid: bool,
}

impl ScanResult {
    /// Create a new empty scan result container.
    pub fn new() -> Self {
        Self {
            seq_name: String::new(),
            seq_len: 0,
            trnas: Vec::new(),
            pseudogenes: Vec::new(),
            total_found: 0,
            total_pseudo: 0,
            isotype_counts: HashMap::new(),
            anticodon_counts: HashMap::new(),
            processing_time: 0.0,
            source: String::new(),
            mode: String::from("Eukaryotic"),
            org_mode: String::from("general"),
            next_id: 1,
            stats_valid: false,
        }
    }

    /// Create a scan result with sequence information.
    pub fn with_seq_info(seq_name: String, seq_len: usize) -> Self {
        Self {
            seq_name,
            seq_len,
            ..Self::new()
        }
    }

    // =========================================================================
    // Setters
    // =========================================================================

    /// Set the sequence name.
    pub fn set_seq_name(&mut self, name: String) {
        self.seq_name = name;
    }

    /// Set the sequence length.
    pub fn set_seq_len(&mut self, len: usize) {
        self.seq_len = len;
    }

    /// Set the source file or identifier.
    pub fn set_source(&mut self, source: String) {
        self.source = source;
    }

    /// Set the search mode.
    pub fn set_mode(&mut self, mode: String) {
        self.mode = mode;
    }

    /// Set the organellar mode.
    pub fn set_org_mode(&mut self, org_mode: String) {
        self.org_mode = org_mode;
    }

    /// Set processing time.
    pub fn set_processing_time(&mut self, time: f64) {
        self.processing_time = time;
    }

    // =========================================================================
    // Adding Results
    // =========================================================================

    /// Add a tRNA hit to the results.
    ///
    /// If the hit is a pseudogene, it will be added to the pseudogene list.
    /// Otherwise, it will be added to the main tRNA list.
    pub fn add_trna(&mut self, mut trna: TrnaHit) {
        // Assign an ID if not set
        if trna.trna_num == 0 {
            trna.trna_num = self.next_id;
            self.next_id += 1;
        } else if trna.trna_num >= self.next_id {
            self.next_id = trna.trna_num + 1;
        }

        // Invalidate stats
        self.stats_valid = false;

        // Add to appropriate list
        if trna.is_pseudo {
            self.pseudogenes.push(trna);
        } else {
            self.trnas.push(trna);
        }
    }

    /// Add a pseudogene hit explicitly.
    ///
    /// This ensures the hit is marked as a pseudogene and added
    /// to the pseudogene list.
    pub fn add_pseudogene(&mut self, mut trna: TrnaHit) {
        trna.is_pseudo = true;

        // Assign an ID if not set
        if trna.trna_num == 0 {
            trna.trna_num = self.next_id;
            self.next_id += 1;
        } else if trna.trna_num >= self.next_id {
            self.next_id = trna.trna_num + 1;
        }

        self.stats_valid = false;
        self.pseudogenes.push(trna);
    }

    /// Add multiple tRNAs at once.
    pub fn add_trnas(&mut self, trnas: Vec<TrnaHit>) {
        for trna in trnas {
            self.add_trna(trna);
        }
    }

    // =========================================================================
    // Filtering
    // =========================================================================

    /// Filter tRNAs by isotype.
    ///
    /// Returns references to all tRNAs matching the given isotype.
    pub fn filter_by_isotype(&self, isotype: &str) -> Vec<&TrnaHit> {
        self.trnas
            .iter()
            .filter(|t| t.isotype == isotype)
            .collect()
    }

    /// Filter tRNAs by minimum score.
    ///
    /// Returns references to all tRNAs with score >= min_score.
    pub fn filter_by_score(&self, min_score: f64) -> Vec<&TrnaHit> {
        self.trnas
            .iter()
            .filter(|t| t.inf_score >= min_score)
            .collect()
    }

    /// Filter tRNAs by strand.
    ///
    /// `strand` should be '+' for forward or '-' for reverse.
    pub fn filter_by_strand(&self, strand: char) -> Vec<&TrnaHit> {
        let forward = strand == '+';
        self.trnas
            .iter()
            .filter(|t| t.forward_strand == forward)
            .collect()
    }

    /// Filter tRNAs by anticodon.
    pub fn filter_by_anticodon(&self, anticodon: &str) -> Vec<&TrnaHit> {
        self.trnas
            .iter()
            .filter(|t| t.anticodon == anticodon)
            .collect()
    }

    /// Filter tRNAs that have introns.
    pub fn filter_with_introns(&self) -> Vec<&TrnaHit> {
        self.trnas.iter().filter(|t| t.has_intron()).collect()
    }

    /// Filter tRNAs within a genomic region.
    pub fn filter_by_region(&self, start: i32, end: i32) -> Vec<&TrnaHit> {
        self.trnas
            .iter()
            .filter(|t| t.start >= start && t.end <= end)
            .collect()
    }

    // =========================================================================
    // Sorting
    // =========================================================================

    /// Sort tRNAs by genomic position.
    ///
    /// Sorts by start position, then by end position for ties.
    pub fn sort_by_position(&mut self) {
        self.trnas.sort_by(|a, b| {
            a.start.cmp(&b.start).then_with(|| a.end.cmp(&b.end))
        });
        self.pseudogenes.sort_by(|a, b| {
            a.start.cmp(&b.start).then_with(|| a.end.cmp(&b.end))
        });
    }

    /// Sort tRNAs by score (highest first).
    pub fn sort_by_score(&mut self) {
        self.trnas.sort_by(|a, b| {
            b.inf_score
                .partial_cmp(&a.inf_score)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        self.pseudogenes.sort_by(|a, b| {
            b.inf_score
                .partial_cmp(&a.inf_score)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }

    /// Sort tRNAs by isotype alphabetically.
    pub fn sort_by_isotype(&mut self) {
        self.trnas.sort_by(|a, b| {
            a.isotype
                .cmp(&b.isotype)
                .then_with(|| a.anticodon.cmp(&b.anticodon))
        });
        self.pseudogenes.sort_by(|a, b| {
            a.isotype
                .cmp(&b.isotype)
                .then_with(|| a.anticodon.cmp(&b.anticodon))
        });
    }

    /// Sort tRNAs by tRNA number (ID).
    pub fn sort_by_id(&mut self) {
        self.trnas.sort_by_key(|t| t.trna_num);
        self.pseudogenes.sort_by_key(|t| t.trna_num);
    }

    /// Sort tRNAs for BED output (by seq_name, then position).
    pub fn sort_for_bed(&mut self) {
        self.trnas.sort_by(|a, b| {
            a.seq_name
                .cmp(&b.seq_name)
                .then_with(|| a.start.cmp(&b.start))
        });
    }

    // =========================================================================
    // Statistics
    // =========================================================================

    /// Calculate statistics from current results.
    ///
    /// Updates `total_found`, `total_pseudo`, `isotype_counts`,
    /// and `anticodon_counts`.
    pub fn calculate_stats(&mut self) {
        self.total_found = self.trnas.len();
        self.total_pseudo = self.pseudogenes.len();

        // Clear and recalculate counts
        self.isotype_counts.clear();
        self.anticodon_counts.clear();

        // Count from tRNAs
        for trna in &self.trnas {
            *self.isotype_counts.entry(trna.isotype.clone()).or_insert(0) += 1;
            *self
                .anticodon_counts
                .entry(trna.anticodon.clone())
                .or_insert(0) += 1;
        }

        self.stats_valid = true;
    }

    /// Get the isotype distribution.
    ///
    /// Returns a map of isotype name to count.
    pub fn get_isotype_distribution(&self) -> HashMap<String, usize> {
        let mut counts = HashMap::new();
        for trna in &self.trnas {
            *counts.entry(trna.isotype.clone()).or_insert(0) += 1;
        }
        counts
    }

    /// Get the anticodon distribution.
    ///
    /// Returns a map of anticodon to count.
    pub fn get_anticodon_distribution(&self) -> HashMap<String, usize> {
        let mut counts = HashMap::new();
        for trna in &self.trnas {
            *counts.entry(trna.anticodon.clone()).or_insert(0) += 1;
        }
        counts
    }

    /// Get count of tRNAs with introns.
    pub fn intron_count(&self) -> usize {
        self.trnas.iter().filter(|t| t.has_intron()).count()
    }

    /// Get total base pairs covered by tRNAs.
    pub fn total_bp_covered(&self) -> i32 {
        self.trnas.iter().map(|t| t.full_length()).sum()
    }

    /// Get average tRNA length.
    pub fn avg_trna_length(&self) -> f64 {
        if self.trnas.is_empty() {
            return 0.0;
        }
        let total: i32 = self.trnas.iter().map(|t| t.full_length()).sum();
        total as f64 / self.trnas.len() as f64
    }

    /// Get count of standard amino acid tRNAs (excluding special types).
    pub fn standard_aa_count(&self) -> usize {
        self.trnas
            .iter()
            .filter(|t| {
                !t.is_pseudo
                    && t.isotype != "SeC"
                    && t.isotype != "Sup"
                    && t.isotype != "Unk"
                    && t.isotype != "Undet"
            })
            .count()
    }

    // =========================================================================
    // Merging
    // =========================================================================

    /// Merge another ScanResult into this one.
    ///
    /// All tRNAs and pseudogenes from `other` are added to this result.
    /// Statistics are invalidated and must be recalculated.
    pub fn merge(&mut self, other: ScanResult) {
        for trna in other.trnas {
            self.add_trna(trna);
        }
        for pseudo in other.pseudogenes {
            self.add_pseudogene(pseudo);
        }
        self.stats_valid = false;
    }

    /// Remove overlapping hits, keeping the higher-scoring one.
    ///
    /// This uses a simple greedy algorithm:
    /// 1. Sort by score (highest first)
    /// 2. For each hit, remove any overlapping hits with lower scores
    pub fn deduplicate(&mut self) {
        self.resolve_overlaps();
    }

    // =========================================================================
    // Overlap Detection
    // =========================================================================

    /// Check if two tRNAs overlap.
    ///
    /// Two tRNAs overlap if their genomic ranges intersect.
    pub fn overlaps(a: &TrnaHit, b: &TrnaHit) -> bool {
        // Must be on same sequence and strand
        if a.seq_name != b.seq_name || a.forward_strand != b.forward_strand {
            return false;
        }

        let a_start = a.start.min(a.end);
        let a_end = a.start.max(a.end);
        let b_start = b.start.min(b.end);
        let b_end = b.start.max(b.end);

        // Check for intersection
        a_start <= b_end && b_start <= a_end
    }

    /// Get the amount of overlap between two tRNAs.
    ///
    /// Returns the number of overlapping base pairs, or 0 if no overlap.
    pub fn overlap_amount(a: &TrnaHit, b: &TrnaHit) -> i64 {
        if !Self::overlaps(a, b) {
            return 0;
        }

        let a_start = a.start.min(a.end) as i64;
        let a_end = a.start.max(a.end) as i64;
        let b_start = b.start.min(b.end) as i64;
        let b_end = b.start.max(b.end) as i64;

        let overlap_start = a_start.max(b_start);
        let overlap_end = a_end.min(b_end);

        overlap_end - overlap_start + 1
    }

    /// Resolve overlapping hits by keeping the higher-scoring one.
    ///
    /// This modifies the tRNA list in place, removing lower-scoring
    /// overlapping hits.
    pub fn resolve_overlaps(&mut self) {
        if self.trnas.len() < 2 {
            return;
        }

        // Sort by score (highest first)
        self.trnas.sort_by(|a, b| {
            b.inf_score
                .partial_cmp(&a.inf_score)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        // Keep track of which hits to retain
        let mut keep = vec![true; self.trnas.len()];

        for i in 0..self.trnas.len() {
            if !keep[i] {
                continue;
            }
            for j in (i + 1)..self.trnas.len() {
                if !keep[j] {
                    continue;
                }
                if Self::overlaps(&self.trnas[i], &self.trnas[j]) {
                    // Keep the higher-scoring one (i, since sorted)
                    keep[j] = false;
                }
            }
        }

        // Retain only the kept hits
        let mut idx = 0;
        self.trnas.retain(|_| {
            let k = keep[idx];
            idx += 1;
            k
        });

        self.stats_valid = false;
    }

    // =========================================================================
    // Output Methods
    // =========================================================================

    /// Write tabular output.
    ///
    /// Corresponds to the main .out file format from tRNAscan-SE.
    pub fn write_output<W: Write>(&self, w: &mut W) -> io::Result<()> {
        // Write header
        self.write_output_header(w)?;

        // Write each tRNA
        for trna in &self.trnas {
            self.write_trna_line(w, trna)?;
        }

        // Write pseudogenes
        for trna in &self.pseudogenes {
            self.write_trna_line(w, trna)?;
        }

        Ok(())
    }

    /// Write tabular output header.
    fn write_output_header<W: Write>(&self, w: &mut W) -> io::Result<()> {
        writeln!(
            w,
            "Sequence\t\ttRNA\tBounds\ttRNA\tAnti\tIntron Bounds\tInf\tIsotype\tIsotype\t"
        )?;
        writeln!(
            w,
            "Name     \ttRNA #\tBegin\tEnd  \tType\tCodon\tBegin\tEnd\tScore\tCM\tScore\tNote"
        )?;
        writeln!(
            w,
            "--------\t------\t-----\t-----\t----\t-----\t-----\t----\t-----\t-------\t-------\t------"
        )?;
        Ok(())
    }

    /// Write a single tRNA as a tabular line.
    fn write_trna_line<W: Write>(&self, w: &mut W, trna: &TrnaHit) -> io::Result<()> {
        let (disp_start, disp_end) = if trna.forward_strand {
            (trna.start, trna.end)
        } else {
            (trna.end, trna.start)
        };

        let intron_start = if trna.intron_start > 0 {
            trna.intron_start.to_string()
        } else {
            "0".to_string()
        };
        let intron_end = if trna.intron_end > 0 {
            trna.intron_end.to_string()
        } else {
            "0".to_string()
        };

        let note = if trna.is_pseudo {
            "pseudo".to_string()
        } else {
            trna.note.clone()
        };

        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.1}\t{}\t{:.1}\t{}",
            trna.seq_name,
            trna.trna_num,
            disp_start,
            disp_end,
            trna.isotype,
            trna.anticodon,
            intron_start,
            intron_end,
            trna.inf_score,
            trna.isotype_cm,
            trna.isotype_score,
            note
        )
    }

    /// Write secondary structure output.
    ///
    /// Corresponds to the .ss file format from tRNAscan-SE.
    pub fn write_ss_output<W: Write>(&self, w: &mut W) -> io::Result<()> {
        for trna in &self.trnas {
            self.write_ss_entry(w, trna)?;
        }
        for trna in &self.pseudogenes {
            self.write_ss_entry(w, trna)?;
        }
        Ok(())
    }

    /// Write a single tRNA's secondary structure entry.
    fn write_ss_entry<W: Write>(&self, w: &mut W, trna: &TrnaHit) -> io::Result<()> {
        let (disp_start, disp_end) = if trna.forward_strand {
            (trna.start, trna.end)
        } else {
            (trna.end, trna.start)
        };

        // Header line
        writeln!(
            w,
            "{}.trna{} ({}-{})\tLength: {} bp",
            trna.seq_name,
            trna.trna_num,
            disp_start,
            disp_end,
            trna.full_length()
        )?;

        // Type line
        writeln!(
            w,
            "Type: {}\tAnticodon: {} at {}-{}\tScore: {:.1}",
            trna.isotype,
            trna.anticodon,
            trna.anticodon_pos_start,
            trna.anticodon_pos_end,
            trna.inf_score
        )?;

        // Intron info
        if trna.has_intron() {
            let rel_start = trna.intron_start - trna.start + 1;
            let rel_end = trna.intron_end - trna.start + 1;
            writeln!(
                w,
                "Possible intron: {}-{} ({}-{})",
                rel_start, rel_end, trna.intron_start, trna.intron_end
            )?;
        }

        // Notes
        let mut note_line = String::new();
        if trna.is_pseudo {
            note_line.push_str("Possible pseudogene");
        }
        if trna.hmm_score > 0.0 || trna.ss_score > 0.0 {
            if !note_line.is_empty() {
                note_line.push_str(": ");
            }
            note_line.push_str(&format!(
                "HMM Sc={:.2}\tSec struct Sc={:.2}",
                trna.hmm_score, trna.ss_score
            ));
        }
        if !note_line.is_empty() {
            writeln!(w, "{}", note_line)?;
        }

        // Ruler
        let seq_len = trna.sequence.len();
        if seq_len > 0 {
            let mut ruler = String::with_capacity(seq_len + 10);
            ruler.push_str("     ");
            for i in 0..seq_len {
                if (i + 1) % 10 == 5 {
                    ruler.push('*');
                } else if (i + 1) % 10 == 0 {
                    ruler.push('|');
                } else {
                    ruler.push(' ');
                }
            }
            writeln!(w, "{}", ruler)?;
        }

        // Sequence
        writeln!(w, "Seq: {}", trna.sequence)?;

        // Structure
        writeln!(w, "Str: {}", trna.secondary_structure)?;

        // Blank separator
        writeln!(w)?;

        Ok(())
    }

    /// Write BED format output.
    ///
    /// BED12 format for genome browser visualization.
    pub fn write_bed_output<W: Write>(&self, w: &mut W) -> io::Result<()> {
        // Sort for BED output (by position)
        let mut sorted_trnas: Vec<_> = self.trnas.iter().collect();
        sorted_trnas.sort_by(|a, b| {
            a.seq_name
                .cmp(&b.seq_name)
                .then_with(|| a.start.cmp(&b.start))
        });

        for trna in sorted_trnas {
            self.write_bed_line(w, trna)?;
        }

        Ok(())
    }

    /// Write a single BED line.
    fn write_bed_line<W: Write>(&self, w: &mut W, trna: &TrnaHit) -> io::Result<()> {
        let (bed_start, bed_end) = trna.bed_coords();
        let strand = trna.strand_char();
        let name = trna.identifier();

        // Convert score to BED score (0-1000 range)
        let score = ((trna.inf_score * 10.0).min(1000.0).max(0.0)) as i32;

        if trna.has_intron() {
            // Two blocks (exon-intron-exon)
            let intron_rel_start = if trna.forward_strand {
                trna.intron_start - trna.start
            } else {
                trna.end - trna.intron_end
            };
            let intron_rel_end = if trna.forward_strand {
                trna.intron_end - trna.start + 1
            } else {
                trna.end - trna.intron_start + 1
            };

            let block1_size = intron_rel_start;
            let block2_size = trna.full_length() - intron_rel_end;

            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t0\t2\t{},{},\t0,{},",
                trna.seq_name,
                bed_start,
                bed_end,
                name,
                score,
                strand,
                bed_start,
                bed_end,
                block1_size,
                block2_size,
                intron_rel_end
            )
        } else {
            // Single block
            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t0\t1\t{},\t0,",
                trna.seq_name,
                bed_start,
                bed_end,
                name,
                score,
                strand,
                bed_start,
                bed_end,
                trna.full_length()
            )
        }
    }

    /// Write GFF3 format output.
    pub fn write_gff_output<W: Write>(&self, w: &mut W) -> io::Result<()> {
        // Header
        writeln!(w, "##gff-version 3")?;

        // Sort for output
        let mut sorted_trnas: Vec<_> = self.trnas.iter().chain(self.pseudogenes.iter()).collect();
        sorted_trnas.sort_by(|a, b| {
            a.seq_name
                .cmp(&b.seq_name)
                .then_with(|| a.start.cmp(&b.start))
        });

        for trna in sorted_trnas {
            self.write_gff_entry(w, trna)?;
        }

        Ok(())
    }

    /// Write a single GFF3 entry.
    fn write_gff_entry<W: Write>(&self, w: &mut W, trna: &TrnaHit) -> io::Result<()> {
        let biotype = if trna.is_pseudo { "pseudogene" } else { "tRNA" };
        let strand = trna.strand_char();
        let id = format!("{}.trna{}", trna.seq_name, trna.trna_num);

        // Main feature line
        writeln!(
            w,
            "{}\ttRNAscan-SE\t{}\t{}\t{}\t{:.1}\t{}\t.\tID={};Name={};isotype={};anticodon={};gene_biotype={};",
            trna.seq_name,
            biotype,
            trna.start,
            trna.end,
            trna.inf_score,
            strand,
            id,
            trna.identifier(),
            trna.isotype,
            trna.anticodon,
            biotype
        )?;

        // Exon features
        if trna.has_intron() {
            // Two exons
            if trna.forward_strand {
                writeln!(
                    w,
                    "{}\ttRNAscan-SE\texon\t{}\t{}\t.\t{}\t.\tID={}.exon1;Parent={};",
                    trna.seq_name,
                    trna.start,
                    trna.intron_start - 1,
                    strand,
                    id,
                    id
                )?;
                writeln!(
                    w,
                    "{}\ttRNAscan-SE\texon\t{}\t{}\t.\t{}\t.\tID={}.exon2;Parent={};",
                    trna.seq_name,
                    trna.intron_end + 1,
                    trna.end,
                    strand,
                    id,
                    id
                )?;
            } else {
                writeln!(
                    w,
                    "{}\ttRNAscan-SE\texon\t{}\t{}\t.\t{}\t.\tID={}.exon1;Parent={};",
                    trna.seq_name,
                    trna.intron_end + 1,
                    trna.end,
                    strand,
                    id,
                    id
                )?;
                writeln!(
                    w,
                    "{}\ttRNAscan-SE\texon\t{}\t{}\t.\t{}\t.\tID={}.exon2;Parent={};",
                    trna.seq_name,
                    trna.start,
                    trna.intron_start - 1,
                    strand,
                    id,
                    id
                )?;
            }
        } else {
            // Single exon
            writeln!(
                w,
                "{}\ttRNAscan-SE\texon\t{}\t{}\t.\t{}\t.\tID={}.exon1;Parent={};",
                trna.seq_name, trna.start, trna.end, strand, id, id
            )?;
        }

        Ok(())
    }

    /// Write statistics summary.
    pub fn write_stats<W: Write>(&self, w: &mut W) -> io::Result<()> {
        writeln!(w)?;
        writeln!(w, "tRNAscan-SE v.2.0 (Rust) scan results")?;
        writeln!(w, "------------------------------------------------------------")?;
        writeln!(w)?;

        if !self.source.is_empty() {
            writeln!(w, "Sequence file(s) to search:        {}", self.source)?;
        }
        writeln!(w, "Search Mode:                       {}", self.mode)?;
        if !self.org_mode.is_empty() && self.org_mode != "general" {
            writeln!(w, "Organellar Mode:                   {}", self.org_mode)?;
        }

        writeln!(w)?;
        writeln!(w, "------------------------------------------------------------")?;
        writeln!(w)?;

        writeln!(w, "Summary Statistics:")?;
        writeln!(w, "-------------------")?;
        writeln!(w)?;

        writeln!(w, "Total tRNAs:                        {}", self.trnas.len())?;
        writeln!(
            w,
            "tRNAs decoding Standard 20 AA:      {}",
            self.standard_aa_count()
        )?;
        writeln!(w, "Predicted pseudogenes:              {}", self.pseudogenes.len())?;
        writeln!(w, "tRNAs with introns:                 {}", self.intron_count())?;
        writeln!(w)?;

        // Isotype distribution
        writeln!(w, "Isotype / Anticodon Counts:")?;
        writeln!(w)?;

        let mut isotypes: Vec<_> = self.get_isotype_distribution().into_iter().collect();
        isotypes.sort_by(|a, b| a.0.cmp(&b.0));

        for (isotype, count) in isotypes {
            writeln!(w, "  {:6}: {}", isotype, count)?;
        }

        writeln!(w)?;

        if self.processing_time > 0.0 {
            writeln!(w, "Processing time: {:.2} seconds", self.processing_time)?;
        }

        Ok(())
    }

    // =========================================================================
    // Iteration
    // =========================================================================

    /// Iterate over all tRNAs (non-pseudogenes).
    pub fn iter(&self) -> impl Iterator<Item = &TrnaHit> {
        self.trnas.iter()
    }

    /// Iterate over all hits including pseudogenes.
    pub fn iter_all(&self) -> impl Iterator<Item = &TrnaHit> {
        self.trnas.iter().chain(self.pseudogenes.iter())
    }

    /// Iterate over pseudogenes only.
    pub fn iter_pseudogenes(&self) -> impl Iterator<Item = &TrnaHit> {
        self.pseudogenes.iter()
    }

    /// Get total count of all hits (tRNAs + pseudogenes).
    pub fn len(&self) -> usize {
        self.trnas.len() + self.pseudogenes.len()
    }

    /// Check if results are empty.
    pub fn is_empty(&self) -> bool {
        self.trnas.is_empty() && self.pseudogenes.is_empty()
    }

    /// Get count of tRNAs only (non-pseudogenes).
    pub fn trna_count(&self) -> usize {
        self.trnas.len()
    }

    /// Get count of pseudogenes.
    pub fn pseudogene_count(&self) -> usize {
        self.pseudogenes.len()
    }

    /// Clear all results.
    pub fn clear(&mut self) {
        self.trnas.clear();
        self.pseudogenes.clear();
        self.isotype_counts.clear();
        self.anticodon_counts.clear();
        self.total_found = 0;
        self.total_pseudo = 0;
        self.next_id = 1;
        self.stats_valid = false;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_test_hit(
        seq_name: &str,
        num: i32,
        start: i32,
        end: i32,
        isotype: &str,
        anticodon: &str,
    ) -> TrnaHit {
        TrnaHit::new(seq_name, num, start, end, isotype, anticodon)
    }

    #[test]
    fn test_new_scan_result() {
        let result = ScanResult::new();
        assert!(result.trnas.is_empty());
        assert!(result.pseudogenes.is_empty());
        assert_eq!(result.total_found, 0);
    }

    #[test]
    fn test_with_seq_info() {
        let result = ScanResult::with_seq_info("chr1".to_string(), 1000000);
        assert_eq!(result.seq_name, "chr1");
        assert_eq!(result.seq_len, 1000000);
    }

    #[test]
    fn test_add_trna() {
        let mut result = ScanResult::new();
        let hit = make_test_hit("seq1", 0, 100, 175, "Leu", "CAA");
        result.add_trna(hit);

        assert_eq!(result.trnas.len(), 1);
        assert_eq!(result.trnas[0].trna_num, 1); // Auto-assigned
    }

    #[test]
    fn test_add_pseudogene() {
        let mut result = ScanResult::new();
        let hit = make_test_hit("seq1", 0, 100, 175, "Leu", "CAA");
        result.add_pseudogene(hit);

        assert_eq!(result.pseudogenes.len(), 1);
        assert!(result.pseudogenes[0].is_pseudo);
        assert!(result.trnas.is_empty());
    }

    #[test]
    fn test_filter_by_isotype() {
        let mut result = ScanResult::new();
        result.add_trna(make_test_hit("seq1", 1, 100, 175, "Leu", "CAA"));
        result.add_trna(make_test_hit("seq1", 2, 200, 275, "Ser", "AGA"));
        result.add_trna(make_test_hit("seq1", 3, 300, 375, "Leu", "TAG"));

        let leu_hits = result.filter_by_isotype("Leu");
        assert_eq!(leu_hits.len(), 2);
    }

    #[test]
    fn test_filter_by_score() {
        let mut result = ScanResult::new();

        let mut hit1 = make_test_hit("seq1", 1, 100, 175, "Leu", "CAA");
        hit1.inf_score = 80.0;
        result.add_trna(hit1);

        let mut hit2 = make_test_hit("seq1", 2, 200, 275, "Ser", "AGA");
        hit2.inf_score = 50.0;
        result.add_trna(hit2);

        let high_score = result.filter_by_score(60.0);
        assert_eq!(high_score.len(), 1);
        assert_eq!(high_score[0].isotype, "Leu");
    }

    #[test]
    fn test_filter_by_strand() {
        let mut result = ScanResult::new();

        let mut hit1 = make_test_hit("seq1", 1, 100, 175, "Leu", "CAA");
        hit1.forward_strand = true;
        result.add_trna(hit1);

        let mut hit2 = make_test_hit("seq1", 2, 275, 200, "Ser", "AGA");
        hit2.forward_strand = false;
        result.add_trna(hit2);

        let forward = result.filter_by_strand('+');
        assert_eq!(forward.len(), 1);

        let reverse = result.filter_by_strand('-');
        assert_eq!(reverse.len(), 1);
    }

    #[test]
    fn test_sort_by_position() {
        let mut result = ScanResult::new();
        result.add_trna(make_test_hit("seq1", 3, 300, 375, "Phe", "GAA"));
        result.add_trna(make_test_hit("seq1", 1, 100, 175, "Leu", "CAA"));
        result.add_trna(make_test_hit("seq1", 2, 200, 275, "Ser", "AGA"));

        result.sort_by_position();

        assert_eq!(result.trnas[0].start, 100);
        assert_eq!(result.trnas[1].start, 200);
        assert_eq!(result.trnas[2].start, 300);
    }

    #[test]
    fn test_sort_by_score() {
        let mut result = ScanResult::new();

        let mut hit1 = make_test_hit("seq1", 1, 100, 175, "Leu", "CAA");
        hit1.inf_score = 50.0;
        result.add_trna(hit1);

        let mut hit2 = make_test_hit("seq1", 2, 200, 275, "Ser", "AGA");
        hit2.inf_score = 80.0;
        result.add_trna(hit2);

        result.sort_by_score();

        assert_eq!(result.trnas[0].inf_score, 80.0); // Highest first
        assert_eq!(result.trnas[1].inf_score, 50.0);
    }

    #[test]
    fn test_calculate_stats() {
        let mut result = ScanResult::new();
        result.add_trna(make_test_hit("seq1", 1, 100, 175, "Leu", "CAA"));
        result.add_trna(make_test_hit("seq1", 2, 200, 275, "Leu", "TAG"));
        result.add_trna(make_test_hit("seq1", 3, 300, 375, "Ser", "AGA"));

        result.calculate_stats();

        assert_eq!(result.total_found, 3);
        assert_eq!(result.isotype_counts.get("Leu"), Some(&2));
        assert_eq!(result.isotype_counts.get("Ser"), Some(&1));
    }

    #[test]
    fn test_get_distributions() {
        let mut result = ScanResult::new();
        result.add_trna(make_test_hit("seq1", 1, 100, 175, "Leu", "CAA"));
        result.add_trna(make_test_hit("seq1", 2, 200, 275, "Leu", "CAA"));
        result.add_trna(make_test_hit("seq1", 3, 300, 375, "Ser", "AGA"));

        let isotype_dist = result.get_isotype_distribution();
        assert_eq!(isotype_dist.get("Leu"), Some(&2));
        assert_eq!(isotype_dist.get("Ser"), Some(&1));

        let anticodon_dist = result.get_anticodon_distribution();
        assert_eq!(anticodon_dist.get("CAA"), Some(&2));
        assert_eq!(anticodon_dist.get("AGA"), Some(&1));
    }

    #[test]
    fn test_overlaps() {
        let mut hit1 = make_test_hit("seq1", 1, 100, 175, "Leu", "CAA");
        hit1.forward_strand = true;

        let mut hit2 = make_test_hit("seq1", 2, 150, 225, "Ser", "AGA");
        hit2.forward_strand = true;

        assert!(ScanResult::overlaps(&hit1, &hit2));

        let mut hit3 = make_test_hit("seq1", 3, 200, 275, "Phe", "GAA");
        hit3.forward_strand = true;

        assert!(!ScanResult::overlaps(&hit1, &hit3));
    }

    #[test]
    fn test_overlap_amount() {
        let mut hit1 = make_test_hit("seq1", 1, 100, 175, "Leu", "CAA");
        hit1.forward_strand = true;

        let mut hit2 = make_test_hit("seq1", 2, 150, 225, "Ser", "AGA");
        hit2.forward_strand = true;

        // Overlap: 150-175 = 26 bp
        assert_eq!(ScanResult::overlap_amount(&hit1, &hit2), 26);
    }

    #[test]
    fn test_resolve_overlaps() {
        let mut result = ScanResult::new();

        let mut hit1 = make_test_hit("seq1", 1, 100, 175, "Leu", "CAA");
        hit1.inf_score = 80.0;
        hit1.forward_strand = true;
        result.add_trna(hit1);

        let mut hit2 = make_test_hit("seq1", 2, 150, 225, "Ser", "AGA");
        hit2.inf_score = 50.0;
        hit2.forward_strand = true;
        result.add_trna(hit2);

        result.resolve_overlaps();

        // Should keep only the higher-scoring hit
        assert_eq!(result.trnas.len(), 1);
        assert_eq!(result.trnas[0].inf_score, 80.0);
    }

    #[test]
    fn test_merge() {
        let mut result1 = ScanResult::new();
        result1.add_trna(make_test_hit("seq1", 1, 100, 175, "Leu", "CAA"));

        let mut result2 = ScanResult::new();
        result2.add_trna(make_test_hit("seq2", 1, 200, 275, "Ser", "AGA"));

        result1.merge(result2);

        assert_eq!(result1.trnas.len(), 2);
    }

    #[test]
    fn test_iteration() {
        let mut result = ScanResult::new();
        result.add_trna(make_test_hit("seq1", 1, 100, 175, "Leu", "CAA"));
        result.add_pseudogene(make_test_hit("seq1", 2, 200, 275, "Ser", "AGA"));

        assert_eq!(result.iter().count(), 1);
        assert_eq!(result.iter_all().count(), 2);
        assert_eq!(result.iter_pseudogenes().count(), 1);
    }

    #[test]
    fn test_write_output() {
        let mut result = ScanResult::new();
        let mut hit = make_test_hit("chr1", 1, 1000, 1075, "Leu", "CAA");
        hit.inf_score = 74.5;
        hit.isotype_cm = "Leu".to_string();
        hit.isotype_score = 80.2;
        result.add_trna(hit);

        let mut output = Vec::new();
        result.write_output(&mut output).unwrap();

        let output_str = String::from_utf8(output).unwrap();
        assert!(output_str.contains("chr1"));
        assert!(output_str.contains("Leu"));
        assert!(output_str.contains("CAA"));
        assert!(output_str.contains("74.5"));
    }

    #[test]
    fn test_write_bed_output() {
        let mut result = ScanResult::new();
        let hit = make_test_hit("chr1", 1, 1000, 1075, "Leu", "CAA");
        result.add_trna(hit);

        let mut output = Vec::new();
        result.write_bed_output(&mut output).unwrap();

        let output_str = String::from_utf8(output).unwrap();
        assert!(output_str.contains("chr1"));
        assert!(output_str.contains("999")); // 0-based start
        assert!(output_str.contains("tRNA1-LeuCAA"));
    }

    #[test]
    fn test_write_gff_output() {
        let mut result = ScanResult::new();
        let hit = make_test_hit("chr1", 1, 1000, 1075, "Leu", "CAA");
        result.add_trna(hit);

        let mut output = Vec::new();
        result.write_gff_output(&mut output).unwrap();

        let output_str = String::from_utf8(output).unwrap();
        assert!(output_str.contains("##gff-version 3"));
        assert!(output_str.contains("tRNAscan-SE"));
        assert!(output_str.contains("tRNA"));
        assert!(output_str.contains("isotype=Leu"));
    }

    #[test]
    fn test_write_stats() {
        let mut result = ScanResult::new();
        result.set_source("test.fa".to_string());
        result.add_trna(make_test_hit("chr1", 1, 1000, 1075, "Leu", "CAA"));
        result.add_trna(make_test_hit("chr1", 2, 2000, 2075, "Leu", "TAG"));
        result.add_trna(make_test_hit("chr1", 3, 3000, 3075, "Ser", "AGA"));

        let mut output = Vec::new();
        result.write_stats(&mut output).unwrap();

        let output_str = String::from_utf8(output).unwrap();
        assert!(output_str.contains("Total tRNAs:"));
        assert!(output_str.contains("3"));
        assert!(output_str.contains("Leu"));
    }

    #[test]
    fn test_clear() {
        let mut result = ScanResult::new();
        result.add_trna(make_test_hit("seq1", 1, 100, 175, "Leu", "CAA"));
        result.add_pseudogene(make_test_hit("seq1", 2, 200, 275, "Ser", "AGA"));

        result.clear();

        assert!(result.is_empty());
        assert_eq!(result.trna_count(), 0);
        assert_eq!(result.pseudogene_count(), 0);
    }
}
