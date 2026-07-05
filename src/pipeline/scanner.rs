//! Main tRNA scanning pipeline for tRNAscan-SE.
//!
//! This module provides the core scanning functionality:
//! - Coordinate first-pass detection and CM-based analysis
//! - Process sequences from FASTA files
//! - Collect and organize tRNA hits

use std::collections::HashMap;
use std::path::Path;

use crate::eufind::pavesi::TrnaInfo;
use crate::isotype::{anticodon_to_isotype, IsotypeScorer};
use crate::squid::sqio::SeqFileReader;
use crate::types::cm::CM;
use crate::types::state::IState;

use super::config::TrnaScanConfig;

/// Origin of a tRNA hit (first-pass method that found it).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HitOrigin {
    /// Found by EuFindtRNA (Pavesi algorithm)
    EuFind,

    /// Found by tRNAscan 1.4
    TrnaScan14,

    /// Found by Infernal first-pass
    Infernal,

    /// Found by both EuFindtRNA and tRNAscan 1.4
    Both,

    /// Found by manual/custom method
    Manual,
}

impl HitOrigin {
    /// Get the display name for this origin.
    pub fn name(&self) -> &'static str {
        match self {
            HitOrigin::EuFind => "EuFind",
            HitOrigin::TrnaScan14 => "tRNAscan",
            HitOrigin::Infernal => "Inf",
            HitOrigin::Both => "Both",
            HitOrigin::Manual => "Manual",
        }
    }

    /// Get the short abbreviation for tabular output.
    pub fn short_name(&self) -> &'static str {
        match self {
            HitOrigin::EuFind => "Eu",
            HitOrigin::TrnaScan14 => "Ts",
            HitOrigin::Infernal => "Inf",
            HitOrigin::Both => "Bo",
            HitOrigin::Manual => "Ma",
        }
    }
}

impl std::fmt::Display for HitOrigin {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.name())
    }
}

/// A single tRNA hit from the scanning pipeline.
#[derive(Debug, Clone)]
pub struct TrnaHit {
    /// Sequence name where tRNA was found
    pub seq_name: String,

    /// tRNA number within sequence (1-based)
    pub trna_num: i32,

    /// Start position in sequence (1-based)
    pub start: i32,

    /// End position in sequence (1-based)
    pub end: i32,

    /// Predicted isotype (e.g., "Leu", "Ser", "SeC")
    pub isotype: String,

    /// Anticodon sequence (3 letters)
    pub anticodon: String,

    /// Intron start position (0 if no intron)
    pub intron_start: i32,

    /// Intron end position (0 if no intron)
    pub intron_end: i32,

    /// Infernal/CM bit score
    pub inf_score: f64,

    /// HMM/first-pass score
    pub hmm_score: f64,

    /// Secondary structure score
    pub ss_score: f64,

    /// Hit origin (which first-pass method found it)
    pub origin: HitOrigin,

    /// Isotype-specific CM model that scored best
    pub isotype_cm: String,

    /// Isotype-specific CM score
    pub isotype_score: f64,

    /// Note field (e.g., "pseudo" for pseudogenes)
    pub note: String,

    /// Whether this is a pseudogene
    pub is_pseudo: bool,

    /// tRNA sequence (with intron if present)
    pub sequence: String,

    /// Secondary structure string (Konings notation)
    pub secondary_structure: String,

    /// Position of anticodon in sequence (1-based start)
    pub anticodon_pos_start: i32,

    /// Position of anticodon in sequence (1-based end)
    pub anticodon_pos_end: i32,

    /// Strand (true = forward, false = reverse)
    pub forward_strand: bool,

    /// Full TrnaInfo from first-pass detection
    pub trna_info: Option<TrnaInfo>,
}

impl TrnaHit {
    /// Create a new tRNA hit with minimal information.
    pub fn new(
        seq_name: impl Into<String>,
        trna_num: i32,
        start: i32,
        end: i32,
        isotype: impl Into<String>,
        anticodon: impl Into<String>,
    ) -> Self {
        let forward_strand = end >= start;
        Self {
            seq_name: seq_name.into(),
            trna_num,
            start,
            end,
            isotype: isotype.into(),
            anticodon: anticodon.into(),
            intron_start: 0,
            intron_end: 0,
            inf_score: 0.0,
            hmm_score: 0.0,
            ss_score: 0.0,
            origin: HitOrigin::Infernal,
            isotype_cm: String::new(),
            isotype_score: 0.0,
            note: String::new(),
            is_pseudo: false,
            sequence: String::new(),
            secondary_structure: String::new(),
            anticodon_pos_start: 0,
            anticodon_pos_end: 0,
            forward_strand,
            trna_info: None,
        }
    }

    /// Get the tRNA length (excluding intron).
    pub fn length(&self) -> i32 {
        let raw_len = (self.end - self.start).abs() + 1;
        if self.intron_start > 0 && self.intron_end > 0 {
            raw_len - (self.intron_end - self.intron_start + 1)
        } else {
            raw_len
        }
    }

    /// Get the full tRNA length (including intron).
    pub fn full_length(&self) -> i32 {
        (self.end - self.start).abs() + 1
    }

    /// Check if this tRNA has an intron.
    pub fn has_intron(&self) -> bool {
        self.intron_start > 0 && self.intron_end > 0
    }

    /// Get BED-format coordinates (0-based, half-open).
    pub fn bed_coords(&self) -> (i32, i32) {
        if self.forward_strand {
            (self.start - 1, self.end)
        } else {
            (self.end - 1, self.start)
        }
    }

    /// Get strand character for output.
    pub fn strand_char(&self) -> char {
        if self.forward_strand {
            '+'
        } else {
            '-'
        }
    }

    /// Generate the tRNA identifier (e.g., "CELF22B7.tRNA1-LeuCAA").
    pub fn identifier(&self) -> String {
        format!(
            "{}.tRNA{}-{}{}",
            self.seq_name, self.trna_num, self.isotype, self.anticodon
        )
    }

    /// Check if isotype and anticodon match expected genetic code.
    pub fn isotype_matches_anticodon(&self) -> bool {
        if let Some(expected) = anticodon_to_isotype(&self.anticodon) {
            expected == self.isotype
        } else {
            false
        }
    }

    /// Set scores from IsotypeScorer.
    pub fn set_isotype_scores(&mut self, scorer: &IsotypeScorer) {
        self.isotype_cm = scorer.cm_best_isotype.clone();
        self.isotype_score = scorer.cm_best_score as f64;
        self.hmm_score = scorer.hmm_score as f64;
        self.ss_score = scorer.secondary_structure_score as f64;
    }
}

/// Results from scanning a complete file or set of sequences.
#[derive(Debug, Clone)]
pub struct ScanResults {
    /// All tRNA hits found
    pub hits: Vec<TrnaHit>,

    /// Number of sequences scanned
    pub sequences_scanned: usize,

    /// Total bases scanned
    pub bases_scanned: usize,

    /// Number of sequences with at least one hit
    pub sequences_with_hits: usize,

    /// Average tRNA length
    pub avg_trna_length: f64,

    /// Count of tRNAs with introns
    pub trnas_with_introns: usize,

    /// Isotype counts: isotype -> count
    pub isotype_counts: HashMap<String, usize>,

    /// Anticodon counts: anticodon -> count
    pub anticodon_counts: HashMap<String, usize>,

    /// Number of predicted pseudogenes
    pub pseudogene_count: usize,

    /// Number of tRNAs with isotype mismatches
    pub mismatch_count: usize,

    /// Number of selenocysteine tRNAs
    pub sec_count: usize,

    /// Number of suppressor tRNAs
    pub sup_count: usize,

    /// Scan start time (for statistics)
    pub start_time: Option<std::time::Instant>,

    /// First-pass scan time in seconds
    pub first_pass_time: f64,

    /// Infernal scan time in seconds
    pub infernal_time: f64,
}

impl ScanResults {
    /// Create new empty results.
    pub fn new() -> Self {
        Self {
            hits: Vec::new(),
            sequences_scanned: 0,
            bases_scanned: 0,
            sequences_with_hits: 0,
            avg_trna_length: 0.0,
            trnas_with_introns: 0,
            isotype_counts: HashMap::new(),
            anticodon_counts: HashMap::new(),
            pseudogene_count: 0,
            mismatch_count: 0,
            sec_count: 0,
            sup_count: 0,
            start_time: None,
            first_pass_time: 0.0,
            infernal_time: 0.0,
        }
    }

    /// Add a hit and update statistics.
    pub fn add_hit(&mut self, hit: TrnaHit) {
        // Update isotype count
        *self.isotype_counts.entry(hit.isotype.clone()).or_insert(0) += 1;

        // Update anticodon count
        *self
            .anticodon_counts
            .entry(hit.anticodon.clone())
            .or_insert(0) += 1;

        // Update special counts
        if hit.is_pseudo {
            self.pseudogene_count += 1;
        }
        if hit.isotype == "SeC" {
            self.sec_count += 1;
        }
        if hit.isotype == "Sup" {
            self.sup_count += 1;
        }
        if hit.has_intron() {
            self.trnas_with_introns += 1;
        }
        if !hit.isotype_matches_anticodon() {
            self.mismatch_count += 1;
        }

        self.hits.push(hit);
    }

    /// Finalize statistics after all hits added.
    pub fn finalize(&mut self) {
        if !self.hits.is_empty() {
            let total_len: i32 = self.hits.iter().map(|h| h.full_length()).sum();
            self.avg_trna_length = total_len as f64 / self.hits.len() as f64;

            // Count unique sequences with hits
            let seq_names: std::collections::HashSet<_> =
                self.hits.iter().map(|h| &h.seq_name).collect();
            self.sequences_with_hits = seq_names.len();
        }
    }

    /// Get total tRNA count (excluding pseudogenes).
    pub fn total_trnas(&self) -> usize {
        self.hits.iter().filter(|h| !h.is_pseudo).count()
    }

    /// Get count of standard 20 AA tRNAs.
    pub fn standard_aa_count(&self) -> usize {
        self.hits
            .iter()
            .filter(|h| {
                !h.is_pseudo
                    && h.isotype != "SeC"
                    && h.isotype != "Sup"
                    && h.isotype != "Unk"
                    && h.isotype != "iMet"
                    && h.isotype != "fMet"
            })
            .count()
    }

    /// Get scan speed in bp/sec.
    pub fn scan_speed(&self) -> f64 {
        let total_time = self.first_pass_time + self.infernal_time;
        if total_time > 0.0 {
            self.bases_scanned as f64 / total_time
        } else {
            0.0
        }
    }

    /// Get isotype-anticodon distribution for statistics.
    pub fn get_isotype_anticodon_table(&self) -> HashMap<String, HashMap<String, usize>> {
        let mut table: HashMap<String, HashMap<String, usize>> = HashMap::new();

        for hit in &self.hits {
            if hit.is_pseudo {
                continue;
            }
            let entry = table.entry(hit.isotype.clone()).or_insert_with(HashMap::new);
            *entry.entry(hit.anticodon.clone()).or_insert(0) += 1;
        }

        table
    }
}

impl Default for ScanResults {
    fn default() -> Self {
        Self::new()
    }
}

/// Main tRNA scanner.
#[derive(Debug)]
#[allow(dead_code)]
pub struct TrnaScanner {
    /// Scanner configuration
    config: TrnaScanConfig,

    /// Covariance model (loaded from file)
    cm: Option<CM>,

    /// Integer state model for Viterbi
    icm: Vec<IState>,
}

impl TrnaScanner {
    /// Create a new scanner with the given configuration.
    pub fn new(config: TrnaScanConfig) -> Self {
        Self {
            config,
            cm: None,
            icm: Vec::new(),
        }
    }

    /// Load covariance model from file.
    pub fn load_cm(&mut self, _path: &Path) -> Result<(), ScanError> {
        // Placeholder: actual CM loading would use save.rs
        // For now, create an empty CM
        self.cm = Some(CM::new(0));
        Ok(())
    }

    /// Scan a single sequence for tRNAs.
    pub fn scan_sequence(&self, seq: &str, seq_name: &str) -> Vec<TrnaHit> {
        let hits = Vec::new();

        // Skip very short sequences
        if seq.len() < 50 {
            return hits;
        }

        // This is a placeholder implementation
        // In production, this would:
        // 1. Run first-pass detection (EuFindtRNA, tRNAscan 1.4, or Infernal)
        // 2. Extract candidate regions
        // 3. Run Infernal CM search on candidates
        // 4. Score with isotype-specific models
        // 5. Determine anticodons and isotypes
        // 6. Build secondary structure

        // For now, just return empty (actual implementation would process sequences)
        let _ = (seq, seq_name);

        hits
    }

    /// Scan all sequences in a FASTA file.
    pub fn scan_file(&self, input_path: &Path) -> Result<ScanResults, ScanError> {
        let mut results = ScanResults::new();
        results.start_time = Some(std::time::Instant::now());

        // Open sequence file
        let mut reader = SeqFileReader::open(input_path).map_err(|e| {
            ScanError::FileError(input_path.to_path_buf(), e.to_string())
        })?;

        // Process each sequence
        let mut trna_num = 0;
        while let Some((seq_data, sqinfo)) = reader.read_seq().map_err(|e| {
            ScanError::FileError(input_path.to_path_buf(), e.to_string())
        })? {
            results.sequences_scanned += 1;
            results.bases_scanned += seq_data.len();

            // Get sequence name
            let seq_name = if !sqinfo.name.is_empty() {
                sqinfo.name.clone()
            } else {
                format!("seq{}", results.sequences_scanned)
            };

            // Convert sequence to string
            let seq_str = String::from_utf8_lossy(&seq_data);

            // Check match pattern if specified
            if let Some(ref pattern) = self.config.match_pattern {
                if !seq_name.contains(pattern) {
                    continue;
                }
            }

            // Scan sequence
            let seq_hits = self.scan_sequence(&seq_str, &seq_name);

            // Add hits with proper numbering
            for mut hit in seq_hits {
                trna_num += 1;
                hit.trna_num = trna_num;
                results.add_hit(hit);
            }

            // Also scan reverse complement if configured
            if self.config.both_strands {
                let rc_seq = reverse_complement(&seq_str);
                let rc_hits = self.scan_sequence(&rc_seq, &seq_name);

                for mut hit in rc_hits {
                    trna_num += 1;
                    hit.trna_num = trna_num;
                    hit.forward_strand = false;
                    // Adjust coordinates for reverse strand
                    let seq_len = seq_data.len() as i32;
                    let new_start = seq_len - hit.end + 1;
                    let new_end = seq_len - hit.start + 1;
                    hit.start = new_start;
                    hit.end = new_end;
                    results.add_hit(hit);
                }
            }
        }

        results.finalize();

        // Record timing
        if let Some(start) = results.start_time {
            let elapsed = start.elapsed().as_secs_f64();
            results.first_pass_time = elapsed * 0.1; // Estimate: 10% first pass
            results.infernal_time = elapsed * 0.9; // Estimate: 90% Infernal
        }

        Ok(results)
    }

    /// Get reference to configuration.
    pub fn config(&self) -> &TrnaScanConfig {
        &self.config
    }

    /// Get mutable reference to configuration.
    pub fn config_mut(&mut self) -> &mut TrnaScanConfig {
        &mut self.config
    }
}

/// Scanner errors.
#[derive(Debug, Clone)]
pub enum ScanError {
    /// File I/O error
    FileError(std::path::PathBuf, String),

    /// No sequences found in input
    NoSequences,

    /// Invalid sequence format
    InvalidFormat(String),

    /// CM model error
    CmError(String),

    /// Memory allocation error
    MemoryError(String),

    /// Interrupted by user
    Interrupted,
}

impl std::fmt::Display for ScanError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ScanError::FileError(path, msg) => {
                write!(f, "File error '{}': {}", path.display(), msg)
            }
            ScanError::NoSequences => {
                write!(f, "No sequences found in input file")
            }
            ScanError::InvalidFormat(msg) => {
                write!(f, "Invalid FASTA format: {}", msg)
            }
            ScanError::CmError(msg) => {
                write!(f, "CM model error: {}", msg)
            }
            ScanError::MemoryError(msg) => {
                write!(f, "Memory allocation failed: {}", msg)
            }
            ScanError::Interrupted => {
                write!(f, "Interrupted by user")
            }
        }
    }
}

impl std::error::Error for ScanError {}

/// Compute reverse complement of a DNA sequence.
fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c.to_ascii_uppercase() {
            'A' => 'T',
            'T' => 'A',
            'U' => 'A',
            'G' => 'C',
            'C' => 'G',
            'N' => 'N',
            _ => 'N',
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_trna_hit_creation() {
        let hit = TrnaHit::new("CELF22B7", 1, 12619, 12738, "Leu", "CAA");

        assert_eq!(hit.seq_name, "CELF22B7");
        assert_eq!(hit.trna_num, 1);
        assert_eq!(hit.start, 12619);
        assert_eq!(hit.end, 12738);
        assert_eq!(hit.isotype, "Leu");
        assert_eq!(hit.anticodon, "CAA");
        assert!(hit.forward_strand);
    }

    #[test]
    fn test_trna_hit_length() {
        let mut hit = TrnaHit::new("seq", 1, 100, 180, "Leu", "CAA");
        assert_eq!(hit.full_length(), 81);

        // With intron
        hit.intron_start = 120;
        hit.intron_end = 150;
        assert_eq!(hit.length(), 50); // 81 - 31 = 50
        assert!(hit.has_intron());
    }

    #[test]
    fn test_trna_hit_identifier() {
        let hit = TrnaHit::new("CELF22B7", 1, 100, 180, "Leu", "CAA");
        assert_eq!(hit.identifier(), "CELF22B7.tRNA1-LeuCAA");
    }

    #[test]
    fn test_trna_hit_bed_coords() {
        // Forward strand
        let hit = TrnaHit::new("seq", 1, 100, 180, "Leu", "CAA");
        let (start, end) = hit.bed_coords();
        assert_eq!(start, 99); // 0-based
        assert_eq!(end, 180); // Half-open

        // Reverse strand
        let mut hit = TrnaHit::new("seq", 1, 180, 100, "Leu", "CAA");
        hit.forward_strand = false;
        let (start, end) = hit.bed_coords();
        assert_eq!(start, 99);
        assert_eq!(end, 180);
    }

    #[test]
    fn test_scan_results() {
        let mut results = ScanResults::new();

        let hit1 = TrnaHit::new("seq1", 1, 100, 180, "Leu", "CAA");
        let hit2 = TrnaHit::new("seq1", 2, 500, 580, "Ser", "AGA");
        let hit3 = TrnaHit::new("seq2", 3, 200, 275, "Phe", "GAA");

        results.add_hit(hit1);
        results.add_hit(hit2);
        results.add_hit(hit3);
        results.finalize();

        assert_eq!(results.hits.len(), 3);
        assert_eq!(results.isotype_counts.get("Leu"), Some(&1));
        assert_eq!(results.isotype_counts.get("Ser"), Some(&1));
        assert_eq!(results.isotype_counts.get("Phe"), Some(&1));
        assert_eq!(results.sequences_with_hits, 2);
    }

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement("ATCG"), "CGAT");
        assert_eq!(reverse_complement("AAAA"), "TTTT");
        assert_eq!(reverse_complement("GCGC"), "GCGC");
    }

    #[test]
    fn test_hit_origin_display() {
        assert_eq!(HitOrigin::Infernal.name(), "Inf");
        assert_eq!(HitOrigin::EuFind.short_name(), "Eu");
        assert_eq!(HitOrigin::TrnaScan14.short_name(), "Ts");
    }

    #[test]
    fn test_scanner_creation() {
        let config = TrnaScanConfig::default();
        let scanner = TrnaScanner::new(config);
        assert!(scanner.cm.is_none());
        assert!(scanner.icm.is_empty());
    }

    #[test]
    fn test_isotype_anticodon_match() {
        let mut hit = TrnaHit::new("seq", 1, 100, 180, "Leu", "CAA");
        assert!(hit.isotype_matches_anticodon());

        hit.isotype = "Ser".to_string();
        assert!(!hit.isotype_matches_anticodon());
    }
}
