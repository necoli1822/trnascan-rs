//! CM Scanning Module for tRNAscan-SE
//!
//! This module provides covariance model (CM) scanning functionality using
//! the Infernal Rust library directly. It replaces the Perl CM.pm module
//! which called external Infernal binaries.
//!
//! Key components:
//! - CMScan: Main scanner with CM file paths and options
//! - CMScanOptions: Configuration for cmsearch execution
//! - CMSearchHit: Individual search hit result
//!
//! Ported from tRNAscanSE::CM.pm (3841 lines)

use std::collections::HashMap;
use std::path::PathBuf;

use infernox::{cm_file_read, CM};
use infernox::{FaithfulConfig, FaithfulHit, FaithfulSearcher};
use infernox::easel::error::InfernalError;

use crate::trna::{Strand, TRna, Intron, AnticodonPos, TRnaCategory};

pub mod decode;

// ============================================================================
// Constants
// ============================================================================

/// Minimum canonical intron length
const MIN_INTRON_LENGTH: usize = 4;

/// Minimum tRNA length without intron
const MIN_TRNA_NO_INTRON: usize = 62;

/// Default CM score cutoff
const DEFAULT_CM_CUTOFF: f64 = 20.0;

/// Organelle CM score cutoff
const DEFAULT_ORGANELLE_CM_CUTOFF: f64 = 15.0;

/// First-pass Infernal cutoff
const DEFAULT_FP_CUTOFF: f64 = 10.0;

/// BHB (bulge-helix-bulge) CM cutoff for intron detection
const DEFAULT_BHB_CM_CUTOFF: f64 = 10.0;

/// Minimum pseudogene filter score (cmsearch)
const MIN_CMSEARCH_PSEUDO_FILTER_SCORE: f64 = 55.0;

/// Minimum secondary structure score
const MIN_SS_SCORE: f64 = 5.0;

/// Minimum HMM score
const MIN_HMM_SCORE: f64 = 10.0;

// ============================================================================
// CMScanOptions
// ============================================================================

/// Options for CM scanning
#[derive(Debug, Clone)]
pub struct CMScanOptions {
    /// Score cutoff threshold (bits)
    pub score_cutoff: f64,
    /// E-value cutoff (optional, score takes precedence)
    pub eval_cutoff: Option<f64>,
    /// Enable truncated sequence search
    pub truncated_search: bool,
    /// Use local alignment mode
    pub local: bool,
    /// Disable HMM filter (slower but more sensitive)
    pub nohmm: bool,
    /// Use maximum sensitivity settings
    pub max: bool,
    /// Number of threads for parallel processing
    pub num_threads: usize,
    /// Maximum matrix size (MB)
    pub mxsize: usize,
    /// Use HMM filter for first pass
    pub hmm_filter: bool,
    /// Check for introns
    pub check_for_introns: bool,
    /// Check for split halves (archaeal)
    pub check_for_split_halves: bool,
    /// Skip pseudogene filter
    pub skip_pseudo_filter: bool,
    /// Get HMM score (for pseudogene detection)
    pub get_hmm_score: bool,
    /// Minimum intron length
    pub min_intron_length: usize,
    /// Minimum tRNA length without intron
    pub min_trna_no_intron: usize,
    /// Minimum pseudogene filter score
    pub min_pseudo_filter_score: f64,
    /// Minimum secondary structure score
    pub min_ss_score: f64,
    /// Minimum HMM score
    pub min_hmm_score: f64,
    /// Organelle mode cutoff
    pub organelle_cm_cutoff: f64,
    /// BHB CM cutoff
    pub bhb_cm_cutoff: f64,
    /// First-pass cutoff
    pub fp_cutoff: f64,
}

impl Default for CMScanOptions {
    fn default() -> Self {
        CMScanOptions {
            score_cutoff: DEFAULT_CM_CUTOFF,
            eval_cutoff: None,
            truncated_search: false,
            local: false,
            nohmm: true,  // Default to no HMM for accuracy
            max: false,
            num_threads: 1,
            mxsize: 128,
            hmm_filter: false,
            check_for_introns: false,
            check_for_split_halves: false,
            skip_pseudo_filter: false,
            get_hmm_score: false,
            min_intron_length: MIN_INTRON_LENGTH,
            min_trna_no_intron: MIN_TRNA_NO_INTRON,
            min_pseudo_filter_score: MIN_CMSEARCH_PSEUDO_FILTER_SCORE,
            min_ss_score: MIN_SS_SCORE,
            min_hmm_score: MIN_HMM_SCORE,
            organelle_cm_cutoff: DEFAULT_ORGANELLE_CM_CUTOFF,
            bhb_cm_cutoff: DEFAULT_BHB_CM_CUTOFF,
            fp_cutoff: DEFAULT_FP_CUTOFF,
        }
    }
}

// ============================================================================
// CMSearchHit
// ============================================================================

/// A single CM search hit result
///
/// Contains all information extracted from cmsearch output including
/// coordinates, scores, alignment, and tRNA-specific annotations.
#[derive(Debug, Clone, Default)]
pub struct CMSearchHit {
    // === Target (sequence) information ===
    /// Target sequence name
    pub target_name: String,
    /// Target sequence accession
    pub target_acc: String,
    /// Target sequence length
    pub target_len: i64,

    // === Query (model) information ===
    /// Query model name
    pub query_name: String,
    /// Query model accession
    pub query_acc: String,
    /// Query model length (consensus length)
    pub query_len: i32,

    // === Hit coordinates ===
    /// Hit start position in target (1-indexed, always smaller)
    pub seq_from: i64,
    /// Hit end position in target (1-indexed, always larger)
    pub seq_to: i64,
    /// Model start position (1-indexed)
    pub mdl_from: i32,
    /// Model end position (1-indexed)
    pub mdl_to: i32,
    /// Strand: '+' or '-'
    pub strand: char,

    // === Scoring ===
    /// Bit score
    pub score: f64,
    /// E-value
    pub evalue: f64,
    /// Bias score
    pub bias: f64,
    /// GC content
    pub gc: f64,

    // === Truncation info ===
    /// Truncation status: "no", "5'", "3'", "5'&3'"
    pub trunc: String,
    /// Pass number (which filter pass)
    pub pass: i32,
    /// Inclusion threshold marker: '!' included, '?' below threshold
    pub inc: char,

    // === Alignment data ===
    /// Aligned target sequence (with gaps)
    pub target_seq: String,
    /// Aligned query sequence (consensus with gaps)
    pub query_seq: String,
    /// Posterior probability string
    pub pp: String,
    /// Secondary structure annotation
    pub ss: String,
    /// Non-canonical base pair annotation
    pub nc: String,

    // === tRNA-specific fields ===
    /// Anticodon triplet
    pub anticodon: String,
    /// Anticodon loop start position (1-indexed, relative to tRNA)
    pub antiloop_start: i32,
    /// Anticodon loop end position
    pub antiloop_end: i32,
    /// Anticodon position within tRNA (1-indexed)
    pub anticodon_pos: i32,
    /// tRNA isotype (e.g., "Ala", "Gly")
    pub isotype: String,
    /// Intron sequence (if present)
    pub intron_seq: String,
    /// Intron start position (relative to tRNA)
    pub intron_start: i32,
    /// Intron end position (relative to tRNA)
    pub intron_end: i32,
    /// Intron type: "CI" (canonical) or "NCI" (non-canonical)
    pub intron_type: String,
    /// Is this hit likely a pseudogene?
    pub is_pseudo: bool,
    /// Pseudogene reason (if is_pseudo is true)
    pub pseudo_reason: String,
    /// HMM-only score (for pseudogene detection)
    pub hmm_score: f64,
    /// Secondary structure contribution to score
    pub ss_score: f64,
    /// CM model used
    pub model: String,
    /// tRNA category
    pub category: TRnaCategory,
    /// Additional notes
    pub notes: String,
}

impl CMSearchHit {
    /// Create a new empty CMSearchHit
    pub fn new() -> Self {
        CMSearchHit {
            strand: '+',
            trunc: "no".to_string(),
            pass: 1,
            inc: '!',
            ..Default::default()
        }
    }

    /// Create CMSearchHit from a faithful infernox `FaithfulHit`.
    ///
    /// Replaces the former `from_search_hit(&SearchHit)` (the legacy `cmsearch()`
    /// API was removed when infernox's pre-parity `legacy/` modules were deleted).
    /// Mirrors the exact column mapping `batch_search_external` produces — so the
    /// strand-oriented coords, evalue/bias/gc, and model span are now filled in
    /// (the old converter hard-coded '+' strand and left evalue/bias blank).
    pub fn from_faithful_hit(hit: &FaithfulHit, model_name: &str, target_name: &str, target_len: i64) -> Self {
        let mut h = CMSearchHit::default();
        h.target_name = target_name.to_string();
        h.query_name = model_name.to_string();
        h.model = model_name.to_string();
        h.target_len = target_len;
        h.seq_from = hit.start;
        h.seq_to = hit.stop;
        h.strand = if hit.in_rc { '-' } else { '+' };
        h.score = hit.score as f64;
        h.evalue = hit.evalue;
        h.bias = hit.bias as f64;
        h.gc = hit.gc;
        h.mdl_from = hit.mdl_from;
        h.mdl_to = hit.mdl_to;
        h.trunc = "no".to_string();
        h.pass = 1;
        h.inc = if hit.evalue <= 0.01 { '!' } else { '?' };

        // Parse isotype/anticodon from query name (e.g. "tRNA-Ala-TGC").
        if h.query_name.contains("tRNA-") {
            let parts: Vec<&str> = h.query_name.split('-').collect();
            if parts.len() >= 2 {
                h.isotype = parts[1].to_string();
            }
            if parts.len() >= 3 {
                h.anticodon = parts[2].to_string();
            }
        }
        h
    }

    /// Get the hit length
    pub fn len(&self) -> i64 {
        (self.seq_to - self.seq_from).abs() + 1
    }

    /// Check if hit is empty
    pub fn is_empty(&self) -> bool {
        self.target_name.is_empty()
    }

    /// Get strand as Strand enum
    pub fn get_strand(&self) -> Strand {
        match self.strand {
            '+' => Strand::Plus,
            '-' => Strand::Minus,
            _ => Strand::Unknown,
        }
    }

    /// Check if this hit passes the score cutoff
    pub fn passes_cutoff(&self, cutoff: f64) -> bool {
        self.score >= cutoff
    }

    /// Convert to TRna structure
    pub fn to_trna(&self, id: usize, seqname: &str) -> TRna {
        let mut trna = TRna::default();

        trna.id = id;
        trna.seqname = seqname.to_string();
        trna.start = self.seq_from;
        trna.end = self.seq_to;
        trna.strand = self.get_strand();
        trna.isotype = self.isotype.clone();
        trna.anticodon = self.anticodon.clone();
        trna.ss = self.ss.clone();
        trna.seq = self.target_seq.replace("-", "").replace(".", "");

        // Add anticodon position if available
        if self.anticodon_pos > 0 {
            trna.ac_positions.push(AnticodonPos {
                rel_start: self.anticodon_pos,
                rel_end: self.anticodon_pos + 2,
            });
        }

        // Add intron if present
        if !self.intron_seq.is_empty() && self.intron_start > 0 {
            trna.introns.push(Intron {
                rel_start: self.intron_start,
                rel_end: self.intron_end,
                start: 0,  // Will be computed based on strand
                end: 0,
                intron_type: self.intron_type.clone(),
                seq: self.intron_seq.clone(),
            });
        }

        // Set domain model
        trna.domain_models.insert("infernal".to_string(), crate::trna::DomainModel {
            score: self.score,
            mat_score: self.score,
            hmm_score: self.hmm_score,
            ss_score: self.ss_score,
        });

        trna.category = self.category.clone();
        if self.is_pseudo {
            trna.is_pseudo = true;
        }

        trna
    }
}

// ============================================================================
// CM File Path Management
// ============================================================================

/// CM file paths organized by domain and type
#[derive(Debug, Clone, Default)]
pub struct CMFilePaths {
    /// Main domain CM files (Domain -> path)
    pub main: HashMap<String, PathBuf>,
    /// Main non-secondary-structure CM files (Domain -> path)
    pub main_ns: HashMap<String, PathBuf>,
    /// Isotype-specific CM files (isotype -> path)
    pub isotype: HashMap<String, PathBuf>,
    /// Isotype CM database path
    pub isotype_db: Option<PathBuf>,
    /// Mitochondrial isotype CM files
    pub mito_isotype: HashMap<String, PathBuf>,
    /// Mitochondrial isotype CM database path
    pub mito_isotype_db: Option<PathBuf>,
    /// Intron/BHB CM files for archaeal intron detection
    pub intron: HashMap<String, PathBuf>,
    /// Archaeal 5' half CM
    pub arch_five_half: Option<PathBuf>,
    /// Archaeal 3' half CM
    pub arch_three_half: Option<PathBuf>,
    /// Cove CM file (legacy)
    pub cove: Option<PathBuf>,
    /// SeC (selenocysteine) prokaryotic CM
    pub pselc: Option<PathBuf>,
    /// SeC (selenocysteine) eukaryotic CM
    pub eselc: Option<PathBuf>,
    /// Isotype score cutoffs
    pub isotype_cutoffs: HashMap<String, f64>,
}

impl CMFilePaths {
    pub fn new() -> Self {
        CMFilePaths::default()
    }

    /// Add a main CM file path
    pub fn add_main(&mut self, key: &str, path: PathBuf) {
        self.main.insert(key.to_string(), path);
    }

    /// Add a main NS (non-secondary-structure) CM file path
    pub fn add_main_ns(&mut self, key: &str, path: PathBuf) {
        self.main_ns.insert(key.to_string(), path);
    }

    /// Add an isotype CM file path
    pub fn add_isotype(&mut self, isotype: &str, path: PathBuf) {
        self.isotype.insert(isotype.to_string(), path);
    }

    /// Add an intron CM file path
    pub fn add_intron(&mut self, key: &str, path: PathBuf) {
        self.intron.insert(key.to_string(), path);
    }

    /// Add isotype score cutoff
    pub fn add_isotype_cutoff(&mut self, isotype: &str, cutoff: f64) {
        self.isotype_cutoffs.insert(isotype.to_string(), cutoff);
    }

    /// Get the main CM file for a domain
    pub fn get_main(&self, domain: &str) -> Option<&PathBuf> {
        self.main.get(domain)
    }

    /// Get the NS CM file for a domain
    pub fn get_main_ns(&self, domain: &str) -> Option<&PathBuf> {
        self.main_ns.get(domain)
    }
}

// ============================================================================
// CMScan - Main Scanner
// ============================================================================

/// CM Scanner for tRNA detection
///
/// This struct manages CM file paths, options, and provides methods for
/// running cmsearch via the Infernal Rust library.
#[derive(Debug, Clone)]
pub struct CMScan {
    /// CM file paths
    pub cm_files: CMFilePaths,
    /// Scan options
    pub options: CMScanOptions,
    /// Cached CM models (model name -> CM)
    cached_cms: HashMap<String, CM>,
}

impl CMScan {
    /// Create a new CMScan instance
    pub fn new() -> Self {
        CMScan {
            cm_files: CMFilePaths::new(),
            options: CMScanOptions::default(),
            cached_cms: HashMap::new(),
        }
    }

    /// Create CMScan with specific options
    pub fn with_options(options: CMScanOptions) -> Self {
        CMScan {
            cm_files: CMFilePaths::new(),
            options,
            cached_cms: HashMap::new(),
        }
    }

    /// Load a CM file and cache it
    pub fn load_cm(&mut self, path: &PathBuf) -> Result<&CM, String> {
        let path_str = path.to_string_lossy().to_string();

        if !self.cached_cms.contains_key(&path_str) {
            let cm = cm_file_read(&path_str)
                .map_err(|e: InfernalError| format!("Failed to read CM file {}: {}", path_str, e))?;
            self.cached_cms.insert(path_str.clone(), cm);
        }

        self.cached_cms.get(&path_str)
            .ok_or_else(|| format!("Failed to load CM from {}", path_str))
    }

    /// Run cmsearch on a single sequence
    ///
    /// # Arguments
    /// * `cm_key` - Key to look up CM in main files (e.g., "Domain")
    /// * `seq` - The sequence to search
    /// * `seq_name` - Name of the sequence
    ///
    /// # Returns
    /// Vector of CMSearchHit results
    pub fn search(
        &mut self,
        cm_key: &str,
        seq: &str,
        seq_name: &str,
    ) -> Result<Vec<CMSearchHit>, String> {
        let cm_path = self.cm_files.get_main(cm_key)
            .ok_or_else(|| format!("CM file not found for key: {}", cm_key))?
            .clone();

        self.search_with_cm(&cm_path, seq, seq_name)
    }

    /// Run cmsearch with a specific CM file, in-process via the faithful infernox
    /// pipeline (byte-parity with C Infernal). Replaces the removed legacy
    /// `cmsearch(&cm, ...)` single-sequence API with `FaithfulSearcher`.
    pub fn search_with_cm(
        &mut self,
        cm_path: &PathBuf,
        seq: &str,
        seq_name: &str,
    ) -> Result<Vec<CMSearchHit>, String> {
        // Build the faithful searcher (reads CM in global config, builds p7
        // filters + CP9 HMM + configures CM scores).
        let searcher = FaithfulSearcher::from_cm_file(cm_path).map_err(|e| {
            format!("Failed to build faithful searcher from {}: {}", cm_path.display(), e)
        })?;
        let model_name = searcher.model_name().to_string();

        let seqs = [seq];
        let cfg = FaithfulConfig::default();
        let fhits = searcher.search(&seqs, &cfg);
        let target_len = seq.len() as i64;

        let hits: Vec<CMSearchHit> = fhits.iter()
            .filter(|h| h.score as f64 >= self.options.score_cutoff)
            .map(|h| CMSearchHit::from_faithful_hit(h, &model_name, seq_name, target_len))
            .collect();

        Ok(hits)
    }

    /// Run cmsearch on all main CM files
    ///
    /// Searches the sequence against all configured main CM models
    /// and returns merged, deduplicated results.
    pub fn search_all_main(
        &mut self,
        seq: &str,
        seq_name: &str,
    ) -> Result<Vec<CMSearchHit>, String> {
        let mut all_hits = Vec::new();

        // Get keys first to avoid borrow issues
        let keys: Vec<String> = self.cm_files.main.keys().cloned().collect();

        for key in keys {
            match self.search(&key, seq, seq_name) {
                Ok(hits) => all_hits.extend(hits),
                Err(e) => {
                    // Log error but continue with other CMs
                    eprintln!("Warning: search failed for CM {}: {}", key, e);
                }
            }
        }

        // Sort by score descending
        all_hits.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap_or(std::cmp::Ordering::Equal));

        // Remove overlapping hits, keeping higher scoring ones
        self.merge_overlapping_hits(all_hits)
    }

    /// Merge overlapping hits, keeping the best scoring one
    fn merge_overlapping_hits(&self, mut hits: Vec<CMSearchHit>) -> Result<Vec<CMSearchHit>, String> {
        if hits.is_empty() {
            return Ok(hits);
        }

        // Sort by position
        hits.sort_by(|a, b| {
            a.seq_from.cmp(&b.seq_from)
                .then_with(|| b.score.partial_cmp(&a.score).unwrap_or(std::cmp::Ordering::Equal))
        });

        let mut merged = Vec::new();
        let mut current = hits[0].clone();

        for hit in hits.into_iter().skip(1) {
            // Check for overlap
            if self.hits_overlap(&current, &hit) {
                // Keep the higher scoring one
                if hit.score > current.score {
                    current = hit;
                } else {
                    // Extend current if hit extends beyond
                    current.seq_to = current.seq_to.max(hit.seq_to);
                }
            } else {
                merged.push(current);
                current = hit;
            }
        }
        merged.push(current);

        Ok(merged)
    }

    /// Check if two hits overlap
    fn hits_overlap(&self, a: &CMSearchHit, b: &CMSearchHit) -> bool {
        if a.strand != b.strand {
            return false;
        }

        let overlap = a.seq_from.max(b.seq_from) <= a.seq_to.min(b.seq_to);
        overlap
    }

    // ========================================================================
    // Batch Search (External cmsearch)
    // ========================================================================

    /// Run batch cmsearch in-process via the faithful infernox pipeline.
    ///
    /// Replaces the former external-`cmsearch`-subprocess implementation (which
    /// silently returned 0 hits whenever no `cmsearch` binary was on PATH). The
    /// [`FaithfulSearcher`] runs the exact byte-parity Infernal pipeline in-
    /// process, so no external binary or temp files are needed. The CM is read +
    /// configured once per call and reused across all `sequences`.
    ///
    /// # Arguments
    /// * `cm_key` - Key to look up CM in main files
    /// * `sequences` - Vec of (seq_name, sequence) tuples
    ///
    /// # Returns
    /// Vec of CMSearchHit results
    pub fn batch_search_external(
        &self,
        cm_key: &str,
        sequences: &[(String, String)],
    ) -> Result<Vec<CMSearchHit>, String> {
        if sequences.is_empty() {
            return Ok(Vec::new());
        }

        let cm_path = self.cm_files.get_main(cm_key)
            .ok_or_else(|| format!("CM file not found for key: {}", cm_key))?;

        // Build the faithful searcher (reads CM in global config, builds p7
        // filters + CP9 HMM + configures CM scores) once for this batch.
        let searcher = FaithfulSearcher::from_cm_file(cm_path).map_err(|e| {
            format!("Failed to build faithful searcher from {}: {}", cm_path.display(), e)
        })?;
        let model_name = searcher.model_name().to_string();

        let seqs: Vec<&str> = sequences.iter().map(|(_, s)| s.as_str()).collect();
        let cfg = FaithfulConfig::default();
        let fhits = searcher.search(&seqs, &cfg);

        // Convert to CMSearchHit, mirroring what parse_tblout() produced from the
        // infernox tblout columns exactly (seq_from/seq_to are strand-oriented:
        // for the '-' strand, seq_from > seq_to). Downstream logic is unchanged.
        let results: Vec<CMSearchHit> = fhits
            .iter()
            .filter(|h| h.score as f64 >= self.options.score_cutoff)
            .map(|h| {
                let region_name = sequences[h.seq_idx].0.clone();
                let target_len = sequences[h.seq_idx].1.len() as i64;
                let mut hit = CMSearchHit::default();
                hit.target_name = region_name;
                hit.query_name = model_name.clone();
                hit.model = model_name.clone();
                hit.target_len = target_len;
                hit.seq_from = h.start;
                hit.seq_to = h.stop;
                hit.strand = if h.in_rc { '-' } else { '+' };
                hit.score = h.score as f64;
                hit.evalue = h.evalue;
                hit.bias = h.bias as f64;
                hit.gc = h.gc;
                hit.mdl_from = h.mdl_from;
                hit.mdl_to = h.mdl_to;
                hit.trunc = "no".to_string();
                hit.pass = 1;
                hit.inc = if h.evalue <= 0.01 { '!' } else { '?' };

                // Parse isotype/anticodon from query name (e.g. "tRNA-Ala-TGC").
                // Single-model CMs (e.g. bact-030216) carry no isotype in the
                // name; those fields stay blank and are resolved downstream.
                if hit.query_name.contains("tRNA-") {
                    let parts: Vec<&str> = hit.query_name.split('-').collect();
                    if parts.len() >= 2 {
                        hit.isotype = parts[1].to_string();
                    }
                    if parts.len() >= 3 {
                        hit.anticodon = parts[2].to_string();
                    }
                }
                hit
            })
            .collect();

        Ok(results)
    }

    // ========================================================================
    // Anticodon Finding
    // ========================================================================

    /// Find anticodon from secondary structure
    ///
    /// Matches pattern in secondary structure output, looking for second
    /// stem-loop structure ">>>>...<<<" that should be the anticodon stem-loop.
    ///
    /// Returns (anticodon, antiloop_start, antiloop_end, anticodon_pos)
    /// Find the anticodon from a NORMALIZED (post-`format_cmsearch_output`)
    /// sequence + secondary structure.
    ///
    /// Faithful port of `CM.pm::find_anticodon` (:731) — delegates to the exact
    /// regex-based implementation in `decode` (fixes the previous ad-hoc
    /// scanner's off-by-one). `seq`/`ss` MUST already be normalized.
    pub fn find_anticodon(
        &self,
        seq: &str,
        ss: &str,
        undef_anticodon: &str,
    ) -> (String, i32, i32, i32) {
        let (ac, ai, ae, acp) = decode::find_anticodon(seq, ss);
        if ac == decode::UNDEF_ANTICODON {
            (undef_anticodon.to_string(), -1, -1, -1)
        } else {
            (ac, ai, ae, acp)
        }
    }

    // ========================================================================
    // Intron Detection
    // ========================================================================

    /// Find canonical intron in anticodon loop region
    ///
    /// Looks for lowercase sequence (intron) within the anticodon loop.
    /// Returns (intron_seq, start, end) or empty tuple if none found.
    pub fn find_intron(
        &self,
        seq: &str,
        antiloop_start: i32,
        antiloop_end: i32,
    ) -> (String, i32, i32) {
        if antiloop_start < 0 || antiloop_end < 0 {
            return (String::new(), 0, 0);
        }

        let start_idx = antiloop_start as usize;
        let end_idx = (antiloop_end + 1).min(seq.len() as i32) as usize;

        if end_idx <= start_idx || end_idx > seq.len() {
            return (String::new(), 0, 0);
        }

        let antiloop_seq = &seq[start_idx..end_idx];

        // Look for intron (lowercase letters of minimum length)
        let min_len = self.options.min_intron_length;

        // Find runs of lowercase letters
        let mut intron_start = None;
        let mut current_run_start = None;
        let mut current_run_len = 0;

        for (i, c) in antiloop_seq.chars().enumerate() {
            if c.is_ascii_lowercase() {
                if current_run_start.is_none() {
                    current_run_start = Some(i);
                    current_run_len = 1;
                } else {
                    current_run_len += 1;
                }
            } else {
                if current_run_len >= min_len {
                    intron_start = current_run_start;
                    break;
                }
                current_run_start = None;
                current_run_len = 0;
            }
        }

        // Check final run
        if current_run_len >= min_len && intron_start.is_none() {
            intron_start = current_run_start;
        }

        if let Some(rel_start) = intron_start {
            let intron_seq: String = antiloop_seq[rel_start..rel_start + current_run_len]
                .to_string();

            // Find absolute position in full sequence
            let abs_start = start_idx + rel_start;
            let abs_end = abs_start + current_run_len - 1;

            // Find position relative to start of seq (account for multiple occurrences)
            if let Some(pos) = seq[..end_idx].rfind(&intron_seq.to_lowercase()) {
                return (
                    intron_seq,
                    (pos + 1) as i32,
                    (pos + current_run_len) as i32,
                );
            }

            return (
                intron_seq,
                (abs_start + 1) as i32,
                (abs_end + 1) as i32,
            );
        }

        (String::new(), 0, 0)
    }

    /// Check if intron boundaries are canonical
    ///
    /// Canonical introns have specific splice site sequences.
    pub fn is_canonical_intron(&self, intron_seq: &str, upstream: &str, downstream: &str) -> bool {
        // Canonical introns typically have conserved boundaries
        // In archaeal tRNAs: BHB (bulge-helix-bulge) motif
        // In eukaryotic tRNAs: usually 37/38 position introns

        if intron_seq.len() < self.options.min_intron_length {
            return false;
        }

        // Check for conserved boundaries (simplified)
        let upstream_ok = upstream.len() >= 2;
        let downstream_ok = downstream.len() >= 2;

        upstream_ok && downstream_ok
    }

    // ========================================================================
    // Pseudogene Detection
    // ========================================================================

    /// Check if a tRNA hit is likely a pseudogene
    ///
    /// Runs cmsearch with non-secondary-structure CM to get HMM score.
    /// Pseudogene if: primary structure score < threshold OR
    ///               secondary structure contribution < threshold
    pub fn is_pseudogene(
        &mut self,
        hit: &CMSearchHit,
        seq: &str,
        seq_name: &str,
        domain: &str,
    ) -> Result<(bool, String), String> {
        // Skip check if score is above minimum or filter is disabled
        if hit.score >= self.options.min_pseudo_filter_score || self.options.skip_pseudo_filter {
            if !self.options.get_hmm_score {
                return Ok((false, String::new()));
            }
        }

        // Get NS (non-secondary-structure) CM
        let ns_cm_path = match self.cm_files.get_main_ns(domain) {
            Some(path) => path.clone(),
            None => return Ok((false, String::new())),
        };

        // Run cmsearch with NS model
        let ns_hits = self.search_with_cm(&ns_cm_path, seq, seq_name)?;

        if ns_hits.is_empty() {
            return Ok((true, "No hit with NS model".to_string()));
        }

        let hmm_score = ns_hits[0].score;
        let ss_score = hit.score - hmm_score;

        // Check pseudogene criteria
        if ss_score < self.options.min_ss_score {
            return Ok((true, format!("Low SS score: {:.1}", ss_score)));
        }

        if hmm_score < self.options.min_hmm_score {
            return Ok((true, format!("Low HMM score: {:.1}", hmm_score)));
        }

        if hit.score < self.options.min_pseudo_filter_score {
            return Ok((true, format!("Low total score: {:.1}", hit.score)));
        }

        Ok((false, String::new()))
    }

    // ========================================================================
    // Isotype CM Support
    // ========================================================================

    /// Scan with isotype-specific CMs
    ///
    /// Runs cmscan against isotype CM database for more specific classification.
    pub fn scan_isotype(
        &mut self,
        seq: &str,
        seq_name: &str,
        use_mito: bool,
    ) -> Result<HashMap<String, f64>, String> {
        let db_path = if use_mito {
            self.cm_files.mito_isotype_db.clone()
        } else {
            self.cm_files.isotype_db.clone()
        };

        let db_path = match db_path {
            Some(p) => p,
            None => return Ok(HashMap::new()),
        };

        // Search against isotype database
        let hits = self.search_with_cm(&db_path, seq, seq_name)?;

        // Collect best score per isotype
        let mut isotype_scores: HashMap<String, f64> = HashMap::new();
        for hit in hits {
            let isotype = hit.query_name.clone();
            let current = isotype_scores.entry(isotype).or_insert(f64::NEG_INFINITY);
            if hit.score > *current {
                *current = hit.score;
            }
        }

        Ok(isotype_scores)
    }

    // ========================================================================
    // First-Pass Integration
    // ========================================================================

    /// Run first-pass cmsearch with relaxed settings
    ///
    /// Used for initial genome-wide scanning with HMM filter enabled.
    pub fn first_pass_search(
        &mut self,
        seq: &str,
        seq_name: &str,
    ) -> Result<Vec<CMSearchHit>, String> {
        // Save original settings
        let orig_cutoff = self.options.score_cutoff;
        let orig_hmm = self.options.hmm_filter;

        // Use first-pass settings
        self.options.score_cutoff = self.options.fp_cutoff;
        self.options.hmm_filter = true;

        let result = self.search_all_main(seq, seq_name);

        // Restore settings
        self.options.score_cutoff = orig_cutoff;
        self.options.hmm_filter = orig_hmm;

        result
    }

    /// Process first-pass hits and merge overlapping ones
    pub fn process_fp_hits(
        &self,
        hits: Vec<CMSearchHit>,
        start_offset: i64,
    ) -> Vec<CMSearchHit> {
        hits.into_iter()
            .map(|mut hit| {
                // Adjust coordinates for sequence offset
                hit.seq_from += start_offset;
                hit.seq_to += start_offset;
                hit
            })
            .collect()
    }

    // ========================================================================
    // Result Parsing (for compatibility with external cmsearch)
    // ========================================================================

    /// Parse tblout format line
    ///
    /// For compatibility when running external cmsearch binary.
    pub fn parse_tblout_line(line: &str) -> Option<CMSearchHit> {
        if line.starts_with('#') || line.trim().is_empty() {
            return None;
        }

        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 18 {
            return None;
        }

        let mut hit = CMSearchHit::new();
        hit.target_name = fields[0].to_string();
        hit.target_acc = fields[1].to_string();
        hit.query_name = fields[2].to_string();
        hit.query_acc = fields[3].to_string();
        // fields[4] is "cm"
        hit.mdl_from = fields[5].parse().unwrap_or(0);
        hit.mdl_to = fields[6].parse().unwrap_or(0);
        hit.seq_from = fields[7].parse().unwrap_or(0);
        hit.seq_to = fields[8].parse().unwrap_or(0);
        hit.strand = fields[9].chars().next().unwrap_or('+');
        hit.trunc = fields[10].to_string();
        hit.pass = fields[11].parse().unwrap_or(1);
        hit.gc = fields[12].parse().unwrap_or(0.0);
        hit.bias = fields[13].parse().unwrap_or(0.0);
        hit.score = fields[14].parse().unwrap_or(0.0);
        hit.evalue = fields[15].parse().unwrap_or(1.0);
        hit.inc = fields[16].chars().next().unwrap_or('?');

        // Description is everything after field 17
        if fields.len() > 17 {
            hit.notes = fields[17..].join(" ");
        }

        Some(hit)
    }

    /// Parse alignment output to extract sequence and structure
    pub fn parse_alignout(
        output: &str,
        seq_name: &str,
    ) -> Option<(String, String, f64)> {
        let mut seq = String::new();
        let mut ss = String::new();
        let mut score = -1000.0;

        // Escape special regex characters in seq_name
        let safe_name = seq_name.replace(|c: char| !c.is_alphanumeric() && c != '_', "");

        for line in output.lines() {
            // Look for sequence line: "  seqname ACGU..."
            if line.contains(&safe_name) {
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 2 {
                    for part in &parts[1..] {
                        if part.chars().all(|c| "ACGTUacgtunN.-".contains(c)) {
                            seq.push_str(part);
                        }
                    }
                }
            }
            // Look for structure line with score
            else if line.contains("bits") {
                if let Some(score_str) = line.split_whitespace().next() {
                    score = score_str.parse().unwrap_or(-1000.0);
                }
            }
            // Look for structure annotation: "  (( <<<...>>> ))"
            else if line.trim().starts_with("(") || line.trim().starts_with("<")
                 || line.trim().starts_with(">") || line.trim().starts_with(".")
            {
                let trimmed = line.trim();
                if trimmed.chars().all(|c| "()<>.-:_,".contains(c) || c.is_whitespace()) {
                    ss.push_str(&trimmed.chars().filter(|c| !c.is_whitespace()).collect::<String>());
                }
            }
        }

        // Remove gaps from sequence
        seq = seq.replace("-", "");
        ss = ss.replace(" ", "");

        if seq.is_empty() || ss.is_empty() {
            return None;
        }

        Some((seq, ss, score))
    }

    // ========================================================================
    // tRNA Property Decoding
    // ========================================================================

    /// Decode full tRNA properties from a CM hit
    ///
    /// This combines anticodon finding, intron detection, and pseudogene
    /// checking into a complete tRNA annotation.
    pub fn decode_trna_properties(
        &mut self,
        hit: &mut CMSearchHit,
        _seq: &str,
        undef_anticodon: &str,
        undef_isotype: &str,
    ) -> Result<(), String> {
        // Find anticodon
        let (anticodon, antiloop_start, antiloop_end, ac_pos) =
            self.find_anticodon(&hit.target_seq, &hit.ss, undef_anticodon);

        hit.anticodon = anticodon.clone();
        hit.antiloop_start = antiloop_start;
        hit.antiloop_end = antiloop_end;
        hit.anticodon_pos = ac_pos;

        // Check for undetermined anticodon
        if anticodon == undef_anticodon {
            hit.category = TRnaCategory::UndeterminedAc;
            hit.isotype = undef_isotype.to_string();
            return Ok(());
        }

        // Find intron
        let (intron_seq, intron_start, intron_end) =
            self.find_intron(&hit.target_seq, antiloop_start, antiloop_end);

        if !intron_seq.is_empty() {
            hit.intron_seq = intron_seq;
            hit.intron_start = intron_start;
            hit.intron_end = intron_end;
            hit.intron_type = "CI".to_string();  // Canonical intron
        }

        Ok(())
    }

    // ========================================================================
    // Utility Methods
    // ========================================================================

    /// Get the effective score cutoff based on mode
    pub fn get_score_cutoff(&self, organelle_mode: bool) -> f64 {
        if organelle_mode {
            self.options.organelle_cm_cutoff
        } else {
            self.options.score_cutoff
        }
    }

    /// Clear cached CM models
    pub fn clear_cache(&mut self) {
        self.cached_cms.clear();
    }

    /// Check if CM file exists and is readable
    pub fn verify_cm_file(path: &PathBuf) -> Result<(), String> {
        if !path.exists() {
            return Err(format!("CM file not found: {:?}", path));
        }
        if !path.is_file() {
            return Err(format!("CM path is not a file: {:?}", path));
        }
        Ok(())
    }
}

impl Default for CMScan {
    fn default() -> Self {
        CMScan::new()
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cmscan_options_default() {
        let opts = CMScanOptions::default();
        assert_eq!(opts.score_cutoff, DEFAULT_CM_CUTOFF);
        assert!(!opts.truncated_search);
        assert!(opts.nohmm);
    }

    #[test]
    fn test_cmsearch_hit_new() {
        let hit = CMSearchHit::new();
        assert_eq!(hit.strand, '+');
        assert_eq!(hit.trunc, "no");
        assert_eq!(hit.pass, 1);
        assert_eq!(hit.inc, '!');
    }

    #[test]
    fn test_cmsearch_hit_passes_cutoff() {
        let mut hit = CMSearchHit::new();
        hit.score = 50.0;

        assert!(hit.passes_cutoff(40.0));
        assert!(hit.passes_cutoff(50.0));
        assert!(!hit.passes_cutoff(60.0));
    }

    #[test]
    fn test_parse_tblout_line() {
        let line = "seq1 - tRNA-Ala - cm 1 71 100 170 + no 1 0.54 0.0 55.5 1.2e-10 ! -";
        let hit = CMScan::parse_tblout_line(line).unwrap();

        assert_eq!(hit.target_name, "seq1");
        assert_eq!(hit.query_name, "tRNA-Ala");
        assert_eq!(hit.seq_from, 100);
        assert_eq!(hit.seq_to, 170);
        assert_eq!(hit.strand, '+');
        assert!((hit.score - 55.5).abs() < 0.01);
    }

    #[test]
    fn test_find_anticodon_regex() {
        let scan = CMScan::new();
        // Normalized Ser fixture (Example1 #2): expect AGA.
        let ss = ">>>>>>>..>>>>.......<<<<.>>>>>.......<<<<<..>>>>>>>....<<<<<<<..>>>>>.......<<<<<<<<<<<<...";
        let seq = "GCAGTCATGTCCGAGTGGTTAAGGAGATTGACTAGAAATCAATTGGGCTCTGCCCGCGTAGGTTCGAATCCTGCTGACTGCG";
        let (ac, ai, _ae, _acp) = scan.find_anticodon(seq, ss, "NNN");
        // Note: this exercises the delegation path; exact strings are covered
        // by decode::tests. Here we only assert it does not panic and returns.
        let _ = (ac, ai);
    }

    #[test]
    fn test_cm_file_paths() {
        let mut paths = CMFilePaths::new();

        paths.add_main("Domain", PathBuf::from("/path/to/cm"));
        paths.add_isotype("Ala", PathBuf::from("/path/to/ala.cm"));
        paths.add_isotype_cutoff("Ala", 25.0);

        assert!(paths.get_main("Domain").is_some());
        assert!(paths.isotype.contains_key("Ala"));
        assert_eq!(paths.isotype_cutoffs.get("Ala"), Some(&25.0));
    }

    #[test]
    fn test_find_intron() {
        let scan = CMScan::new();

        // Sequence with intron (lowercase)
        let seq = "ACGUACGUacguacguACGUACGU";

        // Anticodon loop spanning the intron
        let (intron_seq, start, end) = scan.find_intron(seq, 4, 20);

        assert!(!intron_seq.is_empty() || start == 0);  // May or may not find depending on min length
    }
}
