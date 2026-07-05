//! tRNA representation module for tRNAscan-SE
//!
//! This module provides the core tRNA data structure and associated methods
//! for representing, manipulating, and outputting tRNA gene predictions.
//!
//! Ported from tRNAscanSE::tRNA.pm

use std::collections::HashMap;

// ============================================================================
// Constants
// ============================================================================

/// Type 2 tRNAs by clade (Leu, Ser have longer variable loops)
pub fn get_type2_trnas(clade: &str) -> Option<&'static [&'static str]> {
    match clade {
        "Eukaryota" => Some(&["Leu", "Ser"]),
        "Archaea" => Some(&["Leu", "Ser"]),
        "Bacteria" => Some(&["Leu", "Ser", "Tyr"]),
        _ => None,
    }
}

// ============================================================================
// Supporting Types
// ============================================================================

/// Anticodon position within tRNA sequence
#[derive(Debug, Clone, Default, PartialEq)]
pub struct AnticodonPos {
    /// Relative start position within tRNA
    pub rel_start: i32,
    /// Relative end position within tRNA
    pub rel_end: i32,
}

/// Intron information
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Intron {
    /// Relative start position within tRNA
    pub rel_start: i32,
    /// Relative end position within tRNA
    pub rel_end: i32,
    /// Absolute start position in sequence
    pub start: i64,
    /// Absolute end position in sequence
    pub end: i64,
    /// Intron type (e.g., "CI" for canonical intron)
    pub intron_type: String,
    /// Intron sequence
    pub seq: String,
}

/// Domain model hit information
#[derive(Debug, Clone, Default, PartialEq)]
pub struct DomainModel {
    /// Overall score
    pub score: f64,
    /// Mature score
    pub mat_score: f64,
    /// HMM score
    pub hmm_score: f64,
    /// Secondary structure score
    pub ss_score: f64,
}

/// Multi-model hit information
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ModelHit {
    /// Model type (e.g., "cove", "infernal")
    pub hit_type: String,
    /// Model name
    pub model: String,
    /// Score
    pub score: f64,
    /// Secondary structure
    pub ss: String,
}

/// Non-canonical marker types
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NonCanonical {
    /// Mismatch in stem
    Mismatch,
    /// Non-canonical base pair (wobble allowed)
    NonCanonical,
    /// Insertion
    Insertion,
    /// Deletion
    Deletion,
    /// Normal/canonical
    None,
}

impl Default for NonCanonical {
    fn default() -> Self {
        NonCanonical::None
    }
}

impl NonCanonical {
    pub fn from_char(c: char) -> Self {
        match c {
            'M' => NonCanonical::Mismatch,
            'N' => NonCanonical::NonCanonical,
            'I' => NonCanonical::Insertion,
            'D' => NonCanonical::Deletion,
            _ => NonCanonical::None,
        }
    }

    pub fn to_char(self) -> char {
        match self {
            NonCanonical::Mismatch => 'M',
            NonCanonical::NonCanonical => 'N',
            NonCanonical::Insertion => 'I',
            NonCanonical::Deletion => 'D',
            NonCanonical::None => ' ',
        }
    }
}

/// tRNA category classification
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub enum TRnaCategory {
    #[default]
    Unknown,
    Cytosolic,
    Mitochondrial,
    Organelle,
    NuMt,
    UndeterminedAc,
    MitoIsoConflict,
    MitoNoncanonicalAc,
    MitoAcMislocation,
    MitoMismatchAc,
}

impl TRnaCategory {
    pub fn from_str(s: &str) -> Self {
        match s {
            "cyto" => TRnaCategory::Cytosolic,
            "mito" => TRnaCategory::Mitochondrial,
            "org" => TRnaCategory::Organelle,
            "numt" => TRnaCategory::NuMt,
            "undetermined_ac" => TRnaCategory::UndeterminedAc,
            "mito_iso_conflict" => TRnaCategory::MitoIsoConflict,
            "mito_noncanonical_ac" => TRnaCategory::MitoNoncanonicalAc,
            "mito_ac_mislocation" => TRnaCategory::MitoAcMislocation,
            "mito_mismatch_ac" => TRnaCategory::MitoMismatchAc,
            _ => TRnaCategory::Unknown,
        }
    }

    pub fn as_str(&self) -> &'static str {
        match self {
            TRnaCategory::Unknown => "",
            TRnaCategory::Cytosolic => "cyto",
            TRnaCategory::Mitochondrial => "mito",
            TRnaCategory::Organelle => "org",
            TRnaCategory::NuMt => "numt",
            TRnaCategory::UndeterminedAc => "undetermined_ac",
            TRnaCategory::MitoIsoConflict => "mito_iso_conflict",
            TRnaCategory::MitoNoncanonicalAc => "mito_noncanonical_ac",
            TRnaCategory::MitoAcMislocation => "mito_ac_mislocation",
            TRnaCategory::MitoMismatchAc => "mito_mismatch_ac",
        }
    }
}

/// Strand orientation
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Strand {
    #[default]
    Unknown,
    Plus,
    Minus,
}

impl Strand {
    pub fn from_char(c: char) -> Self {
        match c {
            '+' => Strand::Plus,
            '-' => Strand::Minus,
            _ => Strand::Unknown,
        }
    }

    pub fn from_str(s: &str) -> Self {
        match s {
            "+" => Strand::Plus,
            "-" => Strand::Minus,
            _ => Strand::Unknown,
        }
    }

    pub fn to_char(self) -> char {
        match self {
            Strand::Plus => '+',
            Strand::Minus => '-',
            Strand::Unknown => '.',
        }
    }

    pub fn as_str(&self) -> &'static str {
        match self {
            Strand::Plus => "+",
            Strand::Minus => "-",
            Strand::Unknown => ".",
        }
    }
}

/// Truncation type
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub enum Truncation {
    #[default]
    None,
    FivePrime,
    ThreePrime,
    Both,
}

impl Truncation {
    pub fn from_str(s: &str) -> Self {
        match s {
            "5'" | "5" => Truncation::FivePrime,
            "3'" | "3" => Truncation::ThreePrime,
            "5'3'" | "both" => Truncation::Both,
            _ => Truncation::None,
        }
    }

    pub fn as_str(&self) -> &'static str {
        match self {
            Truncation::None => "",
            Truncation::FivePrime => "5'",
            Truncation::ThreePrime => "3'",
            Truncation::Both => "5'3'",
        }
    }

    pub fn is_truncated(&self) -> bool {
        !matches!(self, Truncation::None)
    }
}

// ============================================================================
// Main tRNA Struct
// ============================================================================

/// Core tRNA representation
#[derive(Debug, Clone, Default)]
pub struct TRna {
    // === Position Information ===
    /// Sequence name (chromosome/contig)
    pub seqname: String,
    /// Ordered sequence name (for sorting)
    pub ordered_seqname: usize,
    /// Start position (always smaller coordinate)
    pub start: i64,
    /// End position (always larger coordinate)
    pub end: i64,
    /// Strand orientation
    pub strand: Strand,

    // === Multi-exon support (for trans-spliced tRNAs) ===
    /// Second exon start
    pub start2: i64,
    /// Second exon end
    pub end2: i64,
    /// Second exon strand
    pub strand2: Strand,
    /// Third exon start
    pub start3: i64,
    /// Third exon end
    pub end3: i64,
    /// Third exon strand
    pub strand3: Strand,

    // === Identification ===
    /// Position/order in results
    pub position: usize,
    /// Internal ID number
    pub id: usize,
    /// tRNAscan-SE ID string (e.g., "chr1.trna1")
    pub trnascan_id: String,
    /// GtRNAdb ID
    pub gtrnadb_id: String,
    /// External database ID
    pub extdb_id: String,

    // === tRNA Identity ===
    /// Isotype (amino acid, e.g., "Ala", "Gly")
    pub isotype: String,
    /// Anticodon sequence (e.g., "TGC")
    pub anticodon: String,
    /// Anticodon positions
    pub ac_positions: Vec<AnticodonPos>,
    /// Clade (e.g., "Bacteria", "Archaea", "Eukaryota")
    pub clade: String,
    /// Model used for identification
    pub model: String,

    // === Introns ===
    /// List of introns
    pub introns: Vec<Intron>,

    // === Scores ===
    /// Primary score
    pub score: f64,
    /// Mature tRNA score
    pub mat_score: f64,
    /// HMM score
    pub hmm_score: f64,
    /// Secondary structure score
    pub ss_score: f64,
    /// Domain-specific model hits
    pub domain_models: HashMap<String, DomainModel>,
    /// Multiple model hits for comparison
    pub multi_models: Vec<ModelHit>,

    // === Sequences ===
    /// Full tRNA sequence (with intron)
    pub seq: String,
    /// Mature tRNA sequence (intron removed)
    pub mat_seq: String,
    /// Secondary structure (dot-bracket notation)
    pub ss: String,
    /// Mature secondary structure
    pub mat_ss: String,
    /// Sprinzl-aligned sequence
    pub sprinzl_align: String,
    /// Sprinzl-aligned secondary structure
    pub sprinzl_ss: String,

    // === Sprinzl Position Mapping ===
    /// Position to Sprinzl position mapping
    pub pos_sprinzl_map: Vec<String>,
    /// Sprinzl position to alignment position mapping
    pub sprinzl_pos_map: HashMap<String, usize>,

    // === Non-canonical positions ===
    /// Non-canonical markers by position
    pub non_canonical: Vec<NonCanonical>,

    // === Classification ===
    /// Is pseudogene
    pub is_pseudo: bool,
    /// Truncation status
    pub trunc: Truncation,
    /// Full truncation label from check_truncation (C: `trunc()`, e.g.
    /// "trunc_start:6" / "trunc_end:3" / "trunc_start:6,trunc_end:3"). Shown in
    /// the main `.out` Note only under `--detail`; struct output shows the
    /// "Possible truncation" text derived from `trunc` (the direction enum).
    pub trunc_label: String,
    /// Category classification
    pub category: TRnaCategory,

    // === Flanking Sequences ===
    /// Upstream sequence
    pub upstream: String,
    /// Downstream sequence
    pub downstream: String,

    // === Metadata ===
    /// Hit source (e.g., "Inf", "Eu", "Cove")
    pub hit_source: String,
    /// Source sequence ID
    pub src_seqid: usize,
    /// Source sequence length
    pub src_seqlen: usize,
    /// Additional notes
    pub note: String,

    // === Isotype-specific CM scan (`--detail`) ===
    /// Highest-scoring isotype CM model basename (e.g. "Ser", "Ile2").
    pub iso_model: String,
    /// Bit score of the highest-scoring isotype CM model.
    pub iso_score: f64,
}

impl TRna {
    /// Create a new empty tRNA
    pub fn new() -> Self {
        TRna::default()
    }

    /// Clear all fields to defaults
    pub fn clear(&mut self) {
        *self = TRna::default();
    }

    // ========================================================================
    // Basic Accessors
    // ========================================================================

    /// Get tRNA length
    pub fn len(&self) -> i64 {
        if self.start <= self.end {
            self.end - self.start + 1
        } else {
            self.start - self.end + 1
        }
    }

    /// Check if tRNA has zero length
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Check if this is a Type 2 tRNA (Leu/Ser with long variable loop)
    pub fn is_type2(&self) -> bool {
        if let Some(type2_list) = get_type2_trnas(&self.clade) {
            type2_list.iter().any(|&iso| iso == self.isotype)
        } else {
            false
        }
    }

    /// Check if truncated
    pub fn is_trunc(&self) -> bool {
        self.trunc.is_truncated()
    }

    /// Check if NuMT (nuclear mitochondrial)
    pub fn is_numt(&self) -> bool {
        matches!(self.category, TRnaCategory::NuMt)
    }

    /// Check if mitochondrial
    pub fn is_mito(&self) -> bool {
        matches!(self.category, TRnaCategory::Mitochondrial)
    }

    /// Check if undetermined anticodon
    pub fn is_undetermined(&self) -> bool {
        matches!(self.category, TRnaCategory::UndeterminedAc)
    }

    /// Check if cytosolic
    pub fn is_cytosolic(&self) -> bool {
        matches!(self.category, TRnaCategory::Cytosolic)
    }

    /// Check if organelle
    pub fn is_organelle(&self) -> bool {
        matches!(self.category, TRnaCategory::Organelle)
    }

    /// Check if mito isotype conflict
    pub fn is_mito_iso_conflict(&self) -> bool {
        matches!(self.category, TRnaCategory::MitoIsoConflict)
    }

    /// Check if mito non-canonical anticodon
    pub fn is_mito_noncanonical_ac(&self) -> bool {
        matches!(self.category, TRnaCategory::MitoNoncanonicalAc)
    }

    /// Check if mito anticodon mislocation
    pub fn is_mito_ac_mislocation(&self) -> bool {
        matches!(self.category, TRnaCategory::MitoAcMislocation)
    }

    /// Check if mito mismatched anticodon
    pub fn is_mito_mismatch_ac(&self) -> bool {
        matches!(self.category, TRnaCategory::MitoMismatchAc)
    }

    // ========================================================================
    // Exon Methods
    // ========================================================================

    /// Get exon start position
    pub fn exon_start(&self, exon: usize) -> i64 {
        match exon {
            1 => self.start,
            2 => self.start2,
            3 => self.start3,
            _ => 0,
        }
    }

    /// Set exon start position
    pub fn set_exon_start(&mut self, exon: usize, start: i64) {
        match exon {
            1 => self.start = start,
            2 => self.start2 = start,
            3 => self.start3 = start,
            _ => {}
        }
    }

    /// Get exon end position
    pub fn exon_end(&self, exon: usize) -> i64 {
        match exon {
            1 => self.end,
            2 => self.end2,
            3 => self.end3,
            _ => 0,
        }
    }

    /// Set exon end position
    pub fn set_exon_end(&mut self, exon: usize, end: i64) {
        match exon {
            1 => self.end = end,
            2 => self.end2 = end,
            3 => self.end3 = end,
            _ => {}
        }
    }

    /// Get exon strand
    pub fn exon_strand(&self, exon: usize) -> Strand {
        match exon {
            1 => self.strand,
            2 => self.strand2,
            3 => self.strand3,
            _ => Strand::Unknown,
        }
    }

    /// Set exon strand
    pub fn set_exon_strand(&mut self, exon: usize, strand: Strand) {
        match exon {
            1 => self.strand = strand,
            2 => self.strand2 = strand,
            3 => self.strand3 = strand,
            _ => {}
        }
    }

    // ========================================================================
    // Anticodon Methods
    // ========================================================================

    /// Add anticodon position
    pub fn add_ac_pos(&mut self, rel_start: i32, rel_end: i32) {
        self.ac_positions.push(AnticodonPos { rel_start, rel_end });
    }

    /// Get anticodon position count
    pub fn get_ac_pos_count(&self) -> usize {
        self.ac_positions.len()
    }

    /// Remove anticodon position at index
    pub fn remove_ac_pos(&mut self, index: usize) {
        if index < self.ac_positions.len() {
            self.ac_positions.remove(index);
        }
    }

    // ========================================================================
    // Intron Methods
    // ========================================================================

    /// Add intron with relative positions (calculates absolute positions)
    pub fn add_rel_intron(&mut self, rel_start: i32, rel_end: i32, intron_type: &str, seq: &str) {
        let (abs_start, abs_end) = match self.strand {
            Strand::Plus => {
                let start = self.start + (rel_start as i64) - 1;
                let end = self.start + (rel_end as i64) - 1;
                (start, end)
            }
            Strand::Minus => {
                let end = self.end - (rel_start as i64) + 1;
                let start = self.end - (rel_end as i64) + 1;
                (start, end)
            }
            Strand::Unknown => (0, 0),
        };
        self.add_intron(rel_start, rel_end, abs_start, abs_end, intron_type, seq);
    }

    /// Add intron with all positions specified
    pub fn add_intron(
        &mut self,
        rel_start: i32,
        rel_end: i32,
        abs_start: i64,
        abs_end: i64,
        intron_type: &str,
        seq: &str,
    ) {
        self.introns.push(Intron {
            rel_start,
            rel_end,
            start: abs_start,
            end: abs_end,
            intron_type: intron_type.to_string(),
            seq: seq.to_uppercase(),
        });
    }

    /// Set intron at index
    pub fn set_intron(&mut self, index: usize, rel_start: i32, rel_end: i32, intron_type: &str, seq: &str) {
        if index < self.introns.len() {
            self.introns[index].rel_start = rel_start;
            self.introns[index].rel_end = rel_end;
            self.introns[index].intron_type = intron_type.to_string();
            self.introns[index].seq = seq.to_string();

            // Update absolute positions
            match self.strand {
                Strand::Plus => {
                    self.introns[index].start = self.start + (rel_start as i64) - 1;
                    self.introns[index].end = self.start + (rel_end as i64) - 1;
                }
                Strand::Minus => {
                    self.introns[index].end = self.end - (rel_start as i64) + 1;
                    self.introns[index].start = self.end - (rel_end as i64) + 1;
                }
                Strand::Unknown => {}
            }
        }
    }

    /// Get intron at index
    pub fn get_intron(&self, index: usize) -> Option<&Intron> {
        self.introns.get(index)
    }

    /// Remove intron at index
    pub fn remove_intron(&mut self, index: usize) {
        if index < self.introns.len() {
            self.introns.remove(index);
        }
    }

    /// Get intron count
    pub fn get_intron_count(&self) -> usize {
        self.introns.len()
    }

    /// Sort introns by relative start position
    pub fn sort_introns(&mut self) {
        self.introns.sort_by_key(|i| i.rel_start);
    }

    /// Merge adjacent introns
    pub fn merge_introns(&mut self) {
        self.sort_introns();
        let mut to_remove = Vec::new();

        for ct in 1..self.introns.len() {
            let prev_end = self.introns[ct - 1].rel_end;
            if prev_end == self.introns[ct].rel_start - 1 {
                // Merge: extend current to include previous
                self.introns[ct].rel_start = self.introns[ct - 1].rel_start;
                match self.strand {
                    Strand::Plus => {
                        self.introns[ct].start = self.introns[ct - 1].start;
                    }
                    Strand::Minus => {
                        self.introns[ct].end = self.introns[ct - 1].end;
                    }
                    _ => {}
                }
                // Preserve CI type
                if self.introns[ct - 1].intron_type == "CI" {
                    self.introns[ct].intron_type = "CI".to_string();
                }
                to_remove.push(ct - 1);
            }
        }

        // Remove merged introns in reverse order
        for idx in to_remove.into_iter().rev() {
            self.introns.remove(idx);
        }
    }

    /// Get intron sequence at index
    pub fn get_intron_seq(&self, index: usize) -> String {
        if index < self.introns.len() {
            let intron = &self.introns[index];
            let start = (intron.rel_start - 1) as usize;
            let len = (intron.rel_end - intron.rel_start + 1) as usize;
            if start + len <= self.seq.len() {
                return self.seq[start..start + len].to_uppercase();
            }
        }
        String::new()
    }

    // ========================================================================
    // Domain Model Methods
    // ========================================================================

    /// Set domain model hit
    pub fn set_domain_model(&mut self, alg: &str, score: f64) {
        self.domain_models.insert(
            alg.to_string(),
            DomainModel {
                score,
                mat_score: 0.0,
                hmm_score: 0.0,
                ss_score: 0.0,
            },
        );
    }

    /// Get domain model
    pub fn get_domain_model(&self, alg: &str) -> Option<&DomainModel> {
        self.domain_models.get(alg)
    }

    /// Update domain model scores
    pub fn update_domain_model(
        &mut self,
        alg: &str,
        score: f64,
        mat_score: f64,
        hmm_score: f64,
        ss_score: f64,
    ) {
        if let Some(model) = self.domain_models.get_mut(alg) {
            model.score = score;
            model.mat_score = mat_score;
            model.hmm_score = hmm_score;
            model.ss_score = ss_score;
        }
    }

    /// Set default scores from domain models
    pub fn set_default_scores(&mut self) {
        if let Some(cove) = self.domain_models.get("cove") {
            self.score = cove.score;
            self.mat_score = cove.mat_score;
        } else if let Some(infernal) = self.domain_models.get("infernal") {
            self.score = infernal.score;
            self.mat_score = infernal.mat_score;
        }

        if let Some(infernal) = self.domain_models.get("infernal") {
            self.hmm_score = infernal.hmm_score;
            self.ss_score = infernal.ss_score;
        } else if let Some(cove) = self.domain_models.get("cove") {
            self.hmm_score = cove.hmm_score;
            self.ss_score = cove.ss_score;
        }
    }

    /// Get default score type
    pub fn get_default_score_type(&self) -> &str {
        if self.domain_models.contains_key("cove") {
            "cove"
        } else if self.domain_models.contains_key("infernal") {
            "infernal"
        } else {
            ""
        }
    }

    // ========================================================================
    // Multi-Model Methods
    // ========================================================================

    /// Add model hit
    pub fn add_model_hit(&mut self, hit_type: &str, model: &str, score: f64, ss: &str) {
        self.multi_models.push(ModelHit {
            hit_type: hit_type.to_string(),
            model: model.to_string(),
            score,
            ss: ss.to_string(),
        });
    }

    /// Sort multi-models by key
    pub fn sort_multi_models(&mut self, key: &str) {
        match key {
            "model" => {
                self.multi_models.sort_by(|a, b| {
                    a.hit_type.cmp(&b.hit_type).then(a.model.cmp(&b.model))
                });
            }
            "score" => {
                self.multi_models.sort_by(|a, b| {
                    b.score.partial_cmp(&a.score).unwrap_or(std::cmp::Ordering::Equal)
                });
            }
            _ => {}
        }
    }

    /// Get model hit by type and name (binary search, assumes sorted by model)
    pub fn get_model_hit(&self, hit_type: &str, model_name: &str) -> Option<&ModelHit> {
        self.multi_models.iter().find(|m| m.hit_type == hit_type && m.model == model_name)
    }

    /// Get highest scoring model
    pub fn get_highest_score_model(&mut self) -> Option<&ModelHit> {
        if self.multi_models.is_empty() {
            return None;
        }
        self.sort_multi_models("score");
        self.multi_models.first()
    }

    // ========================================================================
    // Mature tRNA Methods
    // ========================================================================

    /// Set mature tRNA sequence (remove introns)
    pub fn set_mature_trna(&mut self) {
        let set_ss = self.mat_ss.is_empty();

        if self.mat_seq.is_empty() {
            if self.introns.is_empty() {
                self.mat_seq = self.seq.clone();
                self.mat_ss = self.ss.clone();
            } else {
                let mut mat_seq = String::new();
                let mut mat_ss = String::new();

                for (i, intron) in self.introns.iter().enumerate() {
                    let start = if i == 0 {
                        0
                    } else {
                        self.introns[i - 1].rel_end as usize
                    };
                    let end = (intron.rel_start - 1) as usize;

                    if end > start && end <= self.seq.len() {
                        mat_seq.push_str(&self.seq[start..end]);
                        if set_ss && end <= self.ss.len() {
                            mat_ss.push_str(&self.ss[start..end]);
                        }
                    }
                }

                // Add sequence after last intron
                if let Some(last_intron) = self.introns.last() {
                    let start = last_intron.rel_end as usize;
                    if start < self.seq.len() {
                        mat_seq.push_str(&self.seq[start..]);
                        if set_ss && start < self.ss.len() {
                            mat_ss.push_str(&self.ss[start..]);
                        }
                    }
                }

                self.mat_seq = mat_seq;
                if set_ss {
                    self.mat_ss = mat_ss;
                }
            }
        }
    }

    // ========================================================================
    // Sprinzl Position Methods
    // ========================================================================

    /// Map Sprinzl positions from alignment
    pub fn map_sprinzl_pos(&mut self, ar_sprinzl_pos: &[&str]) {
        self.pos_sprinzl_map.clear();
        self.sprinzl_pos_map.clear();

        let bases: Vec<char> = self.sprinzl_align.chars().collect();
        let mut pos_count = 0;
        let mut ins_count = 0;

        for (i, &base) in bases.iter().enumerate() {
            if matches!(base, 'A' | 'C' | 'G' | 'T' | 'U' | '-') {
                if pos_count < ar_sprinzl_pos.len() {
                    let pos = ar_sprinzl_pos[pos_count].to_string();
                    self.pos_sprinzl_map.push(pos.clone());
                    self.sprinzl_pos_map.insert(pos, i);
                }
                pos_count += 1;
                ins_count = 0;
            } else {
                ins_count += 1;
                if pos_count > 0 && pos_count <= ar_sprinzl_pos.len() {
                    let pos = format!("{}:i{}", ar_sprinzl_pos[pos_count - 1], ins_count);
                    self.pos_sprinzl_map.push(pos.clone());
                    self.sprinzl_pos_map.insert(pos, i);
                }
            }
        }
    }

    /// Get Sprinzl position at relative position
    pub fn get_sprinzl_pos(&self, rel_pos: usize) -> Option<&str> {
        self.pos_sprinzl_map.get(rel_pos).map(|s| s.as_str())
    }

    /// Get relative alignment position from Sprinzl position
    pub fn get_rel_align_pos(&self, sprinzl_pos: &str) -> Option<usize> {
        self.sprinzl_pos_map.get(sprinzl_pos).copied()
    }

    /// Get base at Sprinzl position
    pub fn get_base_at_sprinzl(&self, sprinzl_pos: &str) -> Option<char> {
        self.get_rel_align_pos(sprinzl_pos)
            .and_then(|pos| self.sprinzl_align.chars().nth(pos))
    }

    /// Get base at relative alignment position
    pub fn get_base_at_rel_align_pos(&self, rel_pos: usize) -> Option<char> {
        self.sprinzl_align.chars().nth(rel_pos)
    }

    /// Get sequence at relative alignment positions
    pub fn get_seq_at_rel_align_pos(&self, start: usize, end: usize) -> String {
        if end >= start && end < self.sprinzl_align.len() {
            self.sprinzl_align[start..=end].to_string()
        } else {
            String::new()
        }
    }

    /// Get relative position from absolute coordinate
    pub fn get_rel_pos(&self, coord: i64) -> Option<i32> {
        if coord >= self.start && coord <= self.end {
            let rel_pos = match self.strand {
                Strand::Plus => coord - self.start + 1,
                Strand::Minus => self.end - coord + 1,
                Strand::Unknown => return None,
            };
            Some(rel_pos as i32)
        } else {
            None
        }
    }

    /// Get relative mature position (accounting for introns)
    pub fn get_rel_mature_pos(&self, coord: i64) -> Option<i32> {
        let rel_pos = self.get_rel_pos(coord)?;
        let mut rel_mature_pos = rel_pos;

        for intron in self.introns.iter().rev() {
            if rel_pos > intron.rel_end {
                rel_mature_pos -= intron.rel_end - intron.rel_start + 1;
            }
        }

        Some(rel_mature_pos)
    }

    /// Check if coordinate is within an intron
    pub fn is_intronic(&self, coord: i64) -> bool {
        if let Some(rel_pos) = self.get_rel_pos(coord) {
            self.introns.iter().any(|intron| {
                rel_pos >= intron.rel_start && rel_pos <= intron.rel_end
            })
        } else {
            false
        }
    }

    /// Convert coordinate to Sprinzl position
    pub fn convert_sprinzl_pos(&self, coord: i64) -> Option<String> {
        if self.is_intronic(coord) {
            return None;
        }

        let rel_mature_pos = self.get_rel_mature_pos(coord)?;
        let bases: Vec<char> = self.sprinzl_align.chars().collect();
        let mut ct = 0;

        for (i, &base) in bases.iter().enumerate() {
            if base != '-' {
                ct += 1;
                if ct == rel_mature_pos {
                    return self.get_sprinzl_pos(i).map(|s| s.to_string());
                }
            }
        }

        None
    }

    /// Get sequence fragment by relative positions
    pub fn get_seq_fragment(&self, rel_start: i32, rel_end: i32) -> String {
        if rel_start > 0 && rel_end <= self.len() as i32 {
            let start = (rel_start - 1) as usize;
            let len = (rel_end - rel_start + 1) as usize;
            if start + len <= self.seq.len() {
                return self.seq[start..start + len].to_string();
            }
        }
        String::new()
    }

    // ========================================================================
    // Non-canonical Methods
    // ========================================================================

    /// Add non-canonical marker at position
    pub fn add_non_canonical(&mut self, rel_pos: usize, nc: NonCanonical) {
        if rel_pos >= self.non_canonical.len() {
            self.non_canonical.resize(rel_pos + 1, NonCanonical::None);
        }
        self.non_canonical[rel_pos] = nc;
    }

    /// Get non-canonical marker at Sprinzl position
    pub fn get_non_canonical(&self, sprinzl_pos: &str) -> NonCanonical {
        self.get_rel_align_pos(sprinzl_pos)
            .and_then(|pos| self.non_canonical.get(pos).copied())
            .unwrap_or(NonCanonical::None)
    }

    /// Get mismatch count (pairs counted as 1)
    pub fn get_mismatch_count(&self) -> usize {
        self.non_canonical
            .iter()
            .filter(|&&nc| nc == NonCanonical::Mismatch)
            .count() / 2
    }

    /// Get non-canonical count (pairs counted as 1)
    pub fn get_noncanonical_count(&self) -> usize {
        self.non_canonical
            .iter()
            .filter(|&&nc| nc == NonCanonical::NonCanonical)
            .count() / 2
    }

    /// Check if position is a mismatch
    pub fn is_mismatch(&self, rel_pos: usize) -> bool {
        self.non_canonical.get(rel_pos).copied() == Some(NonCanonical::Mismatch)
    }

    /// Check if position is an insertion
    pub fn is_insertion(&self, rel_pos: usize) -> bool {
        self.non_canonical.get(rel_pos).copied() == Some(NonCanonical::Insertion)
    }

    /// Check if position is a deletion
    pub fn is_deletion(&self, rel_pos: usize) -> bool {
        self.non_canonical.get(rel_pos).copied() == Some(NonCanonical::Deletion)
    }

    /// Check if position is non-canonical
    pub fn is_non_canonical(&self, rel_pos: usize) -> bool {
        self.non_canonical.get(rel_pos).copied() == Some(NonCanonical::NonCanonical)
    }

    /// Get annotated mismatches string
    pub fn get_annotated_mismatches(&self, sprinzl_pairs: &HashMap<String, String>) -> String {
        let mut result = Vec::new();

        for (i, &nc) in self.non_canonical.iter().enumerate() {
            if nc == NonCanonical::Mismatch {
                if let Some(pos1) = self.get_sprinzl_pos(i) {
                    if let Some(pos2) = sprinzl_pairs.get(pos1) {
                        if let (Some(b1), Some(b2)) = (
                            self.get_base_at_rel_align_pos(i),
                            self.get_base_at_sprinzl(pos2),
                        ) {
                            result.push(format!("{}{}:{}{}", b1, pos1, b2, pos2));
                        }
                    }
                }
            }
        }

        result.join(" ")
    }

    /// Get annotated non-canonical positions string
    pub fn get_annotated_noncanonical(&self) -> String {
        let mut result = Vec::new();

        for (i, &nc) in self.non_canonical.iter().enumerate() {
            if nc == NonCanonical::NonCanonical {
                if let (Some(base), Some(pos)) = (
                    self.get_base_at_rel_align_pos(i),
                    self.get_sprinzl_pos(i),
                ) {
                    result.push(format!("{}{}", base, pos));
                }
            }
        }

        result.join(" ")
    }

    /// Check stem mismatches using Sprinzl pair mapping
    pub fn check_stem_mismatches(&mut self, sprinzl_pairs: &HashMap<String, String>) {
        let bases: Vec<char> = self.sprinzl_align.chars().collect();

        // Ensure non_canonical is properly sized
        if self.non_canonical.len() < bases.len() {
            self.non_canonical.resize(bases.len(), NonCanonical::None);
        }

        for i in 0..bases.len() {
            if let Some(pos1) = self.pos_sprinzl_map.get(i) {
                // Skip insertions and position 13 and "e" positions
                if pos1.contains(":i") || pos1 == "13" || pos1.contains('e') {
                    continue;
                }

                if let Some(pos2_str) = sprinzl_pairs.get(pos1) {
                    if let Some(&pos2) = self.sprinzl_pos_map.get(pos2_str) {
                        let b1 = bases[i];
                        let b2 = bases.get(pos2).copied().unwrap_or('-');

                        if !Self::is_valid_base_pair(b1, b2) {
                            self.non_canonical[i] = NonCanonical::Mismatch;
                            if pos2 < self.non_canonical.len() {
                                self.non_canonical[pos2] = NonCanonical::Mismatch;
                            }
                        }
                    }
                }
            }
        }
    }

    /// Check if bases form a valid (WC or wobble) pair
    fn is_valid_base_pair(b1: char, b2: char) -> bool {
        matches!(
            (b1, b2),
            ('A', 'T') | ('T', 'A') | ('A', 'U') | ('U', 'A') |
            ('G', 'C') | ('C', 'G') |
            ('G', 'T') | ('T', 'G') | ('G', 'U') | ('U', 'G') |
            ('-', _) | (_, '-')
        )
    }

    /// Check if bases form a Watson-Crick pair
    pub fn is_wc_base_pair(b1: char, b2: char) -> bool {
        matches!(
            (b1, b2),
            ('A', 'T') | ('T', 'A') | ('A', 'U') | ('U', 'A') |
            ('G', 'C') | ('C', 'G') | ('-', _) | (_, '-')
        )
    }

    /// Check if bases form a G-U wobble pair
    pub fn is_gu_base_pair(b1: char, b2: char) -> bool {
        matches!(
            (b1, b2),
            ('G', 'T') | ('T', 'G') | ('G', 'U') | ('U', 'G')
        )
    }

    // ========================================================================
    // Hit Source Methods
    // ========================================================================

    /// Set from Covels hit
    pub fn set_covels_hit(&mut self, hit_seqname: &str, score: f64, trna_start: i64, trna_end: i64) {
        self.seqname = hit_seqname.to_string();
        self.set_domain_model("cove", score);

        if trna_start < trna_end {
            self.strand = Strand::Plus;
            self.start = trna_start;
            self.end = trna_end;
        } else {
            self.strand = Strand::Minus;
            self.start = trna_end;
            self.end = trna_start;
        }
    }

    /// Set from cmsearch hit
    pub fn set_cmsearch_hit(
        &mut self,
        hit_seqname: &str,
        score: f64,
        strand: &str,
        trna_start: i64,
        trna_end: i64,
        seq: &str,
        ss: &str,
    ) {
        self.seqname = hit_seqname.to_string();
        self.set_domain_model("infernal", score);
        self.strand = Strand::from_str(strand);

        if self.strand == Strand::Plus {
            self.start = trna_start;
            self.end = trna_end;
        } else {
            self.start = trna_end;
            self.end = trna_start;
        }

        self.seq = seq.to_string();
        self.ss = ss.to_string();
    }

    // ========================================================================
    // Output Methods
    // ========================================================================

    /// Generate FASTA format output
    pub fn to_fasta(&self, use_mature: bool) -> String {
        let seq = if use_mature && !self.mat_seq.is_empty() {
            &self.mat_seq
        } else {
            &self.seq
        };

        let header = format!(
            ">{} {}-{} ({}) {} ({}) {}",
            self.seqname,
            self.start,
            self.end,
            self.strand.to_char(),
            self.isotype,
            self.anticodon,
            self.score
        );

        format!("{}\n{}", header, seq)
    }

    /// Generate BED format output (0-based, half-open)
    pub fn to_bed(&self) -> String {
        let score = (self.score * 10.0).round() as i32;
        let score = score.max(0).min(1000);

        format!(
            "{}\t{}\t{}\t{}-{}\t{}\t{}",
            self.seqname,
            self.start - 1,  // BED is 0-based
            self.end,
            self.isotype,
            self.anticodon,
            score,
            self.strand.to_char()
        )
    }

    /// Generate standard output line
    pub fn format_output_line(&self) -> String {
        let intron_info = if self.introns.is_empty() {
            "0\t0".to_string()
        } else {
            format!("{}\t{}", self.introns[0].start, self.introns[0].end)
        };

        let pseudo_str = if self.is_pseudo { "pseudo" } else { "" };

        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{}\t{}",
            self.seqname,
            self.id,
            self.start,
            self.end,
            self.isotype,
            self.anticodon,
            intron_info,
            self.score,
            self.note,
            pseudo_str,
            self.trunc.as_str()
        )
    }

    /// Generate secondary structure output
    pub fn format_ss_output(&self) -> String {
        format!(
            "Seq: {}\nStr: {}",
            if self.mat_seq.is_empty() { &self.seq } else { &self.mat_seq },
            if self.mat_ss.is_empty() { &self.ss } else { &self.mat_ss }
        )
    }

    // ========================================================================
    // Serialization
    // ========================================================================

    /// Convert to JSON string
    pub fn to_json(&self) -> String {
        let mut json = String::from("{\n");

        // Basic fields
        json.push_str(&format!("  \"seqname\": \"{}\",\n", self.seqname));
        json.push_str(&format!("  \"start\": {},\n", self.start));
        json.push_str(&format!("  \"end\": {},\n", self.end));
        json.push_str(&format!("  \"strand\": \"{}\",\n", self.strand.to_char()));
        json.push_str(&format!("  \"isotype\": \"{}\",\n", self.isotype));
        json.push_str(&format!("  \"anticodon\": \"{}\",\n", self.anticodon));
        json.push_str(&format!("  \"score\": {},\n", self.score));
        json.push_str(&format!("  \"is_pseudo\": {},\n", self.is_pseudo));
        json.push_str(&format!("  \"trunc\": \"{}\",\n", self.trunc.as_str()));
        json.push_str(&format!("  \"category\": \"{}\",\n", self.category.as_str()));
        json.push_str(&format!("  \"seq\": \"{}\",\n", self.seq));
        json.push_str(&format!("  \"ss\": \"{}\",\n", self.ss));
        json.push_str(&format!("  \"mat_seq\": \"{}\",\n", self.mat_seq));
        json.push_str(&format!("  \"note\": \"{}\"\n", self.note.replace('"', "\\\"")));

        json.push('}');
        json
    }

    /// Parse from JSON string (basic implementation)
    pub fn from_json(json: &str) -> Result<TRna, String> {
        let mut trna = TRna::new();

        // Simple JSON parsing (for production, use serde_json)
        for line in json.lines() {
            let line = line.trim();
            if let Some(rest) = line.strip_prefix("\"seqname\": \"") {
                trna.seqname = rest.trim_end_matches("\",").trim_end_matches('"').to_string();
            } else if let Some(rest) = line.strip_prefix("\"start\": ") {
                trna.start = rest.trim_end_matches(',').parse().map_err(|e| format!("{}", e))?;
            } else if let Some(rest) = line.strip_prefix("\"end\": ") {
                trna.end = rest.trim_end_matches(',').parse().map_err(|e| format!("{}", e))?;
            } else if let Some(rest) = line.strip_prefix("\"strand\": \"") {
                let c = rest.chars().next().unwrap_or('.');
                trna.strand = Strand::from_char(c);
            } else if let Some(rest) = line.strip_prefix("\"isotype\": \"") {
                trna.isotype = rest.trim_end_matches("\",").trim_end_matches('"').to_string();
            } else if let Some(rest) = line.strip_prefix("\"anticodon\": \"") {
                trna.anticodon = rest.trim_end_matches("\",").trim_end_matches('"').to_string();
            } else if let Some(rest) = line.strip_prefix("\"score\": ") {
                trna.score = rest.trim_end_matches(',').parse().map_err(|e| format!("{}", e))?;
            } else if let Some(rest) = line.strip_prefix("\"is_pseudo\": ") {
                trna.is_pseudo = rest.trim_end_matches(',') == "true";
            } else if let Some(rest) = line.strip_prefix("\"trunc\": \"") {
                trna.trunc = Truncation::from_str(rest.trim_end_matches("\",").trim_end_matches('"'));
            } else if let Some(rest) = line.strip_prefix("\"category\": \"") {
                trna.category = TRnaCategory::from_str(rest.trim_end_matches("\",").trim_end_matches('"'));
            } else if let Some(rest) = line.strip_prefix("\"seq\": \"") {
                trna.seq = rest.trim_end_matches("\",").trim_end_matches('"').to_string();
            } else if let Some(rest) = line.strip_prefix("\"ss\": \"") {
                trna.ss = rest.trim_end_matches("\",").trim_end_matches('"').to_string();
            } else if let Some(rest) = line.strip_prefix("\"mat_seq\": \"") {
                trna.mat_seq = rest.trim_end_matches("\",").trim_end_matches('"').to_string();
            } else if let Some(rest) = line.strip_prefix("\"note\": \"") {
                trna.note = rest.trim_end_matches('"').replace("\\\"", "\"");
            }
        }

        Ok(trna)
    }

    /// Convert to TSV string
    pub fn to_tsv(&self) -> String {
        let intron_starts: Vec<String> = self.introns.iter().map(|i| i.start.to_string()).collect();
        let intron_ends: Vec<String> = self.introns.iter().map(|i| i.end.to_string()).collect();

        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.seqname,
            self.start,
            self.end,
            self.strand.to_char(),
            self.isotype,
            self.anticodon,
            if intron_starts.is_empty() { "0".to_string() } else { intron_starts.join(",") },
            if intron_ends.is_empty() { "0".to_string() } else { intron_ends.join(",") },
            self.score,
            self.is_pseudo,
            self.trunc.as_str(),
            self.category.as_str(),
            self.seq,
            self.note
        )
    }

    /// Parse from TSV line
    pub fn from_tsv(line: &str) -> Result<TRna, String> {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 10 {
            return Err(format!("Expected at least 10 fields, got {}", fields.len()));
        }

        let mut trna = TRna::new();
        trna.seqname = fields[0].to_string();
        trna.start = fields[1].parse().map_err(|e| format!("Invalid start: {}", e))?;
        trna.end = fields[2].parse().map_err(|e| format!("Invalid end: {}", e))?;
        trna.strand = Strand::from_str(fields[3]);
        trna.isotype = fields[4].to_string();
        trna.anticodon = fields[5].to_string();

        // Parse introns
        let intron_starts: Vec<i64> = fields[6]
            .split(',')
            .filter(|s| !s.is_empty() && *s != "0")
            .filter_map(|s| s.parse().ok())
            .collect();
        let intron_ends: Vec<i64> = fields[7]
            .split(',')
            .filter(|s| !s.is_empty() && *s != "0")
            .filter_map(|s| s.parse().ok())
            .collect();

        for (start, end) in intron_starts.into_iter().zip(intron_ends.into_iter()) {
            // Calculate relative positions from absolute
            let (rel_start, rel_end) = match trna.strand {
                Strand::Plus => (
                    (start - trna.start + 1) as i32,
                    (end - trna.start + 1) as i32,
                ),
                Strand::Minus => (
                    (trna.end - end + 1) as i32,
                    (trna.end - start + 1) as i32,
                ),
                Strand::Unknown => (0, 0),
            };
            trna.introns.push(Intron {
                rel_start,
                rel_end,
                start,
                end,
                intron_type: String::new(),
                seq: String::new(),
            });
        }

        trna.score = fields[8].parse().map_err(|e| format!("Invalid score: {}", e))?;
        trna.is_pseudo = fields[9] == "true" || fields[9] == "1";

        if fields.len() > 10 {
            trna.trunc = Truncation::from_str(fields[10]);
        }
        if fields.len() > 11 {
            trna.category = TRnaCategory::from_str(fields[11]);
        }
        if fields.len() > 12 {
            trna.seq = fields[12].to_string();
        }
        if fields.len() > 13 {
            trna.note = fields[13].to_string();
        }

        Ok(trna)
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_trna_new() {
        let trna = TRna::new();
        assert!(trna.seqname.is_empty());
        assert_eq!(trna.start, 0);
        assert_eq!(trna.strand, Strand::Unknown);
    }

    #[test]
    fn test_trna_len() {
        let mut trna = TRna::new();
        trna.start = 100;
        trna.end = 175;
        assert_eq!(trna.len(), 76);
    }

    #[test]
    fn test_strand_conversion() {
        assert_eq!(Strand::from_char('+'), Strand::Plus);
        assert_eq!(Strand::from_char('-'), Strand::Minus);
        assert_eq!(Strand::Plus.to_char(), '+');
        assert_eq!(Strand::Minus.to_char(), '-');
    }

    #[test]
    fn test_truncation() {
        assert_eq!(Truncation::from_str("5'"), Truncation::FivePrime);
        assert_eq!(Truncation::from_str("3'"), Truncation::ThreePrime);
        assert!(Truncation::FivePrime.is_truncated());
        assert!(!Truncation::None.is_truncated());
    }

    #[test]
    fn test_type2_trna() {
        let mut trna = TRna::new();
        trna.clade = "Bacteria".to_string();
        trna.isotype = "Leu".to_string();
        assert!(trna.is_type2());

        trna.isotype = "Ala".to_string();
        assert!(!trna.is_type2());
    }

    #[test]
    fn test_intron_management() {
        let mut trna = TRna::new();
        trna.start = 100;
        trna.end = 175;
        trna.strand = Strand::Plus;

        trna.add_rel_intron(38, 55, "CI", "ATCGATCGATCGATCG");
        assert_eq!(trna.get_intron_count(), 1);

        let intron = trna.get_intron(0).unwrap();
        assert_eq!(intron.rel_start, 38);
        assert_eq!(intron.rel_end, 55);
    }

    #[test]
    fn test_domain_model() {
        let mut trna = TRna::new();
        trna.set_domain_model("infernal", 65.5);

        assert!(trna.get_domain_model("infernal").is_some());
        assert_eq!(trna.get_domain_model("infernal").unwrap().score, 65.5);
        assert!(trna.get_domain_model("cove").is_none());
    }

    #[test]
    fn test_mature_trna() {
        let mut trna = TRna::new();
        trna.seq = "ATCGATCGXXXXXXXXXXATCGATCG".to_string();
        trna.ss = "..((....))........((....))" .to_string();
        trna.start = 1;
        trna.end = 26;
        trna.strand = Strand::Plus;

        // Add intron from position 9-18 (1-indexed)
        trna.add_intron(9, 18, 9, 18, "CI", "XXXXXXXXXX");

        trna.set_mature_trna();
        assert_eq!(trna.mat_seq, "ATCGATCGATCGATCG");
    }

    #[test]
    fn test_base_pair_checks() {
        assert!(TRna::is_wc_base_pair('A', 'T'));
        assert!(TRna::is_wc_base_pair('G', 'C'));
        assert!(!TRna::is_wc_base_pair('G', 'T'));

        assert!(TRna::is_gu_base_pair('G', 'T'));
        assert!(TRna::is_gu_base_pair('G', 'U'));
        assert!(!TRna::is_gu_base_pair('A', 'T'));
    }

    #[test]
    fn test_to_fasta() {
        let mut trna = TRna::new();
        trna.seqname = "chr1".to_string();
        trna.start = 1000;
        trna.end = 1075;
        trna.strand = Strand::Plus;
        trna.isotype = "Ala".to_string();
        trna.anticodon = "TGC".to_string();
        trna.score = 65.5;
        trna.seq = "ATCGATCG".to_string();

        let fasta = trna.to_fasta(false);
        assert!(fasta.starts_with(">chr1 1000-1075 (+) Ala (TGC)"));
        assert!(fasta.contains("ATCGATCG"));
    }

    #[test]
    fn test_to_bed() {
        let mut trna = TRna::new();
        trna.seqname = "chr1".to_string();
        trna.start = 1000;
        trna.end = 1075;
        trna.strand = Strand::Plus;
        trna.isotype = "Ala".to_string();
        trna.anticodon = "TGC".to_string();
        trna.score = 65.5;

        let bed = trna.to_bed();
        assert!(bed.starts_with("chr1\t999\t1075\tAla-TGC"));
    }

    #[test]
    fn test_json_roundtrip() {
        let mut trna = TRna::new();
        trna.seqname = "chr1".to_string();
        trna.start = 1000;
        trna.end = 1075;
        trna.strand = Strand::Plus;
        trna.isotype = "Ala".to_string();
        trna.anticodon = "TGC".to_string();
        trna.score = 65.5;

        let json = trna.to_json();
        let parsed = TRna::from_json(&json).unwrap();

        assert_eq!(parsed.seqname, trna.seqname);
        assert_eq!(parsed.start, trna.start);
        assert_eq!(parsed.strand, trna.strand);
        assert_eq!(parsed.isotype, trna.isotype);
    }

    #[test]
    fn test_tsv_roundtrip() {
        let mut trna = TRna::new();
        trna.seqname = "chr1".to_string();
        trna.start = 1000;
        trna.end = 1075;
        trna.strand = Strand::Plus;
        trna.isotype = "Ala".to_string();
        trna.anticodon = "TGC".to_string();
        trna.score = 65.5;
        trna.is_pseudo = false;
        trna.seq = "ATCGATCG".to_string();

        let tsv = trna.to_tsv();
        let parsed = TRna::from_tsv(&tsv).unwrap();

        assert_eq!(parsed.seqname, trna.seqname);
        assert_eq!(parsed.start, trna.start);
        assert_eq!(parsed.strand, trna.strand);
        assert_eq!(parsed.isotype, trna.isotype);
    }

    #[test]
    fn test_category() {
        let mut trna = TRna::new();
        trna.category = TRnaCategory::Mitochondrial;
        assert!(trna.is_mito());
        assert!(!trna.is_cytosolic());

        trna.category = TRnaCategory::NuMt;
        assert!(trna.is_numt());
    }

    #[test]
    fn test_non_canonical() {
        let mut trna = TRna::new();
        trna.add_non_canonical(5, NonCanonical::Mismatch);
        trna.add_non_canonical(10, NonCanonical::Mismatch);

        assert!(trna.is_mismatch(5));
        assert!(trna.is_mismatch(10));
        assert!(!trna.is_mismatch(7));
        assert_eq!(trna.get_mismatch_count(), 1); // Pairs counted as 1
    }
}
