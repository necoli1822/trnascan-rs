// TrnaScanner: High-level API for tRNA scanning
//
// This module provides the main public interface for tRNAscan-SE.
// It orchestrates the entire pipeline:
// 1. First-pass detection (EufindtRNA for eukaryotes, tRNAscan for prokaryotes)
// 2. Second-pass Infernal cmsearch verification
// 3. Secondary structure prediction
// 4. Isotype classification
// 5. Output generation

use rayon::prelude::*;
use std::cell::RefCell;
use std::io::{Result as IoResult, Write};
use std::path::{Path, PathBuf};

use crate::cm_scan::{CMScan, CMScanOptions, CMSearchHit};
use crate::cm_scan::decode::{
    decode_trna_properties, find_anticodon, find_intron, format_cmsearch_output, get_trna_type,
    AliDisplay, UNDEF_ANTICODON, UNDEF_ISOTYPE,
};
use crate::trna::{AnticodonPos, Intron, TRna, Strand as TStrand, Truncation};
use infernal::{FaithfulConfig, FaithfulSearcher};
use once_cell::sync::Lazy;
use regex::Regex;

/// Bulge-helix-bulge (BHB) intron structure matcher, verbatim from
/// `CM.pm::check_intron_validity` (:1895). Operates on the RAW cmsearch CS line
/// (`<`,`>`,`-`,`_`,`.`) of a BHB-CM hit: group 1 = 5' exon flank, group 2 = the
/// intron (its own helix), group 3 = 3' exon flank.
static BHB_SS_RE: Lazy<Regex> = Lazy::new(|| {
    Regex::new(r"^([<\-.]{11,})(-<[<.]+[_.]{4,}[>.]{9,}-[.]*-)([-.>]+)$").unwrap()
});
use crate::eufind::{run_eufind_scan, EufindOptions, EufindHit};
use crate::isotype::anticodon_to_isotype;
use crate::squid::SqInfo;
use crate::trnascan::{trnascan, ConsensusMatrix, SearchParams, TrnascanHit};
use crate::trnascan::seq_utils::reverse_complement;

/// Scan mode for different organism types
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScanMode {
    /// Eukaryotic mode - uses EuFindtRNA first-pass
    Eukaryotic,
    /// Bacterial mode - uses tRNAscan (Fichant-Burks) first-pass
    Bacterial,
    /// Archaeal mode - uses tRNAscan (Fichant-Burks) first-pass
    Archaeal,
    /// Organellar mode - uses EuFindtRNA first-pass
    Organellar,
    /// General mode - uses both first-pass methods
    General,
    /// Mitochondrial mode - specialized scanning
    Mitochondrial,
}

impl ScanMode {
    /// Convert from char (CLI mode character)
    pub fn from_char(c: char) -> Self {
        match c.to_ascii_uppercase() {
            'B' => ScanMode::Bacterial,
            'A' => ScanMode::Archaeal,
            'O' => ScanMode::Organellar,
            'G' => ScanMode::General,
            'M' => ScanMode::Mitochondrial,
            _ => ScanMode::Eukaryotic,
        }
    }

    /// Get the corresponding CM model filename
    pub fn cm_model_name(&self) -> &'static str {
        match self {
            ScanMode::Bacterial => "TRNAinf-bact.cm",
            ScanMode::Archaeal => "TRNAinf-arch.cm",
            ScanMode::Eukaryotic => "TRNAinf-euk.cm",
            ScanMode::Organellar => "TRNAinf-euk.cm",
            ScanMode::General => "TRNAinf.cm",
            ScanMode::Mitochondrial => "TRNAinf-mito-mammal.cm",
        }
    }

    /// Check if this mode uses tRNAscan (Fichant-Burks) for first-pass
    pub fn uses_trnascan(&self) -> bool {
        matches!(self, ScanMode::Bacterial | ScanMode::Archaeal | ScanMode::General)
    }

    /// Check if this mode uses EuFindtRNA for first-pass
    pub fn uses_eufind(&self) -> bool {
        matches!(self, ScanMode::Eukaryotic | ScanMode::Organellar | ScanMode::General | ScanMode::Mitochondrial)
    }
}

/// Result of a tRNA scan
#[derive(Debug, Clone)]
pub struct TrnaResult {
    pub seq_name: String,
    pub trna_num: usize,
    pub begin: i64,
    pub end: i64,
    pub isotype: String,
    pub anticodon: String,
    pub intron_begin: usize,
    pub intron_end: usize,
    pub inf_score: f64,
    pub hmm_score: f64,
    pub ss_score: f64,
    pub origin: String,
    pub cm_isotype: String,
    pub cm_score: f64,
    pub note: String,
    pub sequence: Vec<u8>,
    pub structure: String,
    pub strand: char,
}

impl TrnaResult {
    /// Create from a CM search hit
    pub fn from_cm_hit(hit: &CMSearchHit, seq_name: &str, offset: usize) -> Self {
        let begin = hit.seq_from + offset as i64;
        let end = hit.seq_to + offset as i64;

        TrnaResult {
            seq_name: seq_name.to_string(),
            trna_num: 0, // Will be assigned later
            begin,
            end,
            isotype: hit.isotype.clone(),
            anticodon: hit.anticodon.clone(),
            intron_begin: if hit.intron_start > 0 { hit.intron_start as usize } else { 0 },
            intron_end: if hit.intron_end > 0 { hit.intron_end as usize } else { 0 },
            inf_score: hit.score,
            hmm_score: hit.hmm_score,
            ss_score: hit.ss_score,
            origin: "Inf".to_string(),
            cm_isotype: hit.isotype.clone(),
            cm_score: hit.score,
            note: hit.notes.clone(),
            sequence: hit.target_seq.replace("-", "").replace(".", "").into_bytes(),
            structure: hit.ss.clone(),
            strand: hit.strand,
        }
    }

    /// Format as tabular output line
    pub fn format_output_line(&self) -> String {
        format!(
            "{:<12}\t{}\t{:<5}\t{:<6}\t{:<4}\t{:<5}\t{}\t{}\t{:.1}\t{:.2}\t{:.2}\t{:<6}\t{:<7}\t{:.1}\t{}",
            self.seq_name,
            self.trna_num,
            self.begin,
            self.end,
            self.isotype,
            self.anticodon,
            if self.intron_begin > 0 {
                format!("{}", self.intron_begin)
            } else {
                "0".to_string()
            },
            if self.intron_end > 0 {
                format!("{}", self.intron_end)
            } else {
                "0".to_string()
            },
            self.inf_score,
            self.hmm_score,
            self.ss_score,
            self.origin,
            self.cm_isotype,
            self.cm_score,
            self.note
        )
    }

    /// Format as secondary structure output
    pub fn format_ss_output(&self) -> String {
        let seq_str = String::from_utf8_lossy(&self.sequence);
        let length = self.sequence.len();

        let mut result = String::new();
        result.push_str(&format!("{}.trna{} ({}-{})\tLength: {} bp\n",
            self.seq_name, self.trna_num, self.begin, self.end, length));
        result.push_str(&format!("Type: {}\tAnticodon: {} at {}-{} ({}-{})\tScore: {:.1}\n",
            self.isotype, self.anticodon,
            34, 36, // Position in tRNA
            self.begin + 33, self.begin + 35, // Genomic position
            self.inf_score));

        if self.intron_begin > 0 {
            result.push_str(&format!("Possible intron: {}-{} ({}-{})\n",
                38, 38 + (self.intron_end - self.intron_begin),
                self.intron_begin, self.intron_end));
        }

        result.push_str(&format!("HMM Sc={:.2}\tSec struct Sc={:.2}\n",
            self.hmm_score, self.ss_score));

        if !self.note.is_empty() {
            result.push_str(&format!("{}\n", self.note));
        }

        // Format sequence and structure with position markers
        let marker_line = (0..=length/10)
            .map(|_i| format!("         *    |"))
            .collect::<String>();
        result.push_str(&format!("{}\n", &marker_line[0..length.min(marker_line.len())]));
        result.push_str(&format!("Seq: {}\n", seq_str));
        result.push_str(&format!("Str: {}\n", self.structure));
        result.push('\n');

        result
    }

    /// Format as BED line
    pub fn format_bed_line(&self, chrom_name: &str) -> String {
        format!("{}\t{}\t{}\t{}.trna{}\t{:.1}\t{}",
            chrom_name,
            self.begin - 1,  // BED is 0-based
            self.end,
            self.seq_name,
            self.trna_num,
            self.inf_score,
            self.strand
        )
    }
}

/// First-pass hit from either tRNAscan or EuFindtRNA
#[derive(Debug, Clone)]
pub struct FirstPassHit {
    pub start: usize,
    pub end: usize,
    pub strand: char,
    pub score: f64,
    pub source: String,
    pub isotype: String,
    pub anticodon: String,
}

impl FirstPassHit {
    /// Create from a tRNAscan hit
    pub fn from_trnascan(hit: &TrnascanHit) -> Self {
        Self {
            start: hit.start.max(1) as usize,
            end: hit.end.max(1) as usize,
            strand: hit.strand,
            score: 0.0, // tRNAscan doesn't provide a score
            source: "tRNAscan".to_string(),
            isotype: hit.isotype.clone(),
            anticodon: hit.anticodon.clone(),
        }
    }

    /// Create from an EuFindtRNA hit
    pub fn from_eufind(hit: &EufindHit) -> Self {
        Self {
            start: hit.start.max(1) as usize,
            end: hit.end.max(1) as usize,
            strand: if hit.start < hit.end { '+' } else { '-' },
            score: hit.score,
            source: "EuFindtRNA".to_string(),
            isotype: hit.isotype.clone(),
            anticodon: hit.anticodon.clone(),
        }
    }

    /// Create from basic coordinates
    pub fn new(start: usize, end: usize, score: f64, source: &str) -> Self {
        Self {
            start,
            end,
            strand: '+',
            score,
            source: source.to_string(),
            isotype: String::new(),
            anticodon: String::new(),
        }
    }
}

/// Main tRNA scanner
pub struct TrnaScanner {
    /// Scan mode (determines first-pass and CM selection)
    mode: ScanMode,
    /// Score cutoff for second-pass
    score_cutoff: f64,
    /// Quiet mode
    quiet: bool,
    /// Verbose mode
    verbose: bool,
    /// Show pseudogenes
    show_pseudogenes: bool,
    /// Accumulated results
    results: Vec<TrnaResult>,
    /// Accumulated faithful results (bacterial/archaeal/general `-B`/`-A`/`-G` path)
    trna_results: Vec<TRna>,

    // Model paths
    models_dir: PathBuf,

    // First-pass components for tRNAscan (Fichant-Burks)
    #[allow(dead_code)]
    tpc_matrix: Option<ConsensusMatrix>,
    #[allow(dead_code)]
    d_matrix: Option<ConsensusMatrix>,
    #[allow(dead_code)]
    search_params: SearchParams,

    // First-pass components for EuFindtRNA
    eufind_options: EufindOptions,

    // Second-pass: Infernal CM scanner
    cm_scanner: Option<CMScan>,

    /// `-H` (get_hmm_score): emit HMM Score + 2'Str Score columns and run the
    /// pseudogene / no-structure rescore for the faithful path.
    get_hmm_score: bool,
    /// `--detail`: emit Isotype CM + Isotype Score columns and the IPD/ISM note.
    detail: bool,
    /// `--no-isotype`: skip the isotype-specific scan and the Met-family
    /// (fMet/iMet/Ile2) Type refinement (C `no_isotype`).
    no_isotype: bool,
    /// Lazily-built no-structure (HMM) + per-isotype CM searchers (M4/M5).
    iso_res: RefCell<Option<IsoResources>>,
    /// Lazily-built Phase I/II Domain (+SeC) searchers. Cached because they are
    /// identical across source sequences — rebuilding per sequence dominated the
    /// `-B` cost on multi-sequence inputs.
    scan_searchers: RefCell<Option<Vec<(&'static str, FaithfulSearcher, bool)>>>,
    /// Lazily-built BHB noncanonical-intron CM searchers (archaeal only:
    /// `Cren-eury-BHB-noncan.cm` + `Thaum-BHB-noncan.cm`). `Some(vec![])` marks
    /// "already tried, none available" so we don't rebuild every sequence.
    bhb_searchers: RefCell<Option<Vec<FaithfulSearcher>>>,
}

/// No-structure (pseudogene/HMM) and per-isotype CM searchers used by the
/// faithful `-B -H --detail` path (spec §2.3 / §2.5). Built once, lazily.
struct IsoResources {
    /// No-structure model `TRNAinf-bact-ns.cm` (HMM Score / 2'Str rescore).
    ns: Option<FaithfulSearcher>,
    /// Per-isotype CMs: (basename without `bact-` prefix, searcher).
    iso: Vec<(String, FaithfulSearcher)>,
}

impl TrnaScanner {
    /// Create a new scanner with the specified mode and score cutoff
    pub fn new(mode: char, score_cutoff: f64) -> Result<Self, String> {
        Self::with_models_dir(mode, score_cutoff, Path::new("models"))
    }

    /// Create a scanner with a specific model path (for compatibility)
    pub fn with_model_path<P: AsRef<Path>>(
        mode: char,
        score_cutoff: f64,
        _model_path: P,
    ) -> Result<Self, String> {
        // For backward compatibility, derive models_dir from model_path
        let models_dir = _model_path.as_ref()
            .parent()
            .unwrap_or(Path::new("models"))
            .to_path_buf();
        Self::with_models_dir(mode, score_cutoff, &models_dir)
    }

    /// Create a new scanner with a models directory
    pub fn with_models_dir<P: AsRef<Path>>(
        mode: char,
        score_cutoff: f64,
        models_dir: P,
    ) -> Result<Self, String> {
        let scan_mode = ScanMode::from_char(mode);
        let models_path = models_dir.as_ref().to_path_buf();

        // Load tRNAscan matrices for prokaryotic modes
        let (tpc_matrix, d_matrix) = if scan_mode.uses_trnascan() {
            (Some(ConsensusMatrix::tpc_signal()), Some(ConsensusMatrix::d_signal()))
        } else {
            (None, None)
        };

        // Configure search parameters based on mode
        let search_params = match scan_mode {
            ScanMode::Bacterial | ScanMode::Archaeal => SearchParams::relaxed(),
            _ => SearchParams::strict(),
        };

        // Initialize CM scanner with main and SeC models
        let cm_model_path = models_path.join(scan_mode.cm_model_name());
        let cm_scanner = if cm_model_path.exists() {
            let mut scanner = CMScan::new();
            scanner.options = CMScanOptions {
                score_cutoff,
                ..CMScanOptions::default()
            };
            // Add the main CM path
            scanner.cm_files.add_main("Domain", cm_model_path);

            // Add SeC-specific CM if it exists (for bacterial/archaeal/eukaryotic)
            let sec_cm_name = match scan_mode {
                ScanMode::Bacterial => "TRNAinf-bact-SeC.cm",
                ScanMode::Archaeal => "TRNAinf-arch-SeC.cm",
                ScanMode::Eukaryotic => "TRNAinf-euk-SeC.cm",
                _ => "",
            };
            if !sec_cm_name.is_empty() {
                let sec_cm_path = models_path.join(sec_cm_name);
                if sec_cm_path.exists() {
                    scanner.cm_files.add_main("SeC", sec_cm_path);
                }
            }

            Some(scanner)
        } else {
            // CM file doesn't exist, will skip second-pass
            None
        };

        Ok(Self {
            mode: scan_mode,
            score_cutoff,
            quiet: false,
            verbose: false,
            show_pseudogenes: false,
            results: Vec::new(),
            trna_results: Vec::new(),
            models_dir: models_path,
            tpc_matrix,
            d_matrix,
            search_params,
            eufind_options: EufindOptions::default(),
            cm_scanner,
            get_hmm_score: false,
            detail: false,
            no_isotype: false,
            iso_res: RefCell::new(None),
            scan_searchers: RefCell::new(None),
            bhb_searchers: RefCell::new(None),
        })
    }

    /// Create a scanner without loading a model (for testing or when model is loaded later)
    pub fn new_without_model(mode: char, score_cutoff: f64) -> Self {
        let scan_mode = ScanMode::from_char(mode);

        let (tpc_matrix, d_matrix) = if scan_mode.uses_trnascan() {
            (Some(ConsensusMatrix::tpc_signal()), Some(ConsensusMatrix::d_signal()))
        } else {
            (None, None)
        };

        Self {
            mode: scan_mode,
            score_cutoff,
            quiet: false,
            verbose: false,
            show_pseudogenes: false,
            results: Vec::new(),
            trna_results: Vec::new(),
            models_dir: PathBuf::from("models"),
            tpc_matrix,
            d_matrix,
            search_params: SearchParams::default(),
            eufind_options: EufindOptions::default(),
            cm_scanner: None,
            get_hmm_score: false,
            detail: false,
            no_isotype: false,
            iso_res: RefCell::new(None),
            scan_searchers: RefCell::new(None),
            bhb_searchers: RefCell::new(None),
        }
    }

    /// Enable the `-H` HMM Score / 2'Str Score columns (faithful path).
    pub fn set_get_hmm_score(&mut self, on: bool) {
        self.get_hmm_score = on;
    }

    /// Enable the `--detail` Isotype CM / Isotype Score columns + IPD note.
    pub fn set_detail(&mut self, on: bool) {
        self.detail = on;
    }

    /// `--no-isotype`: disable the isotype scan + Met-family Type refinement.
    pub fn set_no_isotype(&mut self, on: bool) {
        self.no_isotype = on;
    }

    /// Isotype-specific scanning applies: the mode ships isotype model DBs
    /// (Bact/Arch/Euk) and it is not disabled. Mirrors C `!no_isotype()` — the
    /// Met-family (fMet/iMet/Ile2) refinement runs in the DEFAULT view;
    /// `--detail` only adds the Isotype CM / Score columns + IPD note.
    fn iso_applicable(&self) -> bool {
        !self.no_isotype
            && matches!(
                self.mode,
                ScanMode::Bacterial | ScanMode::Archaeal | ScanMode::Eukaryotic
            )
    }

    /// Set quiet mode
    pub fn set_quiet(&mut self, quiet: bool) {
        self.quiet = quiet;
    }

    /// Set verbose mode
    pub fn set_verbose(&mut self, verbose: bool) {
        self.verbose = verbose;
    }

    /// Set whether to show pseudogenes
    pub fn set_show_pseudogenes(&mut self, show: bool) {
        self.show_pseudogenes = show;
    }

    /// Get number of results
    pub fn result_count(&self) -> usize {
        self.results.len() + self.trna_results.len()
    }

    /// Borrow the faithful (`-B`/`-A`/`-G`) results accumulated across all
    /// scanned sequences. Populated when `uses_faithful()` is true.
    pub fn trna_results(&self) -> &[TRna] {
        &self.trna_results
    }

    /// Borrow the non-faithful (Cove / covels) results.
    pub fn results(&self) -> &[TrnaResult] {
        &self.results
    }

    /// Get current scan mode
    pub fn mode(&self) -> ScanMode {
        self.mode
    }

    /// First-pass scan using mode-appropriate method
    fn first_pass_scan(&self, seq: &[u8], seq_name: &str) -> Vec<FirstPassHit> {
        match self.mode {
            ScanMode::Bacterial | ScanMode::Archaeal => {
                // Use Infernal HMM-enabled first-pass for prokaryotes (like original tRNAscan-SE)
                self.run_infernal_first_pass(seq, seq_name)
            }
            ScanMode::Eukaryotic | ScanMode::Organellar | ScanMode::Mitochondrial => {
                // Use EuFindtRNA for eukaryotes
                self.run_eufind_first_pass(seq, seq_name)
            }
            ScanMode::General => {
                // Use Infernal first-pass
                self.run_infernal_first_pass(seq, seq_name)
            }
        }
    }

    /// Run Infernal HMM-enabled first-pass (like original tRNAscan-SE -B mode)
    fn run_infernal_first_pass(&self, seq: &[u8], _seq_name: &str) -> Vec<FirstPassHit> {
        let mut hits = Vec::new();

        // Whole-sequence ASCII string, searched in-process on both strands by the
        // faithful pipeline (no external cmsearch binary / temp files).
        let seq_str = String::from_utf8_lossy(seq).to_string();

        // Get CM models based on mode
        let cm_models = self.get_first_pass_cm_models();

        for cm_path in cm_models {
            if !cm_path.exists() {
                continue;
            }

            // Build the faithful searcher for this first-pass CM.
            let searcher = match FaithfulSearcher::from_cm_file(&cm_path) {
                Ok(s) => s,
                Err(_) => continue,
            };

            let seqs = [seq_str.as_str()];
            let cfg = FaithfulConfig::default();
            let fhits = searcher.search(&seqs, &cfg);

            for h in fhits {
                // First-pass candidates are strand-normalized to start <= end.
                let (start, end) = if h.start <= h.stop {
                    (h.start as usize, h.stop as usize)
                } else {
                    (h.stop as usize, h.start as usize)
                };

                hits.push(FirstPassHit {
                    start,
                    end,
                    strand: if h.in_rc { '-' } else { '+' },
                    score: h.score as f64,
                    source: "Infernal-FP".to_string(),
                    isotype: String::new(),
                    anticodon: String::new(),
                });
            }
        }

        // Merge overlapping hits
        self.merge_overlapping_hits(hits)
    }

    /// Get CM models for first-pass based on mode
    fn get_first_pass_cm_models(&self) -> Vec<std::path::PathBuf> {
        let mut models = Vec::new();

        // Main domain model
        let main_cm = self.models_dir.join(self.mode.cm_model_name());
        models.push(main_cm);

        // SeC model for bacterial/archaeal
        match self.mode {
            ScanMode::Bacterial => {
                models.push(self.models_dir.join("TRNAinf-bact-SeC.cm"));
            }
            ScanMode::Archaeal => {
                models.push(self.models_dir.join("TRNAinf-arch-SeC.cm"));
            }
            ScanMode::Eukaryotic => {
                models.push(self.models_dir.join("TRNAinf-euk-SeC.cm"));
            }
            _ => {}
        }

        models
    }

    /// Get isotype-specific CM database path based on mode
    fn get_isotype_cm_db(&self) -> Option<std::path::PathBuf> {
        let db_name = match self.mode {
            ScanMode::Bacterial => "TRNAinf-bact-iso",
            ScanMode::Archaeal => "TRNAinf-arch-iso",
            ScanMode::Eukaryotic => "TRNAinf-euk-iso",
            _ => return None,
        };
        let db_path = self.models_dir.join(db_name);
        if db_path.exists() {
            Some(db_path)
        } else {
            None
        }
    }

    /// Run cmscan to determine isotype and anticodon for results
    fn determine_isotypes(&self, results: &mut Vec<TrnaResult>, seq: &[u8], seq_name: &str) {
        use std::io::Write;
        use std::process::Command;
        use std::collections::HashMap;

        let iso_db = match self.get_isotype_cm_db() {
            Some(db) => db,
            None => {
                // Fallback: extract anticodon from sequence and use anticodon_to_isotype
                for result in results.iter_mut() {
                    self.extract_anticodon_from_seq(result, seq);
                }
                return;
            }
        };

        // Find cmscan binary
        let cmscan_path = std::env::var("CMSCAN_PATH")
            .unwrap_or_else(|_| {
                std::env::var("CMSEARCH_PATH")
                    .map(|p| p.replace("cmsearch", "cmscan"))
                    .unwrap_or_else(|_| "cmscan".to_string())
            });

        // Create temp files
        let temp_dir = std::env::temp_dir();
        let pid = std::process::id();
        let fasta_path = temp_dir.join(format!("trnascan_iso_{}.fa", pid));
        let tblout_path = temp_dir.join(format!("trnascan_iso_{}.tblout", pid));

        // Write tRNA sequences to temp FASTA
        {
            let mut fasta_file = match std::fs::File::create(&fasta_path) {
                Ok(f) => f,
                Err(_) => return,
            };

            for (i, result) in results.iter().enumerate() {
                let start = (result.begin.min(result.end) - 1).max(0) as usize;
                let end = (result.begin.max(result.end)) as usize;
                if end <= seq.len() {
                    let trna_seq = &seq[start..end];
                    let _ = writeln!(fasta_file, ">{}.t{}", seq_name, i + 1);
                    let _ = writeln!(fasta_file, "{}", String::from_utf8_lossy(trna_seq));
                }
            }
        }

        // Run cmscan (removed --toponly and --fmt 2 which cause empty output)
        let output = Command::new(&cmscan_path)
            .arg("-g")
            .arg("--mid")
            .arg("--notrunc")
            .arg("--tblout")
            .arg(&tblout_path)
            .arg(&iso_db)
            .arg(&fasta_path)
            .output();

        // Track best hit per tRNA (highest score)
        let mut best_hits: HashMap<usize, (String, f64)> = HashMap::new();

        if let Ok(run_result) = output {
            if run_result.status.success() {
                if let Ok(content) = std::fs::read_to_string(&tblout_path) {
                    for line in content.lines() {
                        if line.starts_with('#') || line.trim().is_empty() {
                            continue;
                        }
                        let fields: Vec<&str> = line.split_whitespace().collect();
                        // tblout format (without --fmt 2):
                        // 0: target name (bact-Ile)
                        // 1: accession
                        // 2: query name (seq.t1)
                        // ...
                        // 14: score
                        if fields.len() >= 15 {
                            let target_name = fields[0]; // e.g., "bact-Ile"
                            let query_name = fields[2];  // e.g., "NC_000913.3.t1"
                            let score: f64 = fields[14].parse().unwrap_or(0.0);

                            // Parse tRNA number from query name
                            if let Some(dot_pos) = query_name.rfind(".t") {
                                if let Ok(idx) = query_name[dot_pos + 2..].parse::<usize>() {
                                    // Keep only best hit (highest score) per tRNA
                                    let is_better = best_hits.get(&idx)
                                        .map(|(_, old_score)| score > *old_score)
                                        .unwrap_or(true);
                                    if is_better {
                                        best_hits.insert(idx, (target_name.to_string(), score));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // Apply best hits to results
        for (idx, (target_name, _score)) in best_hits {
            if idx > 0 && idx <= results.len() {
                let result = &mut results[idx - 1];

                // Parse isotype from target name (e.g., "bact-Ile" -> "Ile")
                let parts: Vec<&str> = target_name.split('-').collect();
                if parts.len() >= 2 {
                    let isotype = parts[1].to_string();
                    if result.isotype.is_empty() || result.isotype == "???" {
                        result.isotype = isotype.clone();
                        result.cm_isotype = isotype;
                    }
                }
            }
        }

        // Extract anticodon from sequence position 34-36 for all results
        for result in results.iter_mut() {
            self.extract_anticodon_from_seq(result, seq);
        }

        // Clean up
        let _ = std::fs::remove_file(&fasta_path);
        let _ = std::fs::remove_file(&tblout_path);
    }

    /// Extract anticodon from tRNA sequence using pattern matching
    /// The anticodon loop has conserved U33 and A37 positions,
    /// so we look for T[NNN]A pattern in the middle region of the tRNA
    fn extract_anticodon_from_seq(&self, result: &mut TrnaResult, seq: &[u8]) {
        if !result.anticodon.is_empty() && result.anticodon != "???" && result.anticodon != "NNN" {
            return; // Already has anticodon
        }

        let start = (result.begin.min(result.end) - 1).max(0) as usize;
        let end = (result.begin.max(result.end)) as usize;

        if end > seq.len() {
            return;
        }

        let trna_seq = &seq[start..end];
        let trna_len = trna_seq.len();

        if trna_len < 60 {
            return; // Too short
        }

        // Search for T[NNN]A pattern in the anticodon loop region (middle third of tRNA)
        // The pattern is U33-anticodon(34-36)-A37 in standard tRNA numbering
        let search_start = trna_len / 4;
        let search_end = (trna_len * 3) / 4;

        let mut best_ac: Option<(usize, Vec<u8>)> = None;

        for i in search_start..search_end.saturating_sub(4) {
            // Look for T at position i (U33) and A at position i+4 (A37)
            if trna_seq[i] == b'T' && trna_seq[i + 4] == b'A' {
                let anticodon = trna_seq[i + 1..i + 4].to_vec();
                // Prefer positions closer to the expected location (~position 33)
                let expected_pos = 32;
                let dist = (i as i32 - expected_pos as i32).abs();
                if best_ac.is_none() || dist < best_ac.as_ref().unwrap().0 as i32 {
                    // Verify this looks like a valid anticodon (all valid nucleotides)
                    let valid = anticodon.iter().all(|&b| b == b'A' || b == b'C' || b == b'G' || b == b'T');
                    if valid {
                        best_ac = Some((dist as usize, anticodon));
                    }
                }
            }
        }

        // If no T...A pattern found, fall back to position-based extraction
        let anticodon = if let Some((_, ac)) = best_ac {
            if result.strand == '-' {
                reverse_complement(&ac)
            } else {
                ac
            }
        } else {
            // Fallback: use position 33-35 (0-indexed)
            let ac_start = if trna_len >= 70 && trna_len <= 95 {
                33
            } else {
                (trna_len as f64 * 0.44) as usize
            };
            let ac_end = (ac_start + 3).min(trna_len);
            if ac_end <= trna_len {
                if result.strand == '-' {
                    reverse_complement(&trna_seq[trna_len - ac_end..trna_len - ac_start])
                } else {
                    trna_seq[ac_start..ac_end].to_vec()
                }
            } else {
                return;
            }
        };

        let ac_str = String::from_utf8_lossy(&anticodon).to_uppercase();
        result.anticodon = ac_str.clone();

        // If isotype is empty, derive from anticodon
        if result.isotype.is_empty() || result.isotype == "???" {
            if let Some(iso) = anticodon_to_isotype(&ac_str) {
                result.isotype = iso.to_string();
                result.cm_isotype = iso.to_string();
            }
        }
    }

    /// Run tRNAscan (Fichant-Burks) first-pass
    #[allow(dead_code)]
    fn run_trnascan_first_pass(&self, seq: &[u8]) -> Vec<FirstPassHit> {
        if let (Some(tpc), Some(d)) = (&self.tpc_matrix, &self.d_matrix) {
            let hits = trnascan(seq, tpc, d, &self.search_params);
            hits.iter().map(|h| FirstPassHit::from_trnascan(h)).collect()
        } else {
            Vec::new()
        }
    }

    /// Run EuFindtRNA first-pass
    fn run_eufind_first_pass(&self, seq: &[u8], seq_name: &str) -> Vec<FirstPassHit> {
        let hits = run_eufind_scan(seq, seq_name, &self.eufind_options);
        hits.iter().map(|h| FirstPassHit::from_eufind(h)).collect()
    }

    /// Merge overlapping first-pass hits
    fn merge_overlapping_hits(&self, mut hits: Vec<FirstPassHit>) -> Vec<FirstPassHit> {
        if hits.is_empty() {
            return hits;
        }

        // Sort by position
        hits.sort_by_key(|h| (h.start, h.end));

        let mut merged = Vec::new();
        let mut current = hits[0].clone();

        for hit in hits.into_iter().skip(1) {
            // Check for overlap
            if hit.start <= current.end + 20 && hit.strand == current.strand {
                // Merge hits, keeping the one with higher score
                current.end = current.end.max(hit.end);
                if hit.score > current.score {
                    current.score = hit.score;
                    current.source = hit.source;
                }
            } else {
                merged.push(current);
                current = hit;
            }
        }
        merged.push(current);

        merged
    }

    /// Find the candidate that matches a CM search hit by region name
    fn find_candidate_for_hit<'a>(
        &self,
        region_name: &str,
        candidate_offsets: &'a [(usize, &FirstPassHit)],
    ) -> Option<(usize, &'a FirstPassHit)> {
        // Parse region name format: "seqname:start-end"
        if let Some(colon_pos) = region_name.rfind(':') {
            let coords = &region_name[colon_pos + 1..];
            if let Some(dash_pos) = coords.find('-') {
                if let Ok(start) = coords[..dash_pos].parse::<usize>() {
                    // Find matching candidate by start position (convert to 0-indexed)
                    let offset_start = start.saturating_sub(1);
                    for (offset, candidate) in candidate_offsets {
                        if *offset == offset_start {
                            return Some((*offset, *candidate));
                        }
                    }
                }
            }
        }
        None
    }

    /// Second-pass scan using Infernal cmsearch
    fn second_pass_scan(
        &mut self,
        candidates: &[FirstPassHit],
        seq: &[u8],
        seq_name: &str,
    ) -> Vec<TrnaResult> {
        let mut results = Vec::new();

        if let Some(ref cm_scanner) = self.cm_scanner {
            // Prepare batch of sequences for cmsearch
            let padding = 20;
            let mut batch_sequences: Vec<(String, String)> = Vec::new();
            let mut candidate_offsets: Vec<(usize, &FirstPassHit)> = Vec::new();

            for candidate in candidates {
                let start = candidate.start.saturating_sub(padding);
                let end = (candidate.end + padding).min(seq.len());

                if start >= end || end > seq.len() {
                    continue;
                }

                let subseq = &seq[start..end];
                let subseq_str = String::from_utf8_lossy(subseq).to_string();
                let region_name = format!("{}:{}-{}", seq_name, start + 1, end);

                batch_sequences.push((region_name, subseq_str));
                candidate_offsets.push((start, candidate));
            }

            // Run batch cmsearch for each CM model (Domain and SeC)
            if !batch_sequences.is_empty() {
                // Search with main Domain CM
                let cm_keys = ["Domain", "SeC"];
                for cm_key in cm_keys {
                    // Check if this CM exists
                    if cm_scanner.cm_files.get_main(cm_key).is_none() {
                        continue;
                    }

                    match cm_scanner.batch_search_external(cm_key, &batch_sequences) {
                        Ok(hits) => {
                            for hit in hits {
                                if hit.score >= self.score_cutoff || self.show_pseudogenes {
                                    // Parse region name to get offset and find matching candidate
                                    if let Some((offset, candidate)) = self.find_candidate_for_hit(
                                        &hit.target_name,
                                        &candidate_offsets,
                                    ) {
                                        let mut result = TrnaResult::from_cm_hit(&hit, seq_name, offset);
                                        result.strand = candidate.strand;

                                        // For SeC model, set isotype appropriately
                                        if cm_key == "SeC" && result.isotype.is_empty() {
                                            result.isotype = "SeC".to_string();
                                            result.cm_isotype = "SeC".to_string();
                                        }

                                        // Inherit isotype/anticodon from first-pass if CM didn't detect
                                        if result.isotype.is_empty() || result.isotype == "???" {
                                            if !candidate.isotype.is_empty() {
                                                result.isotype = candidate.isotype.clone();
                                                result.cm_isotype = candidate.isotype.clone();
                                            }
                                        }
                                        if result.anticodon.is_empty() || result.anticodon == "???" {
                                            if !candidate.anticodon.is_empty() {
                                                result.anticodon = candidate.anticodon.clone();
                                            }
                                        }

                                        // Try to determine isotype from anticodon if still unknown
                                        if (result.isotype.is_empty() || result.isotype == "???")
                                           && !result.anticodon.is_empty()
                                           && result.anticodon != "???" {
                                            if let Some(iso) = anticodon_to_isotype(&result.anticodon) {
                                                result.isotype = iso.to_string();
                                                result.cm_isotype = iso.to_string();
                                            }
                                        }

                                        results.push(result);
                                    }
                                }
                            }
                        }
                        Err(e) => {
                            if self.verbose {
                                eprintln!("Batch CM search error for {}: {}", cm_key, e);
                            }
                        }
                    }
                }
            }
        } else {
            // No CM scanner available, create results from first-pass hits
            // When no CM verification, output all first-pass hits (they've already passed
            // the first-pass thresholds, so score filtering isn't appropriate here)
            for candidate in candidates {
                // For first-pass only mode, include all candidates
                // (tRNAscan doesn't provide a numeric score, so filter by isotype instead)
                let has_valid_isotype = !candidate.isotype.is_empty()
                    && candidate.isotype != "???"
                    && candidate.isotype != "Unk";

                if has_valid_isotype || self.show_pseudogenes {
                    // Extract sequence
                    let start = candidate.start.saturating_sub(1);
                    let end = candidate.end.min(seq.len());
                    let subseq = if start < end { seq[start..end].to_vec() } else { Vec::new() };

                    let isotype = if !candidate.isotype.is_empty() && candidate.isotype != "???" {
                        candidate.isotype.clone()
                    } else if !candidate.anticodon.is_empty() && candidate.anticodon != "???" {
                        anticodon_to_isotype(&candidate.anticodon)
                            .map(|s| s.to_string())
                            .unwrap_or_else(|| "Unk".to_string())
                    } else {
                        "Unk".to_string()
                    };

                    results.push(TrnaResult {
                        seq_name: seq_name.to_string(),
                        trna_num: 0,
                        begin: candidate.start as i64,
                        end: candidate.end as i64,
                        isotype: isotype.clone(),
                        anticodon: if candidate.anticodon.is_empty() { "NNN".to_string() } else { candidate.anticodon.clone() },
                        intron_begin: 0,
                        intron_end: 0,
                        inf_score: 0.0,  // No CM score available
                        hmm_score: 0.0,
                        ss_score: 0.0,
                        origin: candidate.source.clone(),
                        cm_isotype: isotype,
                        cm_score: 0.0,  // No CM score available
                        note: "first-pass".to_string(),
                        sequence: subseq,
                        structure: String::new(),
                        strand: candidate.strand,
                    });
                }
            }
        }

        results
    }

    /// Whether this mode uses the faithful in-process Infernal `-B`/`-A`/`-G`
    /// pipeline (Phase I candidate scan -> Phase II global-nohmm verify -> decode).
    pub fn uses_faithful(&self) -> bool {
        matches!(
            self.mode,
            ScanMode::Bacterial | ScanMode::Archaeal | ScanMode::General
        )
    }

    /// Build the (role, searcher, is_sec_cm) list once for this scan (Domain then
    /// SeC), reused across Phase I (local) and Phase II (global-nohmm) searches.
    fn build_faithful_searchers(&self) -> Vec<(&'static str, FaithfulSearcher, bool)> {
        // Main domain CM + optional SeC CM. Both builds are independent, so run
        // them concurrently (rayon::join) — construction dominates on short inputs.
        let main_cm = self.models_dir.join(self.mode.cm_model_name());
        let sec_name = match self.mode {
            ScanMode::Bacterial => "TRNAinf-bact-SeC.cm",
            ScanMode::Archaeal => "TRNAinf-arch-SeC.cm",
            ScanMode::Eukaryotic => "TRNAinf-euk-SeC.cm",
            _ => "",
        };
        let sec_cm = if sec_name.is_empty() {
            None
        } else {
            let p = self.models_dir.join(sec_name);
            if p.exists() { Some(p) } else { None }
        };

        let (main, sec) = rayon::join(
            || {
                if main_cm.exists() {
                    FaithfulSearcher::from_cm_file(&main_cm).ok()
                } else {
                    None
                }
            },
            || sec_cm.and_then(|p| FaithfulSearcher::from_cm_file(&p).ok()),
        );

        let mut out = Vec::new();
        if let Some(s) = main {
            out.push(("Domain", s, false));
        }
        if let Some(s) = sec {
            out.push(("SeC", s, true));
        }
        out
    }

    /// Build (and lazily cache) the Phase I/II Domain (+SeC) searchers.
    fn ensure_scan_searchers(&self) {
        if self.scan_searchers.borrow().is_some() {
            return;
        }
        let s = self.build_faithful_searchers();
        *self.scan_searchers.borrow_mut() = Some(s);
    }

    /// Faithful `-B` pipeline for one source sequence: returns sorted, numbered
    /// `TRna` records (byte-parity intent with the C/Perl reference 9-col `.out`).
    fn faithful_scan_sequence(&self, seq: &[u8], seqname: &str, seqlen: usize) -> Vec<TRna> {
        self.ensure_scan_searchers();
        let searchers_ref = self.scan_searchers.borrow();
        let searchers = match searchers_ref.as_ref() {
            Some(s) if !s.is_empty() => s,
            _ => return Vec::new(),
        };

        let whole = String::from_utf8_lossy(seq).to_string();

        // ---- Phase I: candidate prescan (C flag-3: `-g --mid --notrunc -T 10`,
        // both strands) ----
        // C runs the first pass in GLOBAL --mid mode with the fixed
        // infernal_fp_cutoff bit-score threshold (10), NOT local/E-value. The
        // first-pass hit boundaries (via --mid) then define the +/-10-padded
        // Phase-II window, so the mode must match C exactly or the candidate
        // window (and thus the final 5'/3' boundary) diverges on aberrant tRNAs.
        let fp_cfg = FaithfulConfig {
            toponly: false,
            global: true,
            mid: true,
            notrunc: true,
            t_cutoff: Some(10.0),
            ..Default::default()
        };
        let mut candidates: Vec<FirstPassHit> = Vec::new();
        for (_role, searcher, _is_sec) in searchers.iter() {
            let fhits = searcher.search(&[whole.as_str()], &fp_cfg);
            for h in fhits {
                let (start, end) = if h.start <= h.stop {
                    (h.start as usize, h.stop as usize)
                } else {
                    (h.stop as usize, h.start as usize)
                };
                // C `process_fp_cmsearch_hits` (CM.pm:3163-3164) stores each
                // first-pass hit as `cmsearch_coord + start_index`, where
                // start_index=1 for a whole contig scanned in one buffer. This +1
                // shifts the candidate window one base toward 3' (the whole tRNA
                // moves +1, then it is padded +/-10). The shift is inert for a
                // well-centred tRNA but decisive at a 5' edge (e.g. a His whose
                // acceptor stem scores higher when position -1 is included): C's
                // window excludes that base, so cmsearch cannot over-extend 5'.
                let start = start + 1;
                let end = std::cmp::min(seqlen, end + 1);
                candidates.push(FirstPassHit {
                    start,
                    end,
                    strand: if h.in_rc { '-' } else { '+' },
                    score: h.score as f64,
                    source: "Infernal-FP".to_string(),
                    isotype: String::new(),
                    anticodon: String::new(),
                });
            }
        }
        let candidates = self.merge_overlapping_hits(candidates);

        // ---- Phase II: per-candidate global-nohmm verify + decode ----
        // Region = first-pass candidate padded by exactly +/-10 (C
        // default_Padding, FpScanResultFile.pm:351-352), clamped to [1, seqlen].
        // This +/-10 window IS the entire Phase-II input in C — there is no extra
        // flank. A wider window lets Phase-II extend the 5'/3' boundary past what
        // C can see on aberrant tRNAs, so it must be exactly 10. Search both
        // strands (bit scores are strand-independent; the C `--toponly` on a
        // strand-oriented extraction yields identical scores/alignments).
        const FLANK: i64 = 10; // C default_Padding (FpScanResultFile.pm:123)
        let flank: i64 = FLANK;
        // C's post-cmsearch 3' boundary trim/extend + rescore (CM.pm:3685-3757)
        // runs for the bacterial/archaeal main models only (euk/general excluded).
        let apply_boundary = matches!(self.mode, ScanMode::Bacterial | ScanMode::Archaeal);
        let mode = self.mode;
        let cfg = FaithfulConfig {
            toponly: false,
            e_report: 10.0,
            global: true,
            nohmm: true,
            ..Default::default()
        };

        // Each candidate's Phase-II verify (nohmm CYK+Inside on its ~250 bp region,
        // × Domain/SeC) is independent; the downstream dedup sorts by score, so the
        // hit order is irrelevant. Fan out over candidates — this is the dominant
        // Phase-II cost on real genomes (hundreds of candidates).
        let score_cutoff = self.score_cutoff;
        let hits: Vec<TRna> = candidates
            .par_iter()
            .flat_map_iter(|cand| {
                let mut local: Vec<TRna> = Vec::new();
                let region_lo = std::cmp::max(1i64, cand.start as i64 - flank);
                let region_hi = std::cmp::min(seqlen as i64, cand.end as i64 + flank);
                if region_lo >= region_hi {
                    return local.into_iter();
                }
                let region = &seq[(region_lo - 1) as usize..region_hi as usize];
                let region_str = String::from_utf8_lossy(region).to_string();

                for (role, searcher, is_sec) in searchers.iter() {
                    let fhits = searcher.search(&[region_str.as_str()], &cfg);
                    for h in fhits {
                        if (h.score as f64) < score_cutoff {
                            continue;
                        }
                        let ali = match &h.alignment {
                            Some(a) => a,
                            None => continue,
                        };

                        // Remap region-frame coords back to genomic.
                        let (gstart, gend, strand) = if h.in_rc {
                            // h.start > h.stop (both region-frame, forward numbering)
                            let gs = region_lo - 1 + h.stop;
                            let ge = region_lo - 1 + h.start;
                            (gs, ge, TStrand::Minus)
                        } else {
                            let gs = region_lo - 1 + h.start;
                            let ge = region_lo - 1 + h.stop;
                            (gs, ge, TStrand::Plus)
                        };

                        let adisp = AliDisplay {
                            aseq: ali.aseq.clone(),
                            ss_cons: ali.csline.clone(),
                            nc: ali.ncline.clone(),
                            model: ali.model.clone(),
                        };
                        let dec = decode_trna_properties(
                            &adisp, role, *is_sec, false, strand, gstart, gend,
                        );

                        let mut t = TRna::new();
                        t.seqname = seqname.to_string();
                        t.strand = strand;
                        t.start = gstart;
                        t.end = gend;
                        t.score = h.score as f64;
                        t.set_domain_model("infernal", h.score as f64);
                        t.isotype = dec.isotype.clone();
                        t.anticodon = dec.anticodon.clone();
                        t.model = role.to_string();
                        t.hit_source = "Inf".to_string();
                        t.src_seqlen = seqlen;
                        if let Some(intron) = dec.intron {
                            t.introns.push(intron);
                        }

                        // C's post-cmsearch 3' boundary trim/extend + rescore
                        // (CM.pm:3685-3757). Adjusts the 3' coord (and rescored
                        // score) BEFORE dedup so coords/score reflect the boundary.
                        // The anticodon/isotype/intron are all on the 5' side and
                        // are unaffected by the 3' adjustment, so `dec` stays valid.
                        // Final C-format seq/ss (`cm_tRNA->seq()`/`ss()`), after
                        // any 3' boundary trim/extend and 5' isotype fixes.
                        let (fin_seq, fin_ss): (Vec<u8>, Vec<u8>) = if apply_boundary {
                            let (mut aseq, mut ass) = Self::apply_boundary_adjust(
                                &mut t,
                                &dec.norm_seq,
                                &dec.norm_ss,
                                seq,
                                searcher,
                            );
                            // C's 5' isotype boundary fixes, run on the post-trim
                            // seq/ss (analyze_with_cmsearch CM.pm:3784-3785).
                            // fix_fMet: bacterial + Met + score>40; fix_His:
                            // archaeal + His + score>35.
                            match mode {
                                ScanMode::Bacterial => Self::fix_fmet(
                                    &mut t, &mut aseq, &mut ass, seq, searcher,
                                ),
                                ScanMode::Archaeal => {
                                    Self::fix_his(&mut t, &mut aseq, &mut ass, searcher)
                                }
                                _ => {}
                            }
                            (aseq, ass)
                        } else {
                            (dec.norm_seq.clone().into_bytes(), dec.norm_ss.clone().into_bytes())
                        };
                        Self::populate_seq_fields(&mut t, fin_seq, fin_ss);

                        local.push(t);
                    }
                }
                local.into_iter()
            })
            .collect();

        // ---- Cross-model / cross-candidate dedup: keep higher score on overlap ----
        let mut deduped = self.dedup_faithful_hits(hits);

        // ---- Output sort + per-sequence numbering (spec 3.4) ----
        Self::sort_faithful(&mut deduped);
        for (i, t) in deduped.iter_mut().enumerate() {
            t.id = i + 1;
        }

        // ---- Noncanonical (BHB) intron re-scan (archaeal only) ----
        // C: tRNAscan-SE.src:441 scan_noncanonical_introns, run BEFORE isotype
        // decoration. Adds a bulge-helix-bulge intron at a noncanonical position
        // that the anticodon-loop-only `find_intron` misses.
        self.scan_noncanonical_introns(&mut deduped, seq);

        // ---- Truncated tRNA search (C: tRNAscan-SE.src:449, CM.pm:2653) ----
        // Re-scan each found tRNA's MATURE sequence with a truncation-allowed
        // search (infernox `notrunc:false` = C cmsearch scan_flag 6 `-g --toponly`,
        // the only pass that drops `--notrunc`). check_truncation (CM.pm:2718)
        // then labels tRNAs whose best (Inside-selected) parse is a truncated-CM
        // alignment (trunc_start:N / trunc_end:N). Runs between the noncanonical-
        // intron rescan and isotype decoration, for euk/bact/arch (SeC uses its CM).
        self.apply_truncation(&mut deduped, seq);

        // ---- M4/M5: HMM Score + 2'Str + Isotype CM/Score + Note ----
        self.decorate_faithful(&mut deduped, seq);

        deduped
    }

    /// Truncated tRNA search (C: CM.pm truncated_tRNA_search 2653 + check_truncation
    /// 2718). Re-scans each found tRNA's MATURE sequence with a truncation-allowed
    /// search (`notrunc:false`, the infernox port of C's `-g --toponly` flag-6) and,
    /// when the best (Inside-selected) parse is a truncated-CM alignment, records the
    /// `trunc_start:N` / `trunc_end:N` label + direction. Only euk/bact/arch run this
    /// (driver src:445); each tRNA is scored against its own model (Domain vs SeC).
    fn apply_truncation(&self, trnas: &mut [TRna], seq: &[u8]) {
        if !matches!(
            self.mode,
            ScanMode::Bacterial | ScanMode::Archaeal | ScanMode::Eukaryotic
        ) {
            return;
        }
        self.ensure_scan_searchers();
        let searchers_ref = self.scan_searchers.borrow();
        let searchers = match searchers_ref.as_ref() {
            Some(s) if !s.is_empty() => s,
            _ => return,
        };
        // flag-6: `-g --toponly` with truncation ENABLED (notrunc:false). The mature
        // sequence is coding-oriented, so the tRNA is on the top strand. Default
        // E-value reporting (<=10) matches C (cm_cutoff > 10 adds no -T).
        let cfg = FaithfulConfig {
            toponly: true,
            e_report: 10.0,
            global: true,
            notrunc: false,
            ..Default::default()
        };
        // CCA 3'-exclusion applies for non-euk / non-general Domain models (CM.pm:2742).
        let is_euk_general = matches!(self.mode, ScanMode::General);

        trnas.par_iter_mut().for_each(|t| {
            let span = Self::faithful_span_seq(t, seq);
            if span.is_empty() {
                return;
            }
            let mature = Self::faithful_mature_seq(t, &span);
            if mature.is_empty() {
                return;
            }
            // Score against the tRNA's own model (Domain / SeC), matching C's
            // per-model write_tRNAs loop.
            let searcher = searchers
                .iter()
                .find(|(role, _, _)| *role == t.model)
                .or_else(|| searchers.iter().find(|(role, _, _)| *role == "Domain"))
                .map(|(_, s, _)| s);
            let searcher = match searcher {
                Some(s) => s,
                None => return,
            };
            let hits = searcher.search(&[mature.as_str()], &cfg);
            // Best-scoring hit (the pipeline reports the highest-scoring parse).
            let best = hits.iter().fold(None::<&infernal::FaithfulHit>, |acc, h| {
                match acc {
                    Some(a) if a.score >= h.score => Some(a),
                    _ => Some(h),
                }
            });
            if let Some(h) = best {
                let label = Self::check_truncation(h, &mature, is_euk_general);
                if !label.is_empty() {
                    let has5 = label.contains("trunc_start");
                    let has3 = label.contains("trunc_end");
                    t.trunc = match (has5, has3) {
                        (true, true) => Truncation::Both,
                        (true, false) => Truncation::FivePrime,
                        (false, true) => Truncation::ThreePrime,
                        (false, false) => Truncation::None,
                    };
                    t.trunc_label = label;
                }
            }
        });
    }

    /// C: check_truncation (CM.pm:2718). Given a truncation-allowed hit on the
    /// mature sequence, build the `trunc_start:N` / `trunc_end:N` label from the
    /// alignment's model span. N = cfrom_emit-cfrom_span (5', `<[N]*`) and
    /// cto_span-cto_emit (3', `*[N]>`). The 3' side carries a CCA exclusion:
    /// a small (<=3) 3'-only overhang on a non-CCA-tailed, non-5'-truncated
    /// bact/arch tRNA is NOT labelled (it is the missing CCA, not a truncation).
    fn check_truncation(h: &infernal::FaithfulHit, mature: &str, is_euk_general: bool) -> String {
        let mut label = String::new();
        if h.trunc.is_empty() || h.trunc == "no" {
            return label;
        }
        let ali = match &h.alignment {
            Some(a) => a,
            None => return label,
        };
        let is5 = h.trunc.contains("5'");
        let is3 = h.trunc.contains("3'");
        // 5' truncation (C: /^\<\[(\d+)\]\*/): N = cfrom_emit - cfrom_span.
        if is5 {
            let n = h.mdl_from - ali.cfrom_span;
            if n > 0 {
                label = format!("trunc_start:{}", n);
            }
        }
        // 3' truncation (C: /\*\[(\d+)\]\>$/): diff = cto_span - cto_emit.
        if is3 {
            let diff = ali.cto_span - h.mdl_to;
            if diff > 0 {
                // CCA exclusion (CM.pm:2740-2744).
                let mat_ends_cca = mature.len() >= 3
                    && mature[mature.len() - 3..].eq_ignore_ascii_case("CCA");
                let skip = diff <= 3
                    && ((!mat_ends_cca && !is5) || is5)
                    && !is_euk_general;
                if !skip {
                    if !label.is_empty() {
                        label.push(',');
                    }
                    label.push_str(&format!("trunc_end:{}", diff));
                }
            }
        }
        label
    }

    /// Isotype-CM reporting cutoff (`isotype_cm_cutoff.bact`, spec §2.5).
    fn isotype_cutoff(&self) -> f64 {
        20.0
    }

    /// Round a bit score to one decimal, matching the `%.1f` cmsearch/cmscan
    /// tblout representation the reference parses and stores. All downstream
    /// arithmetic (2'Str = Inf − HMM, IPD = own − highest) uses these rounded
    /// values, so the `%.2f` columns render e.g. `37.30`, not `37.26`.
    fn round1(x: f64) -> f64 {
        (x * 10.0).round() / 10.0
    }

    /// Build (and lazily cache) the no-structure (HMM) + per-isotype searchers.
    fn ensure_iso_res(&self) {
        if self.iso_res.borrow().is_some() {
            return;
        }
        let res = self.build_iso_resources();
        *self.iso_res.borrow_mut() = Some(res);
    }

    /// Load `TRNAinf-<clade>-ns.cm` (HMM rescore) and split `TRNAinf-<clade>-iso`
    /// into its per-isotype CMs (each `INFERNAL1…//` CM + its paired HMMER3
    /// filter), loading a `FaithfulSearcher` for each.
    fn build_iso_resources(&self) -> IsoResources {
        let (ns_name, iso_name, prefix) = match self.mode {
            ScanMode::Bacterial => ("TRNAinf-bact-ns.cm", "TRNAinf-bact-iso", "bact-"),
            ScanMode::Archaeal => ("TRNAinf-arch-ns.cm", "TRNAinf-arch-iso", "arch-"),
            ScanMode::Eukaryotic => ("TRNAinf-euk-ns.cm", "TRNAinf-euk-iso", "euk-"),
            _ => return IsoResources { ns: None, iso: Vec::new() },
        };

        // No-structure model path (searcher built below, in parallel with the
        // isotype fleet).
        let ns_path = {
            let p = self.models_dir.join(ns_name);
            if p.exists() { Some(p) } else { None }
        };

        // Isotype CM database → split → per-model temp CM files. The split and
        // temp-file writes are cheap I/O; the expensive part is the per-model
        // FaithfulSearcher build, done in parallel below.
        let mut iso_files: Vec<(String, PathBuf)> = Vec::new();
        let iso_path = self.models_dir.join(iso_name);
        if let Ok(text) = std::fs::read_to_string(&iso_path) {
            let tmp_dir = std::env::temp_dir()
                .join(format!("trnascan_iso_{}", std::process::id()));
            let _ = std::fs::create_dir_all(&tmp_dir);

            // Split on each `INFERNAL1…` record start; every segment carries one
            // CM block plus its trailing HMMER3 filter block.
            let mut segments: Vec<String> = Vec::new();
            let mut cur = String::new();
            for line in text.lines() {
                if line.starts_with("INFERNAL1") && !cur.is_empty() {
                    segments.push(std::mem::take(&mut cur));
                }
                cur.push_str(line);
                cur.push('\n');
            }
            if !cur.is_empty() {
                segments.push(cur);
            }

            for seg in segments {
                let name = seg.lines().find_map(|l| {
                    let l = l.trim_start();
                    if l.starts_with("NAME") {
                        l.split_whitespace().nth(1).map(|s| s.to_string())
                    } else {
                        None
                    }
                });
                let name = match name {
                    Some(n) => n,
                    None => continue,
                };
                let basename = name
                    .strip_prefix(prefix)
                    .unwrap_or(&name)
                    .to_string();
                let file = tmp_dir.join(format!("{}.cm", name));
                if std::fs::write(&file, seg.as_bytes()).is_err() {
                    continue;
                }
                iso_files.push((basename, file));
            }
        }

        // Build every searcher in parallel — the ns model concurrently with the
        // isotype fleet, and the isotype models among themselves. Each build
        // (CP9 + global CM + filters + QDB bands + consensus) is the dominant
        // one-time cost of `-H`/`--detail`. `collect` preserves `iso_files`
        // order, so the downstream tie-break fold stays deterministic.
        let (ns, iso) = rayon::join(
            || ns_path.and_then(|p| FaithfulSearcher::from_cm_file(&p).ok()),
            || {
                iso_files
                    .par_iter()
                    .filter_map(|(basename, file)| {
                        FaithfulSearcher::from_cm_file(file)
                            .ok()
                            .map(|s| (basename.clone(), s))
                    })
                    .collect::<Vec<(String, FaithfulSearcher)>>()
            },
        );

        IsoResources { ns, iso }
    }

    /// Strand-corrected genomic span (5'→3' coding sequence, intron included).
    fn faithful_span_seq(t: &TRna, seq: &[u8]) -> String {
        let (s, e) = (t.start as usize, t.end as usize);
        if s < 1 || e > seq.len() || s > e {
            return String::new();
        }
        let sub = &seq[s - 1..e];
        let span = if t.strand == TStrand::Minus {
            reverse_complement(sub)
        } else {
            sub.to_vec()
        };
        String::from_utf8_lossy(&span).to_string()
    }

    /// Mature sequence (intron excised) from the coding span. `rel_start`/
    /// `rel_end` are 1-based within the coding span.
    fn faithful_mature_seq(t: &TRna, span: &str) -> String {
        if let Some(intron) = t.introns.first() {
            let rs = intron.rel_start as usize;
            let re = intron.rel_end as usize;
            if rs >= 1 && re <= span.len() && rs <= re {
                let mut m = String::with_capacity(span.len());
                m.push_str(&span[..rs - 1]);
                m.push_str(&span[re..]);
                return m;
            }
        }
        span.to_string()
    }

    /// Store the final C-format `seq`/`ss` (`cm_tRNA->seq()`/`ss()`), the mature
    /// (intron-excised) `mat_seq`/`mat_ss`, and the anticodon relative positions
    /// onto `t`. The anticodon position is recomputed on the FINAL seq/ss so that
    /// any 3' boundary trim / 5' isotype shift is reflected; genomic mapping in the
    /// struct writer (`t.start + rel_start - 1`) is then self-consistent with `seq`.
    fn populate_seq_fields(t: &mut TRna, seq: Vec<u8>, ss: Vec<u8>) {
        let seq_s = String::from_utf8_lossy(&seq).into_owned();
        let ss_s = String::from_utf8_lossy(&ss).into_owned();

        t.ac_positions.clear();
        if t.anticodon != crate::cm_scan::decode::UNDEF_ANTICODON {
            let (ac, _ai, _ae, ac_pos) = find_anticodon(&seq_s, &ss_s);
            if ac_pos > 0 && ac == t.anticodon {
                t.ac_positions.push(AnticodonPos {
                    rel_start: ac_pos,
                    rel_end: ac_pos + 2,
                });
            }
        }

        // Mature seq/ss = coding seq/ss with the canonical intron excised
        // (rel_start..rel_end, 1-based inclusive). ASCII → byte == char index.
        let splice = |s: &str| -> String {
            if let Some(intron) = t.introns.first() {
                let rs = intron.rel_start as usize;
                let re = intron.rel_end as usize;
                if rs >= 1 && re <= s.len() && rs <= re {
                    let mut m = String::with_capacity(s.len());
                    m.push_str(&s[..rs - 1]);
                    m.push_str(&s[re..]);
                    return m;
                }
            }
            s.to_string()
        };
        t.mat_seq = splice(&seq_s);
        t.mat_ss = splice(&ss_s);
        t.seq = seq_s;
        t.ss = ss_s;
    }

    /// Max bit score over all `-g --nohmm` hits of `seq` under `searcher`
    /// (E-value reporting relaxed so low-scoring rescore hits survive).
    /// Max bit score of any hit of `searcher` on `seq` under `cfg`.
    fn best_score(searcher: &FaithfulSearcher, seq: &str, cfg: &FaithfulConfig) -> Option<f64> {
        if seq.is_empty() {
            return None;
        }
        searcher
            .search(&[seq], cfg)
            .iter()
            .map(|h| h.score as f64)
            .fold(None, |acc, s| Some(acc.map_or(s, |a: f64| a.max(s))))
    }

    /// The 3 coding nucleotides immediately 3' of the mature tRNA end, uppercased,
    /// or `None` if they fall outside the source sequence. `start`/`end` are the
    /// 1-based genomic low/high coords; the 3' end is `end` on the plus strand and
    /// `start` on the minus strand. Used by Rule A (CM.pm:3686) to test for a
    /// genomic `CCA` just past the hit.
    fn next3_coding_3prime(genome: &[u8], start: i64, end: i64, plus: bool) -> Option<[u8; 3]> {
        #[inline]
        fn comp(b: u8) -> u8 {
            match b.to_ascii_uppercase() {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                b'U' => b'A',
                other => other,
            }
        }
        if plus {
            // Genomic 1-based positions end+1, end+2, end+3 -> 0-based end, end+1, end+2.
            let i = end as usize;
            if i + 3 > genome.len() {
                return None;
            }
            Some([
                genome[i].to_ascii_uppercase(),
                genome[i + 1].to_ascii_uppercase(),
                genome[i + 2].to_ascii_uppercase(),
            ])
        } else {
            // Coding continues to lower genomic coords: revcomp of genomic 1-based
            // positions start-1, start-2, start-3 (0-based start-2, start-3, start-4).
            if start < 4 {
                return None;
            }
            Some([
                comp(genome[(start - 2) as usize]),
                comp(genome[(start - 3) as usize]),
                comp(genome[(start - 4) as usize]),
            ])
        }
    }

    /// Faithful port of the CM.pm `analyze_with_cmsearch` 3' boundary trim/extend +
    /// rescore (lines 3685-3757), applied to bacterial/archaeal main-model hits only
    /// (the caller gates on `apply_boundary`; euk/general are excluded there).
    ///
    /// `norm_seq`/`norm_ss` are the C-format `cm_tRNA->seq()`/`ss()` strings
    /// (post-`format_cmsearch_output`: brackets swapped so `<` = a paired *close*,
    /// unpaired = `.`, lowercase = insert-state residue). Mutates `t.start`/`t.end`
    /// (3' coord) and, if any rule fired, re-runs the flag-7 (`-g --max --toponly
    /// --notrunc -T 0`) search on the adjusted mature sequence to refresh `t.score`.
    ///
    /// Returns the final adjusted `(seq, ss)` (C-format) so the caller can feed
    /// them to the 5' isotype fixes (`fix_fMet` / `fix_His`), which C runs on the
    /// post-boundary-trim seq/ss (`analyze_with_cmsearch` CM.pm:3784-3785).
    fn apply_boundary_adjust(
        t: &mut TRna,
        norm_seq: &str,
        norm_ss: &str,
        genome: &[u8],
        searcher: &FaithfulSearcher,
    ) -> (Vec<u8>, Vec<u8>) {
        #[inline]
        fn is_ins(b: u8) -> bool {
            matches!(b, b'a' | b'c' | b'g' | b't' | b'n')
        }
        #[inline]
        fn is_mat(b: u8) -> bool {
            matches!(b, b'A' | b'C' | b'G' | b'T' | b'N')
        }

        let mut seq: Vec<u8> = norm_seq.as_bytes().to_vec();
        let mut ss: Vec<u8> = norm_ss.as_bytes().to_vec();
        let plus = t.strand == TStrand::Plus;
        let mut rescore = false;

        // ---- Rule A: extend +3 (append genomic CCA) — CM.pm:3685-3703 ----
        // Last 4 struct chars NOT all unpaired AND the next 3 genomic coding nt == CCA.
        if ss.len() >= 4 && &ss[ss.len() - 4..] != b"...." {
            if let Some(cca) = Self::next3_coding_3prime(genome, t.start, t.end, plus) {
                if &cca == b"CCA" {
                    if plus {
                        t.end += 3;
                    } else {
                        t.start -= 3;
                    }
                    seq.extend_from_slice(b"CCA");
                    ss.extend_from_slice(b"...");
                    rescore = true;
                }
            }
        }

        let last3_is_cca = |s: &[u8]| s.len() >= 3 && s[s.len() - 3..].eq_ignore_ascii_case(b"CCA");
        let last4_all_dot = ss.len() >= 4 && &ss[ss.len() - 4..] == b"....";

        // ---- Rule B: trim -3 — CM.pm:3705-3722 ----
        // seq[-3:] != CCA AND ss[-4:] == "....".
        if !last3_is_cca(&seq) && last4_all_dot {
            if plus {
                t.end -= 3;
            } else {
                t.start += 3;
            }
            seq.truncate(seq.len() - 3);
            ss.truncate(ss.len() - 3);
            rescore = true;
        }
        // ---- Rule C: trim -4 or -5 — CM.pm:3723-3752 (elsif) ----
        // seq[-3:] == CCA AND one of three stem-remnant patterns at the 3' end.
        else if last3_is_cca(&seq) {
            let n = ss.len();
            let m = seq.len();
            // p1: ss[-6:] == "<....." AND seq[-5:-3] == [acgtn][ACGTN]
            let p1 = n >= 6
                && &ss[n - 6..] == b"<....."
                && m >= 5
                && is_ins(seq[m - 5])
                && is_mat(seq[m - 4]);
            // p2: ss[-7:] == "<......" AND seq[-6:-3] == [ACGTN][acgtn][ACGTN]
            let p2 = n >= 7
                && &ss[n - 7..] == b"<......"
                && m >= 6
                && is_mat(seq[m - 6])
                && is_ins(seq[m - 5])
                && is_mat(seq[m - 4]);
            // p3: ss[-7:] == "<......" AND seq[-6:-3] == [acgtn]{2}[ACGTN]  (trim 5)
            let p3 = n >= 7
                && &ss[n - 7..] == b"<......"
                && m >= 6
                && is_ins(seq[m - 6])
                && is_ins(seq[m - 5])
                && is_mat(seq[m - 4]);
            if p1 || p2 || p3 {
                let trim_len: usize = if p3 { 5 } else { 4 };
                if plus {
                    t.end -= trim_len as i64;
                } else {
                    t.start += trim_len as i64;
                }
                seq.truncate(seq.len() - trim_len);
                // C uppercases the new final base (CM.pm:3749); harmless for scoring.
                if let Some(last) = seq.last_mut() {
                    *last = last.to_ascii_uppercase();
                }
                ss.truncate(ss.len() - trim_len);
                rescore = true;
            }
        }

        // ---- rescore — CM.pm:3754-3757 -> rescore_tRNA -> cmsearch_scoring flag 7 ----
        if rescore {
            let adj: String = seq.iter().map(|&b| b as char).collect();
            if let Some(sc) = Self::ns_max_score(searcher, &adj) {
                t.score = sc;
                t.set_domain_model("infernal", sc);
            }
        }

        (seq, ss)
    }

    /// Faithful port of CM.pm `fix_fMet` (:1352-1433). BACTERIAL only. When the
    /// hit's isotype is `Met` and its (post-boundary-trim) infernal domain score
    /// > 40, C nudges the 5' boundary by ±1 base under two mutually-exclusive
    /// secondary-structure / sequence conditions, then rescores (flag-7 search on
    /// the adjusted mature seq). `seq`/`ss` are the C-format post-boundary seq/ss
    /// returned by `apply_boundary_adjust`. Mutates `t.start`/`t.end`/`t.score`.
    ///
    /// The reported anticodon (a string) is unchanged by a 1-base 5' shift, and
    /// bacterial Met tRNAs carry no intron, so C's `ar_ac_pos` bookkeeping (which
    /// only feeds relative anticodon/intron reporting) has no observable effect on
    /// the 9-column output here; the intron genomic coords, when present, are
    /// stored absolute and unaffected. Only start/end/score are mutated.
    fn fix_fmet(
        t: &mut TRna,
        seq: &mut Vec<u8>,
        ss: &mut Vec<u8>,
        genome: &[u8],
        searcher: &FaithfulSearcher,
    ) {
        let old = t.score;
        if Self::fix_fmet_transform(t, seq, ss, genome) {
            let adj: String = seq.iter().map(|&b| b as char).collect();
            if let Some(sc) = Self::ns_max_score(searcher, &adj) {
                t.score = sc;
                t.set_domain_model("infernal", sc);
            }
            if std::env::var("FIX_DEBUG").is_ok() {
                eprintln!(
                    "FIX_FMET fired: {} {}-{} score {}->{}",
                    t.seqname, t.start, t.end, old, t.score
                );
            }
        }
    }

    /// Pure geometric transform of `fix_fMet` (no rescore): applies the C 5'
    /// boundary edit to `t.start`/`t.end` + `seq`/`ss`, returning `true` if a rule
    /// fired. Split out so the exact C conditions can be unit-tested without a
    /// searcher. Gate (isotype/score) is checked here to mirror C.
    fn fix_fmet_transform(
        t: &mut TRna,
        seq: &mut Vec<u8>,
        ss: &mut Vec<u8>,
        genome: &[u8],
    ) -> bool {
        if t.isotype != "Met" || t.score <= 40.0 {
            return false;
        }
        let plus = t.strand == TStrand::Plus;
        let mut rescore = false;

        fn last3_is_cca(s: &[u8]) -> bool {
            s.len() >= 3 && &s[s.len() - 3..] == b"CCA"
        }
        fn last5(s: &[u8]) -> &[u8] {
            if s.len() >= 5 { &s[s.len() - 5..] } else { s }
        }

        // Outer gate (CM.pm:1363-1364): (seq[-3:]=="CCA" and ss[-5:]==".....") or
        // ss[-5:]=="<<<..".
        let cca_dot5 = last3_is_cca(seq) && last5(ss) == b".....";
        let stem5 = last5(ss) == b"<<<..";
        if !(cca_dot5 || stem5) {
            return false;
        }

        if !ss.is_empty() && ss[0] != b'.' {
            // Branch 1 (CM.pm:1366-1392): prepend the upstream coding base, extend
            // the 5' boundary by 1.
            let can = (plus && t.start > 1)
                || (!plus && (t.end as usize) < t.src_seqlen);
            if can {
                if let Some(b) = Self::upstream_coding_base(genome, t.start, t.end, plus) {
                    let mut ns: Vec<u8> = Vec::with_capacity(seq.len() + 1);
                    ns.push(b);
                    ns.extend_from_slice(seq);
                    *seq = ns;
                    let mut nss: Vec<u8> = Vec::with_capacity(ss.len() + 1);
                    nss.push(b'.');
                    nss.extend_from_slice(ss);
                    *ss = nss;
                    if plus {
                        t.start -= 1;
                    } else {
                        t.end += 1;
                    }
                    rescore = true;
                }
            }
        } else if ss.len() >= 4 && &ss[0..4] == b".>.>" && seq.len() >= 2 && &seq[0..2] == b"CG" {
            // Branch 2 (CM.pm:1393-1425): remove a 5' bulge base, trim 5' by 1, but
            // only when the base 3' of the bulge is the Watson-Crick partner of the
            // 3'-side pos71 base (a genuine extra base, not a real residue).
            // pos71 index depends on which 3' pattern matched (CM.pm:1396-1403).
            let n = ss.len();
            let pos71: Option<u8> = if cca_dot5 {
                // seq[len(ss)-6]
                if n >= 6 { seq.get(n - 6).copied() } else { None }
            } else {
                // stem5: seq[len(ss)-3]
                if n >= 3 { seq.get(n - 3).copied() } else { None }
            };
            if let Some(p71) = pos71 {
                let target = seq[2].to_ascii_uppercase();
                if Self::rev_comp_base(p71.to_ascii_uppercase()) == target {
                    // seq = seq[1] + uc(seq[2]) + seq[3..]  (drop index0, uc old idx2)
                    let mut ns: Vec<u8> = Vec::with_capacity(seq.len() - 1);
                    ns.push(seq[1]);
                    ns.push(seq[2].to_ascii_uppercase());
                    ns.extend_from_slice(&seq[3..]);
                    *seq = ns;
                    // ss = ss[0..2] + ss[3..]  (drop index2)
                    let mut nss: Vec<u8> = Vec::with_capacity(ss.len() - 1);
                    nss.extend_from_slice(&ss[0..2]);
                    nss.extend_from_slice(&ss[3..]);
                    *ss = nss;
                    if plus {
                        t.start += 1;
                    } else {
                        t.end -= 1;
                    }
                    rescore = true;
                }
            }
        }

        rescore
    }

    /// Faithful port of CM.pm `fix_His` (:1436-1479). ARCHAEAL only. When the
    /// hit's isotype is `His` and its (post-boundary-trim) infernal domain score
    /// > 35, and the SS matches `>>>>.>>>.` … `<<<.<<<<.` with a valid base pair
    /// at pos5/pos68, C removes the spurious extra His 5' G bulge (drops one base
    /// from each end), shifts start+1 / end-1, and rescores. `seq`/`ss` are the
    /// C-format post-boundary seq/ss. Mutates `t.start`/`t.end`/`t.score`.
    ///
    /// As with `fix_fMet`, the anticodon string is unchanged and archaeal His
    /// carries no intron, so C's relative `ar_ac_pos` shift has no observable
    /// effect on the 9-column output; only start/end/score are mutated.
    fn fix_his(
        t: &mut TRna,
        seq: &mut Vec<u8>,
        ss: &mut Vec<u8>,
        searcher: &FaithfulSearcher,
    ) {
        let old = t.score;
        if Self::fix_his_transform(t, seq, ss) {
            let adj: String = seq.iter().map(|&b| b as char).collect();
            if let Some(sc) = Self::ns_max_score(searcher, &adj) {
                t.score = sc;
                t.set_domain_model("infernal", sc);
            }
            if std::env::var("FIX_DEBUG").is_ok() {
                eprintln!(
                    "FIX_HIS fired: {} {}-{} score {}->{}",
                    t.seqname, t.start, t.end, old, t.score
                );
            }
        }
    }

    /// Pure geometric transform of `fix_His` (no rescore), returning `true` if the
    /// rule fired. Split out for unit-testing the exact C conditions.
    fn fix_his_transform(t: &mut TRna, seq: &mut Vec<u8>, ss: &mut Vec<u8>) -> bool {
        if t.isotype != "His" || t.score <= 35.0 {
            return false;
        }
        let n = ss.len();
        let m = seq.len();
        // C-format SS 5'/3' gate (CM.pm:1450): ss[0..9]==">>>>.>>>." and
        // ss[-9:]=="<<<.<<<<.". Requires len >= 11 for the mid slice below.
        if m < 11 || n < 11 {
            return false;
        }
        if &ss[0..9] != b">>>>.>>>." || &ss[n - 9..] != b"<<<.<<<<." {
            return false;
        }
        let pos5 = seq[4].to_ascii_uppercase();
        let pos68 = seq[m - 6].to_ascii_uppercase();
        // Base-pair test (CM.pm:1456-1458): AT/TA/GC/CG/GT/TG (T-form, wobble ok).
        let paired = matches!(
            (pos5, pos68),
            (b'A', b'T') | (b'T', b'A') | (b'G', b'C') | (b'C', b'G') | (b'G', b'T') | (b'T', b'G')
        );
        if !paired {
            return false;
        }
        // mid = seq[5 .. m-6]  (length m-11)
        let mid = &seq[5..m - 6];
        // seq = seq[1..4] + pos5 + mid + pos68 + seq[m-5..m-1]  (drop idx0 & last)
        let mut nseq: Vec<u8> = Vec::with_capacity(m - 2);
        nseq.extend_from_slice(&seq[1..4]);
        nseq.push(pos5);
        nseq.extend_from_slice(mid);
        nseq.push(pos68);
        nseq.extend_from_slice(&seq[m - 5..m - 1]);
        // ss = ss[1..4] + ">" + ss[5..n-6] + "<" + ss[n-5..n-2] + "."
        let mut nss: Vec<u8> = Vec::with_capacity(n - 2);
        nss.extend_from_slice(&ss[1..4]);
        nss.push(b'>');
        nss.extend_from_slice(&ss[5..n - 6]);
        nss.push(b'<');
        nss.extend_from_slice(&ss[n - 5..n - 2]);
        nss.push(b'.');
        *seq = nseq;
        *ss = nss;
        t.start += 1;
        t.end -= 1;
        true
    }

    /// Reverse-complement of a single (uppercased) base, T-form (port of
    /// `rev_comp_seq` restricted to one nt). Non-ACGT returns itself.
    #[inline]
    fn rev_comp_base(b: u8) -> u8 {
        match b {
            b'A' => b'T',
            b'T' | b'U' => b'A',
            b'G' => b'C',
            b'C' => b'G',
            other => other,
        }
    }

    /// The single coding base immediately 5' of the mature tRNA start (C's
    /// `substr($trna->upstream(), -1)`), or `None` if it falls outside the source
    /// sequence. On `+` this is genomic `start-1`; on `-` it is the complement of
    /// genomic `end+1` (coding continues to higher genomic coords).
    fn upstream_coding_base(genome: &[u8], start: i64, end: i64, plus: bool) -> Option<u8> {
        if plus {
            if start < 2 {
                return None;
            }
            genome.get((start - 2) as usize).map(|b| b.to_ascii_uppercase())
        } else {
            let i = end as usize; // 0-based index of genomic (end+1)
            genome.get(i).map(|&b| Self::rev_comp_base(b.to_ascii_uppercase()))
        }
    }

    /// M4 pseudogene NS rescore — C flag 7: `-g --max --toponly --notrunc -T 0`.
    /// `--max` (non-banded) is required: `--nohmm` QDB banding drops weak hits
    /// (e.g. MySeq5's 12.7) that C's `-g --max` recovers.
    fn ns_max_score(searcher: &FaithfulSearcher, seq: &str) -> Option<f64> {
        // toponly: the span is extracted in the tRNA's coding orientation, so the
        // hit is always on the top (+) strand — the bottom strand only yields
        // lower spurious hits that never win the max. Scoring the top strand alone
        // is byte-identical and halves the (non-banded --max) work.
        let cfg = FaithfulConfig {
            toponly: true,
            e_report: 1e9,
            global: true,
            max: true,
            t_cutoff: Some(0.0),
            ..Default::default()
        };
        Self::best_score(searcher, seq, &cfg)
    }

    /// M5 isotype scan — C flag 2: `-g --mid --toponly --notrunc` (cmscan).
    /// `--mid` (HMM-banded global final Inside) matches the golden isotype
    /// scores; `--nohmm` (QDB) differs by ~0.1 bit on some models (e.g. His).
    fn iso_mid_score(searcher: &FaithfulSearcher, seq: &str) -> Option<f64> {
        // toponly: the mature span is coding-oriented (see ns_max_score) — top strand
        // holds the isotype hit; scoring one strand is byte-identical and halves work.
        let cfg = FaithfulConfig {
            toponly: true,
            e_report: 1e9,
            global: true,
            mid: true,
            ..Default::default()
        };
        Self::best_score(searcher, seq, &cfg)
    }

    /// M4 (HMM Score / 2'Str) + M5 (Isotype CM / Score / IPD note) decoration.
    fn decorate_faithful(&self, trnas: &mut [TRna], seq: &[u8]) {
        // Isotype refinement runs in the DEFAULT view (C gates it on !no_isotype,
        // not on --detail). In the default 9-col view only Met-family (CAT) tRNAs
        // can change Type, so the isotype scan there is limited to Met hits; a run
        // with no Met tRNA skips the (expensive) iso-searcher build entirely.
        let iso_on = self.iso_applicable();
        let has_met = iso_on && trnas.iter().any(|t| t.isotype == "Met");
        let need_iso = self.detail && iso_on || has_met;
        // Pseudogene filter (C is_pseudo_gene, CM.pm:999) runs for any tRNA with
        // Inf score < 55 even in the default view, so we must not early-return when
        // such candidates exist — they still need the NS rescore + pseudo check.
        let has_pseudo_cand = trnas.iter().any(|t| t.model != "SeC" && t.score < 55.0);
        if !self.get_hmm_score && !need_iso && !has_pseudo_cand {
            return;
        }
        self.ensure_iso_res();
        let res_ref = self.iso_res.borrow();
        let res = match res_ref.as_ref() {
            Some(r) => r,
            None => return,
        };
        let cutoff = self.isotype_cutoff();
        let get_hmm_score = self.get_hmm_score;
        let detail = self.detail;

        // Each tRNA's decoration is independent; the NS rescore (M4) and the 23
        // isotype rescores (M5) are the pipeline's dominant cost (non-banded --max /
        // HMM-banded --mid DP). Fan out over tRNAs, and over the isotype models
        // within each, with rayon. FaithfulSearcher scores via `&self` only, so the
        // shared `res` searchers are Sync; results are order-independent.
        trnas.par_iter_mut().for_each(|t| {
            let span = Self::faithful_span_seq(t, seq);
            if span.is_empty() {
                return;
            }
            let mature = Self::faithful_mature_seq(t, &span);
            let is_sec = t.model == "SeC";

            // ---- M4: HMM Score + 2'Str Score (spec §2.3) + pseudogene filter ----
            // C (is_pseudo_gene, CM.pm:999): for Inf score < 55 the NS rescore +
            // pseudo check runs ALWAYS (even without -H) — the fast-path skip only
            // fires when score >= 55 AND !get_hmm_score. With -H it runs for every
            // tRNA to populate the HMM/2'Str columns. SeC is exempt (score >= 55).
            if !is_sec && (get_hmm_score || t.score < 55.0) {
                let hmm = Self::round1(
                    res.ns
                        .as_ref()
                        .and_then(|ns| Self::ns_max_score(ns, &span))
                        .unwrap_or(0.0),
                );
                let ss = Self::round1(t.score) - hmm;
                // Pseudogene: (ss_score < min_ss_score 5 OR hmm_score < min_hmm_score
                // 10) AND Inf score < min_pseudo_filter_score 55.
                if (ss < 5.0 || hmm < 10.0) && t.score < 55.0 {
                    t.is_pseudo = true;
                }
                // The HMM/2'Str columns are only displayed under -H.
                if get_hmm_score {
                    t.hmm_score = hmm;
                    t.ss_score = ss;
                }
            } else if get_hmm_score && is_sec {
                // SeC hits skip the pseudogene / NS rescore; -H shows 0.00/0.00.
                t.hmm_score = 0.0;
                t.ss_score = 0.0;
            }

            // ---- M5: Isotype scan + Met-family Type refinement (spec §2.5) ----
            // Default view: only Met tRNAs need the scan (that's the only Type that
            // the isotype models can change). --detail: scan every tRNA for the
            // extra Isotype CM / Score columns + IPD note.
            let want_iso = iso_on && (detail || t.isotype == "Met");
            let mut recorded: Vec<(String, f64)> = Vec::new();
            if want_iso {
                // Parallel over the 23 isotype models; `collect` preserves `res.iso`
                // order so the tie-break fold below is identical to the serial loop.
                recorded = res
                    .iso
                    .par_iter()
                    .filter_map(|(name, s)| {
                        Self::iso_mid_score(s, &mature).and_then(|sc| {
                            let sc = Self::round1(sc);
                            if sc >= cutoff {
                                Some((name.clone(), sc))
                            } else {
                                None
                            }
                        })
                    })
                    .collect();
                // Highest-scoring model (ties: keep the first / alphabetically
                // earlier, matching the sorted split order).
                let best = recorded.iter().cloned().fold(
                    None,
                    |acc: Option<(String, f64)>, (m, s)| match acc {
                        Some((_, bs)) if bs >= s => acc,
                        _ => Some((m, s)),
                    },
                );
                if let Some((model, score)) = best {
                    // Met/Ile2 promotion (ScanResult.pm:311-327).
                    if t.isotype == "Met"
                        && (model == "iMet" || model == "fMet" || model == "Ile2")
                    {
                        t.isotype = model.clone();
                    } else if t.isotype == "Met"
                        && model != "Met"
                        && model != "iMet"
                        && model != "fMet"
                    {
                        let met = recorded
                            .iter()
                            .find(|(m, _)| m == "Met")
                            .map(|(_, s)| *s)
                            .unwrap_or(0.0);
                        let ile2 = recorded
                            .iter()
                            .find(|(m, _)| m == "Ile2")
                            .map(|(_, s)| *s)
                            .unwrap_or(0.0);
                        if ile2 > 0.0
                            && met > 0.0
                            && (score - ile2) <= 5.0
                            && (ile2 - met) >= 5.0
                            && t.score > 50.0
                        {
                            t.isotype = "Ile2".to_string();
                        }
                    }
                    // Isotype CM / Score are --detail-only output columns; the
                    // default view already got its Type refinement above.
                    if detail {
                        t.iso_model = model;
                        t.iso_score = score;
                    }
                }
            }

            // ---- Note: pseudo, then IPD/ISM (spec §2.5, ScanResult.pm:758-864) ----
            let mut note = String::new();
            if t.is_pseudo {
                note.push_str("pseudo");
            }
            if detail
                && !t.iso_model.is_empty()
                && !t.isotype.is_empty()
                && t.isotype != "Undet"
                && t.iso_model != t.isotype
            {
                let alias = (t.iso_model == "iMet"
                    || t.iso_model == "fMet"
                    || t.iso_model == "Ile2")
                    && t.isotype == "Met";
                if !alias {
                    let own = recorded
                        .iter()
                        .find(|(m, _)| m == &t.isotype)
                        .map(|(_, s)| *s)
                        .unwrap_or(0.0);
                    if !note.is_empty() {
                        note.push(',');
                    }
                    note.push_str(&format!("IPD:{:.2}", own - t.iso_score));
                }
            }
            // Truncation label (C: construct_tab_output 792-798). The main `.out`
            // Note shows the `trunc_start:N`/`trunc_end:N` label only under
            // --detail; the struct writer shows "Possible truncation" always
            // (derived from `t.trunc`). So the default 9-col output is unchanged.
            if detail && t.is_trunc() && !t.trunc_label.is_empty() {
                if !note.is_empty() {
                    note.push(',');
                }
                note.push_str(&t.trunc_label);
            }
            t.note = note;
        });
    }

    /// Build (and lazily cache) the archaeal BHB noncanonical-intron searchers
    /// (`Cren-eury-BHB-noncan.cm`, `Thaum-BHB-noncan.cm`). Only archaeal mode
    /// ships these models (C: CM.pm:441-446 `nci_cm`). Returns an empty cache
    /// (never `None`) once attempted, so we don't re-scan the models dir.
    fn ensure_bhb_searchers(&self) {
        if self.bhb_searchers.borrow().is_some() {
            return;
        }
        let mut v = Vec::new();
        if self.mode == ScanMode::Archaeal {
            for name in ["Cren-eury-BHB-noncan.cm", "Thaum-BHB-noncan.cm"] {
                let p = self.models_dir.join(name);
                if p.exists() {
                    if let Ok(s) = FaithfulSearcher::from_cm_file(&p) {
                        v.push(s);
                    }
                }
            }
        }
        *self.bhb_searchers.borrow_mut() = Some(v);
    }

    /// Faithful port of `CM.pm::scan_noncanonical_introns` (:1598) — the full
    /// two-round BHB (bulge-helix-bulge) noncanonical intron driver.
    ///
    /// For each already-decoded tRNA, run the BHB intron CM(s) (`-g --max
    /// --toponly --notrunc -T 6.5`, C `scan_flag == 5`) over the tRNA's coding
    /// span padded with 70 nt of coding up/downstream. Round 1 processes every
    /// BHB hit through [`Self::check_intron_validity`]; each accepted intron is
    /// clipped out and its length accumulated. If any intron was accepted (or the
    /// tRNA had a canonical intron but no BHB hit), Round 2 re-searches the
    /// clipped sequence so a second intron whose BHB structure only appears after
    /// the first is removed can be found (C :1673-1709).
    ///
    /// `check_intron_validity` is the acceptance gate: SS-regex parse, mature
    /// re-score with the main Domain (+SeC) CM (`-g --nohmm`, `scan_flag == 0`),
    /// CCA trim, intron-in-precursor location, duplicate / inclusion / overlap
    /// rejection, `hit_overlap` (40) and `score > tRNA.score && mature ≥ 70`.
    ///
    /// After both rounds, if any noncanonical intron was accepted the tRNA is
    /// re-decoded ([`Self::decode_nci_trna_properties`]): all introns are spliced
    /// out and `find_anticodon` is re-run on the clean mature, recovering the true
    /// anticodon / isotype that the intron-garbled loop hid. Otherwise the tRNA is
    /// restored to its pristine pre-scan copy (C `$tRNA_copy`).
    fn scan_noncanonical_introns(&self, trnas: &mut [TRna], seq: &[u8]) {
        self.ensure_bhb_searchers();
        let bhb_ref = self.bhb_searchers.borrow();
        let bhb = match bhb_ref.as_ref() {
            Some(v) if !v.is_empty() => v,
            _ => return,
        };
        // Main Domain (+SeC) searchers for the mature re-score (already built).
        self.ensure_scan_searchers();
        let scan_ref = self.scan_searchers.borrow();
        let main = match scan_ref.as_ref() {
            Some(v) if !v.is_empty() => v,
            _ => return,
        };

        const FLANK: i64 = 70; // C upstream_len / downstream_len
        let seqlen = seq.len() as i64;
        // BHB search config: C scan_flag 5 = `-g --max --toponly --notrunc -T 6.5`.
        let bhb_cfg = FaithfulConfig {
            toponly: true,
            e_report: 1e9,
            global: true,
            max: true,
            t_cutoff: Some(6.5),
            ..Default::default()
        };
        // Mature re-score config: C scan_flag 0 = `-g --nohmm --toponly --notrunc`.
        let main_cfg = FaithfulConfig {
            toponly: true,
            e_report: 1e9,
            global: true,
            nohmm: true,
            ..Default::default()
        };

        for t in trnas.iter_mut() {
            // Coding span (5'→3', intron included) + 70 nt coding flanks
            // (C: tRNA->upstream()/seq()/downstream()).
            let span = Self::faithful_span_seq(t, seq);
            if span.is_empty() {
                continue;
            }
            let (up_seq, dn_seq) = Self::coding_flanks(seq, t, FLANK, seqlen);

            // Pristine copy for the "no noncanonical intron" write path (C
            // `$tRNA_copy`, CM.pm:1638/1751).
            let t_copy = t.clone();

            // Working flank/precursor state; updated on each accepted intron
            // (C: tRNA->upstream()/seq()/downstream() get overwritten, :2090-2092).
            // The precursor is the raw genomic coding span, UPPERCASED (C
            // uppercases input, so `tRNA->seq()` at NCI-scan time is uppercase
            // genomic even in soft-masked regions): the reconstructed mature is
            // thus uppercase, matching C's archaeal struct `Seq:` for introns that
            // do NOT carry an anticodon-loop `find_intron` (the mixed case seen on
            // those that DO comes from re-splicing the clip-rescore mature in
            // `decode_nci_trna_properties`).
            let mut cur_up = up_seq.to_uppercase();
            let mut cur_seq = span.to_uppercase();
            let mut cur_dn = dn_seq.to_uppercase();

            // Clip-space search target (C `$padded_seq`); starts = full padded seq.
            let mut padded_seq = format!("{}{}{}", cur_up, cur_seq, cur_dn);
            let mut prev_len: i64 = 0;
            let mut rnd2 = false;
            let mut add_ci = false;

            // ---- Round 1 (C CM.pm:1656-1670) ----
            let r1 = Self::bhb_hits(bhb, &padded_seq, &bhb_cfg);
            if r1.is_empty() {
                if t.introns.is_empty() {
                    // No BHB hit, no intron: write tRNA unchanged (C :1642-1646).
                    continue;
                }
                // No BHB hit but an existing (canonical) intron: re-search the
                // MATURE sequence in round 2 (C :1648-1652).
                let mat = Self::faithful_mature_seq(t, &cur_seq);
                padded_seq = format!("{}{}{}", cur_up, mat, cur_dn);
                rnd2 = true;
                add_ci = true;
            } else {
                for h in &r1 {
                    let padded_full = format!("{}{}{}", cur_up, cur_seq, cur_dn);
                    if let Some((dup, clip, ilen)) = self.check_intron_validity(
                        t, h, &padded_seq, &padded_full, prev_len, &mut cur_up,
                        &mut cur_seq, &mut cur_dn, main, &main_cfg,
                    ) {
                        padded_seq = clip;
                        prev_len += ilen;
                        rnd2 = true;
                        if dup {
                            add_ci = true;
                        }
                    }
                }
            }

            // ---- Round 2 (C CM.pm:1673-1709): re-search the clipped seq ----
            if rnd2 {
                let r2 = Self::bhb_hits(bhb, &padded_seq, &bhb_cfg);
                prev_len = 0;
                for h in &r2 {
                    let padded_full = format!("{}{}{}", cur_up, cur_seq, cur_dn);
                    if let Some((dup, clip, ilen)) = self.check_intron_validity(
                        t, h, &padded_seq, &padded_full, prev_len, &mut cur_up,
                        &mut cur_seq, &mut cur_dn, main, &main_cfg,
                    ) {
                        padded_seq = clip;
                        prev_len += ilen;
                        if dup {
                            add_ci = true;
                        }
                    }
                }
            }

            // ---- CI/NCI reconciliation (C CM.pm:1711-1752) ----
            let nci_count = t.introns.iter().filter(|i| i.intron_type == "NCI").count();
            if nci_count > 0 {
                let ci_index = t.introns.iter().position(|i| i.intron_type == "CI");
                let mut ci_seq = String::new();
                if let Some(ci) = ci_index {
                    if t.model != "SeC" {
                        if add_ci {
                            ci_seq = t.introns[ci].seq.clone();
                        }
                        t.introns.remove(ci);
                    }
                }
                t.introns.sort_by_key(|i| i.rel_start);
                if add_ci {
                    self.add_canonical_intron(t, &ci_seq, main, &main_cfg);
                }
                self.decode_nci_trna_properties(t);
            } else {
                // No noncanonical intron: restore the untouched tRNA (C :1751).
                *t = t_copy;
            }
        }
    }

    /// `Utils.pm::seg_overlap` (:156). Returns true if the two 1-based inclusive
    /// segments overlap; with `range > 0` an endpoint within `±range` of an
    /// endpoint also counts (used for the tRNA/hit 40-nt overlap test).
    fn seg_overlap(a1: i64, b1: i64, a2: i64, b2: i64, range: i64) -> bool {
        if range == 0 {
            (a1 >= a2 && a1 <= b2)
                || (b1 >= a2 && b1 <= b2)
                || (a2 >= a1 && a2 <= b1)
                || (b2 >= a1 && b2 <= b1)
        } else {
            (a1 >= a2 - range && a1 <= a2 + range)
                || (b1 >= b2 - range && b1 <= b2 + range)
                || (a2 >= a1 - range && a2 <= a1 + range)
                || (b2 >= b1 - range && b2 <= b1 + range)
        }
    }

    /// Run the BHB intron CM(s) over `target` and return the hits after C's
    /// per-merge overlap dedup (`merge_overlapping_hits`, range 0, keep higher
    /// score) sorted into output order (start ascending on the `+` search
    /// strand; C `sort_by_tRNAscanSE_output`).
    fn bhb_hits(
        bhb: &[FaithfulSearcher],
        target: &str,
        cfg: &FaithfulConfig,
    ) -> Vec<infernal::FaithfulHit> {
        let mut cand: Vec<infernal::FaithfulHit> = Vec::new();
        for s in bhb {
            cand.extend(s.search(&[target], cfg));
        }
        cand.sort_by(|a, b| {
            a.start
                .cmp(&b.start)
                .then(b.score.partial_cmp(&a.score).unwrap_or(std::cmp::Ordering::Equal))
        });
        let mut i = cand.len() as isize - 2;
        while i >= 0 {
            let ui = i as usize;
            if Self::seg_overlap(cand[ui].start, cand[ui].stop, cand[ui + 1].start, cand[ui + 1].stop, 0) {
                if cand[ui].score >= cand[ui + 1].score {
                    cand.remove(ui + 1);
                } else {
                    cand.remove(ui);
                }
            }
            i -= 1;
        }
        cand
    }

    /// Faithful port of `CM.pm::check_intron_validity` (:1885). Validates one BHB
    /// hit against the mature re-score gate and, on acceptance, mutates `t`
    /// (boundaries / seq / ss / mat_seq / score / model / introns) and the working
    /// flank state (`cur_up`/`cur_seq`/`cur_dn`). Returns `Some((duplicate,
    /// clip_seq, intron_len))` when the intron is valid (C `$ret_value` true;
    /// `$duplicate` set for a re-found existing intron), or `None` when invalid.
    #[allow(clippy::too_many_arguments)]
    fn check_intron_validity(
        &self,
        t: &mut TRna,
        hit: &infernal::FaithfulHit,
        padded_seq: &str,
        padded_full: &str,
        prev_len: i64,
        cur_up: &mut String,
        cur_seq: &mut String,
        cur_dn: &mut String,
        main: &[(&'static str, FaithfulSearcher, bool)],
        main_cfg: &FaithfulConfig,
    ) -> Option<(bool, String, i64)> {
        // ---- Parse the BHB CS line into (pre, intron, post) (C :1895-1916) ----
        let ali = hit.alignment.as_ref()?;
        let caps = BHB_SS_RE.captures(&ali.csline)?;
        let pre_len = caps.get(1).unwrap().as_str().len();
        let intr_len = caps.get(2).unwrap().as_str().len();

        // cm_intron->seq() with U->T applied (not gap-removed yet, C :1900-1906).
        let full: Vec<char> = ali
            .aseq
            .chars()
            .map(|c| match c {
                'U' => 'T',
                'u' => 't',
                o => o,
            })
            .collect();
        if pre_len + intr_len > full.len() {
            return None;
        }
        let degap = |s: &[char]| -> String { s.iter().filter(|&&c| c != '-').collect() };
        let pre_seq = degap(&full[..pre_len]);
        let intron_seq = degap(&full[pre_len..pre_len + intr_len]);
        let post_seq = degap(&full[pre_len + intr_len..]);
        if intron_seq.is_empty() {
            return None;
        }
        let intron_len = intron_seq.len() as i64;

        // ---- Clip the intron out of the current padded seq (C :1931) ----
        let cm_start = hit.start;
        let cm_end = hit.stop;
        let lo = (cm_start - prev_len + pre_seq.len() as i64 - 1).max(0) as usize;
        let hi = (cm_end - prev_len - post_seq.len() as i64).max(0) as usize;
        if lo > padded_seq.len() || hi > padded_seq.len() || lo > hi {
            return None;
        }
        let mut clip = String::with_capacity(padded_seq.len());
        clip.push_str(&padded_seq[..lo]);
        clip.push_str(&padded_seq[hi..]);

        // ---- Re-score the mature candidate with the main CM(s) (C :1946-1964) ----
        let mut mhits: Vec<(&'static str, infernal::FaithfulHit)> = Vec::new();
        for (role, s, _is_sec) in main.iter() {
            for mh in s.search(&[clip.as_str()], main_cfg) {
                mhits.push((*role, mh));
            }
        }
        Self::merge_main_hits(&mut mhits);

        // No mature hit at all: C leaves $ret_value == 1 (intron still "valid",
        // sequence clipped) but changes nothing on the tRNA (:1961 guard false).
        if mhits.is_empty() {
            return Some((false, clip, intron_len));
        }

        let intron_up = intron_seq.to_uppercase();
        for (role, mh) in &mhits {
            let mut duplicate = false;
            let mut ret_value = true;
            let ad = match &mh.alignment {
                Some(a) => a,
                None => continue,
            };
            let (mut cm_ss, mut cm_seq) = format_cmsearch_output(&ad.csline, &ad.aseq, &ad.ncline);
            let mut cm_end_m = mh.stop;
            let cm_start_m = mh.start;

            // CCA 3' trim (C :1971-1984). The clip is coding-oriented (top strand).
            if cm_seq.len() >= 3 && cm_ss.len() >= 4 {
                let last3 = &cm_seq[cm_seq.len() - 3..];
                let last4 = &cm_ss[cm_ss.len() - 4..];
                if !last3.eq_ignore_ascii_case("CCA") && last4 == "...." {
                    cm_end_m -= 3;
                    cm_seq.truncate(cm_seq.len() - 3);
                    cm_ss.truncate(cm_ss.len() - 3);
                }
            }

            let upstream_len = cm_start_m - 1;
            let downstream_len = clip.len() as i64 - cm_end_m;
            // seq = substr(padded_full, cm_start-1, len - up_len - dn_len)
            //     = padded_full[up_len .. len - dn_len]  (C :1989).
            let fl = padded_full.len() as i64;
            let s0 = upstream_len.max(0);
            let s1 = fl - downstream_len;
            let seq_recon = if s0 >= 0 && s1 >= s0 && (s1 as usize) <= padded_full.len() {
                padded_full[s0 as usize..s1 as usize].to_string()
            } else {
                String::new()
            };
            let seq_up = seq_recon.to_uppercase();

            // Locate the intron in the reconstructed precursor (C :1990).
            let mut intron_start = 0i64;
            let mut intron_end = 0i64;
            if seq_recon.is_empty() {
                ret_value = false;
            } else {
                match seq_up.find(&intron_up) {
                    None => ret_value = false,
                    Some(p) => {
                        intron_start = p as i64 + 1;
                        intron_end = intron_len + p as i64;
                        // Duplicate / inclusion / overlap checks against the tRNA's
                        // existing introns, re-indexed in the new seq (C :2000-2058).
                        for iv in t.introns.iter() {
                            let (rs, re) = Self::adjust_rel(&seq_up, iv);
                            if rs == intron_start && re == intron_end {
                                duplicate = true;
                                break;
                            } else if iv.intron_type == "CI"
                                && iv.seq.len() > 40
                                && rs < intron_start
                                && re > intron_start
                                && rs < intron_end
                                && re > intron_end
                            {
                                ret_value = false;
                                break;
                            } else if iv.intron_type == "CI"
                                && rs == intron_start
                                && Self::seg_overlap(rs, re, intron_start, intron_end, 0)
                            {
                                ret_value = false;
                                break;
                            } else if iv.intron_type == "NCI"
                                && Self::seg_overlap(rs, re, intron_start, intron_end, 0)
                            {
                                ret_value = false;
                                break;
                            }
                        }
                    }
                }
            }

            // ---- Genomic boundary reconstruction (C :2061-2075) ----
            let up_len = cur_up.len() as i64;
            let dn_len = cur_dn.len() as i64;
            let (new_start, new_end) = match t.strand {
                TStrand::Minus => {
                    let downstream_start = t.start - dn_len;
                    let upstream_end = t.end + up_len;
                    let downstream_end = downstream_start + downstream_len - 1;
                    let upstream_start = upstream_end - upstream_len + 1;
                    (downstream_end + 1, upstream_start - 1)
                }
                _ => {
                    let upstream_start = t.start - up_len;
                    let downstream_end = t.end + dn_len;
                    let upstream_end = upstream_start + upstream_len - 1;
                    let downstream_start = downstream_end - downstream_len + 1;
                    (upstream_end + 1, downstream_start - 1)
                }
            };

            let hit_overlap = Self::seg_overlap(t.start, t.end, new_start, new_end, 40);
            if ret_value
                && hit_overlap
                && (mh.score as f64) > t.score
                && cm_seq.len() as i64 >= 70
            {
                // Re-index existing introns in the new seq (C :2080-2089), genomic
                // coords recomputed from the OLD boundaries (order preserved).
                for i in 0..t.introns.len() {
                    let (rs_cur, re_cur, itype, iseq) = {
                        let iv = &t.introns[i];
                        (iv.rel_start, iv.rel_end, iv.intron_type.clone(), iv.seq.clone())
                    };
                    let cur = Self::substr_uc(&seq_up, rs_cur as i64, re_cur as i64);
                    if cur != iseq.to_uppercase() {
                        if let Some(p) = seq_up.find(&iseq.to_uppercase()) {
                            let nrs = p as i64 + 1;
                            let nre = iseq.len() as i64 + p as i64;
                            let (gs, ge) = match t.strand {
                                TStrand::Minus => (t.end - nre + 1, t.end - nrs + 1),
                                _ => (t.start + nrs - 1, t.start + nre - 1),
                            };
                            let iv = &mut t.introns[i];
                            iv.rel_start = nrs as i32;
                            iv.rel_end = nre as i32;
                            iv.intron_type = itype;
                            iv.seq = iseq;
                            iv.start = gs;
                            iv.end = ge;
                        }
                    }
                }

                // Update working flanks / precursor (C :2090-2092).
                *cur_up = clip[..(upstream_len.max(0) as usize).min(clip.len())].to_string();
                *cur_dn = clip[(cm_end_m.max(0) as usize).min(clip.len())..].to_string();
                *cur_seq = seq_recon.clone();

                t.seq = seq_recon;
                t.ss = cm_ss.clone();
                t.mat_seq = cm_seq;
                t.mat_ss = cm_ss;
                t.start = new_start;
                t.end = new_end;
                t.score = mh.score as f64;
                t.set_domain_model("infernal", mh.score as f64);
                t.model = role.to_string();

                if !duplicate {
                    let (gs, ge) = match t.strand {
                        TStrand::Minus => (t.end - intron_end + 1, t.end - intron_start + 1),
                        _ => (t.start + intron_start - 1, t.start + intron_end - 1),
                    };
                    t.introns.push(Intron {
                        rel_start: intron_start as i32,
                        rel_end: intron_end as i32,
                        start: gs,
                        end: ge,
                        intron_type: "NCI".to_string(),
                        seq: intron_up.clone(),
                    });
                }
                return Some((duplicate, clip, intron_len));
            }
        }
        None
    }

    /// C `sort_by_tRNAscanSE_output` + `merge_overlapping_hits(0)` for the main
    /// clip re-score hits: sort by start ascending (then score descending) and
    /// drop lower-scoring hits that overlap an adjacent kept hit.
    fn merge_main_hits(hits: &mut Vec<(&'static str, infernal::FaithfulHit)>) {
        hits.sort_by(|a, b| {
            a.1.start
                .cmp(&b.1.start)
                .then(b.1.score.partial_cmp(&a.1.score).unwrap_or(std::cmp::Ordering::Equal))
        });
        let mut i = hits.len() as isize - 2;
        while i >= 0 {
            let ui = i as usize;
            if Self::seg_overlap(hits[ui].1.start, hits[ui].1.stop, hits[ui + 1].1.start, hits[ui + 1].1.stop, 0) {
                if hits[ui].1.score >= hits[ui + 1].1.score {
                    hits.remove(ui + 1);
                } else {
                    hits.remove(ui);
                }
            }
            i -= 1;
        }
    }

    /// Uppercase `seq[rel_start-1 .. rel_end]` (1-based inclusive), clamped.
    fn substr_uc(seq_up: &str, rel_start: i64, rel_end: i64) -> String {
        let a = (rel_start - 1).max(0) as usize;
        let b = rel_end.max(0) as usize;
        if a >= seq_up.len() || b <= a {
            return String::new();
        }
        seq_up[a..b.min(seq_up.len())].to_string()
    }

    /// Re-index an existing intron against the reconstructed (uppercase) precursor
    /// (C :2011-2021): keep its stored rel coords if the substring still matches,
    /// otherwise relocate by `index()`.
    fn adjust_rel(seq_up: &str, iv: &Intron) -> (i64, i64) {
        let cur = Self::substr_uc(seq_up, iv.rel_start as i64, iv.rel_end as i64);
        if cur == iv.seq.to_uppercase() {
            (iv.rel_start as i64, iv.rel_end as i64)
        } else if let Some(p) = seq_up.find(&iv.seq.to_uppercase()) {
            (p as i64 + 1, iv.seq.len() as i64 + p as i64)
        } else {
            (iv.rel_start as i64, iv.rel_end as i64)
        }
    }

    /// Faithful port of `CM.pm::add_canonical_intron` (:1759). Re-inserts a removed
    /// canonical intron sequence into the precursor and, if the mature changes,
    /// re-scores it with the main CM to refresh ss / mat_seq / model. (Rarely
    /// exercised: only when a BHB hit duplicated an existing canonical intron.)
    fn add_canonical_intron(
        &self,
        t: &mut TRna,
        ci_seq: &str,
        main: &[(&'static str, FaithfulSearcher, bool)],
        main_cfg: &FaithfulConfig,
    ) {
        if ci_seq.is_empty() {
            return;
        }
        let precursor = t.seq.clone();
        let up = precursor.to_uppercase();
        let ci_up = ci_seq.to_uppercase();
        let precursor = if let Some(idx) = up.find(&ci_up) {
            format!(
                "{}{}{}",
                &precursor[..idx],
                ci_seq.to_lowercase(),
                &precursor[idx + ci_seq.len()..]
            )
        } else {
            precursor
        };

        // Splice all introns to form the mature sequence.
        let mut introns = t.introns.clone();
        introns.sort_by_key(|i| i.rel_start);
        let mat_seq = Self::splice_introns(&precursor, &introns);

        if !mat_seq.eq_ignore_ascii_case(&t.mat_seq) {
            let mut mhits: Vec<(&'static str, infernal::FaithfulHit)> = Vec::new();
            for (role, s, _is_sec) in main.iter() {
                for mh in s.search(&[mat_seq.as_str()], main_cfg) {
                    mhits.push((*role, mh));
                }
            }
            Self::merge_main_hits(&mut mhits);
            if let Some((role, mh)) = mhits.first() {
                if let Some(ad) = &mh.alignment {
                    let (mut cm_ss, mut cm_seq) =
                        format_cmsearch_output(&ad.csline, &ad.aseq, &ad.ncline);
                    if cm_seq.len() >= 3 && cm_ss.len() >= 4 {
                        let last3 = &cm_seq[cm_seq.len() - 3..];
                        let last4 = &cm_ss[cm_ss.len() - 4..];
                        if !last3.eq_ignore_ascii_case("CCA") && last4 == "...." {
                            cm_seq.truncate(cm_seq.len() - 3);
                            cm_ss.truncate(cm_ss.len() - 3);
                        }
                    }
                    t.ss = cm_ss.clone();
                    t.mat_seq = cm_seq;
                    t.mat_ss = cm_ss;
                    t.score = mh.score as f64;
                    t.set_domain_model("infernal", mh.score as f64);
                    t.model = role.to_string();
                }
            }
        }
    }

    /// Splice the given (rel_start-sorted) introns out of `precursor`
    /// (C `decode_nci_tRNA_properties`/`add_canonical_intron` splice, 1-based rel).
    fn splice_introns(precursor: &str, introns: &[Intron]) -> String {
        if introns.is_empty() {
            return precursor.to_string();
        }
        let n = precursor.len();
        let mut mat = String::with_capacity(n);
        for (i, iv) in introns.iter().enumerate() {
            let (lo, hi) = if i == 0 {
                (0usize, (iv.rel_start as usize).saturating_sub(1))
            } else {
                let prev_end = introns[i - 1].rel_end as usize; // 1-based inclusive
                let start = prev_end; // 0-based index just after prev intron
                let end = (iv.rel_start as usize).saturating_sub(1);
                (start.min(n), end.min(n))
            };
            if lo <= hi && hi <= n {
                mat.push_str(&precursor[lo..hi]);
            }
        }
        let last_end = introns[introns.len() - 1].rel_end as usize; // 1-based inclusive
        if last_end <= n {
            mat.push_str(&precursor[last_end..]);
        }
        mat
    }

    /// Faithful port of `CM.pm::decode_nci_tRNA_properties` (:1069). Splices all
    /// introns to rebuild the mature sequence, re-runs `find_anticodon` on the
    /// clean mature (recovering the true anticodon that the intron-garbled loop
    /// hid), shifts the anticodon index back into precursor coordinates, and
    /// recomputes the isotype / anticodon positions.
    fn decode_nci_trna_properties(&self, t: &mut TRna) {
        let precursor = t.seq.clone();
        let mut introns = t.introns.clone();
        introns.sort_by_key(|i| i.rel_start);
        if introns.is_empty() {
            return;
        }
        let mut mat_seq = Self::splice_introns(&precursor, &introns);

        // If the spliced mature matches the stored mat_seq, decode on it (C :1099).
        if mat_seq.eq_ignore_ascii_case(&t.mat_seq) {
            t.seq = t.mat_seq.clone();
        }

        let (anticodon, antiloop_index, antiloop_end, mut acodon_index) =
            find_anticodon(&t.seq, &t.ss);
        let mut acodon_index2 = 0i32;
        let mut acodon_end1 = 0i32;
        t.anticodon = anticodon.clone();

        // Shift the anticodon index past each intron -> precursor coords (C :1114).
        for iv in introns.iter() {
            let rs = iv.rel_start;
            let re = iv.rel_end;
            if acodon_index >= rs {
                acodon_index += re - rs + 1;
            } else if acodon_index < rs && acodon_index + 2 >= rs {
                acodon_end1 = rs - 1;
                acodon_index2 = re + 1;
            } else if acodon_index + 2 < rs {
                break;
            }
        }
        t.ac_positions.clear();
        if acodon_index2 == 0 {
            if acodon_index > 0 {
                t.ac_positions.push(AnticodonPos {
                    rel_start: acodon_index,
                    rel_end: acodon_index + 2,
                });
            }
        } else {
            t.ac_positions.push(AnticodonPos {
                rel_start: acodon_index,
                rel_end: acodon_end1,
            });
            t.ac_positions.push(AnticodonPos {
                rel_start: acodon_index2,
                rel_end: (3 - (acodon_end1 - acodon_index + 1)) + acodon_index2 - 1,
            });
        }

        if anticodon == UNDEF_ANTICODON || t.seq == "Error" {
            t.anticodon = UNDEF_ANTICODON.to_string();
            t.isotype = UNDEF_ISOTYPE.to_string();
        } else {
            // Canonical anticodon-loop intron in the (now clean) mature (C :1167).
            let (ci_seq, istart, iend) = find_intron(&t.seq, antiloop_index, antiloop_end);
            if !ci_seq.is_empty() {
                let seqs = t.seq.clone();
                let sss = t.ss.clone();
                let is = istart as usize;
                let ie = iend as usize;
                if is >= 1 && ie <= seqs.len() && is <= ie {
                    // Re-splice the mature to drop the anticodon-loop intron too
                    // (C :1171-1172): mat_seq / mat_ss lose [istart..iend].
                    mat_seq = format!("{}{}", &seqs[..is - 1], &seqs[ie..]);
                    t.mat_ss = format!("{}{}", &sss[..is - 1], &sss[ie..]);
                    t.ss = t.mat_ss.clone();
                }
                let mut ci_start = istart;
                let mut ci_end = iend;
                for iv in introns.iter() {
                    if ci_start > iv.rel_start {
                        ci_start += iv.rel_end - iv.rel_start + 1;
                        ci_end += iv.rel_end - iv.rel_start + 1;
                    } else if ci_end < iv.rel_start {
                        break;
                    }
                }
                let (gs, ge) = match t.strand {
                    TStrand::Minus => (t.end - ci_end as i64 + 1, t.end - ci_start as i64 + 1),
                    _ => (t.start + ci_start as i64 - 1, t.start + ci_end as i64 - 1),
                };
                t.introns.push(Intron {
                    rel_start: ci_start,
                    rel_end: ci_end,
                    start: gs,
                    end: ge,
                    intron_type: "CI".to_string(),
                    seq: ci_seq.to_uppercase(),
                });
                Self::merge_introns(t);
                introns = t.introns.clone();
                introns.sort_by_key(|i| i.rel_start);
            }
            t.isotype = get_trna_type(&t.anticodon, &t.model, t.model == "SeC", false);
        }

        // Reset a noncanonical intron to canonical when it abuts the anticodon and
        // no canonical intron exists (C :1206-1236). Affects only the intron type.
        let has_ci = t.introns.iter().any(|i| i.intron_type == "CI");
        if !has_ci {
            for iv in t.introns.iter_mut() {
                let target_rs = if acodon_index2 == 0 {
                    acodon_index + 2 + 2
                } else {
                    (3 - (acodon_end1 - acodon_index + 1)) + acodon_index2 - 1 + 2
                };
                if target_rs == iv.rel_start && iv.intron_type == "NCI" {
                    iv.intron_type = "CI".to_string();
                    break;
                }
            }
        }

        t.mat_seq = mat_seq;
        t.seq = precursor;
    }

    /// `tRNA.pm::merge_introns` (:602): merge introns that abut (prev.rel_end ==
    /// next.rel_start - 1), keeping the earlier start and promoting to CI.
    fn merge_introns(t: &mut TRna) {
        let strand = t.strand;
        t.introns.sort_by_key(|i| i.rel_start);
        let mut merged: Vec<Intron> = Vec::new();
        for iv in t.introns.drain(..) {
            if let Some(last) = merged.last_mut() {
                if last.rel_end == iv.rel_start - 1 {
                    last.rel_end = iv.rel_end;
                    match strand {
                        TStrand::Minus => last.start = iv.start,
                        _ => last.end = iv.end,
                    }
                    if last.intron_type == "CI" || iv.intron_type == "CI" {
                        last.intron_type = "CI".to_string();
                    }
                    // seq is not re-derived here (C keeps the merged rel span only).
                    continue;
                }
            }
            merged.push(iv);
        }
        t.introns = merged;
    }

    /// The 70-nt coding up/downstream flanks of a tRNA (C `upstream()` /
    /// `downstream()`), clamped to the source sequence. On `+` upstream is the
    /// genomic 5' side; on `-` it is the reverse complement of the genomic 3'
    /// side. Returns `(upstream, downstream)` already in coding orientation.
    fn coding_flanks(seq: &[u8], t: &TRna, flank: i64, seqlen: i64) -> (String, String) {
        let (lo, hi) = (t.start, t.end); // genomic low/high (1-based)
        let s = |a: i64, b: i64| -> &[u8] {
            let a = a.max(1);
            let b = b.min(seqlen);
            if a > b {
                &[]
            } else {
                &seq[(a - 1) as usize..b as usize]
            }
        };
        match t.strand {
            TStrand::Minus => {
                // Coding 5' = revcomp of genomic (hi+1 .. hi+flank).
                let up = reverse_complement(s(hi + 1, hi + flank));
                // Coding 3' = revcomp of genomic (lo-flank .. lo-1).
                let dn = reverse_complement(s(lo - flank, lo - 1));
                (
                    String::from_utf8_lossy(&up).to_string(),
                    String::from_utf8_lossy(&dn).to_string(),
                )
            }
            _ => {
                let up = s(lo - flank, lo - 1);
                let dn = s(hi + 1, hi + flank);
                (
                    String::from_utf8_lossy(up).to_string(),
                    String::from_utf8_lossy(dn).to_string(),
                )
            }
        }
    }

    /// Dedup overlapping same-strand hits (from either model / adjacent
    /// candidate regions), keeping the higher-scoring hit.
    fn dedup_faithful_hits(&self, mut hits: Vec<TRna>) -> Vec<TRna> {
        // Higher score first so the greedy keep retains the best of an overlap set.
        hits.sort_by(|a, b| {
            b.score
                .partial_cmp(&a.score)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let mut kept: Vec<TRna> = Vec::new();
        for h in hits {
            let overlaps = kept.iter().any(|k| {
                k.strand == h.strand
                    && std::cmp::max(k.start, h.start) <= std::cmp::min(k.end, h.end)
            });
            if !overlaps {
                kept.push(h);
            }
        }
        kept
    }

    /// Sort `TRna`s in tRNAscan-SE output order (IntResultFile::sort_by_tRNAscanSE_output):
    /// `+` strand first (ascending start), then `-` strand (descending end).
    fn sort_faithful(hits: &mut [TRna]) {
        hits.sort_by(|a, b| {
            let sa = a.strand == TStrand::Plus;
            let sb = b.strand == TStrand::Plus;
            match (sa, sb) {
                (true, false) => std::cmp::Ordering::Less,
                (false, true) => std::cmp::Ordering::Greater,
                (true, true) => a.start.cmp(&b.start),
                (false, false) => b.end.cmp(&a.end),
            }
        });
    }

    /// Format the 9-column `-B` `.out` (default: no `-H`, no `--detail`) using
    /// column widths frozen from the FIRST tRNA (spec 3.1/3.2/3.3).
    pub fn write_faithful_out<W: Write>(&self, writer: &mut W, brief: bool) -> IoResult<()> {
        let hmm = self.get_hmm_score;
        let detail = self.detail;
        let first = match self.trna_results.first() {
            Some(t) => t,
            None => {
                // No tRNAs: still emit the header (unless brief) with default widths.
                if !brief {
                    write_faithful_header(writer, 8, 1, hmm, detail)?;
                }
                return Ok(());
            }
        };
        let w = std::cmp::max(first.seqname.len() + 1, 8);
        let l = format!("{}", first.src_seqlen).len().max(1);

        if !brief {
            write_faithful_header(writer, w, l, hmm, detail)?;
        }

        for t in &self.trna_results {
            // Strand-oriented Begin/End.
            let (begin, end) = match t.strand {
                TStrand::Minus => (t.end, t.start),
                _ => (t.start, t.end),
            };
            // Intron Begin/End (0/0 if none; strand-swapped for `-`). Multiple
            // introns are comma-joined in coding (rel_start) order (C
            // IntResultFile multi-intron output).
            let (ibeg, iend): (String, String) = if t.introns.is_empty() {
                ("0".to_string(), "0".to_string())
            } else {
                let mut ordered: Vec<&Intron> = t.introns.iter().collect();
                ordered.sort_by_key(|i| i.rel_start);
                let (mut bs, mut es): (Vec<String>, Vec<String>) = (Vec::new(), Vec::new());
                for intron in ordered {
                    let (b, e) = match t.strand {
                        TStrand::Minus => (intron.end, intron.start),
                        _ => (intron.start, intron.end),
                    };
                    bs.push(b.to_string());
                    es.push(e.to_string());
                }
                (bs.join(","), es.join(","))
            };
            // Base 9 columns through Inf Score (raw, no printf).
            write!(
                writer,
                "{:<w$}\t{}\t{:<l$}\t{:<l$}\t{}\t{}\t{}\t{}\t{:.1}",
                t.seqname,
                t.id,
                begin,
                end,
                t.isotype,
                t.anticodon,
                ibeg,
                iend,
                t.score,
                w = w,
                l = l,
            )?;
            // HMM Score + 2'Str Score (`-H`).
            if hmm {
                write!(writer, "\t{:.2}\t{:.2}", t.hmm_score, t.ss_score)?;
            }
            // Isotype CM + Isotype Score (`--detail`).
            if detail {
                write!(writer, "\t{}\t{:.1}", t.iso_model, t.iso_score)?;
            }
            // Note (always the trailing column; leading TAB, then note body).
            writeln!(writer, "\t{}", t.note)?;
        }
        Ok(())
    }

    /// `bed_output` sort order (IntResultFile.pm `sort_by_bed_output`): by genomic
    /// low coord (`start`) ascending, then high coord (`end`). Returns references to
    /// `self.trna_results` in that order (used by BED + GFF, which C sorts identically).
    fn bed_output_order(&self) -> Vec<&TRna> {
        let mut order: Vec<&TRna> = self.trna_results.iter().collect();
        order.sort_by(|a, b| a.start.cmp(&b.start).then(a.end.cmp(&b.end)));
        order
    }

    /// `ScanResult.pm::convert_bed_score` (:1101): `score * 10`, clamped to [0, 1000].
    /// The stored score is already `%.1f`-rounded, so `*10` is integral.
    fn convert_bed_score(score: f64) -> i64 {
        let bs = (Self::round1(score) * 10.0).round() as i64;
        bs.clamp(0, 1000)
    }

    /// `ScanResult.pm::save_allStruct_output` (:408) — secondary-structure (`-f`).
    /// Iterates `trna_results` in id order.
    pub fn write_faithful_struct<W: Write>(&self, writer: &mut W) -> IoResult<()> {
        let ruler: String = "    *    |".repeat(20);
        for t in &self.trna_results {
            let seqlen = t.seq.len();
            let plus = t.strand != TStrand::Minus;
            // Header coords: (start-end) for +, (end-start) for -.
            let (dstart, dend) = if plus { (t.start, t.end) } else { (t.end, t.start) };
            write!(
                writer,
                "{}.trna{} ({}-{})\tLength: {} bp\nType: {}\t",
                t.seqname, t.id, dstart, dend, seqlen, t.isotype
            )?;
            write!(writer, "Anticodon: {} at ", t.anticodon)?;
            if t.anticodon == "NNN" || t.ac_positions.is_empty() {
                write!(writer, "0-0 (0-0)\t")?;
            } else {
                for (i, ap) in t.ac_positions.iter().enumerate() {
                    if i > 0 {
                        write!(writer, ",")?;
                    }
                    write!(writer, "{}-{}", ap.rel_start, ap.rel_end)?;
                }
                let n = t.ac_positions.len();
                for (i, ap) in t.ac_positions.iter().enumerate() {
                    if i == 0 {
                        write!(writer, " (")?;
                    } else {
                        write!(writer, ",")?;
                    }
                    let rs = ap.rel_start as i64;
                    if plus {
                        write!(writer, "{}-{}", rs + t.start - 1, rs + t.start + 1)?;
                    } else {
                        write!(writer, "{}-{}", t.end - rs + 1, t.end - rs - 1)?;
                    }
                    if i == n - 1 {
                        write!(writer, ")\t")?;
                    }
                }
            }
            write!(writer, "Score: {:.1}\n", t.score)?;

            for intron in &t.introns {
                write!(
                    writer,
                    "Possible intron: {}-{} ",
                    intron.rel_start, intron.rel_end
                )?;
                if plus {
                    write!(writer, "({}-{})\n", intron.start, intron.end)?;
                } else {
                    write!(writer, "({}-{})\n", intron.end, intron.start)?;
                }
            }

            // Note (default view: no --mito, no infernal_score): pseudo/trunc, then
            // the -H HMM/2'Str block. Truncation is never labelled (the flag-6
            // truncated-CM search is not ported — see the intentional-gap note in
            // faithful_scan_sequence), so only the pseudogene state can appear.
            let base = if t.is_pseudo && t.is_trunc() {
                "Possible truncation, pseudogene"
            } else if t.is_pseudo {
                "Possible pseudogene"
            } else if t.is_trunc() {
                "Possible truncation"
            } else {
                ""
            };
            let mut line = base.to_string();
            if self.get_hmm_score {
                if !base.is_empty() {
                    line.push_str(": ");
                }
                line.push_str(&format!(
                    "HMM Sc={:.2}\tSec struct Sc={:.2}",
                    t.hmm_score, t.ss_score
                ));
            }
            if !line.is_empty() {
                write!(writer, "{}\n", line)?;
            }

            // C `save_allStruct_output` branches on arch_mode (ScanResult.pm:578).
            if self.mode != ScanMode::Archaeal {
                // Bacterial / general / euk: precursor Seq/Str, ruler = len(seq)-1.
                let take = seqlen.saturating_sub(1).min(ruler.len());
                write!(writer, "     {}\n", &ruler[..take])?;
                write!(writer, "Seq: {}\nStr: {}\n\n", t.seq, t.ss)?;
            } else {
                // Archaeal: mature Seq/Str (ruler = len(mat_seq)-1), plus a bracketed
                // precursor `Pre:` line when the tRNA carries any intron (:585-606).
                let take = t.mat_seq.len().saturating_sub(1).min(ruler.len());
                write!(writer, "     {}\n", &ruler[..take])?;
                write!(writer, "Seq: {}\nStr: {}\n", t.mat_seq, t.mat_ss)?;
                if !t.seq.eq_ignore_ascii_case(&t.mat_seq) {
                    let mut precursor = t.seq.to_uppercase();
                    for intron in &t.introns {
                        let rs = intron.rel_start as usize;
                        let re = intron.rel_end as usize;
                        if rs >= 1 && re <= t.seq.len() && rs <= re {
                            let intron_seq = t.seq[rs - 1..re].to_uppercase();
                            if !intron_seq.is_empty() {
                                if let Some(pos) = precursor.find(&intron_seq) {
                                    precursor.replace_range(
                                        pos..pos + intron_seq.len(),
                                        &format!("[{}]", intron_seq),
                                    );
                                }
                            }
                        }
                    }
                    write!(writer, "Pre: {}\n\n", precursor)?;
                } else {
                    write!(writer, "\n")?;
                }
            }
        }
        Ok(())
    }

    /// `ScanResult.pm::write_bed` (:997) — 12-column BED (`-b`). Sorted by genomic
    /// position (`bed_output`).
    pub fn write_faithful_bed<W: Write>(&self, writer: &mut W) -> IoResult<()> {
        for t in self.bed_output_order() {
            let name = format!("{}.tRNA{}-{}{}", t.seqname, t.id, t.isotype, t.anticodon);
            let strand = t.strand.as_str();
            write!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t0\t{}\t",
                t.seqname,
                t.start - 1,
                t.end,
                name,
                Self::convert_bed_score(t.score),
                strand,
                t.start - 1,
                t.end,
                t.introns.len() + 1,
            )?;
            if t.introns.is_empty() {
                write!(writer, "{},\t0,\n", t.end - t.start + 1)?;
            } else {
                let mut block_sizes = String::new();
                let mut block_starts = String::from("0,");
                if strand == "+" {
                    let mut prev_start: i64 = 1;
                    for intron in &t.introns {
                        block_sizes
                            .push_str(&format!("{},", intron.rel_start as i64 - prev_start));
                        block_starts.push_str(&format!("{},", intron.rel_end));
                        prev_start = intron.rel_end as i64 + 1;
                    }
                    let last = t.introns.last().unwrap();
                    block_sizes.push_str(&format!("{},", t.end - last.end));
                } else {
                    let mut prev_start: i64 = t.seq.len() as i64;
                    for intron in t.introns.iter().rev() {
                        block_sizes
                            .push_str(&format!("{},", prev_start - intron.rel_end as i64));
                        block_starts
                            .push_str(&format!("{},", prev_start - intron.rel_start as i64 + 1));
                        prev_start = intron.rel_start as i64;
                    }
                    let first = t.introns.first().unwrap();
                    block_sizes.push_str(&format!("{},", first.rel_start as i64 - 1));
                }
                write!(writer, "{}\t{}\n", block_sizes, block_starts)?;
            }
        }
        Ok(())
    }

    /// `ScanResult.pm::write_gff` (:1118) — GFF3 (`--gff`). Sorted by genomic
    /// position (`bed_output`); emits a `tRNA` line + one `exon` line per exon.
    pub fn write_faithful_gff<W: Write>(&self, writer: &mut W) -> IoResult<()> {
        write!(writer, "##gff-version 3\n")?;
        for t in self.bed_output_order() {
            let biotype = if t.is_pseudo { "pseudogene" } else { "tRNA" };
            let name = format!("{}.tRNA{}-{}{}", t.seqname, t.id, t.isotype, t.anticodon);
            let strand = t.strand.as_str();
            write!(
                writer,
                "{}\ttRNAscan-SE\t{}\t{}\t{}\t{:.1}\t{}\t.\tID={}.trna{};Name={};isotype={};anticodon={};gene_biotype={};\n",
                t.seqname, biotype, t.start, t.end, t.score, strand,
                t.seqname, t.id, name, t.isotype, t.anticodon, biotype,
            )?;
            if t.introns.is_empty() {
                write!(
                    writer,
                    "{}\ttRNAscan-SE\texon\t{}\t{}\t.\t{}\t.\tID={}.trna{}.exon1;Parent={}.trna{};\n",
                    t.seqname, t.start, t.end, strand, t.seqname, t.id, t.seqname, t.id,
                )?;
            } else if strand == "+" {
                write!(writer, "{}\ttRNAscan-SE\texon\t{}\t", t.seqname, t.start)?;
                for (i, intron) in t.introns.iter().enumerate() {
                    write!(
                        writer,
                        "{}\t.\t{}\t.\tID={}.trna{}.exon{};Parent={}.trna{};\n",
                        intron.start - 1, strand, t.seqname, t.id, i + 1, t.seqname, t.id,
                    )?;
                    write!(
                        writer,
                        "{}\ttRNAscan-SE\texon\t{}\t",
                        t.seqname,
                        intron.end + 1
                    )?;
                }
                write!(
                    writer,
                    "{}\t.\t{}\t.\tID={}.trna{}.exon{};Parent={}.trna{};\n",
                    t.end, strand, t.seqname, t.id, t.introns.len() + 1, t.seqname, t.id,
                )?;
            } else {
                let mut end = t.end;
                for (i, intron) in t.introns.iter().enumerate() {
                    write!(
                        writer,
                        "{}\ttRNAscan-SE\texon\t{}\t{}\t.\t{}\t.\tID={}.trna{}.exon{};Parent={}.trna{};\n",
                        t.seqname, intron.end + 1, end, strand, t.seqname, t.id, i + 1, t.seqname, t.id,
                    )?;
                    end = intron.start - 1;
                }
                write!(
                    writer,
                    "{}\ttRNAscan-SE\texon\t{}\t{}\t.\t{}\t.\tID={}.trna{}.exon{};Parent={}.trna{};\n",
                    t.seqname, t.start, end, strand, t.seqname, t.id, t.introns.len() + 1, t.seqname, t.id,
                )?;
            }
        }
        Ok(())
    }

    /// `ScanResult.pm::write_tRNA_sequence` (:610) — FASTA (`-a`). Iterates
    /// `trna_results` in id order; sequence uppercased, wrapped at 60 columns.
    pub fn write_faithful_fasta<W: Write>(&self, writer: &mut W) -> IoResult<()> {
        for t in &self.trna_results {
            let seqlen = t.seq.len();
            write!(
                writer,
                ">{}.trna{} {}:{}-{} ({}) {} ({}) {} bp Sc: {:.1}",
                t.seqname,
                t.id,
                t.seqname,
                t.start,
                t.end,
                t.strand.as_str(),
                t.isotype,
                t.anticodon,
                seqlen,
                t.score,
            )?;
            if t.is_pseudo {
                write!(writer, " Possible pseudogene\n")?;
            } else {
                write!(writer, "\n")?;
            }
            let up = t.seq.to_uppercase();
            let bytes = up.as_bytes();
            let parts = seqlen / 60;
            let remain = seqlen % 60;
            for j in 0..parts {
                writer.write_all(&bytes[j * 60..j * 60 + 60])?;
                write!(writer, "\n")?;
            }
            if remain > 0 {
                writer.write_all(&bytes[parts * 60..])?;
                write!(writer, "\n")?;
            }
        }
        Ok(())
    }

    /// `Stats.pm::output_summary` (:326) — the deterministic tail of the `-m`
    /// statistics file: the summary counts + the Isotype / Anticodon Counts table.
    /// Driven by `trna_results` (equivalent to C's tab-result parse). The banner /
    /// first-pass / second-pass blocks above it carry timestamps + CPU times and are
    /// inherently non-reproducible; this method emits only the deterministic block.
    pub fn write_faithful_stats<W: Write>(&self, writer: &mut W) -> IoResult<()> {
        // Isotype order + anticodon layout — GeneticCode.pm initialize (:38/:46).
        #[allow(dead_code)]
        const ISOTYPES: [&str; 22] = [
            "Ala", "Gly", "Pro", "Thr", "Val", "Ser", "Arg", "Leu", "Phe", "Asn",
            "Lys", "Asp", "Glu", "His", "Gln", "Ile", "Met", "Tyr", "Supres", "Cys",
            "Trp", "SelCys",
        ];
        const AC_LIST: [(&str, &[&str]); 22] = [
            ("Ala", &["AGC", "GGC", "CGC", "TGC"]),
            ("Gly", &["ACC", "GCC", "CCC", "TCC"]),
            ("Pro", &["AGG", "GGG", "CGG", "TGG"]),
            ("Thr", &["AGT", "GGT", "CGT", "TGT"]),
            ("Val", &["AAC", "GAC", "CAC", "TAC"]),
            ("Ser", &["AGA", "GGA", "CGA", "TGA", "ACT", "GCT"]),
            ("Arg", &["ACG", "GCG", "CCG", "TCG", "CCT", "TCT"]),
            ("Leu", &["AAG", "GAG", "CAG", "TAG", "CAA", "TAA"]),
            ("Phe", &["AAA", "GAA", "&nbsp", "&nbsp"]),
            ("Asn", &["ATT", "GTT", "&nbsp", "&nbsp"]),
            ("Lys", &["&nbsp", "&nbsp", "CTT", "TTT"]),
            ("Asp", &["ATC", "GTC", "&nbsp", "&nbsp"]),
            ("Glu", &["&nbsp", "&nbsp", "CTC", "TTC"]),
            ("His", &["ATG", "GTG", "&nbsp", "&nbsp"]),
            ("Gln", &["&nbsp", "&nbsp", "CTG", "TTG"]),
            ("Ile", &["AAT", "GAT", "CAT", "TAT"]),
            ("Met", &["&nbsp", "&nbsp", "CAT", "&nbsp"]),
            ("Tyr", &["ATA", "GTA", "&nbsp", "&nbsp"]),
            ("Supres", &["&nbsp", "CTA", "TTA", "TCA"]),
            ("Cys", &["ACA", "GCA", "&nbsp", "&nbsp"]),
            ("Trp", &["&nbsp", "&nbsp", "CCA", "&nbsp"]),
            ("SelCys", &["&nbsp", "&nbsp", "&nbsp", "TCA"]),
        ];

        use std::collections::HashMap;
        let mut iso_ar: HashMap<String, usize> = HashMap::new();
        let mut ac_ar: HashMap<String, HashMap<String, usize>> = HashMap::new();
        let mut intron_ac: HashMap<String, HashMap<String, usize>> = HashMap::new();
        let (mut trna_ct, mut selcys_ct, mut pseudo_ct) = (0usize, 0usize, 0usize);
        let (mut mismatch_ct, mut undet_ct, mut stop_sup_ct, mut intron_ct) =
            (0usize, 0usize, 0usize, 0usize);

        for t in &self.trna_results {
            let iso = &t.isotype;
            let ac = &t.anticodon;
            let is_sec = iso.contains("SeC");
            if t.note.contains("pseudo") {
                pseudo_ct += 1;
                *iso_ar.entry("Pseudo".to_string()).or_insert(0) += 1;
            } else if t.note.contains("IPD") {
                mismatch_ct += 1;
            } else if iso == "Undet" {
                undet_ct += 1;
            } else if is_sec {
                selcys_ct += 1;
            } else if iso == "Sup" {
                stop_sup_ct += 1;
            } else {
                trna_ct += 1;
            }

            let key = if is_sec {
                "SelCys".to_string()
            } else if iso == "Sup" {
                "Supres".to_string()
            } else {
                iso.clone()
            };
            *iso_ar.entry(key.clone()).or_insert(0) += 1;
            *ac_ar
                .entry(key.clone())
                .or_default()
                .entry(ac.clone())
                .or_insert(0) += 1;

            let icount = t.introns.len();
            if icount > 0 {
                intron_ct += icount;
                *intron_ac
                    .entry(key)
                    .or_default()
                    .entry(ac.clone())
                    .or_insert(0) += icount;
            }
        }

        let total = trna_ct + selcys_ct + pseudo_ct + mismatch_ct + undet_ct + stop_sup_ct;
        // show_mismatch: (bact|arch|euk) and !no_isotype and --detail (Stats.pm:457).
        let show_mismatch = self.detail
            && !self.no_isotype
            && matches!(
                self.mode,
                ScanMode::Bacterial | ScanMode::Archaeal | ScanMode::Eukaryotic
            );

        write!(writer, "\n")?;
        write!(writer, "tRNAs decoding Standard 20 AA:              {}\n", trna_ct)?;
        write!(writer, "Selenocysteine tRNAs (TCA):                 {}\n", selcys_ct)?;
        write!(writer, "Possible suppressor tRNAs (CTA,TTA,TCA):    {}\n", stop_sup_ct)?;
        write!(writer, "tRNAs with undetermined/unknown isotypes:   {}\n", undet_ct)?;
        if show_mismatch {
            write!(writer, "tRNAs with mismatch isotypes:               {}\n", mismatch_ct)?;
        }
        write!(writer, "Predicted pseudogenes:                      {}\n", pseudo_ct)?;
        write!(writer, "                                            -------\n")?;
        write!(writer, "Total tRNAs:                                {}\n\n", total)?;
        write!(writer, "tRNAs with introns:     \t{}\n\n", intron_ct)?;

        // Intron prefixes, then the closing '|'. C's `defined($intron_ac{aa}{ac})`
        // test autovivifies a key for every isotype, so the `keys > 0` guard is
        // always true → the '|' line is emitted unconditionally.
        for (aa, acs) in AC_LIST.iter() {
            for ac in acs.iter() {
                if let Some(n) = intron_ac.get(*aa).and_then(|m| m.get(*ac)) {
                    write!(writer, "| {}-{}: {} ", aa, ac, n)?;
                }
            }
        }
        write!(writer, "|\n\n")?;

        write!(writer, "Isotype / Anticodon Counts:\n\n")?;
        for (aa, acs) in AC_LIST.iter() {
            let mut label = aa.to_string();
            let mut iso_count = *iso_ar.get(*aa).unwrap_or(&0);
            let iso_cm_count = 0usize; // detail-only isotype-CM tally; 0 in default view.
            if *aa == "Met" {
                if iso_ar.contains_key("iMet") {
                    label = "Met/iMet".to_string();
                    iso_count += iso_ar.get("iMet").copied().unwrap_or(0);
                } else if iso_ar.contains_key("fMet") {
                    label = "Met/fMet".to_string();
                    iso_count += iso_ar.get("fMet").copied().unwrap_or(0);
                }
            } else if *aa == "Ile" {
                iso_count += iso_ar.get("Ile2").copied().unwrap_or(0);
            }

            if show_mismatch {
                write!(writer, "{:<8}: {} ({})\t", label, iso_count, iso_cm_count)?;
            } else {
                write!(writer, "{:<8}: {}\t", label, iso_count)?;
            }

            for ac in acs.iter() {
                if *ac == "&nbsp" {
                    write!(writer, "             ")?;
                } else {
                    let mut count = ac_ar.get(*aa).and_then(|m| m.get(*ac)).copied().unwrap_or(0);
                    if *aa == "Met" {
                        if let Some(m) = ac_ar.get("iMet") {
                            count += m.get(*ac).copied().unwrap_or(0);
                        } else if let Some(m) = ac_ar.get("fMet") {
                            count += m.get(*ac).copied().unwrap_or(0);
                        }
                    } else if *aa == "Ile" {
                        if let Some(m) = ac_ar.get("Ile2") {
                            count += m.get(*ac).copied().unwrap_or(0);
                        }
                    }
                    let cs = if count > 0 { count.to_string() } else { String::new() };
                    write!(writer, "{:>5}: {:<6}", ac, cs)?;
                }
            }
            write!(writer, "\n")?;
        }
        write!(writer, "\n")?;
        Ok(())
    }

    /// Scan a sequence for tRNAs
    pub fn scan_sequence(&mut self, seq: &[u8], sqinfo: &SqInfo) -> Result<(), String> {
        if self.verbose && !self.quiet {
            eprintln!("Scanning {} ({} bp) with {:?} mode", sqinfo.name, seq.len(), self.mode);
        }

        // Faithful in-process Infernal pipeline for bacterial/archaeal/general modes.
        if self.uses_faithful() {
            if seq.len() >= 60 {
                let trnas = self.faithful_scan_sequence(seq, &sqinfo.name, seq.len());
                self.trna_results.extend(trnas);
            }
            return Ok(());
        }

        // Skip very short sequences
        if seq.len() < 60 {
            return Ok(());
        }

        let mut new_results;

        // Check if using Infernal first-pass (scans both strands automatically)
        let uses_infernal_fp = matches!(self.mode, ScanMode::Bacterial | ScanMode::Archaeal | ScanMode::General);

        if uses_infernal_fp {
            // Infernal first-pass scans both strands, no need for reverse complement
            let candidates = self.first_pass_scan(seq, &sqinfo.name);

            if self.verbose && !self.quiet {
                eprintln!("Found {} candidates in first pass (Infernal)", candidates.len());
            }

            // Second-pass: CM verification
            new_results = self.second_pass_scan(&candidates, seq, &sqinfo.name);
        } else {
            // EuFindtRNA-based first-pass - scan both strands separately
            let candidates = self.first_pass_scan(seq, &sqinfo.name);

            if self.verbose && !self.quiet {
                eprintln!("Found {} candidates in first pass (forward)", candidates.len());
            }

            // Second-pass: CM verification
            new_results = self.second_pass_scan(&candidates, seq, &sqinfo.name);

            // Also scan reverse complement for non-Infernal first-pass
            let rc_seq = reverse_complement(seq);
            let rc_candidates = self.first_pass_scan(&rc_seq, &sqinfo.name);

            if self.verbose && !self.quiet {
                eprintln!("Found {} candidates in first pass (reverse)", rc_candidates.len());
            }

            let mut rc_results = self.second_pass_scan(&rc_candidates, &rc_seq, &sqinfo.name);

            // Adjust coordinates for reverse strand
            let seqlen = seq.len() as i64;
            for result in &mut rc_results {
                let tmp_begin = result.begin;
                result.begin = seqlen - result.end + 1;
                result.end = seqlen - tmp_begin + 1;
                result.strand = '-';
            }

            new_results.extend(rc_results);
        }

        // Sort by position
        new_results.sort_by_key(|r| r.begin);

        // Remove duplicates (overlapping hits on same strand)
        new_results = self.deduplicate_results(new_results);

        // Determine isotype and anticodon using cmscan with isotype-specific CMs
        self.determine_isotypes(&mut new_results, seq, &sqinfo.name);

        // Assign tRNA numbers
        let base_num = self.results.iter()
            .filter(|r| r.seq_name == sqinfo.name)
            .count();
        for (i, result) in new_results.iter_mut().enumerate() {
            result.trna_num = base_num + i + 1;
        }

        self.results.extend(new_results);

        Ok(())
    }

    /// Remove duplicate/overlapping results
    fn deduplicate_results(&self, mut results: Vec<TrnaResult>) -> Vec<TrnaResult> {
        if results.is_empty() {
            return results;
        }

        // First, normalize coordinates (ensure begin < end for all)
        for result in &mut results {
            if result.begin > result.end {
                std::mem::swap(&mut result.begin, &mut result.end);
                // Flip strand if we had to swap
                result.strand = if result.strand == '+' { '-' } else { '+' };
            }
        }

        // Sort by position (smallest begin first), then by score (highest first)
        results.sort_by(|a, b| {
            a.begin.cmp(&b.begin)
                .then_with(|| b.inf_score.partial_cmp(&a.inf_score).unwrap_or(std::cmp::Ordering::Equal))
        });

        let mut deduped: Vec<TrnaResult> = Vec::new();

        for result in results {
            // Check if this result overlaps significantly with any existing result
            let overlap_threshold = 20; // Minimum distance to consider as separate

            let overlaps = deduped.iter().any(|existing| {
                let overlap_start = result.begin.max(existing.begin);
                let overlap_end = result.end.min(existing.end);
                let overlap_len = overlap_end - overlap_start;
                overlap_len > overlap_threshold
            });

            if !overlaps {
                deduped.push(result);
            }
        }

        // Re-sort by position for output
        deduped.sort_by_key(|r| r.begin);

        deduped
    }

    /// Write main results table
    pub fn write_results<W: Write>(&self, writer: &mut W) -> IoResult<()> {
        // Results
        for result in &self.results {
            writeln!(writer, "{}", result.format_output_line())?;
        }

        Ok(())
    }

    /// Write secondary structure output
    pub fn write_secondary_structures<W: Write>(&self, writer: &mut W) -> IoResult<()> {
        for result in &self.results {
            write!(writer, "{}", result.format_ss_output())?;
        }
        Ok(())
    }

    /// Write statistics summary
    pub fn write_statistics<W: Write>(&self, writer: &mut W) -> IoResult<()> {
        writeln!(writer, "tRNAscan-SE Statistics")?;
        writeln!(writer, "Total tRNAs: {}", self.results.len())?;

        // Count by isotype
        let mut isotype_counts = std::collections::HashMap::new();
        for result in &self.results {
            *isotype_counts.entry(result.isotype.clone()).or_insert(0) += 1;
        }

        writeln!(writer, "\nIsotype distribution:")?;
        let mut isotypes: Vec<_> = isotype_counts.iter().collect();
        isotypes.sort_by_key(|(k, _)| k.as_str());
        for (isotype, count) in isotypes {
            writeln!(writer, "  {}: {}", isotype, count)?;
        }

        Ok(())
    }

    /// Write BED format output
    pub fn write_bed_format<W: Write>(&self, writer: &mut W) -> IoResult<()> {
        for result in &self.results {
            writeln!(writer, "{}", result.format_bed_line(&result.seq_name))?;
        }
        Ok(())
    }

    /// Write isotype model scores
    pub fn write_isotype_models<W: Write>(&self, writer: &mut W) -> IoResult<()> {
        writeln!(writer, "Sequence\ttRNA#\tIsotype\tCM Score")?;
        for result in &self.results {
            writeln!(writer, "{}\t{}\t{}\t{:.1}",
                result.seq_name, result.trna_num, result.cm_isotype, result.cm_score)?;
        }
        Ok(())
    }
}

/// Write the 9-column `-B` `.out` header (default: no `-H`, no `--detail`).
///
/// Column widths: `w` = max(len(first tRNA seqname)+1, 8); `l` = digits in the
/// first tRNA's source-sequence length (spec 3.1/3.2).
fn write_faithful_header<W: Write>(
    writer: &mut W,
    w: usize,
    l: usize,
    hmm: bool,
    detail: bool,
) -> IoResult<()> {
    // Line 1: base through "Inf", then conditional blocks, then the Note spacer.
    write!(
        writer,
        "{:<w$}\t\t{:<l$}\t{:<l$}\ttRNA\tAnti\tIntron Bounds\tInf",
        "Sequence", "tRNA", "Bounds", w = w, l = l
    )?;
    if hmm {
        write!(writer, "\tHMM\t2'Str")?;
    }
    if detail {
        write!(writer, "\tIsotype\tIsotype")?;
    }
    writeln!(writer, "\t      ")?;

    // Line 2: base through "Score", then conditional blocks, then "Note".
    write!(
        writer,
        "{:<w$}\ttRNA #\t{:<l$}\t{:<l$}\tType\tCodon\tBegin\tEnd\tScore",
        "Name", "Begin", "End", w = w, l = l
    )?;
    if hmm {
        write!(writer, "\tScore\tScore")?;
    }
    if detail {
        write!(writer, "\tCM\tScore")?;
    }
    writeln!(writer, "\tNote")?;

    // Line 3 (separator): Begin underline = 5 dashes, End underline = 6 dashes,
    // Inf Score underline = 6 dashes; conditional HMM/2'Str = 5, Isotype = 7.
    write!(
        writer,
        "{:<w$}\t------\t{:<l$}\t{:<l$}\t----\t-----\t-----\t----\t------",
        "--------", "-----", "------", w = w, l = l
    )?;
    if hmm {
        write!(writer, "\t-----\t-----")?;
    }
    if detail {
        write!(writer, "\t-------\t-------")?;
    }
    writeln!(writer, "\t------")?;
    Ok(())
}

#[cfg(test)]
mod fix_isotype_tests {
    //! Unit tests for the ported CM.pm `fix_fMet` / `fix_His` 5' boundary edits.
    //! These exercise the pure geometric transforms (no rescore) against hand-
    //! worked C `substr` expansions, since neither fix triggered on the available
    //! real genomes (E. coli MG1655, S. acidocaldarius). See report notes.
    use super::*;
    use crate::trna::{Strand, TRna};

    fn met(start: i64, end: i64, plus: bool) -> TRna {
        let mut t = TRna::new();
        t.isotype = "Met".to_string();
        t.score = 45.0; // > 40
        t.strand = if plus { Strand::Plus } else { Strand::Minus };
        t.start = start;
        t.end = end;
        t.src_seqlen = 100000;
        t
    }

    // ---- fix_fMet Branch 1: prepend upstream base, extend 5' ----
    #[test]
    fn fmet_branch1_plus_prepend() {
        // seq[-3:]=="CCA", ss[-5:]==".....", ss[0]!='.' -> prepend upstream base.
        let mut t = met(100, 114, true);
        let mut seq = b"GCACGGATGGCCCCA".to_vec(); // len 15, ends CCA
        let mut ss = b">>>>>>>>>>.....".to_vec(); // ss[0]='>', ss[-5:]='.....'
        // genome[start-2] = genome[98] = 'A'
        let mut genome = vec![b'T'; 200];
        genome[98] = b'A';
        let fired = TrnaScanner::fix_fmet_transform(&mut t, &mut seq, &mut ss, &genome);
        assert!(fired);
        assert_eq!(seq, b"AGCACGGATGGCCCCA".to_vec());
        assert_eq!(ss, b".>>>>>>>>>>.....".to_vec());
        assert_eq!(t.start, 99);
        assert_eq!(t.end, 114);
    }

    #[test]
    fn fmet_branch1_minus_prepend() {
        // Minus strand: upstream coding base = revcomp(genome[end]); end += 1.
        let mut t = met(100, 114, false);
        let mut seq = b"GCACGGATGGCCCCA".to_vec();
        let mut ss = b">>>>>>>>>>.....".to_vec();
        // genome[end] = genome[114] = 'T' -> revcomp -> 'A'
        let mut genome = vec![b'G'; 200];
        genome[114] = b'T';
        let fired = TrnaScanner::fix_fmet_transform(&mut t, &mut seq, &mut ss, &genome);
        assert!(fired);
        assert_eq!(seq, b"AGCACGGATGGCCCCA".to_vec());
        assert_eq!(t.start, 100);
        assert_eq!(t.end, 115);
    }

    // ---- fix_fMet Branch 2: remove 5' bulge base, trim 5' ----
    #[test]
    fn fmet_branch2_stem5_remove() {
        // ss[0]=='.', ss[0..4]=='.>.>', seq[0..2]=='CG', stem5 pattern '<<<..'.
        // pos71 = seq[len(ss)-3] = seq[12] = 'C'; revcomp('C')='G' == uc(seq[2])='G'.
        let mut t = met(200, 214, true);
        let mut seq = b"CGGAAAAAAAAACTT".to_vec(); // len 15
        let mut ss = b".>.>......<<<..".to_vec(); // len 15, ss[-5:]='<<<..'
        let genome = vec![b'T'; 300];
        let fired = TrnaScanner::fix_fmet_transform(&mut t, &mut seq, &mut ss, &genome);
        assert!(fired);
        // seq[1] + uc(seq[2]) + seq[3..] = 'G'+'G'+"AAAAAAAAACTT"
        assert_eq!(seq, b"GGAAAAAAAAACTT".to_vec());
        // ss[0..2] + ss[3..] = ".>" + ">......<<<.."
        assert_eq!(ss, b".>>......<<<..".to_vec());
        assert_eq!(t.start, 201);
        assert_eq!(t.end, 214);
    }

    #[test]
    fn fmet_branch2_no_fire_when_not_partner() {
        // Same as above but seq[2] mismatches revcomp(pos71) -> no edit.
        let mut t = met(200, 214, true);
        let mut seq = b"CGAAAAAAAAAACTT".to_vec(); // seq[2]='A', revcomp('C')='G' != 'A'
        let mut ss = b".>.>......<<<..".to_vec();
        let genome = vec![b'T'; 300];
        let fired = TrnaScanner::fix_fmet_transform(&mut t, &mut seq, &mut ss, &genome);
        assert!(!fired);
        assert_eq!(t.start, 200);
        assert_eq!(t.end, 214);
    }

    #[test]
    fn fmet_no_fire_wrong_isotype_or_score() {
        let mut t = met(100, 114, true);
        t.isotype = "Ala".to_string();
        let mut seq = b"GCACGGATGGCCCCA".to_vec();
        let mut ss = b">>>>>>>>>>.....".to_vec();
        let genome = vec![b'A'; 200];
        assert!(!TrnaScanner::fix_fmet_transform(&mut t, &mut seq, &mut ss, &genome));
        // score gate
        let mut t2 = met(100, 114, true);
        t2.score = 40.0; // not > 40
        let mut seq2 = b"GCACGGATGGCCCCA".to_vec();
        let mut ss2 = b">>>>>>>>>>.....".to_vec();
        assert!(!TrnaScanner::fix_fmet_transform(&mut t2, &mut seq2, &mut ss2, &genome));
    }

    // ---- fix_His: remove extra 5' G bulge, shift start+1/end-1 ----
    #[test]
    fn his_remove_bulge() {
        let mut t = TRna::new();
        t.isotype = "His".to_string();
        t.score = 40.0; // > 35
        t.strand = Strand::Plus;
        t.start = 500;
        t.end = 521;
        let mut seq = b"GATCGAAAAAAAAAAACTTTTT".to_vec(); // len 22, pos5=G, pos68=C
        let mut ss = b">>>>.>>>.....<<<.<<<<.".to_vec(); // len 22
        let fired = TrnaScanner::fix_his_transform(&mut t, &mut seq, &mut ss);
        assert!(fired);
        assert_eq!(seq, b"ATCGAAAAAAAAAAACTTTT".to_vec()); // len 20
        assert_eq!(ss, b">>>>>>>.....<<<<<<<.".to_vec()); // len 20
        assert_eq!(t.start, 501);
        assert_eq!(t.end, 520);
    }

    #[test]
    fn his_no_fire_bad_pair() {
        let mut t = TRna::new();
        t.isotype = "His".to_string();
        t.score = 40.0;
        t.strand = Strand::Plus;
        t.start = 500;
        t.end = 521;
        // pos5=A (idx4), pos68=A (idx16) -> A:A not a valid pair -> no fire.
        let mut seq = b"GATCAAAAAAAAAAAAATTTTT".to_vec();
        let mut ss = b">>>>.>>>.....<<<.<<<<.".to_vec();
        assert!(!TrnaScanner::fix_his_transform(&mut t, &mut seq, &mut ss));
        assert_eq!(t.start, 500);
        assert_eq!(t.end, 521);
    }
}
