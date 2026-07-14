//! Options module for tRNAscan-SE configuration.
//!
//! This module handles all command-line options and configuration
//! for the tRNAscan-SE pipeline, ported from the Perl Options.pm.

use std::fmt;
use std::io::Write;
use std::path::{Path, PathBuf};

/// Error type for options parsing and validation.
#[derive(Debug, Clone)]
pub enum OptionsError {
    /// Invalid search mode combination
    InvalidSearchMode(String),
    /// Invalid organism mode
    InvalidOrgMode(String),
    /// Invalid genetic code
    InvalidGeneticCode(i32),
    /// Missing required file
    MissingFile(String),
    /// Invalid file path
    InvalidPath(PathBuf),
    /// Conflicting options
    ConflictingOptions(String),
    /// Invalid value
    InvalidValue(String),
}

impl fmt::Display for OptionsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            OptionsError::InvalidSearchMode(msg) => write!(f, "Invalid search mode: {}", msg),
            OptionsError::InvalidOrgMode(msg) => write!(f, "Invalid organism mode: {}", msg),
            OptionsError::InvalidGeneticCode(code) => write!(f, "Invalid genetic code: {}", code),
            OptionsError::MissingFile(path) => write!(f, "Missing file: {}", path),
            OptionsError::InvalidPath(path) => write!(f, "Invalid path: {}", path.display()),
            OptionsError::ConflictingOptions(msg) => write!(f, "Conflicting options: {}", msg),
            OptionsError::InvalidValue(msg) => write!(f, "Invalid value: {}", msg),
        }
    }
}

impl std::error::Error for OptionsError {}

/// Search mode for the tRNA scanning pipeline.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SearchMode {
    /// Eukaryotic mode (default)
    #[default]
    Eukaryotic,
    /// Bacterial mode (-B)
    Bacterial,
    /// Archaeal mode (-A)
    Archaeal,
    /// Organellar mode (-O)
    Organellar,
    /// General mode (-G) - combines tRNAs from all 3 domains
    General,
    /// Mitochondrial mode (-M)
    Mitochondrial,
    /// Metagenome mode
    Metagenome,
    /// Nuclear-encoded mitochondrial (numt) mode
    Numt,
    /// Alternate mode
    Alternate,
}

impl SearchMode {
    /// Returns the mode string used internally.
    pub fn as_str(&self) -> &'static str {
        match self {
            SearchMode::Eukaryotic => "euk",
            SearchMode::Bacterial => "bacteria",
            SearchMode::Archaeal => "archaea",
            SearchMode::Organellar => "organelle",
            SearchMode::General => "general",
            SearchMode::Mitochondrial => "mito",
            SearchMode::Metagenome => "metagenome",
            SearchMode::Numt => "numt",
            SearchMode::Alternate => "alt",
        }
    }

    /// Parse from string.
    pub fn from_str(s: &str) -> Option<Self> {
        match s.to_lowercase().as_str() {
            "euk" | "eukaryotic" | "" => Some(SearchMode::Eukaryotic),
            "bacteria" | "bacterial" => Some(SearchMode::Bacterial),
            "archaea" | "archaeal" => Some(SearchMode::Archaeal),
            "organelle" | "organellar" => Some(SearchMode::Organellar),
            "general" => Some(SearchMode::General),
            "mito" | "mitochondrial" => Some(SearchMode::Mitochondrial),
            "metagenome" => Some(SearchMode::Metagenome),
            "numt" => Some(SearchMode::Numt),
            "alt" | "alternate" => Some(SearchMode::Alternate),
            _ => None,
        }
    }

    /// Returns display name for the mode.
    pub fn display_name(&self) -> &'static str {
        match self {
            SearchMode::Eukaryotic => "Eukaryotic",
            SearchMode::Bacterial => "Bacterial",
            SearchMode::Archaeal => "Archaeal",
            SearchMode::Organellar => "Organellar",
            SearchMode::General => "General",
            SearchMode::Mitochondrial => "Mitochondrial",
            SearchMode::Metagenome => "Metagenome",
            SearchMode::Numt => "Numt",
            SearchMode::Alternate => "Alternate",
        }
    }
}

impl fmt::Display for SearchMode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.display_name())
    }
}

/// Mitochondrial model type.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum MitoModel {
    /// No specific mito model
    #[default]
    None,
    /// Mammalian mitochondrial model
    Mammal,
    /// Vertebrate mitochondrial model
    Vertebrate,
}

impl MitoModel {
    /// Returns the model string.
    pub fn as_str(&self) -> &'static str {
        match self {
            MitoModel::None => "",
            MitoModel::Mammal => "mammal",
            MitoModel::Vertebrate => "vert",
        }
    }

    /// Parse from string.
    pub fn from_str(s: &str) -> Option<Self> {
        match s.to_lowercase().as_str() {
            "" => Some(MitoModel::None),
            "mammal" | "mammalian" => Some(MitoModel::Mammal),
            "vert" | "vertebrate" => Some(MitoModel::Vertebrate),
            _ => None,
        }
    }
}

impl fmt::Display for MitoModel {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_str())
    }
}

/// Covariance model mode for second-pass scanning.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum CMMode {
    /// Use Cove for second pass
    #[default]
    Cove,
    /// Use Infernal for second pass
    Infernal,
    /// No covariance model (first-pass only)
    None,
}

impl CMMode {
    /// Returns the mode string.
    pub fn as_str(&self) -> &'static str {
        match self {
            CMMode::Cove => "cove",
            CMMode::Infernal => "infernal",
            CMMode::None => "",
        }
    }

    /// Parse from string.
    pub fn from_str(s: &str) -> Option<Self> {
        match s.to_lowercase().as_str() {
            "cove" => Some(CMMode::Cove),
            "infernal" => Some(CMMode::Infernal),
            "" | "none" => Some(CMMode::None),
            _ => None,
        }
    }

    /// Check if using Cove mode.
    pub fn is_cove(&self) -> bool {
        matches!(self, CMMode::Cove)
    }

    /// Check if using Infernal mode.
    pub fn is_infernal(&self) -> bool {
        matches!(self, CMMode::Infernal)
    }
}

impl fmt::Display for CMMode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_str())
    }
}

/// Main options structure for tRNAscan-SE.
///
/// Contains all configuration options for the tRNA scanning pipeline,
/// including file paths, search modes, thresholds, and output settings.
#[derive(Debug, Clone)]
pub struct Options {
    // Input files
    /// Input FASTA file(s)
    pub fasta_file: PathBuf,
    /// Legacy field for compatibility
    pub fafile: String,
    /// Multiple input files flag
    pub multiple_files: bool,

    // Output files
    /// Main output file ("-" for stdout)
    pub out_file: String,
    /// Results to stdout flag
    pub results_to_stdout: bool,
    /// Secondary structure output file
    pub all_struct_file: String,
    /// BED format output file
    pub bed_file: String,
    /// GFF format output file
    pub gff_file: String,
    /// FASTA output file for predicted tRNAs
    pub output_fasta_file: String,
    /// Isotype-specific results file
    pub isotype_specific_file: String,
    /// Split fragment results file
    pub split_fragment_file: String,
    /// Log file path
    pub log_file: String,
    /// Statistics file path
    pub stats_file: String,
    /// Odd structure file
    pub odd_struct_file: String,
    /// Verbose output file
    pub verb_file: String,
    /// First-pass results file
    pub firstpass_result_file: String,
    /// First-pass flanking file
    pub firstpass_flanking_file: String,
    /// Second-pass intermediate results file
    pub secondpass_int_result_file: String,
    /// Isotype intermediate results file
    pub isotype_int_result_file: String,
    /// Truncated intermediate results file
    pub truncated_int_result_file: String,
    /// False positives file
    pub falsepos_file: String,
    /// Missed sequences file
    pub missed_seq_file: String,
    /// Genetic code file
    pub gc_file: String,

    // Output format flags
    /// ACeDB output format
    pub ace_output: bool,
    /// Brief output (no column headers)
    pub brief_output: bool,
    /// Quiet mode (no credits/options display)
    pub quiet_mode: bool,
    /// Display progress during search
    pub display_progress: bool,
    /// Save progress to log file
    pub save_progress: bool,

    // Save flags
    /// Save statistics
    pub save_stats: bool,
    /// Save odd structures
    pub save_odd_struct: bool,
    /// Save all secondary structures
    pub save_all_struct: bool,
    /// Save verbose output
    pub save_verbose: bool,
    /// Save first-pass results
    pub save_firstpass_res: bool,
    /// Save false positives
    pub save_falsepos: bool,
    /// Save missed sequences
    pub save_missed: bool,
    /// Save source of first-pass hit
    pub save_source: bool,

    // Sequence filtering
    /// Sequence key pattern to match
    pub seq_key: String,
    /// Raw user-input key
    pub raw_seq_key: String,
    /// Start at key instead of matching
    pub start_at_key: bool,

    // Search mode flags
    /// Run tRNAscan first pass
    pub tscan_mode: bool,
    /// Run EufindtRNA first pass
    pub eufind_mode: bool,
    /// Use strict tRNAscan parameters
    pub strict_params: bool,
    /// Run Infernal as first pass
    pub infernal_fp: bool,
    /// Covariance model mode for second pass
    pub cm_mode: CMMode,
    /// Use HMM filter with Infernal
    pub hmm_filter: bool,
    /// Second pass label for display
    pub second_pass_label: String,

    // Organism/search mode
    /// Search mode (euk, bacteria, archaea, etc.)
    pub search_mode: SearchMode,
    /// Disable isotype-specific CM scan
    pub no_isotype: bool,
    /// Mitochondrial model selection
    pub mito_model: MitoModel,
    /// Organellar mode flag (legacy)
    pub org_mode: bool,

    // Genetic code
    /// Use alternate genetic code
    pub alt_gcode: bool,

    // Output options
    /// Output codon instead of anticodon
    pub output_codon: bool,
    /// Display detailed prediction info
    pub detail: bool,
    /// Output Infernal score
    pub infernal_score: bool,

    // Thresholds and parameters
    /// Default padding for first-pass hits
    pub default_padding: usize,
    /// Padding for first-pass hits
    pub padding: usize,
    /// Default max intron + variable loop length
    pub def_max_int_len: usize,
    /// Max intron + variable loop length
    pub max_int_len: usize,

    // File handling
    /// Use previous tRNAscan run results
    pub use_prev_ts_run: bool,
    /// Prompt before overwriting files
    pub prompt_for_overwrite: bool,

    // --- C tRNAscan-SE flags with no prior Rust home (parse + store; some deferred) ---
    /// --mt <model>: mito tRNA models for cytosolic/mito determination
    pub mt_model: String,
    /// -L: legacy search method (tRNAscan + EufindtRNA + COVE)
    pub legacy_mode: bool,
    /// -D --nopseudo: disable pseudogene checking
    pub disable_pseudo: bool,
    /// -U: search using alternate models defined in config file
    pub use_alternate_models: bool,
    /// --nomerge: keep redundant tRNAscan 1.3 hits
    pub nomerge: bool,
    /// -c --conf: configuration file path
    pub conf_file: String,
    /// -p --prefix: default-output-file name prefix
    pub output_prefix: String,
    /// --tmode <mode>: explicit tRNAscan param mode (R or S)
    pub tscan_strictness: String,
    /// --emode <mode>: explicit EufindtRNA param mode (R, N, or S)
    pub eufind_strictness: String,
    /// --iscore <score>: manual EufindtRNA intermediate cutoff score
    pub eufind_intermediate_score: Option<f64>,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            // Input files
            fasta_file: PathBuf::new(),
            fafile: String::new(),
            multiple_files: false,

            // Output files - "-" means stdout
            out_file: "-".to_string(),
            results_to_stdout: true,
            all_struct_file: String::new(),
            bed_file: String::new(),
            gff_file: String::new(),
            output_fasta_file: String::new(),
            isotype_specific_file: String::new(),
            split_fragment_file: String::new(),
            log_file: String::new(),
            stats_file: String::new(),
            odd_struct_file: String::new(),
            verb_file: String::new(),
            firstpass_result_file: String::new(),
            firstpass_flanking_file: String::new(),
            secondpass_int_result_file: String::new(),
            isotype_int_result_file: String::new(),
            truncated_int_result_file: String::new(),
            falsepos_file: String::new(),
            missed_seq_file: String::new(),
            gc_file: String::new(),

            // Output format flags
            ace_output: false,
            brief_output: false,
            quiet_mode: false,
            display_progress: false,
            save_progress: false,

            // Save flags
            save_stats: false,
            save_odd_struct: false,
            save_all_struct: false,
            save_verbose: false,
            save_firstpass_res: false,
            save_falsepos: false,
            save_missed: false,
            save_source: false,

            // Sequence filtering
            seq_key: String::new(),
            raw_seq_key: String::new(),
            start_at_key: false,

            // Search mode flags - defaults match Perl
            tscan_mode: true,
            eufind_mode: true,
            strict_params: true,
            infernal_fp: false,
            cm_mode: CMMode::Cove,
            hmm_filter: false,
            second_pass_label: "Cove".to_string(),

            // Organism/search mode
            search_mode: SearchMode::Eukaryotic,
            no_isotype: false,
            mito_model: MitoModel::None,
            org_mode: false,

            // Genetic code
            alt_gcode: false,

            // Output options
            output_codon: false,
            detail: false,
            infernal_score: false,

            // Thresholds - match Perl defaults
            default_padding: 10,
            padding: 10,
            def_max_int_len: 200,
            max_int_len: 200,

            // File handling
            use_prev_ts_run: false,
            prompt_for_overwrite: true,

            // C flags with no prior Rust home
            mt_model: String::new(),
            legacy_mode: false,
            disable_pseudo: false,
            use_alternate_models: false,
            nomerge: false,
            conf_file: String::new(),
            output_prefix: String::new(),
            tscan_strictness: String::new(),
            eufind_strictness: String::new(),
            eufind_intermediate_score: None,
        }
    }
}

impl Options {
    /// Create a new Options instance with default values.
    pub fn new() -> Self {
        Self::default()
    }

    // Mode query methods

    /// Check if in eukaryotic mode.
    pub fn is_euk_mode(&self) -> bool {
        matches!(self.search_mode, SearchMode::Eukaryotic)
    }

    /// Check if in bacterial mode.
    pub fn is_bact_mode(&self) -> bool {
        matches!(self.search_mode, SearchMode::Bacterial)
    }

    /// Check if in archaeal mode.
    pub fn is_arch_mode(&self) -> bool {
        matches!(self.search_mode, SearchMode::Archaeal)
    }

    /// Check if in mitochondrial mode.
    pub fn is_mito_mode(&self) -> bool {
        matches!(self.search_mode, SearchMode::Mitochondrial)
    }

    /// Check if in mammalian mitochondrial mode.
    pub fn is_mito_mammal_mode(&self) -> bool {
        self.is_mito_mode() && matches!(self.mito_model, MitoModel::Mammal)
    }

    /// Check if in vertebrate mitochondrial mode.
    pub fn is_mito_vert_mode(&self) -> bool {
        self.is_mito_mode() && matches!(self.mito_model, MitoModel::Vertebrate)
    }

    /// Check if in organellar mode.
    pub fn is_organelle_mode(&self) -> bool {
        matches!(self.search_mode, SearchMode::Organellar)
    }

    /// Check if in general mode.
    pub fn is_general_mode(&self) -> bool {
        matches!(self.search_mode, SearchMode::General)
    }

    /// Check if in metagenome mode.
    pub fn is_metagenome_mode(&self) -> bool {
        matches!(self.search_mode, SearchMode::Metagenome)
    }

    /// Check if in numt mode.
    pub fn is_numt_mode(&self) -> bool {
        matches!(self.search_mode, SearchMode::Numt)
    }

    /// Check if in alternate mode.
    pub fn is_alternate_mode(&self) -> bool {
        matches!(self.search_mode, SearchMode::Alternate)
    }

    /// Check if using Cove mode.
    pub fn is_cove_mode(&self) -> bool {
        self.cm_mode.is_cove()
    }

    /// Check if using Infernal mode.
    pub fn is_infernal_mode(&self) -> bool {
        self.cm_mode.is_infernal()
    }

    // Setters

    /// Set the output file.
    pub fn set_out_file(&mut self, path: &str) {
        self.out_file = path.to_string();
        self.results_to_stdout = path == "-";
    }

    /// Set the search mode.
    pub fn set_search_mode(&mut self, mode: SearchMode) {
        self.search_mode = mode;
    }

    /// Set the CM mode for second pass.
    pub fn set_cm_mode(&mut self, mode: CMMode) {
        self.cm_mode = mode;
        self.second_pass_label = match mode {
            CMMode::Cove => "Cove".to_string(),
            CMMode::Infernal => "Infernal".to_string(),
            CMMode::None => String::new(),
        };
    }

    /// Set mitochondrial model.
    pub fn set_mito_model(&mut self, model: MitoModel) {
        self.mito_model = model;
    }

    /// Set padding for first-pass hits.
    pub fn set_padding(&mut self, padding: usize) {
        self.padding = padding;
    }

    /// Set max intron length.
    pub fn set_max_intron_len(&mut self, len: usize) {
        self.max_int_len = len;
    }

    /// Set sequence key for filtering.
    pub fn set_seq_key(&mut self, key: &str) {
        self.raw_seq_key = key.to_string();
        // Convert to regex pattern if needed
        self.seq_key = if key.is_empty() {
            String::new()
        } else {
            key.to_string()
        };
    }

    /// Set the FASTA input file.
    pub fn set_fasta_file(&mut self, path: &Path) {
        self.fasta_file = path.to_path_buf();
        self.fafile = path.to_string_lossy().to_string();
    }

    /// Set the log file.
    pub fn set_log_file(&mut self, path: &str) {
        self.log_file = path.to_string();
        self.save_progress = !path.is_empty();
    }

    /// Set the stats file.
    pub fn set_stats_file(&mut self, path: &str) {
        self.stats_file = path.to_string();
        self.save_stats = !path.is_empty();
    }

    /// Set the secondary structure file.
    pub fn set_struct_file(&mut self, path: &str) {
        self.all_struct_file = path.to_string();
        self.save_all_struct = !path.is_empty();
    }

    /// Set the BED output file.
    pub fn set_bed_file(&mut self, path: &str) {
        self.bed_file = path.to_string();
    }

    /// Set the GFF output file.
    pub fn set_gff_file(&mut self, path: &str) {
        self.gff_file = path.to_string();
    }

    /// Set the FASTA output file.
    pub fn set_output_fasta_file(&mut self, path: &str) {
        self.output_fasta_file = path.to_string();
    }

    /// Set the isotype-specific output file.
    pub fn set_isotype_file(&mut self, path: &str) {
        self.isotype_specific_file = path.to_string();
    }

    /// Set the split fragment output file.
    pub fn set_split_fragment_file(&mut self, path: &str) {
        self.split_fragment_file = path.to_string();
    }

    /// Set the genetic code file.
    pub fn set_gc_file(&mut self, path: &str) {
        self.gc_file = path.to_string();
        self.alt_gcode = !path.is_empty();
    }

    /// Set first-pass results file.
    pub fn set_firstpass_result_file(&mut self, path: &str) {
        self.firstpass_result_file = path.to_string();
        self.save_firstpass_res = !path.is_empty();
    }

    /// Set odd structure file.
    pub fn set_odd_struct_file(&mut self, path: &str) {
        self.odd_struct_file = path.to_string();
        self.save_odd_struct = !path.is_empty();
    }

    /// Set false positives file.
    pub fn set_falsepos_file(&mut self, path: &str) {
        self.falsepos_file = path.to_string();
        self.save_falsepos = !path.is_empty();
    }

    /// Set missed sequences file.
    pub fn set_missed_file(&mut self, path: &str) {
        self.missed_seq_file = path.to_string();
        self.save_missed = !path.is_empty();
    }

    /// Set verbose output file.
    pub fn set_verbose_file(&mut self, path: &str) {
        self.verb_file = path.to_string();
        self.save_verbose = !path.is_empty();
    }

    // Validation

    /// Validate the options configuration.
    pub fn validate(&self) -> Result<(), OptionsError> {
        // Check for conflicting first-pass modes
        if !self.tscan_mode && !self.eufind_mode && !self.infernal_fp {
            if self.cm_mode == CMMode::None {
                return Err(OptionsError::ConflictingOptions(
                    "No search method enabled".to_string(),
                ));
            }
        }

        // Check mito model is only used with mito mode
        if self.mito_model != MitoModel::None && !self.is_mito_mode() && !self.is_euk_mode() {
            return Err(OptionsError::ConflictingOptions(
                "Mito model requires mito or eukaryotic mode".to_string(),
            ));
        }

        // Validate padding
        if self.padding == 0 {
            return Err(OptionsError::InvalidValue(
                "Padding must be greater than 0".to_string(),
            ));
        }

        Ok(())
    }

    /// Resolve relative paths to absolute paths.
    pub fn resolve_paths(&mut self, base_dir: &Path) {
        if !self.fasta_file.as_os_str().is_empty() && self.fasta_file.is_relative() {
            self.fasta_file = base_dir.join(&self.fasta_file);
        }
    }

    /// Get a formatted filename for display (handles "-" for stdout).
    fn format_filename(path: &str) -> &str {
        if path == "-" {
            "Standard output"
        } else {
            path
        }
    }

    /// Display the run options configuration.
    ///
    /// This outputs the current configuration in a format similar to
    /// the original Perl implementation's display_run_options.
    pub fn display_run_options<W: Write>(
        &self,
        writer: &mut W,
        cm_cutoff: f64,
        default_cm_cutoff: f64,
        organelle_cm_cutoff: f64,
        infernal_fp_cutoff: f64,
        temp_dir: &str,
        input_files: &[String],
        check_for_introns: bool,
        check_for_split_halves: bool,
        cm_files: &[String],
        skip_pseudo_filter: bool,
        get_hmm_score: bool,
        keep_tscan_repeats: bool,
        eufind_params: &str,
        eufind_intscore: f64,
        relaxed_param: &str,
        strict_param: &str,
    ) -> std::io::Result<()> {
        let separator = "-".repeat(60);

        writeln!(writer, "{}", separator)?;
        writeln!(
            writer,
            "Sequence file(s) to search:        {}",
            input_files.join(", ")
        )?;

        // Sequence key filtering
        if !self.seq_key.is_empty() && self.seq_key != r"\S*" {
            if self.start_at_key {
                writeln!(
                    writer,
                    "Starting at sequence name:         {}",
                    self.raw_seq_key
                )?;
            } else {
                writeln!(
                    writer,
                    "Search only names matching:        {}",
                    self.raw_seq_key
                )?;
            }
        }

        // Search mode
        write!(writer, "Search Mode:                       ")?;
        let mode_str = match self.search_mode {
            SearchMode::Mitochondrial if self.mito_model == MitoModel::Mammal => {
                "Mitochondrial in mammals"
            }
            SearchMode::Mitochondrial if self.mito_model == MitoModel::Vertebrate => {
                "Mitochondrial in vertebrates"
            }
            _ => self.search_mode.display_name(),
        };
        writeln!(writer, "{}", mode_str)?;

        // Output file
        writeln!(
            writer,
            "Results written to:                {}",
            Self::format_filename(&self.out_file)
        )?;

        // Output format
        write!(writer, "Output format:                     ")?;
        if self.ace_output {
            writeln!(writer, "ACeDB")?;
        } else {
            writeln!(writer, "Tabular")?;
        }

        // Search method
        write!(writer, "Searching with:                    ")?;
        let search_method = self.format_search_method();
        writeln!(writer, "{}", search_method)?;

        // Sensitivity mode for Infernal
        if self.cm_mode == CMMode::Infernal && !self.infernal_fp {
            if self.hmm_filter {
                writeln!(writer, "                                   Fast mode")?;
            } else {
                writeln!(writer, "                                   Maximum sensitivity mode")?;
            }
        }

        // Isotype-specific scan
        if !self.no_isotype
            && (self.is_euk_mode() || self.is_bact_mode() || self.is_arch_mode())
        {
            writeln!(writer, "Isotype-specific model scan:       Yes")?;
        } else if self.no_isotype
            && (self.is_euk_mode() || self.is_bact_mode() || self.is_arch_mode())
        {
            writeln!(writer, "Isotype-specific model scan:       No")?;
        }

        // Mito isotype model
        if self.is_euk_mode() && self.mito_model != MitoModel::None {
            writeln!(
                writer,
                "Mito isotype model scan:           {}",
                self.mito_model
            )?;
        }

        // Intron/split scanning
        if check_for_introns {
            writeln!(writer, "Scan for noncanonical introns")?;
        }
        if check_for_split_halves {
            writeln!(writer, "Scan for fragments of split tRNAs")?;
        }

        // Covariance models
        for (i, cm_file) in cm_files.iter().enumerate() {
            if i == 0 {
                writeln!(writer, "Covariance model:                  {}", cm_file)?;
            } else {
                writeln!(writer, "                                   {}", cm_file)?;
            }
        }

        // CM cutoff
        let show_cutoff = if self.is_mito_mode() || self.is_organelle_mode() {
            cm_cutoff != organelle_cm_cutoff
        } else {
            cm_cutoff != default_cm_cutoff
        };
        if show_cutoff {
            writeln!(writer, "tRNA covariance model search       ")?;
            writeln!(writer, "    cutoff score:                  {}", cm_cutoff)?;
        }

        // Previous run
        if self.use_prev_ts_run {
            writeln!(writer, "Using previous")?;
            writeln!(
                writer,
                "tabular output file:               {}",
                self.firstpass_result_file
            )?;
        }

        // tRNAscan parameters
        if self.tscan_mode {
            write!(writer, "tRNAscan parameters:               ")?;
            if self.strict_params {
                writeln!(writer, "Strict")?;
            } else {
                writeln!(writer, "Relaxed")?;
            }
        }

        // EufindtRNA parameters
        if self.eufind_mode {
            write!(writer, "EufindtRNA parameters:             ")?;
            if eufind_params == relaxed_param {
                writeln!(writer, "Relaxed (Int Cutoff= {})", eufind_intscore)?;
            } else if eufind_params.is_empty() {
                writeln!(writer, "Normal")?;
            } else if eufind_params == strict_param {
                writeln!(writer, "Strict")?;
            } else {
                writeln!(writer, "?")?;
            }
        }

        // Infernal first pass cutoff
        if self.infernal_fp {
            writeln!(
                writer,
                "Infernal first pass cutoff score:  {}",
                infernal_fp_cutoff
            )?;
        }

        // Custom padding
        if self.padding != self.default_padding {
            writeln!(
                writer,
                "First-pass tRNA hit padding:       {} bp",
                self.padding
            )?;
        }

        // Alternate genetic code
        if self.alt_gcode {
            writeln!(
                writer,
                "Alternate transl code used:        from file {}",
                self.gc_file
            )?;
        }

        // Temp directory
        writeln!(writer)?;
        writeln!(writer, "Temporary directory:               {}", temp_dir)?;

        // Secondary structure output
        if self.save_all_struct {
            writeln!(writer, "tRNA secondary structure")?;
            writeln!(
                writer,
                "    predictions saved to:          {}",
                Self::format_filename(&self.all_struct_file)
            )?;
        }

        // BED output
        if !self.bed_file.is_empty() {
            writeln!(
                writer,
                "tRNA predictions saved to:         {}",
                Self::format_filename(&self.bed_file)
            )?;
        }

        // GFF output
        if !self.gff_file.is_empty() {
            writeln!(
                writer,
                "tRNA predictions saved to:         {}",
                Self::format_filename(&self.gff_file)
            )?;
        }

        // FASTA output
        if !self.output_fasta_file.is_empty() {
            writeln!(writer, "Predicted tRNA sequences")?;
            writeln!(
                writer,
                "    saved to:                      {}",
                Self::format_filename(&self.output_fasta_file)
            )?;
        }

        // Isotype-specific output
        if !self.isotype_specific_file.is_empty() {
            writeln!(writer, "Isotype specific")?;
            writeln!(
                writer,
                "    predictions saved to:          {}",
                Self::format_filename(&self.isotype_specific_file)
            )?;
        }

        // Split fragment output
        if !self.split_fragment_file.is_empty() {
            writeln!(writer, "Split tRNA fragment")?;
            writeln!(
                writer,
                "    predictions saved to:          {}",
                Self::format_filename(&self.split_fragment_file)
            )?;
        }

        // Odd structure output
        if self.save_odd_struct {
            writeln!(writer, "Sec structures for tRNAs")?;
            writeln!(
                writer,
                " with no anticodon predictn: {}",
                self.odd_struct_file
            )?;
        }

        // First-pass results
        if self.save_firstpass_res {
            writeln!(
                writer,
                "First-pass results saved in:      {}",
                self.firstpass_result_file
            )?;
        }

        // Log file
        if !self.log_file.is_empty() {
            writeln!(writer, "Search log saved in:               {}", self.log_file)?;
        }

        // Stats file
        if self.save_stats {
            writeln!(
                writer,
                "Search statistics saved in:        {}",
                self.stats_file
            )?;
        }

        // False positives file
        if self.save_falsepos {
            writeln!(
                writer,
                "False positives saved in:          {}",
                self.falsepos_file
            )?;
        }

        // Missed sequences file
        if self.save_missed {
            writeln!(
                writer,
                "Seqs with 0 hits saved in:         {}",
                self.missed_seq_file
            )?;
        }

        // Additional options
        if skip_pseudo_filter || get_hmm_score || keep_tscan_repeats {
            writeln!(writer)?;
        }

        if self.max_int_len != self.def_max_int_len {
            writeln!(
                writer,
                "Max intron + var. length:          {}",
                self.max_int_len
            )?;
        }

        if skip_pseudo_filter {
            writeln!(writer, "Pseudogene checking disabled")?;
        }

        if get_hmm_score {
            writeln!(writer, "Reporting HMM/2' structure score breakdown")?;
        }

        if keep_tscan_repeats {
            writeln!(writer, "Redundant tRNAscan hits not merged")?;
        }

        writeln!(writer, "{}", separator)?;
        writeln!(writer)?;

        Ok(())
    }

    /// Format the search method description.
    fn format_search_method(&self) -> String {
        let second_pass = if self.cm_mode != CMMode::None {
            format!("->{}", self.second_pass_label)
        } else {
            String::new()
        };

        if self.eufind_mode {
            if self.tscan_mode {
                if self.cm_mode != CMMode::None {
                    format!("tRNAscan + EufindtRNA {}", second_pass)
                } else {
                    "tRNAscan + EufindtRNA".to_string()
                }
            } else if self.cm_mode != CMMode::None {
                format!("EufindtRNA{}", second_pass)
            } else {
                "EufindtRNA only".to_string()
            }
        } else if self.tscan_mode {
            if self.cm_mode != CMMode::None {
                format!("tRNAscan{}", second_pass)
            } else {
                "tRNAscan only".to_string()
            }
        } else if self.infernal_fp {
            if self.cm_mode != CMMode::None {
                format!("Infernal First Pass{}", second_pass)
            } else {
                "Infernal First Pass only".to_string()
            }
        } else {
            format!("{} single-pass scan", self.second_pass_label)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_options() {
        let opts = Options::default();
        assert!(opts.tscan_mode);
        assert!(opts.eufind_mode);
        assert!(opts.strict_params);
        assert!(!opts.infernal_fp);
        assert_eq!(opts.cm_mode, CMMode::Cove);
        assert_eq!(opts.search_mode, SearchMode::Eukaryotic);
        assert_eq!(opts.padding, 10);
        assert_eq!(opts.max_int_len, 200);
        assert_eq!(opts.out_file, "-");
        assert!(opts.results_to_stdout);
    }

    #[test]
    fn test_search_mode_parsing() {
        assert_eq!(SearchMode::from_str("euk"), Some(SearchMode::Eukaryotic));
        assert_eq!(SearchMode::from_str("bacteria"), Some(SearchMode::Bacterial));
        assert_eq!(SearchMode::from_str("archaea"), Some(SearchMode::Archaeal));
        assert_eq!(SearchMode::from_str("mito"), Some(SearchMode::Mitochondrial));
        assert_eq!(SearchMode::from_str("general"), Some(SearchMode::General));
        assert_eq!(SearchMode::from_str("invalid"), None);
    }

    #[test]
    fn test_mode_queries() {
        let mut opts = Options::default();

        assert!(opts.is_euk_mode());
        assert!(!opts.is_bact_mode());

        opts.search_mode = SearchMode::Bacterial;
        assert!(opts.is_bact_mode());
        assert!(!opts.is_euk_mode());

        opts.search_mode = SearchMode::Mitochondrial;
        opts.mito_model = MitoModel::Mammal;
        assert!(opts.is_mito_mode());
        assert!(opts.is_mito_mammal_mode());
        assert!(!opts.is_mito_vert_mode());
    }

    #[test]
    fn test_cm_mode() {
        let mut opts = Options::default();

        assert!(opts.is_cove_mode());
        assert!(!opts.is_infernal_mode());

        opts.set_cm_mode(CMMode::Infernal);
        assert!(opts.is_infernal_mode());
        assert!(!opts.is_cove_mode());
        assert_eq!(opts.second_pass_label, "Infernal");
    }

    #[test]
    fn test_set_out_file() {
        let mut opts = Options::default();

        opts.set_out_file("output.txt");
        assert_eq!(opts.out_file, "output.txt");
        assert!(!opts.results_to_stdout);

        opts.set_out_file("-");
        assert_eq!(opts.out_file, "-");
        assert!(opts.results_to_stdout);
    }

    #[test]
    fn test_validation() {
        let opts = Options::default();
        assert!(opts.validate().is_ok());

        let mut bad_opts = Options::default();
        bad_opts.tscan_mode = false;
        bad_opts.eufind_mode = false;
        bad_opts.infernal_fp = false;
        bad_opts.cm_mode = CMMode::None;
        assert!(bad_opts.validate().is_err());
    }

    #[test]
    fn test_mito_model() {
        assert_eq!(MitoModel::from_str("mammal"), Some(MitoModel::Mammal));
        assert_eq!(MitoModel::from_str("vert"), Some(MitoModel::Vertebrate));
        assert_eq!(MitoModel::from_str(""), Some(MitoModel::None));
        assert_eq!(MitoModel::from_str("invalid"), None);
    }
}
