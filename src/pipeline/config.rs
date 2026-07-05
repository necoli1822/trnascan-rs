//! Configuration and command-line options for tRNAscan-SE.
//!
//! This module provides configuration structures and parsing for:
//! - Search modes (Eukaryotic, Bacterial, Archaeal, etc.)
//! - First-pass detection methods
//! - Output file configurations
//! - Score cutoffs and algorithm parameters

use std::path::PathBuf;

/// Search mode for tRNA detection.
///
/// Different organisms have different tRNA characteristics. The search mode
/// determines which covariance models and first-pass methods to use.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SearchMode {
    /// Eukaryotic mode (-E): Uses EuFindtRNA first pass, standard eukaryotic CM
    Eukaryotic,

    /// Bacterial mode (-B): Uses tRNAscan 1.4 first pass, bacterial CM
    Bacterial,

    /// Archaeal mode (-A): Uses tRNAscan 1.4 first pass, archaeal CM
    Archaeal,

    /// Organellar mode (-O): Mitochondrial/chloroplast tRNAs
    Organellar,

    /// Mitochondrial mode (-M): Specific mitochondrial models
    Mitochondrial,

    /// General mode (-G): All domains, more sensitive but slower
    General,

    /// Legacy infernal-only mode (-I): Skips first-pass, uses only Infernal
    InfernalOnly,

    /// Custom mode: User-specified models and parameters
    Custom,
}

impl SearchMode {
    /// Get the default CM model file name for this search mode.
    pub fn default_cm_model(&self) -> &'static str {
        match self {
            SearchMode::Eukaryotic => "TRNAinf-euk.cm",
            SearchMode::Bacterial => "TRNAinf-bact.cm",
            SearchMode::Archaeal => "TRNAinf-arch.cm",
            SearchMode::Organellar => "TRNAinf-mito.cm",
            SearchMode::Mitochondrial => "TRNAinf-mito.cm",
            SearchMode::General => "TRNAinf-all.cm",
            SearchMode::InfernalOnly => "TRNAinf-all.cm",
            SearchMode::Custom => "TRNAinf-all.cm",
        }
    }

    /// Get the short mode name for display.
    pub fn short_name(&self) -> &'static str {
        match self {
            SearchMode::Eukaryotic => "Eukaryotic",
            SearchMode::Bacterial => "Bacterial",
            SearchMode::Archaeal => "Archaeal",
            SearchMode::Organellar => "Organellar",
            SearchMode::Mitochondrial => "Mitochondrial",
            SearchMode::General => "General",
            SearchMode::InfernalOnly => "Infernal",
            SearchMode::Custom => "Custom",
        }
    }

    /// Get the command-line flag for this mode.
    pub fn flag(&self) -> &'static str {
        match self {
            SearchMode::Eukaryotic => "-E",
            SearchMode::Bacterial => "-B",
            SearchMode::Archaeal => "-A",
            SearchMode::Organellar => "-O",
            SearchMode::Mitochondrial => "-M",
            SearchMode::General => "-G",
            SearchMode::InfernalOnly => "-I",
            SearchMode::Custom => "",
        }
    }

    /// Get the default first-pass method for this mode.
    pub fn default_first_pass(&self) -> FirstPassMethod {
        match self {
            SearchMode::Eukaryotic => FirstPassMethod::EuFindtRNA,
            SearchMode::Bacterial | SearchMode::Archaeal => FirstPassMethod::TrnaScan14,
            SearchMode::Organellar | SearchMode::Mitochondrial => FirstPassMethod::TrnaScan14,
            SearchMode::General => FirstPassMethod::Infernal,
            SearchMode::InfernalOnly => FirstPassMethod::Infernal,
            SearchMode::Custom => FirstPassMethod::Infernal,
        }
    }

    /// Get the default score cutoff for this mode.
    pub fn default_score_cutoff(&self) -> f64 {
        match self {
            SearchMode::Eukaryotic => 20.0,
            SearchMode::Bacterial => 20.0,
            SearchMode::Archaeal => 20.0,
            SearchMode::Organellar => 15.0,
            SearchMode::Mitochondrial => 15.0,
            SearchMode::General => 20.0,
            SearchMode::InfernalOnly => 20.0,
            SearchMode::Custom => 20.0,
        }
    }
}

impl std::fmt::Display for SearchMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.short_name())
    }
}

/// First-pass detection method.
///
/// The first pass quickly identifies candidate tRNA regions, which are then
/// analyzed in detail with Infernal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FirstPassMethod {
    /// EuFindtRNA (Pavesi algorithm): Detects A-box and B-box promoter regions
    EuFindtRNA,

    /// tRNAscan 1.4: Profile-based detection using weight matrices
    TrnaScan14,

    /// Infernal direct: Use CM scan as first pass (slower but more sensitive)
    Infernal,

    /// Combined: Run both EuFindtRNA and tRNAscan 1.4
    Combined,
}

impl FirstPassMethod {
    /// Get the display name for this method.
    pub fn name(&self) -> &'static str {
        match self {
            FirstPassMethod::EuFindtRNA => "EuFindtRNA",
            FirstPassMethod::TrnaScan14 => "tRNAscan 1.4",
            FirstPassMethod::Infernal => "Infernal",
            FirstPassMethod::Combined => "EuFindtRNA + tRNAscan 1.4",
        }
    }
}

impl std::fmt::Display for FirstPassMethod {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.name())
    }
}

/// Mitochondrial model variant.
///
/// Different organisms have different mitochondrial genetic codes and
/// tRNA characteristics.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MitoModel {
    /// Vertebrate mitochondrial code
    Vertebrate,

    /// Mammalian mitochondrial code
    Mammal,

    /// Yeast mitochondrial code
    Yeast,

    /// Invertebrate mitochondrial code
    Invertebrate,

    /// Mold/protozoan mitochondrial code
    Mold,

    /// Echinoderm mitochondrial code
    Echinoderm,

    /// Plant mitochondrial code
    Plant,
}

impl MitoModel {
    /// Parse a model name from command-line argument.
    pub fn from_name(name: &str) -> Option<Self> {
        match name.to_lowercase().as_str() {
            "vert" | "vertebrate" => Some(MitoModel::Vertebrate),
            "mammal" | "mammalian" => Some(MitoModel::Mammal),
            "yeast" => Some(MitoModel::Yeast),
            "invert" | "invertebrate" => Some(MitoModel::Invertebrate),
            "mold" | "protozoan" => Some(MitoModel::Mold),
            "echinoderm" => Some(MitoModel::Echinoderm),
            "plant" => Some(MitoModel::Plant),
            _ => None,
        }
    }

    /// Get the display name for this model.
    pub fn name(&self) -> &'static str {
        match self {
            MitoModel::Vertebrate => "vertebrate",
            MitoModel::Mammal => "mammal",
            MitoModel::Yeast => "yeast",
            MitoModel::Invertebrate => "invertebrate",
            MitoModel::Mold => "mold",
            MitoModel::Echinoderm => "echinoderm",
            MitoModel::Plant => "plant",
        }
    }

    /// Get the CM model file name for this variant.
    pub fn cm_model(&self) -> &'static str {
        match self {
            MitoModel::Vertebrate => "TRNAinf-mito-vert.cm",
            MitoModel::Mammal => "TRNAinf-mito-mammal.cm",
            MitoModel::Yeast => "TRNAinf-mito-yeast.cm",
            MitoModel::Invertebrate => "TRNAinf-mito-invert.cm",
            MitoModel::Mold => "TRNAinf-mito-mold.cm",
            MitoModel::Echinoderm => "TRNAinf-mito-echinoderm.cm",
            MitoModel::Plant => "TRNAinf-mito-plant.cm",
        }
    }
}

impl std::fmt::Display for MitoModel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.name())
    }
}

/// Output file configuration.
#[derive(Debug, Clone, Default)]
pub struct OutputConfig {
    /// Main tabular output file (.out)
    pub output_file: Option<PathBuf>,

    /// Secondary structure output file (.ss)
    pub ss_file: Option<PathBuf>,

    /// Statistics output file (.stats)
    pub stats_file: Option<PathBuf>,

    /// BED format output file (.bed)
    pub bed_file: Option<PathBuf>,

    /// Isotype-specific output file (.iso)
    pub iso_file: Option<PathBuf>,

    /// Log file for verbose output (.log)
    pub log_file: Option<PathBuf>,

    /// Prefix for auto-generated output file names
    pub prefix: Option<String>,

    /// Force overwrite existing files
    pub force_overwrite: bool,
}

impl OutputConfig {
    /// Create new output config with no output files.
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the output file prefix.
    pub fn with_prefix(mut self, prefix: impl Into<String>) -> Self {
        self.prefix = Some(prefix.into());
        self
    }

    /// Generate default output file names from prefix.
    pub fn generate_defaults(&mut self) {
        if let Some(ref prefix) = self.prefix {
            if self.output_file.is_none() {
                self.output_file = Some(PathBuf::from(format!("{}-tRNAs.out", prefix)));
            }
            if self.ss_file.is_none() {
                self.ss_file = Some(PathBuf::from(format!("{}-tRNAs.ss", prefix)));
            }
            if self.stats_file.is_none() {
                self.stats_file = Some(PathBuf::from(format!("{}-tRNAs.stats", prefix)));
            }
            if self.bed_file.is_none() {
                self.bed_file = Some(PathBuf::from(format!("{}-tRNAs.bed", prefix)));
            }
            if self.iso_file.is_none() {
                self.iso_file = Some(PathBuf::from(format!("{}-tRNAs.iso", prefix)));
            }
        }
    }
}

/// Main tRNAscan-SE configuration structure.
#[derive(Debug, Clone)]
pub struct TrnaScanConfig {
    /// Search mode (Eukaryotic, Bacterial, etc.)
    pub search_mode: SearchMode,

    /// First-pass detection method
    pub first_pass: FirstPassMethod,

    /// Covariance model file path
    pub cm_model: PathBuf,

    /// Selenocysteine CM model file path (optional)
    pub sec_cm_model: Option<PathBuf>,

    /// First-pass score cutoff
    pub first_pass_cutoff: f64,

    /// Infernal score cutoff
    pub score_cutoff: f64,

    /// Output configuration
    pub output: OutputConfig,

    /// Show secondary structure in output
    pub show_secondary_structure: bool,

    /// Show statistics in output
    pub show_stats: bool,

    /// Show HMM/secondary structure score breakdown (-H)
    pub show_breakdown: bool,

    /// Show hit source origin (-y)
    pub show_hit_source: bool,

    /// Brief output mode (--brief)
    pub brief_output: bool,

    /// Detailed output mode (--detail)
    pub detail_output: bool,

    /// Quiet mode (-q): suppress progress output
    pub quiet: bool,

    /// Disable pseudogene detection (-D)
    pub disable_pseudogene: bool,

    /// Scan both strands (default true)
    pub both_strands: bool,

    /// Mitochondrial model (if applicable)
    pub mito_model: Option<MitoModel>,

    /// Padding bases around tRNA for alignment (-z)
    pub padding: usize,

    /// Temporary directory for intermediate files
    pub temp_dir: PathBuf,

    /// Match pattern filter for sequence names (--match)
    pub match_pattern: Option<String>,

    /// Library path for CM models
    pub lib_path: Option<PathBuf>,

    /// Infernal score cutoff for first pass (if using Infernal first pass)
    pub infernal_first_pass_cutoff: f64,

    /// Enable isotype-specific CM scanning
    pub isotype_specific: bool,

    /// Path to the main covariance model file (overrides default model path)
    pub cm_model_path: Option<PathBuf>,

    /// Path to the models directory (overrides default model directory)
    pub models_dir: Option<PathBuf>,
}

impl TrnaScanConfig {
    /// Create a new configuration with the specified search mode.
    pub fn new(search_mode: SearchMode) -> Self {
        let first_pass = search_mode.default_first_pass();
        let score_cutoff = search_mode.default_score_cutoff();

        Self {
            search_mode,
            first_pass,
            cm_model: PathBuf::from(search_mode.default_cm_model()),
            sec_cm_model: None,
            first_pass_cutoff: 10.0,
            score_cutoff,
            output: OutputConfig::new(),
            show_secondary_structure: false,
            show_stats: false,
            show_breakdown: false,
            show_hit_source: false,
            brief_output: false,
            detail_output: false,
            quiet: false,
            disable_pseudogene: false,
            both_strands: true,
            mito_model: None,
            padding: 8,
            temp_dir: std::env::temp_dir(),
            match_pattern: None,
            lib_path: None,
            infernal_first_pass_cutoff: 10.0,
            isotype_specific: true,
            cm_model_path: None,
            models_dir: None,
        }
    }

    /// Builder: Set the output file.
    pub fn with_output_file(mut self, path: Option<PathBuf>) -> Self {
        self.output.output_file = path;
        self
    }

    /// Builder: Set the secondary structure output file.
    pub fn with_ss_file(mut self, path: Option<PathBuf>) -> Self {
        self.output.ss_file = path;
        self
    }

    /// Builder: Set the statistics output file.
    pub fn with_stats_file(mut self, path: Option<PathBuf>) -> Self {
        self.output.stats_file = path;
        self
    }

    /// Builder: Set the BED output file.
    pub fn with_bed_file(mut self, path: Option<PathBuf>) -> Self {
        self.output.bed_file = path;
        self
    }

    /// Builder: Set the isotype output file.
    pub fn with_iso_file(mut self, path: Option<PathBuf>) -> Self {
        self.output.iso_file = path;
        self
    }

    /// Builder: Set the score cutoff.
    pub fn with_score_cutoff(mut self, cutoff: f64) -> Self {
        self.score_cutoff = cutoff;
        self
    }

    /// Builder: Set the CM model path.
    pub fn with_cm_model(mut self, path: PathBuf) -> Self {
        self.cm_model = path;
        self
    }

    /// Builder: Enable secondary structure output.
    pub fn with_secondary_structure(mut self, enabled: bool) -> Self {
        self.show_secondary_structure = enabled;
        self
    }

    /// Builder: Enable statistics output.
    pub fn with_stats(mut self, enabled: bool) -> Self {
        self.show_stats = enabled;
        self
    }

    /// Builder: Enable breakdown output (-H).
    pub fn with_breakdown(mut self, enabled: bool) -> Self {
        self.show_breakdown = enabled;
        self
    }

    /// Builder: Set quiet mode.
    pub fn with_quiet(mut self, quiet: bool) -> Self {
        self.quiet = quiet;
        self
    }

    /// Builder: Set padding.
    pub fn with_padding(mut self, padding: usize) -> Self {
        self.padding = padding;
        self
    }

    /// Builder: Set first-pass method.
    pub fn with_first_pass(mut self, method: FirstPassMethod) -> Self {
        self.first_pass = method;
        self
    }

    /// Builder: Set mitochondrial model.
    pub fn with_mito_model(mut self, model: MitoModel) -> Self {
        self.mito_model = Some(model);
        self
    }

    /// Get the CM model path, searching in standard locations.
    ///
    /// Priority order:
    /// 1. Explicit cm_model_path (highest priority)
    /// 2. models_dir / cm_model
    /// 3. lib_path / models / cm_model
    /// 4. Standard locations:
    ///    - ./original/lib/models/<model>
    ///    - /usr/local/lib/tRNAscan-SE/models/<model>
    ///    - ~/.trnascan-rs/models/<model>
    /// 5. Current directory / cm_model
    pub fn get_model_path(&self) -> Result<PathBuf, String> {
        // Priority 1: Explicit cm_model_path
        if let Some(ref path) = self.cm_model_path {
            if path.exists() {
                return Ok(path.clone());
            } else {
                return Err(format!(
                    "Specified CM model file not found: {}",
                    path.display()
                ));
            }
        }

        // Priority 2: models_dir / cm_model
        if let Some(ref models_dir) = self.models_dir {
            let path = models_dir.join(&self.cm_model);
            if path.exists() {
                return Ok(path);
            }
        }

        // Priority 3: lib_path / models / cm_model
        if let Some(ref lib_path) = self.lib_path {
            let path = lib_path.join("models").join(&self.cm_model);
            if path.exists() {
                return Ok(path);
            }
        }

        // Priority 4: Standard locations
        let model_name = &self.cm_model;
        let standard_locations = [
            PathBuf::from("original/lib/models").join(model_name),
            PathBuf::from("models").join(model_name),
            PathBuf::from("/usr/local/lib/tRNAscan-SE/models").join(model_name),
        ];

        // Add home directory location if HOME env var is set
        let mut search_paths = standard_locations.to_vec();
        if let Ok(home) = std::env::var("HOME") {
            search_paths.push(PathBuf::from(home).join(".trnascan-rs/models").join(model_name));
        }

        for path in &search_paths {
            if path.exists() {
                return Ok(path.clone());
            }
        }

        // Priority 5: Current directory / cm_model
        if self.cm_model.exists() {
            return Ok(self.cm_model.clone());
        }

        // Not found anywhere
        Err(format!(
            "CM model file '{}' not found in any standard location. Searched:\n{}",
            model_name.display(),
            search_paths
                .iter()
                .map(|p| format!("  - {}", p.display()))
                .collect::<Vec<_>>()
                .join("\n")
        ))
    }

    /// Get the effective CM model path, resolving against lib_path if needed.
    pub fn effective_cm_model(&self) -> PathBuf {
        if self.cm_model.is_absolute() {
            self.cm_model.clone()
        } else if let Some(ref lib_path) = self.lib_path {
            lib_path.join("models").join(&self.cm_model)
        } else {
            self.cm_model.clone()
        }
    }

    /// Validate the configuration for consistency.
    pub fn validate(&self) -> Result<(), ConfigError> {
        // Check for conflicting options
        if self.brief_output && self.detail_output {
            return Err(ConfigError::ConflictingOptions(
                "--brief".to_string(),
                "--detail".to_string(),
            ));
        }

        // Validate score cutoff
        if self.score_cutoff < 0.0 {
            return Err(ConfigError::InvalidScoreCutoff(self.score_cutoff));
        }

        // Validate padding
        if self.padding > 100 {
            return Err(ConfigError::InvalidPadding(self.padding));
        }

        Ok(())
    }
}

impl Default for TrnaScanConfig {
    fn default() -> Self {
        Self::new(SearchMode::Eukaryotic)
    }
}

/// Configuration errors.
#[derive(Debug, Clone)]
pub enum ConfigError {
    /// Missing required input file
    MissingInputFile,

    /// Input file not found
    InputFileNotFound(PathBuf),

    /// Cannot open input file
    CannotOpenInput(PathBuf, String),

    /// Cannot write output file
    CannotWriteOutput(PathBuf, String),

    /// Invalid score cutoff value
    InvalidScoreCutoff(f64),

    /// Invalid padding value
    InvalidPadding(usize),

    /// Conflicting options specified
    ConflictingOptions(String, String),

    /// Multiple search modes specified
    MultipleSearchModes,

    /// Unknown command-line option
    UnknownOption(String),

    /// Invalid mitochondrial model
    InvalidMitoModel(String),

    /// CM model file not found
    CmModelNotFound(PathBuf),

    /// Invalid match pattern (regex error)
    InvalidMatchPattern(String),

    /// Permission denied
    PermissionDenied(PathBuf),

    /// Temp directory not writable
    TempDirNotWritable(PathBuf),
}

impl std::fmt::Display for ConfigError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ConfigError::MissingInputFile => {
                write!(f, "No input file specified")
            }
            ConfigError::InputFileNotFound(path) => {
                write!(f, "Cannot open input file '{}'", path.display())
            }
            ConfigError::CannotOpenInput(path, msg) => {
                write!(f, "Cannot read input file '{}': {}", path.display(), msg)
            }
            ConfigError::CannotWriteOutput(path, msg) => {
                write!(f, "Cannot write to output file '{}': {}", path.display(), msg)
            }
            ConfigError::InvalidScoreCutoff(val) => {
                write!(f, "Invalid score cutoff: {}", val)
            }
            ConfigError::InvalidPadding(val) => {
                write!(f, "Invalid padding value: {}", val)
            }
            ConfigError::ConflictingOptions(opt1, opt2) => {
                write!(
                    f,
                    "Conflicting options: {} and {} cannot be used together",
                    opt1, opt2
                )
            }
            ConfigError::MultipleSearchModes => {
                write!(f, "Multiple search modes specified (only one allowed)")
            }
            ConfigError::UnknownOption(opt) => {
                write!(f, "Unknown option '{}'", opt)
            }
            ConfigError::InvalidMitoModel(model) => {
                write!(
                    f,
                    "Invalid mitochondrial model: {}\nValid options: mammal, vert",
                    model
                )
            }
            ConfigError::CmModelNotFound(path) => {
                write!(
                    f,
                    "Unable to open '{}' covariance model file",
                    path.display()
                )
            }
            ConfigError::InvalidMatchPattern(pattern) => {
                write!(f, "Invalid match pattern: {}", pattern)
            }
            ConfigError::PermissionDenied(path) => {
                write!(f, "Permission denied: {}", path.display())
            }
            ConfigError::TempDirNotWritable(path) => {
                write!(f, "Cannot write to temporary directory: {}", path.display())
            }
        }
    }
}

impl std::error::Error for ConfigError {}

/// Find the default CM model file in standard locations.
///
/// Searches for TRNA2.cm (or other common models) in:
/// - ./original/lib/models/
/// - ./models/
/// - /usr/local/lib/tRNAscan-SE/models/
/// - ~/.trnascan-rs/models/
///
/// Returns the first path that exists, or None if not found.
pub fn find_default_model() -> Option<PathBuf> {
    let model_names = [
        "TRNA2.cm",
        "TRNAinf-euk.cm",
        "TRNAinf-bact.cm",
        "TRNAinf-arch.cm",
        "TRNAinf-all.cm",
    ];

    let mut possible_paths = vec![
        PathBuf::from("original/lib/models"),
        PathBuf::from("models"),
        PathBuf::from("/usr/local/lib/tRNAscan-SE/models"),
    ];

    // Add home directory location if HOME env var is set
    if let Ok(home) = std::env::var("HOME") {
        possible_paths.push(PathBuf::from(home).join(".trnascan-rs/models"));
    }

    // Search for each model name in each directory
    for dir in &possible_paths {
        for model_name in &model_names {
            let path = dir.join(model_name);
            if path.exists() {
                return Some(path);
            }
        }
    }

    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_search_mode_defaults() {
        let euk = SearchMode::Eukaryotic;
        assert_eq!(euk.short_name(), "Eukaryotic");
        assert_eq!(euk.flag(), "-E");
        assert_eq!(euk.default_first_pass(), FirstPassMethod::EuFindtRNA);
        assert_eq!(euk.default_score_cutoff(), 20.0);

        let bact = SearchMode::Bacterial;
        assert_eq!(bact.short_name(), "Bacterial");
        assert_eq!(bact.flag(), "-B");
        assert_eq!(bact.default_first_pass(), FirstPassMethod::TrnaScan14);
    }

    #[test]
    fn test_config_builder() {
        let config = TrnaScanConfig::new(SearchMode::Eukaryotic)
            .with_score_cutoff(30.0)
            .with_padding(16)
            .with_quiet(true);

        assert_eq!(config.score_cutoff, 30.0);
        assert_eq!(config.padding, 16);
        assert!(config.quiet);
    }

    #[test]
    fn test_config_validation() {
        // Valid config
        let config = TrnaScanConfig::new(SearchMode::Eukaryotic);
        assert!(config.validate().is_ok());

        // Conflicting options
        let mut config = TrnaScanConfig::new(SearchMode::Eukaryotic);
        config.brief_output = true;
        config.detail_output = true;
        assert!(matches!(
            config.validate(),
            Err(ConfigError::ConflictingOptions(_, _))
        ));

        // Negative score cutoff
        let config = TrnaScanConfig::new(SearchMode::Eukaryotic).with_score_cutoff(-10.0);
        assert!(matches!(
            config.validate(),
            Err(ConfigError::InvalidScoreCutoff(_))
        ));
    }

    #[test]
    fn test_mito_model_parsing() {
        assert_eq!(MitoModel::from_name("mammal"), Some(MitoModel::Mammal));
        assert_eq!(MitoModel::from_name("vert"), Some(MitoModel::Vertebrate));
        assert_eq!(MitoModel::from_name("YEAST"), Some(MitoModel::Yeast));
        assert_eq!(MitoModel::from_name("invalid"), None);
    }

    #[test]
    fn test_output_config_defaults() {
        let mut output = OutputConfig::new().with_prefix("test");
        output.generate_defaults();

        assert_eq!(
            output.output_file,
            Some(PathBuf::from("test-tRNAs.out"))
        );
        assert_eq!(
            output.ss_file,
            Some(PathBuf::from("test-tRNAs.ss"))
        );
        assert_eq!(
            output.stats_file,
            Some(PathBuf::from("test-tRNAs.stats"))
        );
    }

    #[test]
    fn test_first_pass_method_display() {
        assert_eq!(FirstPassMethod::EuFindtRNA.name(), "EuFindtRNA");
        assert_eq!(FirstPassMethod::TrnaScan14.name(), "tRNAscan 1.4");
        assert_eq!(FirstPassMethod::Infernal.name(), "Infernal");
    }

    #[test]
    fn test_model_path_priority() {
        use std::fs;
        use std::io::Write;

        // Create a temporary test directory
        let temp_dir = std::env::temp_dir().join("trnascan_test_models");
        let _ = fs::create_dir_all(&temp_dir);

        // Create a test model file
        let test_model = temp_dir.join("test.cm");
        let mut file = fs::File::create(&test_model).unwrap();
        file.write_all(b"test model content").unwrap();

        // Test with explicit cm_model_path (highest priority)
        let mut config = TrnaScanConfig::new(SearchMode::Eukaryotic);
        config.cm_model_path = Some(test_model.clone());

        let result = config.get_model_path();
        assert!(result.is_ok());
        assert_eq!(result.unwrap(), test_model);

        // Test with models_dir
        config.cm_model_path = None;
        config.models_dir = Some(temp_dir.clone());
        config.cm_model = PathBuf::from("test.cm");

        let result = config.get_model_path();
        assert!(result.is_ok());
        assert_eq!(result.unwrap(), test_model);

        // Cleanup
        let _ = fs::remove_file(&test_model);
        let _ = fs::remove_dir(&temp_dir);
    }

    #[test]
    fn test_model_not_found_error() {
        let mut config = TrnaScanConfig::new(SearchMode::Eukaryotic);
        config.cm_model = PathBuf::from("nonexistent_model.cm");

        let result = config.get_model_path();
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("not found"));
    }

    #[test]
    fn test_find_default_model() {
        // This test will pass if any default model exists in standard locations
        // or fail silently if none are found (which is expected in test environments)
        let result = find_default_model();
        // Just ensure the function runs without panicking
        if let Some(path) = result {
            assert!(path.exists());
        }
    }
}
