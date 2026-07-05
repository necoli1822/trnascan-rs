//! tRNAscan-SE - tRNA detection in genomic sequences
//!
//! This is the main CLI binary for tRNAscan-SE, providing a complete
//! interface for tRNA detection using the combined tRNAscan/EufindtRNA
//! first-pass with Cove/Infernal covariance model verification.

use clap::{Parser, ArgAction, ValueEnum};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::process;
use std::time::Instant;

use trnascan_rs::core::TrnaScanner;
use trnascan_rs::log::LogFile;
use trnascan_rs::options::{Options, SearchMode, CMMode, MitoModel};
use trnascan_rs::squid::SeqFileReader;
use trnascan_rs::stats::ScanStats;

/// Search mode for organism type
#[derive(Debug, Clone, Copy, ValueEnum, Default)]
pub enum OrgMode {
    /// Eukaryotic search mode (default)
    #[default]
    Eukaryotic,
    /// Bacterial search mode
    Bacterial,
    /// Archaeal search mode
    Archaeal,
    /// General search mode (all domains)
    General,
    /// Organellar search mode
    Organellar,
}

impl From<OrgMode> for SearchMode {
    fn from(mode: OrgMode) -> Self {
        match mode {
            OrgMode::Eukaryotic => SearchMode::Eukaryotic,
            OrgMode::Bacterial => SearchMode::Bacterial,
            OrgMode::Archaeal => SearchMode::Archaeal,
            OrgMode::General => SearchMode::General,
            OrgMode::Organellar => SearchMode::Organellar,
        }
    }
}

#[derive(Parser)]
#[command(name = "trnascan-rs")]
#[command(version)]
#[command(about = "tRNAscan-SE: A program for tRNA detection and analysis")]
#[command(long_about = "
tRNAscan-SE identifies transfer RNA genes in genomic DNA or RNA sequences.
It combines two tRNA detection algorithms: tRNAscan and EufindtRNA, with
covariance model analysis (Cove or Infernal) for verification and refinement.

USAGE:
    trnascan-rs [OPTIONS] <INPUT>...

SEARCH MODES:
    By default, tRNAscan-SE uses eukaryotic search mode. Use -B for bacterial,
    -A for archaeal, -G for general (all domains), or -O for organellar mode.

EXAMPLES:
    trnascan-rs genome.fa                    # Eukaryotic mode (default)
    trnascan-rs -B bacteria.fa -o output.txt # Bacterial mode with output file
    trnascan-rs -A archaea.fa -f struct.ss   # Archaeal mode with structure output
    trnascan-rs -G -b output.bed input.fa    # General mode with BED output

For more information, see: https://github.com/UCSC-LoweLab/tRNAscan-SE
")]
struct Args {
    /// Input sequence file(s) (FASTA, GenBank, EMBL formats supported)
    #[arg(required = true)]
    input: Vec<PathBuf>,

    // === Search Mode Flags ===

    /// General search mode - combines models from all 3 phylogenetic domains
    #[arg(short = 'G', long = "general", action = ArgAction::SetTrue,
          conflicts_with_all = ["bacterial", "archaeal", "organellar", "mitochondrial"])]
    general: bool,

    /// Bacterial search mode - uses bacterial-specific CM and parameters
    #[arg(short = 'B', long = "bacterial", action = ArgAction::SetTrue,
          conflicts_with_all = ["general", "archaeal", "organellar", "mitochondrial"])]
    bacterial: bool,

    /// Archaeal search mode - uses archaeal-specific CM and parameters
    #[arg(short = 'A', long = "archaeal", action = ArgAction::SetTrue,
          conflicts_with_all = ["general", "bacterial", "organellar", "mitochondrial"])]
    archaeal: bool,

    /// Organellar (chloroplast/mitochondrial) search mode
    #[arg(short = 'O', long = "organellar", action = ArgAction::SetTrue,
          conflicts_with_all = ["general", "bacterial", "archaeal", "mitochondrial"])]
    organellar: bool,

    /// Mitochondrial search mode with specified model (mammal, vert)
    #[arg(short = 'M', long = "mito", value_name = "MODEL",
          conflicts_with_all = ["general", "bacterial", "archaeal", "organellar"])]
    mitochondrial: Option<String>,

    // === Output Files ===

    /// Main output file (default: stdout)
    #[arg(short = 'o', long = "output", value_name = "FILE")]
    output: Option<PathBuf>,

    /// Secondary structure output file
    #[arg(short = 'f', long = "struct", value_name = "FILE")]
    ss_file: Option<PathBuf>,

    /// BED format output file
    #[arg(short = 'b', long = "bed", value_name = "FILE")]
    bed_file: Option<PathBuf>,

    /// Statistics output file
    #[arg(short = 'm', long = "stats", value_name = "FILE")]
    stats_file: Option<PathBuf>,

    /// Isotype-specific scores output file
    #[arg(short = 'j', long = "iso", value_name = "FILE")]
    iso_file: Option<PathBuf>,

    /// Log file for progress and debugging
    #[arg(short = 'y', long = "log", value_name = "FILE")]
    log_file: Option<PathBuf>,

    /// FASTA output file for predicted tRNA sequences
    #[arg(short = 'a', long = "fasta-out", value_name = "FILE")]
    fasta_out: Option<PathBuf>,

    /// GFF3 format output file
    #[arg(long = "gff", value_name = "FILE")]
    gff_file: Option<PathBuf>,

    // === Score Thresholds ===

    /// Covariance model score cutoff (default: 20.0)
    #[arg(short = 'X', long = "cutoff", default_value = "20.0")]
    cutoff: f64,

    /// Infernal first-pass score cutoff (for -I mode)
    #[arg(long = "inf-cutoff", default_value = "10.0")]
    infernal_cutoff: f64,

    // === Search Options ===

    /// Show pseudogenes (tRNAs with score below cutoff)
    #[arg(long = "pseudo", action = ArgAction::SetTrue)]
    pseudo: bool,

    /// Search top strand only
    #[arg(short = 'H', long = "top-strand", action = ArgAction::SetTrue)]
    top_strand: bool,

    /// Disable pseudogene checking
    #[arg(long = "nopseudo", action = ArgAction::SetTrue)]
    no_pseudo: bool,

    /// First-pass padding around hits (default: 10)
    #[arg(short = 'L', long = "pad", default_value = "10")]
    padding: usize,

    /// Maximum intron + variable loop length (default: 200)
    #[arg(long = "max-intron", default_value = "200")]
    max_intron: usize,

    // === First-pass Options ===

    /// Use tRNAscan only for first pass (no EufindtRNA)
    #[arg(short = 'T', long = "tscan-only", action = ArgAction::SetTrue)]
    tscan_only: bool,

    /// Use EufindtRNA only for first pass (no tRNAscan)
    #[arg(short = 'E', long = "eufind-only", action = ArgAction::SetTrue)]
    eufind_only: bool,

    /// Use Infernal for first pass (instead of tRNAscan/EufindtRNA)
    #[arg(short = 'I', long = "infernal-fp", action = ArgAction::SetTrue)]
    infernal_fp: bool,

    /// Relaxed tRNAscan parameters (more sensitive, more false positives)
    #[arg(short = 'r', long = "relaxed", action = ArgAction::SetTrue)]
    relaxed: bool,

    // === Second-pass Options ===

    /// Use Infernal instead of Cove for second pass
    #[arg(short = 'C', long = "infernal", action = ArgAction::SetTrue)]
    use_infernal: bool,

    /// Skip covariance model search (first-pass only)
    #[arg(short = 'N', long = "no-cove", action = ArgAction::SetTrue)]
    no_cove: bool,

    /// Disable isotype-specific CM scanning
    #[arg(long = "no-isotype", action = ArgAction::SetTrue)]
    no_isotype: bool,

    // === Output Format Options ===

    /// Quiet mode - suppress progress display
    #[arg(short = 'q', long = "quiet", action = ArgAction::SetTrue)]
    quiet: bool,

    /// Brief output format (no headers)
    #[arg(long = "brief", action = ArgAction::SetTrue)]
    brief: bool,

    /// Show progress during search
    #[arg(long = "progress", action = ArgAction::SetTrue)]
    progress: bool,

    /// Verbose output
    #[arg(short = 'v', long = "verbose", action = ArgAction::SetTrue)]
    verbose: bool,

    /// Output codon instead of anticodon
    #[arg(long = "codon", action = ArgAction::SetTrue)]
    output_codon: bool,

    /// Show detailed prediction info
    #[arg(short = 'D', long = "detail", action = ArgAction::SetTrue)]
    detail: bool,

    // === Model Paths ===

    /// Path to main covariance model file
    #[arg(short = 'c', long = "cm", value_name = "FILE")]
    cm_model: Option<PathBuf>,

    /// Path to models directory
    #[arg(long = "models-dir", value_name = "DIR")]
    models_dir: Option<PathBuf>,

    /// Use alternate genetic code from file
    #[arg(short = 'g', long = "gc", value_name = "FILE")]
    gc_file: Option<PathBuf>,

    // === Threading ===

    /// Worker threads for the parallel CM search (0 = auto: min(cores, 8)).
    /// Peak memory scales ~linearly with this (a Phase-II alignment DP matrix per
    /// worker); raise for more speed on large machines at higher memory cost.
    #[arg(long = "threads", visible_alias = "thread", short = 't', default_value = "0")]
    threads: usize,

    // === Legacy/Debug Options ===

    /// Read first-pass results from previous run
    #[arg(long = "prev-run", value_name = "FILE")]
    prev_run: Option<PathBuf>,

    /// Save first-pass results to file
    #[arg(long = "save-fp", value_name = "FILE")]
    save_firstpass: Option<PathBuf>,

    /// Save false positives to file
    #[arg(long = "falsepos", value_name = "FILE")]
    falsepos_file: Option<PathBuf>,
}

/// Build Options struct from CLI arguments
fn build_options(args: &Args) -> Options {
    let mut opts = Options::new();

    // Set search mode based on flags
    let search_mode = if args.general {
        SearchMode::General
    } else if args.bacterial {
        SearchMode::Bacterial
    } else if args.archaeal {
        SearchMode::Archaeal
    } else if args.organellar {
        SearchMode::Organellar
    } else if args.mitochondrial.is_some() {
        SearchMode::Mitochondrial
    } else {
        SearchMode::Eukaryotic
    };
    opts.set_search_mode(search_mode);

    // Set mito model if specified
    if let Some(ref model) = args.mitochondrial {
        let mito_model = match model.to_lowercase().as_str() {
            "mammal" | "mammalian" => MitoModel::Mammal,
            "vert" | "vertebrate" => MitoModel::Vertebrate,
            _ => {
                eprintln!("Warning: Unknown mito model '{}', using default", model);
                MitoModel::None
            }
        };
        opts.set_mito_model(mito_model);
    }

    // Set CM mode
    if args.no_cove {
        opts.set_cm_mode(CMMode::None);
    } else if args.use_infernal {
        opts.set_cm_mode(CMMode::Infernal);
    } else {
        opts.set_cm_mode(CMMode::Cove);
    }

    // First-pass mode configuration
    if args.tscan_only {
        opts.tscan_mode = true;
        opts.eufind_mode = false;
    } else if args.eufind_only {
        opts.tscan_mode = false;
        opts.eufind_mode = true;
    } else if args.infernal_fp {
        opts.tscan_mode = false;
        opts.eufind_mode = false;
        opts.infernal_fp = true;
    } else {
        // Default: use both tRNAscan and EufindtRNA
        opts.tscan_mode = true;
        opts.eufind_mode = true;
    }

    // Parameters
    opts.strict_params = !args.relaxed;
    opts.set_padding(args.padding);
    opts.set_max_intron_len(args.max_intron);

    // Output options
    opts.quiet_mode = args.quiet;
    opts.brief_output = args.brief;
    opts.display_progress = args.progress;
    opts.output_codon = args.output_codon;
    opts.detail = args.detail;
    opts.no_isotype = args.no_isotype;

    // Set output files
    if let Some(ref path) = args.output {
        opts.set_out_file(&path.to_string_lossy());
    }
    if let Some(ref path) = args.ss_file {
        opts.set_struct_file(&path.to_string_lossy());
    }
    if let Some(ref path) = args.bed_file {
        opts.set_bed_file(&path.to_string_lossy());
    }
    if let Some(ref path) = args.gff_file {
        opts.set_gff_file(&path.to_string_lossy());
    }
    if let Some(ref path) = args.stats_file {
        opts.set_stats_file(&path.to_string_lossy());
    }
    if let Some(ref path) = args.log_file {
        opts.set_log_file(&path.to_string_lossy());
    }
    if let Some(ref path) = args.iso_file {
        opts.set_isotype_file(&path.to_string_lossy());
    }
    if let Some(ref path) = args.fasta_out {
        opts.set_output_fasta_file(&path.to_string_lossy());
    }
    if let Some(ref path) = args.gc_file {
        opts.set_gc_file(&path.to_string_lossy());
    }
    if let Some(ref path) = args.save_firstpass {
        opts.set_firstpass_result_file(&path.to_string_lossy());
    }
    if let Some(ref path) = args.falsepos_file {
        opts.set_falsepos_file(&path.to_string_lossy());
    }
    if let Some(ref path) = args.prev_run {
        opts.use_prev_ts_run = true;
        opts.firstpass_result_file = path.to_string_lossy().to_string();
    }

    opts
}

/// Display credits and version info
fn display_credits(verbose: bool) {
    if verbose {
        eprintln!("trnascan-rs v{} (Rust port of tRNAscan-SE 2.0)", env!("CARGO_PKG_VERSION"));
        eprintln!("Based on tRNAscan-SE by Todd Lowe and Sean Eddy");
        eprintln!("Rust port integrating Infernal/Easel Rust crates");
        eprintln!();
    }
}

/// Format name for detected sequence format
fn format_name(format: i32) -> &'static str {
    match format {
        1 => "FASTA/Pearson",
        2 => "EMBL",
        3 => "GenBank",
        4 => "PIR",
        5 => "IG",
        6 => "NBRF",
        7 => "Extended FASTA",
        _ => "Unknown"
    }
}

/// Main entry point
fn main() {
    let args = Args::parse();
    let start_time = Instant::now();

    // Size the global rayon pool. Default (0) caps at 8 threads to keep peak
    // memory bounded — the parallel Phase-II verify holds a non-banded alignment
    // DP matrix per worker (~10-15 MB), so 128 idle cores would otherwise balloon
    // RSS to >1.5 GB on a large genome for little extra speed. Explicit --threads
    // overrides. Output is byte-identical regardless of thread count.
    let nthreads = if args.threads > 0 {
        args.threads
    } else {
        std::thread::available_parallelism()
            .map(|n| n.get().min(8))
            .unwrap_or(1)
    };
    let _ = rayon::ThreadPoolBuilder::new()
        .num_threads(nthreads)
        .build_global();

    // Validate input files exist
    for input_path in &args.input {
        if !input_path.exists() {
            eprintln!("Error: Input file '{}' not found", input_path.display());
            process::exit(1);
        }
    }

    // Build options from args
    let opts = build_options(&args);

    // Display credits if not quiet
    display_credits(args.verbose && !args.quiet);

    // Initialize log file if specified
    let mut log = if let Some(ref log_path) = args.log_file {
        match LogFile::new(Some(log_path.as_path()), args.quiet) {
            Ok(mut l) => {
                let _ = l.initialize("tRNAscan-SE", Some(&format!("{:?}", std::env::args().collect::<Vec<_>>())));
                Some(l)
            }
            Err(e) => {
                eprintln!("Warning: Could not open log file: {}", e);
                None
            }
        }
    } else {
        None
    };

    // Initialize statistics
    let mut stats = ScanStats::new();
    stats.start_fp_timer();

    // Determine search mode character
    let mode_char = if args.bacterial {
        'B'
    } else if args.archaeal {
        'A'
    } else if args.organellar {
        'O'
    } else if args.general {
        'G'
    } else if args.mitochondrial.is_some() {
        'M'
    } else {
        'E' // Eukaryotic (default)
    };

    let mode_name = match mode_char {
        'B' => "Bacterial",
        'A' => "Archaeal",
        'O' => "Organellar",
        'G' => "General",
        'M' => "Mitochondrial",
        _ => "Eukaryotic",
    };

    // Display run options if verbose
    if args.verbose && !args.quiet {
        eprintln!("Search mode: {} ({})", mode_name, mode_char);
        eprintln!("Score cutoff: {}", args.cutoff);
        if let Some(ref models_dir) = args.models_dir {
            eprintln!("Models directory: {}", models_dir.display());
        }
        eprintln!("Input files: {:?}", args.input);
        eprintln!();
    }

    // Initialize scanner with appropriate models directory
    let mut scanner = if let Some(ref models_dir) = args.models_dir {
        // Use explicit models directory
        match TrnaScanner::with_models_dir(mode_char, args.cutoff, models_dir) {
            Ok(s) => s,
            Err(e) => {
                eprintln!("Error initializing scanner with models from '{}': {}", models_dir.display(), e);
                if let Some(ref mut l) = log {
                    let _ = l.error(&format!("Failed to load models: {}", e));
                }
                process::exit(1);
            }
        }
    } else if let Some(ref cm_path) = args.cm_model {
        // Legacy: derive models_dir from CM file path
        match TrnaScanner::with_model_path(mode_char, args.cutoff, cm_path) {
            Ok(s) => s,
            Err(e) => {
                eprintln!("Error initializing scanner with model '{}': {}", cm_path.display(), e);
                if let Some(ref mut l) = log {
                    let _ = l.error(&format!("Failed to load model: {}", e));
                }
                process::exit(1);
            }
        }
    } else {
        // Default: use built-in models path
        match TrnaScanner::new(mode_char, args.cutoff) {
            Ok(s) => s,
            Err(e) => {
                eprintln!("Error initializing scanner: {}", e);
                if let Some(ref mut l) = log {
                    let _ = l.error(&format!("Failed to initialize scanner: {}", e));
                }
                process::exit(1);
            }
        }
    };

    // Configure scanner
    scanner.set_quiet(args.quiet);
    scanner.set_verbose(args.verbose);
    scanner.set_show_pseudogenes(args.pseudo);
    // `-H` → HMM Score / 2'Str Score columns; `--detail` → Isotype columns + note.
    scanner.set_get_hmm_score(args.top_strand);
    scanner.set_detail(args.detail);
    scanner.set_no_isotype(args.no_isotype);

    // Process each input file
    let mut total_seq_count = 0;
    let mut total_base_count = 0;

    for input_path in &args.input {
        if args.verbose && !args.quiet {
            eprintln!("Processing: {}", input_path.display());
        }
        if let Some(ref mut l) = log {
            let _ = l.status(&format!("Processing: {}", input_path.display()));
        }

        // Open input file
        let mut seq_reader = match SeqFileReader::open(input_path) {
            Ok(r) => r,
            Err(e) => {
                eprintln!("Error opening input file '{}': {}", input_path.display(), e);
                if let Some(ref mut l) = log {
                    let _ = l.error(&format!("Failed to open file: {}", e));
                }
                continue;
            }
        };

        if args.verbose && !args.quiet {
            eprintln!("Input format detected: {}", format_name(seq_reader.format()));
        }

        // Process sequences
        let mut file_seq_count = 0;
        loop {
            match seq_reader.read_seq() {
                Ok(Some((seq, sqinfo))) => {
                    file_seq_count += 1;
                    total_seq_count += 1;
                    total_base_count += seq.len();

                    stats.inc_numscanned();
                    stats.increment_first_pass_base_ct(seq.len());

                    if args.progress && !args.quiet {
                        eprint!("\rScanning sequence {}: {} ({} bp)     ",
                            total_seq_count, sqinfo.name, sqinfo.len);
                    }

                    if let Err(e) = scanner.scan_sequence(&seq, &sqinfo) {
                        eprintln!("\nWarning: Error scanning {}: {}", sqinfo.name, e);
                        if let Some(ref mut l) = log {
                            let _ = l.warning(&format!("Scan error for {}: {}", sqinfo.name, e));
                        }
                    }

                    // Update stats if tRNAs found
                    if scanner.result_count() > stats.trnatotal() {
                        let new_trnas = scanner.result_count() - stats.trnatotal();
                        stats.increment_trnatotal(new_trnas);
                        stats.inc_seqs_hit();
                    }
                }
                Ok(None) => break,  // EOF
                Err(e) => {
                    eprintln!("Error reading sequence from '{}': {}", input_path.display(), e);
                    if let Some(ref mut l) = log {
                        let _ = l.error(&format!("Read error: {}", e));
                    }
                    break;
                }
            }
        }

        if args.verbose && !args.quiet {
            eprintln!("\nProcessed {} sequences from {}", file_seq_count, input_path.display());
        }
    }

    // Clear progress line
    if args.progress && !args.quiet {
        eprintln!();
    }

    // End timers
    stats.end_fp_timer();
    stats.start_sp_timer();
    stats.end_sp_timer();

    // Write outputs
    if let Err(e) = write_outputs(&scanner, &args, &opts, &stats) {
        eprintln!("Error writing output: {}", e);
        if let Some(ref mut l) = log {
            let _ = l.error(&format!("Output error: {}", e));
        }
        process::exit(1);
    }

    // Summary
    let elapsed = start_time.elapsed();
    if !args.quiet {
        eprintln!();
        eprintln!("Scan complete.");
        eprintln!("Sequences scanned: {}", total_seq_count);
        eprintln!("Total bases: {}", total_base_count);
        eprintln!("tRNAs found: {}", scanner.result_count());
        eprintln!("Time elapsed: {:.2}s", elapsed.as_secs_f64());
    }

    if let Some(ref mut l) = log {
        let _ = l.status(&format!("Scan complete. Found {} tRNAs in {} sequences",
            scanner.result_count(), total_seq_count));
        let _ = l.finish();
    }
}

/// Write all output files
fn write_outputs(
    scanner: &TrnaScanner,
    args: &Args,
    _opts: &Options,
    stats: &ScanStats,
) -> std::io::Result<()> {
    // Main output
    if let Some(ref path) = args.output {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        if scanner.uses_faithful() {
            scanner.write_faithful_out(&mut writer, args.brief)?;
        } else {
            if !args.brief {
                write_header(&mut writer)?;
            }
            scanner.write_results(&mut writer)?;
        }
    } else if !args.quiet {
        let stdout = std::io::stdout();
        let mut writer = stdout.lock();
        if scanner.uses_faithful() {
            scanner.write_faithful_out(&mut writer, args.brief)?;
        } else {
            if !args.brief {
                write_header(&mut writer)?;
            }
            scanner.write_results(&mut writer)?;
        }
    }

    // Secondary structure output (-f)
    if let Some(ref path) = args.ss_file {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        if scanner.uses_faithful() {
            scanner.write_faithful_struct(&mut writer)?;
        } else {
            scanner.write_secondary_structures(&mut writer)?;
        }
    }

    // Statistics output (-m). The banner / first-pass / second-pass blocks carry
    // timestamps + CPU times and cannot byte-match C; only the summary tail is
    // deterministic (faithful path uses the trna_results-driven writer).
    if let Some(ref path) = args.stats_file {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        stats.save_firstpass_stats(&mut writer)?;
        if scanner.uses_faithful() {
            scanner.write_faithful_stats(&mut writer)?;
        } else {
            scanner.write_statistics(&mut writer)?;
        }
    }

    // BED format output (-b)
    if let Some(ref path) = args.bed_file {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        if scanner.uses_faithful() {
            scanner.write_faithful_bed(&mut writer)?;
        } else {
            scanner.write_bed_format(&mut writer)?;
        }
    }

    // GFF3 output (--gff)
    if let Some(ref path) = args.gff_file {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        if scanner.uses_faithful() {
            scanner.write_faithful_gff(&mut writer)?;
        }
    }

    // FASTA output (-a)
    if let Some(ref path) = args.fasta_out {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        if scanner.uses_faithful() {
            scanner.write_faithful_fasta(&mut writer)?;
        }
    }

    // Isotype model output
    if let Some(ref path) = args.iso_file {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        scanner.write_isotype_models(&mut writer)?;
    }

    Ok(())
}

/// Write output header
fn write_header<W: Write>(writer: &mut W) -> std::io::Result<()> {
    writeln!(writer, "Sequence\t\ttRNA\tBounds\ttRNA\tAnti\tIntron Bounds\tInf\tHMM\t2'Str\tHit\tIsotype\tIsotype\t      ")?;
    writeln!(writer, "Name    \ttRNA #\tBegin\tEnd\tType\tCodon\tBegin\tEnd\tScore\tScore\tScore\tOrigin\tCM\tScore\tNote")?;
    writeln!(writer, "--------\t------\t-----\t------\t----\t-----\t-----\t----\t------\t-----\t-----\t------\t-------\t-------\t------")?;
    Ok(())
}
