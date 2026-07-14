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
#[command(version = env!("CARGO_PKG_VERSION"))]
#[command(about = "trnascan-rs: A program for tRNA detection and analysis")]
#[command(long_about = "
trnascan-rs identifies transfer RNA genes in genomic DNA or RNA sequences.
It combines two tRNA detection algorithms: tRNAscan and EufindtRNA, with
covariance model analysis (Cove or Infernal) for verification and refinement.

USAGE:
    trnascan-rs [OPTIONS] <INPUT>...

SEARCH MODES:
    By default, trnascan-rs uses eukaryotic search mode. Use -B for bacterial,
    -A for archaeal, -G for general (all domains), or -O for organellar mode.

EXAMPLES:
    trnascan-rs genome.fa                    # Eukaryotic mode (default)
    trnascan-rs -B bacteria.fa -o output.txt # Bacterial mode with output file
    trnascan-rs -A archaea.fa -f struct.ss   # Archaeal mode with structure output
    trnascan-rs -G -b output.bed input.fa    # General mode with BED output

For more information, see: https://github.com/necoli1822/trnascan-rs
")]
struct Args {
    /// Input sequence file(s) (FASTA, GenBank, EMBL formats supported)
    #[arg(required = true)]
    input: Vec<PathBuf>,

    // === Search Mode Options ===

    // C tRNAscan-SE GetOptions: "euk|E" — search for eukaryotic tRNAs (default)
    /// search for eukaryotic tRNAs (default)
    #[arg(short = 'E', long = "euk", action = ArgAction::SetTrue,
          conflicts_with_all = ["general", "bacterial", "archaeal", "organellar", "mitochondrial"])]
    eukaryotic: bool,

    // C tRNAscan-SE GetOptions: "bact|B" — search for bacterial tRNAs
    /// search for bacterial tRNAs
    #[arg(short = 'B', long = "bact", action = ArgAction::SetTrue,
          conflicts_with_all = ["eukaryotic", "general", "archaeal", "organellar", "mitochondrial"])]
    bacterial: bool,

    // C tRNAscan-SE GetOptions: "arch|A" — search for archaeal tRNAs
    /// search for archaeal tRNAs
    #[arg(short = 'A', long = "arch", action = ArgAction::SetTrue,
          conflicts_with_all = ["eukaryotic", "general", "bacterial", "organellar", "mitochondrial"])]
    archaeal: bool,

    // C tRNAscan-SE GetOptions: "mito|M=s" — search for mitochondrial tRNAs (mammal, vert)
    /// search for mitochondrial tRNAs (options: mammal, vert)
    #[arg(short = 'M', long = "mito", value_name = "model",
          conflicts_with_all = ["eukaryotic", "general", "bacterial", "archaeal", "organellar"])]
    mitochondrial: Option<String>,

    // C tRNAscan-SE GetOptions: "organ|O" — search for other organellar tRNAs
    /// search for other organellar tRNAs
    #[arg(short = 'O', long = "organ", action = ArgAction::SetTrue,
          conflicts_with_all = ["eukaryotic", "general", "bacterial", "archaeal", "mitochondrial"])]
    organellar: bool,

    // C tRNAscan-SE GetOptions: "general|G" — use general tRNA model (all 3 domains)
    /// use general tRNA model (cytosolic tRNAs from all 3 domains included)
    #[arg(short = 'G', long = "general", action = ArgAction::SetTrue,
          conflicts_with_all = ["eukaryotic", "bacterial", "archaeal", "organellar", "mitochondrial"])]
    general: bool,

    // C tRNAscan-SE -h: --mt <model> : use mito tRNA models for cytosolic/mito determination
    /// use mito tRNA models for cytosolic/mito determination
    #[arg(long = "mt", value_name = "model")]
    mt: Option<String>,

    // C tRNAscan-SE GetOptions: "inf|I" — search using Infernal
    /// search using Infernal
    #[arg(short = 'I', long = "inf", action = ArgAction::SetTrue)]
    infernal: bool,

    // C tRNAscan-SE -h: --max : maximum sensitivity (Infernal, no hmm filter)
    /// maximum sensitivity mode - search using Infernal without hmm filter (very slow)
    #[arg(long = "max", action = ArgAction::SetTrue)]
    max_sensitivity: bool,

    // C tRNAscan-SE -h: --mid : fast scan mode - Infernal with mid-level hmm filter
    /// fast scan mode - search using Infernal with mid-level strictness of hmm filter
    #[arg(long = "mid", action = ArgAction::SetTrue)]
    mid: bool,

    // C tRNAscan-SE GetOptions: "legacy|L" — search using the legacy method (tRNAscan, EufindtRNA, and COVE)
    /// search using the legacy method (tRNAscan, EufindtRNA, and COVE)
    #[arg(short = 'L', long = "legacy", action = ArgAction::SetTrue)]
    legacy: bool,

    // C tRNAscan-SE -h: -C --cove : search using COVE analysis only (legacy, extremely slow)
    /// search using COVE analysis only (legacy, extremely slow)
    #[arg(short = 'C', long = "cove", action = ArgAction::SetTrue)]
    cove: bool,

    // C tRNAscan-SE -h: -H --breakdown : show breakdown of primary/secondary structure scores
    /// show breakdown of primary and secondary structure components to CM bit scores
    #[arg(short = 'H', long = "breakdown", action = ArgAction::SetTrue)]
    breakdown: bool,

    // C tRNAscan-SE -h: -D --nopseudo : disable pseudogene checking
    /// disable pseudogene checking
    #[arg(short = 'D', long = "nopseudo", action = ArgAction::SetTrue)]
    nopseudo: bool,

    // === Output options ===

    // C tRNAscan-SE -h: -o --output <file> : save final results in <file>
    /// save final results in <file>
    #[arg(short = 'o', long = "output", value_name = "file")]
    output: Option<PathBuf>,

    // C tRNAscan-SE -h: -f --struct <file> : save tRNA secondary structures to <file>
    /// save tRNA secondary structures to <file>
    #[arg(short = 'f', long = "struct", value_name = "file")]
    ss_file: Option<PathBuf>,

    // C tRNAscan-SE -h: -s --isospecific <file> : save isotype-specific model results
    /// save results using isotype-specific models in <file>
    #[arg(short = 's', long = "isospecific", value_name = "file")]
    isospecific_file: Option<PathBuf>,

    // C tRNAscan-SE -h: -m --stats <file> : save statistics summary for run
    /// save statistics summary for run in <file>
    #[arg(short = 'm', long = "stats", value_name = "file")]
    stats_file: Option<PathBuf>,

    // C tRNAscan-SE -h: -b --bed <file> : save results in BED file format
    /// save results in BED file format of <file>
    #[arg(short = 'b', long = "bed", value_name = "file")]
    bed_file: Option<PathBuf>,

    // C tRNAscan-SE -h: -j --gff <file> : save results in GFF3 file format
    /// save results in GFF3 file format of <file>
    #[arg(short = 'j', long = "gff", value_name = "file")]
    gff_file: Option<PathBuf>,

    // C tRNAscan-SE -h: -a --fasta <file> : save predicted tRNA sequences in FASTA format
    /// save predicted tRNA sequences in FASTA file format of <file>
    #[arg(short = 'a', long = "fasta", value_name = "file")]
    fasta: Option<PathBuf>,

    // C tRNAscan-SE -h: -l --log <file> : save log of program progress
    /// save log of program progress in <file>
    #[arg(short = 'l', long = "log", value_name = "file")]
    log_file: Option<PathBuf>,

    // C tRNAscan-SE -h: --detail : display prediction outputs in detailed view
    /// display prediction outputs in detailed view
    #[arg(long = "detail", action = ArgAction::SetTrue)]
    detail: bool,

    // C tRNAscan-SE -h: --brief : brief output format (no column headers)
    /// brief output format (no column headers)
    #[arg(long = "brief", action = ArgAction::SetTrue)]
    brief: bool,

    // C tRNAscan-SE GetOptions: "acedb" — ACeDB output format (tRNAscan-SE:881
    // `$opts->ace_output(1)`).
    /// use ACeDB output format instead of the default tabular output
    #[arg(long = "acedb", action = ArgAction::SetTrue)]
    acedb: bool,

    // C tRNAscan-SE -h: -p --prefix <label> : use <label> prefix for default output file names
    /// use <label> prefix for all default output file names
    #[arg(short = 'p', long = "prefix", value_name = "label")]
    prefix: Option<String>,

    // C tRNAscan-SE -h: -d --progress : display program progress messages
    /// display program progress messages
    #[arg(short = 'd', long = "progress", action = ArgAction::SetTrue)]
    progress: bool,

    // C tRNAscan-SE -h: -q --quiet : quiet mode (credits & run option selections suppressed)
    /// quiet mode (credits & run option selections suppressed)
    #[arg(short = 'q', long = "quiet", action = ArgAction::SetTrue)]
    quiet: bool,

    // C tRNAscan-SE -h: -y --hitsrc : show origin of hits (Ts/Eu/Bo/Inf)
    /// show origin of hits (Ts=tRNAscan 1.4, Eu=EufindtRNA, Bo=Both, Inf=Infernal)
    #[arg(short = 'y', long = "hitsrc", action = ArgAction::SetTrue)]
    hitsrc: bool,

    // === Specify Alternate Cutoffs / Data Files ===

    // C tRNAscan-SE -h: -X --score <score> : cutoff score (bits) for reporting tRNAs (default=20)
    /// set cutoff score (in bits) for reporting tRNAs (default=20)
    #[arg(short = 'X', long = "score", value_name = "score", default_value = "20.0")]
    cutoff: f64,

    // C tRNAscan-SE -h: -g --gencode <file> : use alternate genetic codes for tRNA type
    /// use alternate genetic codes specified in <file> for determining tRNA type
    #[arg(short = 'g', long = "gencode", value_name = "file")]
    gencode_file: Option<PathBuf>,

    // C tRNAscan-SE GetOptions: "isocm|S=s" — toggle isotype-specific CM scanning
    // (on|off), overriding the per-search-mode default (tRNAscan-SE:955-1349).
    /// turn isotype-specific model scanning on or off (values: on, off)
    #[arg(short = 'S', long = "isocm", value_name = "on|off")]
    isocm: Option<String>,

    // C tRNAscan-SE -h: -z --pad <number> : padding when passing first-pass bounds to CM (default=8)
    /// use <number> nucleotides padding when passing first-pass tRNA bounds to CM analysis (default=8)
    #[arg(short = 'z', long = "pad", value_name = "number", default_value = "8")]
    padding: usize,

    // C tRNAscan-SE GetOptions: "len=i" : max length of tRNA intron+variable region (legacy, default=116)
    /// set max length of tRNA intron+variable region for legacy search mode (default=116bp)
    #[arg(long = "len", value_name = "length", default_value = "116")]
    max_intron: usize,

    // === Misc Options ===

    // C tRNAscan-SE -h: -c --conf <file> : tRNAscan-SE configuration file
    /// tRNAscan-SE configuration file (default: tRNAscan-SE.conf)
    #[arg(short = 'c', long = "conf", value_name = "file")]
    conf_file: Option<PathBuf>,

    // C tRNAscan-SE -h: -Q --forceow : do not prompt before overwriting result files
    /// do not prompt user before overwriting pre-existing result files (batch processing)
    #[arg(short = 'Q', long = "forceow", action = ArgAction::SetTrue)]
    forceow: bool,

    // C tRNAscan-SE -h: --match <EXPR> : search only sequences with names matching <EXPR>
    /// search only sequences with names matching <EXPR> string
    #[arg(long = "match", value_name = "EXPR")]
    match_expr: Option<String>,

    // C tRNAscan-SE -h: --search <EXPR> : start search at sequence with name matching <EXPR>
    /// start search at sequence with name matching <EXPR> string and continue to end
    #[arg(long = "search", value_name = "EXPR")]
    search_expr: Option<String>,

    // === Special Advanced Options (for testing & special purposes) ===

    // C tRNAscan-SE GetOptions: "alt|U" — search for tRNAs with alternate models defined in configuration file
    /// search for tRNAs with alternate models defined in configuration file
    #[arg(short = 'U', long = "alt", action = ArgAction::SetTrue)]
    alternate_models: bool,

    // C tRNAscan-SE -h: -t --tscan : search using tRNAscan only (defaults to strict params)
    /// search using tRNAscan only (defaults to strict params)
    #[arg(short = 't', long = "tscan", action = ArgAction::SetTrue)]
    tscan: bool,

    // C tRNAscan-SE -h: --tmode <mode> : explicitly set tRNAscan params (R=relaxed, S=strict)
    /// explicitly set tRNAscan params, where <mode>=R or S
    #[arg(long = "tmode", value_name = "mode")]
    tmode: Option<String>,

    // C tRNAscan-SE -h: -v --verbose <file> : save verbose tRNAscan 1.3 output to <file>
    /// save verbose tRNAscan 1.3 output to <file>
    #[arg(short = 'v', long = "verbose", value_name = "file")]
    verbose_file: Option<PathBuf>,

    // C tRNAscan-SE -h: --nomerge : keep redundant tRNAscan 1.3 hits
    /// keep redundant tRNAscan 1.3 hits (don't filter multiple predictions per tRNA)
    #[arg(long = "nomerge", action = ArgAction::SetTrue)]
    nomerge: bool,

    // C tRNAscan-SE -h: -e --eufind : search using EufindtRNA only
    /// search using Eukaryotic tRNA finder (EufindtRNA) only
    #[arg(short = 'e', long = "eufind", action = ArgAction::SetTrue)]
    eufind: bool,

    // C tRNAscan-SE -h: --emode <mode> : explicitly set EufindtRNA params (R, N, or S)
    /// explicitly set EufindtRNA params, where <mode>=R, N, or S
    #[arg(long = "emode", value_name = "mode")]
    emode: Option<String>,

    // C tRNAscan-SE -h: --iscore <score> : manually set intermediate cutoff for EufindtRNA
    /// manually set "intermediate" cutoff score for EufindtRNA
    #[arg(long = "iscore", value_name = "score")]
    iscore: Option<f64>,

    // C tRNAscan-SE -h: -r --fsres <file> : save first-pass scan results in tabular format
    /// save first-pass scan results from EufindtRNA, tRNAscan, or Infernal hmm in <file>
    #[arg(short = 'r', long = "fsres", value_name = "file")]
    fsres_file: Option<PathBuf>,

    // C tRNAscan-SE -h: -F --falsepos <file> : save first-pass candidates found to be false positives
    /// save first-pass candidate tRNAs later found to be false positives in <file>
    #[arg(short = 'F', long = "falsepos", value_name = "file")]
    falsepos_file: Option<PathBuf>,

    // C tRNAscan-SE -h: --missed <file> : save seqs with no tRNA prediction
    /// save all seqs that do NOT have at least one tRNA prediction (aka "missed" seqs)
    #[arg(long = "missed", value_name = "file")]
    missed_file: Option<PathBuf>,

    // C tRNAscan-SE GetOptions: "w=s" — save tRNAs with odd (uncallable) secondary
    // structure to <file> (tRNAscan-SE:1729 `$opts->save_odd_struct(1)` +
    // `odd_struct_file`).
    /// save tRNAs with odd (uncallable) secondary structures in <file>
    #[arg(short = 'w', value_name = "file")]
    odd_struct_file: Option<PathBuf>,

    // C tRNAscan-SE GetOptions: "Y" — write a <prefix>.pid file recording the run's
    // process id (tRNAscan-SE:1784).
    /// write a <prefix>.pid file recording the process id of this run
    #[arg(short = 'Y', action = ArgAction::SetTrue)]
    write_pid: bool,

    // C tRNAscan-SE -h: --thread <number> : number of threads used for running infernal
    /// number of threads used for running infernal (default is to use available threads)
    #[arg(long = "thread", value_name = "number", default_value = "0")]
    threads: usize,

    // === Rust-only options (no C tRNAscan-SE equivalent; long-only to avoid collisions) ===

    /// [Rust-only] Path to models directory (used by the test harness)
    #[arg(long = "models-dir", value_name = "DIR")]
    models_dir: Option<PathBuf>,

    /// [Rust-only] Path to main covariance model file
    #[arg(long = "cm", value_name = "FILE")]
    cm_model: Option<PathBuf>,

    /// [Rust-only] Infernal first-pass score cutoff
    #[arg(long = "inf-cutoff", default_value = "10.0")]
    infernal_cutoff: f64,

    /// [Rust-only] Read first-pass results from a previous run
    #[arg(long = "prev-run", value_name = "FILE")]
    prev_run: Option<PathBuf>,

    /// [Rust-only] Save first-pass results to file
    #[arg(long = "save-fp", value_name = "FILE")]
    save_firstpass: Option<PathBuf>,

    /// [Rust-only] Disable isotype-specific CM scanning
    #[arg(long = "no-isotype", action = ArgAction::SetTrue)]
    no_isotype: bool,

    /// [Rust-only] Skip covariance model search (first-pass only)
    #[arg(long = "no-cove", action = ArgAction::SetTrue)]
    no_cove: bool,

    /// [Rust-only] Show pseudogenes (tRNAs with score below cutoff)
    #[arg(long = "pseudo", action = ArgAction::SetTrue)]
    pseudo: bool,

    /// [Rust-only] Relaxed tRNAscan parameters
    #[arg(long = "relaxed", action = ArgAction::SetTrue)]
    relaxed: bool,

    // C tRNAscan-SE GetOptions: "codons" — output the codon corresponding to each
    // predicted anticodon (C tRNAscan-SE:916 `$opts->output_codon(1)`).
    /// output a tRNA's corresponding codon in place of its anticodon
    #[arg(long = "codons", action = ArgAction::SetTrue)]
    output_codon: bool,

    /// [Rust-only] Extra diagnostic output on stderr
    #[arg(long = "rs-verbose", action = ArgAction::SetTrue)]
    rs_verbose: bool,
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

    // --mt <model>: mito tRNA models for cytosolic/mito determination.
    if let Some(ref model) = args.mt {
        opts.mt_model = model.clone();
    }

    // Set CM / search method mode.
    // C tRNAscan-SE -h: -C --cove (COVE only), -I (Infernal), --max / --mid (Infernal),
    //                   -L (legacy = tRNAscan+EufindtRNA+COVE), --no-cove (Rust-only).
    if args.no_cove {
        opts.set_cm_mode(CMMode::None);
    } else if args.cove {
        opts.set_cm_mode(CMMode::Cove);
    } else if args.legacy {
        opts.legacy_mode = true;
        opts.set_cm_mode(CMMode::Cove);
    } else if args.infernal || args.max_sensitivity || args.mid {
        opts.set_cm_mode(CMMode::Infernal);
    }

    // Infernal hmm-filter strictness: --max = no filter; --mid = mid-level filter.
    if args.max_sensitivity {
        opts.hmm_filter = false;
    } else if args.mid {
        opts.hmm_filter = true;
    }

    // First-pass mode configuration.
    // C tRNAscan-SE -h: -t --tscan (tRNAscan only), -e --eufind (EufindtRNA only).
    if args.tscan {
        opts.tscan_mode = true;
        opts.eufind_mode = false;
    } else if args.eufind {
        opts.tscan_mode = false;
        opts.eufind_mode = true;
    } else if args.infernal || args.max_sensitivity || args.mid {
        opts.tscan_mode = false;
        opts.eufind_mode = false;
        opts.infernal_fp = true;
    } else {
        // Default: use both tRNAscan and EufindtRNA
        opts.tscan_mode = true;
        opts.eufind_mode = true;
    }

    // tRNAscan / EufindtRNA explicit param modes (parses; strictness stored).
    if let Some(ref m) = args.tmode {
        opts.tscan_strictness = m.clone();
    }
    if let Some(ref m) = args.emode {
        opts.eufind_strictness = m.clone();
    }
    if let Some(score) = args.iscore {
        opts.eufind_intermediate_score = Some(score);
    }

    // -U : alternate models; --nomerge : keep redundant tRNAscan hits; -D : nopseudo.
    opts.use_alternate_models = args.alternate_models;
    opts.nomerge = args.nomerge;
    opts.disable_pseudo = args.nopseudo;

    // -y --hitsrc : show origin of hits.
    opts.save_source = args.hitsrc;

    // -c --conf : configuration file (parses; behavior deferred).
    if let Some(ref path) = args.conf_file {
        opts.conf_file = path.to_string_lossy().to_string();
    }

    // -p --prefix : default-output-file prefix (parses; behavior deferred).
    if let Some(ref label) = args.prefix {
        opts.output_prefix = label.clone();
    }

    // --match / --search : sequence-name filtering.
    if let Some(ref expr) = args.match_expr {
        opts.set_seq_key(expr);
    }
    if let Some(ref expr) = args.search_expr {
        opts.set_seq_key(expr);
        opts.start_at_key = true;
    }

    // -Q --forceow : do not prompt before overwriting result files.
    if args.forceow {
        opts.prompt_for_overwrite = false;
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

    // -S/--isocm <on|off>: toggle isotype-specific CM scanning, overriding the
    // per-mode default (C tRNAscan-SE:955-1349). `off` forces no_isotype; `on` is
    // valid only for euk/bact/arch (already on) and FATAL for -G/-O/-M/-U; any
    // other value is FATAL (tRNAscan-SE:955-957).
    if let Some(ref v) = args.isocm {
        match v.as_str() {
            "off" => opts.no_isotype = true,
            "on" => {
                if args.general
                    || args.organellar
                    || args.mitochondrial.is_some()
                    || args.alternate_models
                {
                    eprintln!(
                        "FATAL: Conflicting search options have been selected. \
                         The selected search mode cannot be combined with --isocm on."
                    );
                    std::process::exit(1);
                }
            }
            _ => {
                eprintln!(
                    "FATAL: Invalid value for --isocm. Please use on or off or \
                     leave out the option for default setting"
                );
                std::process::exit(1);
            }
        }
    }

    // --acedb: ACeDB output format (C tRNAscan-SE:881 `$opts->ace_output(1)`).
    opts.ace_output = args.acedb;

    // -w <file>: save tRNAs with odd (uncallable) secondary structure
    // (C tRNAscan-SE:1729 `$opts->save_odd_struct(1)` + `odd_struct_file`).
    if let Some(ref path) = args.odd_struct_file {
        opts.set_odd_struct_file(&path.to_string_lossy());
    }

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
    if let Some(ref path) = args.isospecific_file {
        opts.set_isotype_file(&path.to_string_lossy());
    }
    if let Some(ref path) = args.fasta {
        opts.set_output_fasta_file(&path.to_string_lossy());
    }
    if let Some(ref path) = args.gencode_file {
        opts.set_gc_file(&path.to_string_lossy());
    }
    // -r --fsres and --save-fp (Rust-only alias) both target the first-pass result file.
    if let Some(ref path) = args.fsres_file {
        opts.set_firstpass_result_file(&path.to_string_lossy());
    }
    if let Some(ref path) = args.save_firstpass {
        opts.set_firstpass_result_file(&path.to_string_lossy());
    }
    if let Some(ref path) = args.falsepos_file {
        opts.set_falsepos_file(&path.to_string_lossy());
    }
    if let Some(ref path) = args.missed_file {
        opts.set_missed_file(&path.to_string_lossy());
    }
    if let Some(ref path) = args.verbose_file {
        opts.set_verbose_file(&path.to_string_lossy());
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
        eprintln!("trnascan-rs v{} (Rust reimplementation)", env!("CARGO_PKG_VERSION"));
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

    // -Y: write a `<prefix>.pid` file recording this run's process id
    // (C tRNAscan-SE:1784-1790: `print TESTF "PID=$$\n"`). `$fafile` is the -p
    // prefix when given, else the first input path.
    if args.write_pid {
        let fafile = args
            .prefix
            .clone()
            .unwrap_or_else(|| args.input[0].to_string_lossy().to_string());
        let _ = std::fs::write(format!("{fafile}.pid"), format!("PID={}\n", process::id()));
    }

    // Build options from args
    let opts = build_options(&args);

    // Display credits if not quiet
    display_credits(args.rs_verbose && !args.quiet);

    // Initialize log file if specified
    let mut log = if let Some(ref log_path) = args.log_file {
        match LogFile::new(Some(log_path.as_path()), args.quiet) {
            Ok(mut l) => {
                let _ = l.initialize("trnascan-rs", Some(&format!("{:?}", std::env::args().collect::<Vec<_>>())));
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
    if args.rs_verbose && !args.quiet {
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
    scanner.set_verbose(args.rs_verbose);
    scanner.set_show_pseudogenes(args.pseudo);
    // `-H`/`--breakdown` → HMM Score / 2'Str Score columns; `--detail` → Isotype columns + note.
    scanner.set_get_hmm_score(args.breakdown);
    scanner.set_detail(args.detail);
    // Effective no_isotype = Rust-only `--no-isotype` OR `-S/--isocm off` (C
    // tRNAscan-SE:1012/1074 set `no_isotype(1)` on `--isocm off` for -B/-A). The
    // per-mode defaults for -G/-O/-M (isotype always off) are handled by those
    // modes' own faithful paths.
    scanner.set_no_isotype(args.no_isotype || args.isocm.as_deref() == Some("off"));
    // `-D`/`--nopseudo`: disable pseudogene checking (C tRNAscan-SE:910
    // `$cm->skip_pseudo_filter(1)`) → suppress the `pseudo` note in the faithful path.
    scanner.set_disable_pseudo(args.nopseudo);
    // `--codons`: report codon instead of anticodon in the faithful `.out`/`-f`.
    scanner.set_output_codon(args.output_codon);
    // `-s`/`--isospecific`: retain the full per-model isotype score vector for `.iso`.
    scanner.set_iso_output(args.isospecific_file.is_some());
    // C tRNAscan-SE:1514-1531: only $opt_inf sets infernal_fp(1). The legacy
    // opts (-e EufindtRNA / -t tRNAscan / -L legacy / -C COVE) leave infernal_fp=0,
    // so they bypass the DEFAULT Infernal first-pass and use the heuristic path.
    scanner.set_legacy_first_pass(args.eufind || args.tscan || args.legacy || args.cove);

    // `-M <model>`: select the mito CM set (vert/mammal) + force cm_cutoff=15.
    // C driver (tRNAscan-SE:1164-1205) rejects anything but "mammal"/"vert".
    if let Some(ref model) = args.mitochondrial {
        let lc = model.to_lowercase();
        match lc.as_str() {
            "mammal" | "mammalian" => scanner.set_mito_model_name("mammal"),
            "vert" | "vertebrate" => scanner.set_mito_model_name("vert"),
            _ => {
                eprintln!(
                    "FATAL: Invalid mitochondrial tRNA option. Only mammal or vert can be used."
                );
                process::exit(1);
            }
        }
    }

    // Process each input file
    let mut total_seq_count = 0;
    let mut total_base_count = 0;

    for input_path in &args.input {
        if args.rs_verbose && !args.quiet {
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

        if args.rs_verbose && !args.quiet {
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

        if args.rs_verbose && !args.quiet {
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

    // Reconcile the `.stats` counters with the faithful pipeline's real first-pass /
    // second-pass accounting (C Stats.pm). During scanning `trnatotal` tracked the
    // running *confirmed* count (for the seqs-hit delta); the `.stats` "tRNAs
    // predicted" / "Candidate tRNAs read" fields are the FIRST-PASS candidate count.
    if scanner.uses_faithful() {
        stats.set_trnatotal(scanner.fp_candidate_ct());
        stats.set_fpass_trna_base_ct(scanner.fp_candidate_bases());
        stats.set_secpass_base_ct(scanner.sp_scanned_bases());
        stats.set_total_secpass_ct(scanner.result_count());
    }

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
            // --acedb routes the main output through the ACeDB-format writer.
            if args.acedb {
                scanner.write_faithful_acedb(&mut writer)?;
            } else {
                scanner.write_faithful_out(&mut writer, args.brief)?;
            }
        } else {
            // C ScanResult.pm:362 writes the column header lazily, INSIDE the
            // per-result loop — zero results => empty `.out` (e.g. `-O` on an
            // organelle genome with no CM hits). Match that here.
            if !args.brief && !scanner.results().is_empty() {
                write_header(&mut writer)?;
            }
            scanner.write_results(&mut writer)?;
        }
    } else if !args.quiet {
        let stdout = std::io::stdout();
        let mut writer = stdout.lock();
        if scanner.uses_faithful() {
            if args.acedb {
                scanner.write_faithful_acedb(&mut writer)?;
            } else {
                scanner.write_faithful_out(&mut writer, args.brief)?;
            }
        } else {
            if !args.brief && !scanner.results().is_empty() {
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

    // -w: odd-struct output (tRNAs with an uncallable anticodon).
    if let Some(ref path) = args.odd_struct_file {
        if scanner.uses_faithful() {
            let file = File::create(path)?;
            let mut writer = BufWriter::new(file);
            scanner.write_faithful_odd_struct(&mut writer)?;
        }
    }

    // Statistics output (-m). The banner / first-pass / second-pass blocks carry
    // timestamps + CPU times and cannot byte-match C; only the summary tail is
    // deterministic (faithful path uses the trna_results-driven writer).
    if let Some(ref path) = args.stats_file {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        // Preamble: banner + Started + run-parameter block (C tRNAscan-SE.src:178 +
        // Options.pm display_run_options). Host/date + absolute paths are runtime.
        write_stats_preamble(&mut writer, scanner, args)?;
        stats.save_firstpass_stats(&mut writer)?;
        if scanner.uses_faithful() {
            // Second-pass "Infernal Stats" block + "Overall scan speed" (C
            // Stats.pm save_final_stats, minus the isotype-count summary), then the
            // isotype/anticodon summary from the faithful writer.
            stats.save_secondpass_stats(&mut writer)?;
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
    if let Some(ref path) = args.fasta {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        if scanner.uses_faithful() {
            scanner.write_faithful_fasta(&mut writer)?;
        }
    }

    // Isotype-specific model output (-s / --isospecific)
    if let Some(ref path) = args.isospecific_file {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        if scanner.uses_faithful() {
            scanner.write_faithful_iso(&mut writer)?;
        } else {
            scanner.write_isotype_models(&mut writer)?;
        }
    }

    Ok(())
}

/// Write the `.stats` preamble: the run banner (`tRNAscan-SE v.<ver> (<rel>) scan
/// results (on host <H>)` / `Started: <date>`) plus the run-parameter block. Faithful
/// port of C tRNAscan-SE.src:178 (`$stats->write_line(...)`) + Options.pm
/// `display_run_options`. Host/date + absolute file paths are runtime values
/// (normalized in parity checks); the labels, spacing, and block structure match C.
fn write_stats_preamble<W: Write>(
    writer: &mut W,
    scanner: &TrnaScanner,
    args: &Args,
) -> std::io::Result<()> {
    // Self-identification banner: trnascan-rs reports its own name and crate
    // version here (it no longer emits C tRNAscan-SE 2.0.12's banner string).
    let host = std::process::Command::new("hostname")
        .output()
        .ok()
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .unwrap_or_else(|| "localhost".to_string());
    let now = chrono::Local::now().format("%a %b %d %H:%M:%S %Z %Y");

    // Banner: "\ntrnascan-rs v.<version> scan results (on host <host>)\nStarted: <date>".
    // (Mirrors C tRNAscan-SE's banner layout, but self-identifies as trnascan-rs.)
    writeln!(
        writer,
        "\ntrnascan-rs v.{} scan results (on host {})",
        env!("CARGO_PKG_VERSION"), host
    )?;
    writeln!(writer, "Started: {}\n", now)?;

    // ---- Run-parameter block (Options.pm display_run_options) ----
    let sep = "-".repeat(60);
    writeln!(writer, "{}", sep)?;

    let inputs: Vec<String> = args
        .input
        .iter()
        .map(|p| p.display().to_string())
        .collect();
    writeln!(writer, "Sequence file(s) to search:        {}", inputs.join(", "))?;

    let mode_name = if args.bacterial {
        "Bacterial"
    } else if args.archaeal {
        "Archaeal"
    } else if args.organellar {
        "Organellar"
    } else if args.general {
        "General"
    } else if args.mitochondrial.is_some() {
        "Mitochondrial"
    } else {
        "Eukaryotic"
    };
    writeln!(writer, "Search Mode:                       {}", mode_name)?;

    match &args.output {
        Some(p) => writeln!(writer, "Results written to:                {}", p.display())?,
        None => writeln!(writer, "Results written to:                Standard output")?,
    }
    writeln!(writer, "Output format:                     Tabular")?;
    writeln!(writer, "Searching with:                    Infernal First Pass->Infernal")?;

    // Isotype-specific model scan (Options.pm:786-793): only euk/bact/arch.
    if args.bacterial || args.archaeal {
        if args.no_isotype {
            writeln!(writer, "Isotype-specific model scan:       No")?;
        } else {
            writeln!(writer, "Isotype-specific model scan:       Yes")?;
        }
    }

    // Covariance model(s): first line labeled, the rest 35-space-indented (Options.pm:806-818).
    for (i, cm) in scanner.covariance_model_paths().iter().enumerate() {
        if i == 0 {
            writeln!(writer, "Covariance model:                  {}", cm.display())?;
        } else {
            writeln!(writer, "                                   {}", cm.display())?;
        }
    }

    // Infernal first-pass cutoff (Options.pm:882; infernal_fp_cutoff = 10).
    writeln!(writer, "Infernal first pass cutoff score:  10")?;

    // Temporary directory + output-file listing (Options.pm:894-981).
    writeln!(writer, "\nTemporary directory:               {}", std::env::temp_dir().display())?;
    if let Some(p) = &args.ss_file {
        writeln!(writer, "tRNA secondary structure")?;
        writeln!(writer, "    predictions saved to:          {}", p.display())?;
    }
    if let Some(p) = &args.bed_file {
        writeln!(writer, "tRNA predictions saved to:         {}", p.display())?;
    }
    if let Some(p) = &args.gff_file {
        writeln!(writer, "tRNA predictions saved to:         {}", p.display())?;
    }
    if let Some(p) = &args.fasta {
        writeln!(writer, "Predicted tRNA sequences")?;
        writeln!(writer, "    saved to:                      {}", p.display())?;
    }
    if let Some(p) = &args.isospecific_file {
        writeln!(writer, "Isotype specific")?;
        writeln!(writer, "    predictions saved to:          {}", p.display())?;
    }
    if let Some(p) = &args.stats_file {
        writeln!(writer, "Search statistics saved in:        {}", p.display())?;
    }
    // -H breakdown reporting (Options.pm:1003-1006, preceded by a blank line).
    if args.breakdown {
        writeln!(writer)?;
        writeln!(writer, "Reporting HMM/2' structure score breakdown")?;
    }
    writeln!(writer, "{}\n", sep)?;

    Ok(())
}

/// Write output header
fn write_header<W: Write>(writer: &mut W) -> std::io::Result<()> {
    writeln!(writer, "Sequence\t\ttRNA\tBounds\ttRNA\tAnti\tIntron Bounds\tInf\tHMM\t2'Str\tHit\tIsotype\tIsotype\t      ")?;
    writeln!(writer, "Name    \ttRNA #\tBegin\tEnd\tType\tCodon\tBegin\tEnd\tScore\tScore\tScore\tOrigin\tCM\tScore\tNote")?;
    writeln!(writer, "--------\t------\t-----\t------\t----\t-----\t-----\t----\t------\t-----\t-----\t------\t-------\t-------\t------")?;
    Ok(())
}
