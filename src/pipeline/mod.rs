//! Main tRNA scanning pipeline for tRNAscan-SE.
//!
//! This module provides the complete tRNA detection and analysis pipeline,
//! including:
//! - Configuration and command-line option handling
//! - Main scanning pipeline coordinating first-pass and CM-based detection
//! - Output formatting in multiple formats (.out, .ss, .stats, .bed, .iso)
//!
//! # Example Usage
//!
//! ```no_run
//! use trnascan_rs::pipeline::{TrnaScanConfig, TrnaScanner, SearchMode};
//! use std::path::PathBuf;
//!
//! // Create configuration
//! let config = TrnaScanConfig::new(SearchMode::Eukaryotic)
//!     .with_output_file(Some(PathBuf::from("results.out")))
//!     .with_stats_file(Some(PathBuf::from("results.stats")));
//!
//! // Create scanner and process file
//! let scanner = TrnaScanner::new(config);
//! let results = scanner.scan_file(&PathBuf::from("input.fa")).unwrap();
//!
//! // Access results
//! for hit in &results.hits {
//!     println!("{}: {} at {}-{}", hit.seq_name, hit.isotype, hit.start, hit.end);
//! }
//! ```

pub mod config;
pub mod output;
pub mod scanner;

// Re-export commonly used items
pub use config::{
    FirstPassMethod, MitoModel, OutputConfig, SearchMode, TrnaScanConfig,
};
pub use output::{
    format_bed_line, format_iso_line, format_ss_entry, format_tabular_header,
    format_tabular_line, BedFormatter, IsoFormatter, SsFormatter, StatsFormatter,
    TabularFormatter,
};
pub use scanner::{HitOrigin, ScanResults, TrnaHit, TrnaScanner};
