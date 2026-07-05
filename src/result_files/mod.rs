//! Result file modules for tRNAscan-SE
//!
//! This module provides file I/O for various result file formats used
//! in the tRNAscan-SE pipeline:
//!
//! - **IntResultFile**: Intermediate results with full tRNA details
//! - **FpScanResultFile**: First-pass scan results (EufindtRNA/tRNAscan)
//! - **CMscanResultFile**: Covariance model scan results (cmscan/cmsearch)
//! - **MultiResultFile**: Multi-model scoring results
//! - **ResultFileReader**: Generic reader for tRNAscan-SE output files
//!
//! Ported from tRNAscanSE::*ResultFile.pm modules

use std::fs::File;
use std::io::{self, BufRead, BufReader, Seek, SeekFrom};
use std::path::Path;
use thiserror::Error;

// ============================================================================
// Error Types
// ============================================================================

#[derive(Debug, Error)]
pub enum ResultError {
    #[error("I/O error: {0}")]
    Io(#[from] io::Error),

    #[error("Parse error: {0}")]
    Parse(String),

    #[error("Invalid format: {0}")]
    InvalidFormat(String),

    #[error("File not found: {0}")]
    FileNotFound(String),

    #[error("Invalid record at line {line}: {msg}")]
    InvalidRecord { line: usize, msg: String },
}

// ============================================================================
// Common Traits
// ============================================================================

/// Trait for result file writers/readers
pub trait ResultFile {
    /// Open a result file at the given path
    fn open(path: &Path) -> Result<Self, ResultError> where Self: Sized;

    /// Write file header if needed
    fn write_header(&mut self) -> io::Result<()>;

    /// Close the file (flushes buffers)
    fn close(&mut self) -> io::Result<()>;
}

/// Trait for result records that can be serialized/deserialized
pub trait ResultRecord {
    /// Convert this record to a tab-delimited line
    fn to_line(&self) -> String;

    /// Parse a tab-delimited line into this record
    fn from_line(line: &str) -> Result<Self, ResultError> where Self: Sized;
}

// ============================================================================
// Result Format Detection
// ============================================================================

/// Result file formats supported by tRNAscan-SE
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ResultFormat {
    /// Tabular output (.out)
    Tabular,
    /// Secondary structure (.ss)
    SecondaryStructure,
    /// BED format (.bed)
    Bed,
    /// GFF3 format (.gff)
    Gff,
    /// Intermediate format (internal)
    Intermediate,
    /// First-pass format (internal)
    FirstPass,
    /// CM scan format (internal)
    CMScan,
    /// Unknown format
    Unknown,
}

impl ResultFormat {
    /// Detect format from file path extension
    pub fn from_path(path: &Path) -> Self {
        match path.extension().and_then(|s| s.to_str()) {
            Some("out") => ResultFormat::Tabular,
            Some("ss") => ResultFormat::SecondaryStructure,
            Some("bed") => ResultFormat::Bed,
            Some("gff") | Some("gff3") => ResultFormat::Gff,
            Some("int") => ResultFormat::Intermediate,
            Some("fp") => ResultFormat::FirstPass,
            Some("cm") => ResultFormat::CMScan,
            _ => ResultFormat::Unknown,
        }
    }

    /// Detect format from file contents
    pub fn from_contents(reader: &mut BufReader<File>) -> io::Result<Self> {
        let mut first_line = String::new();
        reader.read_line(&mut first_line)?;
        reader.seek(SeekFrom::Start(0))?;

        if first_line.starts_with("##gff-version") {
            Ok(ResultFormat::Gff)
        } else if first_line.starts_with("Sequence\t\ttRNA") {
            Ok(ResultFormat::Tabular)
        } else if first_line.contains(".trna") && first_line.contains("Length:") {
            Ok(ResultFormat::SecondaryStructure)
        } else if first_line.starts_with("seqname\tordered_seqname\tid") {
            Ok(ResultFormat::Intermediate)
        } else {
            Ok(ResultFormat::Unknown)
        }
    }
}

// ============================================================================
// File Index Record (for seek-based access)
// ============================================================================

/// Index record for fast file seeking
#[derive(Debug, Clone)]
pub struct IndexRecord {
    /// File position (byte offset)
    pub file_pos: u64,
    /// Additional indexing data (format-specific)
    pub data: Vec<String>,
}

impl IndexRecord {
    pub fn new(file_pos: u64) -> Self {
        Self {
            file_pos,
            data: Vec::new(),
        }
    }

    pub fn with_data(file_pos: u64, data: Vec<String>) -> Self {
        Self { file_pos, data }
    }
}

// ============================================================================
// Helper Functions
// ============================================================================

/// Parse intron string format: "rel_start-rel_end,start-end,type,seq;"
pub fn parse_intron(s: &str) -> Result<Vec<(i32, i32, i64, i64, String, String)>, ResultError> {
    let mut introns = Vec::new();

    for part in s.split(';') {
        if part.is_empty() {
            continue;
        }

        let fields: Vec<&str> = part.split(',').collect();
        if fields.len() < 4 {
            return Err(ResultError::Parse(format!("Invalid intron format: {}", s)));
        }

        // Parse rel_start-rel_end
        let rel_range: Vec<&str> = fields[0].split('-').collect();
        if rel_range.len() != 2 {
            return Err(ResultError::Parse(format!("Invalid relative range: {}", fields[0])));
        }
        let rel_start: i32 = rel_range[0].parse()
            .map_err(|_| ResultError::Parse(format!("Invalid rel_start: {}", rel_range[0])))?;
        let rel_end: i32 = rel_range[1].parse()
            .map_err(|_| ResultError::Parse(format!("Invalid rel_end: {}", rel_range[1])))?;

        // Parse start-end
        let abs_range: Vec<&str> = fields[1].split('-').collect();
        if abs_range.len() != 2 {
            return Err(ResultError::Parse(format!("Invalid absolute range: {}", fields[1])));
        }
        let start: i64 = abs_range[0].parse()
            .map_err(|_| ResultError::Parse(format!("Invalid start: {}", abs_range[0])))?;
        let end: i64 = abs_range[1].parse()
            .map_err(|_| ResultError::Parse(format!("Invalid end: {}", abs_range[1])))?;

        let intron_type = fields[2].to_string();
        let seq = fields[3].to_string();

        introns.push((rel_start, rel_end, start, end, intron_type, seq));
    }

    Ok(introns)
}

/// Format introns to string
pub fn format_introns(introns: &[(i32, i32, i64, i64, String, String)]) -> String {
    introns.iter()
        .map(|(rs, re, s, e, t, seq)| format!("{}-{},{}-{},{},{}", rs, re, s, e, t, seq))
        .collect::<Vec<_>>()
        .join(";")
}

/// Parse anticodon position string format: "start-end;"
pub fn parse_ac_pos(s: &str) -> Result<Vec<(i32, i32)>, ResultError> {
    let mut positions = Vec::new();

    for part in s.split(';') {
        if part.is_empty() {
            continue;
        }

        let range: Vec<&str> = part.split('-').collect();
        if range.len() != 2 {
            return Err(ResultError::Parse(format!("Invalid anticodon position: {}", part)));
        }

        let start: i32 = range[0].parse()
            .map_err(|_| ResultError::Parse(format!("Invalid ac_pos start: {}", range[0])))?;
        let end: i32 = range[1].parse()
            .map_err(|_| ResultError::Parse(format!("Invalid ac_pos end: {}", range[1])))?;

        positions.push((start, end));
    }

    Ok(positions)
}

/// Format anticodon positions to string
pub fn format_ac_pos(positions: &[(i32, i32)]) -> String {
    positions.iter()
        .map(|(s, e)| format!("{}-{}", s, e))
        .collect::<Vec<_>>()
        .join(";")
}

/// Parse optional float field (empty string -> None)
pub fn parse_opt_f64(s: &str) -> Option<f64> {
    if s.is_empty() {
        None
    } else {
        s.parse().ok()
    }
}

/// Parse optional int field (empty string -> None)
pub fn parse_opt_i32(s: &str) -> Option<i32> {
    if s.is_empty() {
        None
    } else {
        s.parse().ok()
    }
}

// ============================================================================
// Submodules
// ============================================================================

pub mod int_result;
pub mod fp_scan_result;
pub mod cm_scan_result;
pub mod multi_result;
pub mod reader;

pub use int_result::IntResultFile;
pub use fp_scan_result::FpScanResultFile;
pub use cm_scan_result::CMscanResultFile;
pub use multi_result::MultiResultFile;
pub use reader::ResultFileReader;
