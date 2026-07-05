//! Generic result file reader module
//!
//! Corresponds to tRNAscanSE::ResultFileReader.pm
//!
//! Provides functions to read various tRNAscan-SE output formats:
//! - Tabular output (.out)
//! - Secondary structure output (.ss)
//! - Sprinzl position maps
//! - Non-canonical feature files

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

use super::{ResultError, ResultFormat};

/// Generic result file reader
pub struct ResultFileReader {
    #[allow(dead_code)]
    path: std::path::PathBuf,
    reader: BufReader<File>,
    format: ResultFormat,
    line_num: usize,
}

impl ResultFileReader {
    /// Open a result file for reading
    pub fn open(path: &Path) -> Result<Self, ResultError> {
        let file = File::open(path)?;
        let mut reader = BufReader::new(file);

        // Detect format from contents
        let format = ResultFormat::from_contents(&mut reader)?;

        Ok(Self {
            path: path.to_path_buf(),
            reader,
            format,
            line_num: 0,
        })
    }

    /// Get the detected format
    pub fn format(&self) -> ResultFormat {
        self.format
    }

    /// Read next line
    fn next_line(&mut self) -> io::Result<Option<String>> {
        let mut line = String::new();
        match self.reader.read_line(&mut line) {
            Ok(0) => Ok(None),
            Ok(_) => {
                self.line_num += 1;
                Ok(Some(line.trim_end().to_string()))
            }
            Err(e) => Err(e),
        }
    }

    /// Read all tabular output records
    pub fn read_tabular_output(&mut self) -> Result<Vec<TabularRecord>, ResultError> {
        let mut records = Vec::new();
        let mut header = TabularHeader::default();

        while let Some(line) = self.next_line()? {
            if line.is_empty() {
                continue;
            }

            let columns: Vec<String> = line.split('\t')
                .map(|s| s.trim().to_string())
                .collect();

            // Parse headers
            if line.starts_with("Sequence") {
                header.parse_line1(&columns);
            } else if line.starts_with("Name") {
                header.parse_line2(&columns);
            } else if line.starts_with("---") {
                continue;
            } else {
                // Parse data line
                if let Ok(record) = TabularRecord::from_columns(&columns, &header) {
                    records.push(record);
                }
            }
        }

        Ok(records)
    }

    /// Read all secondary structure records
    pub fn read_ss_output(&mut self) -> Result<Vec<SSRecord>, ResultError> {
        let mut records = Vec::new();
        let mut current: Option<SSRecord> = None;

        while let Some(line) = self.next_line()? {
            if line.is_empty() {
                // End of record
                if let Some(record) = current.take() {
                    records.push(record);
                }
                continue;
            }

            // Header line: seqname.trna# (coords) Length: N bp
            if let Some((id_part, rest)) = line.split_once(' ') {
                if rest.contains("Length:") {
                    if let Some(record) = current.take() {
                        records.push(record);
                    }
                    current = Some(SSRecord::new(id_part.to_string()));
                    continue;
                }
            }

            if let Some(ref mut record) = current {
                // Type line
                if line.starts_with("Type:") {
                    record.parse_type_line(&line);
                }
                // Score line
                else if line.starts_with("HMM Sc=") || line.starts_with("pseudogene:") {
                    record.parse_score_line(&line);
                }
                // Intron line
                else if line.starts_with("Possible intron:") {
                    record.parse_intron_line(&line);
                }
                // Sequence line
                else if line.starts_with("Seq:") {
                    record.seq = line.strip_prefix("Seq:").unwrap_or("").trim().to_string();
                }
                // Structure line
                else if line.starts_with("Str:") {
                    record.ss = line.strip_prefix("Str:").unwrap_or("").trim().to_string();
                }
                // Pre-tRNA line
                else if line.starts_with("Pre:") || line.starts_with("PRE:") || line.starts_with("BHB:") {
                    let pre_seq = line.split(':').nth(1).unwrap_or("").trim();
                    let pre_seq_clean = pre_seq.replace(&['[', ']'][..], "");
                    record.mat_seq = record.seq.clone();
                    record.seq = pre_seq_clean;
                }
            }
        }

        // Add final record
        if let Some(record) = current {
            records.push(record);
        }

        Ok(records)
    }

    /// Read GFF3 output
    pub fn read_gff_output(&mut self) -> Result<Vec<GFFRecord>, ResultError> {
        let mut records = Vec::new();

        while let Some(line) = self.next_line()? {
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            if let Ok(record) = GFFRecord::from_line(&line) {
                records.push(record);
            }
        }

        Ok(records)
    }

    /// Read BED output
    pub fn read_bed_output(&mut self) -> Result<Vec<BEDRecord>, ResultError> {
        let mut records = Vec::new();

        while let Some(line) = self.next_line()? {
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            if let Ok(record) = BEDRecord::from_line(&line) {
                records.push(record);
            }
        }

        Ok(records)
    }
}

// ============================================================================
// Tabular Output Records
// ============================================================================

#[derive(Debug, Default)]
struct TabularHeader {
    seqname: usize,
    trna_id: usize,
    start: usize,
    end: usize,
    isotype: usize,
    anticodon: usize,
    intron_start: usize,
    intron_end: usize,
    score: usize,
    isotype_cm: Option<usize>,
    isotype_score: Option<usize>,
    note: Option<usize>,
}

impl TabularHeader {
    fn parse_line1(&mut self, columns: &[String]) {
        for (i, col) in columns.iter().enumerate() {
            match col.as_str() {
                "Sequence" => self.seqname = i,
                "Anti" => self.anticodon = i,
                "Intron" => {
                    self.intron_start = i;
                    self.intron_end = i + 1;
                }
                "Inf" | "Cove" => self.score = i,
                _ => {}
            }
        }
    }

    fn parse_line2(&mut self, columns: &[String]) {
        for (i, col) in columns.iter().enumerate() {
            match col.as_str() {
                "tRNA#" | "tRNA #" => self.trna_id = i,
                "Begin" => {
                    self.start = i;
                    self.end = i + 1;
                }
                "Type" => self.isotype = i,
                "CM" => {
                    self.isotype_cm = Some(i);
                    self.isotype_score = Some(i + 1);
                }
                "Note" => self.note = Some(i),
                _ => {}
            }
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct TabularRecord {
    pub seqname: String,
    pub trna_num: i32,
    pub start: i64,
    pub end: i64,
    pub isotype: String,
    pub anticodon: String,
    pub intron_start: i64,
    pub intron_end: i64,
    pub score: f64,
    pub isotype_cm: String,
    pub isotype_score: f64,
    pub note: String,
    pub strand: String,
}

impl TabularRecord {
    fn from_columns(columns: &[String], header: &TabularHeader) -> Result<Self, ResultError> {
        let start: i64 = columns.get(header.start)
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);
        let end: i64 = columns.get(header.end)
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);

        let strand = if start < end {
            "+".to_string()
        } else {
            "-".to_string()
        };

        Ok(Self {
            seqname: columns.get(header.seqname).cloned().unwrap_or_default(),
            trna_num: columns.get(header.trna_id)
                .and_then(|s| s.parse().ok())
                .unwrap_or(0),
            start,
            end,
            isotype: columns.get(header.isotype).cloned().unwrap_or_default(),
            anticodon: columns.get(header.anticodon).cloned().unwrap_or_default(),
            intron_start: columns.get(header.intron_start)
                .and_then(|s| s.parse().ok())
                .unwrap_or(0),
            intron_end: columns.get(header.intron_end)
                .and_then(|s| s.parse().ok())
                .unwrap_or(0),
            score: columns.get(header.score)
                .and_then(|s| s.parse().ok())
                .unwrap_or(0.0),
            isotype_cm: header.isotype_cm
                .and_then(|idx| columns.get(idx))
                .cloned()
                .unwrap_or_default(),
            isotype_score: header.isotype_score
                .and_then(|idx| columns.get(idx))
                .and_then(|s| s.parse().ok())
                .unwrap_or(0.0),
            note: header.note
                .and_then(|idx| columns.get(idx))
                .cloned()
                .unwrap_or_default(),
            strand,
        })
    }
}

// ============================================================================
// Secondary Structure Records
// ============================================================================

#[derive(Debug, Clone, Default)]
pub struct SSRecord {
    pub trnascan_id: String,
    pub seqname: String,
    pub start: i64,
    pub end: i64,
    pub strand: String,
    pub isotype: String,
    pub anticodon: String,
    pub anticodon_pos: Vec<(i32, i32)>,
    pub score: f64,
    pub hmm_score: f64,
    pub ss_score: f64,
    pub is_pseudo: bool,
    pub introns: Vec<(i32, i32, i64, i64)>,
    pub seq: String,
    pub mat_seq: String,
    pub ss: String,
}

impl SSRecord {
    fn new(id: String) -> Self {
        let seqname = if let Some(pos) = id.rfind('.') {
            id[..pos].to_string()
        } else {
            id.clone()
        };

        Self {
            trnascan_id: id,
            seqname,
            ..Default::default()
        }
    }

    fn parse_type_line(&mut self, line: &str) {
        // Type: Ala    Anticodon: TGC at 34-36    Score: 77.4
        let parts: Vec<&str> = line.split_whitespace().collect();
        for (i, part) in parts.iter().enumerate() {
            match *part {
                "Type:" => {
                    if let Some(isotype) = parts.get(i + 1) {
                        self.isotype = isotype.to_string();
                    }
                }
                "Anticodon:" => {
                    if let Some(ac) = parts.get(i + 1) {
                        self.anticodon = ac.to_string();
                    }
                    if let Some(pos) = parts.get(i + 3) {
                        if let Some((s, e)) = pos.split_once('-') {
                            if let (Ok(start), Ok(end)) = (s.parse(), e.parse()) {
                                self.anticodon_pos.push((start, end));
                            }
                        }
                    }
                }
                "Score:" => {
                    if let Some(score) = parts.get(i + 1) {
                        self.score = score.parse().unwrap_or(0.0);
                    }
                }
                _ => {}
            }
        }
    }

    fn parse_score_line(&mut self, line: &str) {
        if line.contains("pseudogene") {
            self.is_pseudo = true;
        }

        // HMM Sc=75.2  Sec struct Sc=2.2
        if let Some(hmm_part) = line.split("HMM Sc=").nth(1) {
            if let Some(hmm_str) = hmm_part.split_whitespace().next() {
                self.hmm_score = hmm_str.parse().unwrap_or(0.0);
            }
        }
        if let Some(ss_part) = line.split("Sec struct Sc=").nth(1) {
            if let Some(ss_str) = ss_part.split_whitespace().next() {
                self.ss_score = ss_str.parse().unwrap_or(0.0);
            }
        }
    }

    fn parse_intron_line(&mut self, line: &str) {
        // Possible intron: 37-38 (1234-1235)
        if let Some(rest) = line.strip_prefix("Possible intron:") {
            let parts: Vec<&str> = rest.trim().split_whitespace().collect();
            if parts.len() >= 2 {
                // Parse relative positions
                if let Some((rs, re)) = parts[0].split_once('-') {
                    if let (Ok(rel_start), Ok(rel_end)) = (rs.parse(), re.parse()) {
                        // Parse absolute positions
                        let abs_part = parts[1].trim_matches(|c| c == '(' || c == ')');
                        if let Some((s, e)) = abs_part.split_once('-') {
                            if let (Ok(start), Ok(end)) = (s.parse(), e.parse()) {
                                self.introns.push((rel_start, rel_end, start, end));
                            }
                        }
                    }
                }
            }
        }
    }
}

// ============================================================================
// GFF Records
// ============================================================================

#[derive(Debug, Clone, Default)]
pub struct GFFRecord {
    pub seqname: String,
    pub source: String,
    pub feature_type: String,
    pub start: i64,
    pub end: i64,
    pub score: f64,
    pub strand: String,
    pub frame: String,
    pub attributes: HashMap<String, String>,
}

impl GFFRecord {
    fn from_line(line: &str) -> Result<Self, ResultError> {
        let columns: Vec<&str> = line.split('\t').collect();
        if columns.len() < 9 {
            return Err(ResultError::Parse(format!("Invalid GFF line: {}", line)));
        }

        let mut attributes = HashMap::new();
        for pair in columns[8].split(';') {
            if let Some((key, val)) = pair.split_once('=') {
                attributes.insert(key.to_string(), val.to_string());
            }
        }

        Ok(Self {
            seqname: columns[0].to_string(),
            source: columns[1].to_string(),
            feature_type: columns[2].to_string(),
            start: columns[3].parse().unwrap_or(0),
            end: columns[4].parse().unwrap_or(0),
            score: columns[5].parse().unwrap_or(0.0),
            strand: columns[6].to_string(),
            frame: columns[7].to_string(),
            attributes,
        })
    }
}

// ============================================================================
// BED Records
// ============================================================================

#[derive(Debug, Clone, Default)]
pub struct BEDRecord {
    pub chrom: String,
    pub chrom_start: i64,
    pub chrom_end: i64,
    pub name: String,
    pub score: i32,
    pub strand: String,
    pub thick_start: i64,
    pub thick_end: i64,
    pub item_rgb: String,
    pub block_count: i32,
    pub block_sizes: Vec<i32>,
    pub block_starts: Vec<i32>,
}

impl BEDRecord {
    fn from_line(line: &str) -> Result<Self, ResultError> {
        let columns: Vec<&str> = line.split('\t').collect();
        if columns.len() < 12 {
            return Err(ResultError::Parse(format!("Invalid BED line: {}", line)));
        }

        let block_sizes: Vec<i32> = columns[10]
            .trim_end_matches(',')
            .split(',')
            .filter_map(|s| s.parse().ok())
            .collect();

        let block_starts: Vec<i32> = columns[11]
            .trim_end_matches(',')
            .split(',')
            .filter_map(|s| s.parse().ok())
            .collect();

        Ok(Self {
            chrom: columns[0].to_string(),
            chrom_start: columns[1].parse().unwrap_or(0),
            chrom_end: columns[2].parse().unwrap_or(0),
            name: columns[3].to_string(),
            score: columns[4].parse().unwrap_or(0),
            strand: columns[5].to_string(),
            thick_start: columns[6].parse().unwrap_or(0),
            thick_end: columns[7].parse().unwrap_or(0),
            item_rgb: columns[8].to_string(),
            block_count: columns[9].parse().unwrap_or(0),
            block_sizes,
            block_starts,
        })
    }
}
