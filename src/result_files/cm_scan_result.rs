//! CM scan result file module
//!
//! Corresponds to tRNAscanSE::CMscanResultFile.pm
//!
//! Parses output from cmsearch/cmscan (Infernal covariance model search).
//! Extracts hit information including:
//! - Target sequence name and bounds
//! - Score and E-value
//! - Strand orientation
//! - Secondary structure alignment
//! - Model and target sequence alignment

use std::fs::File;
use std::io::{self, BufRead, BufReader, Seek, SeekFrom};
use std::path::{Path, PathBuf};

use super::{ResultError, ResultFile, IndexRecord};

/// CM scan result file reader
pub struct CMscanResultFile {
    path: PathBuf,
    reader: Option<BufReader<File>>,
    indexes: Vec<IndexRecord>,
}

impl CMscanResultFile {
    pub fn new(path: &Path, _file_type: &str) -> Self {
        Self {
            path: path.to_path_buf(),
            reader: None,
            indexes: Vec::new(),
        }
    }

    /// Open for reading
    pub fn open_read(&mut self) -> io::Result<()> {
        let file = File::open(&self.path)?;
        self.reader = Some(BufReader::new(file));
        Ok(())
    }

    /// Sort and index cmsearch records
    pub fn sort_cmsearch_records(&mut self) -> io::Result<()> {
        if self.reader.is_none() {
            self.open_read()?;
        }

        self.indexes.clear();

        if let Some(ref mut r) = self.reader {
            let mut file_pos = 0u64;
            let mut in_hits = false;
            let mut current_seq = String::new();
            let mut current_pos = 0u64;

            for line_result in r.lines() {
                let line = line_result?;

                if line.contains("Hit alignments:") {
                    in_hits = true;
                    file_pos += line.len() as u64 + 1;
                    continue;
                }

                if in_hits {
                    // Match sequence header: >> seqname
                    if line.starts_with(">>") {
                        if let Some(seq) = line.strip_prefix(">>").map(|s| s.trim()) {
                            current_seq = seq.to_string();
                            current_pos = file_pos;
                        }
                    }
                    // Match hit line with score/bounds
                    else if let Some(record) = Self::parse_hit_line(&line, &current_seq, current_pos) {
                        self.indexes.push(record);
                    }
                }

                file_pos += line.len() as u64 + 1;
            }

            // Sort by sequence name, then position
            self.sort_indexes();
        }

        Ok(())
    }

    /// Parse a cmsearch hit line
    fn parse_hit_line(line: &str, seq_name: &str, file_pos: u64) -> Option<IndexRecord> {
        // Format: (rank) model e-val score bias mdl_from mdl_to seq_from seq_to strand ...
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 13 {
            return None;
        }

        // Check if starts with (number)
        if !parts[0].starts_with('(') {
            return None;
        }

        let score: f64 = parts[2].parse().ok()?;
        let mut start: i64 = parts[7].parse().ok()?;
        let mut end: i64 = parts[8].parse().ok()?;
        let strand = parts[9];
        let trunc = parts[12];

        // Swap start/end for reverse strand
        if strand == "-" {
            std::mem::swap(&mut start, &mut end);
        }

        let trunc_str = if trunc == "no" {
            String::new()
        } else {
            trunc.to_string()
        };

        Some(IndexRecord::with_data(
            file_pos,
            vec![
                seq_name.to_string(),
                start.to_string(),
                end.to_string(),
                strand.to_string(),
                score.to_string(),
                trunc_str,
            ],
        ))
    }

    /// Sort indexes by tRNAscanSE output order
    fn sort_indexes(&mut self) {
        self.indexes.sort_by(|a, b| {
            let a_strand = a.data.get(3).map(|s| s.as_str()).unwrap_or("+");
            let b_strand = b.data.get(3).map(|s| s.as_str()).unwrap_or("+");
            let a_start: i64 = a.data.get(1).and_then(|s| s.parse().ok()).unwrap_or(0);
            let b_start: i64 = b.data.get(1).and_then(|s| s.parse().ok()).unwrap_or(0);
            let a_end: i64 = a.data.get(2).and_then(|s| s.parse().ok()).unwrap_or(0);
            let b_end: i64 = b.data.get(2).and_then(|s| s.parse().ok()).unwrap_or(0);
            let a_score: f64 = a.data.get(4).and_then(|s| s.parse().ok()).unwrap_or(0.0);
            let b_score: f64 = b.data.get(4).and_then(|s| s.parse().ok()).unwrap_or(0.0);

            // Sort by seqname, then strand, then position, then score
            a.data[0].cmp(&b.data[0])
                .then_with(|| {
                    if a_strand == b_strand && a_strand == "+" {
                        a_strand.cmp(b_strand)
                            .then_with(|| a_start.cmp(&b_start))
                            .then_with(|| b_score.partial_cmp(&a_score).unwrap_or(std::cmp::Ordering::Equal))
                    } else if a_strand == b_strand && a_strand == "-" {
                        b_strand.cmp(a_strand)
                            .then_with(|| b_end.cmp(&a_end))
                            .then_with(|| b_score.partial_cmp(&a_score).unwrap_or(std::cmp::Ordering::Equal))
                    } else {
                        a_strand.cmp(b_strand)
                    }
                })
        });
    }

    /// Get cmsearch record at file position
    pub fn get_record(&mut self, file_pos: u64, format: bool) -> Result<CMsearchRecord, ResultError> {
        if self.reader.is_none() {
            self.open_read()?;
        }

        if let Some(ref mut r) = self.reader {
            r.seek(SeekFrom::Start(file_pos))?;

            let mut seq_name = String::new();
            let mut ss = String::new();
            let mut seq = String::new();
            let mut model = String::new();
            let mut nc = String::new();

            // Collect lines into a vector to avoid borrowing conflicts
            let mut lines = Vec::new();
            for line_result in r.lines() {
                let line = line_result?;
                if line.contains("Internal CM pipeline") {
                    break;
                }
                if !seq_name.is_empty() && line.starts_with(">>") {
                    break;
                }
                lines.push(line);
            }

            // Process collected lines
            let mut i = 0;
            while i < lines.len() {
                let line = &lines[i];

                if line.starts_with(">>") {
                    if seq_name.is_empty() {
                        seq_name = line.strip_prefix(">>")
                            .map(|s| s.trim().to_string())
                            .unwrap_or_default();
                    }
                }
                // Parse NC line
                else if line.ends_with(" NC") {
                    let nc_line = line.strip_suffix(" NC").unwrap_or("");
                    nc = nc_line.to_string();
                }
                // Parse CS (consensus structure) line
                else if let Some(cs_match) = line.strip_suffix(" CS") {
                    if let Some(structure) = Self::extract_structure(cs_match) {
                        ss.push_str(&structure);

                        // Read next 4 lines: model seq, spacer, target seq, PP
                        if i + 1 < lines.len() {
                            if let Some(m) = Self::extract_sequence(&lines[i + 1]) {
                                model.push_str(&m);
                            }
                        }
                        // Skip spacer (i+2)
                        if i + 3 < lines.len() {
                            let seq_line_clean = lines[i + 3].replace("*[0]*", "-----");
                            if let Some(s) = Self::extract_sequence(&seq_line_clean) {
                                seq.push_str(&s);
                                // Extract NC segment
                                let seq_len = s.len();
                                if nc.len() >= seq_len {
                                    let nc_segment = &nc[nc.len() - seq_len..];
                                    nc = nc_segment.to_string();
                                }
                            }
                        }
                        // Skip PP line (i+4)
                        i += 4;
                    }
                }
                i += 1;
            }

            if format {
                let (formatted_ss, formatted_seq) = Self::format_output(&ss, &seq, &nc);
                return Ok(CMsearchRecord {
                    seq_name,
                    ss: formatted_ss,
                    seq: formatted_seq,
                    model,
                    nc,
                });
            }

            Ok(CMsearchRecord {
                seq_name,
                ss,
                seq,
                model,
                nc,
            })
        } else {
            Err(ResultError::Io(io::Error::new(io::ErrorKind::Other, "Reader not open")))
        }
    }

    fn extract_structure(line: &str) -> Option<String> {
        let parts: Vec<&str> = line.split_whitespace().collect();
        // Structure is the last field before CS marker
        if parts.len() >= 2 {
            Some(parts[parts.len() - 1].to_string())
        } else {
            None
        }
    }

    fn extract_sequence(line: &str) -> Option<String> {
        let parts: Vec<&str> = line.split_whitespace().collect();
        // Sequence is between position numbers: seqname start SEQ end
        if parts.len() >= 4 {
            Some(parts[parts.len() - 2].to_string())
        } else {
            None
        }
    }

    /// Format cmsearch output (convert structure notation, remove gaps)
    fn format_output(ss: &str, seq: &str, nc: &str) -> (String, String) {
        let mut formatted_seq = seq.replace("U", "T").replace("u", "t");
        let mut formatted_ss = ss.to_string();

        // Fix mismatches based on NC line
        formatted_ss = Self::fix_mismatch_ss(&formatted_ss, &formatted_seq, nc);

        // Mark gaps with * and then remove them
        let mut seq_chars: Vec<char> = formatted_seq.chars().collect();
        let mut ss_chars: Vec<char> = formatted_ss.chars().collect();

        for i in 0..seq_chars.len() {
            if seq_chars[i] == '-' {
                seq_chars[i] = '*';
                if i < ss_chars.len() {
                    ss_chars[i] = '*';
                }
            }
        }

        formatted_seq = seq_chars.iter().filter(|&&c| c != '*').collect();
        formatted_ss = ss_chars.iter().filter(|&&c| c != '*').collect();

        // Convert structure notation
        formatted_ss = formatted_ss
            .replace(&[',', '_', '-', ':'][..], ".")
            .replace('>', "@")
            .replace(')', "@")
            .replace('<', ">")
            .replace('(', ">")
            .replace('@', "<");

        // Pad structure to match sequence length
        while formatted_ss.len() < formatted_seq.len() {
            formatted_ss.push('.');
        }

        (formatted_ss, formatted_seq)
    }

    /// Fix mismatches in secondary structure based on NC (non-canonical) markers
    fn fix_mismatch_ss(ss: &str, seq: &str, nc: &str) -> String {
        let mut result = ss.to_string();
        let seq_bytes = seq.as_bytes();
        let ss_bytes = unsafe { result.as_bytes_mut() };
        let nc_bytes = nc.as_bytes();

        // Mark non-canonical positions as unpaired
        for i in 0..nc_bytes.len().min(ss_bytes.len()) {
            if nc_bytes[i] == b'v' {
                ss_bytes[i] = b'.';
            }
        }

        // Check base pairing validity
        let mut left_stack: Vec<usize> = Vec::new();
        let mut pairs: Vec<(usize, usize)> = Vec::new();

        for (i, &c) in ss_bytes.iter().enumerate() {
            if c == b'<' || c == b'(' {
                left_stack.push(i);
            } else if c == b'>' || c == b')' {
                if let Some(left_pos) = left_stack.pop() {
                    pairs.push((left_pos, i));
                }
            }
        }

        // Validate pairs
        for (left, right) in pairs {
            if left < seq_bytes.len() && right < seq_bytes.len() {
                let left_base = seq_bytes[left].to_ascii_uppercase();
                let right_base = seq_bytes[right].to_ascii_uppercase();

                let valid = match (left_base, right_base) {
                    (b'A', b'U') | (b'A', b'T') |
                    (b'U', b'A') | (b'U', b'G') |
                    (b'T', b'A') | (b'T', b'G') |
                    (b'G', b'C') | (b'G', b'U') | (b'G', b'T') |
                    (b'C', b'G') => true,
                    _ => false,
                };

                if !valid || left_base == b'-' || right_base == b'-' {
                    ss_bytes[left] = b'.';
                    ss_bytes[right] = b'.';
                }
            }
        }

        String::from_utf8(ss_bytes.to_vec()).unwrap_or_else(|_| ss.to_string())
    }

    /// Read tabular format output (--tblout)
    pub fn get_next_tab_hits(&mut self) -> Result<Vec<CMsearchTabRecord>, ResultError> {
        if self.reader.is_none() {
            self.open_read()?;
        }

        let mut hits = Vec::new();
        let mut last_seq = String::new();

        if let Some(ref mut r) = self.reader {
            let file_pos = r.stream_position()?;

            for line_result in r.lines() {
                let line = line_result?;

                // Skip comment lines
                if line.starts_with('#') {
                    continue;
                }

                let columns: Vec<&str> = line.split_whitespace().collect();
                if columns.len() < 17 {
                    continue;
                }

                let seq_name = columns[3].to_string();

                // If we've moved to a new sequence and have hits, stop
                if seq_name != last_seq {
                    if !hits.is_empty() {
                        r.seek(SeekFrom::Start(file_pos))?;
                        break;
                    }
                    last_seq = seq_name.clone();
                }

                // Parse hit
                let mut start: i64 = columns[9].parse().unwrap_or(0);
                let mut end: i64 = columns[10].parse().unwrap_or(0);
                let strand = columns[11];

                if strand == "-" {
                    std::mem::swap(&mut start, &mut end);
                }

                hits.push(CMsearchTabRecord {
                    model: columns[1].to_string(),
                    seq_name,
                    start,
                    end,
                    strand: strand.to_string(),
                    score: columns[16].parse().unwrap_or(0.0),
                });
            }
        }

        Ok(hits)
    }

    /// Get indexes
    pub fn get_indexes(&self) -> &[IndexRecord] {
        &self.indexes
    }

    /// Clear indexes
    pub fn clear_indexes(&mut self) {
        self.indexes.clear();
    }
}

impl ResultFile for CMscanResultFile {
    fn open(path: &Path) -> Result<Self, ResultError> {
        let mut file = Self::new(path, "cmsearch");
        file.open_read()?;
        Ok(file)
    }

    fn write_header(&mut self) -> io::Result<()> {
        // CM scan files are read-only (generated by Infernal)
        Ok(())
    }

    fn close(&mut self) -> io::Result<()> {
        self.reader = None;
        Ok(())
    }
}

// ============================================================================
// CM Search Records
// ============================================================================

/// CM search alignment record
#[derive(Debug, Clone, Default)]
pub struct CMsearchRecord {
    pub seq_name: String,
    pub ss: String,
    pub seq: String,
    pub model: String,
    pub nc: String,
}

/// CM search tabular record
#[derive(Debug, Clone, Default)]
pub struct CMsearchTabRecord {
    pub model: String,
    pub seq_name: String,
    pub start: i64,
    pub end: i64,
    pub strand: String,
    pub score: f64,
}
