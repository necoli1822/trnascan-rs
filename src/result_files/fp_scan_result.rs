//! First-pass scan result file module
//!
//! Corresponds to tRNAscanSE::FpScanResultFile.pm
//!
//! This file format stores first-pass scanning results from
//! EufindtRNA, tRNAscan-SE first pass, or Infernal first pass.
//!
//! Format: Tab-delimited with columns:
//! Seqname, tRNA#, Begin, End, Type, Codon, SeqID, SeqLen, Score, Model

use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use super::{ResultError, ResultFile, IndexRecord};

/// First-pass scan result file
pub struct FpScanResultFile {
    path: PathBuf,
    writer: Option<BufWriter<File>>,
    reader: Option<BufReader<File>>,
    indexes: Vec<IndexRecord>,
    current_seq_idx: usize,
    current_record: usize,
}

impl FpScanResultFile {
    pub fn new(path: &Path) -> Self {
        Self {
            path: path.to_path_buf(),
            writer: None,
            reader: None,
            indexes: Vec::new(),
            current_seq_idx: 0,
            current_record: 0,
        }
    }

    /// Open for writing
    pub fn open_write(&mut self) -> io::Result<()> {
        let file = File::create(&self.path)?;
        self.writer = Some(BufWriter::new(file));
        Ok(())
    }

    /// Open for appending
    pub fn open_append(&mut self) -> io::Result<()> {
        let file = File::options().append(true).open(&self.path)?;
        self.writer = Some(BufWriter::new(file));
        Ok(())
    }

    /// Open for reading
    pub fn open_read(&mut self) -> io::Result<()> {
        let file = File::open(&self.path)?;
        self.reader = Some(BufReader::new(file));
        Ok(())
    }

    /// Initialize first-pass result file with header
    pub fn init_with_header(&mut self) -> io::Result<()> {
        if self.writer.is_none() {
            self.open_append()?;
        }

        if let Some(ref mut w) = self.writer {
            writeln!(w, "Sequence\t\ttRNA Bounds\ttRNA\tAnti\t")?;
            writeln!(w, "Name     \ttRNA #\tBegin\tEnd\tType\tCodon\tSeqID\tSeqLen\tScore")?;
            writeln!(w, "--------\t------\t-----\t---\t----\t-----\t-----\t------\t-----")?;
        }
        Ok(())
    }

    /// Write a first-pass hit
    pub fn write_hit(&mut self, record: &FpScanRecord) -> io::Result<()> {
        if let Some(ref mut w) = self.writer {
            writeln!(w, "{}", record.to_line())?;
            Ok(())
        } else {
            Err(io::Error::new(io::ErrorKind::Other, "File not open for writing"))
        }
    }

    /// Index results by sequence name
    pub fn index_results(&mut self) -> io::Result<bool> {
        if self.reader.is_none() {
            self.open_read()?;
        }

        self.indexes.clear();
        let mut seqinfo_flag = false;

        if let Some(ref mut r) = self.reader {
            let mut file_pos = 0u64;
            let mut prev_seqname = String::new();
            let mut line_ct = 0usize;

            for line_result in r.lines() {
                let line = line_result?;

                // Check for column header with SeqID/SeqLen
                if line.contains("Type\tCodon\tSeqID\tSeqLen") {
                    seqinfo_flag = true;
                    file_pos += line.len() as u64 + 1;
                    continue;
                }

                // Skip other header lines
                if line.starts_with("Sequence") || line.starts_with("Name") || line.starts_with("---") {
                    file_pos += line.len() as u64 + 1;
                    continue;
                }

                // Parse data line
                let columns: Vec<&str> = line.split_whitespace().collect();
                if columns.len() >= 9 {
                    let seqname = columns[0].to_string();

                    if seqname != prev_seqname {
                        // New sequence - save previous count
                        if !prev_seqname.is_empty() && !self.indexes.is_empty() {
                            let last_idx = self.indexes.len() - 1;
                            self.indexes[last_idx].data.push(line_ct.to_string());
                        }

                        // Start new index
                        self.indexes.push(IndexRecord::with_data(
                            file_pos,
                            vec![seqname.clone()],
                        ));

                        prev_seqname = seqname;
                        line_ct = 0;
                    }
                    line_ct += 1;
                }

                file_pos += line.len() as u64 + 1;
            }

            // Save final count
            if !prev_seqname.is_empty() && !self.indexes.is_empty() {
                let last_idx = self.indexes.len() - 1;
                self.indexes[last_idx].data.push(line_ct.to_string());
            }
        }

        Ok(seqinfo_flag)
    }

    /// Get next tRNA candidate from indexed file
    pub fn get_next_candidate(&mut self, seq_idx: usize) -> Result<Option<FpScanRecord>, ResultError> {
        if seq_idx >= self.indexes.len() {
            return Ok(None);
        }

        // Reset to new sequence if needed
        if self.current_seq_idx != seq_idx {
            self.current_seq_idx = seq_idx;
            self.current_record = 0;

            if self.reader.is_none() {
                self.open_read()?;
            }
        }

        let record = &self.indexes[seq_idx];
        let max_records: usize = record.data.get(1)
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);

        if self.current_record >= max_records {
            return Ok(None);
        }

        // Read next line from current position
        if let Some(ref mut r) = self.reader {
            let mut line = String::new();
            r.read_line(&mut line)?;

            self.current_record += 1;

            if line.trim().is_empty() {
                return Ok(None);
            }

            Ok(Some(FpScanRecord::from_line(&line)?))
        } else {
            Ok(None)
        }
    }

    /// Read all records
    pub fn read_all(&mut self) -> Result<Vec<FpScanRecord>, ResultError> {
        if self.reader.is_none() {
            self.open_read()?;
        }

        let mut records = Vec::new();
        if let Some(ref mut r) = self.reader {
            for line_result in r.lines() {
                let line = line_result?;

                // Skip headers
                if line.starts_with("Sequence") || line.starts_with("Name") || line.starts_with("---") {
                    continue;
                }

                if line.trim().is_empty() {
                    continue;
                }

                match FpScanRecord::from_line(&line) {
                    Ok(record) => records.push(record),
                    Err(e) => {
                        eprintln!("Warning: Skipping invalid record: {}", e);
                        continue;
                    }
                }
            }
        }

        Ok(records)
    }

    /// Reset current sequence tracking
    pub fn reset_current_seq(&mut self) {
        self.current_seq_idx = 0;
        self.current_record = 0;
    }

    /// Get hit count across all sequences
    pub fn get_hit_count(&self) -> usize {
        self.indexes.iter()
            .filter_map(|idx| idx.data.get(1))
            .filter_map(|s| s.parse::<usize>().ok())
            .sum()
    }

    /// Clear indexes
    pub fn clear_indexes(&mut self) {
        self.indexes.clear();
        self.reset_current_seq();
    }
}

impl ResultFile for FpScanResultFile {
    fn open(path: &Path) -> Result<Self, ResultError> {
        let mut file = Self::new(path);
        file.open_read()?;
        Ok(file)
    }

    fn write_header(&mut self) -> io::Result<()> {
        self.init_with_header()
    }

    fn close(&mut self) -> io::Result<()> {
        if let Some(ref mut w) = self.writer {
            w.flush()?;
        }
        self.writer = None;
        self.reader = None;
        Ok(())
    }
}

// ============================================================================
// First-Pass Scan Record
// ============================================================================

/// A single first-pass hit record
#[derive(Debug, Clone, Default)]
pub struct FpScanRecord {
    pub seqname: String,
    pub trna_num: i32,
    pub start: i64,
    pub end: i64,
    pub isotype: String,
    pub anticodon: String,
    pub seq_id: i32,
    pub seq_len: i64,
    pub score: f64,
    pub model: String,
    pub hit_source: String,
}

impl FpScanRecord {
    pub fn to_line(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{}",
            self.seqname,
            self.trna_num,
            self.start,
            self.end,
            self.isotype,
            self.anticodon,
            self.seq_id,
            self.seq_len,
            self.score,
            self.model
        )
    }

    pub fn from_line(line: &str) -> Result<Self, ResultError> {
        let columns: Vec<&str> = line.split_whitespace().collect();
        if columns.len() < 9 {
            return Err(ResultError::Parse(format!(
                "Expected at least 9 columns, got {}", columns.len()
            )));
        }

        Ok(Self {
            seqname: columns[0].to_string(),
            trna_num: columns[1].parse().unwrap_or(0),
            start: columns[2].parse().unwrap_or(0),
            end: columns[3].parse().unwrap_or(0),
            isotype: columns[4].to_string(),
            anticodon: columns[5].to_string(),
            seq_id: columns[6].parse().unwrap_or(0),
            seq_len: columns[7].parse().unwrap_or(0),
            score: columns[8].parse().unwrap_or(0.0),
            model: columns.get(9).unwrap_or(&"").to_string(),
            hit_source: columns.get(10).unwrap_or(&"").to_string(),
        })
    }
}
