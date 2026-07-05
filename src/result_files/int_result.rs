//! Intermediate result file module
//!
//! Corresponds to tRNAscanSE::IntResultFile.pm
//!
//! This file format stores detailed tRNA information between pipeline stages:
//! - Full tRNA structure data
//! - All score components (HMM, SS, mat, cove, infernal)
//! - Intron and anticodon positions
//! - Sequences and secondary structure
//!
//! Format: Tab-delimited with header line containing 40+ columns

use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write, Seek, SeekFrom};
use std::path::{Path, PathBuf};

use super::{ResultError, ResultFile, IndexRecord, parse_intron, format_introns, parse_ac_pos, format_ac_pos, parse_opt_f64};

/// Intermediate result file writer/reader
pub struct IntResultFile {
    path: PathBuf,
    writer: Option<BufWriter<File>>,
    reader: Option<BufReader<File>>,
    indexes: Vec<IndexRecord>,
}

impl IntResultFile {
    /// Create a new intermediate result file
    pub fn new(path: &Path) -> Self {
        Self {
            path: path.to_path_buf(),
            writer: None,
            reader: None,
            indexes: Vec::new(),
        }
    }

    /// Open for writing
    pub fn open_write(&mut self) -> io::Result<()> {
        let file = File::create(&self.path)?;
        self.writer = Some(BufWriter::new(file));
        self.write_header()?;
        Ok(())
    }

    /// Open for reading
    pub fn open_read(&mut self) -> io::Result<()> {
        let file = File::open(&self.path)?;
        self.reader = Some(BufReader::new(file));
        Ok(())
    }

    /// Open for appending
    pub fn open_append(&mut self) -> io::Result<()> {
        let file = File::options().append(true).open(&self.path)?;
        self.writer = Some(BufWriter::new(file));
        Ok(())
    }

    /// Write a tRNA record
    pub fn write_trna(&mut self, record: &IntResultRecord) -> io::Result<()> {
        if let Some(ref mut w) = self.writer {
            writeln!(w, "{}", record.to_line())?;
            Ok(())
        } else {
            Err(io::Error::new(io::ErrorKind::Other, "File not open for writing"))
        }
    }

    /// Read all tRNA records
    pub fn read_all(&mut self) -> Result<Vec<IntResultRecord>, ResultError> {
        if self.reader.is_none() {
            self.open_read()?;
        }

        let mut records = Vec::new();
        if let Some(ref mut r) = self.reader {
            let mut line_num = 0;
            for line_result in r.lines() {
                line_num += 1;
                let line = line_result?;

                // Skip header
                if line_num == 1 {
                    continue;
                }

                match IntResultRecord::from_line(&line) {
                    Ok(record) => records.push(record),
                    Err(e) => {
                        eprintln!("Warning: Skipping invalid record at line {}: {}", line_num, e);
                        continue;
                    }
                }
            }
        }

        Ok(records)
    }

    /// Build index by tRNAscan_id
    pub fn build_index(&mut self) -> io::Result<()> {
        if self.reader.is_none() {
            self.open_read()?;
        }

        self.indexes.clear();

        if let Some(ref mut r) = self.reader {
            let mut file_pos = 0u64;
            let mut line_num = 0;

            for line_result in r.lines() {
                line_num += 1;
                let line = line_result?;

                // Skip header
                if line_num == 1 {
                    file_pos += line.len() as u64 + 1; // +1 for newline
                    continue;
                }

                let columns: Vec<&str> = line.split('\t').collect();
                if columns.len() >= 4 {
                    let trna_id = columns[3].to_string();
                    self.indexes.push(IndexRecord::with_data(file_pos, vec![trna_id]));
                }

                file_pos += line.len() as u64 + 1;
            }

            // Sort by tRNAscan_id
            self.indexes.sort_by(|a, b| a.data[0].cmp(&b.data[0]));
        }

        Ok(())
    }

    /// Get record at file position
    pub fn get_at(&mut self, file_pos: u64) -> Result<IntResultRecord, ResultError> {
        if self.reader.is_none() {
            self.open_read()?;
        }

        if let Some(ref mut r) = self.reader {
            r.seek(SeekFrom::Start(file_pos))?;
            let mut line = String::new();
            r.read_line(&mut line)?;
            IntResultRecord::from_line(&line)
        } else {
            Err(ResultError::Io(io::Error::new(io::ErrorKind::Other, "Reader not open")))
        }
    }

    /// Clear indexes
    pub fn clear_indexes(&mut self) {
        self.indexes.clear();
    }

    /// Get index count
    pub fn index_count(&self) -> usize {
        self.indexes.len()
    }
}

impl ResultFile for IntResultFile {
    fn open(path: &Path) -> Result<Self, ResultError> {
        let mut file = Self::new(path);
        file.open_read()?;
        Ok(file)
    }

    fn write_header(&mut self) -> io::Result<()> {
        if let Some(ref mut w) = self.writer {
            writeln!(w, "seqname\tordered_seqname\tid\ttRNAscan_id\tstart\tend\tstrand\t\
                         start2\tend2\tstrand2\tstart3\tend3\tstrand3\tposition\t\
                         isotype\tanticodon\tac_pos\tintron\tmodel\tscore\tmat_score\t\
                         hmm_score\tss_score\tcove_score\tcove_mat_score\tcove_hmm_score\t\
                         cove_ss_score\tinf_score\tinf_mat_score\tinf_hmm_score\t\
                         inf_ss_score\tpseudo\ttrunc\tcategory\thit_source\t\
                         src_seqid\tsrc_seq_len\tseq\tmat_seq\tss\tmat_ss\tnote")?;
        }
        Ok(())
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
// Intermediate Result Record
// ============================================================================

/// A single tRNA record from intermediate results
#[derive(Debug, Clone, Default)]
pub struct IntResultRecord {
    pub seqname: String,
    pub ordered_seqname: String,
    pub id: i32,
    pub trnascan_id: String,
    pub start: i64,
    pub end: i64,
    pub strand: String,
    // Exon 2 (for introns)
    pub start2: i64,
    pub end2: i64,
    pub strand2: String,
    // Exon 3 (for rare cases with 2 introns)
    pub start3: i64,
    pub end3: i64,
    pub strand3: String,
    pub position: String,
    pub isotype: String,
    pub anticodon: String,
    pub ac_pos: Vec<(i32, i32)>,
    pub introns: Vec<(i32, i32, i64, i64, String, String)>,
    pub model: String,
    // Scores
    pub score: f64,
    pub mat_score: f64,
    pub hmm_score: f64,
    pub ss_score: f64,
    // Cove scores (optional)
    pub cove_score: Option<f64>,
    pub cove_mat_score: Option<f64>,
    pub cove_hmm_score: Option<f64>,
    pub cove_ss_score: Option<f64>,
    // Infernal scores (optional)
    pub inf_score: Option<f64>,
    pub inf_mat_score: Option<f64>,
    pub inf_hmm_score: Option<f64>,
    pub inf_ss_score: Option<f64>,
    pub pseudo: i32,
    pub trunc: String,
    pub category: String,
    pub hit_source: String,
    pub src_seqid: i32,
    pub src_seqlen: i64,
    pub seq: String,
    pub mat_seq: String,
    pub ss: String,
    pub mat_ss: String,
    pub note: String,
}

impl IntResultRecord {
    pub fn to_line(&self) -> String {
        let ac_pos_str = format_ac_pos(&self.ac_pos);
        let intron_str = format_introns(&self.introns);

        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t\
             {}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.seqname, self.ordered_seqname, self.id, self.trnascan_id,
            self.start, self.end, self.strand,
            self.start2, self.end2, self.strand2,
            self.start3, self.end3, self.strand3,
            self.position, self.isotype, self.anticodon,
            ac_pos_str, intron_str, self.model,
            self.score, self.mat_score, self.hmm_score, self.ss_score,
            self.cove_score.map_or(String::new(), |s| s.to_string()),
            self.cove_mat_score.map_or(String::new(), |s| s.to_string()),
            self.cove_hmm_score.map_or(String::new(), |s| s.to_string()),
            self.cove_ss_score.map_or(String::new(), |s| s.to_string()),
            self.inf_score.map_or(String::new(), |s| s.to_string()),
            self.inf_mat_score.map_or(String::new(), |s| s.to_string()),
            self.inf_hmm_score.map_or(String::new(), |s| s.to_string()),
            self.inf_ss_score.map_or(String::new(), |s| s.to_string()),
            self.pseudo, self.trunc, self.category, self.hit_source,
            self.src_seqid, self.src_seqlen,
            self.seq, self.mat_seq, self.ss, self.mat_ss, self.note
        )
    }

    pub fn from_line(line: &str) -> Result<Self, ResultError> {
        let columns: Vec<&str> = line.split('\t').collect();
        if columns.len() < 42 {
            return Err(ResultError::Parse(format!(
                "Expected at least 42 columns, got {}", columns.len()
            )));
        }

        let ac_pos = if !columns[16].is_empty() {
            parse_ac_pos(columns[16])?
        } else {
            Vec::new()
        };

        let introns = if !columns[17].is_empty() {
            parse_intron(columns[17])?
        } else {
            Vec::new()
        };

        Ok(Self {
            seqname: columns[0].to_string(),
            ordered_seqname: columns[1].to_string(),
            id: columns[2].parse().unwrap_or(0),
            trnascan_id: columns[3].to_string(),
            start: columns[4].parse().unwrap_or(0),
            end: columns[5].parse().unwrap_or(0),
            strand: columns[6].to_string(),
            start2: columns[7].parse().unwrap_or(0),
            end2: columns[8].parse().unwrap_or(0),
            strand2: columns[9].to_string(),
            start3: columns[10].parse().unwrap_or(0),
            end3: columns[11].parse().unwrap_or(0),
            strand3: columns[12].to_string(),
            position: columns[13].to_string(),
            isotype: columns[14].to_string(),
            anticodon: columns[15].to_string(),
            ac_pos,
            introns,
            model: columns[18].to_string(),
            score: columns[19].parse().unwrap_or(0.0),
            mat_score: columns[20].parse().unwrap_or(0.0),
            hmm_score: columns[21].parse().unwrap_or(0.0),
            ss_score: columns[22].parse().unwrap_or(0.0),
            cove_score: parse_opt_f64(columns[23]),
            cove_mat_score: parse_opt_f64(columns[24]),
            cove_hmm_score: parse_opt_f64(columns[25]),
            cove_ss_score: parse_opt_f64(columns[26]),
            inf_score: parse_opt_f64(columns[27]),
            inf_mat_score: parse_opt_f64(columns[28]),
            inf_hmm_score: parse_opt_f64(columns[29]),
            inf_ss_score: parse_opt_f64(columns[30]),
            pseudo: columns[31].parse().unwrap_or(0),
            trunc: columns[32].to_string(),
            category: columns[33].to_string(),
            hit_source: columns[34].to_string(),
            src_seqid: columns[35].parse().unwrap_or(0),
            src_seqlen: columns[36].parse().unwrap_or(0),
            seq: columns[37].to_string(),
            mat_seq: columns.get(38).unwrap_or(&"").to_string(),
            ss: columns.get(39).unwrap_or(&"").to_string(),
            mat_ss: columns.get(40).unwrap_or(&"").to_string(),
            note: columns.get(41).unwrap_or(&"").to_string(),
        })
    }
}
