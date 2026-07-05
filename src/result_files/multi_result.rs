//! Multi-model result file module
//!
//! Corresponds to tRNAscanSE::MultiResultFile.pm
//!
//! Stores results from scanning with multiple covariance models.
//! Each tRNA has scores from different models (e.g., cyto vs mito models).
//!
//! Format: Tab-delimited with first line containing model names

use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write, Seek, SeekFrom};
use std::path::{Path, PathBuf};

use super::{ResultError, ResultFile, IndexRecord};

/// Multi-model result file
pub struct MultiResultFile {
    path: PathBuf,
    writer: Option<BufWriter<File>>,
    reader: Option<BufReader<File>>,
    models: Vec<String>,
    indexes: Vec<IndexRecord>,
}

impl MultiResultFile {
    pub fn new(path: &Path) -> Self {
        Self {
            path: path.to_path_buf(),
            writer: None,
            reader: None,
            models: Vec::new(),
            indexes: Vec::new(),
        }
    }

    /// Open for writing
    pub fn open_write(&mut self) -> io::Result<()> {
        let file = File::create(&self.path)?;
        self.writer = Some(BufWriter::new(file));
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

    /// Write a line to the file
    pub fn write_line(&mut self, line: &str) -> io::Result<()> {
        if let Some(ref mut w) = self.writer {
            writeln!(w, "{}", line)?;
            Ok(())
        } else {
            Err(io::Error::new(io::ErrorKind::Other, "File not open for writing"))
        }
    }

    /// Read a line from the file
    pub fn get_line(&mut self) -> io::Result<Option<String>> {
        if let Some(ref mut r) = self.reader {
            let mut line = String::new();
            match r.read_line(&mut line) {
                Ok(0) => Ok(None),
                Ok(_) => {
                    line = line.trim_end().to_string();
                    Ok(Some(line))
                }
                Err(e) => Err(e),
            }
        } else {
            Err(io::Error::new(io::ErrorKind::Other, "File not open for reading"))
        }
    }

    /// Read model names from first line
    pub fn read_models(&mut self) -> io::Result<()> {
        if let Some(line) = self.get_line()? {
            let columns: Vec<&str> = line.split('\t').collect();
            self.models.clear();
            // Skip first column (tRNAscan_id header)
            for col in columns.iter().skip(1) {
                self.models.push(col.to_string());
            }
        }
        Ok(())
    }

    /// Get model names
    pub fn get_models(&self) -> &[String] {
        &self.models
    }

    /// Get next record
    pub fn get_next_record(&mut self) -> io::Result<Option<(String, String)>> {
        if let Some(line) = self.get_line()? {
            let columns: Vec<&str> = line.split('\t').collect();
            if !columns.is_empty() {
                let trna_id = columns[0].to_string();
                return Ok(Some((trna_id, line)));
            }
        }
        Ok(None)
    }

    /// Parse model hits from a record line
    pub fn parse_record(&self, line: &str) -> Result<Vec<MultiModelHit>, ResultError> {
        let columns: Vec<&str> = line.split('\t').collect();
        if columns.is_empty() {
            return Ok(Vec::new());
        }

        let mut hits = Vec::new();

        // Skip first column (tRNAscan_id)
        for (i, score_str) in columns.iter().skip(1).enumerate() {
            if score_str.is_empty() {
                continue;
            }

            if i >= self.models.len() {
                break;
            }

            let score: f64 = score_str.parse()
                .map_err(|_| ResultError::Parse(format!("Invalid score: {}", score_str)))?;

            let model_name = &self.models[i];
            let (domain, model) = if let Some(stripped) = model_name.strip_prefix("mito_") {
                ("mito", stripped)
            } else {
                ("cyto", model_name.as_str())
            };

            hits.push(MultiModelHit {
                domain: domain.to_string(),
                model: model.to_string(),
                score,
            });
        }

        Ok(hits)
    }

    /// Get record at file position
    pub fn get_at(&mut self, file_pos: u64) -> Result<Vec<MultiModelHit>, ResultError> {
        if self.reader.is_none() {
            self.open_read()?;
        }

        if let Some(ref mut r) = self.reader {
            r.seek(SeekFrom::Start(file_pos))?;
            if let Some(line) = self.get_line()? {
                return self.parse_record(&line);
            }
        }

        Ok(Vec::new())
    }

    /// Build index by tRNAscan_id
    pub fn build_index(&mut self) -> io::Result<()> {
        if self.reader.is_none() {
            self.open_read()?;
        }

        self.indexes.clear();

        if let Some(ref mut r) = self.reader {
            // Skip header line
            let mut line = String::new();
            let mut file_pos = 0u64;
            r.read_line(&mut line)?;
            file_pos += line.len() as u64;

            // Read records
            for line_result in r.lines() {
                let line = line_result?;
                let columns: Vec<&str> = line.split('\t').collect();

                if !columns.is_empty() {
                    let trna_id = columns[0].to_string();
                    self.indexes.push(IndexRecord::with_data(
                        file_pos,
                        vec![trna_id],
                    ));
                }

                file_pos += line.len() as u64 + 1;
            }

            // Sort by tRNAscan_id
            self.indexes.sort_by(|a, b| a.data[0].cmp(&b.data[0]));
        }

        Ok(())
    }

    /// Binary search for tRNAscan_id
    pub fn bsearch_id(&self, id: &str) -> Option<usize> {
        let mut left = 0;
        let mut right = self.indexes.len();

        while left < right {
            let mid = (left + right) / 2;
            match self.indexes[mid].data[0].as_str().cmp(id) {
                std::cmp::Ordering::Less => left = mid + 1,
                std::cmp::Ordering::Greater => right = mid,
                std::cmp::Ordering::Equal => return Some(mid),
            }
        }

        None
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

impl ResultFile for MultiResultFile {
    fn open(path: &Path) -> Result<Self, ResultError> {
        let mut file = Self::new(path);
        file.open_read()?;
        file.read_models()?;
        Ok(file)
    }

    fn write_header(&mut self) -> io::Result<()> {
        // Header is written dynamically with model names
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
// Multi-Model Hit Record
// ============================================================================

/// A model hit from multi-model scanning
#[derive(Debug, Clone)]
pub struct MultiModelHit {
    pub domain: String,  // "cyto" or "mito"
    pub model: String,   // model name
    pub score: f64,      // model score
}

impl Default for MultiModelHit {
    fn default() -> Self {
        Self {
            domain: String::new(),
            model: String::new(),
            score: 0.0,
        }
    }
}
