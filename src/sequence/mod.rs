//! Sequence handling module for tRNAscan-SE
//!
//! This module provides sequence reading, manipulation, and writing capabilities,
//! with support for FASTA format files and efficient handling of large sequences.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Errors that can occur during sequence operations
#[derive(Debug)]
pub enum SeqError {
    Io(std::io::Error),
    Parse(String),
    InvalidFormat(String),
    NotFound(String),
}

impl std::fmt::Display for SeqError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SeqError::Io(e) => write!(f, "I/O error: {}", e),
            SeqError::Parse(s) => write!(f, "Parse error: {}", s),
            SeqError::InvalidFormat(s) => write!(f, "Invalid format: {}", s),
            SeqError::NotFound(s) => write!(f, "Not found: {}", s),
        }
    }
}

impl std::error::Error for SeqError {}

impl From<std::io::Error> for SeqError {
    fn from(error: std::io::Error) -> Self {
        SeqError::Io(error)
    }
}

/// Sequence format
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SeqFormat {
    Fasta,
    Raw,
}

/// A sequence with associated metadata
#[derive(Debug, Clone)]
pub struct Sequence {
    pub name: String,
    pub description: String,
    pub data: Vec<u8>,
    pub length: usize,
    pub format: SeqFormat,
}

impl Sequence {
    /// Create a new sequence
    pub fn new(name: &str, data: Vec<u8>) -> Self {
        let length = data.len();
        Sequence {
            name: name.to_string(),
            description: String::new(),
            data,
            length,
            format: SeqFormat::Raw,
        }
    }

    /// Create a sequence from FASTA format string
    pub fn from_fasta(content: &str) -> Result<Self, SeqError> {
        let mut lines = content.lines();

        let header = lines.next()
            .ok_or_else(|| SeqError::Parse("Empty FASTA content".to_string()))?;

        if !header.starts_with('>') {
            return Err(SeqError::InvalidFormat(
                "FASTA must start with '>'".to_string()
            ));
        }

        let header = &header[1..]; // Remove '>'
        let parts: Vec<&str> = header.splitn(2, ' ').collect();
        let name = parts[0].to_string();
        let description = parts.get(1).unwrap_or(&"").to_string();

        let mut data = Vec::new();
        for line in lines {
            let line = line.trim();
            if !line.is_empty() && !line.starts_with('>') {
                // Strip whitespace and digits
                for ch in line.chars() {
                    if ch.is_ascii_alphabetic() {
                        data.push(ch.to_ascii_uppercase() as u8);
                    }
                }
            }
        }

        let length = data.len();
        Ok(Sequence {
            name,
            description,
            data,
            length,
            format: SeqFormat::Fasta,
        })
    }

    /// Read a sequence from a file
    pub fn from_file(path: &Path) -> Result<Self, SeqError> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let format = detect_format(path);

        match format {
            SeqFormat::Fasta => {
                let mut content = String::new();
                for line in reader.lines() {
                    content.push_str(&line?);
                    content.push('\n');
                }
                Self::from_fasta(&content)
            }
            SeqFormat::Raw => {
                let mut data = Vec::new();
                for line in reader.lines() {
                    let line = line?;
                    for ch in line.chars() {
                        if ch.is_ascii_alphabetic() {
                            data.push(ch.to_ascii_uppercase() as u8);
                        }
                    }
                }
                let length = data.len();
                Ok(Sequence {
                    name: path.file_stem()
                        .and_then(|s| s.to_str())
                        .unwrap_or("unknown")
                        .to_string(),
                    description: String::new(),
                    data,
                    length,
                    format: SeqFormat::Raw,
                })
            }
        }
    }

    /// Extract a subsequence
    pub fn subseq(&self, start: usize, end: usize) -> Vec<u8> {
        if start >= self.length || end > self.length || start >= end {
            return Vec::new();
        }
        self.data[start..end].to_vec()
    }

    /// Get reverse complement of the sequence
    pub fn reverse_complement(&self) -> Self {
        let mut rc_data = Vec::with_capacity(self.length);

        for &base in self.data.iter().rev() {
            let complement = match base {
                b'A' => b'T',
                b'T' => b'A',
                b'G' => b'C',
                b'C' => b'G',
                b'U' => b'A',
                b'N' => b'N',
                _ => base,
            };
            rc_data.push(complement);
        }

        Sequence {
            name: format!("{}_rc", self.name),
            description: format!("Reverse complement of {}", self.description),
            data: rc_data,
            length: self.length,
            format: self.format,
        }
    }

    /// Calculate GC content as a fraction (0.0 to 1.0)
    pub fn gc_content(&self) -> f64 {
        if self.length == 0 {
            return 0.0;
        }

        let gc_count = self.data.iter()
            .filter(|&&b| b == b'G' || b == b'C')
            .count();

        gc_count as f64 / self.length as f64
    }

    /// Convert to RNA (T -> U)
    pub fn to_rna(&self) -> Self {
        let rna_data: Vec<u8> = self.data.iter()
            .map(|&b| if b == b'T' { b'U' } else { b })
            .collect();

        Sequence {
            name: self.name.clone(),
            description: self.description.clone(),
            data: rna_data,
            length: self.length,
            format: self.format,
        }
    }

    /// Convert to DNA (U -> T)
    pub fn to_dna(&self) -> Self {
        let dna_data: Vec<u8> = self.data.iter()
            .map(|&b| if b == b'U' { b'T' } else { b })
            .collect();

        Sequence {
            name: self.name.clone(),
            description: self.description.clone(),
            data: dna_data,
            length: self.length,
            format: self.format,
        }
    }

    /// Clean sequence by removing non-ACGT bases (convert to N)
    pub fn clean(&mut self) {
        for base in &mut self.data {
            if !matches!(*base, b'A' | b'C' | b'G' | b'T' | b'U' | b'N') {
                *base = b'N';
            }
        }
    }

    /// Convert sequence to FASTA format (single line)
    pub fn to_fasta(&self) -> String {
        let seq_str = String::from_utf8_lossy(&self.data);
        if self.description.is_empty() {
            format!(">{}\n{}\n", self.name, seq_str)
        } else {
            format!(">{} {}\n{}\n", self.name, self.description, seq_str)
        }
    }

    /// Convert sequence to FASTA format with line wrapping
    pub fn to_fasta_wrapped(&self, width: usize) -> String {
        let mut result = String::new();

        if self.description.is_empty() {
            result.push_str(&format!(">{}\n", self.name));
        } else {
            result.push_str(&format!(">{} {}\n", self.name, self.description));
        }

        let seq_str = String::from_utf8_lossy(&self.data);
        let mut pos = 0;
        while pos < self.length {
            let end = (pos + width).min(self.length);
            result.push_str(&seq_str[pos..end]);
            result.push('\n');
            pos = end;
        }

        result
    }
}

/// Detect the format of a sequence file
pub fn detect_format(path: &Path) -> SeqFormat {
    if let Ok(file) = File::open(path) {
        let mut reader = BufReader::new(file);
        let mut first_line = String::new();
        if reader.read_line(&mut first_line).is_ok() {
            if first_line.starts_with('>') {
                return SeqFormat::Fasta;
            }
        }
    }
    SeqFormat::Raw
}

/// Read all sequences from a file
pub fn read_sequences(path: &Path) -> Result<Vec<Sequence>, SeqError> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let format = detect_format(path);

    match format {
        SeqFormat::Fasta => {
            let mut sequences = Vec::new();
            let mut current_name = String::new();
            let mut current_desc = String::new();
            let mut current_data = Vec::new();

            for line in reader.lines() {
                let line = line?;
                let trimmed = line.trim();

                if trimmed.starts_with('>') {
                    // Save previous sequence if any
                    if !current_name.is_empty() {
                        let length = current_data.len();
                        sequences.push(Sequence {
                            name: current_name,
                            description: current_desc,
                            data: current_data,
                            length,
                            format: SeqFormat::Fasta,
                        });
                    }

                    // Start new sequence
                    let header = &trimmed[1..];
                    let parts: Vec<&str> = header.splitn(2, ' ').collect();
                    current_name = parts[0].to_string();
                    current_desc = parts.get(1).unwrap_or(&"").to_string();
                    current_data = Vec::new();
                } else if !trimmed.is_empty() {
                    // Accumulate sequence data
                    for ch in trimmed.chars() {
                        if ch.is_ascii_alphabetic() {
                            current_data.push(ch.to_ascii_uppercase() as u8);
                        }
                    }
                }
            }

            // Save last sequence
            if !current_name.is_empty() {
                let length = current_data.len();
                sequences.push(Sequence {
                    name: current_name,
                    description: current_desc,
                    data: current_data,
                    length,
                    format: SeqFormat::Fasta,
                });
            }

            Ok(sequences)
        }
        SeqFormat::Raw => {
            // For raw format, treat entire file as single sequence
            let seq = Sequence::from_file(path)?;
            Ok(vec![seq])
        }
    }
}

/// Read a single sequence from a file
pub fn read_sequence(path: &Path) -> Result<Sequence, SeqError> {
    Sequence::from_file(path)
}

/// Iterator for reading sequences from large files
pub struct SequenceReader {
    reader: BufReader<File>,
    saved_line: Option<String>,
    seq_id: usize,
    finished: bool,
}

impl SequenceReader {
    /// Create a new sequence reader
    pub fn new(path: &Path) -> Result<Self, SeqError> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);

        Ok(SequenceReader {
            reader,
            saved_line: None,
            seq_id: 0,
            finished: false,
        })
    }
}

impl Iterator for SequenceReader {
    type Item = Result<Sequence, SeqError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }

        let mut name = String::new();
        let mut description = String::new();
        let mut data = Vec::new();
        let mut found_header = false;

        // Process saved line from previous iteration
        if let Some(line) = self.saved_line.take() {
            if line.starts_with('>') {
                let header = &line[1..];
                let parts: Vec<&str> = header.splitn(2, ' ').collect();
                name = parts[0].trim().to_string();
                description = parts.get(1).unwrap_or(&"").trim().to_string();
                found_header = true;
                self.seq_id += 1;
            }
        }

        // Read lines until we find a header or EOF
        loop {
            let mut line = String::new();
            match self.reader.read_line(&mut line) {
                Ok(0) => {
                    // EOF
                    self.finished = true;
                    if found_header {
                        let length = data.len();
                        return Some(Ok(Sequence {
                            name,
                            description,
                            data,
                            length,
                            format: SeqFormat::Fasta,
                        }));
                    }
                    return None;
                }
                Ok(_) => {
                    let trimmed = line.trim();
                    if trimmed.starts_with('>') {
                        if !found_header {
                            // First header
                            let header = &trimmed[1..];
                            let parts: Vec<&str> = header.splitn(2, ' ').collect();
                            name = parts[0].to_string();
                            description = parts.get(1).unwrap_or(&"").to_string();
                            found_header = true;
                            self.seq_id += 1;
                        } else {
                            // Next sequence header - save it and return current sequence
                            self.saved_line = Some(trimmed.to_string());
                            let length = data.len();
                            return Some(Ok(Sequence {
                                name,
                                description,
                                data,
                                length,
                                format: SeqFormat::Fasta,
                            }));
                        }
                    } else if found_header && !trimmed.is_empty() {
                        // Accumulate sequence data
                        for ch in trimmed.chars() {
                            if ch.is_ascii_alphabetic() {
                                data.push(ch.to_ascii_uppercase() as u8);
                            }
                        }
                    }
                }
                Err(e) => {
                    self.finished = true;
                    return Some(Err(SeqError::Io(e)));
                }
            }
        }
    }
}

/// Reader for handling large sequences with buffering
#[allow(dead_code)]
pub struct BufferedSequenceReader {
    reader: BufReader<File>,
    max_buffer_size: usize,
    buffer_overlap: usize,
    seq_name_map: HashMap<String, usize>,
    current_seq_id: usize,
}

impl BufferedSequenceReader {
    /// Create a new buffered sequence reader
    pub fn new(path: &Path) -> Result<Self, SeqError> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);

        Ok(BufferedSequenceReader {
            reader,
            max_buffer_size: 50_000_000, // 50 MB
            buffer_overlap: 200,
            seq_name_map: HashMap::new(),
            current_seq_id: 0,
        })
    }

    /// Set maximum buffer size for reading sequences
    pub fn set_max_buffer_size(&mut self, size: usize) {
        self.max_buffer_size = size;
    }

    /// Set overlap size between buffers
    pub fn set_buffer_overlap(&mut self, overlap: usize) {
        self.buffer_overlap = overlap;
    }
}

/// Reverse complement a sequence (raw function for utility)
pub fn rev_comp_seq(seq: &[u8]) -> Vec<u8> {
    let mut rc = Vec::with_capacity(seq.len());

    for &base in seq.iter().rev() {
        let complement = match base {
            b'A' => b'T',
            b'T' => b'A',
            b'G' => b'C',
            b'C' => b'G',
            b'U' => b'A',
            b'N' => b'N',
            _ => base,
        };
        rc.push(complement);
    }

    rc
}

/// Write sequences to a FASTA file
pub fn write_fasta(path: &Path, sequences: &[Sequence]) -> Result<(), SeqError> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    for seq in sequences {
        writer.write_all(seq.to_fasta_wrapped(60).as_bytes())?;
    }

    writer.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sequence_new() {
        let seq = Sequence::new("test", vec![b'A', b'C', b'G', b'T']);
        assert_eq!(seq.name, "test");
        assert_eq!(seq.length, 4);
        assert_eq!(seq.data, vec![b'A', b'C', b'G', b'T']);
    }

    #[test]
    fn test_from_fasta() {
        let fasta = ">seq1 description\nACGTACGT\nACGT\n";
        let seq = Sequence::from_fasta(fasta).unwrap();
        assert_eq!(seq.name, "seq1");
        assert_eq!(seq.description, "description");
        assert_eq!(seq.length, 12);
    }

    #[test]
    fn test_reverse_complement() {
        let seq = Sequence::new("test", vec![b'A', b'C', b'G', b'T']);
        let rc = seq.reverse_complement();
        assert_eq!(rc.data, vec![b'A', b'C', b'G', b'T']);

        let seq2 = Sequence::new("test2", vec![b'A', b'T', b'G', b'C']);
        let rc2 = seq2.reverse_complement();
        assert_eq!(rc2.data, vec![b'G', b'C', b'A', b'T']);
    }

    #[test]
    fn test_gc_content() {
        let seq = Sequence::new("test", vec![b'G', b'C', b'A', b'T']);
        assert!((seq.gc_content() - 0.5).abs() < 0.001);

        let seq2 = Sequence::new("test2", vec![b'G', b'G', b'C', b'C']);
        assert!((seq2.gc_content() - 1.0).abs() < 0.001);
    }

    #[test]
    fn test_subseq() {
        let seq = Sequence::new("test", vec![b'A', b'C', b'G', b'T', b'A', b'C']);
        let subseq = seq.subseq(1, 4);
        assert_eq!(subseq, vec![b'C', b'G', b'T']);
    }

    #[test]
    fn test_to_rna() {
        let seq = Sequence::new("test", vec![b'A', b'T', b'G', b'C']);
        let rna = seq.to_rna();
        assert_eq!(rna.data, vec![b'A', b'U', b'G', b'C']);
    }

    #[test]
    fn test_to_dna() {
        let seq = Sequence::new("test", vec![b'A', b'U', b'G', b'C']);
        let dna = seq.to_dna();
        assert_eq!(dna.data, vec![b'A', b'T', b'G', b'C']);
    }

    #[test]
    fn test_clean() {
        let mut seq = Sequence::new("test", vec![b'A', b'X', b'G', b'Z', b'T']);
        seq.clean();
        assert_eq!(seq.data, vec![b'A', b'N', b'G', b'N', b'T']);
    }

    #[test]
    fn test_to_fasta() {
        let seq = Sequence::new("test", vec![b'A', b'C', b'G', b'T']);
        let fasta = seq.to_fasta();
        assert_eq!(fasta, ">test\nACGT\n");
    }
}
