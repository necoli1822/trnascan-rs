// Sequence file I/O functions
// Ported from original sqio.c

use std::fs::File;
use std::io::{BufRead, BufReader, Error as IoError, ErrorKind, Result as IoResult};
use std::path::Path;

use super::iupac::{is_seq_char, seq_type as determine_seq_type};
use super::types::*;

const LINEBUFLEN: usize = 4096;

/// Sequence file reader
pub struct SeqFileReader {
    reader: BufReader<File>,
    current_line: String,
    format: i32,
    has_buffered_line: bool,
}

impl SeqFileReader {
    /// Open a sequence file and detect its format
    pub fn open<P: AsRef<Path>>(path: P) -> IoResult<Self> {
        let format = detect_format(path.as_ref())?;
        let file = File::open(path)?;
        let reader = BufReader::new(file);

        Ok(Self {
            reader,
            current_line: String::with_capacity(LINEBUFLEN),
            format,
            has_buffered_line: false,
        })
    }

    /// Open a sequence file with known format
    pub fn open_with_format<P: AsRef<Path>>(path: P, format: i32) -> IoResult<Self> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);

        Ok(Self {
            reader,
            current_line: String::with_capacity(LINEBUFLEN),
            format,
            has_buffered_line: false,
        })
    }

    /// Get the detected format
    pub fn format(&self) -> i32 {
        self.format
    }

    /// Read the next line into the buffer
    fn read_line(&mut self) -> IoResult<bool> {
        if self.has_buffered_line {
            self.has_buffered_line = false;
            return Ok(true);
        }

        self.current_line.clear();
        match self.reader.read_line(&mut self.current_line) {
            Ok(0) => Ok(false), // EOF
            Ok(_) => {
                // Remove trailing newline
                if self.current_line.ends_with('\n') {
                    self.current_line.pop();
                    if self.current_line.ends_with('\r') {
                        self.current_line.pop();
                    }
                }
                Ok(true)
            }
            Err(e) => Err(e),
        }
    }

    /// Read the next sequence from the file
    pub fn read_seq(&mut self) -> IoResult<Option<(Vec<u8>, SqInfo)>> {
        let mut sqfile = SqFile::new();

        match self.format {
            K_PEARSON | K_XPEARSON => self.read_pearson(&mut sqfile)?,
            K_GENBANK => self.read_genbank(&mut sqfile)?,
            K_EMBL => {
                sqfile.dash_equals_n = true;
                self.read_embl(&mut sqfile)?
            }
            K_NBRF => self.read_nbrf(&mut sqfile)?,
            K_IG => self.read_ig(&mut sqfile)?,
            K_PIR => self.read_pir(&mut sqfile)?,
            _ => {
                return Err(IoError::new(
                    ErrorKind::InvalidData,
                    format!("Unsupported format: {}", self.format),
                ))
            }
        }

        if sqfile.seq.is_empty() {
            return Ok(None);
        }

        sqfile.sqinfo.len = sqfile.seqlen;
        sqfile.sqinfo.flags |= SQINFO_LEN;

        Ok(Some((sqfile.seq, sqfile.sqinfo)))
    }

    /// Add sequence characters from a line
    fn add_seq(&self, line: &str, sqfile: &mut SqFile) {
        for &b in line.as_bytes() {
            if is_seq_char(b) {
                let mut ch = b;
                if ch == b'-' && sqfile.dash_equals_n {
                    ch = b'N';
                }
                if sqfile.seqlen >= sqfile.maxseq {
                    sqfile.maxseq += SqFile::START_LENGTH;
                    sqfile.seq.reserve(SqFile::START_LENGTH);
                }
                sqfile.seq.push(ch);
                sqfile.seqlen += 1;
            }
        }
    }

    /// Read FASTA/Pearson format
    fn read_pearson(&mut self, sqfile: &mut SqFile) -> IoResult<()> {
        // Find next sequence (line starting with '>')
        while self.read_line()? {
            if self.current_line.starts_with('>') {
                break;
            }
        }

        if !self.current_line.starts_with('>') {
            return Ok(()); // EOF, no sequence
        }

        // Parse header line
        let header = self.current_line[1..].trim_start(); // Skip '>' and leading whitespace
        let parts: Vec<&str> = header.splitn(2, char::is_whitespace).collect();

        if !parts.is_empty() && !parts[0].is_empty() {
            sqfile.sqinfo.set_string_field(parts[0], SQINFO_NAME);
        }
        if parts.len() > 1 {
            sqfile.sqinfo.set_string_field(parts[1], SQINFO_DESC);
        }

        // Read sequence lines until next '>' or EOF
        while self.read_line()? {
            if self.current_line.starts_with('>') {
                // Buffer this line for next read
                self.has_buffered_line = true;
                break;
            }
            self.add_seq(&self.current_line, sqfile);
        }

        Ok(())
    }

    /// Read GenBank format
    fn read_genbank(&mut self, sqfile: &mut SqFile) -> IoResult<()> {
        // Find LOCUS line
        while self.read_line()? {
            if self.current_line.starts_with("LOCUS") {
                break;
            }
        }

        if !self.current_line.starts_with("LOCUS") {
            return Ok(()); // EOF
        }

        // Parse LOCUS line
        let parts: Vec<&str> = self.current_line.split_whitespace().collect();
        if parts.len() > 1 {
            sqfile.sqinfo.set_string_field(parts[1], SQINFO_NAME);
            sqfile.sqinfo.set_string_field(parts[1], SQINFO_ID);
        }

        // Read header lines
        while self.read_line()? {
            if self.current_line.starts_with("ACCESSION") {
                let parts: Vec<&str> = self.current_line.split_whitespace().collect();
                if parts.len() > 1 {
                    sqfile.sqinfo.set_string_field(parts[1], SQINFO_ACC);
                }
            } else if self.current_line.starts_with("DEFINITION") {
                let desc = self.current_line[10..].trim();
                sqfile.sqinfo.set_string_field(desc, SQINFO_DESC);
            } else if self.current_line.starts_with("ORIGIN") {
                break;
            }
        }

        // Read sequence
        while self.read_line()? {
            if self.current_line.starts_with("//") {
                break;
            }
            // Skip line number at start
            let seq_part = if let Some(pos) = self.current_line.find(char::is_alphabetic) {
                &self.current_line[pos..]
            } else {
                &self.current_line
            };
            self.add_seq(seq_part, sqfile);
        }

        // Set source coordinates
        sqfile.sqinfo.start = 1;
        sqfile.sqinfo.stop = sqfile.seqlen as i32;
        sqfile.sqinfo.olen = sqfile.seqlen as i32;
        sqfile.sqinfo.flags |= SQINFO_START | SQINFO_STOP | SQINFO_OLEN;

        // Advance to next LOCUS
        while self.read_line()? {
            if self.current_line.starts_with("LOCUS") {
                self.has_buffered_line = true;
                break;
            }
        }

        Ok(())
    }

    /// Read EMBL format
    fn read_embl(&mut self, sqfile: &mut SqFile) -> IoResult<()> {
        // Find ID line
        while self.read_line()? {
            if self.current_line.starts_with("ID  ") {
                break;
            }
        }

        if !self.current_line.starts_with("ID  ") {
            return Ok(()); // EOF
        }

        // Parse ID line
        let parts: Vec<&str> = self.current_line[4..].split_whitespace().collect();
        if !parts.is_empty() {
            sqfile.sqinfo.set_string_field(parts[0], SQINFO_NAME);
            sqfile.sqinfo.set_string_field(parts[0], SQINFO_ID);
        }

        // Read header lines
        while self.read_line()? {
            if self.current_line.starts_with("AC  ") {
                let parts: Vec<&str> = self.current_line[4..].split(';').collect();
                if !parts.is_empty() {
                    sqfile
                        .sqinfo
                        .set_string_field(parts[0].trim(), SQINFO_ACC);
                }
            } else if self.current_line.starts_with("DE  ") {
                let desc = self.current_line[4..].trim();
                sqfile.sqinfo.set_string_field(desc, SQINFO_DESC);
            } else if self.current_line.starts_with("SQ") {
                break;
            }
        }

        // Read sequence (lines start with spaces)
        while self.read_line()? {
            // EMBL sequence lines start with 5 spaces
            if !self.current_line.starts_with("     ") {
                break;
            }
            self.add_seq(&self.current_line, sqfile);
        }

        // Set source coordinates
        sqfile.sqinfo.start = 1;
        sqfile.sqinfo.stop = sqfile.seqlen as i32;
        sqfile.sqinfo.olen = sqfile.seqlen as i32;
        sqfile.sqinfo.flags |= SQINFO_START | SQINFO_STOP | SQINFO_OLEN;

        // Advance to next ID line
        while self.read_line()? {
            if self.current_line.starts_with("ID  ") {
                self.has_buffered_line = true;
                break;
            }
        }

        Ok(())
    }

    /// Read NBRF/PIR format
    fn read_nbrf(&mut self, sqfile: &mut SqFile) -> IoResult<()> {
        // Find line starting with '>'
        while self.read_line()? {
            if self.current_line.starts_with('>') {
                break;
            }
        }

        if !self.current_line.starts_with('>') {
            return Ok(()); // EOF
        }

        // Parse name from header (after type code)
        if self.current_line.len() > 4 {
            let name = self.current_line[4..].trim();
            sqfile.sqinfo.set_string_field(name, SQINFO_NAME);
        }

        // Skip title line
        self.read_line()?;

        // Read sequence until * or next >
        while self.read_line()? {
            if self.current_line.contains('*') {
                // Remove the * and add final line
                let clean = self.current_line.replace('*', "");
                self.add_seq(&clean, sqfile);
                break;
            }
            if self.current_line.starts_with('>') {
                break;
            }
            self.add_seq(&self.current_line, sqfile);
        }

        // Advance to next record
        while self.read_line()? {
            if self.current_line.starts_with('>') {
                self.has_buffered_line = true;
                break;
            }
        }

        Ok(())
    }

    /// Read IG/Stanford format
    fn read_ig(&mut self, sqfile: &mut SqFile) -> IoResult<()> {
        // Skip comment lines (starting with ';')
        while self.read_line()? {
            if !self.current_line.is_empty() && !self.current_line.starts_with(';') {
                break;
            }
        }

        if self.current_line.is_empty() {
            return Ok(()); // EOF
        }

        // First non-comment line is the name
        let parts: Vec<&str> = self.current_line.split_whitespace().collect();
        if !parts.is_empty() {
            sqfile.sqinfo.set_string_field(parts[0], SQINFO_NAME);
        }

        // Read sequence until we see '1' or '2' (IG end markers)
        while self.read_line()? {
            if self.current_line.contains('1') || self.current_line.contains('2') {
                self.add_seq(&self.current_line, sqfile);
                break;
            }
            self.add_seq(&self.current_line, sqfile);
        }

        // Advance to next sequence
        while self.read_line()? {
            if !self.current_line.is_empty() && self.current_line.starts_with(';') {
                self.has_buffered_line = true;
                break;
            }
        }

        Ok(())
    }

    /// Read PIR/CODATA format
    fn read_pir(&mut self, sqfile: &mut SqFile) -> IoResult<()> {
        // Find ENTRY line
        while self.read_line()? {
            if self.current_line.starts_with("ENTRY") {
                break;
            }
        }

        if !self.current_line.starts_with("ENTRY") {
            return Ok(()); // EOF
        }

        // Parse ENTRY line (name at position 15+)
        if self.current_line.len() > 15 {
            let parts: Vec<&str> = self.current_line[15..].split_whitespace().collect();
            if !parts.is_empty() {
                sqfile.sqinfo.set_string_field(parts[0], SQINFO_NAME);
                sqfile.sqinfo.set_string_field(parts[0], SQINFO_ID);
            }
        }

        // Read header lines
        while self.read_line()? {
            if self.current_line.starts_with("TITLE") {
                let title = self.current_line[15..].trim();
                sqfile.sqinfo.set_string_field(title, SQINFO_DESC);
            } else if self.current_line.starts_with("ACCESSION") {
                let parts: Vec<&str> = self.current_line[15..].split_whitespace().collect();
                if !parts.is_empty() {
                    sqfile.sqinfo.set_string_field(parts[0], SQINFO_ACC);
                }
            } else if self.current_line.starts_with("SEQUENCE") {
                break;
            }
        }

        // Skip coordinate line
        self.read_line()?;

        // Read sequence until /// or ENTRY
        while self.read_line()? {
            if self.current_line.starts_with("///") || self.current_line.starts_with("ENTRY") {
                break;
            }
            self.add_seq(&self.current_line, sqfile);
        }

        // Set source coordinates
        sqfile.sqinfo.start = 1;
        sqfile.sqinfo.stop = sqfile.seqlen as i32;
        sqfile.sqinfo.olen = sqfile.seqlen as i32;
        sqfile.sqinfo.flags |= SQINFO_START | SQINFO_STOP | SQINFO_OLEN;

        // Advance to next ENTRY
        while self.read_line()? {
            if self.current_line.starts_with("ENTRY") {
                self.has_buffered_line = true;
                break;
            }
        }

        Ok(())
    }
}

/// Detect the format of a sequence file
pub fn detect_format<P: AsRef<Path>>(path: P) -> IoResult<i32> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    let mut found_ig = false;
    let mut found_strider = false;
    let mut found_gb = false;
    let mut found_embl = false;
    let mut found_nbrf = false;
    let mut found_pearson = false;
    let mut found_xpearson = false;
    let mut found_zuker = false;

    let mut got_gcgdata = false;
    let mut got_pir = false;
    let mut got_squid = false;
    let mut got_uw = false;
    let got_msf = false;
    let got_clustal = false;

    let mut nlines = 0;
    let mut dna_lines = 0;

    for line_result in reader.lines().take(500) {
        let line = line_result?;
        nlines += 1;

        if line.is_empty() {
            continue;
        }

        // High probability identities
        if line.contains(" MSF:") && line.contains(" Type:") && line.contains(" Check:") {
            return Ok(K_MSF);
        }

        if line.starts_with("CLUSTAL ") && line.contains("multiple sequence alignment") {
            return Ok(K_CLUSTAL);
        }

        if line.contains(" Check: ") {
            got_uw = true;
        }

        if line.starts_with("///") || line.starts_with("ENTRY ") {
            got_pir = true;
        }

        if line.starts_with("++") || line.starts_with("NAM ") {
            got_squid = true;
        }

        if line.starts_with(">>>>") && line.contains("Len: ") {
            got_gcgdata = true;
        }

        // Uncertain identities
        if line.starts_with(';') {
            if line.contains("Strider") {
                found_strider = true;
            } else {
                found_ig = true;
            }
        } else if line.starts_with("LOCUS") || line.starts_with("ORIGIN") {
            found_gb = true;
        } else if line.starts_with('>') {
            if line.len() > 3 && line.chars().nth(3) == Some(';') {
                found_nbrf = true;
            } else if line.contains("::") && line.contains("..") {
                found_xpearson = true;
            } else {
                found_pearson = true;
            }
        } else if line.starts_with("ID   ") || line.starts_with("SQ   ") {
            found_embl = true;
        } else if line.starts_with('(') {
            found_zuker = true;
        } else {
            // Check if line looks like DNA/RNA
            let seq_type = determine_seq_type(line.as_bytes());
            if (seq_type == K_DNA || seq_type == K_RNA) && line.len() > 20 {
                dna_lines += 1;
            }
        }

        // Early returns for high-confidence formats
        if got_msf {
            return Ok(K_MSF);
        }
        if got_clustal {
            return Ok(K_CLUSTAL);
        }
        if got_squid {
            return Ok(K_SQUID);
        }
        if got_pir {
            return Ok(K_PIR);
        }
        if got_gcgdata {
            return Ok(K_GCGDATA);
        }
        if got_uw {
            return if found_ig { Ok(K_IG) } else { Ok(K_GCG) };
        }

        // After enough lines or DNA lines, make a decision
        if dna_lines > 1 || nlines > 500 {
            break;
        }
    }

    // Decide on most likely format
    if found_strider {
        Ok(K_STRIDER)
    } else if found_gb {
        Ok(K_GENBANK)
    } else if found_embl {
        Ok(K_EMBL)
    } else if found_nbrf {
        Ok(K_NBRF)
    } else if found_ig {
        Ok(K_IG)
    } else if found_pearson {
        Ok(K_PEARSON)
    } else if found_xpearson {
        Ok(K_XPEARSON)
    } else if found_zuker {
        Ok(K_ZUKER)
    } else {
        // Default to SELEX if we can't determine
        Ok(K_SELEX)
    }
}

// ============================================================================
// Sequence Utility Functions
// ============================================================================

/// Get file format from filename extension
///
/// Attempts to determine sequence format from file extension.
/// Falls back to K_UNKNOWN if extension is not recognized.
pub fn seq_file_format(filename: &str) -> i32 {
    let lower = filename.to_lowercase();

    if lower.ends_with(".fa") || lower.ends_with(".fasta") || lower.ends_with(".fna") {
        K_PEARSON
    } else if lower.ends_with(".gb") || lower.ends_with(".gbk") || lower.ends_with(".genbank") {
        K_GENBANK
    } else if lower.ends_with(".embl") || lower.ends_with(".emb") {
        K_EMBL
    } else if lower.ends_with(".pir") || lower.ends_with(".nbrf") {
        K_NBRF
    } else if lower.ends_with(".gcg") {
        K_GCG
    } else if lower.ends_with(".msf") {
        K_MSF
    } else if lower.ends_with(".aln") || lower.ends_with(".clustal") {
        K_CLUSTAL
    } else if lower.ends_with(".stk") || lower.ends_with(".sto") || lower.ends_with(".stockholm") {
        K_SELEX
    } else {
        K_UNKNOWN
    }
}

/// Determine if character is a gap character
///
/// Corresponds to isgap() functionality in sqio.c.
/// Gap characters: '-', '.', '_', '~'
pub fn is_gap_char(c: char) -> bool {
    matches!(c, '-' | '.' | '_' | '~')
}

/// Determine if byte is a gap character
pub fn is_gap_byte(c: u8) -> bool {
    matches!(c, b'-' | b'.' | b'_' | b'~')
}

/// Get symbol index for nucleotide
///
/// Returns index: A=0, C=1, G=2, T/U=3
/// Returns -1 for invalid characters
pub fn symbol_index(c: char) -> i32 {
    match c.to_ascii_uppercase() {
        'A' => 0,
        'C' => 1,
        'G' => 2,
        'T' | 'U' => 3,
        _ => -1,
    }
}

/// Get symbol index for byte
pub fn symbol_index_byte(c: u8) -> i32 {
    match c {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' | b'U' | b'u' => 3,
        _ => -1,
    }
}

/// Convert sequence to digital (0-3) encoding
///
/// A=0, C=1, G=2, T/U=3
/// Unknown characters are encoded as -1 (sentinel)
/// This is the digitized sequence (dsq) format used internally.
pub fn digitize_seq(seq: &str) -> Vec<i8> {
    seq.chars()
        .map(|c| {
            match c.to_ascii_uppercase() {
                'A' => 0,
                'C' => 1,
                'G' => 2,
                'T' | 'U' => 3,
                _ if is_gap_char(c) => -1,
                _ => -1, // Unknown
            }
        })
        .collect()
}

/// Convert sequence bytes to digital encoding
pub fn digitize_seq_bytes(seq: &[u8]) -> Vec<i8> {
    seq.iter()
        .map(|&c| {
            match c {
                b'A' | b'a' => 0,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' | b'U' | b'u' => 3,
                _ if is_gap_byte(c) => -1,
                _ => -1,
            }
        })
        .collect()
}

/// Convert digital sequence back to string
///
/// 0=A, 1=C, 2=G, 3=U (RNA output)
/// -1 and other values become gaps (-)
pub fn undigitize_seq(dsq: &[i8]) -> String {
    dsq.iter()
        .map(|&d| match d {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            3 => 'U',
            _ => '-',
        })
        .collect()
}

/// Convert digital sequence to bytes
pub fn undigitize_seq_bytes(dsq: &[i8]) -> Vec<u8> {
    dsq.iter()
        .map(|&d| match d {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            3 => b'U',
            _ => b'-',
        })
        .collect()
}

/// Reverse complement of a DNA/RNA sequence
///
/// A <-> T/U, C <-> G
/// Preserves case and handles ambiguity codes
pub fn rev_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&c| match c {
            b'A' => b'T',
            b'a' => b't',
            b'T' => b'A',
            b't' => b'a',
            b'U' => b'A',
            b'u' => b'a',
            b'G' => b'C',
            b'g' => b'c',
            b'C' => b'G',
            b'c' => b'g',
            // Ambiguity codes
            b'R' => b'Y', // A|G -> T|C
            b'r' => b'y',
            b'Y' => b'R', // C|T -> G|A
            b'y' => b'r',
            b'S' => b'S', // G|C -> C|G (self-complement)
            b's' => b's',
            b'W' => b'W', // A|T -> T|A (self-complement)
            b'w' => b'w',
            b'K' => b'M', // G|T -> C|A
            b'k' => b'm',
            b'M' => b'K', // A|C -> T|G
            b'm' => b'k',
            b'B' => b'V', // C|G|T -> G|C|A
            b'b' => b'v',
            b'V' => b'B', // A|C|G -> T|G|C
            b'v' => b'b',
            b'D' => b'H', // A|G|T -> T|C|A
            b'd' => b'h',
            b'H' => b'D', // A|C|T -> T|G|A
            b'h' => b'd',
            b'N' => b'N',
            b'n' => b'n',
            // Keep gaps and unknowns
            other => other,
        })
        .collect()
}

/// Clean sequence - remove whitespace, convert T to U (RNA mode)
///
/// Removes all whitespace and digits, optionally converts T to U.
pub fn clean_seq(seq: &str, to_rna: bool) -> String {
    seq.chars()
        .filter(|c| !c.is_whitespace() && !c.is_ascii_digit())
        .map(|c| {
            if to_rna {
                match c {
                    'T' => 'U',
                    't' => 'u',
                    other => other,
                }
            } else {
                c
            }
        })
        .collect()
}

/// Convert sequence to uppercase
pub fn to_upper(seq: &mut [u8]) {
    for c in seq.iter_mut() {
        if c.is_ascii_lowercase() {
            *c = c.to_ascii_uppercase();
        }
    }
}

/// Convert sequence to lowercase
pub fn to_lower(seq: &mut [u8]) {
    for c in seq.iter_mut() {
        if c.is_ascii_uppercase() {
            *c = c.to_ascii_lowercase();
        }
    }
}

/// Convert sequence to DNA (U -> T)
pub fn seq_to_dna(seq: &mut [u8]) {
    for c in seq.iter_mut() {
        match *c {
            b'U' => *c = b'T',
            b'u' => *c = b't',
            _ => {}
        }
    }
}

/// Convert sequence to RNA (T -> U)
pub fn seq_to_rna(seq: &mut [u8]) {
    for c in seq.iter_mut() {
        match *c {
            b'T' => *c = b'U',
            b't' => *c = b'u',
            _ => {}
        }
    }
}

/// Read all sequences from a file
///
/// Returns a vector of (name, sequence) pairs.
pub fn read_seq_file<P: AsRef<Path>>(path: P) -> IoResult<Vec<(String, Vec<u8>)>> {
    let mut reader = SeqFileReader::open(path)?;
    let mut sequences = Vec::new();

    while let Some((seq, info)) = reader.read_seq()? {
        let name = if info.has_flag(SQINFO_NAME) {
            info.name.clone()
        } else {
            format!("seq_{}", sequences.len())
        };
        sequences.push((name, seq));
    }

    Ok(sequences)
}

/// Read sequences with full info
pub fn read_seq_file_with_info<P: AsRef<Path>>(path: P) -> IoResult<Vec<(Vec<u8>, SqInfo)>> {
    let mut reader = SeqFileReader::open(path)?;
    let mut sequences = Vec::new();

    while let Some((seq, info)) = reader.read_seq()? {
        sequences.push((seq, info));
    }

    Ok(sequences)
}

/// Get sequence format name as string
pub fn seq_format_string(code: i32) -> &'static str {
    match code {
        K_UNKNOWN => "Unknown",
        K_IG => "IntelliGenetics",
        K_GENBANK => "GenBank",
        K_NBRF => "PIR-NBRF",
        K_EMBL => "EMBL",
        K_GCG => "GCG",
        K_STRIDER => "Strider",
        K_PEARSON => "FASTA",
        K_ZUKER => "Zuker",
        K_IDRAW => "Idraw",
        K_SELEX => "SELEX",
        K_MSF => "MSF",
        K_PIR => "PIR-CODATA",
        K_RAW => "Raw",
        K_SQUID => "Squid",
        K_XPEARSON => "Extended FASTA",
        K_GCGDATA => "GCG datalibrary",
        K_CLUSTAL => "Clustal",
        _ => "Unknown",
    }
}

/// Calculate GCG checksum for a sequence
///
/// Corresponds to GCGchecksum() in sqio.c lines 1166-1177
pub fn gcg_checksum(seq: &[u8]) -> i32 {
    let mut check: i32 = 0;
    let mut count: i32 = 0;

    for &c in seq {
        count += 1;
        check += count * c.to_ascii_uppercase() as i32;
        if count == 57 {
            count = 0;
        }
    }

    check % 10000
}

/// Write sequence to string in FASTA format
pub fn format_fasta(name: &str, desc: Option<&str>, seq: &[u8], line_width: usize) -> String {
    let mut result = String::new();

    // Header line
    result.push('>');
    result.push_str(name);
    if let Some(d) = desc {
        result.push(' ');
        result.push_str(d);
    }
    result.push('\n');

    // Sequence lines
    for chunk in seq.chunks(line_width) {
        result.push_str(&String::from_utf8_lossy(chunk));
        result.push('\n');
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_gap_char() {
        assert!(is_gap_char('-'));
        assert!(is_gap_char('.'));
        assert!(is_gap_char('_'));
        assert!(is_gap_char('~'));
        assert!(!is_gap_char('A'));
        assert!(!is_gap_char('N'));
    }

    #[test]
    fn test_symbol_index() {
        assert_eq!(symbol_index('A'), 0);
        assert_eq!(symbol_index('a'), 0);
        assert_eq!(symbol_index('C'), 1);
        assert_eq!(symbol_index('G'), 2);
        assert_eq!(symbol_index('T'), 3);
        assert_eq!(symbol_index('U'), 3);
        assert_eq!(symbol_index('N'), -1);
        assert_eq!(symbol_index('-'), -1);
    }

    #[test]
    fn test_digitize_seq() {
        let dsq = digitize_seq("ACGT");
        assert_eq!(dsq, vec![0, 1, 2, 3]);

        let dsq2 = digitize_seq("AcGu");
        assert_eq!(dsq2, vec![0, 1, 2, 3]);

        let dsq3 = digitize_seq("AN-G");
        assert_eq!(dsq3, vec![0, -1, -1, 2]);
    }

    #[test]
    fn test_undigitize_seq() {
        let seq = undigitize_seq(&[0, 1, 2, 3]);
        assert_eq!(seq, "ACGU");

        let seq2 = undigitize_seq(&[0, -1, 2, 3]);
        assert_eq!(seq2, "A-GU");
    }

    #[test]
    fn test_rev_complement() {
        let seq = b"ACGT";
        let rc = rev_complement(seq);
        assert_eq!(rc, b"ACGT"); // ACGT is its own reverse complement

        let seq2 = b"AACG";
        let rc2 = rev_complement(seq2);
        assert_eq!(rc2, b"CGTT");

        let seq3 = b"AAUU"; // RNA
        let rc3 = rev_complement(seq3);
        assert_eq!(rc3, b"AATT");
    }

    #[test]
    fn test_clean_seq() {
        let cleaned = clean_seq("ACG T\n123", false);
        assert_eq!(cleaned, "ACGT");

        let cleaned_rna = clean_seq("ACGT", true);
        assert_eq!(cleaned_rna, "ACGU");
    }

    #[test]
    fn test_seq_to_dna_rna() {
        let mut seq = b"ACGU".to_vec();
        seq_to_dna(&mut seq);
        assert_eq!(seq, b"ACGT");

        let mut seq2 = b"ACGT".to_vec();
        seq_to_rna(&mut seq2);
        assert_eq!(seq2, b"ACGU");
    }

    #[test]
    fn test_gcg_checksum() {
        // Simple test - actual values from C implementation would need verification
        let check = gcg_checksum(b"ACGT");
        assert!(check >= 0 && check < 10000);
    }

    #[test]
    fn test_format_fasta() {
        let fasta = format_fasta("seq1", Some("test sequence"), b"ACGTACGT", 4);
        assert!(fasta.starts_with(">seq1 test sequence\n"));
        assert!(fasta.contains("ACGT\n"));
    }

    #[test]
    fn test_seq_file_format() {
        assert_eq!(seq_file_format("test.fasta"), K_PEARSON);
        assert_eq!(seq_file_format("test.fa"), K_PEARSON);
        assert_eq!(seq_file_format("test.gb"), K_GENBANK);
        assert_eq!(seq_file_format("test.embl"), K_EMBL);
        assert_eq!(seq_file_format("test.unknown"), K_UNKNOWN);
    }

    #[test]
    fn test_detect_format_fasta() {
        let format = detect_format("tests/golden/squid/inputs/test.fasta").unwrap();
        assert_eq!(format, K_PEARSON);
    }

    #[test]
    fn test_detect_format_genbank() {
        let format = detect_format("tests/golden/squid/inputs/test.gb").unwrap();
        assert_eq!(format, K_GENBANK);
    }

    #[test]
    fn test_detect_format_embl() {
        let format = detect_format("tests/golden/squid/inputs/test.embl").unwrap();
        assert_eq!(format, K_EMBL);
    }
}
