/// Utility functions for tRNAscan-SE
///
/// This module provides general utility functions including:
/// - File handling utilities
/// - String manipulation
/// - Sequence utilities (reverse complement, etc.)
/// - tRNA-specific helpers

use std::collections::HashMap;
use std::fs;
use std::io::{self, Write};
use std::path::{Path, PathBuf};

lazy_static::lazy_static! {
    /// Complement map for nucleotides
    static ref COMP_MAP: HashMap<u8, u8> = {
        let mut m = HashMap::new();
        m.insert(b'A', b'T'); m.insert(b'a', b't');
        m.insert(b'T', b'A'); m.insert(b't', b'a');
        m.insert(b'U', b'A'); m.insert(b'u', b'a');
        m.insert(b'G', b'C'); m.insert(b'g', b'c');
        m.insert(b'C', b'G'); m.insert(b'c', b'g');
        m.insert(b'Y', b'R'); m.insert(b'y', b'r');
        m.insert(b'R', b'Y'); m.insert(b'r', b'y');
        m.insert(b'S', b'S'); m.insert(b's', b's');
        m.insert(b'W', b'W'); m.insert(b'w', b'w');
        m.insert(b'M', b'K'); m.insert(b'm', b'k');
        m.insert(b'K', b'M'); m.insert(b'k', b'm');
        m.insert(b'B', b'V'); m.insert(b'b', b'v');
        m.insert(b'V', b'B'); m.insert(b'v', b'b');
        m.insert(b'H', b'D'); m.insert(b'h', b'd');
        m.insert(b'D', b'H'); m.insert(b'd', b'h');
        m.insert(b'N', b'N'); m.insert(b'n', b'n');
        m.insert(b'X', b'X'); m.insert(b'x', b'x');
        m.insert(b'?', b'?');
        m.insert(b'-', b'-');
        m
    };
}

// ============================================================================
// File utilities
// ============================================================================

/// Create a temporary file with given prefix and suffix
pub fn create_temp_file(prefix: &str, suffix: &str) -> io::Result<PathBuf> {
    let temp_dir = std::env::temp_dir();
    let pid = std::process::id();
    let filename = format!("{}{}{}", prefix, pid, suffix);
    Ok(temp_dir.join(filename))
}

/// Ensure a directory exists, creating it if necessary
pub fn ensure_dir(path: &Path) -> io::Result<()> {
    if !path.exists() {
        fs::create_dir_all(path)?;
    }
    Ok(())
}

/// Check if a file exists
pub fn file_exists(path: &Path) -> bool {
    path.exists() && path.is_file()
}

/// Copy a file from source to destination
pub fn copy_file(src: &Path, dst: &Path) -> io::Result<u64> {
    fs::copy(src, dst)
}

/// Print filename, converting "-" to "Standard output"
pub fn print_filename(fname: &str) -> String {
    if fname == "-" {
        "Standard output".to_string()
    } else {
        fname.to_string()
    }
}

// ============================================================================
// String utilities
// ============================================================================

/// Trim whitespace from both ends of a string
pub fn trim_whitespace(s: &str) -> String {
    s.trim().to_string()
}

/// Alignment for padding
pub enum Alignment {
    Left,
    Right,
    Center,
}

/// Pad string to given width with alignment
pub fn pad_string(s: &str, width: usize, align: Alignment) -> String {
    let len = s.len();
    if len >= width {
        return s.to_string();
    }

    let padding = width - len;
    match align {
        Alignment::Left => format!("{}{}", s, " ".repeat(padding)),
        Alignment::Right => format!("{}{}", " ".repeat(padding), s),
        Alignment::Center => {
            let left_pad = padding / 2;
            let right_pad = padding - left_pad;
            format!("{}{}{}", " ".repeat(left_pad), s, " ".repeat(right_pad))
        }
    }
}

/// Pad number with leading zeros
pub fn pad_num(num: usize, width: usize) -> String {
    format!("{:0width$}", num, width = width)
}

/// Truncate string to maximum length
pub fn truncate_string(s: &str, max_len: usize) -> String {
    if s.len() <= max_len {
        s.to_string()
    } else {
        s[..max_len].to_string()
    }
}

// ============================================================================
// Sequence utilities
// ============================================================================

/// Reverse complement a DNA/RNA sequence
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| *COMP_MAP.get(&b).unwrap_or(&b))
        .collect()
}

/// Complement a DNA/RNA sequence (without reversing)
pub fn complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|&b| *COMP_MAP.get(&b).unwrap_or(&b))
        .collect()
}

/// Get minimum of two values
pub fn min<T: Ord>(a: T, b: T) -> T {
    if a < b { a } else { b }
}

/// Get maximum of two values
pub fn max<T: Ord>(a: T, b: T) -> T {
    if a > b { a } else { b }
}

/// Check if two segments overlap
pub fn seg_overlap(seg1_a: i64, seg1_b: i64, seg2_a: i64, seg2_b: i64, range: i64) -> bool {
    if range == 0 {
        (seg1_a >= seg2_a && seg1_a <= seg2_b)
            || (seg1_b >= seg2_a && seg1_b <= seg2_b)
            || (seg2_a >= seg1_a && seg2_a <= seg1_b)
            || (seg2_b >= seg1_a && seg2_b <= seg1_b)
    } else {
        (seg1_a >= seg2_a - range && seg1_a <= seg2_a + range)
            || (seg1_b >= seg2_b - range && seg1_b <= seg2_b + range)
            || (seg2_a >= seg1_a - range && seg2_a <= seg1_a + range)
            || (seg2_b >= seg1_b - range && seg2_b <= seg1_b + range)
    }
}

// ============================================================================
// tRNA-specific utilities
// ============================================================================

/// Format anticodon for display (uppercase)
pub fn format_anticodon(ac: &str) -> String {
    ac.to_uppercase()
}

/// Validate anticodon (should be 3 nucleotides)
pub fn validate_anticodon(ac: &str) -> bool {
    if ac.len() != 3 {
        return false;
    }
    ac.bytes().all(|b| matches!(b, b'A'..=b'Z' | b'a'..=b'z'))
}

/// Convert amino acid single-letter code to isotype name
pub fn aa_code_to_isotype(code: &str) -> String {
    match code.to_uppercase().as_str() {
        "A" => "Ala",
        "C" => "Cys",
        "D" => "Asp",
        "E" => "Glu",
        "F" => "Phe",
        "G" => "Gly",
        "H" => "His",
        "I" => "Ile",
        "K" => "Lys",
        "L" => "Leu",
        "M" => "Met",
        "N" => "Asn",
        "P" => "Pro",
        "Q" => "Gln",
        "R" => "Arg",
        "S" => "Ser",
        "T" => "Thr",
        "V" => "Val",
        "W" => "Trp",
        "Y" => "Tyr",
        "*" => "Ter",
        "U" => "Sec",
        "O" => "Pyl",
        _ => "Unk",
    }
    .to_string()
}

/// Convert isotype name to amino acid single-letter code
pub fn isotype_to_aa_code(isotype: &str) -> String {
    match isotype {
        "Ala" => "A",
        "Cys" => "C",
        "Asp" => "D",
        "Glu" => "E",
        "Phe" => "F",
        "Gly" => "G",
        "His" => "H",
        "Ile" => "I",
        "Lys" => "K",
        "Leu" => "L",
        "Met" => "M",
        "Asn" => "N",
        "Pro" => "P",
        "Gln" => "Q",
        "Arg" => "R",
        "Ser" => "S",
        "Thr" => "T",
        "Val" => "V",
        "Trp" => "W",
        "Tyr" => "Y",
        "Ter" | "Sup" => "*",
        "Sec" => "U",
        "Pyl" => "O",
        _ => "?",
    }
    .to_string()
}

/// Check if a process exit status indicates an error
pub fn check_exit_status(prog_name: &str, seq_name: &str, status: i32) -> Result<(), String> {
    if status != 0 {
        Err(format!(
            "{} could not complete successfully for {}. \
             Possible memory allocation problem or missing file. (Exit code={})",
            prog_name, seq_name, status
        ))
    } else {
        Ok(())
    }
}

/// Write to file or stdout
pub fn write_output(path: &str, content: &str) -> io::Result<()> {
    if path == "-" {
        io::stdout().write_all(content.as_bytes())?;
    } else {
        fs::write(path, content)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        let seq = b"ACGT";
        let rc = reverse_complement(seq);
        assert_eq!(&rc, b"ACGT");

        let seq2 = b"AAAA";
        let rc2 = reverse_complement(seq2);
        assert_eq!(&rc2, b"TTTT");
    }

    #[test]
    fn test_complement() {
        let seq = b"ACGT";
        let c = complement(seq);
        assert_eq!(&c, b"TGCA");
    }

    #[test]
    fn test_pad_string() {
        assert_eq!(pad_string("test", 10, Alignment::Right), "      test");
        assert_eq!(pad_string("test", 10, Alignment::Left), "test      ");
    }

    #[test]
    fn test_pad_num() {
        assert_eq!(pad_num(42, 5), "00042");
        assert_eq!(pad_num(1, 3), "001");
    }

    #[test]
    fn test_seg_overlap() {
        assert!(seg_overlap(10, 20, 15, 25, 0));
        assert!(seg_overlap(10, 20, 5, 15, 0));
        assert!(!seg_overlap(10, 20, 25, 30, 0));
    }

    #[test]
    fn test_validate_anticodon() {
        assert!(validate_anticodon("AAA"));
        assert!(validate_anticodon("ttt"));
        assert!(!validate_anticodon("AA"));
        assert!(!validate_anticodon("AAAA"));
    }

    #[test]
    fn test_aa_code_conversion() {
        assert_eq!(aa_code_to_isotype("A"), "Ala");
        assert_eq!(aa_code_to_isotype("F"), "Phe");
        assert_eq!(isotype_to_aa_code("Ala"), "A");
        assert_eq!(isotype_to_aa_code("Phe"), "F");
    }
}
