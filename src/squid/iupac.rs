// IUPAC nucleotide and amino acid codes
// Ported from original SQUID library

/// Valid amino acid characters
pub const AMINOS: &[u8] = b"ABCDEFGHIKLMNPQRSTVWXYZ*";

/// Primary nucleotide characters
pub const PRIME_NUC: &[u8] = b"ACGTUN";

/// Protein-only characters (not in nucleotide sequences)
pub const PROT_ONLY: &[u8] = b"EFIPQZ";

/// All allowed symbols in sequences
pub const ALL_SYMBOLS: &[u8] = b"_.-*?<>{}[]()!@#$%^&=+;:'|`~\"\\";

/// Check if character is a valid sequence character
#[inline]
pub fn is_seq_char(c: u8) -> bool {
    c.is_ascii_alphabetic() || ALL_SYMBOLS.contains(&c)
}

/// Check if character is a gap character
#[inline]
pub fn is_gap(c: u8) -> bool {
    c == b'-' || c == b'.' || c == b'_'
}

/// Determine sequence type from content
pub fn seq_type(seq: &[u8]) -> i32 {
    use crate::squid::types::{K_AMINO, K_DNA, K_OTHER_SEQ, K_RNA};

    if seq.is_empty() {
        return K_OTHER_SEQ;
    }

    let mut saw_t = false;
    let mut saw_u = false;
    let mut saw_prot = false;
    let mut total = 0;
    let mut nuclike = 0;

    for &c in seq {
        let upper = c.to_ascii_uppercase();
        if !c.is_ascii_alphabetic() {
            continue;
        }

        total += 1;

        // Check for protein-specific residues
        if PROT_ONLY.contains(&upper) {
            saw_prot = true;
        }

        // Check for T vs U (DNA vs RNA)
        if upper == b'T' {
            saw_t = true;
        }
        if upper == b'U' {
            saw_u = true;
        }

        // Count nucleotide-like characters
        if PRIME_NUC.contains(&upper) {
            nuclike += 1;
        }
    }

    // If we saw protein-only chars, it's protein
    if saw_prot {
        return K_AMINO;
    }

    // If mostly nucleotide-like characters
    if total > 0 && nuclike * 2 > total {
        // Decide between DNA and RNA
        if saw_u && !saw_t {
            return K_RNA;
        } else if saw_t && !saw_u {
            return K_DNA;
        } else if saw_t || saw_u {
            // Mixed or default to DNA
            return K_DNA;
        } else {
            // No T or U, assume DNA
            return K_DNA;
        }
    }

    // Check if it could be amino acid
    if total > 0 {
        let mut aa_count = 0;
        for &c in seq {
            let upper = c.to_ascii_uppercase();
            if AMINOS.contains(&upper) {
                aa_count += 1;
            }
        }
        if aa_count * 2 > total {
            return K_AMINO;
        }
    }

    K_OTHER_SEQ
}

/// Convert DNA sequence to RNA (T -> U)
pub fn to_rna(seq: &mut [u8]) {
    for c in seq.iter_mut() {
        if *c == b'T' {
            *c = b'U';
        } else if *c == b't' {
            *c = b'u';
        }
    }
}

/// Convert RNA sequence to DNA (U -> T)
pub fn to_dna(seq: &mut [u8]) {
    for c in seq.iter_mut() {
        if *c == b'U' {
            *c = b'T';
        } else if *c == b'u' {
            *c = b't';
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::squid::types::{K_AMINO, K_DNA, K_OTHER_SEQ, K_RNA};

    #[test]
    fn test_is_seq_char() {
        assert!(is_seq_char(b'A'));
        assert!(is_seq_char(b'z'));
        assert!(is_seq_char(b'-'));
        assert!(is_seq_char(b'*'));
        assert!(!is_seq_char(b' '));
        assert!(!is_seq_char(b'\n'));
    }

    #[test]
    fn test_is_gap() {
        assert!(is_gap(b'-'));
        assert!(is_gap(b'.'));
        assert!(is_gap(b'_'));
        assert!(!is_gap(b'A'));
        assert!(!is_gap(b' '));
    }

    #[test]
    fn test_seq_type_dna() {
        let dna = b"ACGTACGTACGT";
        assert_eq!(seq_type(dna), K_DNA);
    }

    #[test]
    fn test_seq_type_rna() {
        let rna = b"ACGUACGUACGU";
        assert_eq!(seq_type(rna), K_RNA);
    }

    #[test]
    fn test_seq_type_protein() {
        let protein = b"MEILPQRSTV";
        assert_eq!(seq_type(protein), K_AMINO);
    }

    #[test]
    fn test_seq_type_empty() {
        assert_eq!(seq_type(b""), K_OTHER_SEQ);
    }

    #[test]
    fn test_to_rna() {
        let mut seq = b"ACGTACGT".to_vec();
        to_rna(&mut seq);
        assert_eq!(seq, b"ACGUACGU");
    }

    #[test]
    fn test_to_dna() {
        let mut seq = b"ACGUACGU".to_vec();
        to_dna(&mut seq);
        assert_eq!(seq, b"ACGTACGT");
    }
}
