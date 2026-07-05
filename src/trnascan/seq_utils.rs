//! Sequence utility functions for tRNAscan.
//!
//! This module provides utilities for sequence manipulation including
//! reverse complement, anticodon encoding, and amino acid lookup.

/// Compute the reverse complement of a DNA sequence.
///
/// This implements the compstrand function from trnascan.c (lines 1726-1773).
/// Complementation rules: A<->T, C<->G, other->N
///
/// # Arguments
/// * `seq` - Input DNA sequence as bytes (lowercase a, c, g, t)
///
/// # Returns
/// The reverse complement sequence
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'a' | b'A' => b't',
            b't' | b'T' => b'a',
            b'c' | b'C' => b'g',
            b'g' | b'G' => b'c',
            _ => b'n',
        })
        .collect()
}

/// Encode an anticodon triplet as a number 1-65.
///
/// This implements the codage function from trnascan.c (lines 1782-1821).
/// Encoding is alphabetical: AAA=1, AAC=2, ..., TTT=64.
/// If any base is ambiguous (not A, C, G, T), returns 65.
///
/// # Arguments
/// * `anticodon` - 3-byte anticodon sequence (lowercase)
///
/// # Returns
/// Encoded number 1-64, or 65 for indeterminate
pub fn encode_anticodon(anticodon: &[u8]) -> i32 {
    if anticodon.len() < 3 {
        return 65;
    }

    let mut num = 1;
    let mut iba = 1;
    let mut has_ambiguous = false;

    // Process in reverse order (positions 2, 1, 0)
    for i in (0..3).rev() {
        let j = match anticodon[i] {
            b'a' | b'A' => 0,
            b'c' | b'C' => 1,
            b'g' | b'G' => 2,
            b't' | b'T' => 3,
            _ => {
                has_ambiguous = true;
                0
            }
        };

        if has_ambiguous {
            break;
        }

        num += j * iba;
        iba *= 4;
    }

    if has_ambiguous {
        65
    } else {
        num
    }
}

/// Amino acid lookup table indexed by anticodon code (1-65).
///
/// This table gives the correspondence between the anticodon code
/// and the amino acid. Index 0-63 correspond to codes 1-64,
/// index 64 corresponds to code 65 (indeterminate).
const AMINO_ACID_TABLE: [&str; 65] = [
    "Phe", "Val", "Leu", "Ile", // AAA, AAC, AAG, AAT -> codes 1-4
    "Cys", "Trp", "Arg", "Ser", // ACA, ACC, ACG, ACT -> codes 5-8
    "Ser", "Ala", "Pro", "Thr", // AGA, AGC, AGG, AGT -> codes 9-12
    "Tyr", "Asp", "His", "Asn", // ATA, ATC, ATG, ATT -> codes 13-16
    "Leu", "Val", "Leu", "Met", // CAA, CAC, CAG, CAT -> codes 17-20
    "Trp", "Gly", "Arg", "Arg", // CCA, CCC, CCG, CCT -> codes 21-24
    "Ser", "Ala", "Pro", "Thr", // CGA, CGC, CGG, CGT -> codes 25-28
    "Sup", "Glu", "Gln", "Lys", // CTA, CTC, CTG, CTT -> codes 29-32
    "Phe", "Val", "Leu", "Ile", // GAA, GAC, GAG, GAT -> codes 33-36
    "Cys", "Gly", "Arg", "Ser", // GCA, GCC, GCG, GCT -> codes 37-40
    "Ser", "Ala", "Pro", "Thr", // GGA, GGC, GGG, GGT -> codes 41-44
    "Tyr", "Asp", "His", "Asn", // GTA, GTC, GTG, GTT -> codes 45-48
    "Leu", "Val", "Leu", "Ile", // TAA, TAC, TAG, TAT -> codes 49-52
    "Sup", "Gly", "Arg", "Arg", // TCA, TCC, TCG, TCT -> codes 53-56
    "Ser", "Ala", "Pro", "Thr", // TGA, TGC, TGG, TGT -> codes 57-60
    "Sup", "Glu", "Gln", "Lys", // TTA, TTC, TTG, TTT -> codes 61-64
    "Ind",                      // Code 65 (indeterminate)
];

/// Get amino acid three-letter code from anticodon encoding.
///
/// This implements the corresaa function from trnascan.c (lines 1826-1850).
///
/// # Arguments
/// * `num` - Anticodon code (1-65) from encode_anticodon
///
/// # Returns
/// Three-letter amino acid code, "Sup" for stop codons, "Ind" for indeterminate
pub fn anticodon_to_aa(num: i32) -> &'static str {
    if num < 1 || num > 65 {
        return "Ind";
    }
    AMINO_ACID_TABLE[(num - 1) as usize]
}

/// Convert a base character to index (A=0, C=1, G=2, T=3, other=-1).
#[inline]
pub fn base_to_index(b: u8) -> i32 {
    match b {
        b'a' | b'A' => 0,
        b'c' | b'C' => 1,
        b'g' | b'G' => 2,
        b't' | b'T' => 3,
        _ => -1,
    }
}

/// Check if a byte is a valid nucleotide (A, C, G, T).
#[inline]
pub fn is_nucleotide(b: u8) -> bool {
    matches!(b, b'a' | b'A' | b'c' | b'C' | b'g' | b'G' | b't' | b'T')
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement_simple() {
        let seq = b"acgt";
        let rc = reverse_complement(seq);
        assert_eq!(rc, b"acgt"); // acgt reversed complement is acgt
    }

    #[test]
    fn test_reverse_complement_longer() {
        let seq = b"aaaaaa";
        let rc = reverse_complement(seq);
        assert_eq!(rc, b"tttttt");
    }

    #[test]
    fn test_reverse_complement_mixed() {
        let seq = b"atcgatcg";
        let rc = reverse_complement(seq);
        assert_eq!(rc, b"cgatcgat");
    }

    #[test]
    fn test_reverse_complement_ambiguous() {
        let seq = b"acngt";
        let rc = reverse_complement(seq);
        assert_eq!(rc, b"acngt");
    }

    #[test]
    fn test_encode_anticodon_aaa() {
        assert_eq!(encode_anticodon(b"aaa"), 1);
    }

    #[test]
    fn test_encode_anticodon_aac() {
        assert_eq!(encode_anticodon(b"aac"), 2);
    }

    #[test]
    fn test_encode_anticodon_ttt() {
        // TTT should be 64 (all T's)
        // T=3, so 3*16 + 3*4 + 3*1 + 1 = 48 + 12 + 3 + 1 = 64
        assert_eq!(encode_anticodon(b"ttt"), 64);
    }

    #[test]
    fn test_encode_anticodon_ambiguous() {
        assert_eq!(encode_anticodon(b"ann"), 65);
        assert_eq!(encode_anticodon(b"acn"), 65);
    }

    #[test]
    fn test_anticodon_to_aa() {
        assert_eq!(anticodon_to_aa(1), "Phe");  // AAA
        assert_eq!(anticodon_to_aa(20), "Met"); // CAT
        assert_eq!(anticodon_to_aa(65), "Ind"); // Indeterminate
    }

    #[test]
    fn test_anticodon_to_aa_stop_codons() {
        assert_eq!(anticodon_to_aa(29), "Sup"); // CTA
        assert_eq!(anticodon_to_aa(53), "Sup"); // TCA
        assert_eq!(anticodon_to_aa(61), "Sup"); // TTA
    }

    #[test]
    fn test_base_to_index() {
        assert_eq!(base_to_index(b'a'), 0);
        assert_eq!(base_to_index(b'A'), 0);
        assert_eq!(base_to_index(b'c'), 1);
        assert_eq!(base_to_index(b'g'), 2);
        assert_eq!(base_to_index(b't'), 3);
        assert_eq!(base_to_index(b'n'), -1);
    }

    #[test]
    fn test_is_nucleotide() {
        assert!(is_nucleotide(b'a'));
        assert!(is_nucleotide(b'A'));
        assert!(is_nucleotide(b'c'));
        assert!(is_nucleotide(b'g'));
        assert!(is_nucleotide(b't'));
        assert!(!is_nucleotide(b'n'));
        assert!(!is_nucleotide(b'N'));
    }
}
