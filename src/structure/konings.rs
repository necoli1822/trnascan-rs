//! Konings/Hogeweg secondary structure notation
//!
//! Implementation of secondary structure functions from konings.c
//! Reference: Konings and Hogeweg, J. Mol. Biol. 207:597-614 1989

use crate::types::constants::*;
use crate::types::trace::Trace;

/// Check if two RNA bases are complementary
///
/// Watson-Crick pairs: A-U, U-A, G-C, C-G
/// Wobble pairs (if allow_gu): G-U, U-G
///
/// Based on IsRNAComplement() from konings.c lines 224-246
///
/// # Arguments
/// * `base1` - First base (case insensitive, T treated as U)
/// * `base2` - Second base (case insensitive, T treated as U)
/// * `allow_gu` - Whether to allow G-U wobble pairs
///
/// # Returns
/// `true` if bases can pair, `false` otherwise
pub fn is_rna_complement(base1: char, base2: char, allow_gu: bool) -> bool {
    // Normalize to uppercase and convert T to U
    let mut b1 = base1.to_ascii_uppercase();
    let mut b2 = base2.to_ascii_uppercase();

    if b1 == 'T' {
        b1 = 'U';
    }
    if b2 == 'T' {
        b2 = 'U';
    }

    // Watson-Crick pairs
    match (b1, b2) {
        ('A', 'U') | ('U', 'A') | ('G', 'C') | ('C', 'G') => true,
        ('G', 'U') | ('U', 'G') if allow_gu => true,
        _ => false,
    }
}

/// Convert traceback tree to Konings/Hogeweg secondary structure string
///
/// Based on Trace2KHS() from konings.c lines 172-222
///
/// Walks the trace tree and marks paired positions with '>' (left) and '<' (right).
/// Unpaired positions are marked with '.'.
///
/// # Arguments
/// * `trace` - Traceback tree from alignment
/// * `seq` - Sequence being aligned (0-indexed)
/// * `rlen` - Length of sequence
/// * `watsoncrick` - If true, only annotate Watson-Crick canonical pairs
///
/// # Returns
/// Secondary structure string of length rlen
pub fn trace2khs(trace: &Trace, seq: &[u8], rlen: usize, watsoncrick: bool) -> String {
    // Initialize structure string with dots (unpaired)
    let mut ss = vec![b'.'; rlen];

    // Use a stack to traverse the trace tree (simulating the C code's tracestack)
    let mut stack: Vec<&Trace> = Vec::new();

    // Start with the left child of root (skip root BEGIN state)
    if let Some(ref left) = trace.nxtl {
        stack.push(left);
    }

    // Process trace tree using stack
    while let Some(curr) = stack.pop() {
        // Check if this is a MATP (paired) state
        if curr.trace_type == U_MATP_ST {
            // Convert emitl and emitr from 1-based to 0-based indices
            let left_idx = (curr.emitl - 1) as usize;
            let right_idx = (curr.emitr - 1) as usize;

            // Check bounds
            if left_idx < rlen && right_idx < rlen {
                // If watsoncrick is true, verify bases can pair
                let can_pair = if watsoncrick {
                    is_rna_complement(seq[left_idx] as char, seq[right_idx] as char, true)
                } else {
                    true // Accept all MATP pairs if not requiring Watson-Crick
                };

                if can_pair {
                    ss[left_idx] = b'>';
                    ss[right_idx] = b'<';
                }
            }
        }

        // Push children onto stack (right first, then left, for correct traversal order)
        if let Some(ref right) = curr.nxtr {
            stack.push(right);
        }
        if let Some(ref left) = curr.nxtl {
            stack.push(left);
        }
    }

    // Convert to string
    String::from_utf8(ss).expect("Invalid UTF-8 in structure string")
}

/// Convert Konings/Hogeweg structure string to connect table
///
/// Based on KHS2ct() from konings.c lines 249-322
///
/// The connect table ct[i] = j means position i pairs with position j (0-indexed).
/// Unpaired positions have ct[i] = -1.
///
/// Structure notation:
/// - '>' and '<' for base pairs
/// - 'A'-'Z' and 'a'-'z' for pseudoknot pairs (if allow_pseudoknots)
/// - Other characters are ignored (unpaired)
///
/// # Arguments
/// * `ss` - Secondary structure string
/// * `allow_pseudoknots` - Whether to parse pseudoknot notation
///
/// # Returns
/// `Ok(ct)` where ct[i] is the pairing partner of position i (or -1 if unpaired)
/// `Err(msg)` if structure string is inconsistent
pub fn khs2ct(ss: &str, allow_pseudoknots: bool) -> Result<Vec<i32>, String> {
    let len = ss.len();
    let mut ct = vec![-1i32; len];

    // Stack for each pairing level (0 = main structure, 1-26 = pseudoknots A-Z)
    let mut stacks: Vec<Vec<usize>> = vec![Vec::new(); 27];

    let chars: Vec<char> = ss.chars().collect();

    for (pos, &ch) in chars.iter().enumerate() {
        // Skip high-value characters (bulletproof against buggy ctype.h)
        if ch as u32 > 127 {
            return Err("Invalid character in structure string".to_string());
        }

        if ch == '>' {
            // Left side of pair: push onto stack 0
            stacks[0].push(pos);
        } else if ch == '<' {
            // Right side of pair: pop from stack 0 and create pairing
            if let Some(pair) = stacks[0].pop() {
                ct[pos] = pair as i32;
                ct[pair] = pos as i32;
            } else {
                return Err("Unmatched '<' in structure string".to_string());
            }
        } else if allow_pseudoknots && ch.is_ascii_uppercase() {
            // Pseudoknot left side: push onto corresponding stack
            let stack_idx = (ch as usize) - ('A' as usize) + 1;
            if stack_idx < 27 {
                stacks[stack_idx].push(pos);
            }
        } else if allow_pseudoknots && ch.is_ascii_lowercase() {
            // Pseudoknot right side: pop from corresponding stack
            let stack_idx = (ch as usize) - ('a' as usize) + 1;
            if stack_idx < 27 {
                if let Some(pair) = stacks[stack_idx].pop() {
                    ct[pos] = pair as i32;
                    ct[pair] = pos as i32;
                } else {
                    return Err(format!("Unmatched '{}' in structure string", ch));
                }
            }
        } else if allow_pseudoknots && !is_gap(ch) && ch != '.' && ch != ' ' {
            // Bad character in pseudoknot mode
            return Err(format!("Unexpected character '{}' in structure string", ch));
        }
        // Otherwise ignore the character (unpaired position)
    }

    // Check that all stacks are empty (all pairs matched)
    for (i, stack) in stacks.iter().enumerate() {
        if !stack.is_empty() {
            if i == 0 {
                return Err("Unmatched '>' in structure string".to_string());
            } else {
                let ch = ('A' as u8 + (i - 1) as u8) as char;
                return Err(format!("Unmatched '{}' in structure string", ch));
            }
        }
    }

    Ok(ct)
}

/// Check if character is a gap character
fn is_gap(ch: char) -> bool {
    ch == '-' || ch == '_' || ch == '~' || ch == ' '
}

/// Alignment state types (from structs.h)
#[derive(Debug, Clone, Copy, PartialEq)]
#[allow(dead_code)]
pub enum AlignType {
    Begin,
    Bifurc,
    Del,
    InsL,
    InsR,
    MatL,
    MatR,
    MatP,
}

/// Alignment structure element
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct AlignNode {
    pub align_type: AlignType,
    pub nodeidx: i32,
    pub pos: i32,        // Position in sequence (-1 if deleted)
    pub sym: char,       // Symbol at this position
    pub ss: char,        // Secondary structure annotation
}

/// Convert alignment to Konings/Hogeweg structure string
///
/// Based on Align2kh() from konings.c lines 25-103
///
/// Converts an alignment list into a secondary structure string.
/// Symbols used: '>' and '<' for base pairs (MATP), '.' for other positions.
///
/// Also returns the "aligned" sequence where deleted consensus positions
/// are '-' and insert positions are lowercase.
///
/// # Arguments
/// * `ali` - Alignment list
///
/// # Returns
/// `(aligned_sequence, structure_string)`
pub fn align2kh(ali: &[AlignNode]) -> (String, String) {
    let len = ali.len();
    let mut aseq = Vec::with_capacity(len);
    let mut khseq = Vec::with_capacity(len);

    for node in ali {
        match node.align_type {
            AlignType::Begin | AlignType::Bifurc => {
                // These shouldn't appear in alignment list
                continue;
            }
            AlignType::Del => {
                khseq.push(' ');
                aseq.push('-');
            }
            AlignType::InsL | AlignType::InsR => {
                khseq.push(node.ss);
                aseq.push(node.sym.to_ascii_lowercase());
            }
            AlignType::MatL | AlignType::MatR | AlignType::MatP => {
                khseq.push(node.ss);
                aseq.push(node.sym.to_ascii_uppercase());
            }
        }
    }

    (
        aseq.into_iter().collect(),
        khseq.into_iter().collect(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_rna_complement_watson_crick() {
        // Watson-Crick pairs
        assert!(is_rna_complement('A', 'U', false));
        assert!(is_rna_complement('U', 'A', false));
        assert!(is_rna_complement('G', 'C', false));
        assert!(is_rna_complement('C', 'G', false));

        // Case insensitive
        assert!(is_rna_complement('a', 'u', false));
        assert!(is_rna_complement('g', 'c', false));

        // T treated as U
        assert!(is_rna_complement('A', 'T', false));
        assert!(is_rna_complement('T', 'A', false));

        // Non-pairs
        assert!(!is_rna_complement('A', 'A', false));
        assert!(!is_rna_complement('A', 'G', false));
        assert!(!is_rna_complement('C', 'U', false));
    }

    #[test]
    fn test_is_rna_complement_wobble() {
        // G-U wobble pairs (only when allowed)
        assert!(!is_rna_complement('G', 'U', false));
        assert!(!is_rna_complement('U', 'G', false));

        assert!(is_rna_complement('G', 'U', true));
        assert!(is_rna_complement('U', 'G', true));

        // Case insensitive
        assert!(is_rna_complement('g', 'u', true));

        // G-T wobble (T as U)
        assert!(is_rna_complement('G', 'T', true));
        assert!(is_rna_complement('T', 'G', true));
    }

    #[test]
    fn test_khs2ct_simple() {
        // Simple hairpin: >>>...<<<
        let ss = ">>>...<<<";
        let ct = khs2ct(ss, false).unwrap();

        assert_eq!(ct.len(), 9);
        assert_eq!(ct[0], 8); // First > pairs with last <
        assert_eq!(ct[8], 0);
        assert_eq!(ct[1], 7);
        assert_eq!(ct[7], 1);
        assert_eq!(ct[2], 6);
        assert_eq!(ct[6], 2);

        // Unpaired positions
        assert_eq!(ct[3], -1);
        assert_eq!(ct[4], -1);
        assert_eq!(ct[5], -1);
    }

    #[test]
    fn test_khs2ct_nested() {
        // Nested structure: >>..>>..<<..<<
        let ss = ">>..>>..<<..<<";
        let ct = khs2ct(ss, false).unwrap();

        // Outer pairs
        assert_eq!(ct[0], 13);
        assert_eq!(ct[13], 0);
        assert_eq!(ct[1], 12);
        assert_eq!(ct[12], 1);

        // Inner pairs
        assert_eq!(ct[4], 9);
        assert_eq!(ct[9], 4);
        assert_eq!(ct[5], 8);
        assert_eq!(ct[8], 5);
    }

    #[test]
    fn test_khs2ct_errors() {
        // Unmatched > (more > than <)
        assert!(khs2ct(">>.<", false).is_err());

        // Unmatched < (more < than >)
        assert!(khs2ct("><<<", false).is_err());

        // Empty string is OK
        assert!(khs2ct("", false).is_ok());
    }

    #[test]
    fn test_khs2ct_pseudoknots() {
        // Pseudoknot: >AAA>...<aaa<
        // Main structure: > >   < <
        // Pseudoknot:      AAA aaa (stack is LIFO, so 3-9, 2-10, 1-11)
        let ss = ">AAA>...<aaa<";
        let ct = khs2ct(ss, true).unwrap();

        // Main structure pairs
        assert_eq!(ct[0], 12);
        assert_eq!(ct[12], 0);
        assert_eq!(ct[4], 8);
        assert_eq!(ct[8], 4);

        // Pseudoknot pairs (LIFO order)
        assert_eq!(ct[3], 9); // Last A pushed pairs with first a popped
        assert_eq!(ct[9], 3);
        assert_eq!(ct[2], 10);
        assert_eq!(ct[10], 2);
        assert_eq!(ct[1], 11); // First A pushed pairs with last a popped
        assert_eq!(ct[11], 1);
    }

    #[test]
    fn test_khs2ct_ignore_pseudoknots() {
        // Same string, but pseudoknots disabled
        let ss = ">AAA>...<aaa<";
        let ct = khs2ct(ss, false).unwrap();

        // Only main structure pairs
        assert_eq!(ct[0], 12);
        assert_eq!(ct[4], 8);

        // Pseudoknot positions are unpaired
        assert_eq!(ct[1], -1);
        assert_eq!(ct[2], -1);
        assert_eq!(ct[3], -1);
        assert_eq!(ct[9], -1);
        assert_eq!(ct[10], -1);
        assert_eq!(ct[11], -1);
    }

    #[test]
    fn test_trace2khs_simple() {
        use crate::types::trace::Trace;

        // Build simple trace tree with one MATP pair
        // Sequence: AUGC (positions 1-4, 0-indexed in seq array)
        let seq = b"AUGC";

        // MATP at positions 1 and 4 (1-based, so indices 0 and 3)
        let trace = Trace::with_left(
            0, 0, 0, U_BEGIN_ST,
            Trace::leaf(1, 4, 1, U_MATP_ST),
        );

        let ss = trace2khs(&trace, seq, 4, false);
        assert_eq!(ss, ">..<");
    }

    #[test]
    fn test_trace2khs_watson_crick_only() {
        use crate::types::trace::Trace;

        // Sequence with non-Watson-Crick pair: AGUC
        let seq = b"AGUC";

        // MATP at positions 1 and 4 (A-C, not a valid pair)
        let trace = Trace::with_left(
            0, 0, 0, U_BEGIN_ST,
            Trace::leaf(1, 4, 1, U_MATP_ST),
        );

        // Without watsoncrick check, pair is marked
        let ss1 = trace2khs(&trace, seq, 4, false);
        assert_eq!(ss1, ">..<");

        // With watsoncrick check, pair is not marked (A-C invalid)
        let ss2 = trace2khs(&trace, seq, 4, true);
        assert_eq!(ss2, "....");
    }

    #[test]
    fn test_trace2khs_multiple_pairs() {
        use crate::types::trace::Trace;

        // Sequence: AUGCAU (valid pairs: A-U at 0-5, U-A at 1-4)
        let seq = b"AUGCAU";

        // Build tree with two MATP nodes
        let left = Trace::leaf(1, 6, 1, U_MATP_ST);  // A-U at positions 1,6
        let right = Trace::leaf(2, 5, 2, U_MATP_ST); // U-A at positions 2,5
        let root = Trace::with_left(
            0, 0, 0, U_BEGIN_ST,
            Trace::bifurc(1, 6, 3, left, right),
        );

        let ss = trace2khs(&root, seq, 6, true);
        assert_eq!(ss, ">>..<<");
    }

    #[test]
    fn test_is_gap() {
        assert!(is_gap('-'));
        assert!(is_gap('_'));
        assert!(is_gap('~'));
        assert!(is_gap(' '));

        assert!(!is_gap('A'));
        assert!(!is_gap('.'));
        assert!(!is_gap('>'));
    }
}
