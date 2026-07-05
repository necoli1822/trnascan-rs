// Phase 7: Secondary structure tests
// Tests for converting tRNA alignments to secondary structure notation

mod common;

use trnascan_rs::structure::{is_rna_complement, khs2ct};

#[test]
fn test_is_rna_complement_watson_crick() {
    // Watson-Crick pairs (allow_gu=FALSE)
    assert!(is_rna_complement('A', 'U', false));
    assert!(is_rna_complement('U', 'A', false));
    assert!(is_rna_complement('G', 'C', false));
    assert!(is_rna_complement('C', 'G', false));
    assert!(!is_rna_complement('G', 'U', false), "G-U should be FALSE without wobble");

    // Non-complementary pairs
    assert!(!is_rna_complement('A', 'A', false));
    assert!(!is_rna_complement('A', 'G', false));
    assert!(!is_rna_complement('C', 'C', false));
}

#[test]
fn test_is_rna_complement_wobble() {
    // With GU wobble pairs (allow_gu=TRUE)
    assert!(is_rna_complement('A', 'U', true));
    assert!(is_rna_complement('G', 'C', true));
    assert!(is_rna_complement('G', 'U', true), "G-U should be TRUE with wobble");
    assert!(is_rna_complement('U', 'G', true), "U-G should be TRUE with wobble");
}

#[test]
fn test_is_rna_complement_case_insensitive() {
    // Case insensitive
    assert!(is_rna_complement('a', 'u', false));
    assert!(is_rna_complement('G', 'c', false));
}

#[test]
fn test_is_rna_complement_dna_compatibility() {
    // DNA compatibility (T converted to U)
    assert!(is_rna_complement('A', 'T', false));
    assert!(is_rna_complement('T', 'A', false));
    assert!(is_rna_complement('G', 'T', true), "G-T should be TRUE with wobble (T as U)");
}

#[test]
fn test_khs2ct_simple_hairpin() {
    // Test case 1: Simple hairpin
    // Structure: >>>....<<<
    let ss = ">>>....<<<";
    let ct = khs2ct(ss, false).expect("Failed to parse simple hairpin");

    assert_eq!(ct.len(), 10);
    assert_eq!(ct[0], 9);
    assert_eq!(ct[1], 8);
    assert_eq!(ct[2], 7);
    assert_eq!(ct[3], -1); // Unpaired
    assert_eq!(ct[4], -1); // Unpaired
    assert_eq!(ct[5], -1); // Unpaired
    assert_eq!(ct[6], -1); // Unpaired
    assert_eq!(ct[7], 2);
    assert_eq!(ct[8], 1);
    assert_eq!(ct[9], 0);
}

#[test]
fn test_khs2ct_pseudoknot() {
    // Test case 3: Pseudoknot structure
    // Structure: >>>AAA<<<aaa
    let ss = ">>>AAA<<<aaa";
    let ct = khs2ct(ss, true).expect("Failed to parse pseudoknot");

    assert_eq!(ct.len(), 12);

    // Main structure pairs (stack 0)
    assert_eq!(ct[0], 8);
    assert_eq!(ct[1], 7);
    assert_eq!(ct[2], 6);
    assert_eq!(ct[6], 2);
    assert_eq!(ct[7], 1);
    assert_eq!(ct[8], 0);

    // Pseudoknot pairs (stack 1 for 'A'/'a')
    // Note: LIFO order, so last A pushed pairs with first a popped
    assert_eq!(ct[3], 11);
    assert_eq!(ct[4], 10);
    assert_eq!(ct[5], 9);
    assert_eq!(ct[9], 5);
    assert_eq!(ct[10], 4);
    assert_eq!(ct[11], 3);
}

#[test]
fn test_khs2ct_errors() {
    // Unmatched > (more > than <)
    assert!(khs2ct(">>.<", false).is_err(), "Should fail on unmatched >");

    // Unmatched < (more < than >)
    assert!(khs2ct("><<<", false).is_err(), "Should fail on unmatched <");
}

#[test]
#[ignore] // Integration test - requires traceback data
fn test_trace2khs_trna_sequence() {
    // Test Trace2KHS with actual tRNA sequence
    // This requires building a proper trace tree from Viterbi alignment

    // Test sequence (76 nt):
    // GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA

    // Expected Watson-Crick only (watsoncrick=TRUE):
    // ..>........>...........<.........................>..............<........<..

    // Expected all pairs (watsoncrick=FALSE):
    // >>>.......>>...........<<..>>.............<<.....>>............<<........<<<

    // TODO: Build trace tree from golden file and verify output
}

#[test]
#[ignore] // Pending alignment module implementation
fn test_structure_standard_trna() {
    // Test structure conversion for standard cloverleaf tRNA

    // TODO: When alignment module is implemented:
    // - Load Example1 tRNA results
    // - Convert to secondary structure
    // - Verify acceptor stem pairing
    // - Verify D-stem/loop structure
    // - Verify anticodon stem/loop structure
    // - Verify T-stem/loop structure
}

#[test]
#[ignore] // Future feature
fn test_structure_with_introns() {
    // Test structure conversion for tRNAs with introns

    // TODO: When structure module is extended:
    // - Handle tRNAs with intron sequences
    // - Verify intron positions marked correctly
    // - Check structure excludes intron bases
}

#[test]
#[ignore] // Future feature
fn test_structure_unusual_trna() {
    // Test structure for unusual tRNA types (e.g., missing D-loop)

    // TODO: When structure module is extended:
    // - Convert tRNA-Sec or tRNA-Ile with non-standard structure
    // - Verify handling of variable loops
    // - Check shortened or extended stems
}
