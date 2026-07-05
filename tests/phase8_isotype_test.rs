// Phase 8: Isotype classification tests
// Tests for determining tRNA isotype (amino acid type) from anticodon and sequence

mod common;

use trnascan_rs::isotype::{anticodon_to_isotype, get_all_anticodons_for_isotype};

#[test]
fn test_isotype_from_anticodon() {
    // Test basic isotype determination from anticodon
    assert_eq!(anticodon_to_isotype("CAT"), Some("Met"));
    assert_eq!(anticodon_to_isotype("AAG"), Some("Leu"));
    assert_eq!(anticodon_to_isotype("TGT"), Some("Thr"));
    assert_eq!(anticodon_to_isotype("CGA"), Some("Ser"));
    assert_eq!(anticodon_to_isotype("GAA"), Some("Phe"));
    assert_eq!(anticodon_to_isotype("CTT"), Some("Lys"));
    assert_eq!(anticodon_to_isotype("TCA"), Some("SeC"));

    // Unknown anticodon
    assert_eq!(anticodon_to_isotype("XYZ"), None);
}

#[test]
fn test_isotype_all_standard() {
    // Test all 20 standard amino acids plus SeC

    // Alanine
    assert_eq!(anticodon_to_isotype("AGC"), Some("Ala"));
    assert_eq!(anticodon_to_isotype("TGC"), Some("Ala"));

    // Arginine
    assert_eq!(anticodon_to_isotype("ACG"), Some("Arg"));
    assert_eq!(anticodon_to_isotype("TCT"), Some("Arg"));

    // Asparagine
    assert_eq!(anticodon_to_isotype("GTT"), Some("Asn"));

    // Aspartate
    assert_eq!(anticodon_to_isotype("GTC"), Some("Asp"));

    // Cysteine
    assert_eq!(anticodon_to_isotype("GCA"), Some("Cys"));

    // Glutamine
    assert_eq!(anticodon_to_isotype("CTG"), Some("Gln"));
    assert_eq!(anticodon_to_isotype("TTG"), Some("Gln"));

    // Glutamate
    assert_eq!(anticodon_to_isotype("CTC"), Some("Glu"));
    assert_eq!(anticodon_to_isotype("TTC"), Some("Glu"));

    // Glycine
    assert_eq!(anticodon_to_isotype("ACC"), Some("Gly"));
    assert_eq!(anticodon_to_isotype("TCC"), Some("Gly"));

    // Histidine
    assert_eq!(anticodon_to_isotype("GTG"), Some("His"));

    // Isoleucine
    assert_eq!(anticodon_to_isotype("AAT"), Some("Ile"));
    assert_eq!(anticodon_to_isotype("GAT"), Some("Ile"));

    // Leucine
    assert_eq!(anticodon_to_isotype("CAA"), Some("Leu"));
    assert_eq!(anticodon_to_isotype("AAG"), Some("Leu"));

    // Lysine
    assert_eq!(anticodon_to_isotype("CTT"), Some("Lys"));
    assert_eq!(anticodon_to_isotype("TTT"), Some("Lys"));

    // Methionine
    assert_eq!(anticodon_to_isotype("CAT"), Some("Met"));

    // Phenylalanine
    assert_eq!(anticodon_to_isotype("GAA"), Some("Phe"));

    // Proline
    assert_eq!(anticodon_to_isotype("AGG"), Some("Pro"));
    assert_eq!(anticodon_to_isotype("TGG"), Some("Pro"));

    // Serine
    assert_eq!(anticodon_to_isotype("ACT"), Some("Ser"));
    assert_eq!(anticodon_to_isotype("CGA"), Some("Ser"));

    // Threonine
    assert_eq!(anticodon_to_isotype("AGT"), Some("Thr"));
    assert_eq!(anticodon_to_isotype("TGT"), Some("Thr"));

    // Tryptophan
    assert_eq!(anticodon_to_isotype("CCA"), Some("Trp"));

    // Tyrosine
    assert_eq!(anticodon_to_isotype("GTA"), Some("Tyr"));

    // Valine
    assert_eq!(anticodon_to_isotype("AAC"), Some("Val"));
    assert_eq!(anticodon_to_isotype("TAC"), Some("Val"));

    // Selenocysteine
    assert_eq!(anticodon_to_isotype("TCA"), Some("SeC"));
}

#[test]
#[ignore] // Genetic code variants not yet implemented
fn test_isotype_genetic_code_variants() {
    // Test isotype determination with different genetic codes
    // Current implementation only supports standard genetic code
    // Future work: support mitochondrial and other genetic code variants
}

#[test]
fn test_isotype_ambiguous_anticodons() {
    // Test handling of ambiguous or modified anticodons
    // Unknown anticodons return None
    assert_eq!(anticodon_to_isotype("NNN"), None);
    assert_eq!(anticodon_to_isotype("XYZ"), None);
    assert_eq!(anticodon_to_isotype(""), None);

    // Wobble position variants should still work if defined
    assert_eq!(anticodon_to_isotype("AAG"), Some("Leu")); // Wobble
    assert_eq!(anticodon_to_isotype("TCT"), Some("Arg")); // Wobble
}

#[test]
fn test_isotype_output_formatting() {
    // Test isotype output in various formats
    // Currently returns three-letter code
    let isotype = anticodon_to_isotype("CAT").unwrap();
    assert_eq!(isotype, "Met");
    assert_eq!(isotype.len(), 3); // Three-letter code

    let isotype = anticodon_to_isotype("TCA").unwrap();
    assert_eq!(isotype, "SeC");
    assert_eq!(isotype.len(), 3); // Three-letter code
}

#[test]
fn test_initiator_methionine() {
    // Test detection of initiator tRNA-Met vs elongator tRNA-Met
    // Both use CAT anticodon - differentiation requires structure analysis

    // Basic anticodon mapping gives Met
    assert_eq!(anticodon_to_isotype("CAT"), Some("Met"));

    // Get all anticodons for iMet
    let imet_anticodons = get_all_anticodons_for_isotype("iMet");
    assert!(!imet_anticodons.is_empty());
    assert!(imet_anticodons.contains(&"CAT"));

    // Get all anticodons for fMet
    let fmet_anticodons = get_all_anticodons_for_isotype("fMet");
    assert!(!fmet_anticodons.is_empty());
    assert!(fmet_anticodons.contains(&"CAT"));
}

#[test]
fn test_suppressor_trnas() {
    // Test handling of suppressor tRNAs (amber, ochre, opal)

    // Amber suppressor (UAG stop codon)
    assert_eq!(anticodon_to_isotype("CTA"), Some("Sup"));

    // Ochre suppressor (UAA stop codon)
    assert_eq!(anticodon_to_isotype("TTA"), Some("Sup"));

    // Opal suppressor (UGA stop codon) - also SeC
    assert_eq!(anticodon_to_isotype("TCA"), Some("SeC"));

    // Get all suppressor anticodons
    let sup_anticodons = get_all_anticodons_for_isotype("Sup");
    assert!(!sup_anticodons.is_empty());
    assert!(sup_anticodons.contains(&"CTA"));
    assert!(sup_anticodons.contains(&"TTA"));
}

#[test]
#[ignore] // Pseudogene detection requires full scorer implementation
fn test_pseudo_trnas() {
    // Test classification of pseudogenes
    // Requires IsotypeScorer with CM scoring to detect pseudogenes
    // based on low scores and structural anomalies
}

#[test]
fn test_isotype_from_example_files() {
    // Test isotype classification on Example1 and Example2 data
    // Based on golden file: tests/golden/isotype/isotype_assignments.txt

    // CELF22B7.trna1: AAG -> Leu
    assert_eq!(anticodon_to_isotype("AAG"), Some("Leu"));

    // CELF22B7.trna2: CGA -> Ser
    assert_eq!(anticodon_to_isotype("CGA"), Some("Ser"));

    // CELF22B7.trna3: GAA -> Phe
    assert_eq!(anticodon_to_isotype("GAA"), Some("Phe"));

    // CELF22B7.trna4: GAA -> Phe
    assert_eq!(anticodon_to_isotype("GAA"), Some("Phe"));

    // CELF22B7.trna5: TGG -> Pro
    assert_eq!(anticodon_to_isotype("TGG"), Some("Pro"));

    // MySeq1.trna1: TGT -> Thr
    assert_eq!(anticodon_to_isotype("TGT"), Some("Thr"));

    // MySeq2.trna1: TCT -> Arg
    assert_eq!(anticodon_to_isotype("TCT"), Some("Arg"));

    // MySeq3.trna1: CGA -> Ser
    assert_eq!(anticodon_to_isotype("CGA"), Some("Ser"));

    // MySeq4.trna1: AAG -> Leu
    assert_eq!(anticodon_to_isotype("AAG"), Some("Leu"));

    // MySeq5.trna1: TCA -> SeC
    assert_eq!(anticodon_to_isotype("TCA"), Some("SeC"));

    // MySeq6.trna1: CTT -> Lys
    assert_eq!(anticodon_to_isotype("CTT"), Some("Lys"));
}
