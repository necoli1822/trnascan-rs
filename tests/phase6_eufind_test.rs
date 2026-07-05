// Phase 6: EufindtRNA tests
// Tests for tRNA detection heuristics (GetBbox, GetBestABox, etc.)

mod common;

use common::{load_golden_values, float_eq};
use trnascan_rs::eufind::{get_bbox, int_encode_seq, BBOX_CUTOFF};

#[test]
fn test_getbbox_scores() {
    // Test GetBbox scoring function for B-box detection
    // Golden file: golden/eufind/bbox_scores.txt
    let golden = load_golden_values("tests/golden/eufind/bbox_scores.txt");

    // Test 1: tRNA-Phe sequence
    let seq = "TCGCTGTTAGTTACCATCGCACGGATGGCCGAGTGGTCTAAGGCGCCAGACTCAAGCGAAATGCTTGCCTCATGCTCGAGGTCGACTGGGTGTTCTGGTACTCGTATGGGTGCGTGGGTTCGAATCCCACTTCGTGCAG";
    let iseq = int_encode_seq(seq);

    // Expected B-box at position 116 with score -1.9170
    let score = get_bbox(&iseq, 116);
    assert!(
        float_eq(score as f64, -1.917, 0.001),
        "B-box score at pos 116: expected -1.917, got {}",
        score
    );
    assert!(score > BBOX_CUTOFF, "B-box score should pass cutoff");

    assert!(!golden.is_empty(), "Golden file should contain B-box scores");
}

#[test]
fn test_getbestabox_scores() {
    // Test GetBestABox for A-box detection and scoring
    // Golden file: golden/eufind/abox_scores.txt
    let golden = load_golden_values("tests/golden/eufind/abox_scores.txt");

    use trnascan_rs::eufind::{get_best_abox, AB_BOX_DIST_RANGE};

    // tRNA-Phe sequence
    let seq = "TCGCTGTTAGTTACCATCGCACGGATGGCCGAGTGGTCTAAGGCGCCAGACTCAAGCGAAATGCTTGCCTCATGCTCGAGGTCGACTGGGTGTTCTGGTACTCGTATGGGTGCGTGGGTTCGAATCCCACTTCGTGCAG";
    let iseq = int_encode_seq(seq);
    let seq_bytes = seq.as_bytes();

    // B-box at position 116, end at 126
    let bbox_st = 116;
    let max_ab_dist = AB_BOX_DIST_RANGE;

    let (abox_st, abox_end, _abox_gap, abox_sc, abdist_sc) =
        get_best_abox(&iseq, seq_bytes, bbox_st, max_ab_dist, -1);

    // Expected values from golden file:
    // Position: 24-43, A-box score: -13.764, AB distance: 72, AB dist score: -5.442
    assert_eq!(abox_st, 24, "A-box start position");
    assert_eq!(abox_end, 43, "A-box end position");
    assert!(
        float_eq(abox_sc as f64, -13.764, 0.01),
        "A-box score: expected -13.764, got {}",
        abox_sc
    );
    assert!(
        float_eq(abdist_sc as f64, -5.442, 0.01),
        "AB dist score: expected -5.442, got {}",
        abdist_sc
    );

    // Verify AB distance
    let ab_dist = bbox_st - abox_end - 1;
    assert_eq!(ab_dist, 72, "AB distance should be 72");

    assert!(!golden.is_empty(), "Golden file should contain A-box scores");
}

#[test]
fn test_trna_detection_full() {
    // Test complete tRNA detection pipeline
    // Golden file: golden/eufind/trna_detection.txt
    let golden = load_golden_values("tests/golden/eufind/trna_detection.txt");

    use trnascan_rs::eufind::{get_best_abox, get_best_trx_term, TrnaInfo, AB_BOX_DIST_RANGE};

    // tRNA-Phe sequence
    let seq = "TCGCTGTTAGTTACCATCGCACGGATGGCCGAGTGGTCTAAGGCGCCAGACTCAAGCGAAATGCTTGCCTCATGCTCGAGGTCGACTGGGTGTTCTGGTACTCGTATGGGTGCGTGGGTTCGAATCCCACTTCGTGCAG";
    let iseq = int_encode_seq(seq);
    let seq_bytes = seq.as_bytes();
    let mut trna = TrnaInfo::new();

    // Step 1: B-box detection (at position 116)
    trna.bbox_st = 116;
    trna.bbox_end = 126;
    trna.bbox_sc = get_bbox(&iseq, 116);

    // Step 2: A-box detection
    let (abox_st, abox_end, _abox_gap, abox_sc, abdist_sc) =
        get_best_abox(&iseq, seq_bytes, trna.bbox_st, AB_BOX_DIST_RANGE, -1);
    trna.abox_st = abox_st;
    trna.abox_end = abox_end;
    trna.abox_sc = abox_sc;
    trna.abdist_sc = abdist_sc;

    // Step 3: Termination signal
    let (term_st, term_sc) = get_best_trx_term(seq_bytes, trna.bbox_end, seq.len());
    trna.term_st = term_st;
    trna.term_sc = term_sc;

    // Step 4: Calculate total score
    trna.tot_sc = trna.abox_sc + trna.bbox_sc + trna.abdist_sc + trna.term_sc;

    // Verify against golden values
    assert_eq!(trna.abox_st, 24, "A-box start");
    assert_eq!(trna.abox_end, 43, "A-box end");
    assert_eq!(trna.bbox_st, 116, "B-box start");
    assert_eq!(trna.bbox_end, 126, "B-box end");
    assert_eq!(trna.term_st, -1, "No TTTT termination found");

    assert!(
        float_eq(trna.abox_sc as f64, -13.764, 0.01),
        "A-box score: expected -13.764, got {}",
        trna.abox_sc
    );
    assert!(
        float_eq(trna.bbox_sc as f64, -1.917, 0.01),
        "B-box score: expected -1.917, got {}",
        trna.bbox_sc
    );
    assert!(
        float_eq(trna.abdist_sc as f64, -5.442, 0.01),
        "AB dist score: expected -5.442, got {}",
        trna.abdist_sc
    );
    assert!(
        float_eq(trna.term_sc as f64, -0.55, 0.01),
        "Term score: expected -0.55, got {}",
        trna.term_sc
    );
    assert!(
        float_eq(trna.tot_sc as f64, -21.673, 0.01),
        "Total score: expected -21.673, got {}",
        trna.tot_sc
    );

    assert!(!golden.is_empty(), "Golden file should contain tRNA detection results");
}

#[test]
#[ignore] // Pending eufind module implementation
fn test_astem_scoring() {
    // Test A-stem scoring function

    // TODO: When eufind module is implemented:
    // - Score known A-stem sequences
    // - Verify base pairing contributions
    // - Test with mismatches and wobble pairs
}

#[test]
#[ignore] // Pending eufind module implementation
fn test_dstem_scoring() {
    // Test D-stem scoring function

    // TODO: When eufind module is implemented:
    // - Score known D-stem sequences
    // - Verify stem length detection
    // - Test variable D-loop sizes
}

#[test]
#[ignore] // Pending eufind module implementation
fn test_anticodon_stem_scoring() {
    // Test anticodon stem scoring

    // TODO: When eufind module is implemented:
    // - Score anticodon stem regions
    // - Verify anticodon loop detection
    // - Test with various anticodon sequences
}

#[test]
#[ignore] // Pending eufind module implementation
fn test_tstem_scoring() {
    // Test T-stem scoring function

    // TODO: When eufind module is implemented:
    // - Score T-stem regions
    // - Verify T-loop consensus matches
    // - Test variable loop sizes
}

#[test]
#[ignore] // Pending eufind module implementation
fn test_spacer_distance_calculation() {
    // Test spacer distance calculations between boxes

    // TODO: When eufind module is implemented:
    // - Calculate A-box to B-box spacer
    // - Verify distance constraints
    // - Test boundary cases (min/max spacer lengths)
}

#[test]
#[ignore] // Pending eufind module implementation
fn test_orientation_detection() {
    // Test detection of tRNA orientation (forward/reverse strand)

    // TODO: When eufind module is implemented:
    // - Detect tRNAs on forward strand
    // - Detect tRNAs on reverse strand
    // - Verify coordinates are correctly oriented
}

#[test]
#[ignore] // Pending eufind module implementation
fn test_false_positive_filtering() {
    // Test filtering of false positive tRNA predictions

    // TODO: When eufind module is implemented:
    // - Run on sequences known to produce false positives
    // - Verify score thresholds filter appropriately
    // - Compare with golden expected detections
}
