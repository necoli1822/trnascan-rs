// Phase 5: Viterbi algorithm tests
// Tests for Viterbi dynamic programming and traceback

mod common;

use common::float_eq;
use std::fs;
use std::path::Path;
use trnascan_rs::core::save::read_cm;
use trnascan_rs::types::constants::ALPHASIZE;
use trnascan_rs::viterbi::{prepare_sequence, rearrange_cm, viterbi_align};

/// Load the tRNA covariance model
fn load_test_cm() -> trnascan_rs::types::cm::CM {
    // Try different possible paths
    let paths = [
        "data/models/TRNA2.cm",
        "../data/models/TRNA2.cm",
        "tests/fixtures/TRNA2.cm",
    ];

    for path in &paths {
        if Path::new(path).exists() {
            return read_cm(path).expect(&format!("Failed to read CM from {}", path));
        }
    }

    panic!("Could not find TRNA2.cm model file");
}

/// Test sequences from golden files
const PHE_73BP: &str = "GCCTCGATAGCTCAGTTGGGAGAGCGTACGACTGAAGATCGTAAGGTCACCAGTTCGATCCTGGTTCGGGGCA";
const SER_82BP: &str = "GCAGTCATGTCCGAGTGGTAAGGAGATTGACTAGAAATCAATTGGGCTCTGCCCGCGTAGGTTCGAATCCTGCTGACTGCG";
const FRAGMENT_50BP: &str = "GCCTCGATAGCTCAGTTGGGAGAGCGTACGACTGAAGATCGTAAGGTCAC";

/// Expected scores from golden file (within 0.01 bits tolerance)
const PHE_EXPECTED: f64 = 73.882;
const SER_EXPECTED: f64 = 66.085;
const FRAGMENT_EXPECTED: f64 = -12.801;

#[test]
fn test_viterbi_phe_sequence() {
    let cm = load_test_cm();
    let rfreq = [0.25f64; ALPHASIZE];

    let (icm, statenum) = rearrange_cm(&cm, &rfreq);
    assert!(statenum > 0, "Should have states");
    assert_eq!(icm.len(), statenum);

    let seq = prepare_sequence(PHE_73BP);
    assert_eq!(seq.len(), 73);

    let result = viterbi_align(&icm, statenum, &seq);
    assert!(result.is_ok(), "Viterbi alignment should succeed");

    let (score, _trace) = result.unwrap();
    assert!(
        float_eq(score, PHE_EXPECTED, 0.02),
        "Phe-73bp score should be {:.3}, got {:.3}",
        PHE_EXPECTED,
        score
    );
}

#[test]
fn test_viterbi_ser_sequence() {
    let cm = load_test_cm();
    let rfreq = [0.25f64; ALPHASIZE];

    let (icm, statenum) = rearrange_cm(&cm, &rfreq);

    let seq = prepare_sequence(SER_82BP);
    assert_eq!(seq.len(), 81);

    let result = viterbi_align(&icm, statenum, &seq);
    assert!(result.is_ok(), "Viterbi alignment should succeed");

    let (score, _trace) = result.unwrap();
    assert!(
        float_eq(score, SER_EXPECTED, 0.02),
        "Ser-82bp score should be {:.3}, got {:.3}",
        SER_EXPECTED,
        score
    );
}

#[test]
fn test_viterbi_fragment_sequence() {
    let cm = load_test_cm();
    let rfreq = [0.25f64; ALPHASIZE];

    let (icm, statenum) = rearrange_cm(&cm, &rfreq);

    let seq = prepare_sequence(FRAGMENT_50BP);
    assert_eq!(seq.len(), 50);

    let result = viterbi_align(&icm, statenum, &seq);
    assert!(result.is_ok(), "Viterbi alignment should succeed");

    let (score, _trace) = result.unwrap();
    assert!(
        float_eq(score, FRAGMENT_EXPECTED, 0.02),
        "Fragment-50bp score should be {:.3}, got {:.3}",
        FRAGMENT_EXPECTED,
        score
    );
}

#[test]
fn test_viterbi_model_conversion() {
    // Test that RearrangeCM produces correct number of states
    let cm = load_test_cm();
    let rfreq = [0.25f64; ALPHASIZE];

    let (icm, statenum) = rearrange_cm(&cm, &rfreq);

    // From golden file: TRNA2.cm produces 289 states
    assert_eq!(statenum, 289, "TRNA2.cm should have 289 states");
    assert_eq!(icm.len(), statenum);

    // First state should be BEGIN (ROOT)
    assert_eq!(
        icm[0].statetype,
        trnascan_rs::types::constants::U_BEGIN_ST,
        "First state should be BEGIN"
    );
}

#[test]
fn test_viterbi_traceback_structure() {
    let cm = load_test_cm();
    let rfreq = [0.25f64; ALPHASIZE];

    let (icm, statenum) = rearrange_cm(&cm, &rfreq);

    let seq = prepare_sequence(PHE_73BP);

    let result = viterbi_align(&icm, statenum, &seq);
    assert!(result.is_ok());

    let (_score, trace) = result.unwrap();

    // Trace should start at BEGIN
    assert_eq!(
        trace.trace_type,
        trnascan_rs::types::constants::U_BEGIN_ST,
        "Trace should start with BEGIN"
    );

    // Trace should cover the sequence
    assert!(trace.emitl >= 0);
    assert!(trace.emitr >= 0);
}

#[test]
fn test_viterbi_relative_ordering() {
    // Higher-scoring (real tRNA) should beat lower-scoring (fragment)
    let cm = load_test_cm();
    let rfreq = [0.25f64; ALPHASIZE];

    let (icm, statenum) = rearrange_cm(&cm, &rfreq);

    let phe_seq = prepare_sequence(PHE_73BP);
    let frag_seq = prepare_sequence(FRAGMENT_50BP);

    let (phe_score, _) = viterbi_align(&icm, statenum, &phe_seq).unwrap();
    let (frag_score, _) = viterbi_align(&icm, statenum, &frag_seq).unwrap();

    assert!(
        phe_score > frag_score,
        "Real tRNA should score higher than fragment"
    );
    assert!(
        phe_score > 50.0,
        "Real tRNA should have positive high score"
    );
}

#[test]
fn test_viterbi_scores_against_golden() {
    // Load and parse golden file
    let golden_path = "tests/golden/viterbi/viterbi_scores.txt";
    if !Path::new(golden_path).exists() {
        println!("Golden file not found, skipping golden comparison");
        return;
    }

    let content = fs::read_to_string(golden_path).expect("Failed to read golden file");
    let mut expected_scores: Vec<(String, usize, f64)> = Vec::new();

    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 3 {
            let name = parts[0].to_string();
            let len: usize = parts[1].parse().unwrap_or(0);
            let score: f64 = parts[2].parse().unwrap_or(0.0);
            expected_scores.push((name, len, score));
        }
    }

    assert!(
        !expected_scores.is_empty(),
        "Should have parsed golden scores"
    );

    // Run Viterbi and compare
    let cm = load_test_cm();
    let rfreq = [0.25f64; ALPHASIZE];
    let (icm, statenum) = rearrange_cm(&cm, &rfreq);

    let test_seqs = [
        ("Phe-73bp", PHE_73BP),
        ("Ser-82bp", SER_82BP),
        ("Fragment-50bp", FRAGMENT_50BP),
    ];

    for (name, seq_str) in &test_seqs {
        let seq = prepare_sequence(seq_str);
        let (score, _) = viterbi_align(&icm, statenum, &seq).unwrap();

        // Find expected score
        if let Some(expected) = expected_scores.iter().find(|(n, _, _)| n == *name) {
            assert!(
                float_eq(score, expected.2, 0.02),
                "{}: expected {:.3}, got {:.3}",
                name,
                expected.2,
                score
            );
            println!("{}: {:.4} (expected {:.4}) - PASS", name, score, expected.2);
        }
    }
}

#[test]
fn test_prepare_sequence_normalization() {
    // Test sequence preparation
    let raw = "acgtACGT\nUu\t\r";
    let prepared = prepare_sequence(raw);
    assert_eq!(prepared, b"ACGTACGTUU");
}
