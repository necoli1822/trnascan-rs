mod common;

use common::load_golden_values;
use trnascan_rs::util::sre_math::{sre_srandom, sre_random};

#[test]
fn test_random_sequence_1000() {
    // Load golden values from C implementation
    let golden_file = "tests/golden/random_sequence/random_1000.txt";
    let lines = load_golden_values(golden_file);

    // Parse expected random values from format "random[N] = X.XXXXXXXXXX"
    let mut expected_values = Vec::new();
    for line in lines {
        let parts: Vec<&str> = line.split(" = ").collect();
        if parts.len() == 2 {
            if let Ok(value) = parts[1].parse::<f64>() {
                expected_values.push(value);
            }
        }
    }

    // Verify we have exactly 1000 values
    assert_eq!(
        expected_values.len(),
        1000,
        "Expected 1000 random values, got {}",
        expected_values.len()
    );

    // Initialize with seed 666 (as used in golden generation)
    sre_srandom(666);

    // Generate 1000 random values and compare
    for (i, expected) in expected_values.iter().enumerate() {
        let actual = sre_random();
        let diff = (actual - expected).abs();

        // Allow small floating point tolerance
        assert!(
            diff < 1e-9,
            "Random value mismatch at index {}: expected {:.10}, got {:.10}, diff {:.10}",
            i, expected, actual, diff
        );
    }

    println!("✓ All 1000 random values match exactly!");
}

#[test]
fn test_random_sequence_basic() {
    // Load golden values
    let golden_file = "tests/golden/random_sequence/random_basic.txt";

    // Check if file exists, if not skip this test
    if std::fs::metadata(golden_file).is_err() {
        println!("Skipping test_random_sequence_basic - golden file not found");
        return;
    }

    let lines = load_golden_values(golden_file);

    // Parse seed and expected values
    // Format: first line is seed, remaining lines are expected random values
    if lines.is_empty() {
        panic!("No data in golden file");
    }

    let seed: i64 = lines[0].parse().expect("Failed to parse seed");
    let mut expected_values = Vec::new();

    for line in &lines[1..] {
        let parts: Vec<&str> = line.split(" = ").collect();
        if parts.len() == 2 {
            if let Ok(value) = parts[1].parse::<f64>() {
                expected_values.push(value);
            }
        } else if let Ok(value) = line.parse::<f64>() {
            expected_values.push(value);
        }
    }

    assert!(!expected_values.is_empty(), "No random values in golden file");

    // Initialize with specified seed
    sre_srandom(seed);

    // Generate and verify values
    for (i, expected) in expected_values.iter().enumerate() {
        let actual = sre_random();
        let diff = (actual - expected).abs();

        assert!(
            diff < 1e-9,
            "Random value mismatch at index {}: expected {:.10}, got {:.10}, diff {:.10}",
            i, expected, actual, diff
        );
    }

    println!("✓ {} random values match for seed {}", expected_values.len(), seed);
}

#[test]
fn test_random_critical_values() {
    // Test the critical values mentioned in requirements
    sre_srandom(666);

    // Collect all 1000 values
    let mut values = Vec::new();
    for _ in 0..1000 {
        values.push(sre_random());
    }

    let r0 = values[0];
    let r1 = values[1];
    let r999 = values[999];

    assert!((r0 - 0.5948935151).abs() < 1e-9, "random[0] mismatch: got {:.10}", r0);
    assert!((r1 - 0.4373151660).abs() < 1e-9, "random[1] mismatch: got {:.10}", r1);
    assert!((r999 - 0.1448028088).abs() < 1e-9, "random[999] mismatch: got {:.10}", r999);

    println!("✓ Critical random values verified");
}
