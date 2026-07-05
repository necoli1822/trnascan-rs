mod common;

use common::load_golden_values;
use trnascan_rs::types::constants::ilog2;

#[test]
fn test_ilog2_values() {
    // Load golden values from C implementation
    let golden_file = "tests/golden/constants/ilog2_values.txt";
    let lines = load_golden_values(golden_file);

    // Parse expected values from format: "ILOG2(value) = result"
    let mut test_cases = Vec::new();
    for line in lines {
        // Format: "ILOG2(0.5) = -1000"
        if line.starts_with("ILOG2(") {
            // Extract value between ( and )
            let start = line.find('(').unwrap() + 1;
            let end = line.find(')').unwrap();
            let value_str = &line[start..end];

            // Extract result after =
            let equals_pos = line.find('=').unwrap();
            let result_str = line[equals_pos + 1..].trim();

            let input: f64 = value_str.parse().expect("Failed to parse input");
            let expected: i32 = result_str.parse().expect("Failed to parse expected");
            test_cases.push((input, expected));
        }
    }

    // Verify we have data
    assert!(!test_cases.is_empty(), "No test data loaded from golden file");

    // Test actual implementation against golden values
    for (input, expected) in test_cases {
        let result = ilog2(input);
        assert_eq!(result, expected, "ilog2({}) failed: expected {}, got {}", input, expected, result);
    }
}

#[test]
#[ignore] // Implementation doesn't exist yet
fn test_srandom_init_values() {
    // Load golden values
    let golden_file = "tests/golden/constants/srandom_init_values.txt";
    let lines = load_golden_values(golden_file);

    // Parse expected seed and first few random values
    // Format: seed value1 value2 value3 ...
    let mut test_cases = Vec::new();
    for line in lines {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 2 {
            let seed: u32 = parts[0].parse().expect("Failed to parse seed");
            let values: Vec<i32> = parts[1..]
                .iter()
                .map(|v| v.parse().expect("Failed to parse random value"))
                .collect();
            test_cases.push((seed, values));
        }
    }

    // Verify we have data
    assert!(!test_cases.is_empty(), "No test data loaded from golden file");

    // Placeholder for actual implementation test
    // When implementation exists, replace this with:
    // for (seed, expected_values) in test_cases {
    //     srandom(seed);
    //     let mut actual_values = Vec::new();
    //     for _ in 0..expected_values.len() {
    //         actual_values.push(random());
    //     }
    //     assert_eq!(actual_values, expected_values,
    //                "srandom({}) sequence mismatch", seed);
    // }

    println!("Loaded {} srandom test cases from golden file", test_cases.len());
}
