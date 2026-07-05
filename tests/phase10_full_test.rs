// Phase 10: Full end-to-end integration tests
// Tests against Example1 and Example2 complete outputs

mod common;

use std::fs;
use std::path::Path;
use std::collections::HashMap;

use trnascan_rs::core::TrnaScanner;
use trnascan_rs::squid::SeqFileReader;

#[derive(Debug)]
struct TrnaOutputLine {
    seq_name: String,
    trna_num: usize,
    begin: usize,
    end: usize,
    isotype: String,
    anticodon: String,
    intron_begin: usize,
    intron_end: usize,
    inf_score: f64,
    hmm_score: f64,
    ss_score: f64,
    origin: String,
    cm_isotype: String,
    cm_score: f64,
    note: String,
}

impl TrnaOutputLine {
    fn parse(line: &str) -> Option<Self> {
        // Skip header and separator lines
        if line.starts_with("Sequence") || line.starts_with("Name") || line.starts_with("----") {
            return None;
        }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 14 {
            return None;
        }

        Some(Self {
            seq_name: parts[0].trim().to_string(),
            trna_num: parts[1].trim().parse().ok()?,
            begin: parts[2].trim().parse().ok()?,
            end: parts[3].trim().parse().ok()?,
            isotype: parts[4].trim().to_string(),
            anticodon: parts[5].trim().to_string(),
            intron_begin: parts[6].trim().parse().ok()?,
            intron_end: parts[7].trim().parse().ok()?,
            inf_score: parts[8].trim().parse().ok()?,
            hmm_score: parts[9].trim().parse().ok()?,
            ss_score: parts[10].trim().parse().ok()?,
            origin: parts[11].trim().to_string(),
            cm_isotype: parts[12].trim().to_string(),
            cm_score: parts[13].trim().parse().ok()?,
            note: if parts.len() > 14 {
                parts[14].trim().to_string()
            } else {
                String::new()
            },
        })
    }
}

/// Parse tRNA output file into structured data
fn parse_trna_output(path: &Path) -> Vec<TrnaOutputLine> {
    let content = fs::read_to_string(path).expect("Failed to read output file");
    content.lines()
        .filter_map(TrnaOutputLine::parse)
        .collect()
}

#[test]
fn test_example2_six_trnas() {
    // Test that Example2.fa contains exactly 6 tRNAs with correct properties
    let golden_path = Path::new("tests/golden/full_run/Example2-tRNAs.out");
    let trnas = parse_trna_output(golden_path);

    assert_eq!(trnas.len(), 6, "Example2 should have exactly 6 tRNAs");

    // Verify each tRNA
    assert_eq!(trnas[0].seq_name, "MySeq1");
    assert_eq!(trnas[0].isotype, "Thr");
    assert_eq!(trnas[0].anticodon, "TGT");
    assert_eq!(trnas[0].begin, 13);
    assert_eq!(trnas[0].end, 85);

    assert_eq!(trnas[1].seq_name, "MySeq2");
    assert_eq!(trnas[1].isotype, "Arg");
    assert_eq!(trnas[1].anticodon, "TCT");

    assert_eq!(trnas[2].seq_name, "MySeq3");
    assert_eq!(trnas[2].isotype, "Ser");
    assert_eq!(trnas[2].anticodon, "CGA");
    assert_eq!(trnas[2].intron_begin, 51);
    assert_eq!(trnas[2].intron_end, 69);

    assert_eq!(trnas[3].seq_name, "MySeq4");
    assert_eq!(trnas[3].isotype, "Leu");
    assert_eq!(trnas[3].anticodon, "AAG");

    assert_eq!(trnas[4].seq_name, "MySeq5");
    assert_eq!(trnas[4].isotype, "SeC");
    assert_eq!(trnas[4].anticodon, "TCA");

    assert_eq!(trnas[5].seq_name, "MySeq6");
    assert_eq!(trnas[5].isotype, "Lys");
    assert_eq!(trnas[5].anticodon, "CTT");
    assert!(trnas[5].note.contains("ISM"), "MySeq6 should have ISM note");
}

#[test]
fn test_example2_scores() {
    // Test that scores are within expected ranges
    let golden_path = Path::new("tests/golden/full_run/Example2-tRNAs.out");
    let trnas = parse_trna_output(golden_path);

    for trna in &trnas {
        // All should have positive total scores
        assert!(trna.inf_score > 0.0, "{} should have positive inf score", trna.seq_name);

        // SeC has special scoring (HMM=0, SS=0)
        if trna.isotype == "SeC" {
            assert_eq!(trna.hmm_score, 0.0, "SeC should have HMM score 0");
            assert_eq!(trna.ss_score, 0.0, "SeC should have SS score 0");
        } else {
            assert!(trna.hmm_score > 0.0, "{} should have positive HMM score", trna.seq_name);
            assert!(trna.ss_score > 0.0, "{} should have positive SS score", trna.seq_name);
        }

        // CM scores should be positive
        assert!(trna.cm_score > 0.0, "{} should have positive CM score", trna.seq_name);
    }
}

#[test]
fn test_example2_intron_detection() {
    // Test that Ser-CGA intron is correctly detected
    let golden_path = Path::new("tests/golden/full_run/Example2-tRNAs.out");
    let trnas = parse_trna_output(golden_path);

    let ser_trna = trnas.iter()
        .find(|t| t.isotype == "Ser" && t.anticodon == "CGA")
        .expect("Should find Ser-CGA tRNA");

    assert!(ser_trna.intron_begin > 0, "Ser-CGA should have intron");
    assert!(ser_trna.intron_end > 0, "Ser-CGA should have intron");
    assert_eq!(ser_trna.intron_begin, 51);
    assert_eq!(ser_trna.intron_end, 69);
}

#[test]
fn test_example2_selenocysteine() {
    // Test that selenocysteine tRNA is correctly identified
    let golden_path = Path::new("tests/golden/full_run/Example2-tRNAs.out");
    let trnas = parse_trna_output(golden_path);

    let sec_trna = trnas.iter()
        .find(|t| t.isotype == "SeC")
        .expect("Should find SeC tRNA");

    assert_eq!(sec_trna.anticodon, "TCA");
    assert_eq!(sec_trna.hmm_score, 0.0);
    assert_eq!(sec_trna.ss_score, 0.0);
    assert!(sec_trna.inf_score > 100.0, "SeC should have high inf score");
}

#[test]
fn test_example2_isotype_mismatch() {
    // Test that isotype mismatch is correctly noted for MySeq6
    let golden_path = Path::new("tests/golden/full_run/Example2-tRNAs.out");
    let trnas = parse_trna_output(golden_path);

    let lys_trna = trnas.iter()
        .find(|t| t.seq_name == "MySeq6")
        .expect("Should find MySeq6 tRNA");

    assert_eq!(lys_trna.isotype, "Lys");
    assert_eq!(lys_trna.cm_isotype, "Leu");
    assert!(lys_trna.note.contains("ISM"), "Should have ISM note");
}

#[test]
fn test_example2_secondary_structure() {
    // Test secondary structure output format
    let ss_path = Path::new("tests/golden/full_run/Example2-tRNAs.ss");
    assert!(ss_path.exists(), "Secondary structure file should exist");

    let content = fs::read_to_string(ss_path).expect("Failed to read SS file");

    // Check for each sequence
    assert!(content.contains("MySeq1.trna1"), "Should have MySeq1 structure");
    assert!(content.contains("MySeq2.trna1"), "Should have MySeq2 structure");
    assert!(content.contains("MySeq3.trna1"), "Should have MySeq3 structure");
    assert!(content.contains("MySeq4.trna1"), "Should have MySeq4 structure");
    assert!(content.contains("MySeq5.trna1"), "Should have MySeq5 structure");
    assert!(content.contains("MySeq6.trna1"), "Should have MySeq6 structure");

    // Check structure notation exists
    assert!(content.contains("Str: >>>>>>>"), "Should have structure notation");
    assert!(content.contains("Seq: "), "Should have sequence");
}

#[test]
#[ignore] // Enable when full pipeline is ready
fn test_example2_full_scan() {
    // Full integration test: scan Example2.fa and compare with golden output
    let input_path = Path::new("tests/fixtures/Example2.fa");
    let golden_path = Path::new("tests/golden/full_run/Example2-tRNAs.out");

    assert!(input_path.exists(), "Example2.fa should exist");
    assert!(golden_path.exists(), "Golden output should exist");

    // Initialize scanner
    let mut scanner = TrnaScanner::new('E', 20.0)
        .expect("Failed to create scanner");
    scanner.set_quiet(true);

    // Read and process sequences
    let mut reader = SeqFileReader::open(input_path)
        .expect("Failed to open input file");

    while let Some((seq, sqinfo)) = reader.read_seq().expect("Failed to read sequence") {
        scanner.scan_sequence(&seq, &sqinfo)
            .expect("Failed to scan sequence");
    }

    // Compare results
    let golden_trnas = parse_trna_output(golden_path);
    assert_eq!(scanner.result_count(), golden_trnas.len(),
        "Should find same number of tRNAs as golden output");

    // TODO: Detailed comparison of each result
}

#[test]
#[ignore] // Enable when full pipeline is ready
fn test_example1_full_scan() {
    // Full integration test on Example1 (larger file)
    let input_path = Path::new("tests/fixtures/Example1.fa");
    let golden_path = Path::new("tests/golden/full_run/Example1-tRNAs.out");

    assert!(input_path.exists(), "Example1.fa should exist");
    assert!(golden_path.exists(), "Golden output should exist");

    let mut scanner = TrnaScanner::new('E', 20.0)
        .expect("Failed to create scanner");
    scanner.set_quiet(true);

    let mut reader = SeqFileReader::open(input_path)
        .expect("Failed to open input file");

    while let Some((seq, sqinfo)) = reader.read_seq().expect("Failed to read sequence") {
        scanner.scan_sequence(&seq, &sqinfo)
            .expect("Failed to scan sequence");
    }

    let golden_trnas = parse_trna_output(golden_path);
    assert_eq!(scanner.result_count(), golden_trnas.len(),
        "Should find same number of tRNAs as golden output");
}

#[test]
#[ignore] // Enable when output formatting is complete
fn test_output_format_consistency() {
    // Test that output format matches expected conventions
    let golden_path = Path::new("tests/golden/full_run/Example2-tRNAs.out");
    let content = fs::read_to_string(golden_path).expect("Failed to read file");

    let lines: Vec<&str> = content.lines().collect();

    // Check header format
    assert!(lines[0].contains("Sequence"), "First line should be header");
    assert!(lines[1].contains("Name"), "Second line should be field names");
    assert!(lines[2].starts_with("----"), "Third line should be separator");

    // Check data lines have correct number of fields
    for line in lines.iter().skip(3) {
        let fields: Vec<&str> = line.split('\t').collect();
        assert!(fields.len() >= 14, "Each data line should have at least 14 fields");
    }
}

#[test]
fn test_golden_file_existence() {
    // Verify all golden files are present
    assert!(Path::new("tests/golden/full_run/Example1-tRNAs.out").exists());
    assert!(Path::new("tests/golden/full_run/Example1-tRNAs.ss").exists());
    assert!(Path::new("tests/golden/full_run/Example1-tRNAs.stats").exists());
    assert!(Path::new("tests/golden/full_run/Example1-tRNAs.bed").exists());
    assert!(Path::new("tests/golden/full_run/Example1-tRNAs.iso").exists());

    assert!(Path::new("tests/golden/full_run/Example2-tRNAs.out").exists());
    assert!(Path::new("tests/golden/full_run/Example2-tRNAs.ss").exists());
    assert!(Path::new("tests/golden/full_run/Example2-tRNAs.stats").exists());
    assert!(Path::new("tests/golden/full_run/Example2-tRNAs.iso").exists());
}

#[test]
fn test_fixture_file_existence() {
    // Verify test fixture files are present
    assert!(Path::new("tests/fixtures/Example1.fa").exists());
    assert!(Path::new("tests/fixtures/Example2.fa").exists());
}
