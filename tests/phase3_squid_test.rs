// Phase 3: SQUID library tests
// Tests for sequence format detection and SQINFO field parsing

use trnascan_rs::squid::{detect_format, K_EMBL, K_GENBANK, K_PEARSON, K_SELEX};

#[test]
fn test_format_detection() {
    // Test automatic sequence format detection (FASTA, GenBank, EMBL, etc.)

    // Test FASTA format detection
    let format = detect_format("tests/golden/squid/inputs/test.fasta").unwrap();
    assert_eq!(format, K_PEARSON, "Expected FASTA/Pearson format (7)");

    // Test GenBank format detection
    let format = detect_format("tests/golden/squid/inputs/test.gb").unwrap();
    assert_eq!(format, K_GENBANK, "Expected GenBank format (2)");

    // Test EMBL format detection
    let format = detect_format("tests/golden/squid/inputs/test.embl").unwrap();
    assert_eq!(format, K_EMBL, "Expected EMBL format (4)");

    // Test raw format (should default to SELEX)
    let format = detect_format("tests/golden/squid/inputs/test.raw").unwrap();
    assert_eq!(format, K_SELEX, "Expected SELEX format (10) for raw file");
}

#[test]
fn test_sqinfo_fields() {
    use trnascan_rs::squid::{SeqFileReader, SQINFO_ACC, SQINFO_DESC, SQINFO_ID, SQINFO_NAME, SQINFO_OLEN, SQINFO_START, SQINFO_STOP};

    // Test reading FASTA file
    let mut reader = SeqFileReader::open("tests/golden/squid/inputs/test.fasta").unwrap();

    // Read first sequence
    let result = reader.read_seq().unwrap();
    assert!(result.is_some(), "Should read first sequence");
    let (seq, sqinfo) = result.unwrap();

    // Check flags and values
    assert_eq!(sqinfo.name, "seq1");
    assert_eq!(sqinfo.desc, "First test sequence");
    assert_eq!(sqinfo.len, 64);
    assert_eq!(seq.len() as usize, 64);
    assert!(sqinfo.has_flag(SQINFO_NAME));
    assert!(sqinfo.has_flag(SQINFO_DESC));

    // Test reading GenBank file with more fields
    let mut reader = SeqFileReader::open("tests/golden/squid/inputs/test.gb").unwrap();

    let result = reader.read_seq().unwrap();
    assert!(result.is_some(), "Should read GenBank sequence");
    let (_seq, sqinfo) = result.unwrap();

    // Check all expected fields
    assert_eq!(sqinfo.name, "TESTSEQ1");
    assert_eq!(sqinfo.id, "TESTSEQ1");
    assert_eq!(sqinfo.acc, "TEST001");
    assert!(sqinfo.desc.contains("Test GenBank sequence"));
    assert_eq!(sqinfo.len, 64);
    assert_eq!(sqinfo.start, 1);
    assert_eq!(sqinfo.stop, 64);
    assert_eq!(sqinfo.olen, 64);

    assert!(sqinfo.has_flag(SQINFO_NAME));
    assert!(sqinfo.has_flag(SQINFO_ID));
    assert!(sqinfo.has_flag(SQINFO_ACC));
    assert!(sqinfo.has_flag(SQINFO_DESC));
    assert!(sqinfo.has_flag(SQINFO_START));
    assert!(sqinfo.has_flag(SQINFO_STOP));
    assert!(sqinfo.has_flag(SQINFO_OLEN));
}

#[test]
fn test_sequence_reading() {
    use trnascan_rs::squid::SeqFileReader;

    // Test reading multiple sequences from FASTA
    let mut reader = SeqFileReader::open("tests/golden/squid/inputs/test.fasta").unwrap();

    let mut count = 0;
    while let Some((seq, sqinfo)) = reader.read_seq().unwrap() {
        count += 1;
        assert!(!seq.is_empty(), "Sequence should not be empty");
        assert!(!sqinfo.name.is_empty(), "Name should not be empty");
        assert_eq!(sqinfo.len, seq.len(), "Length should match sequence");
    }
    assert_eq!(count, 3, "Should read 3 sequences from test.fasta");

    // Test reading from GenBank
    let mut reader = SeqFileReader::open("tests/golden/squid/inputs/test.gb").unwrap();

    let mut count = 0;
    while let Some((seq, sqinfo)) = reader.read_seq().unwrap() {
        count += 1;
        // GenBank sequences should be lowercase in the test file
        assert!(!seq.is_empty(), "Sequence should not be empty");
        assert!(!sqinfo.name.is_empty(), "Name should not be empty");
        assert!(!sqinfo.acc.is_empty(), "Accession should not be empty");
    }
    assert_eq!(count, 2, "Should read 2 sequences from test.gb");

    // Test reading from EMBL
    let mut reader = SeqFileReader::open("tests/golden/squid/inputs/test.embl").unwrap();

    let mut count = 0;
    while let Some((_seq, sqinfo)) = reader.read_seq().unwrap() {
        count += 1;
        assert!(!sqinfo.name.is_empty(), "Name should not be empty");
        assert!(!sqinfo.id.is_empty(), "ID should not be empty");
    }
    assert_eq!(count, 2, "Should read 2 sequences from test.embl");
}

#[test]
fn test_sqinfo_initialization() {
    use trnascan_rs::squid::{SqInfo, SQINFO_ACC, SQINFO_DESC, SQINFO_ID, SQINFO_NAME, SQINFO_OLEN, SQINFO_START, SQINFO_STOP};

    // Create empty SQINFO
    let mut sqinfo = SqInfo::new();
    assert_eq!(sqinfo.flags, 0, "Flags should be 0 initially");

    // Set various fields
    sqinfo.set_string_field("TestSeq", SQINFO_NAME);
    assert_eq!(sqinfo.name, "TestSeq");
    assert!(sqinfo.has_flag(SQINFO_NAME));

    sqinfo.set_string_field("TEST001", SQINFO_ID);
    assert_eq!(sqinfo.id, "TEST001");
    assert!(sqinfo.has_flag(SQINFO_ID));

    sqinfo.set_string_field("ACC123", SQINFO_ACC);
    assert_eq!(sqinfo.acc, "ACC123");
    assert!(sqinfo.has_flag(SQINFO_ACC));

    sqinfo.set_string_field("Test description", SQINFO_DESC);
    assert_eq!(sqinfo.desc, "Test description");
    assert!(sqinfo.has_flag(SQINFO_DESC));

    // Set integer fields
    sqinfo.set_int_field(1, SQINFO_START);
    assert_eq!(sqinfo.start, 1);
    assert!(sqinfo.has_flag(SQINFO_START));

    sqinfo.set_int_field(100, SQINFO_STOP);
    assert_eq!(sqinfo.stop, 100);
    assert!(sqinfo.has_flag(SQINFO_STOP));

    sqinfo.set_int_field(100, SQINFO_OLEN);
    assert_eq!(sqinfo.olen, 100);
    assert!(sqinfo.has_flag(SQINFO_OLEN));

    // Test that setting to 0 doesn't set flag
    let mut sqinfo2 = SqInfo::new();
    sqinfo2.set_int_field(0, SQINFO_START);
    assert!(!sqinfo2.has_flag(SQINFO_START), "Zero value shouldn't set flag");
}

#[test]
fn test_iupac_functions() {
    use trnascan_rs::squid::{is_gap, is_seq_char, seq_type, to_dna, to_rna, K_DNA, K_RNA};

    // Test character validation
    assert!(is_seq_char(b'A'));
    assert!(is_seq_char(b'G'));
    assert!(is_seq_char(b'-'));
    assert!(!is_seq_char(b' '));

    // Test gap detection
    assert!(is_gap(b'-'));
    assert!(is_gap(b'.'));
    assert!(!is_gap(b'A'));

    // Test sequence type detection
    let dna_seq = b"ACGTACGT";
    assert_eq!(seq_type(dna_seq), K_DNA);

    let rna_seq = b"ACGUACGU";
    assert_eq!(seq_type(rna_seq), K_RNA);

    // Test RNA/DNA conversion
    let mut seq = b"ACGTACGT".to_vec();
    to_rna(&mut seq);
    assert_eq!(seq, b"ACGUACGU");

    let mut seq = b"ACGUACGU".to_vec();
    to_dna(&mut seq);
    assert_eq!(seq, b"ACGTACGT");
}
