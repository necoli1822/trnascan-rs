# SQUID Library Golden Files - Phase 3

This directory contains golden test files for the SQUID library sequence I/O functions, specifically testing format detection and SQINFO structure parsing.

## Purpose

These golden files serve as reference outputs for testing the Rust reimplementation of SQUID library functions in Phase 3 of the tRNAscan-SE Rust port.

## Files

### Test Generator
- **gen_squid.c** - C program that generates golden output using original SQUID library
  - Tests `SeqfileFormat()` for format detection
  - Tests `ReadSeq()` for parsing sequences
  - Tests SQINFO structure field extraction

### Input Files (inputs/)
- **test.fasta** - FASTA format with 3 sequences
- **test.gb** - GenBank format with 2 sequences
- **test.embl** - EMBL format with 2 sequences
- **test.raw** - Raw unformatted sequence

### Golden Outputs
- **format_detection.txt** - Results of `SeqfileFormat()` calls
- **sqinfo_fields.txt** - Results of `ReadSeq()` and SQINFO field extraction

## Key Test Coverage

### Format Detection
Tests that `SeqfileFormat()` correctly identifies:
- **kPearson (7)** - FASTA format
- **kGenBank (2)** - GenBank format
- **kEMBL (4)** - EMBL format
- **kSelex (10)** - Selex format (raw file detected as this)

### SQINFO Fields
Tests extraction of SQINFO structure fields with correct flags:

**FASTA Format (flags 0x41 or 0x49):**
- SQINFO_NAME (0x01) - sequence name
- SQINFO_DESC (0x08) - description (if present)
- SQINFO_LEN (0x40) - sequence length

**GenBank/EMBL Format (flags 0x27f):**
- SQINFO_NAME (0x01)
- SQINFO_ID (0x02)
- SQINFO_ACC (0x04)
- SQINFO_DESC (0x08)
- SQINFO_START (0x10)
- SQINFO_STOP (0x20)
- SQINFO_LEN (0x40)
- SQINFO_OLEN (0x200)

## Building and Running

```bash
# Build the test generator
cd original
gcc -I src -o ../tests/golden/squid/gen_squid \
    ../tests/golden/squid/gen_squid.c \
    src/sqio.c src/sre_string.c src/sqerror.c src/iupac.c \
    src/selex.c src/interleaved.c src/msf.c src/sre_ctype.c \
    src/sre_math.c src/alignio.c src/types.c -lm

# Run the generator
cd ..
./tests/golden/squid/gen_squid
```

## Usage in Rust Tests

The Rust implementation should:
1. Parse the same input files
2. Generate equivalent output
3. Compare against these golden files

Example test structure:
```rust
#[test]
fn test_format_detection() {
    let format = detect_format("tests/golden/squid/inputs/test.fasta");
    assert_eq!(format, SequenceFormat::Pearson);
}

#[test]
fn test_sqinfo_parsing() {
    let sqinfo = read_sequence("tests/golden/squid/inputs/test.gb", Format::GenBank);
    assert_eq!(sqinfo.name, "TESTSEQ1");
    assert_eq!(sqinfo.len, 64);
    // ... more assertions
}
```

## Notes

- Flag values are in hexadecimal (e.g., 0x49 = 0x01 | 0x08 | 0x40 = NAME | DESC | LEN)
- The raw file is detected as kSelex (format 10) rather than kRaw (13)
- GenBank and EMBL formats extract more metadata than FASTA
- Sequence content case is preserved from input files
