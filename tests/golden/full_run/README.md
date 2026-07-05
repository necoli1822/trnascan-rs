# tRNAscan-SE Phase 9 Golden Files - Full Pipeline Test Cases

## Overview

This directory contains golden reference files for validating the complete tRNAscan-SE pipeline in the Rust reimplementation. These files represent the expected output from running tRNAscan-SE on various input sequences.

## Directory Structure

```
tests/golden/
├── full_run/           # Complete pipeline outputs
│   ├── Example1-*      # Eukaryotic test case outputs
│   └── Example2-*      # Additional eukaryotic test case outputs
└── pipeline/           # Pipeline validation test cases
    ├── input_validation.txt    # Input format handling tests
    ├── cmdline_args.txt        # Command-line argument tests
    └── error_cases.txt         # Error handling scenarios
```

## Test Cases

### Test Case 1: Example1 (Eukaryotic Mode)

**Input File:** `tests/fixtures/Example1.fa`
- **Source:** Single eukaryotic sequence (CELF22B7)
- **Length:** 40,222 bp
- **Expected tRNAs:** 5

**Command Used:**
```bash
tRNAscan-SE -E -H \
  -o Example1-tRNAs.out \
  -f Example1-tRNAs.ss \
  -m Example1-tRNAs.stats \
  -b Example1-tRNAs.bed \
  -s Example1-tRNAs.iso \
  Example1.fa
```

**Output Files:**
- `Example1-tRNAs.out` - Tabular results with all tRNA predictions
- `Example1-tRNAs.ss` - Secondary structure predictions
- `Example1-tRNAs.stats` - Search statistics and summary
- `Example1-tRNAs.bed` - BED format output
- `Example1-tRNAs.iso` - Isotype-specific model scores

**Critical Fields in Example1-tRNAs.out:**
| Field | Description | Example |
|-------|-------------|---------|
| Sequence Name | Input sequence identifier | CELF22B7 |
| tRNA # | Sequential tRNA number | 1, 2, 3... |
| Begin | Start position (1-indexed) | 12619 |
| End | End position (1-indexed) | 12738 |
| Type | tRNA isotype | Leu, Ser, Phe, Pro |
| Codon | Anticodon sequence | CAA, AGA, GAA, CGG |
| Intron Begin | Intron start (0 if none) | 12657 (tRNA1), 0 (others) |
| Intron End | Intron end (0 if none) | 12692 (tRNA1), 0 (others) |
| Inf Score | Infernal score (bits) | 74.2, 81.6, 82.5, 71.5 |
| HMM Score | HMM component score | 51.20, 47.50, 56.60, 48.20 |
| 2'Str Score | Secondary structure score | 23.00, 34.10, 25.90, 23.30 |
| Hit Origin | Detection method | Inf (Infernal) |
| Isotype CM | Best isotype model | Leu, Ser, Phe, Pro |
| Isotype Score | Isotype model score | 119.9, 125.0, 112.1, 113.0 |

**tRNA Predictions:**
1. **tRNA1** - Leu-CAA (12619-12738, 120 bp)
   - **Has intron:** 39-74 (12657-12692)
   - Inf Score: 74.2, HMM: 51.20, 2'Str: 23.00
   - Isotype: Leu (119.9)

2. **tRNA2** - Ser-AGA (19480-19561, 82 bp)
   - No intron
   - Inf Score: 81.6, HMM: 47.50, 2'Str: 34.10
   - Isotype: Ser (125.0)

3. **tRNA3** - Phe-GAA (26367-26439, 73 bp)
   - No intron
   - Inf Score: 82.5, HMM: 56.60, 2'Str: 25.90
   - Isotype: Phe (112.1)

4. **tRNA4** - Phe-GAA (26992-26920, 73 bp)
   - **Reverse strand** (end < begin)
   - No intron
   - Inf Score: 82.5, HMM: 56.60, 2'Str: 25.90
   - Isotype: Phe (112.1)
   - **Note:** Identical sequence to tRNA3 (reverse complement)

5. **tRNA5** - Pro-CGG (23765-23694, 72 bp)
   - **Reverse strand**
   - No intron
   - Inf Score: 71.5, HMM: 48.20, 2'Str: 23.30
   - Isotype: Pro (113.0)

---

### Test Case 2: Example2 (Eukaryotic Mode)

**Input File:** `tests/fixtures/Example2.fa`
- **Source:** Multiple synthetic eukaryotic sequences (MySeq1-6)
- **Expected tRNAs:** 6

**Command Used:**
```bash
tRNAscan-SE -E -H \
  -o Example2-tRNAs.out \
  -f Example2-tRNAs.ss \
  -m Example2-tRNAs.stats \
  -s Example2-tRNAs.iso \
  Example2.fa
```

**Output Files:**
- `Example2-tRNAs.out` - Tabular results
- `Example2-tRNAs.ss` - Secondary structures
- `Example2-tRNAs.stats` - Statistics
- `Example2-tRNAs.iso` - Isotype scores

**tRNA Predictions:**
1. **MySeq1.tRNA1** - Thr-TGT (13-85, 73 bp)
   - Inf: 78.0, HMM: 54.80, 2'Str: 23.20
   - Isotype: Thr (93.1)

2. **MySeq2.tRNA1** - Arg-TCT (6-79, 74 bp)
   - Inf: 75.1, HMM: 56.60, 2'Str: 18.50
   - Isotype: Arg (89.3)

3. **MySeq3.tRNA1** - Ser-CGA (14-114, 101 bp)
   - **Has intron:** 51-69
   - Inf: 71.8, HMM: 49.10, 2'Str: 22.70
   - Isotype: Ser (118.3)

4. **MySeq4.tRNA1** - Leu-AAG (6-88, 83 bp)
   - Inf: 65.0, HMM: 43.90, 2'Str: 21.10
   - Isotype: Leu (92.2)

5. **MySeq5.tRNA1** - SeC-TCA (3-89, 87 bp)
   - **Selenocysteine tRNA**
   - Inf: 146.9, HMM: 0.00, 2'Str: 0.00
   - Isotype: SeC (146.9)
   - **Note:** Special case with unique scoring

6. **MySeq6.tRNA1** - Lys-CTT (7-92, 86 bp)
   - Inf: 72.1, HMM: 40.60, 2'Str: 31.50
   - **Isotype mismatch:** Predicted as Lys but best model is Leu (75.7)
   - **Note:** ISM (-73.90) indicates isotype mismatch

---

## Output File Formats

### 1. Tabular Output (.out)

Standard tab-delimited format with the following columns:

```
Sequence    tRNA    Bounds        tRNA  Anti   Intron Bounds   Inf    HMM    2'Str  Hit     Isotype  Isotype
Name        #       Begin  End    Type  Codon  Begin  End      Score  Score  Score  Origin  CM       Score    Note
```

**Key Points:**
- Header spans 3 lines
- Positions are 1-indexed
- Reverse strand: End < Begin
- Intron bounds: 0,0 if no intron
- Hit Origin: Inf (Infernal), Ts (tRNAscan), Eu (EufindtRNA), Bo (Both)

### 2. Secondary Structure (.ss)

Contains:
- tRNA name and coordinates
- Length in base pairs
- Type and anticodon with positions
- Score breakdown (Inf, HMM, 2'Str)
- Intron annotation (if present)
- Sequence and structure notation

**Structure Notation:**
- `>` - 5' side of helix
- `<` - 3' side of helix
- `.` - unpaired base
- Anticodon bases shown in lowercase

**Example:**
```
CELF22B7.trna1 (12619-12738)	Length: 120 bp
Type: Leu	Anticodon: CAA at 35-37 (12653-12655)	Score: 74.2
Possible intron: 39-74 (12657-12692)
HMM Sc=51.20	Sec struct Sc=23.00
Seq: GCACGGATGGCCGAGTGGTctAAGGCGCCAGACTCAAGcgaaatgcttgcctcatgctcgaggtcgactgggtgTTCTGGTACTCGTATGGGTGCGTGGGTTCGAATCCCACTTCGTGCA
Str: >>>>>>>..>>>...........<<<.>>>>>...........................................<<<<<.>>>>....<<<<..>>>>>.......<<<<<<<<<<<<.
```

### 3. Statistics (.stats)

Contains:
- Command-line parameters
- Search mode and models used
- First-pass statistics
- Infernal statistics
- Overall summary counts
- Isotype/anticodon distribution

**Key Metrics:**
- Sequences read
- Bases read (×2 for both strands)
- tRNAs predicted
- Average tRNA length
- Scan speed (bp/sec)
- tRNAs with introns
- Isotype counts

### 4. BED Format (.bed)

BED12 format for genome browsers:

```
chrom  chromStart  chromEnd  name              score  strand  thickStart  thickEnd  itemRgb  blockCount  blockSizes  blockStarts
```

**Example:**
```
CELF22B7  12618  12738  CELF22B7.tRNA1-LeuCAA  742  +  12618  12738  0  2  38,46,  0,74,
```

**Notes:**
- 0-indexed positions (chromStart is begin-1)
- Score is Inf score × 10 (integer)
- Introns represented as multi-block features
- No intron: blockCount=1
- With intron: blockCount=2

### 5. Isotype Scores (.iso)

Tab-delimited matrix of scores for all isotype models:

```
tRNAscanID  Anticodon_predicted_isotype  Ala  Arg  Asn  ...  Val  iMet
```

**Notes:**
- One row per tRNA
- Columns for each amino acid isotype
- -999 indicates model not tested
- Highest score usually matches predicted type
- Useful for identifying mismatches

---

## Pipeline Validation Test Cases

### Input Validation (`pipeline/input_validation.txt`)

Test cases covering:
- Valid FASTA formats (single/multiple sequences)
- Empty files
- Missing headers
- Whitespace handling
- Invalid characters
- Long sequence names
- Special characters
- Case sensitivity
- Multi-line sequences

### Command-Line Arguments (`pipeline/cmdline_args.txt`)

Test cases covering:
- Search mode options (-E, -B, -A, -M, -O, -G)
- Output file options (-o, -f, -b, -m, -s)
- Score cutoffs (-X)
- Breakdown mode (-H)
- Pseudogene detection (-D)
- Output formats (--detail, --brief)
- Default file generation (?)
- Prefix option (-p)
- Padding (-z)
- Force overwrite (-Q)
- Invalid/conflicting options

### Error Cases (`pipeline/error_cases.txt`)

Test cases covering:
- Empty/corrupted input files
- Permission errors
- Disk space issues
- Invalid parameter values
- Missing model files
- Memory allocation failures
- Signal interrupts
- Configuration errors
- Numerical overflow

---

## Critical Output Fields for Validation

When implementing the Rust version, ensure these fields are preserved exactly:

### Tabular Output (.out)
| Field | Type | Notes |
|-------|------|-------|
| Sequence Name | String | Must match input exactly |
| tRNA # | Integer | Sequential per sequence |
| Begin | Integer | 1-indexed, may be > End for reverse |
| End | Integer | 1-indexed |
| Type | String | 3-letter amino acid code or SeC |
| Codon | String | 3-letter anticodon |
| Intron Begin | Integer | 0 if none, relative to sequence start |
| Intron End | Integer | 0 if none |
| Inf Score | Float | Format: %.1f (one decimal) |
| HMM Score | Float | Format: %.2f (two decimals) |
| 2'Str Score | Float | Format: %.2f (two decimals) |
| Hit Origin | String | Inf, Ts, Eu, Bo |
| Isotype CM | String | Best matching isotype model |
| Isotype Score | Float | Format: %.1f |
| Note | String | Optional (ISM, pseudogene, etc.) |

### Secondary Structure (.ss)
- Coordinate format: `(begin-end)` with spaces
- Anticodon position: `at X-Y (genomic_X-genomic_Y)`
- Intron annotation: `Possible intron: X-Y (genomic_X-genomic_Y)`
- Score format: `Inf Sc=%.1f  HMM Sc=%.2f  Sec struct Sc=%.2f`
- Structure notation must align with sequence (same length)

### Statistics (.stats)
- Timestamps format: `Day Mon DD HH:MM:SS TZ YYYY`
- Counts must be accurate and consistent
- Speed calculations: bases/time
- Isotype table format must be preserved

### BED Format (.bed)
- 0-indexed positions (adjust from 1-indexed output)
- Score = Inf score × 10, rounded to integer
- Name format: `seqname.tRNAX-TypeCodon`
- Block representation for introns

### Isotype Scores (.iso)
- Tab-delimited matrix
- -999 for non-tested models
- Floating-point precision consistent

---

## Version Information

**Original tRNAscan-SE Version:** 2.0rc1 (April 2017)
- Golden files generated on: Mon Apr 10 01:02:12 PDT 2017
- Host: pismo.soe.ucsc.edu
- User: pchan

**Models Used:**
- Eukaryotic: TRNAinf-euk.cm
- Selenocysteine: TRNAinf-euk-SeC.cm
- Isotype-specific: TRNAinf-euk-iso

**Search Parameters:**
- Mode: Infernal First Pass → Infernal
- First pass cutoff: 10 bits
- Final cutoff: 20 bits (default)
- Isotype-specific scan: Yes
- Pseudogene checking: Yes
- Padding: 8 bp (default)

---

## Usage in Rust Implementation

### Running Validation Tests

```bash
# Test against Example1
cargo test --test integration_test -- example1_full_pipeline

# Test against Example2
cargo test --test integration_test -- example2_full_pipeline

# Test all output formats
cargo test --test integration_test -- output_formats

# Test error handling
cargo test --test integration_test -- error_cases
```

### Comparing Outputs

```bash
# Compare tabular output
diff -u tests/golden/full_run/Example1-tRNAs.out output/Example1-tRNAs.out

# Compare secondary structure
diff -u tests/golden/full_run/Example1-tRNAs.ss output/Example1-tRNAs.ss

# Compare statistics (note: timestamps and hosts will differ)
diff -I "Started:" -I "ended:" -I "host" \
  tests/golden/full_run/Example1-tRNAs.stats output/Example1-tRNAs.stats
```

### Field-by-Field Validation

When testing, validate each field independently:

```rust
#[test]
fn test_trna_coordinates() {
    let results = parse_output("Example1-tRNAs.out");
    assert_eq!(results[0].begin, 12619);
    assert_eq!(results[0].end, 12738);
    assert_eq!(results[0].strand, Strand::Forward);
}

#[test]
fn test_intron_detection() {
    let results = parse_output("Example1-tRNAs.out");
    assert_eq!(results[0].intron_begin, 12657);
    assert_eq!(results[0].intron_end, 12692);
}

#[test]
fn test_score_precision() {
    let results = parse_output("Example1-tRNAs.out");
    assert_eq!(format!("{:.1}", results[0].inf_score), "74.2");
    assert_eq!(format!("{:.2}", results[0].hmm_score), "51.20");
}
```

---

## Known Issues and Special Cases

### 1. Reverse Strand Representation
- Coordinates: End < Begin (e.g., 23765-23694)
- Anticodon coordinates are also reversed
- Intron coordinates follow same convention

### 2. Isotype Mismatches
- Example: MySeq6.tRNA1 predicted as Lys but scores higher for Leu
- Indicated by "ISM" note with score difference
- Should still output predicted type, not best-scoring type

### 3. Selenocysteine tRNAs
- Type: SeC, Anticodon: TCA
- May have 0.00 for HMM and 2'Str scores
- Uses specialized model (TRNAinf-euk-SeC.cm)

### 4. Floating-Point Precision
- Inf Score: 1 decimal place
- HMM Score: 2 decimal places
- 2'Str Score: 2 decimal places
- Isotype Score: 1 decimal place

### 5. BED Coordinate Conversion
- BED uses 0-indexed, half-open intervals
- chromStart = Begin - 1
- chromEnd = End (unchanged)
- For reverse strand: chromStart = End - 1, chromEnd = Begin

---

## References

- tRNAscan-SE 2.0 Paper: Chan & Lowe (2019) *Nucleic Acids Research*
- Infernal: Nawrocki & Eddy (2013) *Bioinformatics*
- BED Format: UCSC Genome Browser format specification
- tRNA Structure: Sprinzl Database conventions

---

## Maintenance

When adding new test cases:

1. Document the exact command used
2. Include input file characteristics
3. List expected outputs and counts
4. Note any special features (introns, pseudogenes, etc.)
5. Verify all output files are consistent
6. Update this README with new test case details

---

## Contact

For questions about golden files or test cases:
- GitHub Issues: https://github.com/bactars/tRNAscan-SE
- Original tRNAscan-SE: http://lowelab.ucsc.edu/tRNAscan-SE/
