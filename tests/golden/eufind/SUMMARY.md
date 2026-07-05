# EuFindtRNA Golden Files - Summary

## Overview

This directory contains golden reference files for Phase 6 (EuFindtRNA) of the tRNAscan-SE pipeline. These files document the expected behavior of the Pavesi algorithm for detecting eukaryotic tRNA genes via transcriptional control regions.

## Generated Files

### Core Golden Files

1. **bbox_scores.txt** (535 bytes)
   - B-box detection scores
   - Documents scanning behavior and score calculation
   - Shows position 116 detected with score -1.9170

2. **abox_scores.txt** (540 bytes)
   - A-box detection results
   - Shows best A-box at position 24-43
   - Documents AB distance scoring (72 bp → -5.4420)

3. **trna_detection.txt** (1017 bytes)
   - Complete tRNA detection workflow
   - Full TRNA_TYPE structure output
   - Shows all scoring components and final boundaries

4. **weight_matrices.txt** (5.3 KB)
   - Complete Abox_Mat[6][21] and Bbox_Mat[6][11] matrices
   - AB distance lookup table (7 entries)
   - B-to-Term distance lookup table (9 entries)

## Test Programs

### gen_eufind.c (9.8 KB)
Main test program that exercises:
- `IntEncodeSeq()` - DNA to integer encoding
- `GetBbox()` - B-box scanning with position weight matrix
- `GetBestABox()` - A-box detection with gap handling
- `Get_ABdist_weight()` - AB distance scoring
- `GetBestTrxTerm()` - Termination signal detection

**Test sequence:** tRNA-Phe from *C. elegans* (139 bp)
```
TCGCTGTTAGTTACCATCGCACGGATGGCCGAGTGGTCTAAGGCGCCAGA
CTCAAGCGAAATGCTTGCCTCATGCTCGAGGTCGACTGGGTGTTCTGGTA
CTCGTATGGGTGCGTGGGTTCGAATCCCACTTCGTGCAG
```

### gen_weight_matrices.c (3.5 KB)
Extracts scoring matrices from pavesi.c source code.

## Build System

**Makefile** (595 bytes)
```bash
make          # Build and generate all golden files
make clean    # Remove generated files
```

## Key Scoring Parameters

| Parameter | Value | Meaning |
|-----------|-------|---------|
| BBOX_CUTOFF | -14.14 | Minimum B-box score to accept |
| BBOX_START_IDX | 45 | Start scanning from position 45 |
| ABOX_LEN | 21 | A-box length (with gaps) |
| BBOX_LEN | 11 | B-box length |
| MIN_AB_BOX_DIST | 24 | Minimum A-to-B distance |
| AB_BOX_DIST_RANGE | 116 | Maximum search range |
| MAX_TERM_SEARCH | 133 | Maximum B-to-Term distance |
| TOT_SCORE_THRESH | -31.8 | Total score threshold |

## Test Results

### B-box Detection
- Position: 116
- Score: -1.9170
- Range: 116-126

### A-box Detection
- Position: 24-43
- A-box score: -13.7640
- AB distance: 72 bp
- AB score: -5.4420
- Combined internal: -21.1230

### Termination Signal
- Not found (returns -1)
- Penalty score: -0.5500

### Total Score
- **Final: -21.6730**
- Threshold: -31.8
- **Result: PASS** (above threshold)

### tRNA Boundaries
- Start: 19 (5 bp upstream of A-box)
- End: 138 (sequence end, no terminator)
- Length: 120 bp

## Matrix Structure

### A-box Matrix (Abox_Mat)
- Dimensions: 6 rows × 21 positions
- Rows:
  - 0: A scores (log-odds)
  - 1: C scores
  - 2: G scores
  - 3: T scores
  - 4: Gap penalties
  - 5: Ambiguous base (best of A/C/G/T)
- Special positions:
  - Position 10: Gap at tRNA position 17
  - Positions 14-15: Gap at tRNA positions 20a/20b (4 configurations tested)

### B-box Matrix (Bbox_Mat)
- Dimensions: 6 rows × 11 positions
- Same row structure as A-box
- No gap positions

### Distance Scoring Tables

**AB Distance (ABDistIdx_Mat/ABDistSc_Mat):**
| Distance ≤ | Score |
|------------|-------|
| 30 | -0.46 |
| 36 | -1.83 |
| 42 | -2.35 |
| 48 | -3.24 |
| 54 | -4.06 |
| 60 | -3.83 |
| 66 | -4.75 |

**B-Term Distance (BTermDistIdx_Mat/BTermDistSc_Mat):**
| Distance ≤ | Score |
|------------|-------|
| 17 | -0.54 |
| 23 | -1.40 |
| 29 | -2.80 |
| 35 | -3.36 |
| 41 | -3.24 |
| 47 | -5.44 |
| 53 | -5.44 |
| 59 | -4.06 |
| 100 | -5.44 |

## Usage for Rust Implementation

These files provide comprehensive test data for validating:

1. **Sequence Encoding**
   - Verify DNA → integer conversion (A=0, C=1, G=2, T=3, N=5)

2. **Matrix Scoring**
   - Test weight matrix lookups match exactly
   - Verify accumulation of position scores

3. **B-box Scanning**
   - Check scanning starts at correct position
   - Verify score threshold comparison
   - Test position advancement logic

4. **A-box Detection**
   - Verify gap handling (4 configurations)
   - Test AB distance calculation
   - Check best score selection

5. **Distance Scoring**
   - Verify lookup table interpolation
   - Test boundary conditions

6. **Integration Test**
   - Full tRNA detection should match all intermediate values
   - Final structure should match exactly

## Validation Criteria

- Floating-point scores: tolerance ≤ 0.0001
- Position indices: must match exactly
- Structure fields: must match exactly

## References

Pavesi, A., Conterio, F., Bolchi, A., Dieci, G., & Ottonello, S. (1994).
"Identification of new eukaryotic tRNA genes in genomic DNA databases by a
multistep weight matrix analysis of transcriptional control regions."
*Nucleic Acids Research*, 22(7), 1247-1256.

## Notes

- The test sequence comes from *C. elegans* Example1.fa lines 254-256
- This is a real tRNA-Phe that tRNAscan-SE successfully detects
- The sequence does not contain a TTTT terminator, so termination score is a penalty
- Despite no terminator, the total score (-21.67) is well above the threshold (-31.8)
