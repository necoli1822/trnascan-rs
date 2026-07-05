# EuFindtRNA Golden Files

This directory contains golden test files for Phase 6 (EuFindtRNA) of the tRNAscan-SE pipeline.

## Overview

EuFindtRNA implements the Pavesi algorithm for detecting eukaryotic tRNA genes by finding:
- **A-box**: Internal promoter element (position -31 to -51 relative to transcription start)
- **B-box**: Internal promoter element (position +52 to +62)
- **Termination signal**: Run of 4+ T's (TTTT) downstream of B-box

## Test Sequence

The primary test sequence is a real tRNA-Phe from *C. elegans* (from Example1.fa):

```
TCGCTGTTAGTTACCATCGCACGGATGGCCGAGTGGTCTAAGGCGCCAGA
CTCAAGCGAAATGCTTGCCTCATGCTCGAGGTCGACTGGGTGTTCTGGTA
CTCGTATGGGTGCGTGGGTTCGAATCCCACTTCGTGCAG
```

Length: 139 bp

## Golden Files

### 1. bbox_scores.txt

Tests the `GetBbox()` function which scans for B-box patterns using a position weight matrix.

**Key parameters:**
- B-box length: 11 bp (BBOX_LEN)
- Score cutoff: -14.14 (BBOX_CUTOFF)
- Start scan position: 45 (BBOX_START_IDX)

**Test results:**
- Position 116: score = -1.9170 (FOUND)
- B-box region: positions 116-126

### 2. abox_scores.txt

Tests the `GetBestABox()` function which finds the best-scoring A-box given a B-box position.

**Key parameters:**
- A-box length: 21 bp (ABOX_LEN)
- Min AB distance: 24 bp (MIN_AB_BOX_DIST)
- AB distance range: 116 bp (AB_BOX_DIST_RANGE)
- Gap handling: Searches 4 gap configurations at positions 20a/20b

**Test results:**
- A-box position: 24-43
- A-box score: -13.7640
- AB distance: 72 bp (between A-box end and B-box start)
- AB distance score: -5.4420
- Combined internal score: -21.1230

### 3. trna_detection.txt

Tests the complete tRNA detection workflow:

**Step 1: B-box Detection**
- Position: 116
- Score: -1.9170

**Step 2: A-box Detection**
- Position: 24-43
- A-box score: -13.7640
- AB distance: 72 bp
- AB distance score: -5.4420

**Step 3: Termination Signal**
- Position: -1 (not found in sequence)
- Term score: -0.5500 (penalty for missing terminator)

**Step 4: Total Score**
- Total: -21.6730
- Components:
  - A-box: -13.7640
  - B-box: -1.9170
  - AB distance: -5.4420
  - Termination: -0.5500

**Step 5: tRNA Boundaries**
- Start: 19 (5 bp upstream of A-box)
- End: 138 (end of sequence since no terminator found)

## Score Interpretation

In the Pavesi scoring system:
- **Higher (less negative) scores are better**
- Scores close to 0 indicate perfect matches to the weight matrices
- Scores below the cutoff (-14.14 for B-box) are rejected
- Total score threshold for detection: -31.8 (TOT_SCORE_THRESH)

## Key Functions

1. **IntEncodeSeq()**: Converts DNA sequence to integer encoding (A=0, C=1, G=2, T=3, ambiguous=5)
2. **GetBbox()**: Scans for B-box using weight matrix scoring
3. **GetBestABox()**: Finds best A-box position considering gap variations
4. **Get_ABdist_weight()**: Scores A-box to B-box distance
5. **GetBestTrxTerm()**: Finds TTTT termination signal

## Weight Matrices

The scoring uses two position weight matrices:
- **Abox_Mat[6][21]**: 21 positions with gap at position 17
- **Bbox_Mat[6][11]**: 11 positions

Each matrix has 6 rows:
- Rows 0-3: A, C, G, T scores (log-odds)
- Row 4: Gap penalty
- Row 5: Ambiguous base score (best of A/C/G/T)

## Compilation

```bash
gcc -I../../original/src \
    gen_eufind.c \
    ../../original/src/pavesi.o \
    ../../original/src/revcomp.o \
    ../../original/src/iupac.o \
    ../../original/src/sqerror.o \
    ../../original/src/sre_ctype.o \
    ../../original/src/sre_string.o \
    ../../original/src/sre_math.o \
    -o gen_eufind -lm
```

## References

Pavesi, A., Conterio, F., Bolchi, A., Dieci, G., & Ottonello, S. (1994).
"Identification of new eukaryotic tRNA genes in genomic DNA databases by a
multistep weight matrix analysis of transcriptional control regions."
*Nucleic Acids Research*, 22(7), 1247-1256.

## Usage for Rust Implementation

These golden files should be used to validate:
1. Integer sequence encoding matches C implementation
2. B-box scoring produces identical results
3. A-box scoring with gap handling is correct
4. AB distance scoring follows lookup table
5. Termination signal detection works as expected
6. Overall tRNA structure calculation is accurate

Floating-point comparison should allow for small rounding differences (< 0.0001).
