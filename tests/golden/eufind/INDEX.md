# EuFindtRNA Golden Files - Index

## Quick Start

```bash
# Regenerate all golden files
make clean && make

# Verify golden files
./verify.sh
```

## File Organization

### Documentation
- **README.md** - Detailed explanation of test methodology and results
- **SUMMARY.md** - Comprehensive summary with all key values and parameters
- **INDEX.md** - This file (navigation and quick reference)

### Golden Data Files
- **bbox_scores.txt** - B-box detection scores at each position
- **abox_scores.txt** - A-box detection results and AB distance scores
- **trna_detection.txt** - Complete tRNA structure after full detection
- **weight_matrices.txt** - All scoring matrices and lookup tables

### Source Code
- **gen_eufind.c** - Main golden file generator (9.9 KB)
- **gen_weight_matrices.c** - Matrix extractor (3.5 KB)

### Build Tools
- **Makefile** - Build system for regenerating golden files
- **verify.sh** - Verification script for data integrity

### Compiled Binaries
- **gen_eufind** - Executable for main tests
- **gen_weight_matrices** - Executable for matrix extraction

## Test Sequence

**Type:** tRNA-Phe from *Caenorhabditis elegans*
**Source:** Example1.fa (lines 254-256)
**Length:** 139 bp

```
TCGCTGTTAGTTACCATCGCACGGATGGCCGAGTGGTCTAAGGCGCCAGA
CTCAAGCGAAATGCTTGCCTCATGCTCGAGGTCGACTGGGTGTTCTGGTA
CTCGTATGGGTGCGTGGGTTCGAATCCCACTTCGTGCAG
```

## Key Test Results (Quick Reference)

| Metric | Value |
|--------|-------|
| B-box position | 116 |
| B-box score | -1.9170 |
| B-box range | 116-126 |
| A-box position | 24-43 |
| A-box score | -13.7640 |
| AB distance | 72 bp |
| AB dist score | -5.4420 |
| Term position | -1 (not found) |
| Term score | -0.5500 (penalty) |
| **Total score** | **-21.6730** |
| Score threshold | -31.8 |
| **Detection result** | **PASS** |
| tRNA start | 19 |
| tRNA end | 138 |
| tRNA length | 120 bp |

## Function Coverage

### Tested Functions
- [x] `IntEncodeSeq()` - DNA to integer encoding
- [x] `GetBbox()` - B-box scanning
- [x] `GetBestABox()` - A-box detection with gaps
- [x] `Get_ABdist_weight()` - AB distance scoring
- [x] `GetBestTrxTerm()` - Termination signal detection
- [x] `Init_tRNA()` - Structure initialization

### Scoring Matrices Covered
- [x] Abox_Mat[6][21] - A-box position weight matrix
- [x] Bbox_Mat[6][11] - B-box position weight matrix
- [x] ABDistIdx_Mat[7] + ABDistSc_Mat[7] - AB distance lookup
- [x] BTermDistIdx_Mat[9] + BTermDistSc_Mat[9] - B-Term distance lookup

## Usage in Rust Testing

### Unit Tests
1. Test integer encoding with known sequence
2. Test B-box score calculation at position 116
3. Test A-box score calculation at position 24
4. Test AB distance score lookup (72 bp → -5.4420)
5. Test termination signal detection (missing → -0.5500)

### Integration Test
Run full detection on test sequence and verify:
- All intermediate scores match golden values
- Final TRNA_TYPE structure matches exactly
- Detection passes threshold check

### Floating-Point Tolerance
Use `assert_approx_eq!` with epsilon = 0.0001 for score comparisons.

## Regeneration

To regenerate golden files after code changes:

```bash
cd tests/golden/eufind
make clean
make
./verify.sh
```

## Dependencies

### C Source Files Required
- `pavesi.c` / `pavesi.h` - Main EuFindtRNA implementation
- `eufind_const.h` - Constants and structures
- `revcomp.c` - Reverse complement function
- `iupac.c` - IUPAC code tables
- `sqerror.c` - Error handling
- `sre_ctype.c` - Character utilities
- `sre_string.c` - String utilities
- `sre_math.c` - Math utilities

### Object Files Linked
Located in `../../../original/src/`:
- pavesi.o
- revcomp.o
- iupac.o
- sqerror.o
- sre_ctype.o
- sre_string.o
- sre_math.o

## Version Information

**Created:** March 4, 2026
**tRNAscan-SE Version:** 2.0 branch
**Algorithm:** Pavesi et al. (1994)

## References

Pavesi, A., Conterio, F., Bolchi, A., Dieci, G., & Ottonello, S. (1994).
"Identification of new eukaryotic tRNA genes in genomic DNA databases by a
multistep weight matrix analysis of transcriptional control regions."
*Nucleic Acids Research*, 22(7), 1247-1256.
https://doi.org/10.1093/nar/22.7.1247
