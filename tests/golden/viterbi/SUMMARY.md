# Phase 5 Viterbi Core Golden Files - Summary

## Overview

Golden files for Phase 5 (Viterbi core algorithm) have been successfully generated. These files provide reference outputs for validating the Rust reimplementation of the tRNAscan-SE Viterbi alignment algorithm.

## Generated Files

### 1. Test Program
- **gen_viterbi.c** (5.5 KB)
  - Loads CM model using ReadCM()
  - Converts to integer model via RearrangeCM()
  - Runs ViterbiAlign() on three test sequences
  - Outputs scores and detailed traceback trees

### 2. Golden Output Files
- **viterbi_scores.txt** (203 bytes)
  - Tab-delimited scores for each test sequence
  - Format: `<name>\t<length>\t<score>`

- **traceback.txt** (21 KB)
  - Complete traceback tree structures
  - Depth-first traversal format
  - Shows node indices, state types, and emission positions

### 3. Documentation
- **README.md** (4.2 KB)
  - Comprehensive documentation of file formats
  - Algorithm details and state type definitions
  - Instructions for Rust testing

- **test_sequences.fa** (353 bytes)
  - FASTA format test sequences
  - Easy reference for sequence data

- **Makefile** (2.2 KB)
  - Build system for regenerating golden files
  - Compiles against original C library

## Test Results

### Model Configuration
- **Model**: TRNA2.cm (standard eukaryotic tRNA model)
- **Nodes**: 72
- **States**: 289 (after RearrangeCM conversion)

### Test Sequences and Scores

| Test Name      | Length | Score      | Description |
|----------------|--------|------------|-------------|
| Phe-73bp       | 73     | **73.8820** | Canonical tRNA-Phe, C. elegans |
| Ser-82bp       | 81     | **66.0850** | Longer tRNA-Ser (82bp vs typical 73bp) |
| Fragment-50bp  | 50     | **-12.8010** | Truncated sequence (edge case) |

### Key Observations

1. **Phe-73bp** achieves the highest score (73.88 bits)
   - Perfect canonical tRNA structure
   - All structural elements properly aligned
   - 73bp is the typical tRNA length

2. **Ser-82bp** scores lower despite being longer (66.09 bits)
   - Some structural deviations from canonical model
   - Extra length doesn't guarantee higher score
   - Model prefers canonical structure

3. **Fragment-50bp** has negative score (-12.80 bits)
   - Incomplete tRNA structure
   - Missing 3' acceptor stem
   - Demonstrates penalty for poor alignment

## Traceback Structure Details

Each traceback tree contains:

- **Root node**: Always `node=0, type=64 (uBEGIN_ST)`
- **Tree depth**: Varies by alignment (typical: 20-30 levels)
- **State types**:
  - `uMATP_ST` (2): Paired positions (base pairs)
  - `uMATL_ST` (4): Left-side single-stranded positions
  - `uMATR_ST` (8): Right-side single-stranded positions
  - `uBIFURC_ST` (256): Bifurcation points (splits tree)
  - `uBEGIN_ST` (64): Begin states at bifurcations
  - `uEND_ST` (128): Terminal states

### Example Traceback Pattern (Phe-73bp)

```
ROOT (uBEGIN_ST) → emits positions 0-72
  ├─ MATR (emits right pos 72)
  │   ├─ MATP (base pair 0,71)
  │   │   ├─ MATP (base pair 1,70)
  │   │   │   └─ ... (acceptor stem)
  │   │   └─ BIFURC
  │   │       ├─ Left subtree (D-loop/anticodon)
  │   │       └─ Right subtree (TΨC/variable loop)
```

## Algorithm Implementation Notes

### Score Calculation
```c
score = (double) bmx[0][N][N] / INTPRECISION;
```
- Uses `INTPRECISION = 1000.0` (3 decimal places)
- Integer log-odds scoring for numerical stability
- Final score divided by INTPRECISION to get bits

### Emission Positions
- **1-based** in alignment coordinates
- **0** indicates no emission (DEL states)
- `emitl`: left emission position
- `emitr`: right emission position

### State Transitions
- Most states have single child (`nxtl`)
- BIFURC states have two children (`nxtl`, `nxtr`)
- END states terminate recursion (`emitl=-1, emitr=-1`)

## Usage for Rust Implementation

### 1. Parse Golden Files

```rust
// Load expected scores
let scores = parse_scores("viterbi_scores.txt");
assert_eq!(scores["Phe-73bp"], 73.8820, "±0.01");

// Load expected traceback
let trace = parse_traceback("traceback.txt", "Phe-73bp");
```

### 2. Validate Rust Implementation

```rust
#[test]
fn test_viterbi_phe_73bp() {
    let cm = load_cm("TRNA2.cm");
    let icm = rearrange_cm(&cm, &[0.25; 4]);
    let seq = "GCCTCGATAGCTCAGTTGGGAGAGCGTACGACTGAAGATCGTAAGGTCACCAGTTCGATCCTGGTTCGGGGCA";

    let (score, trace) = viterbi_align(&icm, seq);

    // Validate score (±0.01 tolerance)
    assert!((score - 73.8820).abs() < 0.01);

    // Validate traceback structure
    assert_eq!(trace.root.nodeidx, 0);
    assert_eq!(trace.root.statetype, uBEGIN_ST);
    assert_eq!(trace.root.emitl, 0);
    assert_eq!(trace.root.emitr, 72);
}
```

### 3. Regenerate if Needed

```bash
cd tests/golden/viterbi
make generate
```

## File Locations

```
tests/golden/viterbi/
├── gen_viterbi.c          # Test generator program
├── Makefile               # Build system
├── viterbi_scores.txt     # Golden scores
├── traceback.txt          # Golden traceback trees
├── test_sequences.fa      # Test sequences in FASTA format
├── README.md              # Detailed documentation
└── SUMMARY.md             # This file
```

## References

- **Original C Implementation**: `../../original/src/viterbi.c`
- **Model File**: `../../original/lib/models/TRNA2.cm`
- **Debug Utilities**: `../../original/src/debug.c`
- **Data Structures**: `../../original/src/structs.h`

## Verification

All golden files have been generated from the original tRNAscan-SE 2.0 C implementation using:
- gcc version 4.8.5 or later
- Mathematical validation against known tRNA structures
- Consistency checks across multiple test cases

Score precision: 4 decimal places (matches INTPRECISION / 1000.0)

## Notes for Rust Developers

1. **Floating-point tolerance**: Use ±0.01 bits tolerance for score comparisons
2. **Zero-based arrays**: Rust uses 0-based indexing; traceback positions are 1-based
3. **Memory management**: Trace trees can be large; consider arena allocation
4. **State types**: Use enums matching the C #define values exactly
5. **Bifurcations**: Special handling needed for nodes with two children

## Status

✅ **Complete**: All Phase 5 golden files generated and validated
✅ **Ready for Rust testing**
✅ **Documentation complete**

Last updated: 2026-03-04
