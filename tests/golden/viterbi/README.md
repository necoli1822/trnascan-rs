# Viterbi Alignment Golden Files

This directory contains golden files for testing the Viterbi alignment algorithm implementation.

## Files

- **gen_viterbi.c** - C test program that generates golden output
- **Makefile** - Build system for compiling and running the generator
- **viterbi_scores.txt** - Expected alignment scores for test sequences
- **traceback.txt** - Expected traceback tree structures

## Test Sequences

Three test sequences are used:

1. **Phe-73bp** - A 73bp tRNA-Phe sequence from C. elegans
   - Canonical tRNA structure
   - Score: ~73.88 bits

2. **Ser-82bp** - An 82bp tRNA-Ser sequence from C. elegans
   - Slightly longer than canonical (81bp vs 73bp)
   - Score: ~66.09 bits

3. **Fragment-50bp** - A 50bp fragment (partial tRNA)
   - Truncated sequence for edge case testing
   - Score: ~-12.80 bits (negative, as expected for incomplete structure)

## Model

All tests use the **TRNA2.cm** model (standard eukaryotic tRNA covariance model):
- 72 nodes
- 289 states after conversion to integer model
- Located at: `../../original/lib/models/TRNA2.cm`

## Algorithm Details

The Viterbi algorithm implementation tests:

1. **RearrangeCM()** - Converts floating-point CM to integer-based istate_s model
   - Uses INTPRECISION = 1000.0 for 3 decimal places
   - Log-odds scoring with integer precision

2. **ViterbiAlign()** - Performs alignment with dynamic programming
   - Returns score: `bmx[0][N][N] / INTPRECISION`
   - Returns traceback tree structure

3. **Traceback Tree** - Binary tree structure with:
   - `nodeidx` - Node index in the model (0-71)
   - `type` - State type (uMATP_ST, uMATL_ST, uMATR_ST, etc.)
   - `emitl` - Left emission position (1-based, 0 if none)
   - `emitr` - Right emission position (1-based, 0 if none)
   - `nxtl` - Left child pointer
   - `nxtr` - Right child pointer (for BIFURC states)

## State Types (Unique State Type Flags)

- `uDEL_ST` (1) - Deletion state
- `uMATP_ST` (2) - Match pair state (emits both left and right)
- `uMATL_ST` (4) - Match left state
- `uMATR_ST` (8) - Match right state
- `uINSL_ST` (16) - Insert left state
- `uINSR_ST` (32) - Insert right state
- `uBEGIN_ST` (64) - Begin state
- `uEND_ST` (128) - End state
- `uBIFURC_ST` (256) - Bifurcation state

## Output Format

### viterbi_scores.txt

```
<test_name> <sequence_length> <score>
```

Example:
```
Phe-73bp	73	73.8820
```

### traceback.txt

Depth-first traversal of the traceback tree:

```
node=<nodeidx> type=<type> emitl=<emitl> emitr=<emitr> (<statetype_name>)
```

Indentation shows tree depth. Example:

```
node=0 type=64 emitl=0 emitr=72 (uBEGIN_ST)
  node=1 type=8 emitl=0 emitr=80 (uMATR_ST)
    node=2 type=2 emitl=0 emitr=79 (uMATP_ST)
      ...
```

## Building and Running

### Compile the generator:

```bash
make gen_viterbi
```

### Generate new golden files:

```bash
make generate
# Or manually:
./gen_viterbi ../../../original/lib/models/TRNA2.cm
```

### Clean:

```bash
make clean
```

## Using in Rust Tests

For Rust implementation testing, parse these golden files:

1. Parse `viterbi_scores.txt` to get expected scores
2. Parse `traceback.txt` to verify traceback structure
3. Compare Rust implementation outputs against these values
4. Tolerance: ±0.01 bits for scores (due to floating-point precision)

Example Rust test structure:

```rust
#[test]
fn test_viterbi_phe_73bp() {
    let seq = "GCCTCGATAGCTCAGTTGGGAGAGCGTACGACTGAAGATCGTAAGGTCACCAGTTCGATCCTGGTTCGGGGCA";
    let (score, trace) = viterbi_align(icm, seq);

    // Check score
    assert!((score - 73.8820).abs() < 0.01);

    // Check traceback structure
    assert_eq!(trace.root.nodeidx, 0);
    assert_eq!(trace.root.statetype, uBEGIN_ST);
    // ... more assertions
}
```

## Notes

- Sequence positions in traceback are 0-based in the tree but represent 1-based alignment positions
- END states have `emitl=-1, emitr=-1` to indicate no emission
- BIFURC states have both `nxtl` and `nxtr` children; other states only use `nxtl`
- The score is computed from the final BEGIN state at `bmx[0][N][N]`

## References

- Original implementation: `../../original/src/viterbi.c`
- Debug utilities: `../../original/src/debug.c` (PrintTrace function)
- Model format: `../../original/src/save.c` (ReadCM function)
