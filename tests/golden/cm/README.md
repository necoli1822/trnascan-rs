# CM Structure Golden Files

This directory contains golden reference files for Phase 4 (CM structure parsing) of the tRNAscan-SE Rust reimplementation.

## Generated Files

These golden files were generated from `TRNA2.cm` (72 nodes, 289 states after RearrangeCM):

- **cm_structure.txt** (82 lines) - High-level CM structure overview
  - Number of nodes: 72
  - Constants: STATETYPES=6, ALPHASIZE=4, NODETYPES=7
  - Per-node: type, nxt, nxt2 connections

- **node_details.txt** (1371 lines) - Detailed per-node information
  - Node type and connections
  - Transition matrix [6x6]
  - MATP emissions [4x4]
  - INSL/INSR/MATL/MATR emissions [4]

- **istate_dump.txt** (2321 lines) - Integer state array after RearrangeCM
  - Total states: 289
  - Per-state: nodeidx, statetype, offset, connectnum, tmx[6], emit[4 or 16]
  - Verifies istate_s has nodeidx (not bifr)

## Source Model

**Input:** `/mnt/DAS/sunju/programme/bactars/tRNAscan-SE/original/lib/models/TRNA2.cm`
- Format: Binary CM format (v20magic = 0xe3edb2b0)
- Size: 47K
- Purpose: Standard tRNA covariance model

## Generator Program

**Source:** `gen_cm.c`
- Reads CM using `ReadCM()` from original codebase
- Converts to integer states using `RearrangeCM()`
- Dumps all structure fields for verification

**Build:**
```bash
cd /mnt/DAS/sunju/programme/bactars/tRNAscan-SE/tests/golden/cm
gcc -I../../original/src -I../../original -o gen_cm gen_cm.c \
  ../../original/src/model.o \
  ../../original/src/save.o \
  ../../original/src/sqio.o \
  ../../original/src/iupac.o \
  ../../original/src/misc.o \
  ../../original/src/sqerror.o \
  ../../original/src/sre_math.o \
  ../../original/src/structs.o \
  ../../original/src/sre_string.o \
  ../../original/src/sre_ctype.o \
  ../../original/src/stack.o \
  ../../original/src/selex.o \
  ../../original/src/msf.o \
  ../../original/src/interleaved.o \
  ../../original/src/alignio.o \
  ../../original/src/types.o \
  -lm
```

**Run:**
```bash
./gen_cm /path/to/model.cm
```

## Critical Verification Points

### Structure Verification

1. **Constants (cm_structure.txt)**
   - v20magic = 0xe3edb2b0 (binary format magic)
   - STATETYPES = 6
   - ALPHASIZE = 4
   - NODETYPES = 7

2. **Node Types (node_details.txt)**
   - 0 = BIFURC_NODE (also END_NODE)
   - 1 = MATP_NODE (paired match)
   - 2 = MATL_NODE (left match)
   - 3 = MATR_NODE (right match)
   - 4 = BEGINL_NODE
   - 5 = BEGINR_NODE
   - 6 = ROOT_NODE

3. **State Types (istate_dump.txt)**
   - Unique flags (not indices):
     - 0x01 = uDEL_ST
     - 0x02 = uMATP_ST
     - 0x04 = uMATL_ST
     - 0x08 = uMATR_ST
     - 0x10 = uINSL_ST
     - 0x20 = uINSR_ST
     - 0x40 = uBEGIN_ST
     - 0x80 = uEND_ST
     - 0x100 = uBIFURC_ST

### Field Differences

**istate_s** (integer states for alignment):
- Has: nodeidx, statetype, offset, connectnum, tmx[6], emit[]
- **Does NOT have:** bifr

**pstate_s** (probability states):
- Has: nodeidx, statetype, offset, connectnum, **bifr**, tmx[6], emit[]
- bifr only used for BIFURC states

### Emission Arrays

- **MATP_ST:** emit[16] (4x4 paired emissions)
- **MATL_ST, MATR_ST, INSL_ST, INSR_ST:** emit[4] (single position)
- **DEL_ST, BEGIN_ST, END_ST, BIFURC_ST:** No emissions

### Transition Order

State transitions rearranged in istate_s:
**Order:** INSL, INSR, DEL, MATP, MATL, MATR

(Different from node_s order: DEL, MATP, MATL, MATR, INSL, INSR)

## Usage in Rust Tests

These golden files should be used to verify:

1. **CM Parser** - Correctly reads binary CM files
2. **Structure Conversion** - Properly converts cm_s to istate_s
3. **Field Values** - Exact match on transition/emission probabilities
4. **Log-Odds Conversion** - Correct integer conversion from doubles

Example test:
```rust
#[test]
fn test_cm_structure_matches_golden() {
    let cm = read_cm("tests/models/TRNA2.cm").unwrap();
    assert_eq!(cm.nodes, 72);

    // Verify node 0 (ROOT_NODE)
    assert_eq!(cm.nd[0].type_, 6); // ROOT_NODE
    assert_eq!(cm.nd[0].nxt, 1);
    assert_eq!(cm.nd[0].nxt2, -1);

    // Compare against node_details.txt
    assert_approx_eq!(cm.nd[0].tmx[0][0], 0.008460);
    // ... more assertions
}
```

## Regeneration

To regenerate golden files (if source model changes):
```bash
cd /mnt/DAS/sunju/programme/bactars/tRNAscan-SE/tests/golden/cm
./gen_cm ../../original/lib/models/TRNA2.cm
```

## Date Generated

2026-03-04
