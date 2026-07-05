# Phase 9 Golden Files - Quick Reference Card

## 📁 File Locations

```
tests/golden/
├── full_run/           # Golden reference outputs
│   └── README.md       # Complete documentation
├── pipeline/           # Validation test cases
│   ├── input_validation.txt      # 10 input tests
│   ├── cmdline_args.txt          # 25 argument tests
│   ├── error_cases.txt           # 20 error tests
│   ├── FIELD_VALIDATION.md       # Field specifications
│   └── QUICK_REFERENCE.md        # ← You are here
└── PHASE9_SUMMARY.md   # Overall status
```

---

## 🎯 Golden Reference Files

### Example1 (CELF22B7 - 5 tRNAs)
- `Example1-tRNAs.out` - Tabular output
- `Example1-tRNAs.ss` - Secondary structures
- `Example1-tRNAs.stats` - Statistics
- `Example1-tRNAs.bed` - BED format
- `Example1-tRNAs.iso` - Isotype scores

### Example2 (MySeq1-6 - 6 tRNAs)
- `Example2-tRNAs.out` - Tabular output
- `Example2-tRNAs.ss` - Secondary structures
- `Example2-tRNAs.stats` - Statistics
- `Example2-tRNAs.iso` - Isotype scores

---

## 📊 Test Coverage Summary

| Category | Count | File |
|----------|-------|------|
| Input validation tests | 10 | `pipeline/input_validation.txt` |
| Command-line arg tests | 25 | `pipeline/cmdline_args.txt` |
| Error handling tests | 20 | `pipeline/error_cases.txt` |
| **Total test scenarios** | **55** | |
| Golden tRNA examples | 11 | Example1 (5) + Example2 (6) |
| Output formats documented | 5 | .out, .ss, .stats, .bed, .iso |

---

## 🔑 Critical Output Fields

### Tabular Output (.out) - 15 Fields

| # | Field | Type | Format | Example |
|---|-------|------|--------|---------|
| 1 | Sequence Name | String | Exact match | `CELF22B7` |
| 2 | tRNA # | Integer | Sequential | `1` |
| 3 | Begin | Integer | 1-indexed | `12619` |
| 4 | End | Integer | 1-indexed | `12738` |
| 5 | Type | 3-char | AA code | `Leu` |
| 6 | Codon | 3-char | Uppercase | `CAA` |
| 7 | Intron Begin | Integer | 0 if none | `12657` |
| 8 | Intron End | Integer | 0 if none | `12692` |
| 9 | Inf Score | Float | `%.1f` | `74.2` |
| 10 | HMM Score | Float | `%.2f` | `51.20` |
| 11 | 2'Str Score | Float | `%.2f` | `23.00` |
| 12 | Hit Origin | String | Inf/Ts/Eu/Bo | `Inf` |
| 13 | Isotype CM | 3-char | AA code | `Leu` |
| 14 | Isotype Score | Float | `%.1f` | `119.9` |
| 15 | Note | String | Optional | `ISM (-73.90)` |

---

## 🧪 Special Test Cases

### Forward Strand (Begin < End)
```
CELF22B7  1  12619  12738  Leu  CAA  ...
```

### Reverse Strand (Begin > End)
```
CELF22B7  5  23765  23694  Pro  CGG  ...
```

### With Intron
```
CELF22B7  1  12619  12738  Leu  CAA  12657  12692  ...
```

### Selenocysteine (SeC)
```
MySeq5  1  3  89  SeC  TCA  0  0  146.9  0.00  0.00  Inf  SeC  146.9
```

### Isotype Mismatch (ISM)
```
MySeq6  1  7  92  Lys  CTT  0  0  72.1  40.60  31.50  Inf  Leu  75.7  ISM (-73.90)
```

---

## 📏 Format Precision Rules

| Field | Precision | Rust Format | Example |
|-------|-----------|-------------|---------|
| Inf Score | 1 decimal | `{:.1}` | 74.2 |
| HMM Score | 2 decimals | `{:.2}` | 51.20 |
| 2'Str Score | 2 decimals | `{:.2}` | 23.00 |
| Isotype Score | 1 decimal | `{:.1}` | 119.9 |

---

## 🧮 BED Coordinate Conversion

### Forward Strand (No Intron)
```
Input:  Begin=19480, End=19561
BED:    chromStart=19479, chromEnd=19561
```

### Forward Strand (With Intron)
```
Input:  Begin=12619, End=12738, IntronBegin=12657, IntronEnd=12692
Exon1:  12619-12656 (38 bp)
Exon2:  12693-12738 (46 bp)
BED:    blockCount=2, blockSizes=38,46, blockStarts=0,74
```

### Reverse Strand
```
Input:  Begin=23765, End=23694
BED:    chromStart=23693, chromEnd=23765
```

---

## ✅ Validation Checklist

### Before Claiming Implementation Complete:

- [ ] All 15 tabular output fields match exactly
- [ ] Floating-point precision correct (1 or 2 decimals)
- [ ] Reverse strand coordinates handled (Begin > End)
- [ ] Intron coordinates preserved (or 0,0)
- [ ] BED format uses 0-indexed coordinates
- [ ] SeC tRNAs output correctly
- [ ] ISM note added for mismatches
- [ ] Secondary structure format matches
- [ ] Statistics counts accurate
- [ ] All 5 output formats generated

### Comparison Commands:
```bash
# Tabular
diff -u tests/golden/full_run/Example1-tRNAs.out output/Example1-tRNAs.out

# Secondary structure
diff -u tests/golden/full_run/Example1-tRNAs.ss output/Example1-tRNAs.ss

# BED
diff -u tests/golden/full_run/Example1-tRNAs.bed output/Example1-tRNAs.bed

# Statistics (ignore timestamps)
diff -I "Started:" -I "ended:" -I "host" \
  tests/golden/full_run/Example1-tRNAs.stats output/Example1-tRNAs.stats
```

---

## 📖 Documentation Index

| Document | Purpose | Lines |
|----------|---------|-------|
| `full_run/README.md` | Complete golden file documentation | 484 |
| `pipeline/FIELD_VALIDATION.md` | Exact field specifications | 609 |
| `pipeline/input_validation.txt` | Input handling tests | 97 |
| `pipeline/cmdline_args.txt` | Argument parsing tests | 228 |
| `pipeline/error_cases.txt` | Error handling tests | 263 |
| `PHASE9_SUMMARY.md` | Overall status and next steps | 427 |
| **Total** | | **2,108** |

---

## 🚀 Quick Start for Rust Implementation

### 1. Read the Documentation
```bash
# Start here
cat tests/golden/PHASE9_SUMMARY.md

# Complete reference
cat tests/golden/full_run/README.md

# Field specifications
cat tests/golden/pipeline/FIELD_VALIDATION.md
```

### 2. Run Your Implementation
```bash
# Example1
cargo run -- -E -o output/Example1-tRNAs.out tests/fixtures/Example1.fa

# Example2
cargo run -- -E -o output/Example2-tRNAs.out tests/fixtures/Example2.fa
```

### 3. Compare Outputs
```bash
# Quick check
diff tests/golden/full_run/Example1-tRNAs.out output/Example1-tRNAs.out

# Detailed field-by-field comparison
./tests/scripts/validate_output.sh output/Example1-tRNAs.out
```

### 4. Run Test Suite
```bash
cargo test --test integration_test -- phase9
```

---

## 🎓 Key Concepts

### Coordinate Systems
- **Tabular/SS:** 1-indexed (Begin=1 is first base)
- **BED:** 0-indexed, half-open (chromStart=0 is first base)

### Strand Convention
- **Forward:** Begin < End
- **Reverse:** Begin > End (coordinates still in genomic orientation)

### Score Components
```
Inf Score ≈ HMM Score + 2'Str Score
```

### Isotype Prediction
- **Type:** Predicted isotype (from anticodon)
- **Isotype CM:** Best-scoring isotype model
- **ISM Note:** Added when Type ≠ Isotype CM

---

## 💡 Common Pitfalls

❌ **Wrong:** Using 0-indexed coordinates in tabular output
✅ **Right:** 1-indexed (Begin=1 is first base)

❌ **Wrong:** Truncating scores (51.2 instead of 51.20)
✅ **Right:** Using correct precision (%.2f for HMM)

❌ **Wrong:** Not handling reverse strand
✅ **Right:** Begin > End for reverse strand

❌ **Wrong:** Missing trailing commas in BED
✅ **Right:** `blockSizes=38,46,` (note trailing comma)

❌ **Wrong:** Using spaces instead of tabs
✅ **Right:** All fields tab-separated

---

## 📞 Need More Details?

- **Field formats:** See `FIELD_VALIDATION.md`
- **Complete examples:** See `full_run/README.md`
- **Test scenarios:** See `input_validation.txt`, `cmdline_args.txt`, `error_cases.txt`
- **Overall status:** See `PHASE9_SUMMARY.md`

---

## 📊 Status: ✅ COMPLETE

All Phase 9 golden files and documentation are ready for Rust implementation validation.

**Next Step:** Implement Phase 9 pipeline in Rust and validate against these golden files.
