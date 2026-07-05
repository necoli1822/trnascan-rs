# Phase 8 Golden Files Index

## File Overview

| File | Purpose | Lines | Description |
|------|---------|-------|-------------|
| **README.md** | Overview | 85 | Introduction to isotype scoring, process, and references |
| **IMPLEMENTATION_GUIDE.md** | Implementation | 390 | Complete Rust implementation strategy with code examples |
| **isotype_assignments.txt** | Test Data | 26 | Final isotype assignments for all 11 test tRNAs |
| **isotype_scores.txt** | Test Data | 34 | Full score matrix (11 tRNAs × 23 isotypes) |
| **anticodon_table.txt** | Reference | 122 | Standard genetic code anticodon→isotype mapping |
| **cm_model_info.txt** | Reference | 137 | CM model specifications and scoring details |
| **test_cases.md** | Testing | 286 | Detailed test cases with validation criteria |

## Quick Start

### For Developers
1. **Read first**: `README.md` - Understand what isotype scoring is
2. **Implement**: `IMPLEMENTATION_GUIDE.md` - Follow Rust implementation steps
3. **Test**: `test_cases.md` - Validate against these test cases

### For Testing
```bash
# Run your isotype scorer
./bactars isotype \
  --input tests/data/example_trnas.fa \
  --cm-db original/lib/models/TRNAinf-bact-iso \
  --output test_output.iso

# Validate against golden files
diff test_output.iso tests/golden/isotype/isotype_scores.txt
```

### For Reference
- **anticodon_table.txt**: Look up expected isotype for an anticodon
- **cm_model_info.txt**: Check CM model specifications and scoring details

## Key Test Sequences

### Standard Cases
1. **CELF22B7.trna1**: Leucine (Leu) - High confidence, score 119.9
2. **MySeq1.trna1**: Threonine (Thr) - Medium confidence, score 93.1
3. **MySeq2.trna1**: Arginine (Arg) - High confidence, score 89.3

### Special Cases
4. **CELF22B7.trna2**: Serine (Ser) - Type II tRNA with variable arm, score 125.0
5. **MySeq5.trna1**: Selenocysteine (SeC) - 21st amino acid, score 146.9 ⚠️
6. **MySeq6.trna1**: Lysine (Lys) - Low confidence, unusual score pattern ⚠️

### Edge Cases
7. **CELF22B7.trna3/4**: Phenylalanine (Phe) - Duplicate detection test
8. **CELF22B7.trna3**: Met cross-scoring - Phe vs Met discrimination

## Implementation Checklist

- [ ] Set up Infernal cmsearch invocation
- [ ] Parse tabular output format (--tblout)
- [ ] Extract bit scores for all 23 isotypes
- [ ] Handle score -999 (invalid/cannot score)
- [ ] Select highest-scoring model
- [ ] Calculate confidence from score difference
- [ ] Validate against anticodon prediction
- [ ] Special handling for SeC (Selenocysteine)
- [ ] Special handling for iMet/fMet (Initiator Met)
- [ ] Flag low-confidence assignments
- [ ] Generate output in multiple formats (TSV, JSON, human-readable)

## Validation Criteria

### Isotype Assignment (100% accuracy required)
All 11 test tRNAs must match golden file assignments:
```
CELF22B7.trna1 → Leu ✓
CELF22B7.trna2 → Ser ✓
CELF22B7.trna3 → Phe ✓
CELF22B7.trna4 → Phe ✓
CELF22B7.trna5 → Pro ✓
MySeq1.trna1 → Thr ✓
MySeq2.trna1 → Arg ✓
MySeq3.trna1 → Ser ✓
MySeq4.trna1 → Leu ✓
MySeq5.trna1 → SeC ✓
MySeq6.trna1 → Lys ✓
```

### Score Accuracy (±0.1 bits tolerance)
Top scores must match within 0.1 bits:
- CELF22B7.trna1: 119.9 ± 0.1
- CELF22B7.trna2: 125.0 ± 0.1
- MySeq5.trna1: 146.9 ± 0.1
- etc.

### Confidence Flagging
- High confidence: Score diff > 30 bits (8/11 tRNAs)
- Medium confidence: Score diff 10-30 bits (2/11 tRNAs)
- Low confidence: Score diff < 10 bits (0/11 tRNAs)
- Very low: Negative diff (1/11 tRNAs - MySeq6)

## Dependencies

### Required Software
- **Infernal 1.1+**: cmsearch command must be in PATH
- **CM Database**: `original/lib/models/TRNAinf-bact-iso` (1.2 MB)

### Optional Software
- **SeC Model**: `original/lib/models/TRNAinf-bact-SeC.cm` (for Selenocysteine)

## Common Issues

### 1. cmsearch not found
```
Error: cmsearch command not found
Solution: Install Infernal: conda install -c bioconda infernal
```

### 2. CM database not found
```
Error: Cannot open TRNAinf-bact-iso
Solution: Check path in config, ensure database files are present
```

### 3. Score mismatch
```
Error: Expected Leu=119.9, got Leu=119.8
Solution: Check Infernal version (must be 1.1.1+), verify CM database is correct version
```

### 4. Wrong isotype assigned
```
Error: Expected Leu, got Ser
Solution: Verify score parsing, check if all 23 models were scored
```

## References

### Papers
1. **tRNAscan-SE 2.0**: Chan PP, Lowe TM (2021) *Nucleic Acids Research*
   - doi: 10.1093/nar/gkab688

2. **Infernal**: Nawrocki EP, Eddy SR (2013) *Bioinformatics*
   - doi: 10.1093/bioinformatics/btt509

3. **Covariance Models**: Eddy SR, Durbin R (1994) *Nucleic Acids Research*
   - doi: 10.1093/nar/22.11.2079

### Web Resources
- **tRNAscan-SE**: http://lowelab.ucsc.edu/tRNAscan-SE/
- **Infernal**: http://eddylab.org/infernal/
- **tRNA Database**: http://gtrnadb.ucsc.edu/

### Source Code
- **Original tRNAscan-SE**: `original/tRNAscan-SE.src` (Perl)
- **Isotype Perl Module**: `original/lib/tRNAscanSE/ScanResult.pm`
- **CM Module**: `original/lib/tRNAscanSE/CM.pm`

## Change Log

### 2026-03-04
- Initial creation of Phase 8 golden files
- Extracted from tRNAscan-SE 2.0.12 output
- 11 test sequences covering all major cases
- 7 documentation files created

---

**Total Files**: 7 documentation + reference files
**Total Lines**: 1,080 lines of documentation
**Test Coverage**: 11 tRNAs, 23 isotypes, 253 individual scores
**Implementation Time Estimate**: 2-3 days for complete Rust implementation
