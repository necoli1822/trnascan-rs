# Isotype Scoring Test Cases

## Purpose
Validate that the Rust implementation correctly:
1. Invokes Infernal's cmsearch with the correct CM database
2. Parses isotype-specific bit scores from cmsearch output
3. Selects the highest-scoring model as the predicted isotype
4. Handles special cases (SeC, iMet, score ties)
5. Reports accurate score differences for confidence estimation

## Test Data Source
- **Example1.fa**: CELF22B7 sequence (5 tRNAs)
- **Example2.fa**: Multiple sequences (6 tRNAs, including SeC)
- **Golden files**: Extracted from tRNAscan-SE 2.0.12 output

---

## Test Case 1: Standard Isotype (Leucine)

### Input
- **Sequence**: CELF22B7.trna1
- **Anticodon**: AAG
- **Expected**: Leu

### Expected Scores (Top 5)
```
Leu: 119.9  ← Winner
Ser: 69.0
Arg: 55.3
Thr: 52.1
Cys: 51.2
```

### Validation
- ✓ Leu has highest score
- ✓ Score difference (119.9 - 69.0 = 50.9) indicates high confidence
- ✓ Anticodon AAG matches Leu codon family (CUU/CUC/CUA/CUG/CUU/CUC)

---

## Test Case 2: Type II tRNA (Serine with Variable Arm)

### Input
- **Sequence**: CELF22B7.trna2
- **Anticodon**: CGA
- **Expected**: Ser

### Expected Scores (Top 5)
```
Ser: 125.0  ← Winner
Leu: 72.1
Cys: 56.8
His: 45.0
Thr: 43.9
```

### Validation
- ✓ Ser has highest score
- ✓ High score (125.0) reflects good structural match to Ser variable arm
- ✓ Leu also scores well (72.1) because it also has a long variable arm
- ✓ Score difference (125.0 - 72.1 = 52.9) is decisive

---

## Test Case 3: Phenylalanine (High Met Cross-Score)

### Input
- **Sequence**: CELF22B7.trna3 and CELF22B7.trna4 (identical)
- **Anticodon**: GAA
- **Expected**: Phe

### Expected Scores (Top 5)
```
Phe: 112.1  ← Winner
Met: 85.0
Thr: 84.2
Tyr: 81.3
Ile: 73.6
```

### Validation
- ✓ Phe has highest score
- ✓ Met scores high (85.0) but not high enough
- ✓ Score difference (112.1 - 85.0 = 27.1) is moderate
- ✓ Variable loop structure distinguishes Phe from Met
- ✓ Both trna3 and trna4 have identical scores (duplicate detection test)

---

## Test Case 4: Special Case - Selenocysteine (SeC)

### Input
- **Sequence**: MySeq5.trna1
- **Anticodon**: TCA
- **Expected**: SeC

### Expected Scores
```
SeC: 146.9  ← Winner (extremely high!)
Asn: 5.4
Asp: -2.1
Lys: -2.0
Gly: -2.0
[Most others: -999 or negative]
```

### Validation
- ✓ SeC has exceptionally high score (146.9)
- ✓ Massive score difference (146.9 - 5.4 = 141.5)
- ✓ Most other models cannot score (-999)
- ✓ Recognizes unique SeC structure:
  - 8 bp acceptor stem (vs. 7 bp standard)
  - 6 bp D-stem (vs. 4 bp standard)
- ✓ Anticodon TCA matches UGA stop codon suppression

---

## Test Case 5: Low Score but Correct (Lysine Outlier)

### Input
- **Sequence**: MySeq6.trna1
- **Anticodon**: CTT
- **Expected**: Lys

### Expected Scores (Top 5)
```
Ser: 72.4
Leu: 75.7
His: 69.9
Ile: 55.8
Tyr: 54.8
...
Lys: 1.8  ← Winner (but very low!)
```

### Validation
- ✓ Lys score is positive (1.8) but very low
- ✓ Other models score higher but are incorrect
- ⚠️ Negative score difference (1.8 - 72.4 = -70.6) signals low confidence
- ⚠️ May indicate:
  - Degenerate tRNA structure
  - Potential pseudogene
  - Unusual modifications
  - Need for manual review

**Expected Behavior**: Flag this tRNA for manual inspection due to low confidence.

---

## Test Case 6: Met vs. iMet Discrimination

### Context
- **Met**: Elongator methionine (normal tRNA)
- **iMet**: Initiator methionine (translation start)
- **Structural difference**: 3 mismatches in acceptor stem

### Test Data
From the score matrices, observe:
- Met and iMet have different scores for the same sequence
- Example: CELF22B7.trna3
  - Met: 85.0
  - iMet: -17.7

### Validation
- ✓ iMet scores very low (-17.7) for non-initiator tRNAs
- ✓ Models correctly distinguish elongator vs. initiator forms
- ✓ Special handling may be needed for bacterial fMet

---

## Test Case 7: Score -999 Handling

### Purpose
Verify correct handling of models that cannot score a sequence.

### Examples
From CELF22B7.trna2:
```
Pro: -999  (No variable loop in Pro model)
iMet: -999 (Not an initiator structure)
```

From MySeq5.trna1 (SeC):
```
Ala: -999
Arg: -999
Glu: -999
Met: -999
Phe: -999
Thr: -999
Val: -999
...
```

### Validation
- ✓ Score -999 is recognized as "model cannot score"
- ✓ These models are excluded from isotype selection
- ✓ Only scored models participate in ranking

---

## Edge Cases to Test

### 1. All Models Score -999
- **Action**: Report error, cannot determine isotype
- **Reason**: Invalid tRNA structure

### 2. Tie Scores
- **Example**: Two models have identical top scores
- **Action**: Use tie-breaking rules:
  1. Prefer anticodon-predicted isotype
  2. Prefer simpler isotype (Met over iMet)
  3. Report both as equally likely

### 3. Anticodon Mismatch
- **Example**: Anticodon says Ala, but Leu model scores highest
- **Action**: Report mismatch, flag for review
- **Reason**: May indicate:
  - Anticodon editing
  - Sequencing error
  - Novel tRNA function

### 4. Multiple SeC Models
- **Action**: If both SeC and Sec models exist, use higher score
- **Reason**: Naming inconsistency (SeC vs Sec vs Sel)

---

## Performance Requirements

### Timing
- **Single tRNA**: < 1 second (including cmsearch invocation)
- **Batch (11 tRNAs)**: < 10 seconds

### Memory
- **CM database**: ~1-2 MB loaded once
- **Per tRNA**: < 10 KB for score storage

### Accuracy
- **Isotype assignment**: 100% match with tRNAscan-SE 2.0.12
- **Score precision**: Within 0.1 bits of reference

---

## Validation Script

```bash
# Run isotype scoring on test data
./bactars isotype \
  --input tests/golden/full_run/Example1-tRNAs.fa \
  --cm-db original/lib/models/TRNAinf-bact-iso \
  --output test_output.iso

# Compare with golden file
diff test_output.iso tests/golden/isotype/isotype_scores.txt

# Check assignments
./scripts/check_isotype_assignments.sh \
  test_output.iso \
  tests/golden/isotype/isotype_assignments.txt
```

---

## Expected Output Format

### Isotype Scores (Tab-separated)
```
tRNAscanID	Anticodon_predicted_isotype	Ala	Arg	Asn	...	Val	iMet
CELF22B7.trna1	Leu	37.4	55.3	10.1	...	28.8	-999
```

### Assignment Report
```
CELF22B7.trna1: Leu (score=119.9, runner-up=Ser:69.0, diff=50.9, confidence=HIGH)
CELF22B7.trna2: Ser (score=125.0, runner-up=Leu:72.1, diff=52.9, confidence=HIGH)
MySeq6.trna1: Lys (score=1.8, runner-up=Ser:72.4, diff=-70.6, confidence=LOW) ⚠️
```

---

## References

1. **tRNAscan-SE Paper**: Chan PP, Lowe TM (2021) Nucleic Acids Res.
2. **Infernal Manual**: http://eddylab.org/infernal/Userguide.pdf
3. **tRNA Database**: http://trna.ucsc.edu/
