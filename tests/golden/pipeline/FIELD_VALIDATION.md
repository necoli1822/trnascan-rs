# Critical Output Field Validation

## Purpose

This document specifies exact formatting and validation rules for all output fields in the tRNAscan-SE pipeline. Use this as a reference when implementing output formatting in the Rust version.

---

## Tabular Output (.out) Field Validation

### Field 1: Sequence Name
- **Type:** String
- **Constraints:**
  - Must exactly match sequence name from input FASTA
  - Preserve all characters (including special chars)
  - No truncation
- **Padding:** Left-aligned, minimum 8 characters
- **Test Cases:**
  ```
  CELF22B7    → "CELF22B7 "
  MySeq1      → "MySeq1  "
  seq|x:y/z   → "seq|x:y/z"
  ```

### Field 2: tRNA #
- **Type:** Integer
- **Constraints:**
  - Sequential numbering starting at 1
  - Resets for each new sequence
  - Always positive
- **Format:** Right-aligned, tab-separated
- **Test Cases:**
  ```
  First tRNA in sequence  → 1
  Fifth tRNA in sequence  → 5
  ```

### Field 3: Begin
- **Type:** Integer (1-indexed position)
- **Constraints:**
  - Must be ≥ 1
  - For forward strand: Begin < End
  - For reverse strand: Begin > End
- **Format:** Right-aligned, tab-separated
- **Test Cases:**
  ```
  Forward: 12619 < 12738  → Begin=12619
  Reverse: 23765 > 23694  → Begin=23765
  ```

### Field 4: End
- **Type:** Integer (1-indexed position)
- **Constraints:**
  - Must be ≥ 1
  - Absolute value of |End - Begin| + 1 = tRNA length
- **Format:** Right-aligned, tab-separated
- **Test Cases:**
  ```
  Forward: Begin=12619, End=12738 → Length=120
  Reverse: Begin=23765, End=23694 → Length=72
  ```

### Field 5: Type (tRNA Isotype)
- **Type:** String (3-letter amino acid code or SeC)
- **Constraints:**
  - Valid codes: Ala, Arg, Asn, Asp, Cys, Gln, Glu, Gly, His, Ile, Leu, Lys, Met, Phe, Pro, Ser, Thr, Trp, Tyr, Val, SeC, Sup, Undet
  - SeC = Selenocysteine
  - Sup = Suppressor
  - Undet = Undetermined
- **Format:** Left-aligned, 4 characters wide
- **Test Cases:**
  ```
  Leucine        → "Leu "
  Selenocysteine → "SeC "
  Undetermined   → "Undet"
  ```

### Field 6: Codon (Anticodon)
- **Type:** String (3-letter nucleotide sequence)
- **Constraints:**
  - Valid characters: A, C, G, T (uppercase)
  - Always 3 characters
  - ??? for undetermined
- **Format:** Left-aligned, 5 characters wide
- **Test Cases:**
  ```
  CAA → "CAA  "
  AGA → "AGA  "
  ??? → "???  "
  ```

### Field 7: Intron Begin
- **Type:** Integer (1-indexed, relative to sequence start)
- **Constraints:**
  - 0 if no intron
  - If non-zero: Must be > Begin (for forward) or < Begin (for reverse)
- **Format:** Right-aligned, tab-separated
- **Test Cases:**
  ```
  No intron           → 0
  Intron at 12657     → 12657
  ```

### Field 8: Intron End
- **Type:** Integer (1-indexed)
- **Constraints:**
  - 0 if no intron
  - If non-zero: Must be paired with non-zero Intron Begin
  - For forward: Intron End > Intron Begin
- **Format:** Right-aligned, tab-separated
- **Test Cases:**
  ```
  No intron           → 0
  Intron 12657-12692  → 12692
  ```

### Field 9: Inf Score (Infernal Score)
- **Type:** Float (bits)
- **Format:** `%.1f` (exactly 1 decimal place)
- **Constraints:**
  - Typically 0.0 to 200.0
  - Can be negative for poor matches
- **Test Cases:**
  ```
  74.2     → "74.2"
  146.9    → "146.9"
  8.3      → "8.3"
  -5.1     → "-5.1"
  ```

### Field 10: HMM Score
- **Type:** Float (bits)
- **Format:** `%.2f` (exactly 2 decimal places)
- **Constraints:**
  - Typically positive
  - Can be 0.00 for special cases (SeC tRNA)
- **Test Cases:**
  ```
  51.20    → "51.20"
  0.00     → "0.00"
  43.90    → "43.90"
  123.45   → "123.45"
  ```

### Field 11: 2'Str Score (Secondary Structure Score)
- **Type:** Float (bits)
- **Format:** `%.2f` (exactly 2 decimal places)
- **Constraints:**
  - Typically positive
  - Can be 0.00 for special cases
  - Inf Score ≈ HMM Score + 2'Str Score (approximately)
- **Test Cases:**
  ```
  23.00    → "23.00"
  34.10    → "34.10"
  0.00     → "0.00"
  ```

### Field 12: Hit Origin
- **Type:** String (3-letter code)
- **Constraints:**
  - Valid codes: Inf, Ts, Eu, Bo, Cove
  - Inf = Infernal
  - Ts = tRNAscan 1.4
  - Eu = EufindtRNA
  - Bo = Both Ts and Eu
  - Cove = COVE (legacy)
- **Format:** Left-aligned, minimum 3 characters
- **Test Cases:**
  ```
  Infernal detection    → "Inf"
  Legacy tRNAscan       → "Ts "
  EufindtRNA            → "Eu "
  ```

### Field 13: Isotype CM (Best Isotype Model)
- **Type:** String (3-letter amino acid code)
- **Constraints:**
  - Same codes as Type field
  - Usually matches Type
  - May differ for isotype mismatches
- **Format:** Left-aligned, minimum 3 characters
- **Test Cases:**
  ```
  Match: Type=Leu, Isotype=Leu     → "Leu"
  Mismatch: Type=Lys, Isotype=Leu  → "Leu"
  ```

### Field 14: Isotype Score
- **Type:** Float (bits)
- **Format:** `%.1f` (exactly 1 decimal place)
- **Constraints:**
  - Score from isotype-specific covariance model
  - Usually higher than Inf Score
- **Test Cases:**
  ```
  119.9    → "119.9"
  125.0    → "125.0"
  -73.90   → "-73.9" (for ISM note)
  ```

### Field 15: Note
- **Type:** String (optional)
- **Constraints:**
  - Empty if no special annotation
  - Common values:
    - `ISM (X.XX)` - Isotype mismatch with score
    - `pseudo` - Predicted pseudogene
  - Tab before note if present
- **Test Cases:**
  ```
  Normal tRNA           → ""
  Isotype mismatch      → "ISM (-73.90)"
  Pseudogene            → "pseudo"
  ```

---

## Secondary Structure (.ss) Field Validation

### Header Line
- **Format:** `SeqName.trnaX (Begin-End)\tLength: N bp`
- **Test Cases:**
  ```
  CELF22B7.trna1 (12619-12738)	Length: 120 bp
  MySeq3.trna1 (14-114)	Length: 101 bp
  ```

### Type/Anticodon Line
- **Format:** `Type: XXX\tAnticodon: YYY at A-B (C-D)\tScore: S.S`
- **Constraints:**
  - A-B: positions relative to tRNA (1-indexed)
  - C-D: genomic positions (1-indexed)
- **Test Cases:**
  ```
  Type: Leu	Anticodon: CAA at 35-37 (12653-12655)	Score: 74.2
  Type: Ser	Anticodon: AGA at 34-36 (19513-19515)	Score: 81.6
  ```

### Intron Annotation Line (Optional)
- **Format:** `Possible intron: A-B (C-D)`
- **Only present if intron detected**
- **Test Cases:**
  ```
  Possible intron: 39-74 (12657-12692)
  Possible intron: 51-69 (64-82)
  ```

### Score Breakdown Line
- **Format:** `HMM Sc=X.XX\tSec struct Sc=Y.YY`
- **Test Cases:**
  ```
  HMM Sc=51.20	Sec struct Sc=23.00
  HMM Sc=0.00	Sec struct Sc=0.00
  ```

### Ruler Line
- **Format:** Repeating ` *    |` pattern (6 chars per unit)
- **Length:** Matches sequence length
- **Test Cases:**
  ```
  For 73bp: "         *    |    *    |    *    |    *    |    *    |    *    |    *    |  "
  For 82bp: "         *    |    *    |    *    |    *    |    *    |    *    |    *    |    *    | "
  ```

### Sequence Line
- **Format:** `Seq: NNNN...`
- **Constraints:**
  - Uppercase nucleotides
  - Anticodon in lowercase
  - Intron (if present) in lowercase
- **Test Cases:**
  ```
  No intron:  "Seq: GCAGTCATGTCCGAGTGGTtAAGGAGATTGACTAGAAATCAATTGGGCTCTGCCCGCGTAGGTTCGAATCCTGCTGACTGCG"
  With intron: "Seq: GCACGGATGGCCGAGTGGTctAAGGCGCCAGACTCAAGcgaaatgcttgcctcatgctcgaggtcgactgggtgTTCTGGTACTCGTATGGGTGCGTGGGTTCGAATCCCACTTCGTGCA"
  ```

### Structure Line
- **Format:** `Str: >>>...<<<`
- **Constraints:**
  - Same length as sequence
  - `>` = 5' side of helix
  - `<` = 3' side of helix
  - `.` = unpaired
- **Test Cases:**
  ```
  "Str: >>>>>>>..>>>..........<<<.>>>>>.......<<<<<.>>>>...<<<<..>>>>>.......<<<<<<<<<<<<."
  "Str: >>>>>>>..>>>...........<<<.>>>>>...........................................<<<<<.>>>>....<<<<..>>>>>.......<<<<<<<<<<<<."
  ```

---

## Statistics (.stats) Field Validation

### Header Section
- **Format:**
  ```
  tRNAscan-SE v.X.X (Month Year) scan results (on host hostname)
  Started: Day Mon DD HH:MM:SS TZ YYYY
  ```
- **Note:** Hostname and timestamp will differ - these should be excluded from validation

### Parameters Section
- **Required fields:**
  - Sequence file(s) to search
  - Search Mode
  - Results written to
  - Output format
  - Searching with
  - Isotype-specific model scan
  - Covariance model(s)
  - Cutoff score
  - Temporary directory
  - Output file paths

### First-pass Stats
- **Required fields:**
  - Sequences read (integer)
  - Seqs w/at least 1 hit (integer)
  - Bases read (integer with note "x2 for both strands")
  - Bases in tRNAs (integer)
  - tRNAs predicted (integer)
  - Av. tRNA length (integer)
  - Script CPU time (float, s)
  - Scan CPU time (float, s)
  - Scan speed (float, Kbp/sec or bp/sec)

### Infernal Stats
- **Required fields:**
  - Candidate tRNAs read
  - Infernal-confirmed tRNAs
  - Bases scanned by Infernal
  - % seq scanned by Infernal
  - CPU times
  - Scan speed

### Summary Counts
- **Format:**
  ```
  tRNAs decoding Standard 20 AA:              N
  Selenocysteine tRNAs (TCA):                 N
  Possible suppressor tRNAs (CTA,TTA):        N
  tRNAs with undetermined/unknown isotypes:   N
  tRNAs with mismatch isotypes:               N
  Predicted pseudogenes:                      N
                                              -------
  Total tRNAs:                                N
  ```

### Intron Summary
- **Format:**
  ```
  tRNAs with introns:     	N

  | Type-Codon: count | ...
  ```

### Isotype Table
- **Format:** Aligned columns showing isotype and anticodon counts
- **Structure:**
  ```
  Isotype : total (with_introns)  Codon1: N  Codon2: N  ...
  ```

---

## BED Format (.bed) Field Validation

### BED12 Fields

1. **chrom** - Sequence name (string)
2. **chromStart** - 0-indexed start (Begin - 1)
3. **chromEnd** - End position (unchanged for forward, Begin for reverse)
4. **name** - Format: `SeqName.tRNAX-TypeCodon`
5. **score** - Inf Score × 10, rounded to integer
6. **strand** - `+` or `-`
7. **thickStart** - Same as chromStart
8. **thickEnd** - Same as chromEnd
9. **itemRgb** - 0 (or color code)
10. **blockCount** - 1 (no intron) or 2 (with intron)
11. **blockSizes** - Comma-separated sizes
12. **blockStarts** - Comma-separated starts (relative to chromStart)

### Test Cases

**No Intron (Forward):**
```
Input:  Begin=19480, End=19561, Strand=+, Length=82, Score=81.6
Output: CELF22B7	19479	19561	CELF22B7.tRNA2-SerAGA	816	+	19479	19561	0	1	82,	0,
```

**With Intron (Forward):**
```
Input:  Begin=12619, End=12738, Intron=12657-12692, Score=74.2
        Exon1 Length: 12657-12619 = 38
        Exon2 Length: 12738-12692 = 46
Output: CELF22B7	12618	12738	CELF22B7.tRNA1-LeuCAA	742	+	12618	12738	0	2	38,46,	0,74,
```

**No Intron (Reverse):**
```
Input:  Begin=23765, End=23694, Strand=-, Length=72, Score=71.5
Output: CELF22B7	23693	23765	CELF22B7.tRNA5-ProCGG	715	-	23693	23765	0	1	72,	0,
```

---

## Isotype Scores (.iso) Field Validation

### Header
- **Format:**
  ```
  tRNAscanID	Anticodon_predicted_isotype	Ala	Arg	Asn	Asp	Cys	Gln	Glu	Gly	His	Ile	Leu	Lys	Met	Phe	Pro	SeC	Ser	Thr	Trp	Tyr	Val	iMet
  ```
- **All fields tab-separated**

### Data Rows
- **Format:** `SeqName.trnaX\tType\tScore1\tScore2\t...`
- **Score Format:** Float with 1 decimal place
- **Special Value:** `-999` for models not tested
- **Test Cases:**
  ```
  CELF22B7.trna1	Leu	37.4	55.3	10.1	...	119.9	...	-999
  CELF22B7.trna5	Pro	74.8	50.4	13.0	...	113.0	...	-999
  ```

---

## Validation Checklist

When implementing output formatting, verify:

### ✓ Tabular Output (.out)
- [ ] Header is exactly 3 lines
- [ ] All fields tab-separated
- [ ] Sequence name matches input exactly
- [ ] Begin/End coordinates are 1-indexed
- [ ] Reverse strand: Begin > End
- [ ] Intron coordinates: 0,0 if none
- [ ] Inf Score: 1 decimal (`%.1f`)
- [ ] HMM/2'Str Scores: 2 decimals (`%.2f`)
- [ ] Type codes are valid 3-letter codes
- [ ] Anticodon is 3 uppercase letters or ???
- [ ] Hit Origin is valid code (Inf, Ts, Eu, Bo)
- [ ] Isotype Score: 1 decimal
- [ ] Note field present only if needed

### ✓ Secondary Structure (.ss)
- [ ] Coordinate format: `(Begin-End)`
- [ ] Length calculation correct
- [ ] Anticodon positions (both relative and genomic)
- [ ] Intron annotation if present
- [ ] Score breakdown: HMM and Sec struct
- [ ] Ruler line matches sequence length
- [ ] Sequence: anticodon in lowercase
- [ ] Sequence: intron in lowercase (if present)
- [ ] Structure notation matches sequence length
- [ ] Structure notation valid: `>`, `<`, `.`

### ✓ Statistics (.stats)
- [ ] All required sections present
- [ ] Parameter values match command line
- [ ] Counts are accurate
- [ ] Summary totals match individual counts
- [ ] Isotype table formatted correctly
- [ ] Intron summary correct

### ✓ BED Format (.bed)
- [ ] Positions are 0-indexed (chromStart = Begin - 1)
- [ ] Score = Inf Score × 10 (integer)
- [ ] Name format: `SeqName.tRNAX-TypeCodon`
- [ ] Strand: `+` or `-`
- [ ] Block counts correct (1 or 2)
- [ ] Block sizes and starts correct for introns
- [ ] All fields tab-separated
- [ ] Trailing commas on block fields

### ✓ Isotype Scores (.iso)
- [ ] Header row with all isotypes
- [ ] Tab-separated values
- [ ] tRNAscanID format: `SeqName.trnaX`
- [ ] Predicted isotype in column 2
- [ ] Score format: 1 decimal or -999
- [ ] All isotypes covered

---

## Numerical Precision Rules

### Floating-Point Comparisons
When comparing floating-point values:
- **Inf Score:** `|expected - actual| < 0.05` (tolerance: 0.05 bits)
- **HMM Score:** `|expected - actual| < 0.005` (tolerance: 0.005 bits)
- **2'Str Score:** `|expected - actual| < 0.005`
- **Isotype Score:** `|expected - actual| < 0.05`

### Rounding Rules
- **All scores:** Round half-up (0.5 → 1)
- **BED score:** `floor(Inf Score × 10 + 0.5)`

### Display Format
Use Rust's format strings:
```rust
format!("{:.1}", inf_score)      // Inf Score: 74.2
format!("{:.2}", hmm_score)      // HMM Score: 51.20
format!("{:.2}", ss_score)       // 2'Str Score: 23.00
format!("{:.1}", isotype_score)  // Isotype: 119.9
```

---

## Edge Cases

### 1. Reverse Strand Coordinates
- **Input:** Begin=23765, End=23694
- **Validation:** Begin > End indicates reverse strand
- **BED conversion:** chromStart=23693, chromEnd=23765

### 2. Selenocysteine tRNAs
- **Type:** SeC
- **Anticodon:** TCA
- **Scores:** May have HMM=0.00, 2'Str=0.00, Inf>0
- **No pseudogene check**

### 3. Isotype Mismatches
- **Type ≠ Isotype CM**
- **Note field:** `ISM (score_difference)`
- **Example:** Type=Lys, Isotype CM=Leu, Note=ISM (-73.90)

### 4. Very Long Introns
- **Length > 100 bp**
- **Still valid if detected**
- **Check block sizes in BED format**

### 5. Overlapping tRNAs
- **Same sequence can have multiple tRNAs**
- **May overlap on different strands**
- **Each gets unique tRNA #**

### 6. tRNAs at Sequence Ends
- **Begin=1 is valid**
- **End=sequence_length is valid**

---

## Testing Strategy

### Unit Tests
Test each field formatter independently:
```rust
#[test]
fn test_format_inf_score() {
    assert_eq!(format_inf_score(74.23), "74.2");
    assert_eq!(format_inf_score(146.89), "146.9");
    assert_eq!(format_inf_score(-5.12), "-5.1");
}
```

### Integration Tests
Compare full output against golden files:
```rust
#[test]
fn test_example1_output() {
    let expected = read_golden("Example1-tRNAs.out");
    let actual = run_pipeline("Example1.fa");
    assert_outputs_match(expected, actual);
}
```

### Field-Level Validation
Parse and compare field by field:
```rust
#[test]
fn test_example1_coordinates() {
    let results = parse_output("Example1-tRNAs.out");
    assert_eq!(results[0].sequence_name, "CELF22B7");
    assert_eq!(results[0].trna_num, 1);
    assert_eq!(results[0].begin, 12619);
    assert_eq!(results[0].end, 12738);
    // ... etc
}
```

---

## Common Formatting Errors to Avoid

1. **Wrong decimal places:** HMM score as `51.2` instead of `51.20`
2. **0-indexed positions:** Using 0-indexed when 1-indexed required
3. **Missing tabs:** Using spaces instead of tabs
4. **Wrong strand logic:** Not reversing coordinates for reverse strand
5. **BED coordinate errors:** Not adjusting to 0-indexed
6. **Score precision:** Truncating instead of rounding
7. **Missing trailing commas:** BED block fields need trailing commas
8. **Case errors:** Anticodon in wrong case in sequence line
9. **Structure length mismatch:** Structure notation not same length as sequence
10. **Header formatting:** Not preserving exact header format

---

## Summary

This validation document ensures that the Rust implementation produces output that is:
- **Bit-identical** to the original tRNAscan-SE output
- **Format-compliant** with all standard specifications
- **Test-verified** against golden reference files
- **Edge-case-robust** for all known special cases

All output formatting must pass validation against these specifications before considering Phase 9 complete.
