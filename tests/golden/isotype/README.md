# Phase 8: Isotype Scoring Golden Files

## Overview

Phase 8 involves scoring tRNA candidates against **isotype-specific Covariance Models (CMs)** to determine the best-matching amino acid type. This is done using Infernal's `cmsearch` with pre-calibrated CM models.

## What is Isotype Scoring?

After initial tRNA detection (Phase 1-7), tRNAscan-SE performs a second round of scoring using **isotype-specific models**:

1. **General tRNA model** (Phase 7) identifies the sequence as a tRNA
2. **Isotype-specific models** (Phase 8) determine which amino acid it carries

Each tRNA is scored against ~23 different isotype models (Ala, Arg, Asn, Asp, Cys, Gln, Glu, Gly, His, Ile, Leu, Lys, Met, Phe, Pro, SeC, Ser, Thr, Trp, Tyr, Val, iMet, fMet)

## CM Model Locations

For bacterial mode (`-B`):
- **Main isotype models**: `lib/models/TRNAinf-bact-iso` (concatenated CM database)
- **Individual models**: Named like `bact-Ala`, `bact-Arg`, etc.
- **Selenocysteine model**: `lib/models/TRNAinf-bact-SeC.cm`

Model format: **Infernal 1.1** CM format

## Scoring Process

1. Run `cmsearch` with each isotype CM against the tRNA sequence
2. Extract the **bit score** for each model
3. The **highest-scoring model** determines the isotype
4. Score `-999` means the model could not score the sequence (structural mismatch)

## Output Format

### isotype_assignments.txt
Maps sequence ID to the final isotype determination.

### isotype_scores.txt
Full score matrix showing all isotype model scores for each tRNA.

### anticodon_table.txt
Reference table mapping anticodons to expected isotypes based on the genetic code.

## Example Sequences

From Example1.fa and Example2.fa:

1. **CELF22B7.trna1**: Leu (AAG anticodon) - Score 119.9
2. **CELF22B7.trna2**: Ser (CGA anticodon) - Score 125.0
3. **CELF22B7.trna3/4**: Phe (GAA anticodon) - Score 112.1
4. **CELF22B7.trna5**: Pro (TGG anticodon) - Score 113.0
5. **MySeq1.trna1**: Thr (TGT anticodon) - Score 93.1
6. **MySeq2.trna1**: Arg (TCT anticodon) - Score 89.3
7. **MySeq3.trna1**: Ser (CGA anticodon) - Score 118.3
8. **MySeq4.trna1**: Leu (AAG anticodon) - Score 92.2
9. **MySeq5.trna1**: SeC (TCA anticodon) - Score 146.9 (Selenocysteine!)
10. **MySeq6.trna1**: Lys (CTT anticodon) - Score 1.8 (but top score for Lys)

## Implementation Notes

### Rust Implementation Strategy

1. **Call Infernal's cmsearch** via `std::process::Command`
   - Path to CM database from config
   - Parse tabular output format (`--tblout`)

2. **Score Parsing**
   - Extract bit scores from Infernal output
   - Handle `-999` for invalid models
   - Track the top-scoring model

3. **Isotype Assignment**
   - Select model with highest bit score
   - Apply any tie-breaking rules (e.g., Met vs iMet)
   - Handle special cases like SeC (Selenocysteine)

4. **Validation**
   - Compare assigned isotype with anticodon prediction
   - Flag mismatches (e.g., anticodon says Ala but model says Leu)
   - Report confidence based on score differences

## References

- Infernal User Guide: http://eddylab.org/infernal/
- tRNAscan-SE Paper: doi.org/10.1093/nar/gkab688
- Covariance Models: doi.org/10.1093/bioinformatics/10.3.269
