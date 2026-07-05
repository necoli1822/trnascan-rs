# 🧬 Phase 8: Isotype Scoring - START HERE

Welcome to the Phase 8 (Isotype Scoring) golden files for tRNAscan-SE Rust reimplementation!

## 📚 Quick Navigation

### New to Isotype Scoring?
👉 **Start with**: [`README.md`](README.md)
- What is isotype scoring?
- How does it work?
- Why is it important?

### Ready to Implement?
👉 **Follow**: [`IMPLEMENTATION_GUIDE.md`](IMPLEMENTATION_GUIDE.md)
- Complete Rust implementation strategy
- Code examples and patterns
- Error handling and optimization

### Need Reference Data?
👉 **Check**:
- [`anticodon_table.txt`](anticodon_table.txt) - Genetic code mapping
- [`cm_model_info.txt`](cm_model_info.txt) - CM model specifications

### Writing Tests?
👉 **Use**:
- [`test_cases.md`](test_cases.md) - Detailed test scenarios
- [`isotype_scores.txt`](isotype_scores.txt) - Golden score matrix
- [`isotype_assignments.txt`](isotype_assignments.txt) - Expected assignments

### Looking for Overview?
👉 **See**: [`INDEX.md`](INDEX.md)
- Complete file index
- Implementation checklist
- Validation criteria

---

## ⚡ Quick Start (2 Minutes)

### What You'll Build
An **isotype scorer** that:
1. Takes a tRNA sequence
2. Scores it against 23 isotype-specific CM models using Infernal
3. Selects the best-matching amino acid type
4. Reports confidence and validates against anticodon

### Example
```bash
Input:  CELF22B7.trna1 (AAG anticodon)
Output: Leucine (Leu), score=119.9, confidence=HIGH
```

### Test Data
- **11 tRNAs** from Example1.fa and Example2.fa
- **253 scores** (11 tRNAs × 23 isotypes)
- **7 special cases** including SeC (Selenocysteine)

---

## 🎯 Implementation Steps

1. **Invoke Infernal's cmsearch** with CM database
2. **Parse tabular output** to extract bit scores
3. **Select top-scoring model** as predicted isotype
4. **Calculate confidence** from score differences
5. **Validate** against anticodon prediction

---

## ✅ Success Criteria

Your implementation passes if:
- ✓ All 11 test tRNAs match expected isotypes
- ✓ Scores match within ±0.1 bits
- ✓ SeC (Selenocysteine) detected correctly
- ✓ Low-confidence cases flagged
- ✓ Outputs match golden files

---

## 📦 What's Included

| File | Size | Purpose |
|------|------|---------|
| 00-START-HERE.md | You are here | Navigation guide |
| INDEX.md | 5.4K | Master index |
| README.md | 3.1K | Introduction |
| IMPLEMENTATION_GUIDE.md | 9.9K | **Main implementation guide** |
| test_cases.md | 6.5K | Test scenarios |
| isotype_scores.txt | 2.6K | **Golden score matrix** |
| isotype_assignments.txt | 1.0K | **Expected results** |
| anticodon_table.txt | 2.2K | Genetic code reference |
| cm_model_info.txt | 4.0K | Model specifications |
| VERIFICATION.txt | - | Completion report |

**Total**: 10 files, ~40KB documentation

---

## 🔬 Key Concepts

### What is an Isotype?
The **amino acid specificity** of a tRNA. Examples: Leucine (Leu), Serine (Ser), Phenylalanine (Phe).

### What is a CM?
A **Covariance Model** - a probabilistic model that captures both sequence and structure of RNA families.

### What is cmsearch?
Infernal's tool to search sequences against CM databases. Like BLAST but for RNA structure.

### What is a bit score?
A log-odds score indicating how well a sequence matches a model:
- **>100 bits**: Excellent match
- **60-100**: Good match  
- **30-60**: Marginal match
- **<30**: Poor match
- **-999**: Cannot score (structure mismatch)

---

## 🚀 Next Steps

1. **Read** [`README.md`](README.md) (5 min)
2. **Study** [`IMPLEMENTATION_GUIDE.md`](IMPLEMENTATION_GUIDE.md) (30 min)
3. **Implement** IsotypeScorer struct (2-3 days)
4. **Test** against golden files (1 day)

---

## 💡 Pro Tips

- **Use parallel processing** for multiple tRNAs
- **Cache CM database** in memory if possible
- **Flag low-confidence** assignments for review
- **Handle SeC specially** - it's the 21st amino acid!

---

## 📞 Need Help?

- **Infernal Manual**: http://eddylab.org/infernal/Userguide.pdf
- **tRNAscan-SE Paper**: doi.org/10.1093/nar/gkab688
- **CM Theory**: Eddy & Durbin (1994) NAR 22:2079

---

## ✨ Fun Facts

- **23 isotypes**: 20 standard amino acids + SeC + iMet + fMet
- **SeC = 21st amino acid**: Selenocysteine, inserted at UGA stop codons
- **Type II tRNAs**: Ser and Leu have long variable arms (13-21 bp)
- **Highest score**: 146.9 bits (MySeq5.trna1 SeC)
- **Most unusual**: MySeq6.trna1 Lys (score=1.8, but correct!)

---

**Happy Coding! 🦀**

*Generated 2026-03-04 for tRNAscan-SE Rust Reimplementation*
