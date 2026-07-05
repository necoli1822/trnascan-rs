# Phase 8 Implementation Guide: Isotype Scoring

## Overview

Phase 8 determines the **amino acid specificity** (isotype) of each tRNA candidate by scoring it against 23 isotype-specific Covariance Models using Infernal's `cmsearch`.

## Key Components

### 1. Input
- tRNA candidates from Phase 7 (CM screening)
- Isotype-specific CM database: `TRNAinf-bact-iso`

### 2. Process
```
For each tRNA candidate:
  1. Extract sequence
  2. Run cmsearch against all 23 isotype models
  3. Parse bit scores for each model
  4. Select highest-scoring model as predicted isotype
  5. Calculate confidence (score difference from runner-up)
```

### 3. Output
- Isotype assignment (e.g., "Leu", "Ser", "SeC")
- Score matrix showing all 23 model scores
- Confidence metric

## Rust Implementation Strategy

### Dependencies
```toml
[dependencies]
# For running cmsearch
subprocess = "0.2"

# For parsing tabular output
csv = "1.2"

# For score comparison
ordered-float = "3.0"
```

### Core Structure

```rust
// src/isotype/mod.rs

pub struct IsotypeScorer {
    cm_database_path: PathBuf,
    temp_dir: TempDir,
}

pub struct IsotypeScore {
    pub isotype: String,
    pub score: f64,
}

pub struct IsotypeResult {
    pub trna_id: String,
    pub predicted_isotype: String,
    pub top_score: f64,
    pub runner_up_isotype: String,
    pub runner_up_score: f64,
    pub scores: HashMap<String, f64>,  // All 23 scores
    pub confidence: Confidence,
}

pub enum Confidence {
    High,      // Score diff > 30 bits
    Medium,    // Score diff 10-30 bits
    Low,       // Score diff < 10 bits
    VeryLow,   // Runner-up higher (negative diff)
}

impl IsotypeScorer {
    pub fn new(cm_database_path: PathBuf) -> Result<Self>;

    pub fn score_trna(&self, trna: &Trna) -> Result<IsotypeResult>;

    fn run_cmsearch(&self, sequence: &str) -> Result<Vec<IsotypeScore>>;

    fn parse_cmsearch_output(&self, output: &str) -> Result<Vec<IsotypeScore>>;

    fn select_top_isotype(&self, scores: Vec<IsotypeScore>) -> IsotypeResult;
}
```

### Step 1: Run cmsearch

```rust
fn run_cmsearch(&self, sequence: &str) -> Result<Vec<IsotypeScore>> {
    // Write sequence to temporary FASTA file
    let temp_fa = self.temp_dir.path().join("candidate.fa");
    write_fasta(&temp_fa, sequence)?;

    // Run cmsearch
    let output_tbl = self.temp_dir.path().join("output.tbl");
    let status = Command::new("cmsearch")
        .arg("--tblout").arg(&output_tbl)
        .arg("--noali")          // Don't output alignments
        .arg("--cpu").arg("1")   // Single CPU
        .arg(&self.cm_database_path)
        .arg(&temp_fa)
        .output()?;

    if !status.status.success() {
        return Err(Error::CmsearchFailed);
    }

    // Parse output
    self.parse_cmsearch_output(&read_to_string(&output_tbl)?)
}
```

### Step 2: Parse cmsearch Output

```rust
fn parse_cmsearch_output(&self, output: &str) -> Result<Vec<IsotypeScore>> {
    let mut scores = Vec::new();

    for line in output.lines() {
        // Skip comments
        if line.starts_with('#') {
            continue;
        }

        // Parse tab-delimited output
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 15 {
            continue;
        }

        // Extract fields:
        // [2] = query name (model name, e.g., "bact-Leu")
        // [14] = bit score
        let model_name = fields[2];
        let isotype = model_name.strip_prefix("bact-")
            .unwrap_or(model_name)
            .to_string();

        let score = fields[14].parse::<f64>()?;

        scores.push(IsotypeScore { isotype, score });
    }

    Ok(scores)
}
```

### Step 3: Select Top Isotype

```rust
fn select_top_isotype(&self, mut scores: Vec<IsotypeScore>) -> IsotypeResult {
    // Sort by score (descending)
    scores.sort_by(|a, b| {
        b.score.partial_cmp(&a.score).unwrap_or(Ordering::Equal)
    });

    // Handle -999 scores (invalid)
    let valid_scores: Vec<_> = scores.iter()
        .filter(|s| s.score > -900.0)
        .collect();

    if valid_scores.is_empty() {
        return IsotypeResult::invalid();
    }

    let top = &valid_scores[0];
    let runner_up = valid_scores.get(1);

    let score_diff = if let Some(ru) = runner_up {
        top.score - ru.score
    } else {
        top.score  // Only one valid score
    };

    let confidence = match score_diff {
        d if d > 30.0 => Confidence::High,
        d if d > 10.0 => Confidence::Medium,
        d if d > 0.0 => Confidence::Low,
        _ => Confidence::VeryLow,
    };

    IsotypeResult {
        predicted_isotype: top.isotype.clone(),
        top_score: top.score,
        runner_up_isotype: runner_up.map(|s| s.isotype.clone()).unwrap_or_default(),
        runner_up_score: runner_up.map(|s| s.score).unwrap_or(-999.0),
        scores: scores.into_iter().map(|s| (s.isotype, s.score)).collect(),
        confidence,
    }
}
```

### Step 4: Validation

```rust
fn validate_with_anticodon(&self, result: &IsotypeResult, anticodon: &str)
    -> ValidationResult
{
    let expected_isotype = ANTICODON_TABLE.get(anticodon);

    match expected_isotype {
        Some(exp) if exp == &result.predicted_isotype => {
            ValidationResult::Match
        },
        Some(exp) => {
            ValidationResult::Mismatch {
                expected: exp.clone(),
                predicted: result.predicted_isotype.clone(),
            }
        },
        None => ValidationResult::UnknownAnticodon,
    }
}
```

## Testing Strategy

### Unit Tests
```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_cmsearch_output() {
        let output = include_str!("../../tests/data/cmsearch_output.tbl");
        let scores = parse_cmsearch_output(output).unwrap();
        assert_eq!(scores.len(), 23);

        // Find Leu score
        let leu = scores.iter().find(|s| s.isotype == "Leu").unwrap();
        assert!((leu.score - 119.9).abs() < 0.1);
    }

    #[test]
    fn test_select_top_isotype_leu() {
        let scores = vec![
            IsotypeScore { isotype: "Leu".into(), score: 119.9 },
            IsotypeScore { isotype: "Ser".into(), score: 69.0 },
            IsotypeScore { isotype: "Arg".into(), score: 55.3 },
        ];

        let result = select_top_isotype(scores);
        assert_eq!(result.predicted_isotype, "Leu");
        assert_eq!(result.confidence, Confidence::High);
    }

    #[test]
    fn test_sec_detection() {
        // MySeq5.trna1 should be SeC with very high score
        let scores = vec![
            IsotypeScore { isotype: "SeC".into(), score: 146.9 },
            IsotypeScore { isotype: "Asn".into(), score: 5.4 },
            // ... rest
        ];

        let result = select_top_isotype(scores);
        assert_eq!(result.predicted_isotype, "SeC");
        assert!(result.top_score > 140.0);
    }

    #[test]
    fn test_low_confidence_lys() {
        // MySeq6.trna1 - Lys with very low score
        let scores = vec![
            IsotypeScore { isotype: "Ser".into(), score: 72.4 },
            IsotypeScore { isotype: "Leu".into(), score: 75.7 },
            IsotypeScore { isotype: "Lys".into(), score: 1.8 },
        ];

        let result = select_top_isotype(scores);
        // Top scorer is Leu, not Lys!
        assert_eq!(result.predicted_isotype, "Leu");

        // This reveals a bug in the original test case!
        // Need to check how tRNAscan-SE actually handles this.
    }
}
```

### Integration Tests
```rust
#[test]
fn test_example1_isotype_scoring() {
    let scorer = IsotypeScorer::new("original/lib/models/TRNAinf-bact-iso".into()).unwrap();

    // Load tRNAs from Example1
    let trnas = load_test_trnas("tests/data/Example1.fa");

    // Score each
    for trna in trnas {
        let result = scorer.score_trna(&trna).unwrap();

        // Load expected from golden file
        let expected = load_expected_isotype(&trna.id);

        assert_eq!(result.predicted_isotype, expected.isotype);
        assert!((result.top_score - expected.score).abs() < 0.1);
    }
}
```

## Output Formats

### 1. Score Matrix (TSV)
```
tRNAscanID	Anticodon_predicted_isotype	Ala	Arg	...	Val	iMet
CELF22B7.trna1	Leu	37.4	55.3	...	28.8	-999
```

### 2. Assignment Report (Human-readable)
```
CELF22B7.trna1: Leu (score=119.9, confidence=HIGH, anticodon=AAG)
  Top scores: Leu:119.9, Ser:69.0, Arg:55.3
  Score difference: 50.9 bits
  Validation: ✓ Anticodon matches expected isotype
```

### 3. JSON Output (Machine-readable)
```json
{
  "trna_id": "CELF22B7.trna1",
  "predicted_isotype": "Leu",
  "top_score": 119.9,
  "runner_up": "Ser",
  "runner_up_score": 69.0,
  "confidence": "High",
  "scores": {
    "Ala": 37.4,
    "Arg": 55.3,
    ...
  },
  "anticodon_validation": "match"
}
```

## Error Handling

```rust
#[derive(Debug, thiserror::Error)]
pub enum IsotypeError {
    #[error("cmsearch failed: {0}")]
    CmsearchFailed(String),

    #[error("CM database not found: {0}")]
    DatabaseNotFound(PathBuf),

    #[error("No valid isotype scores")]
    NoValidScores,

    #[error("Failed to parse cmsearch output: {0}")]
    ParseError(String),
}
```

## Performance Optimization

1. **Parallel Processing**: Score multiple tRNAs in parallel
2. **CM Database Caching**: Keep database in memory if possible
3. **Batch Processing**: Process multiple sequences in one cmsearch call

## Special Cases

### 1. Selenocysteine (SeC)
- Use separate `TRNAinf-bact-SeC.cm` model
- Very high scores (>140 bits) typical
- Anticodon TCA, suppresses UGA stop codon

### 2. Initiator Methionine (iMet/fMet)
- Distinguish from elongator Met
- Check acceptor stem mismatches
- Bacterial: May be formylated (fMet)

### 3. Type II tRNAs (Ser, Leu)
- Long variable arm (13-21 bp)
- Higher structural complexity
- May cross-score with each other

### 4. Low-Confidence Assignments
- Flag if score difference < 10 bits
- Flag if top score < 30 bits
- Suggest manual review

## References

- **Infernal**: http://eddylab.org/infernal/
- **tRNAscan-SE**: http://trna.ucsc.edu/
- **CM Theory**: Eddy & Durbin (1994) NAR 22:2079
