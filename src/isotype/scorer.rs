//! Isotype scoring logic for tRNA identification.
//!
//! This module provides the core isotype scoring functionality, including:
//! - Scoring tRNAs against covariance models for each isotype
//! - Determining the best isotype match and confidence scores
//! - Detecting isotype mismatches between anticodon prediction and CM scores

use crate::eufind::pavesi::TrnaInfo;
use super::anticodon::anticodon_to_isotype;
use std::collections::HashMap;

/// Isotype scorer structure containing prediction results.
#[derive(Debug, Clone)]
pub struct IsotypeScorer {
    /// The anticodon sequence (3'->5')
    pub anticodon: String,

    /// The predicted isotype based on anticodon
    pub predicted_isotype: String,

    /// The overall isotype confidence score
    pub isotype_score: f32,

    /// The HMM/CM bit score for the best matching isotype model
    pub hmm_score: f32,

    /// Secondary structure match score component
    pub secondary_structure_score: f32,

    /// The best matching isotype from CM scoring (may differ from predicted)
    pub cm_best_isotype: String,

    /// The CM bit score for the best isotype
    pub cm_best_score: f32,

    /// The second-best isotype from CM scoring
    pub cm_runner_up_isotype: String,

    /// The CM bit score for the runner-up isotype
    pub cm_runner_up_score: f32,

    /// Score difference between best and runner-up
    pub score_difference: f32,

    /// All isotype CM scores (isotype -> bit score)
    pub all_cm_scores: HashMap<String, f32>,
}

impl IsotypeScorer {
    /// Score a tRNA against all isotype covariance models.
    ///
    /// This determines the most likely isotype based on:
    /// 1. Anticodon sequence (genetic code mapping)
    /// 2. Structural fit to CM models for each isotype
    /// 3. Confidence based on score separation
    ///
    /// # Arguments
    /// * `trna_info` - The tRNA information structure
    /// * `anticodon` - The 3-letter anticodon sequence
    /// * `cm_scores` - Map of isotype names to CM bit scores
    ///
    /// # Returns
    /// An `IsotypeScorer` with all scoring results
    pub fn score_isotype(
        trna_info: &TrnaInfo,
        anticodon: &str,
        cm_scores: &HashMap<String, f32>,
    ) -> Self {
        // Determine predicted isotype from anticodon
        let predicted_isotype = anticodon_to_isotype(anticodon)
            .unwrap_or("Unk")
            .to_string();

        // Special handling for initiator methionine
        let predicted_isotype = if predicted_isotype == "Met" {
            // Check if this is an initiator Met based on structure
            // (in practice, this would examine D-loop and other structural features)
            // For now, default to Met unless context suggests otherwise
            predicted_isotype
        } else {
            predicted_isotype
        };

        // Find best and runner-up isotypes from CM scores
        let mut sorted_scores: Vec<_> = cm_scores.iter().collect();
        sorted_scores.sort_by(|a, b| b.1.partial_cmp(a.1).unwrap());

        let (cm_best_isotype, cm_best_score) = if !sorted_scores.is_empty() {
            (sorted_scores[0].0.clone(), *sorted_scores[0].1)
        } else {
            (predicted_isotype.clone(), -999.0)
        };

        let (cm_runner_up_isotype, cm_runner_up_score) = if sorted_scores.len() > 1 {
            (sorted_scores[1].0.clone(), *sorted_scores[1].1)
        } else {
            ("None".to_string(), -999.0)
        };

        let score_difference = cm_best_score - cm_runner_up_score;

        // Calculate overall isotype score
        // This combines the CM score with anticodon concordance
        let isotype_score = if cm_best_isotype == predicted_isotype {
            // Bonus for agreement between anticodon and CM
            cm_best_score + 5.0
        } else {
            // Penalty for disagreement
            cm_best_score - 10.0
        };

        // Use the total score from trna_info as HMM score
        let hmm_score = trna_info.tot_sc;

        // Secondary structure score (simplified - use A-box and B-box scores)
        let secondary_structure_score = trna_info.abox_sc + trna_info.bbox_sc;

        IsotypeScorer {
            anticodon: anticodon.to_string(),
            predicted_isotype,
            isotype_score,
            hmm_score,
            secondary_structure_score,
            cm_best_isotype,
            cm_best_score,
            cm_runner_up_isotype,
            cm_runner_up_score,
            score_difference,
            all_cm_scores: cm_scores.clone(),
        }
    }

    /// Create a default scorer for cases where CM scoring is unavailable.
    pub fn default_scorer(anticodon: &str) -> Self {
        let predicted_isotype = anticodon_to_isotype(anticodon)
            .unwrap_or("Unk")
            .to_string();

        IsotypeScorer {
            anticodon: anticodon.to_string(),
            predicted_isotype: predicted_isotype.clone(),
            isotype_score: -999.0,
            hmm_score: -999.0,
            secondary_structure_score: -999.0,
            cm_best_isotype: predicted_isotype,
            cm_best_score: -999.0,
            cm_runner_up_isotype: "None".to_string(),
            cm_runner_up_score: -999.0,
            score_difference: 0.0,
            all_cm_scores: HashMap::new(),
        }
    }

    /// Check if this is a high-confidence isotype assignment.
    ///
    /// High confidence requires:
    /// - CM score > threshold (typically 40 bits)
    /// - Score difference > threshold (typically 20 bits)
    /// - Agreement between anticodon and CM best match
    pub fn is_high_confidence(&self, min_score: f32, min_difference: f32) -> bool {
        self.cm_best_score >= min_score
            && self.score_difference >= min_difference
            && self.cm_best_isotype == self.predicted_isotype
    }

    /// Get the final isotype call.
    ///
    /// This returns the most likely isotype, considering both anticodon
    /// prediction and CM scoring.
    pub fn final_isotype(&self) -> &str {
        // If CM scoring is available and has high confidence, use it
        if self.cm_best_score > 40.0 && self.score_difference > 15.0 {
            &self.cm_best_isotype
        } else {
            // Otherwise, trust the anticodon prediction
            &self.predicted_isotype
        }
    }

    /// Check if this tRNA is a potential pseudogene.
    ///
    /// Indicators include:
    /// - Very low CM scores (< 20 bits)
    /// - Negative score difference (runner-up better than prediction)
    /// - Poor secondary structure scores
    pub fn is_potential_pseudogene(&self) -> bool {
        self.cm_best_score < 20.0
            || self.score_difference < -30.0
            || self.secondary_structure_score < -50.0
    }
}

/// Check for isotype mismatch (ISM) between anticodon prediction and CM scoring.
///
/// An isotype mismatch is detected when:
/// 1. The anticodon-predicted isotype differs from the best CM match
/// 2. The score difference is significant (typically > 20 bits)
/// 3. Neither is a special case (SeC, Sup, iMet, etc.)
///
/// # Arguments
/// * `predicted_isotype` - The isotype predicted from the anticodon
/// * `cm_isotype` - The best-scoring isotype from CM analysis
/// * `score_diff` - The score difference between best and runner-up
///
/// # Returns
/// * `true` if a significant mismatch is detected
pub fn check_isotype_mismatch(
    predicted_isotype: &str,
    cm_isotype: &str,
    score_diff: f32,
) -> bool {
    // No mismatch if they agree
    if predicted_isotype == cm_isotype {
        return false;
    }

    // No mismatch if score difference is too small
    if score_diff < 20.0 {
        return false;
    }

    // Special cases that may legitimately differ
    match (predicted_isotype, cm_isotype) {
        // Met/iMet/fMet can be confused
        ("Met", "iMet") | ("iMet", "Met") | ("fMet", "Met") | ("Met", "fMet") => false,

        // SeC may be confused with Sup
        ("SeC", "Sup") | ("Sup", "SeC") => false,

        // Ser and Leu both have variable arms and may cross-score
        ("Ser", "Leu") | ("Leu", "Ser") if score_diff < 30.0 => false,

        // Otherwise, this is a genuine mismatch
        _ => true,
    }
}

/// Calculate isotype-specific scoring thresholds.
///
/// Different isotypes have different typical score ranges due to
/// structural differences (e.g., variable arm length).
///
/// # Arguments
/// * `isotype` - The amino acid isotype
///
/// # Returns
/// * `(min_score, min_difference)` - Recommended thresholds
pub fn get_isotype_thresholds(isotype: &str) -> (f32, f32) {
    match isotype {
        // SeC has very distinctive structure -> high thresholds
        "SeC" => (100.0, 80.0),

        // Ser and Leu have variable arms -> moderate thresholds
        "Ser" | "Leu" => (60.0, 25.0),

        // Standard tRNAs
        _ => (40.0, 20.0),
    }
}

/// Parse CM scores from a result string or file.
///
/// This is a placeholder for the actual CM scoring integration.
/// In production, this would interface with Infernal cmscan results.
///
/// # Arguments
/// * `isotype` - The isotype to score against
/// * `sequence` - The tRNA sequence
///
/// # Returns
/// * The CM bit score
pub fn score_cm_model(_isotype: &str, _sequence: &str) -> f32 {
    // Placeholder: In production, this would run cmscan
    // For now, return a dummy score
    -999.0
}

/// Score all isotype models for a tRNA sequence.
///
/// # Arguments
/// * `sequence` - The tRNA sequence
///
/// # Returns
/// * HashMap of isotype -> bit score
pub fn score_all_isotypes(sequence: &str) -> HashMap<String, f32> {
    let mut scores = HashMap::new();

    // List of all standard isotypes
    let isotypes = vec![
        "Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly",
        "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "SeC",
        "Ser", "Thr", "Trp", "Tyr", "Val", "iMet",
    ];

    for isotype in isotypes {
        scores.insert(isotype.to_string(), score_cm_model(isotype, sequence));
    }

    scores
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_isotype_mismatch_detection() {
        // Same isotype -> no mismatch
        assert!(!check_isotype_mismatch("Leu", "Leu", 50.0));

        // Different isotypes, high score -> mismatch
        assert!(check_isotype_mismatch("Leu", "Ser", 30.0));

        // Different isotypes, low score -> no mismatch (ambiguous)
        assert!(!check_isotype_mismatch("Leu", "Ser", 10.0));

        // Met/iMet -> no mismatch (related)
        assert!(!check_isotype_mismatch("Met", "iMet", 50.0));

        // SeC/Sup -> no mismatch (both suppress UGA)
        assert!(!check_isotype_mismatch("SeC", "Sup", 50.0));
    }

    #[test]
    fn test_isotype_thresholds() {
        let (sec_score, sec_diff) = get_isotype_thresholds("SeC");
        assert_eq!(sec_score, 100.0);
        assert_eq!(sec_diff, 80.0);

        let (ser_score, ser_diff) = get_isotype_thresholds("Ser");
        assert_eq!(ser_score, 60.0);
        assert_eq!(ser_diff, 25.0);

        let (ala_score, ala_diff) = get_isotype_thresholds("Ala");
        assert_eq!(ala_score, 40.0);
        assert_eq!(ala_diff, 20.0);
    }

    #[test]
    fn test_default_scorer() {
        let scorer = IsotypeScorer::default_scorer("AAG");
        assert_eq!(scorer.anticodon, "AAG");
        assert_eq!(scorer.predicted_isotype, "Leu");
        assert_eq!(scorer.cm_best_score, -999.0);
    }

    #[test]
    fn test_confidence_levels() {
        let mut scorer = IsotypeScorer::default_scorer("AAG");
        scorer.cm_best_score = 80.0;
        scorer.score_difference = 25.0;
        scorer.cm_best_isotype = "Leu".to_string();

        assert!(scorer.is_high_confidence(40.0, 20.0));
        assert_eq!(scorer.final_isotype(), "Leu");
    }

    #[test]
    fn test_pseudogene_detection() {
        let mut scorer = IsotypeScorer::default_scorer("AAG");
        scorer.cm_best_score = 15.0; // Very low
        scorer.score_difference = -40.0; // Negative
        scorer.secondary_structure_score = -60.0; // Poor

        assert!(scorer.is_potential_pseudogene());
    }
}
