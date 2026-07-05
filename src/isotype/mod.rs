//! Isotype scoring module for tRNAscan-SE.
//!
//! This module provides functionality for:
//! - Mapping anticodons to amino acid isotypes using the standard genetic code
//! - Scoring tRNAs against covariance models for each isotype
//! - Detecting isotype mismatches between anticodon prediction and CM scoring
//! - Handling special isotypes (SeC, iMet/fMet, Sup)
//!
//! # Example Usage
//!
//! ```
//! use trnascan_rs::isotype::{anticodon_to_isotype, IsotypeScorer};
//!
//! // Basic anticodon mapping
//! let isotype = anticodon_to_isotype("AAG");
//! assert_eq!(isotype, Some("Leu"));
//!
//! // Full isotype scoring (requires TrnaInfo and CM scores)
//! // let scorer = IsotypeScorer::score_isotype(&trna_info, "AAG", &cm_scores);
//! // println!("Predicted: {}", scorer.predicted_isotype);
//! // println!("CM best: {}", scorer.cm_best_isotype);
//! ```

pub mod anticodon;
pub mod scorer;

// Re-export commonly used items
pub use anticodon::{
    anticodon_to_isotype,
    get_all_anticodons_for_isotype,
    is_selenocysteine,
    is_suppressor,
    has_long_variable_arm,
};

pub use scorer::{
    IsotypeScorer,
    check_isotype_mismatch,
    get_isotype_thresholds,
    score_all_isotypes,
};
