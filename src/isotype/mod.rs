//! Anticodon → isotype mapping for tRNAscan-SE.
//!
//! Maps anticodons to amino-acid isotypes using the standard genetic code and
//! handles the special isotypes (SeC, iMet/fMet, Sup). The actual per-isotype
//! covariance-model scoring lives in `core::scanner` (in-process infernox
//! cmscan), not here.
//!
//! # Example Usage
//!
//! ```
//! use trnascan_rs::isotype::anticodon_to_isotype;
//!
//! let isotype = anticodon_to_isotype("AAG");
//! assert_eq!(isotype, Some("Leu"));
//! ```

pub mod anticodon;

// Re-export commonly used items
pub use anticodon::{
    anticodon_to_isotype,
    get_all_anticodons_for_isotype,
    is_selenocysteine,
    is_suppressor,
    has_long_variable_arm,
};
