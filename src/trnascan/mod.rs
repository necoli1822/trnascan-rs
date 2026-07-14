//! Sequence utilities shared by the tRNAscan-SE port.
//!
//! Historically this module also hosted a Rust port of the legacy tRNAscan-1.4
//! Fichant-Burks signal first-pass. That first-pass is never reached — every
//! supported mode (euk/bact/arch/organ/mito/general) runs the Infernal
//! covariance-model first-pass instead (tRNAscan-SE:1514-1531 `infernal_fp=1`) —
//! so the dead scaffolding was removed. Only the live sequence helpers remain.

pub mod seq_utils;

pub use seq_utils::{anticodon_to_aa, encode_anticodon, reverse_complement};
