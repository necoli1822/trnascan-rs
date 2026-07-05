//! Anticodon to isotype mapping based on the standard genetic code.
//!
//! This module provides functions to map 3-letter anticodons (written 3'->5')
//! to their corresponding amino acid isotypes based on the standard genetic code.

use std::collections::HashMap;
use once_cell::sync::Lazy;

/// Static anticodon to isotype mapping table.
///
/// Anticodons are written 3'->5' as they appear in the tRNA sequence.
/// The actual codon they recognize is the reverse complement.
static ANTICODON_MAP: Lazy<HashMap<&'static str, &'static str>> = Lazy::new(|| {
    let mut map = HashMap::new();

    // Alanine (Ala)
    map.insert("AGC", "Ala");
    map.insert("CGC", "Ala");
    map.insert("TGC", "Ala");
    map.insert("GGC", "Ala");

    // Arginine (Arg)
    map.insert("ACG", "Arg");
    map.insert("CCG", "Arg");
    map.insert("TCG", "Arg");
    map.insert("GCG", "Arg");
    map.insert("CCT", "Arg");
    map.insert("TCT", "Arg"); // Wobble

    // Asparagine (Asn)
    map.insert("GTT", "Asn");

    // Aspartate (Asp)
    map.insert("GTC", "Asp");

    // Cysteine (Cys)
    map.insert("GCA", "Cys");

    // Glutamine (Gln)
    map.insert("CTG", "Gln");
    map.insert("TTG", "Gln");

    // Glutamate (Glu)
    map.insert("CTC", "Glu");
    map.insert("TTC", "Glu");

    // Glycine (Gly)
    map.insert("ACC", "Gly");
    map.insert("GCC", "Gly");
    map.insert("TCC", "Gly");
    map.insert("CCC", "Gly");

    // Histidine (His)
    map.insert("GTG", "His");

    // Isoleucine (Ile)
    map.insert("AAT", "Ile");
    map.insert("GAT", "Ile");
    map.insert("TAT", "Ile");

    // Leucine (Leu)
    map.insert("CAA", "Leu");
    map.insert("TAA", "Leu");
    map.insert("CAG", "Leu");
    map.insert("TAG", "Leu");
    map.insert("AAG", "Leu"); // Wobble
    map.insert("GAG", "Leu");

    // Lysine (Lys)
    map.insert("CTT", "Lys");
    map.insert("TTT", "Lys");

    // Methionine (Met) / Initiator Methionine (iMet/fMet)
    // Note: CAT can be Met, iMet, or fMet - determined by structure
    map.insert("CAT", "Met");

    // Phenylalanine (Phe)
    map.insert("GAA", "Phe");

    // Proline (Pro)
    map.insert("AGG", "Pro");
    map.insert("CGG", "Pro");
    map.insert("TGG", "Pro");
    map.insert("GGG", "Pro");

    // Selenocysteine (SeC) - special case, also a Sup
    // Note: TCA can be SeC or Sup depending on context
    map.insert("TCA", "SeC");

    // Serine (Ser)
    map.insert("ACT", "Ser");
    map.insert("GCT", "Ser");
    map.insert("TGA", "Ser");
    map.insert("CGA", "Ser");
    map.insert("AGA", "Ser");
    map.insert("GGA", "Ser");

    // Threonine (Thr)
    map.insert("AGT", "Thr");
    map.insert("CGT", "Thr");
    map.insert("TGT", "Thr");
    map.insert("GGT", "Thr");

    // Tryptophan (Trp)
    map.insert("CCA", "Trp");

    // Tyrosine (Tyr)
    map.insert("GTA", "Tyr");

    // Valine (Val)
    map.insert("AAC", "Val");
    map.insert("CAC", "Val");
    map.insert("TAC", "Val");
    map.insert("GAC", "Val");

    // Suppressor tRNAs (stop codon suppressors)
    // Note: These may override other assignments
    // TTA -> UAA (Ochre suppressor)
    // CTA -> UAG (Amber suppressor)
    // TCA -> UGA (Opal suppressor or SeC)
    map.insert("TTA", "Sup");
    map.insert("CTA", "Sup");
    // TCA already mapped to SeC above

    map
});

/// Reverse mapping: isotype to all possible anticodons.
static ISOTYPE_TO_ANTICODONS: Lazy<HashMap<&'static str, Vec<&'static str>>> = Lazy::new(|| {
    let mut map: HashMap<&'static str, Vec<&'static str>> = HashMap::new();

    for (&anticodon, &isotype) in ANTICODON_MAP.iter() {
        map.entry(isotype).or_insert_with(Vec::new).push(anticodon);
    }

    // Add special cases that don't have standard anticodons in the main map
    map.insert("iMet", vec!["CAT"]);
    map.insert("fMet", vec!["CAT"]);

    map
});

/// Maps a 3-letter anticodon (3'->5') to its corresponding amino acid isotype.
///
/// # Arguments
/// * `anticodon` - The 3-letter anticodon sequence (e.g., "TGT", "CGA")
///
/// # Returns
/// * `Some(&str)` - The isotype name (e.g., "Thr", "Ser", "SeC")
/// * `None` - If the anticodon is not recognized
///
/// # Examples
/// ```
/// use trnascan_rs::isotype::anticodon_to_isotype;
/// assert_eq!(anticodon_to_isotype("TGT"), Some("Thr"));
/// assert_eq!(anticodon_to_isotype("CGA"), Some("Ser"));
/// assert_eq!(anticodon_to_isotype("TCA"), Some("SeC"));
/// assert_eq!(anticodon_to_isotype("XYZ"), None);
/// ```
pub fn anticodon_to_isotype(anticodon: &str) -> Option<&'static str> {
    ANTICODON_MAP.get(anticodon).copied()
}

/// Returns all anticodons that decode to the specified isotype.
///
/// # Arguments
/// * `isotype` - The amino acid isotype (e.g., "Leu", "Ser", "Met")
///
/// # Returns
/// * A vector of all anticodon sequences for this isotype
/// * An empty vector if the isotype is not recognized
///
/// # Examples
/// ```
/// use trnascan_rs::isotype::get_all_anticodons_for_isotype;
/// let leu_anticodons = get_all_anticodons_for_isotype("Leu");
/// assert!(leu_anticodons.contains(&"AAG"));
/// assert!(leu_anticodons.contains(&"CAG"));
/// ```
pub fn get_all_anticodons_for_isotype(isotype: &str) -> Vec<&'static str> {
    ISOTYPE_TO_ANTICODONS
        .get(isotype)
        .map(|v: &Vec<&'static str>| v.clone())
        .unwrap_or_default()
}

/// Checks if an anticodon is a known suppressor tRNA.
///
/// Suppressor tRNAs recognize stop codons (UAA, UAG, UGA).
///
/// # Arguments
/// * `anticodon` - The 3-letter anticodon sequence
///
/// # Returns
/// * `true` if the anticodon corresponds to a suppressor tRNA
pub fn is_suppressor(anticodon: &str) -> bool {
    matches!(anticodon, "TTA" | "CTA" | "TCA")
}

/// Checks if an anticodon corresponds to Selenocysteine (SeC).
///
/// # Arguments
/// * `anticodon` - The 3-letter anticodon sequence
///
/// # Returns
/// * `true` if the anticodon is TCA (SeC/Opal suppressor)
pub fn is_selenocysteine(anticodon: &str) -> bool {
    anticodon == "TCA"
}

/// Checks if an isotype has a long variable arm (Type II tRNA).
///
/// Serine and Leucine tRNAs have long variable arms (13-21 bp).
///
/// # Arguments
/// * `isotype` - The amino acid isotype
///
/// # Returns
/// * `true` if the isotype has a long variable arm
pub fn has_long_variable_arm(isotype: &str) -> bool {
    matches!(isotype, "Ser" | "Leu")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_anticodon_mapping() {
        assert_eq!(anticodon_to_isotype("TGT"), Some("Thr"));
        assert_eq!(anticodon_to_isotype("CGA"), Some("Ser"));
        assert_eq!(anticodon_to_isotype("AAG"), Some("Leu"));
        assert_eq!(anticodon_to_isotype("GAA"), Some("Phe"));
        assert_eq!(anticodon_to_isotype("CTT"), Some("Lys"));
    }

    #[test]
    fn test_special_isotypes() {
        assert_eq!(anticodon_to_isotype("TCA"), Some("SeC"));
        assert_eq!(anticodon_to_isotype("CAT"), Some("Met"));
        assert_eq!(anticodon_to_isotype("TTA"), Some("Sup"));
        assert_eq!(anticodon_to_isotype("CTA"), Some("Sup"));
    }

    #[test]
    fn test_unknown_anticodon() {
        assert_eq!(anticodon_to_isotype("XYZ"), None);
        assert_eq!(anticodon_to_isotype(""), None);
    }

    #[test]
    fn test_get_anticodons_for_isotype() {
        let leu_anticodons = get_all_anticodons_for_isotype("Leu");
        assert_eq!(leu_anticodons.len(), 6);
        assert!(leu_anticodons.contains(&"AAG"));
        assert!(leu_anticodons.contains(&"CAG"));
        assert!(leu_anticodons.contains(&"TAG"));

        let phe_anticodons = get_all_anticodons_for_isotype("Phe");
        assert_eq!(phe_anticodons.len(), 1);
        assert!(phe_anticodons.contains(&"GAA"));
    }

    #[test]
    fn test_suppressor_detection() {
        assert!(is_suppressor("TTA"));
        assert!(is_suppressor("CTA"));
        assert!(is_suppressor("TCA"));
        assert!(!is_suppressor("AAG"));
    }

    #[test]
    fn test_selenocysteine_detection() {
        assert!(is_selenocysteine("TCA"));
        assert!(!is_selenocysteine("CTA"));
    }

    #[test]
    fn test_variable_arm_detection() {
        assert!(has_long_variable_arm("Ser"));
        assert!(has_long_variable_arm("Leu"));
        assert!(!has_long_variable_arm("Ala"));
        assert!(!has_long_variable_arm("Met"));
    }
}
