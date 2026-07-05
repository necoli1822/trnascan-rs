//! Genetic code module for tRNAscan-SE
//!
//! This module provides genetic code translation tables and tRNA-related
//! functions for converting between codons, anticodons, and amino acids.
//!
//! Ported from tRNAscanSE::GeneticCode.pm

use std::collections::HashMap;

/// Genetic code error types
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GeneticCodeError {
    InvalidCodeId(i32),
    InvalidCodon(String),
    InvalidAnticodon(String),
}

impl std::fmt::Display for GeneticCodeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            GeneticCodeError::InvalidCodeId(id) => write!(f, "Invalid genetic code ID: {}", id),
            GeneticCodeError::InvalidCodon(c) => write!(f, "Invalid codon: {}", c),
            GeneticCodeError::InvalidAnticodon(a) => write!(f, "Invalid anticodon: {}", a),
        }
    }
}

impl std::error::Error for GeneticCodeError {}

// Standard genetic code IDs from NCBI
pub const STANDARD_CODE: i32 = 1;
pub const VERTEBRATE_MITO: i32 = 2;
pub const YEAST_MITO: i32 = 3;
pub const MOLD_MITO: i32 = 4;
pub const INVERTEBRATE_MITO: i32 = 5;
pub const CILIATE_NUCLEAR: i32 = 6;
pub const ECHINODERM_MITO: i32 = 9;
pub const EUPLOTID_NUCLEAR: i32 = 10;
pub const BACTERIAL: i32 = 11;
pub const ALTERNATIVE_YEAST: i32 = 12;
pub const ASCIDIAN_MITO: i32 = 13;
pub const FLATWORM_MITO: i32 = 14;
pub const CHLOROPHYCEAN_MITO: i32 = 16;
pub const TREMATODE_MITO: i32 = 21;
pub const THRAUSTOCHYTRIUM_MITO: i32 = 23;
pub const PTEROBRANCH_MITO: i32 = 24;
pub const SR1_GRACILIBACTERIA: i32 = 25;

/// Genetic code table
#[derive(Debug, Clone)]
pub struct GeneticCode {
    pub id: i32,
    pub name: String,
    pub start_codons: Vec<String>,
    pub stop_codons: Vec<String>,
    table: HashMap<String, char>,
    one_letter_map: HashMap<String, char>,
    /// Anticodon to isotype mapping (e.g., "TGC" -> "Ala")
    anticodon_to_isotype: HashMap<String, String>,
}

impl GeneticCode {
    /// Create a new genetic code from ID
    ///
    /// # Arguments
    /// * `id` - NCBI genetic code ID
    ///
    /// # Returns
    /// Result containing GeneticCode or error
    pub fn new(id: i32) -> Result<Self, GeneticCodeError> {
        match id {
            STANDARD_CODE | BACTERIAL => Ok(Self::standard_code()),
            VERTEBRATE_MITO => Ok(Self::vertebrate_mito()),
            _ => Err(GeneticCodeError::InvalidCodeId(id)),
        }
    }

    /// Standard genetic code (ID 1)
    fn standard_code() -> Self {
        let mut table = HashMap::new();
        let mut one_letter_map = HashMap::new();

        // Build standard codon table
        let codons = [
            // GCN -> Ala (A)
            ("GCA", 'A'), ("GCC", 'A'), ("GCG", 'A'), ("GCT", 'A'),
            // TGY -> Cys (C)
            ("TGC", 'C'), ("TGT", 'C'),
            // GAY -> Asp (D)
            ("GAC", 'D'), ("GAT", 'D'),
            // GAR -> Glu (E)
            ("GAA", 'E'), ("GAG", 'E'),
            // TTY -> Phe (F)
            ("TTC", 'F'), ("TTT", 'F'),
            // GGN -> Gly (G)
            ("GGA", 'G'), ("GGC", 'G'), ("GGG", 'G'), ("GGT", 'G'),
            // CAY -> His (H)
            ("CAC", 'H'), ("CAT", 'H'),
            // ATH -> Ile (I)
            ("ATA", 'I'), ("ATC", 'I'), ("ATT", 'I'),
            // AAR -> Lys (K)
            ("AAA", 'K'), ("AAG", 'K'),
            // TTR, CTN -> Leu (L)
            ("TTA", 'L'), ("TTG", 'L'), ("CTA", 'L'), ("CTC", 'L'), ("CTG", 'L'), ("CTT", 'L'),
            // ATG -> Met (M)
            ("ATG", 'M'),
            // AAY -> Asn (N)
            ("AAC", 'N'), ("AAT", 'N'),
            // CCN -> Pro (P)
            ("CCA", 'P'), ("CCC", 'P'), ("CCG", 'P'), ("CCT", 'P'),
            // CAR -> Gln (Q)
            ("CAA", 'Q'), ("CAG", 'Q'),
            // AGR, CGN -> Arg (R)
            ("AGA", 'R'), ("AGG", 'R'), ("CGA", 'R'), ("CGC", 'R'), ("CGG", 'R'), ("CGT", 'R'),
            // AGY, TCN -> Ser (S)
            ("AGC", 'S'), ("AGT", 'S'), ("TCA", 'S'), ("TCC", 'S'), ("TCG", 'S'), ("TCT", 'S'),
            // ACN -> Thr (T)
            ("ACA", 'T'), ("ACC", 'T'), ("ACG", 'T'), ("ACT", 'T'),
            // GTN -> Val (V)
            ("GTA", 'V'), ("GTC", 'V'), ("GTG", 'V'), ("GTT", 'V'),
            // TGG -> Trp (W)
            ("TGG", 'W'),
            // TAY -> Tyr (Y)
            ("TAC", 'Y'), ("TAT", 'Y'),
            // TAR -> Stop (*)
            ("TAA", '*'), ("TAG", '*'),
            // TGA -> Stop (*) or SeC (U)
            ("TGA", '*'),
        ];

        for (codon, aa) in &codons {
            table.insert(codon.to_string(), *aa);
        }

        // Three-letter amino acid codes to one-letter
        let aa_codes = [
            ("Ala", 'A'), ("Cys", 'C'), ("Asp", 'D'), ("Glu", 'E'),
            ("Phe", 'F'), ("Gly", 'G'), ("His", 'H'), ("Ile", 'I'),
            ("Lys", 'K'), ("Leu", 'L'), ("Met", 'M'), ("Asn", 'N'),
            ("Pro", 'P'), ("Gln", 'Q'), ("Arg", 'R'), ("Ser", 'S'),
            ("Thr", 'T'), ("Val", 'V'), ("Trp", 'W'), ("Tyr", 'Y'),
            ("Sup", '?'), ("Supres", '?'), ("SeC", 'Z'), ("SelCys", 'Z'),
            ("Undet", '?'), ("???", '?'),
        ];

        for (name, letter) in &aa_codes {
            one_letter_map.insert(name.to_string(), *letter);
        }

        // Anticodon to isotype mapping (standard code)
        let mut anticodon_to_isotype = HashMap::new();
        let ac_map = [
            // Ala
            ("AGC", "Ala"), ("GGC", "Ala"), ("CGC", "Ala"), ("TGC", "Ala"),
            // Gly
            ("ACC", "Gly"), ("GCC", "Gly"), ("CCC", "Gly"), ("TCC", "Gly"),
            // Pro
            ("AGG", "Pro"), ("GGG", "Pro"), ("CGG", "Pro"), ("TGG", "Pro"),
            // Thr
            ("AGT", "Thr"), ("GGT", "Thr"), ("CGT", "Thr"), ("TGT", "Thr"),
            // Val
            ("AAC", "Val"), ("GAC", "Val"), ("CAC", "Val"), ("TAC", "Val"),
            // Ser
            ("AGA", "Ser"), ("GGA", "Ser"), ("CGA", "Ser"), ("TGA", "Ser"),
            ("ACT", "Ser"), ("GCT", "Ser"),
            // Arg
            ("ACG", "Arg"), ("GCG", "Arg"), ("CCG", "Arg"), ("TCG", "Arg"),
            ("CCT", "Arg"), ("TCT", "Arg"),
            // Leu
            ("AAG", "Leu"), ("GAG", "Leu"), ("CAG", "Leu"), ("TAG", "Leu"),
            ("CAA", "Leu"), ("TAA", "Leu"),
            // Phe
            ("AAA", "Phe"), ("GAA", "Phe"),
            // Asn
            ("ATT", "Asn"), ("GTT", "Asn"),
            // Lys
            ("CTT", "Lys"), ("TTT", "Lys"),
            // Asp
            ("ATC", "Asp"), ("GTC", "Asp"),
            // Glu
            ("CTC", "Glu"), ("TTC", "Glu"),
            // His
            ("ATG", "His"), ("GTG", "His"),
            // Gln
            ("CTG", "Gln"), ("TTG", "Gln"),
            // Tyr
            ("ATA", "Tyr"), ("GTA", "Tyr"),
            // Stop suppressors
            ("CTA", "Supres"), ("TTA", "Supres"),
            // Ile
            ("AAT", "Ile"), ("GAT", "Ile"), ("CAT", "Ile"), ("TAT", "Ile"),
            // Met (CAT also used for Ile, but can be iMet)
            // Cys
            ("ACA", "Cys"), ("GCA", "Cys"),
            // Trp
            ("CCA", "Trp"),
            // SelCys (selenocysteine)
            ("TCA", "SelCys"),
        ];

        for (ac, iso) in &ac_map {
            anticodon_to_isotype.insert(ac.to_string(), iso.to_string());
        }

        GeneticCode {
            id: STANDARD_CODE,
            name: "Standard".to_string(),
            start_codons: vec!["ATG".to_string()],
            stop_codons: vec!["TAA".to_string(), "TAG".to_string(), "TGA".to_string()],
            table,
            one_letter_map,
            anticodon_to_isotype,
        }
    }

    /// Vertebrate mitochondrial genetic code (ID 2)
    fn vertebrate_mito() -> Self {
        let mut gc = Self::standard_code();
        gc.id = VERTEBRATE_MITO;
        gc.name = "Vertebrate Mitochondrial".to_string();

        // Modifications from standard code
        // AGA, AGG -> Stop (not Arg)
        gc.table.insert("AGA".to_string(), '*');
        gc.table.insert("AGG".to_string(), '*');
        // ATA -> Met (not Ile)
        gc.table.insert("ATA".to_string(), 'M');
        // TGA -> Trp (not Stop)
        gc.table.insert("TGA".to_string(), 'W');

        gc.start_codons = vec![
            "ATG".to_string(),
            "ATA".to_string(),
            "ATT".to_string(),
            "ATC".to_string(),
        ];
        gc.stop_codons = vec![
            "TAA".to_string(),
            "TAG".to_string(),
            "AGA".to_string(),
            "AGG".to_string(),
        ];

        // Update anticodon mapping for vertebrate mito
        let vert_mito_ac = [
            ("TGC", "Ala"), ("TCC", "Gly"), ("TGG", "Pro"), ("TGT", "Thr"), ("TAC", "Val"),
            ("TGA", "Ser"), ("GCT", "Ser"), ("TCG", "Arg"), ("TAG", "Leu"), ("TAA", "Leu"),
            ("GAA", "Phe"), ("GTT", "Asn"), ("TTT", "Lys"), ("GTC", "Asp"), ("TTC", "Glu"),
            ("GTG", "His"), ("TTG", "Gln"), ("GTA", "Tyr"),
            ("GAT", "Ile"), ("TAT", "Met"), ("CAT", "Met"),
            ("GCA", "Cys"), ("TCA", "Trp"), ("GCC", "Asp"),
        ];

        gc.anticodon_to_isotype.clear();
        for (ac, iso) in &vert_mito_ac {
            gc.anticodon_to_isotype.insert(ac.to_string(), iso.to_string());
        }

        gc
    }

    /// Translate a codon to amino acid
    ///
    /// # Arguments
    /// * `codon` - 3-letter codon sequence
    ///
    /// # Returns
    /// Some(amino acid) or None if invalid
    pub fn translate(&self, codon: &str) -> Option<char> {
        let codon_upper = codon.to_uppercase();
        self.table.get(&codon_upper).copied()
    }

    /// Check if codon is a start codon
    pub fn is_start_codon(&self, codon: &str) -> bool {
        let codon_upper = codon.to_uppercase();
        self.start_codons.contains(&codon_upper)
    }

    /// Check if codon is a stop codon
    pub fn is_stop_codon(&self, codon: &str) -> bool {
        let codon_upper = codon.to_uppercase();
        self.stop_codons.contains(&codon_upper)
    }

    /// Get isotype from anticodon
    ///
    /// # Arguments
    /// * `anticodon` - 3-letter anticodon sequence
    ///
    /// # Returns
    /// Isotype name (e.g., "Ala", "Gly") or "Undet" if unknown
    pub fn isotype_from_anticodon(&self, anticodon: &str) -> String {
        let ac_upper = anticodon.to_uppercase();
        self.anticodon_to_isotype
            .get(&ac_upper)
            .cloned()
            .unwrap_or_else(|| "Undet".to_string())
    }

    /// Convert anticodon to corresponding codon (reverse complement)
    ///
    /// # Arguments
    /// * `anticodon` - Anticodon sequence
    ///
    /// # Returns
    /// Codon sequence
    pub fn anticodon_to_codon(&self, anticodon: &str) -> String {
        anticodon
            .chars()
            .rev()
            .map(|c| match c {
                'A' | 'a' => 'T',
                'T' | 't' | 'U' | 'u' => 'A',
                'G' | 'g' => 'C',
                'C' | 'c' => 'G',
                _ => 'N',
            })
            .collect()
    }

    /// Get amino acid from anticodon (via codon translation)
    ///
    /// # Arguments
    /// * `anticodon` - Anticodon sequence
    ///
    /// # Returns
    /// Some(amino acid) or None if invalid
    pub fn anticodon_to_aa(&self, anticodon: &str) -> Option<char> {
        let codon = self.anticodon_to_codon(anticodon);
        self.translate(&codon)
    }

    /// Get all possible amino acids from anticodon considering wobble
    ///
    /// # Arguments
    /// * `anticodon` - Anticodon sequence
    ///
    /// # Returns
    /// Vector of possible amino acids
    pub fn get_wobble_aa(&self, anticodon: &str) -> Vec<char> {
        let mut aas = Vec::new();

        // Get base codon
        if let Some(aa) = self.anticodon_to_aa(anticodon) {
            aas.push(aa);
        }

        // Try wobble variants at first position (5' of anticodon)
        let wobble_bases = ['A', 'G', 'C', 'T', 'U', 'I'];
        for base in &wobble_bases {
            let mut wobble_ac = anticodon.to_string();
            if !wobble_ac.is_empty() {
                wobble_ac.replace_range(0..1, &base.to_string());
                if let Some(aa) = self.anticodon_to_aa(&wobble_ac) {
                    if !aas.contains(&aa) {
                        aas.push(aa);
                    }
                }
            }
        }

        aas
    }

    /// Check if anticodon is a stop suppressor
    ///
    /// # Arguments
    /// * `anticodon` - Anticodon sequence
    ///
    /// # Returns
    /// True if anticodon suppresses stop codons
    pub fn is_suppressor(&self, anticodon: &str) -> bool {
        let codon = self.anticodon_to_codon(anticodon);
        self.is_stop_codon(&codon)
    }

    /// Check if anticodon codes for selenocysteine
    pub fn is_selenocysteine(&self, anticodon: &str) -> bool {
        let ac_upper = anticodon.to_uppercase();
        ac_upper == "TCA" // Selenocysteine anticodon
    }

    /// Check if anticodon codes for pyrrolysine
    pub fn is_pyrrolysine(&self, anticodon: &str) -> bool {
        let ac_upper = anticodon.to_uppercase();
        ac_upper == "CTA" // Pyrrolysine anticodon (TAG stop codon suppressor)
    }

    /// Convert three-letter amino acid code to one-letter
    pub fn three_to_one(&self, three_letter: &str) -> Option<char> {
        self.one_letter_map.get(three_letter).copied()
    }
}

/// Get genetic code by ID
///
/// # Arguments
/// * `id` - NCBI genetic code ID
///
/// # Returns
/// Result containing GeneticCode or error
pub fn get_genetic_code(id: i32) -> Result<GeneticCode, GeneticCodeError> {
    GeneticCode::new(id)
}

/// List all available genetic codes
///
/// # Returns
/// Vector of (id, name) tuples
pub fn list_genetic_codes() -> Vec<(i32, &'static str)> {
    vec![
        (STANDARD_CODE, "Standard / Universal"),
        (VERTEBRATE_MITO, "Vertebrate Mitochondrial"),
        (YEAST_MITO, "Yeast Mitochondrial"),
        (MOLD_MITO, "Mold/Protozoan/Coelenterate Mitochondrial"),
        (INVERTEBRATE_MITO, "Invertebrate Mitochondrial"),
        (CILIATE_NUCLEAR, "Ciliate/Dasycladacean/Hexamita Nuclear"),
        (ECHINODERM_MITO, "Echinoderm/Flatworm Mitochondrial"),
        (EUPLOTID_NUCLEAR, "Euplotid Nuclear"),
        (BACTERIAL, "Bacterial/Archaeal/Plant Plastid"),
        (ALTERNATIVE_YEAST, "Alternative Yeast Nuclear"),
        (ASCIDIAN_MITO, "Ascidian Mitochondrial"),
        (FLATWORM_MITO, "Alternative Flatworm Mitochondrial"),
        (CHLOROPHYCEAN_MITO, "Chlorophycean Mitochondrial"),
        (TREMATODE_MITO, "Trematode Mitochondrial"),
        (THRAUSTOCHYTRIUM_MITO, "Thraustochytrium Mitochondrial"),
        (PTEROBRANCH_MITO, "Pterobranchia Mitochondrial"),
        (SR1_GRACILIBACTERIA, "Candidate Division SR1 and Gracilibacteria"),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_standard_code() {
        let gc = GeneticCode::new(STANDARD_CODE).unwrap();

        // Test translation
        assert_eq!(gc.translate("ATG"), Some('M'));
        assert_eq!(gc.translate("GCA"), Some('A'));
        assert_eq!(gc.translate("TAA"), Some('*'));

        // Test start/stop codons
        assert!(gc.is_start_codon("ATG"));
        assert!(gc.is_stop_codon("TAA"));
        assert!(gc.is_stop_codon("TAG"));
        assert!(gc.is_stop_codon("TGA"));
    }

    #[test]
    fn test_vertebrate_mito() {
        let gc = GeneticCode::new(VERTEBRATE_MITO).unwrap();

        // Test mito-specific changes
        assert_eq!(gc.translate("AGA"), Some('*')); // Stop in mito
        assert_eq!(gc.translate("ATA"), Some('M')); // Met in mito
        assert_eq!(gc.translate("TGA"), Some('W')); // Trp in mito
    }

    #[test]
    fn test_anticodon_to_codon() {
        let gc = GeneticCode::new(STANDARD_CODE).unwrap();

        // TGC anticodon -> GCA codon (Ala)
        assert_eq!(gc.anticodon_to_codon("TGC"), "GCA");
        // CAT anticodon -> ATG codon (Met)
        assert_eq!(gc.anticodon_to_codon("CAT"), "ATG");
    }

    #[test]
    fn test_isotype_from_anticodon() {
        let gc = GeneticCode::new(STANDARD_CODE).unwrap();

        assert_eq!(gc.isotype_from_anticodon("TGC"), "Ala");
        assert_eq!(gc.isotype_from_anticodon("CAT"), "Ile");
        assert_eq!(gc.isotype_from_anticodon("CCA"), "Trp");
    }

    #[test]
    fn test_suppressor() {
        let gc = GeneticCode::new(STANDARD_CODE).unwrap();

        // CTA anticodon reads TAG stop codon
        assert!(gc.is_suppressor("CTA"));
        // TGC anticodon reads GCA (Ala)
        assert!(!gc.is_suppressor("TGC"));
    }

    #[test]
    fn test_selenocysteine() {
        let gc = GeneticCode::new(STANDARD_CODE).unwrap();

        assert!(gc.is_selenocysteine("TCA"));
        assert!(!gc.is_selenocysteine("CCA"));
    }
}
