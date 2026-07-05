/// Sprinzl position numbering system for tRNA
///
/// Standard tRNA positions numbered 1-76 with insertions (17a, 20a, 20b)
/// and variable loop extensions (e1-e27).

use std::collections::HashMap;

/// Sprinzl position constants
pub const SPRINZL_POSITIONS: &[&str] = &[
    "1", "2", "3", "4", "5", "6", "7",
    "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "17a", "18", "19", "20", "20a", "20b", "21",
    "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38",
    "39", "40", "41", "42", "43", "44", "45",
    "e11", "e12", "e13", "e14", "e15", "e16", "e17",
    "e1", "e2", "e3", "e4", "e5",
    "e27", "e26", "e25", "e24", "e23", "e22", "e21",
    "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58",
    "59", "60", "61", "62", "63", "64", "65", "66", "67", "68", "69", "70", "71",
    "72", "73", "74", "75", "76"
];

/// Sprinzl position manager
pub struct SprinzlPositions {
    /// Base pairing partners
    pairs: HashMap<String, String>,
    /// Reverse pairing map
    rev_pairs: HashMap<String, String>,
    /// Stem position codes (A1-A7, D1-D4, C1-C5, V1-V7, T1-T5)
    stem_pos: HashMap<String, String>,
    /// Secondary structure region codes
    ss_pos: HashMap<String, String>,
    /// Region descriptions
    regions: HashMap<String, String>,
    /// Universal/conserved bases at specific positions
    universal: HashMap<String, char>,
}

impl SprinzlPositions {
    pub fn new() -> Self {
        let mut pairs = HashMap::new();
        let mut rev_pairs = HashMap::new();

        // Acceptor stem pairs
        for i in 1..=7 {
            pairs.insert(i.to_string(), (73-i).to_string());
            rev_pairs.insert((73-i).to_string(), i.to_string());
        }

        // D-arm pairs
        pairs.insert("10".to_string(), "25".to_string());
        pairs.insert("11".to_string(), "24".to_string());
        pairs.insert("12".to_string(), "23".to_string());
        pairs.insert("13".to_string(), "22".to_string());
        rev_pairs.insert("25".to_string(), "10".to_string());
        rev_pairs.insert("24".to_string(), "11".to_string());
        rev_pairs.insert("23".to_string(), "12".to_string());
        rev_pairs.insert("22".to_string(), "13".to_string());

        // Anticodon stem pairs
        pairs.insert("27".to_string(), "43".to_string());
        pairs.insert("28".to_string(), "42".to_string());
        pairs.insert("29".to_string(), "41".to_string());
        pairs.insert("30".to_string(), "40".to_string());
        pairs.insert("31".to_string(), "39".to_string());
        rev_pairs.insert("43".to_string(), "27".to_string());
        rev_pairs.insert("42".to_string(), "28".to_string());
        rev_pairs.insert("41".to_string(), "29".to_string());
        rev_pairs.insert("40".to_string(), "30".to_string());
        rev_pairs.insert("39".to_string(), "31".to_string());

        // Variable stem pairs
        pairs.insert("e11".to_string(), "e21".to_string());
        pairs.insert("e12".to_string(), "e22".to_string());
        pairs.insert("e13".to_string(), "e23".to_string());
        pairs.insert("e14".to_string(), "e24".to_string());
        pairs.insert("e15".to_string(), "e25".to_string());
        pairs.insert("e16".to_string(), "e26".to_string());
        pairs.insert("e17".to_string(), "e27".to_string());
        rev_pairs.insert("e21".to_string(), "e11".to_string());
        rev_pairs.insert("e22".to_string(), "e12".to_string());
        rev_pairs.insert("e23".to_string(), "e13".to_string());
        rev_pairs.insert("e24".to_string(), "e14".to_string());
        rev_pairs.insert("e25".to_string(), "e15".to_string());
        rev_pairs.insert("e26".to_string(), "e16".to_string());
        rev_pairs.insert("e27".to_string(), "e17".to_string());

        // T-arm pairs
        pairs.insert("49".to_string(), "65".to_string());
        pairs.insert("50".to_string(), "64".to_string());
        pairs.insert("51".to_string(), "63".to_string());
        pairs.insert("52".to_string(), "62".to_string());
        pairs.insert("53".to_string(), "61".to_string());
        rev_pairs.insert("65".to_string(), "49".to_string());
        rev_pairs.insert("64".to_string(), "50".to_string());
        rev_pairs.insert("63".to_string(), "51".to_string());
        rev_pairs.insert("62".to_string(), "52".to_string());
        rev_pairs.insert("61".to_string(), "53".to_string());

        // Stem position codes
        let mut stem_pos = HashMap::new();
        for i in 1..=7 {
            stem_pos.insert(i.to_string(), format!("A{}", i));
        }
        for i in 10..=13 {
            stem_pos.insert(i.to_string(), format!("D{}", i-9));
        }
        for i in 27..=31 {
            stem_pos.insert(i.to_string(), format!("C{}", i-26));
        }
        for i in 1..=7 {
            stem_pos.insert(format!("e1{}", i), format!("V{}", i));
        }
        for i in 49..=53 {
            stem_pos.insert(i.to_string(), format!("T{}", i-48));
        }

        // Secondary structure positions
        let mut ss_pos = HashMap::new();
        for i in 1..=7 {
            ss_pos.insert(i.to_string(), "5P1".to_string());
            ss_pos.insert((73-i).to_string(), "3P1".to_string());
        }
        ss_pos.insert("8".to_string(), "L1".to_string());
        ss_pos.insert("9".to_string(), "L1".to_string());

        for i in 10..=13 {
            ss_pos.insert(i.to_string(), "5P2".to_string());
            ss_pos.insert((35-i).to_string(), "3P2".to_string());
        }

        for pos in &["14", "15", "16", "17", "17a", "18", "19", "20", "20a", "20b", "21"] {
            ss_pos.insert(pos.to_string(), "L2".to_string());
        }
        ss_pos.insert("26".to_string(), "L3".to_string());

        for i in 27..=31 {
            ss_pos.insert(i.to_string(), "5P3".to_string());
            ss_pos.insert((70-i).to_string(), "3P3".to_string());
        }

        for i in 32..=38 {
            ss_pos.insert(i.to_string(), "L4".to_string());
        }

        for i in 44..=48 {
            ss_pos.insert(i.to_string(), "L5".to_string());
        }

        for pos in &["e11", "e12", "e13", "e14", "e15", "e16", "e17",
                     "e1", "e2", "e3", "e4", "e5",
                     "e21", "e22", "e23", "e24", "e25", "e26", "e27"] {
            ss_pos.insert(pos.to_string(), "P4".to_string());
        }

        for i in 49..=53 {
            ss_pos.insert(i.to_string(), "5P5".to_string());
            ss_pos.insert((114-i).to_string(), "3P5".to_string());
        }

        for i in 54..=60 {
            ss_pos.insert(i.to_string(), "L6".to_string());
        }

        for i in 73..=76 {
            ss_pos.insert(i.to_string(), "L7".to_string());
        }

        // Region descriptions
        let mut regions = HashMap::new();
        regions.insert("5P1".to_string(), "5p Acceptor Stem".to_string());
        regions.insert("3P1".to_string(), "3p Acceptor Stem".to_string());
        regions.insert("L1".to_string(), "Acceptor-D-arm-linker".to_string());
        regions.insert("5P2".to_string(), "5p D-arm".to_string());
        regions.insert("3P2".to_string(), "3p D-arm".to_string());
        regions.insert("L2".to_string(), "D-loop".to_string());
        regions.insert("L3".to_string(), "D-arm-Anticodon-linker".to_string());
        regions.insert("5P3".to_string(), "5p Anticodon Stem".to_string());
        regions.insert("3P3".to_string(), "3p Anticodon Stem".to_string());
        regions.insert("L4".to_string(), "Anticodon Loop".to_string());
        regions.insert("L5".to_string(), "Variable Loop".to_string());
        regions.insert("P4".to_string(), "Variable Stem".to_string());
        regions.insert("5P5".to_string(), "5p T-arm".to_string());
        regions.insert("3P5".to_string(), "3p T-arm".to_string());
        regions.insert("L6".to_string(), "T-Psi-C Loop".to_string());
        regions.insert("L7".to_string(), "3p end".to_string());

        // Universal/conserved bases
        let mut universal = HashMap::new();
        universal.insert("8".to_string(), 'T');
        universal.insert("14".to_string(), 'A');
        universal.insert("18".to_string(), 'G');
        universal.insert("19".to_string(), 'G');
        universal.insert("21".to_string(), 'A');
        universal.insert("33".to_string(), 'T');
        universal.insert("53".to_string(), 'G');
        universal.insert("54".to_string(), 'T');
        universal.insert("55".to_string(), 'T');
        universal.insert("56".to_string(), 'C');
        universal.insert("58".to_string(), 'A');
        universal.insert("74".to_string(), 'C');
        universal.insert("75".to_string(), 'C');
        universal.insert("76".to_string(), 'A');

        Self {
            pairs,
            rev_pairs,
            stem_pos,
            ss_pos,
            regions,
            universal,
        }
    }

    /// Get pairing partner for a Sprinzl position
    pub fn get_pair(&self, pos: &str) -> Option<&String> {
        self.pairs.get(pos)
    }

    /// Get reverse pairing partner
    pub fn get_rev_pair(&self, pos: &str) -> Option<&String> {
        self.rev_pairs.get(pos)
    }

    /// Get stem position code (e.g., "A1", "D2")
    pub fn get_stem_code(&self, pos: &str) -> Option<&String> {
        self.stem_pos.get(pos)
    }

    /// Get secondary structure region for position
    pub fn get_region(&self, pos: &str) -> Option<&String> {
        // Handle intron notation (e.g., "32:i1")
        let clean_pos = if let Some(idx) = pos.find(":i") {
            &pos[..idx]
        } else {
            pos
        };
        self.ss_pos.get(clean_pos)
    }

    /// Get region description
    pub fn get_region_description(&self, region: &str) -> Option<&String> {
        self.regions.get(region)
    }

    /// Get position region description directly
    pub fn get_pos_description(&self, pos: &str) -> Option<String> {
        if let Some(region) = self.get_region(pos) {
            self.get_region_description(region).cloned()
        } else {
            None
        }
    }

    /// Get expected universal/conserved base at position
    pub fn get_universal_base(&self, pos: &str) -> Option<char> {
        self.universal.get(pos).copied()
    }

    /// Check if position is in a stem region
    pub fn is_stem_position(&self, pos: &str) -> bool {
        self.pairs.contains_key(pos) || self.rev_pairs.contains_key(pos)
    }

    /// Check if position is in a loop region
    pub fn is_loop_position(&self, pos: &str) -> bool {
        if let Some(region) = self.get_region(pos) {
            region.starts_with('L')
        } else {
            false
        }
    }
}

impl Default for SprinzlPositions {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sprinzl_positions() {
        let sp = SprinzlPositions::new();

        // Test acceptor stem pairing
        assert_eq!(sp.get_pair("1").unwrap(), "72");
        assert_eq!(sp.get_pair("7").unwrap(), "66");

        // Test D-arm pairing
        assert_eq!(sp.get_pair("10").unwrap(), "25");
        assert_eq!(sp.get_pair("13").unwrap(), "22");

        // Test reverse pairing
        assert_eq!(sp.get_rev_pair("72").unwrap(), "1");
    }

    #[test]
    fn test_stem_codes() {
        let sp = SprinzlPositions::new();

        assert_eq!(sp.get_stem_code("1").unwrap(), "A1");
        assert_eq!(sp.get_stem_code("10").unwrap(), "D1");
        assert_eq!(sp.get_stem_code("27").unwrap(), "C1");
    }

    #[test]
    fn test_regions() {
        let sp = SprinzlPositions::new();

        assert_eq!(sp.get_region("1").unwrap(), "5P1");
        assert_eq!(sp.get_region("14").unwrap(), "L2");
        assert_eq!(sp.get_region("34").unwrap(), "L4");
    }

    #[test]
    fn test_universal_bases() {
        let sp = SprinzlPositions::new();

        assert_eq!(sp.get_universal_base("8"), Some('T'));
        assert_eq!(sp.get_universal_base("18"), Some('G'));
        assert_eq!(sp.get_universal_base("76"), Some('A'));
    }

    #[test]
    fn test_position_types() {
        let sp = SprinzlPositions::new();

        assert!(sp.is_stem_position("1"));
        assert!(sp.is_loop_position("14"));
        assert!(!sp.is_loop_position("1"));
    }
}
