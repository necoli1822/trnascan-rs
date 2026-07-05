/// Secondary structure utilities for tRNA analysis
///
/// This module provides functions for parsing and validating tRNA secondary structures,
/// extracting stem regions, and working with dot-bracket notation.

/// Stem region in secondary structure
#[derive(Debug, Clone)]
pub struct Stem {
    pub start_left: usize,
    pub end_left: usize,
    pub start_right: usize,
    pub end_right: usize,
}

impl Stem {
    pub fn length(&self) -> usize {
        self.end_left - self.start_left + 1
    }
}

/// Validation result for tRNA structure
#[derive(Debug, Default)]
pub struct ValidationResult {
    pub is_valid_trna: bool,
    pub is_valid_acceptor: bool,
    pub is_valid_darm: bool,
    pub is_valid_anticodon: bool,
    pub is_valid_variable: bool,
    pub is_valid_tstem: bool,
}

/// Secondary structure representation
pub struct SecondaryStructure {
    pub sequence: Vec<u8>,
    pub structure: String,  // Dot-bracket notation
}

impl SecondaryStructure {
    pub fn new(seq: &[u8], ss: &str) -> Self {
        Self {
            sequence: seq.to_vec(),
            structure: ss.to_string(),
        }
    }

    /// Parse stems from secondary structure notation
    /// Returns (stems, mismatches_per_stem)
    pub fn get_stems(&self) -> (Vec<Stem>, Vec<usize>) {
        let ss = self.structure.as_bytes();
        let mut left_positions = Vec::new();
        let mut right_positions = Vec::new();
        let mut pairs: Vec<Option<usize>> = Vec::new();

        // Find base pairs
        for (pos, &ch) in ss.iter().enumerate() {
            if ch == b'>' {
                left_positions.push(pos);
                pairs.push(None);
            } else if ch == b'<' {
                right_positions.push(pos);

                // Find matching left position
                let mut left_index = left_positions.len() - 1;
                loop {
                    if pairs[left_index].is_none() {
                        pairs[left_index] = Some(right_positions.len() - 1);
                        break;
                    }
                    if left_index == 0 {
                        break;
                    }
                    left_index -= 1;
                }
            }
        }

        // Group pairs into stems
        let mut stems = Vec::new();
        if pairs.is_empty() {
            return (stems, Vec::new());
        }

        let mut start_left_idx = 0;
        let mut last_right_idx = pairs[0];

        for (left_idx, &right_idx_opt) in pairs.iter().enumerate().skip(1) {
            if let Some(right_idx) = right_idx_opt {
                if let Some(last_right) = last_right_idx {
                    if right_idx != last_right - 1 {
                        // New stem starts
                        if let Some(end_right) = last_right_idx {
                            stems.push(Stem {
                                start_left: left_positions[start_left_idx],
                                end_left: left_positions[left_idx - 1],
                                start_right: right_positions[end_right],
                                end_right: right_positions[pairs[start_left_idx].unwrap()],
                            });
                        }
                        start_left_idx = left_idx;
                    }
                }
                last_right_idx = Some(right_idx);
            }
        }

        // Add final stem
        if let Some(last_right) = last_right_idx {
            stems.push(Stem {
                start_left: left_positions[start_left_idx],
                end_left: left_positions[pairs.len() - 1],
                start_right: right_positions[last_right],
                end_right: right_positions[pairs[start_left_idx].unwrap()],
            });
        }

        // Count mismatches in each stem
        let mismatches = stems.iter().map(|stem| {
            let mut count = 0;
            let mut left_idx = stem.end_left;
            let mut right_idx = stem.start_right;

            while left_idx >= stem.start_left && right_idx <= stem.end_right {
                if ss[left_idx] == b'.' && ss[right_idx] == b'.' {
                    count += 1;
                    right_idx += 1;
                }
                if ss[left_idx] == b'.' {
                    left_idx = left_idx.saturating_sub(1);
                } else if ss[right_idx] == b'.' {
                    right_idx += 1;
                } else {
                    left_idx = left_idx.saturating_sub(1);
                    right_idx += 1;
                }

                if left_idx == 0 {
                    break;
                }
            }
            count
        }).collect();

        (stems, mismatches)
    }

    /// Validate tRNA structure
    pub fn validate(&self, canonical_intron_len: usize) -> ValidationResult {
        let mut result = ValidationResult {
            is_valid_trna: true,
            is_valid_acceptor: true,
            is_valid_darm: true,
            is_valid_anticodon: true,
            is_valid_variable: true,
            is_valid_tstem: true,
        };

        let (stems, mismatches) = self.get_stems();
        let total_mismatches: usize = mismatches.iter().sum();

        if total_mismatches > 1 || stems.len() < 4 || stems.len() > 5 {
            result.is_valid_trna = false;
        }

        let ss_len = self.structure.len();

        match stems.len() {
            5 => {
                result.is_valid_acceptor = mismatches[0] <= 1 && stems[0].length() == 7;
                result.is_valid_darm = mismatches[1] <= 1 && stems[1].length() >= 3;
                result.is_valid_anticodon = mismatches[2] <= 1 && stems[2].length() == 5;
                result.is_valid_variable = mismatches[3] <= 1 && stems[3].length() >= 2;
                result.is_valid_tstem = mismatches[4] <= 1 && stems[4].length() == 5;
                result.is_valid_trna = result.is_valid_acceptor
                    && result.is_valid_darm
                    && result.is_valid_anticodon
                    && result.is_valid_variable
                    && result.is_valid_tstem
                    && (ss_len - canonical_intron_len <= 90);
            }
            4 => {
                result.is_valid_variable = false;
                result.is_valid_acceptor = mismatches[0] <= 1 && stems[0].length() == 7;
                result.is_valid_darm = mismatches[1] <= 1 && stems[1].length() >= 3;
                result.is_valid_anticodon = mismatches[2] <= 1 && stems[2].length() == 5;
                result.is_valid_tstem = mismatches[3] <= 1 && stems[3].length() == 5;
                result.is_valid_trna = result.is_valid_acceptor
                    && result.is_valid_darm
                    && result.is_valid_anticodon
                    && result.is_valid_tstem
                    && (ss_len - canonical_intron_len <= 80);
            }
            3 => {
                result.is_valid_variable = false;
                result.is_valid_acceptor = mismatches[0] <= 1 && stems[0].length() == 7;
                if mismatches[1] == 0 {
                    result.is_valid_darm = mismatches[1] <= 1 && stems[1].length() >= 3;
                    result.is_valid_anticodon = mismatches[2] <= 1 && stems[2].length() == 5;
                    result.is_valid_tstem = false;
                } else if mismatches[2] == 0 {
                    result.is_valid_darm = false;
                    result.is_valid_anticodon = mismatches[1] <= 1 && stems[1].length() == 5;
                    result.is_valid_tstem = mismatches[2] <= 1 && stems[2].length() == 5;
                } else {
                    result.is_valid_acceptor = false;
                    result.is_valid_darm = false;
                    result.is_valid_anticodon = false;
                    result.is_valid_tstem = false;
                }
                result.is_valid_trna = false;
            }
            _ => {
                result.is_valid_trna = false;
                result.is_valid_acceptor = false;
                result.is_valid_darm = false;
                result.is_valid_anticodon = false;
                result.is_valid_variable = false;
                result.is_valid_tstem = false;
            }
        }

        result
    }

    /// Get acceptor stem half (5' or 3')
    pub fn get_acceptor_half(&self, half: &str) -> Vec<u8> {
        match half {
            "5h" => {
                if self.sequence.len() >= 7 {
                    self.sequence[..7].to_vec()
                } else {
                    Vec::new()
                }
            }
            "3h" => {
                let ss_len = self.structure.len();
                let seq_len = self.sequence.len();

                if ss_len >= 12 && &self.structure[ss_len - 12..] == "<..........." {
                    if seq_len >= 11 {
                        self.sequence[seq_len - 11..seq_len - 4].to_vec()
                    } else {
                        Vec::new()
                    }
                } else {
                    if seq_len >= 8 {
                        self.sequence[seq_len - 8..seq_len - 1].to_vec()
                    } else {
                        Vec::new()
                    }
                }
            }
            _ => Vec::new(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_validate_structure() {
        let seq = b"GCCGCGGTAGTCTAGTGGTTAGAATACACCCTGATAACGGTGAGGTCGGTGGTTCGATTCCGCTCCGTGGCA";
        let ss = ">>>>>>>..>>>>........<<<<.>>>>>.......<<<<<....>>>>>.......<<<<<<<<<<<<.";

        let sec_struct = SecondaryStructure::new(seq, ss);
        let result = sec_struct.validate(0);

        assert!(result.is_valid_trna);
    }
}
