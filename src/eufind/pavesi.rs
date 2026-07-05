// Pavesi algorithm for finding eukaryotic tRNA transcriptional control regions
// Based on original pavesi.c implementation

use super::constants::*;

/// tRNA information structure
/// Corresponds to TRNA_TYPE in original eufind_const.h
#[derive(Debug, Clone)]
pub struct TrnaInfo {
    pub iso_type: String,    // Isotype (5 chars)
    pub acodon: String,      // Anticodon (4 chars)
    pub start: i32,          // tRNA start position
    pub end: i32,            // tRNA end position
    pub abox_st: i32,        // A-box start position
    pub abox_end: i32,       // A-box end position
    pub abox_gap: i32,       // A-box gap
    pub bbox_st: i32,        // B-box start position
    pub bbox_end: i32,       // B-box end position
    pub term_st: i32,        // Termination signal start position
    pub acodon_idx: i32,     // Anticodon index
    pub intron: i32,         // Intron position
    pub idno: i32,           // ID number
    pub tot_sc: f32,         // Total score
    pub abox_sc: f32,        // A-box score
    pub bbox_sc: f32,        // B-box score
    pub abdist_sc: f32,      // A-B box distance score
    pub term_sc: f32,        // Termination signal score
}

impl TrnaInfo {
    /// Initialize a new TrnaInfo structure with default values
    pub fn new() -> Self {
        TrnaInfo {
            iso_type: "???".to_string(),
            acodon: "???".to_string(),
            start: 0,
            end: 0,
            abox_st: 0,
            abox_end: 0,
            abox_gap: 0,
            bbox_st: 0,
            bbox_end: 0,
            term_st: 0,
            acodon_idx: 0,
            intron: 0,
            idno: 0,
            tot_sc: -1000.0,
            abox_sc: -1000.0,
            bbox_sc: -1000.0,
            abdist_sc: -100.0,
            term_sc: -100.0,
        }
    }
}

impl Default for TrnaInfo {
    fn default() -> Self {
        Self::new()
    }
}

/// Encode DNA sequence as integers (A=0, C=1, G=2, T=3, other=5)
///
/// # Arguments
/// * `seq` - DNA sequence string
///
/// # Returns
/// Vector of integer-encoded bases
pub fn int_encode_seq(seq: &str) -> Vec<i32> {
    seq.bytes()
        .map(|c| match c {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' | b'U' | b'u' => 3,
            _ => 5, // Ambiguous base
        })
        .collect()
}

/// Score a B-box at a specific position
///
/// # Arguments
/// * `iseq` - Integer-encoded sequence
/// * `pos` - Position to score
///
/// # Returns
/// B-box score
pub fn get_bbox(iseq: &[i32], pos: usize) -> f32 {
    let mut score = 0.0;
    for j in 0..BBOX_LEN {
        let base_idx = iseq[pos + j] as usize;
        score += BBOX_MAT[base_idx][j];
    }
    score
}

/// Find best A-box given a B-box position
///
/// # Arguments
/// * `iseq` - Integer-encoded sequence
/// * `seq` - Original sequence string (for gap detection)
/// * `bbox_st` - B-box start position
/// * `max_ab_dist` - Maximum A-B box distance
/// * `prev_abox_st` - Previous A-box start position (for searching multiple)
///
/// # Returns
/// Tuple of (abox_st, abox_end, abox_gap, abox_sc, abdist_sc)
pub fn get_best_abox(
    iseq: &[i32],
    seq: &[u8],
    bbox_st: i32,
    max_ab_dist: usize,
    prev_abox_st: i32,
) -> (i32, i32, i32, f32, f32) {
    let bbox_st = bbox_st as usize;
    let startidx = std::cmp::max(
        std::cmp::max(0, bbox_st as i32 - max_ab_dist as i32 - ABOX_LEN as i32),
        prev_abox_st + 2,
    ) as usize;
    let endidx = std::cmp::max(0, bbox_st as i32 - MIN_AB_BOX_DIST as i32 - ABOX_LEN as i32 + 4) as usize;

    let mut best_abox_st = 0;
    let mut best_abox_end = 0;
    let mut best_abox_sc = -1000.0;
    let mut best_abdist_sc = -1000.0;
    let mut best_gap = 0;

    for i in startidx..endidx {
        // Score positions 7-16 (0-9 in 0-indexed)
        let mut sc1 = 0.0;
        for j in 0..10 {
            sc1 += ABOX_MAT[iseq[i + j] as usize][j];
        }

        // Score gap at position 17 by looking for conserved 'G' at position 18
        let (sc2_base, offset1) = if seq[i + 10] == b'G' {
            (ABOX_MAT[GAP_ROW][10], 1)
        } else {
            (0.0, 0)
        };

        // Score positions 18-20
        let mut sc2 = sc2_base;
        for j in 10..(14 - offset1) {
            sc2 += ABOX_MAT[iseq[i + j] as usize][j + offset1];
        }

        // Try all 4 possible gap configurations at positions 20a & 20b
        for gapct in 0..4 {
            let j = 14 - offset1;
            let offset2;
            let mut sc3;

            match gapct {
                0 => {
                    // No gap
                    sc3 = ABOX_MAT[iseq[i + j] as usize][j + offset1]
                        + ABOX_MAT[iseq[i + j + 1] as usize][j + offset1 + 1];
                    offset2 = 0;
                }
                1 => {
                    // Gap after first base
                    offset2 = 1;
                    sc3 = ABOX_MAT[GAP_ROW][j + offset1]
                        + ABOX_MAT[iseq[i + j] as usize][j + offset1 + offset2];
                }
                2 => {
                    // Gap after second base
                    offset2 = 1;
                    sc3 = ABOX_MAT[iseq[i + j] as usize][j + offset1]
                        + ABOX_MAT[GAP_ROW][j + offset1 + 1];
                }
                3 => {
                    // Two gaps
                    offset2 = 2;
                    sc3 = ABOX_MAT[GAP_ROW][j + offset1]
                        + ABOX_MAT[GAP_ROW][j + offset1 + 1];
                }
                _ => unreachable!(),
            }

            // Score remaining positions
            let mut jj = j + if gapct == 0 { 2 } else { 1 };
            while (jj + offset1 + offset2) < ABOX_LEN {
                sc3 += ABOX_MAT[iseq[i + jj] as usize][jj + offset1 + offset2];
                jj += 1;
            }

            let abox_end = (i + ABOX_LEN - offset1 - offset2 - 1) as i32;
            let abdist_sc = get_abdist_weight(bbox_st as i32 - abox_end - 1);
            let total_sc = sc1 + sc2 + sc3 + abdist_sc;

            if total_sc > (best_abox_sc + best_abdist_sc) {
                best_abox_st = i as i32;
                best_abox_end = abox_end;
                best_abox_sc = sc1 + sc2 + sc3;
                best_abdist_sc = abdist_sc;
                best_gap = gapct;
            }
        }
    }

    (best_abox_st, best_abox_end, best_gap, best_abox_sc, best_abdist_sc)
}

/// Get weight for A-B box distance
///
/// # Arguments
/// * `ab_dist` - Distance between A-box end and B-box start
///
/// # Returns
/// Distance score
pub fn get_abdist_weight(ab_dist: i32) -> f32 {
    if ab_dist < MIN_AB_BOX_DIST as i32 {
        return MAX_PENALTY;
    }

    for ct in 0..ABDIST_MAT_SIZE {
        if ab_dist <= AB_DIST_IDX_MAT[ct] {
            return AB_DIST_SC_MAT[ct];
        }
    }

    MAX_PENALTY
}

/// Find best transcription termination signal (TTTT)
///
/// # Arguments
/// * `seq` - DNA sequence (bytes)
/// * `bbox_end` - B-box end position
/// * `seqlen` - Sequence length
///
/// # Returns
/// Tuple of (term_st, term_sc)
pub fn get_best_trx_term(seq: &[u8], bbox_end: i32, seqlen: usize) -> (i32, f32) {
    let startidx = (bbox_end + MIN_BTERM_DIST as i32 - 1) as usize;
    let endidx = std::cmp::min(startidx + MAX_TERM_SEARCH, seqlen - 4);

    for i in startidx..endidx {
        if seq[i] == b'T' && seq[i + 1] == b'T' && seq[i + 2] == b'T' && seq[i + 3] == b'T' {
            let bterm_dist = i as i32 - bbox_end - 1;
            let mut ct = 0;
            while ct < BTERM_MAT_SIZE && BTERM_DIST_IDX_MAT[ct] < bterm_dist {
                ct += 1;
            }
            let term_sc = if ct < BTERM_MAT_SIZE {
                BTERM_DIST_SC_MAT[ct]
            } else {
                MAX_PENALTY
            };
            return (i as i32, term_sc);
        }
    }

    // No termination signal found
    if endidx == seqlen - 4 {
        // At end of sequence, use threshold-based score
        let term_sc = TOT_SCORE_THRESH - INT_SCORE_THRESH;
        (-1, term_sc)
    } else {
        // Not at end and no signal found
        (-1, MAX_PENALTY)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_int_encode_seq() {
        let seq = "ACGTACGT";
        let encoded = int_encode_seq(seq);
        assert_eq!(encoded, vec![0, 1, 2, 3, 0, 1, 2, 3]);

        // Test ambiguous bases
        let seq_amb = "ACGTN";
        let encoded_amb = int_encode_seq(seq_amb);
        assert_eq!(encoded_amb, vec![0, 1, 2, 3, 5]);
    }

    #[test]
    fn test_trna_info_new() {
        let trna = TrnaInfo::new();
        assert_eq!(trna.iso_type, "???");
        assert_eq!(trna.acodon, "???");
        assert_eq!(trna.start, 0);
        assert_eq!(trna.tot_sc, -1000.0);
    }

    #[test]
    fn test_get_abdist_weight() {
        // Test below minimum
        assert_eq!(get_abdist_weight(20), MAX_PENALTY);

        // Test at specific thresholds
        assert_eq!(get_abdist_weight(30), -0.46);
        assert_eq!(get_abdist_weight(35), -1.83);
        assert_eq!(get_abdist_weight(40), -2.35);

        // Test above maximum
        assert_eq!(get_abdist_weight(70), MAX_PENALTY);
    }
}
