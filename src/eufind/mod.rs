// EuFindtRNA module - Eukaryotic tRNA finder
// Based on Pavesi et al. algorithm (NAR 22:1247-56, 1994)

pub mod constants;
pub mod pavesi;

// Re-export commonly used items
pub use constants::*;
pub use pavesi::*;

use crate::trna::{TRna, Strand};

/// EuFind scanning options
#[derive(Debug, Clone)]
pub struct EufindOptions {
    /// Maximum intron length to search
    pub max_int_len: usize,
    /// Internal score threshold
    pub int_score_thresh: f32,
    /// Total score threshold
    pub tot_score_thresh: f32,
    /// Start index for scanning
    pub start_index: usize,
    /// A-B box distance range
    pub max_ab_box_dist: usize,
}

impl Default for EufindOptions {
    fn default() -> Self {
        EufindOptions {
            max_int_len: 116,
            int_score_thresh: INT_SCORE_THRESH,
            tot_score_thresh: TOT_SCORE_THRESH,
            start_index: 0,
            max_ab_box_dist: AB_BOX_DIST_RANGE,
        }
    }
}

/// EuFind hit result
#[derive(Debug, Clone)]
pub struct EufindHit {
    pub seqname: String,
    pub start: i64,
    pub end: i64,
    pub isotype: String,
    pub anticodon: String,
    pub score: f64,
    pub strand: Strand,
    pub abox_st: i64,
    pub abox_end: i64,
    pub bbox_st: i64,
    pub bbox_end: i64,
    pub term_st: i64,
}

/// Run EuFind scan on a sequence
///
/// # Arguments
/// * `seq` - DNA sequence to scan
/// * `options` - Scanning options
///
/// # Returns
/// Vector of EuFind hits
pub fn run_eufind_scan(
    seq: &[u8],
    seqname: &str,
    options: &EufindOptions,
) -> Vec<EufindHit> {
    let mut hits = Vec::new();
    let seqlen = seq.len();

    if seqlen < 100 {
        return hits;
    }

    // Encode sequence as integers
    let iseq = int_encode_seq(
        std::str::from_utf8(seq).unwrap_or("")
    );

    // Scan for B-boxes
    let mut pos = BBOX_START_IDX;
    while pos + BBOX_LEN < seqlen {
        let bbox_sc = get_bbox(&iseq, pos);

        if bbox_sc >= BBOX_CUTOFF {
            // Found a B-box, now find best A-box
            let bbox_st = pos as i32;
            let bbox_end = (pos + BBOX_LEN - 1) as i32;

            let (abox_st, abox_end, _abox_gap, abox_sc, abdist_sc) =
                get_best_abox(&iseq, seq, bbox_st, options.max_ab_box_dist, -1);

            // Find termination signal
            let (term_st, term_sc) = get_best_trx_term(seq, bbox_end, seqlen);

            // Calculate total score
            let int_sc = abox_sc + bbox_sc + abdist_sc;
            let tot_sc = int_sc + term_sc;

            if int_sc >= options.int_score_thresh && tot_sc >= options.tot_score_thresh {
                // Determine strand based on A-box and B-box positions
                let strand = if abox_st < bbox_st {
                    Strand::Plus
                } else {
                    Strand::Minus
                };

                let (start, end) = if strand == Strand::Plus {
                    (abox_st as i64, if term_st > 0 { (term_st + 3) as i64 } else { bbox_end as i64 })
                } else {
                    (if term_st > 0 { (term_st + 3) as i64 } else { bbox_end as i64 }, abox_st as i64)
                };

                hits.push(EufindHit {
                    seqname: seqname.to_string(),
                    start,
                    end,
                    isotype: "???".to_string(),
                    anticodon: "???".to_string(),
                    score: tot_sc as f64,
                    strand,
                    abox_st: abox_st as i64,
                    abox_end: abox_end as i64,
                    bbox_st: bbox_st as i64,
                    bbox_end: bbox_end as i64,
                    term_st: term_st as i64,
                });
            }
        }

        pos += 1;
    }

    hits
}

/// Parse EuFind output from external program
///
/// # Arguments
/// * `content` - Output text from eufindtRNA
///
/// # Returns
/// Vector of EuFind hits
pub fn parse_eufind_output(content: &str) -> Vec<EufindHit> {
    let mut hits = Vec::new();

    for line in content.lines() {
        // Parse format: seqname trnact start end isotype anticodon intron score
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 9 {
            let seqname = parts[0].to_string();
            let start = parts[2].parse::<i64>().unwrap_or(0);
            let end = parts[3].parse::<i64>().unwrap_or(0);
            let isotype = parts[4].to_string();
            let anticodon = parts[5].to_string();
            let score = parts[8].parse::<f64>().unwrap_or(0.0);

            let strand = if start < end {
                Strand::Plus
            } else {
                Strand::Minus
            };

            hits.push(EufindHit {
                seqname,
                start,
                end,
                isotype,
                anticodon,
                score,
                strand,
                abox_st: 0,
                abox_end: 0,
                bbox_st: 0,
                bbox_end: 0,
                term_st: 0,
            });
        }
    }

    hits
}

/// Convert EuFind hit to TRna structure
///
/// # Arguments
/// * `hit` - EuFind hit
///
/// # Returns
/// TRna structure
pub fn eufind_to_trna(hit: &EufindHit) -> TRna {
    let mut trna = TRna::default();

    trna.seqname = hit.seqname.clone();
    trna.start = hit.start;
    trna.end = hit.end;
    trna.strand = hit.strand;
    trna.isotype = hit.isotype.clone();
    trna.anticodon = hit.anticodon.clone();
    trna.score = hit.score;

    trna
}

/// Check if two segments overlap
///
/// # Arguments
/// * `start1` - Start of first segment
/// * `end1` - End of first segment
/// * `start2` - Start of second segment
/// * `end2` - End of second segment
///
/// # Returns
/// True if segments overlap
pub fn segments_overlap(start1: i64, end1: i64, start2: i64, end2: i64) -> bool {
    !(end1 < start2 || end2 < start1)
}

/// Merge overlapping hits
///
/// # Arguments
/// * `hits` - Vector of hits to merge
///
/// # Returns
/// Vector of merged hits
pub fn merge_overlapping_hits(hits: &mut Vec<EufindHit>) -> Vec<EufindHit> {
    if hits.is_empty() {
        return Vec::new();
    }

    // Sort by position
    hits.sort_by_key(|h| (h.start, h.end));

    let mut merged = Vec::new();
    let mut current = hits[0].clone();

    for hit in hits.iter().skip(1) {
        if hit.strand == current.strand &&
           segments_overlap(current.start, current.end, hit.start, hit.end) {
            // Merge hits
            current.start = current.start.min(hit.start);
            current.end = current.end.max(hit.end);
            current.score = current.score.max(hit.score);
            current.isotype = hit.isotype.clone();
            current.anticodon = hit.anticodon.clone();
        } else {
            merged.push(current);
            current = hit.clone();
        }
    }
    merged.push(current);

    merged
}
