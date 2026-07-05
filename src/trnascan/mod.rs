//! tRNAscan first-pass scanner (Fichant-Burks algorithm).
//!
//! This module implements the tRNAscan algorithm as described in:
//! Fichant and Burks, J. Mol. Biol. (1991) 220:659-671.
//!
//! The algorithm scans sequences for tRNA genes by detecting
//! characteristic signal sequences (T-Psi-C and D signals) and
//! verifying stem-loop structures.

pub mod detect;
pub mod matrix;
pub mod params;
pub mod seq_utils;
pub mod signals;

pub use detect::{basepairing, read_signal, scoring, SignalResult};
pub use matrix::{load_consensus_matrix, ConsensusMatrix};
pub use params::{SearchParams, MAX_INTRON_LEN, MIN_SEQ_LEN, MIN_VAR_LOOP};
pub use seq_utils::{anticodon_to_aa, encode_anticodon, reverse_complement};

/// A tRNA gene hit from tRNAscan first-pass.
#[derive(Debug, Clone)]
pub struct TrnascanHit {
    /// Start position (1-based, on original strand)
    pub start: i64,
    /// End position (1-based, on original strand)
    pub end: i64,
    /// Strand: '+' for forward, '-' for reverse complement
    pub strand: char,
    /// Amino acid isotype (e.g., "Phe", "Leu", "Met")
    pub isotype: String,
    /// Anticodon sequence (3 nucleotides)
    pub anticodon: String,
    /// Position of D signal
    pub d_signal_pos: i64,
    /// Position of T-Psi-C signal
    pub tpc_signal_pos: i64,
    /// Whether an intron was detected
    pub has_intron: bool,
    /// Intron start position (0 if no intron)
    pub intron_start: i64,
    /// Intron end position (0 if no intron)
    pub intron_end: i64,
}

impl TrnascanHit {
    /// Create a new hit with default values.
    pub fn new() -> Self {
        Self {
            start: 0,
            end: 0,
            strand: '+',
            isotype: String::new(),
            anticodon: String::new(),
            d_signal_pos: 0,
            tpc_signal_pos: 0,
            has_intron: false,
            intron_start: 0,
            intron_end: 0,
        }
    }
}

impl Default for TrnascanHit {
    fn default() -> Self {
        Self::new()
    }
}

/// Scan a sequence for tRNA genes using the Fichant-Burks algorithm.
///
/// This is the main entry point for tRNA scanning. It searches both
/// strands of the input sequence for potential tRNA genes.
///
/// # Arguments
/// * `seq` - The DNA sequence (lowercase preferred)
/// * `tpc_matrix` - T-Psi-C signal consensus matrix
/// * `d_matrix` - D signal consensus matrix
/// * `params` - Search parameters controlling sensitivity
///
/// # Returns
/// Vector of tRNA hits found in the sequence
pub fn trnascan(
    seq: &[u8],
    tpc_matrix: &ConsensusMatrix,
    d_matrix: &ConsensusMatrix,
    params: &SearchParams,
) -> Vec<TrnascanHit> {
    let mut hits = Vec::new();
    let seqlen = seq.len();

    // Skip sequences shorter than minimum
    if seqlen < MIN_SEQ_LEN {
        return hits;
    }

    // Convert to lowercase for consistent comparison
    let seq_lower: Vec<u8> = seq.iter().map(|&b| b.to_ascii_lowercase()).collect();

    // Scan forward strand
    let forward_hits = scan_strand(&seq_lower, tpc_matrix, d_matrix, params, false);
    hits.extend(forward_hits);

    // Scan reverse complement strand
    let rc_seq = reverse_complement(&seq_lower);
    let reverse_hits = scan_strand(&rc_seq, tpc_matrix, d_matrix, params, true);

    // Convert reverse strand coordinates
    for mut hit in reverse_hits {
        // Flip coordinates for reverse strand
        let orig_start = seqlen as i64 - hit.end + 1;
        let orig_end = seqlen as i64 - hit.start + 1;
        hit.start = orig_start;
        hit.end = orig_end;

        // Also flip intron coordinates if present
        if hit.has_intron {
            let orig_intron_start = seqlen as i64 - hit.intron_end + 1;
            let orig_intron_end = seqlen as i64 - hit.intron_start + 1;
            hit.intron_start = orig_intron_start;
            hit.intron_end = orig_intron_end;
        }

        // Also flip signal positions
        let orig_d_pos = seqlen as i64 - hit.d_signal_pos + 1;
        let orig_tpc_pos = seqlen as i64 - hit.tpc_signal_pos + 1;
        hit.d_signal_pos = orig_d_pos;
        hit.tpc_signal_pos = orig_tpc_pos;

        hits.push(hit);
    }

    hits
}

/// Scan a single strand for tRNA genes.
fn scan_strand(
    seq: &[u8],
    tpc_matrix: &ConsensusMatrix,
    d_matrix: &ConsensusMatrix,
    params: &SearchParams,
    is_reverse: bool,
) -> Vec<TrnascanHit> {
    let mut hits = Vec::new();
    let seqlen = seq.len();

    if seqlen < MIN_SEQ_LEN {
        return hits;
    }

    // Search for T-Psi-C signal starting at position 44 (1-indexed)
    // In 0-indexed: start at 43
    let start_pos = 43;
    let end_pos = if seqlen > 23 { seqlen - 23 } else { return hits };

    for pos in start_pos..end_pos {
        let mut score = 0;

        // Check for T-Psi-C signal
        let tpc_result = read_signal(seq, pos, tpc_matrix, params.tpc_inv);

        if !tpc_result.found {
            continue;
        }

        // Increment score if enough invariant bases found
        if tpc_result.ninv >= tpc_matrix.ktot - 1 {
            score += 1;
        }

        // Score the T-Psi-C signal
        let (tpc_passes, _tpc_score) = scoring(
            &tpc_result.weight,
            tpc_matrix.lsig,
            tpc_matrix.maxtot,
            tpc_matrix.ktot,
            params.tpc_sig_thresh,
            tpc_result.ninv,
        );

        if !tpc_passes {
            continue;
        }

        // Check T-Psi-C stem (5 base pairs, lpair=16)
        // The stem starts at ptr+1 in C code, which is pos+1 here
        if pos + 17 >= seqlen {
            continue;
        }

        let ncomp = basepairing(&seq[pos + 1..], 5, 16);

        if ncomp >= params.tpc_incsg {
            score += 1;
        }

        if ncomp < params.tpc_keep {
            continue;
        }

        let first_score = score;

        // Search for D signal upstream (37-120 bp before TPC)
        let (d_begin, d_end) = if pos <= 127 {
            (7usize, pos.saturating_sub(37))
        } else {
            (
                pos.saturating_sub((MAX_INTRON_LEN + 60) as usize),
                pos.saturating_sub(37),
            )
        };

        for d_pos in d_begin..=d_end {
            let d_result = read_signal(seq, d_pos, d_matrix, params.d_inv);

            if !d_result.found {
                continue;
            }

            let mut current_score = first_score;

            // Increment score if all D invariants found
            if d_result.ninv >= d_matrix.ktot {
                current_score += 1;
            }

            // Score the D signal
            let (d_passes, _d_score) = scoring(
                &d_result.weight,
                d_matrix.lsig,
                d_matrix.maxtot,
                d_matrix.ktot,
                params.d_sig_thresh,
                d_result.ninv,
            );

            if !d_passes {
                continue;
            }

            // Check D stem with variable loop lengths (lpair 14-18)
            // ptr3 = ptr1 + 2 in C code
            let d_stem_start = d_pos + 2;
            if d_stem_start + 18 >= seqlen {
                continue;
            }

            // Try different D loop sizes
            let mut best_d_pairs = 0;
            let mut best_lpair = 14;

            for lpair in 14..=18 {
                if d_stem_start + lpair + 2 >= seqlen {
                    continue;
                }
                let d_ncomp = basepairing(&seq[d_stem_start..], 3, lpair as i32);
                if d_ncomp > best_d_pairs {
                    best_d_pairs = d_ncomp;
                    best_lpair = lpair;
                }
            }

            // Continue with following_search logic
            if let Some(hit) = following_search(
                seq,
                seqlen,
                pos,
                d_pos,
                best_lpair as i32,
                best_d_pairs == 3,
                current_score,
                params,
                is_reverse,
            ) {
                hits.push(hit);
            }
        }
    }

    hits
}

/// Complete the tRNA search after finding TPC and D signals.
///
/// This implements the following_search function from trnascan.c.
fn following_search(
    seq: &[u8],
    seqlen: usize,
    tpc_pos: usize,
    d_pos: usize,
    d_lpair: i32,
    d_has_3_pairs: bool,
    initial_score: i32,
    params: &SearchParams,
    is_reverse: bool,
) -> Option<TrnascanHit> {
    let mut score = initial_score;

    // If D arm has 3 base pairings, increment score
    if d_has_3_pairs {
        score += 1;
    }

    // Check aminoacyl stem
    // npair1 = 7
    // lpair1 = pos - pos1 + 8 + 23
    // ptr2 = ptr1 - 7 (start of aa stem)
    let lpair1 = (tpc_pos as i64 - d_pos as i64 + 8 + 23) as i32;

    if d_pos < 7 {
        return None;
    }

    let aa_start = d_pos - 7;
    if aa_start + lpair1 as usize + 2 > seqlen {
        return None;
    }

    let aa_ncomp = basepairing(&seq[aa_start..], 7, lpair1);

    if aa_ncomp >= params.aa_incsg {
        score += 1;
    }

    if aa_ncomp < params.aa_keep {
        return None;
    }

    // Look for anticodon stem if score is high enough
    if score < params.look_for_acloop_sg {
        return None;
    }

    // Anticodon stem position: ptr4 = ptr3 + lpair + 2
    // ptr3 = d_pos + 2
    let ac_start = d_pos + 2 + d_lpair as usize + 2;
    if ac_start + 20 > seqlen {
        return None;
    }

    // Check anticodon stem (5 pairs, lpair2=16 for no intron)
    let ac_ncomp = basepairing(&seq[ac_start..], 5, 16);

    let mut final_score = score;
    let mut lpair2 = 16;

    if ac_ncomp >= params.acloop_min {
        // Check for invariant T before anticodon
        if ac_start + 6 < seqlen && seq[ac_start + 6] == b't' {
            final_score += 1;
        }

        // Check if score passes threshold
        if final_score >= params.sg_cutoff {
            // Extract anticodon
            if ac_start + 9 < seqlen {
                let anticodon: String =
                    String::from_utf8_lossy(&seq[ac_start + 7..ac_start + 10]).to_string();
                let ac_code = encode_anticodon(&seq[ac_start + 7..ac_start + 10]);
                let isotype = anticodon_to_aa(ac_code).to_string();

                let hit = create_hit(
                    seq,
                    seqlen,
                    tpc_pos,
                    d_pos,
                    d_lpair,
                    lpair1,
                    lpair2,
                    is_reverse,
                    isotype,
                    anticodon,
                    false,
                    0,
                    0,
                );

                return Some(hit);
            }
        }
    }

    // Try to find anticodon with intron
    let pos6 = tpc_pos as i32 - (d_pos as i32 + d_lpair + 3);

    if pos6 >= MIN_VAR_LOOP {
        // Search for anticodon stem with intron
        for lpair2_try in (MIN_VAR_LOOP - 4) as usize..pos6 as usize {
            let mut score2 = score;

            if ac_start + lpair2_try >= seqlen {
                continue;
            }

            let ac_ncomp2 = basepairing(&seq[ac_start..], 5, lpair2_try as i32);

            // Variable loop size check
            let pos4 = tpc_pos as i32 - d_pos as i32 - d_lpair - lpair2_try as i32 - 4;

            if ac_ncomp2 >= 4 && pos4 >= 3 {
                // Check for A or G before intron
                let ptr5_pos = ac_start + 10;
                if ptr5_pos < seqlen && (seq[ptr5_pos] == b'a' || seq[ptr5_pos] == b'g') {
                    // Check invariant T before anticodon
                    if ac_start + 6 < seqlen && seq[ac_start + 6] == b't' {
                        score2 += 1;
                    }

                    if score2 >= params.sg_cutoff {
                        lpair2 = lpair2_try as i32;

                        // Calculate intron positions (1-based)
                        let intron_start_1based = d_pos as i64 + d_lpair as i64 + 15 + 1;
                        let intron_end_1based = d_pos as i64 + d_lpair as i64 + lpair2 as i64 - 2 + 1;

                        // Extract anticodon
                        if ac_start + 9 < seqlen {
                            let anticodon: String =
                                String::from_utf8_lossy(&seq[ac_start + 7..ac_start + 10])
                                    .to_string();
                            let ac_code = encode_anticodon(&seq[ac_start + 7..ac_start + 10]);
                            let isotype = anticodon_to_aa(ac_code).to_string();

                            let hit = create_hit(
                                seq,
                                seqlen,
                                tpc_pos,
                                d_pos,
                                d_lpair,
                                lpair1,
                                lpair2,
                                is_reverse,
                                isotype,
                                anticodon,
                                true,
                                intron_start_1based,
                                intron_end_1based,
                            );

                            return Some(hit);
                        }
                    }
                }
            }
        }
    }

    None
}

/// Create a TrnascanHit from detected positions.
fn create_hit(
    _seq: &[u8],
    _seqlen: usize,
    tpc_pos: usize,
    d_pos: usize,
    _d_lpair: i32,
    lpair1: i32,
    _lpair2: i32,
    is_reverse: bool,
    isotype: String,
    anticodon: String,
    has_intron: bool,
    intron_start: i64,
    intron_end: i64,
) -> TrnascanHit {
    // Start position: d_pos - 7 (0-indexed) -> d_pos - 6 (1-indexed)
    // End position: d_pos - 7 + lpair1 + 1 (0-indexed) -> d_pos - 6 + lpair1 + 1 (1-indexed)
    let start = (d_pos as i64 - 7) + 1; // Convert to 1-based
    let end = start + lpair1 as i64;

    TrnascanHit {
        start,
        end,
        strand: if is_reverse { '-' } else { '+' },
        isotype,
        anticodon,
        d_signal_pos: (d_pos as i64) + 1,
        tpc_signal_pos: (tpc_pos as i64) + 1,
        has_intron,
        intron_start,
        intron_end,
    }
}

/// Convenience function to scan with default matrices.
pub fn trnascan_default(seq: &[u8], params: &SearchParams) -> Vec<TrnascanHit> {
    let tpc_matrix = ConsensusMatrix::tpc_signal();
    let d_matrix = ConsensusMatrix::d_signal();
    trnascan(seq, &tpc_matrix, &d_matrix, params)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_trnascan_hit_default() {
        let hit = TrnascanHit::default();
        assert_eq!(hit.start, 0);
        assert_eq!(hit.strand, '+');
        assert!(!hit.has_intron);
    }

    #[test]
    fn test_short_sequence() {
        let seq = b"acgtacgt";
        let params = SearchParams::relaxed();
        let hits = trnascan_default(seq, &params);
        assert!(hits.is_empty());
    }

    #[test]
    fn test_matrices_loaded() {
        let tpc = ConsensusMatrix::tpc_signal();
        let d = ConsensusMatrix::d_signal();

        assert_eq!(tpc.lsig, 15);
        assert_eq!(d.lsig, 8);
        assert!(tpc.maxtot > 0.0);
        assert!(d.maxtot > 0.0);
    }

    // A proper test would need a real tRNA sequence
    // This is a placeholder showing the API
    #[test]
    fn test_scan_random_sequence() {
        // Random sequence unlikely to contain tRNA
        let seq: Vec<u8> = (0..200)
            .map(|i| match i % 4 {
                0 => b'a',
                1 => b'c',
                2 => b'g',
                _ => b't',
            })
            .collect();

        let params = SearchParams::strict();
        let hits = trnascan_default(&seq, &params);

        // Random sequence should not have tRNA-like structure
        // (though it's possible to get false positives)
        // This test just ensures the function runs without panicking
        assert!(hits.len() < 10); // Sanity check
    }
}
