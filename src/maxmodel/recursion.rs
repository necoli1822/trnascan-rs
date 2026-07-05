//! Recursion calculations for MaxModelMaker algorithm
//!
//! This module implements the core DP recursion from maxmodelmaker.c
//! recurse_maxmx function (lines 435-813).

use crate::maxmodel::emission::{dot_score_2d, pair_emissioncost};
use crate::maxmodel::matrix::MaxMx;
use crate::maxmodel::prior::Prior;
use crate::maxmodel::transition::*;
use crate::maxmodel::types::*;
use crate::types::constants::*;

/// Recursion calculations of the maximum likelihood CM construction algorithm
///
/// This is the heart of the MaxModelMaker algorithm. It fills in the scoring
/// matrix by considering all possible ways to explain each subsequence (i,j).
///
/// From maxmodelmaker.c recurse_maxmx (lines 435-813)
///
/// # Arguments
/// * `aseqs_t` - Transposed alignment [1..alen+1][0..nseq-1]
/// * `weights` - Weights on sequences
/// * `prior` - Prior probability distributions
/// * `mscore` - Singlet match emission scores
/// * `gapcount` - Weighted gap counts per column
/// * `gapthresh` - Fractional occupancy threshold (scaled by nseq)
/// * `mmx` - Scoring matrix to fill
pub fn recurse_maxmx(
    aseqs_t: &[Vec<i8>],
    weights: &[f32],
    prior: &Prior,
    mscore: &[i32],
    gapcount: &[f64],
    gapthresh: f64,
    mmx: &mut MaxMx,
) {
    let alen = mmx.alen;
    let nseq = weights.len();

    // Scale gapthresh to be comparable to counts in gapcount array
    let gapthresh_scaled = gapthresh * nseq as f64;

    // Allocate insertion accumulators
    // insr_accum[j2][idx] = count of insertions between j2 and current j for seq idx
    // insl_accum[i2][idx] = count of insertions between current i and i2 for seq idx
    let mut insr_accum: Vec<Vec<i32>> = vec![vec![0; nseq]; alen + 1];
    let mut insl_accum: Vec<Vec<i32>> = vec![vec![0; nseq]; alen + 2];

    // Initialize insr_accum
    for j2 in 1..=alen {
        for idx in 0..nseq {
            insr_accum[j2][idx] = 0;
        }
    }
    for idx in 0..nseq {
        if aseqs_t[1][idx] >= 0 {
            insr_accum[0][idx] = 1;
        }
    }

    // Main recursion: for each row j
    for j in 2..=alen {
        // Initialize insl_accum for this row
        for i2 in 1..=(j + 1) {
            for idx in 0..nseq {
                insl_accum[i2][idx] = 0;
            }
        }
        for idx in 0..nseq {
            if aseqs_t[j][idx] >= 0 {
                insl_accum[j + 1][idx] += 1;
            }
        }

        // For each column i (from j-1 down to 1)
        for i in (1..j).rev() {
            // === BIFURC ===
            // Explain i,j as sum of i,mid,BEGINL + mid+1,j,BEGINR
            for mid in i..=j {
                let sc = mmx.get(mid, i).sc[maxmx_node::BEGINL]
                    + mmx.get(j, mid + 1).sc[maxmx_node::BEGINR];
                if sc > mmx.get(j, i).sc[maxmx_node::BIFURC] {
                    mmx.get_mut(j, i).sc[maxmx_node::BIFURC] = sc;
                    mmx.get_mut(j, i).bifurc_mid = mid as i16;
                }
            }

            // === MATP ===
            // Score subsequence i,j given that i,j are emitted by MATP
            let mut matp_done = false;
            for j2 in (i..j).rev() {
                if (j - j2 - 1) > MAXINSERT {
                    break;
                }

                for i2 in (i + 1)..=(j2 + 1) {
                    if (i2 - i - 1) > MAXINSERT {
                        break;
                    }

                    // Early pruning: check if any target score could improve
                    let curr_matp = mmx.get(j, i).sc[maxmx_node::MATP];
                    if mmx.get(j2, i2).sc[maxmx_node::MATP] < curr_matp
                        && mmx.get(j2, i2).sc[maxmx_node::MATL] < curr_matp
                        && mmx.get(j2, i2).sc[maxmx_node::MATR] < curr_matp
                        && mmx.get(j2, i2).sc[maxmx_node::BIFURC] < curr_matp
                    {
                        continue;
                    }

                    let mut tmaster = zero_transtable();
                    frommatp_transtable(
                        aseqs_t,
                        weights,
                        i,
                        j,
                        i2,
                        j2,
                        &insl_accum[i2],
                        &insr_accum[j2],
                        &mut tmaster,
                    );

                    for tonode in 0..4 {
                        if i2 > j2 && tonode != maxmx_node::BIFURC {
                            continue;
                        }
                        if i2 == j2 && tonode == maxmx_node::MATP {
                            continue;
                        }
                        if mmx.get(j2, i2).sc[tonode] < mmx.get(j, i).sc[maxmx_node::MATP] {
                            continue;
                        }

                        let mut tcounts = zero_transtable();
                        match tonode {
                            maxmx_node::MATP => to_matp_transtable(&tmaster, &mut tcounts),
                            maxmx_node::MATL => to_matl_transtable(&tmaster, &mut tcounts),
                            maxmx_node::MATR => to_matr_transtable(&tmaster, &mut tcounts),
                            maxmx_node::BIFURC => to_bifurc_transtable(&tmaster, &mut tcounts),
                            _ => continue,
                        }

                        let mut tmx = [[0.0; STATETYPES]; STATETYPES];
                        copy_transtable(&mut tmx, &tcounts);
                        prior.probify_transition_matrix(&mut tmx, MATP_NODE, tonode);

                        let sc = dot_score_2d(&tcounts, &tmx) + mmx.get(j2, i2).sc[tonode];

                        if sc > mmx.get(j, i).sc[maxmx_node::MATP] {
                            mmx.get_mut(j, i).sc[maxmx_node::MATP] = sc;
                            mmx.get_mut(j, i).matp_i2 = i2 as i16;
                            mmx.get_mut(j, i).matp_j2 = j2 as i16;
                            mmx.get_mut(j, i).matp_ftype = tonode as u8;
                        }
                    }

                    if gapcount[i2] <= gapthresh_scaled {
                        break;
                    }
                }

                if gapcount[j2] <= gapthresh_scaled {
                    matp_done = true;
                    break;
                }
            }

            // Add pair emission cost
            if !matp_done || mmx.get(j, i).sc[maxmx_node::MATP] > NEGINFINITY {
                mmx.get_mut(j, i).sc[maxmx_node::MATP] +=
                    pair_emissioncost(&aseqs_t[i], &aseqs_t[j], weights, prior);
            }

            // === MATR ===
            // Account for i,j by emitting j and connecting to some i,j2
            for j2 in (i.saturating_sub(1)..j).rev() {
                let curr_matr = mmx.get(j, i).sc[maxmx_node::MATR];
                if mmx.get(j2, i).sc[maxmx_node::MATP] < curr_matr
                    && mmx.get(j2, i).sc[maxmx_node::MATL] < curr_matr
                    && mmx.get(j2, i).sc[maxmx_node::MATR] < curr_matr
                    && mmx.get(j2, i).sc[maxmx_node::BIFURC] < curr_matr
                {
                    continue;
                }

                let mut tmaster = zero_transtable();
                frommatr_transtable(aseqs_t, weights, i, j, j2, &insr_accum[j2], &mut tmaster);

                for tonode in 0..4 {
                    if i > j2 && tonode != maxmx_node::BIFURC {
                        continue;
                    }
                    if i == j2 && tonode == maxmx_node::MATP {
                        continue;
                    }
                    if mmx.get(j2, i).sc[tonode] < mmx.get(j, i).sc[maxmx_node::MATR] {
                        continue;
                    }

                    let mut tcounts = zero_transtable();
                    match tonode {
                        maxmx_node::MATP => to_matp_transtable(&tmaster, &mut tcounts),
                        maxmx_node::MATL => to_matl_transtable(&tmaster, &mut tcounts),
                        maxmx_node::MATR => to_matr_transtable(&tmaster, &mut tcounts),
                        maxmx_node::BIFURC => to_bifurc_transtable(&tmaster, &mut tcounts),
                        _ => continue,
                    }

                    let mut tmx = [[0.0; STATETYPES]; STATETYPES];
                    copy_transtable(&mut tmx, &tcounts);
                    prior.probify_transition_matrix(&mut tmx, MATR_NODE, tonode);

                    let sc = dot_score_2d(&tcounts, &tmx) + mmx.get(j2, i).sc[tonode];

                    if sc > mmx.get(j, i).sc[maxmx_node::MATR] {
                        mmx.get_mut(j, i).sc[maxmx_node::MATR] = sc;
                        mmx.get_mut(j, i).matr_j2 = j2 as i16;
                        mmx.get_mut(j, i).matr_ftype = tonode as u8;
                    }
                }

                if gapcount[j2] <= gapthresh_scaled {
                    break;
                }
            }
            mmx.get_mut(j, i).sc[maxmx_node::MATR] += mscore[j];

            // === MATL ===
            // Account for i,j by emitting i and connecting to some i2,j
            for i2 in (i + 1)..=(j + 1) {
                let curr_matl = mmx.get(j, i).sc[maxmx_node::MATL];
                if mmx.get(j, i2).sc[maxmx_node::MATP] < curr_matl
                    && mmx.get(j, i2).sc[maxmx_node::MATL] < curr_matl
                    && mmx.get(j, i2).sc[maxmx_node::MATR] < curr_matl
                    && mmx.get(j, i2).sc[maxmx_node::BIFURC] < curr_matl
                {
                    continue;
                }

                let mut tmaster = zero_transtable();
                frommatl_transtable(aseqs_t, weights, i, j, i2, &insl_accum[i2], &mut tmaster);

                for tonode in 0..4 {
                    if i2 > j && tonode != maxmx_node::BIFURC {
                        continue;
                    }
                    if i2 == j && tonode == maxmx_node::MATP {
                        continue;
                    }
                    if mmx.get(j, i2).sc[tonode] < mmx.get(j, i).sc[maxmx_node::MATL] {
                        continue;
                    }

                    let mut tcounts = zero_transtable();
                    match tonode {
                        maxmx_node::MATP => to_matp_transtable(&tmaster, &mut tcounts),
                        maxmx_node::MATL => to_matl_transtable(&tmaster, &mut tcounts),
                        maxmx_node::MATR => to_matr_transtable(&tmaster, &mut tcounts),
                        maxmx_node::BIFURC => to_bifurc_transtable(&tmaster, &mut tcounts),
                        _ => continue,
                    }

                    let mut tmx = [[0.0; STATETYPES]; STATETYPES];
                    copy_transtable(&mut tmx, &tcounts);
                    prior.probify_transition_matrix(&mut tmx, MATL_NODE, tonode);

                    let sc = dot_score_2d(&tcounts, &tmx) + mmx.get(j, i2).sc[tonode];

                    if sc > mmx.get(j, i).sc[maxmx_node::MATL] {
                        mmx.get_mut(j, i).sc[maxmx_node::MATL] = sc;
                        mmx.get_mut(j, i).matl_i2 = i2 as i16;
                        mmx.get_mut(j, i).matl_ftype = tonode as u8;
                    }
                }

                if gapcount[i2] <= gapthresh_scaled {
                    break;
                }
            }
            mmx.get_mut(j, i).sc[maxmx_node::MATL] += mscore[i];

            // Bump insl_accum: add column i to horizontal accumulator
            for i2 in (i + 1)..=(j + 1) {
                for idx in 0..nseq {
                    if aseqs_t[i][idx] >= 0 {
                        insl_accum[i2][idx] += 1;
                    }
                }
            }

            // === BEGINR ===
            // Has INSL state, so can connect to any i2,j inclusive of i,j
            for i2 in i..=(j + 1) {
                let curr_begr = mmx.get(j, i).sc[maxmx_node::BEGINR];
                if mmx.get(j, i2).sc[maxmx_node::MATP] < curr_begr
                    && mmx.get(j, i2).sc[maxmx_node::MATL] < curr_begr
                    && mmx.get(j, i2).sc[maxmx_node::MATR] < curr_begr
                    && mmx.get(j, i2).sc[maxmx_node::BIFURC] < curr_begr
                {
                    continue;
                }

                let mut tmaster = zero_transtable();
                frombeginr_transtable(aseqs_t, weights, j, i2, &insl_accum[i2], &mut tmaster);

                for tonode in 0..4 {
                    if i2 > j && tonode != maxmx_node::BIFURC {
                        continue;
                    }
                    if i2 == j && tonode == maxmx_node::MATP {
                        continue;
                    }
                    if mmx.get(j, i2).sc[tonode] < mmx.get(j, i).sc[maxmx_node::BEGINR] {
                        continue;
                    }

                    let mut tcounts = zero_transtable();
                    match tonode {
                        maxmx_node::MATP => to_matp_transtable(&tmaster, &mut tcounts),
                        maxmx_node::MATL => to_matl_transtable(&tmaster, &mut tcounts),
                        maxmx_node::MATR => to_matr_transtable(&tmaster, &mut tcounts),
                        maxmx_node::BIFURC => to_bifurc_transtable(&tmaster, &mut tcounts),
                        _ => continue,
                    }

                    let mut tmx = [[0.0; STATETYPES]; STATETYPES];
                    copy_transtable(&mut tmx, &tcounts);
                    prior.probify_transition_matrix(&mut tmx, BEGINR_NODE, tonode);

                    let sc = dot_score_2d(&tcounts, &tmx) + mmx.get(j, i2).sc[tonode];

                    if sc > mmx.get(j, i).sc[maxmx_node::BEGINR] {
                        mmx.get_mut(j, i).sc[maxmx_node::BEGINR] = sc;
                        mmx.get_mut(j, i).begr_i2 = i2 as i16;
                        mmx.get_mut(j, i).begr_ftype = tonode as u8;
                    }
                }

                if gapcount[i2] <= gapthresh_scaled {
                    break;
                }
            }

            // === BEGINL ===
            // No inserts, must connect to i,j
            let mut tmaster = zero_transtable();
            frombeginl_transtable(aseqs_t, weights, i, j, &mut tmaster);

            for tonode in 0..4 {
                let curr_begl = mmx.get(j, i).sc[maxmx_node::BEGINL];
                if mmx.get(j, i).sc[maxmx_node::MATP] < curr_begl
                    && mmx.get(j, i).sc[maxmx_node::MATL] < curr_begl
                    && mmx.get(j, i).sc[maxmx_node::MATR] < curr_begl
                    && mmx.get(j, i).sc[maxmx_node::BIFURC] < curr_begl
                {
                    continue;
                }

                let mut tcounts = zero_transtable();
                match tonode {
                    maxmx_node::MATP => to_matp_transtable(&tmaster, &mut tcounts),
                    maxmx_node::MATL => to_matl_transtable(&tmaster, &mut tcounts),
                    maxmx_node::MATR => to_matr_transtable(&tmaster, &mut tcounts),
                    maxmx_node::BIFURC => to_bifurc_transtable(&tmaster, &mut tcounts),
                    _ => continue,
                }

                let mut tmx = [[0.0; STATETYPES]; STATETYPES];
                copy_transtable(&mut tmx, &tcounts);
                prior.probify_transition_matrix(&mut tmx, BEGINL_NODE, tonode);

                let sc = dot_score_2d(&tcounts, &tmx) + mmx.get(j, i).sc[tonode];

                if sc > mmx.get(j, i).sc[maxmx_node::BEGINL] {
                    mmx.get_mut(j, i).sc[maxmx_node::BEGINL] = sc;
                    mmx.get_mut(j, i).begl_ftype = tonode as u8;
                }
            }
        } // end loop over columns i

        // Bump insr_accum: add row j to vertical accumulator
        for j2 in 0..j {
            for idx in 0..nseq {
                if aseqs_t[j][idx] >= 0 {
                    insr_accum[j2][idx] += 1;
                }
            }
        }
    } // end loop over rows j

    // === Termination: ROOT ===
    // ROOT can connect anywhere. Store in mmx[alen][0] using MATP_NODE slot
    for j2 in (0..=alen).rev() {
        for i2 in 1..=(j2 + 1) {
            let curr_root = mmx.get(alen, 0).sc[maxmx_node::MATP];
            if mmx.get(j2, i2).sc[maxmx_node::MATP] < curr_root
                && mmx.get(j2, i2).sc[maxmx_node::MATL] < curr_root
                && mmx.get(j2, i2).sc[maxmx_node::MATR] < curr_root
                && mmx.get(j2, i2).sc[maxmx_node::BIFURC] < curr_root
            {
                continue;
            }

            let mut tmaster = zero_transtable();
            fromroot_transtable(
                aseqs_t,
                weights,
                i2,
                j2,
                &insl_accum[i2],
                &insr_accum[j2],
                &mut tmaster,
            );

            for tonode in 0..4 {
                if i2 > j2 && tonode != maxmx_node::BIFURC {
                    continue;
                }
                if i2 == j2 && tonode == maxmx_node::MATP {
                    continue;
                }
                if mmx.get(j2, i2).sc[tonode] < mmx.get(alen, 0).sc[maxmx_node::MATP] {
                    continue;
                }

                let mut tcounts = zero_transtable();
                match tonode {
                    maxmx_node::MATP => to_matp_transtable(&tmaster, &mut tcounts),
                    maxmx_node::MATL => to_matl_transtable(&tmaster, &mut tcounts),
                    maxmx_node::MATR => to_matr_transtable(&tmaster, &mut tcounts),
                    maxmx_node::BIFURC => to_bifurc_transtable(&tmaster, &mut tcounts),
                    _ => continue,
                }

                let mut tmx = [[0.0; STATETYPES]; STATETYPES];
                copy_transtable(&mut tmx, &tcounts);
                prior.probify_transition_matrix(&mut tmx, ROOT_NODE, tonode);

                let sc = dot_score_2d(&tcounts, &tmx) + mmx.get(j2, i2).sc[tonode];

                if sc > mmx.get(alen, 0).sc[maxmx_node::MATP] {
                    mmx.get_mut(alen, 0).sc[maxmx_node::MATP] = sc;
                    mmx.get_mut(alen, 0).matp_i2 = i2 as i16;
                    mmx.get_mut(alen, 0).matp_j2 = j2 as i16;
                    mmx.get_mut(alen, 0).matp_ftype = tonode as u8;
                }
            }

            if gapcount[i2] <= gapthresh_scaled {
                break;
            }
        }

        if gapcount[j2] <= gapthresh_scaled {
            break;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::maxmodel::emission::{is_gap, singlet_emissions, symbol_index, transpose_alignment};

    #[test]
    fn test_recurse_maxmx_simple() {
        // Simple 2-sequence, 4-column alignment
        let seq1: &[u8] = b"ACGT";
        let seq2: &[u8] = b"ACGT";
        let aseqs: Vec<&[u8]> = vec![seq1, seq2];
        let weights = vec![1.0f32, 1.0];
        let prior = Prior::new();

        let aseqs_t = transpose_alignment(&aseqs, 4, is_gap, symbol_index);
        let (mscore, gapcount) = singlet_emissions(&aseqs_t, &weights, &prior);

        let mut mmx = MaxMx::new(4);
        mmx.initialize(2, &prior, &mscore, &gapcount);

        // Run recursion
        recurse_maxmx(&aseqs_t, &weights, &prior, &mscore, &gapcount, 0.5, &mut mmx);

        // Check that ROOT score was calculated (stored in mmx[alen][0].sc[MATP])
        let root_score = mmx.get(4, 0).sc[maxmx_node::MATP];
        assert!(root_score > NEGINFINITY);
    }
}
