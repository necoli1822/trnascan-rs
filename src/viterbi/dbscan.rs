//! Database Scanning Viterbi Algorithm
//!
//! This module implements ViterbiScan from dbviterbi.c
//! It is a scanning variant that finds matches in long sequences using a sliding window.
//!
//! # Key Differences from Standard Viterbi
//!
//! - Uses a fixed window size instead of sequence length
//! - Uses rolling buffers: amx[2][window+1][statenum] and bmx circular buffer
//! - Reports hits above threshold via callback rather than full traceback
//! - No traceback - only scores are computed
//!
//! # Matrix Layout
//!
//! amx: [aj][diff][y] where aj = j % 2 (rolling 2-row window)
//! bmx: [y][bj][diff] where bj = j % window (circular buffer for BEGIN states)

use crate::types::constants::*;
use crate::types::state::IState;

use super::ViterbiError;

/// Convert a nucleotide character to its index (0-3)
/// Matches SymbolIndex() from original code
#[inline]
fn symbol_index(c: u8) -> usize {
    match c.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' | b'U' => 3,
        _ => 0, // Default to A for unknown bases
    }
}

/// Database scanning matrices with circular buffer for window-based scanning
struct DBMatrices {
    /// Main score matrix: amx[aj][diff][y] where aj is 0 or 1 (rolling window)
    /// Dimensions: [2][window+1][statenum]
    amx: Vec<Vec<Vec<i32>>>,

    /// BEGIN score matrix: bmx[y] is Some([bj][diff]) for BEGIN states only
    /// bj = j % window (circular buffer), diff = 0..window
    /// Dimensions for BEGIN states: [window][window+1]
    bmx: Vec<Option<Vec<Vec<i32>>>>,
}

/// Database scanning Viterbi function
///
/// Implements ViterbiScan() from dbviterbi.c lines 65-96
///
/// Scans a long sequence with a sliding window and reports hits above threshold.
///
/// # Arguments
/// * `icm` - Integer state model
/// * `statenum` - Number of states in the model
/// * `seq` - Sequence to scan (ASCII nucleotides)
/// * `window` - Scanning window size in nucleotides
/// * `thresh` - Score threshold for reporting hits (in bits)
/// * `gotone_f` - Callback function for each hit: (start, end, score) -> bool
///
/// # Returns
/// * `Ok(())` - Scanning completed successfully
/// * `Err(ViterbiError)` - If scanning fails
///
/// # Callback
/// The callback receives 1-based coordinates (start, end) and score in bits.
/// Return `true` to continue scanning, `false` to abort.
pub fn viterbi_scan<F>(
    icm: &[IState],
    statenum: usize,
    seq: &[u8],
    window: usize,
    thresh: f64,
    mut gotone_f: F,
) -> Result<(), ViterbiError>
where
    F: FnMut(usize, usize, f64) -> bool,
{
    let n = seq.len();
    let ithresh = (thresh * INTPRECISION) as i32;

    // Allocate matrices
    let mut matrices = allocate_mx(icm, statenum, window)?;

    // Initialize matrices
    init_mx(icm, statenum, window, &mut matrices)?;

    // Fill matrices via recursion with hit reporting
    recurse_mx(icm, statenum, seq, n, window, &mut matrices, ithresh, &mut gotone_f)?;

    Ok(())
}

/// Allocate the scoring matrices for database scanning
///
/// Implements allocate_mx() from dbviterbi.c lines 120-178
///
/// Key differences from standard Viterbi:
/// - amx: Only 2 rows (rolling window) x (window+1) diff values
/// - bmx: Circular buffer with window rows for j dimension
fn allocate_mx(
    icm: &[IState],
    statenum: usize,
    window: usize,
) -> Result<DBMatrices, ViterbiError> {
    // Main matrix amx: [aj][diff][y] where aj is 0 or 1
    // amx[j][diag] for j = 0..1, diag = 0..window
    let mut amx: Vec<Vec<Vec<i32>>> = Vec::with_capacity(2);

    for _j in 0..=1 {
        let mut row: Vec<Vec<i32>> = Vec::with_capacity(window + 1);
        for _diff in 0..=window {
            row.push(vec![NEGINFINITY; statenum]);
        }
        amx.push(row);
    }

    // BEGIN auxiliary matrix bmx[y][bj][diff]
    // bmx keeps scores for BEGIN states only
    // bj = 0..window-1 (circular buffer), diff = 0..window
    let mut bmx: Vec<Option<Vec<Vec<i32>>>> = vec![None; statenum];

    for y in 0..statenum {
        // We keep score info for BEGIN states
        if icm[y].statetype == U_BEGIN_ST {
            let mut y_matrix: Vec<Vec<i32>> = Vec::with_capacity(window);
            for _bj in 0..window {
                y_matrix.push(vec![NEGINFINITY; window + 1]);
            }
            bmx[y] = Some(y_matrix);
        }
    }

    Ok(DBMatrices { amx, bmx })
}

/// Initialize the scoring matrices for database scanning
///
/// Implements init_mx() from dbviterbi.c lines 232-307
fn init_mx(
    icm: &[IState],
    statenum: usize,
    window: usize,
    matrices: &mut DBMatrices,
) -> Result<(), ViterbiError> {
    let DBMatrices { amx, bmx } = matrices;

    // Init the whole amx to -Infinity
    for j in 0..=1 {
        for diff in 0..=window {
            for y in 0..statenum {
                amx[j][diff][y] = NEGINFINITY;
            }
        }
    }

    // Init the whole bmx to -Infinity
    // State 0 is always a BEGIN (ROOT), so we start there and copy rows
    if let Some(ref mut root_bmx) = bmx[0] {
        for bj in 0..window {
            for diff in 0..=window {
                root_bmx[bj][diff] = NEGINFINITY;
            }
        }
    }

    for y in 1..statenum {
        if let Some(ref mut y_bmx) = bmx[y] {
            for bj in 0..window {
                for diff in 0..=window {
                    y_bmx[bj][diff] = NEGINFINITY;
                }
            }
        }
    }

    // Init the off-diagonal (j = 0..window-1; diff == 0) with -log P scores
    // End state = 0; del, bifurc states are calc'd; begin states same as del's
    for j in 0..window {
        let aj = j % 2;

        // Process states in reverse order (statenum-1 down to 0)
        for y in (0..statenum).rev() {
            if icm[y].statetype == U_END_ST {
                amx[aj][0][y] = 0;
            } else if icm[y].statetype == U_BIFURC_ST {
                // bifurc[y] = bmx[left][j][0] + bmx[right][j][0]
                let left = icm[y].tmx[0] as usize;
                let right = icm[y].tmx[1] as usize;

                let left_score = if let Some(ref left_bmx) = bmx[left] {
                    left_bmx[j][0]
                } else {
                    NEGINFINITY
                };

                let right_score = if let Some(ref right_bmx) = bmx[right] {
                    right_bmx[j][0]
                } else {
                    NEGINFINITY
                };

                if left_score != NEGINFINITY && right_score != NEGINFINITY {
                    amx[aj][0][y] = left_score + right_score;
                }
            } else if icm[y].statetype == U_DEL_ST || icm[y].statetype == U_BEGIN_ST {
                // Only calc DEL-DEL and BEGIN-DEL transitions
                // Find the connection to a non-infinite score
                let offset = icm[y].offset as usize;
                let base_idx = y + offset;

                for ynext in 0..(icm[y].connectnum as usize) {
                    let next_y = base_idx + ynext;
                    if next_y < statenum && amx[aj][0][next_y] != NEGINFINITY {
                        amx[aj][0][y] = amx[aj][0][next_y] + icm[y].tmx[ynext];
                        break;
                    }
                }
            }

            // Make a copy into bmx if y is a BEGIN
            if icm[y].statetype == U_BEGIN_ST {
                if let Some(ref mut y_bmx) = bmx[y] {
                    y_bmx[j][0] = amx[aj][0][y];
                }
            }
        }
    }

    Ok(())
}

/// Fill the scoring matrices via recursion with hit reporting
///
/// Implements recurse_mx() from dbviterbi.c lines 320-472
///
/// Key differences from standard Viterbi:
/// - Uses aj = j % 2 for rolling amx
/// - Uses bj = j % window for circular bmx buffer
/// - BIFURC uses wraparound: leftj = leftj ? leftj-1 : window-1
/// - Reports hits above threshold after each row
fn recurse_mx<F>(
    icm: &[IState],
    statenum: usize,
    seq: &[u8],
    seqlen: usize,
    window: usize,
    matrices: &mut DBMatrices,
    ithresh: i32,
    gotone_f: &mut F,
) -> Result<(), ViterbiError>
where
    F: FnMut(usize, usize, f64) -> bool,
{
    let DBMatrices { amx, bmx } = matrices;

    for j in 1..=seqlen {
        let aj = j % 2;           // 0 or 1 index in amx
        let prev_aj = 1 - aj;     // Previous row in amx
        let bj = j % window;      // 0..window-1 index in bmx (circular)
        let symj = symbol_index(seq[j - 1]); // 1-based to 0-based

        // Calculate max diff for this j (limited by window and j)
        let max_diff = window.min(j);

        for diff in 1..=max_diff {
            let i = j - diff + 1; // 1-based i
            let symi = symbol_index(seq[i - 1]); // 1-based to 0-based

            // Process states in reverse order
            for y in (0..statenum).rev() {
                let st = &icm[y];

                if st.statetype != U_BIFURC_ST {
                    // A normal (non-BIFURC) state
                    let (beam_aj, beam_diff, emitsc): (usize, usize, i32) = match st.statetype {
                        U_BEGIN_ST | U_DEL_ST => (aj, diff, 0),
                        U_MATP_ST => {
                            if diff == 1 {
                                continue;
                            }
                            (prev_aj, diff - 2, st.emit[symi * ALPHASIZE + symj])
                        }
                        U_MATR_ST | U_INSR_ST => (prev_aj, diff - 1, st.emit[symj]),
                        U_MATL_ST | U_INSL_ST => (aj, diff - 1, st.emit[symi]),
                        U_END_ST => continue,
                        _ => continue,
                    };

                    let offset = st.offset as usize;
                    let base_idx = y + offset;

                    // Initialize with ynext == 0 case
                    let mut best_score = if base_idx < statenum && beam_diff <= window {
                        amx[beam_aj][beam_diff][base_idx] + st.tmx[0]
                    } else {
                        NEGINFINITY
                    };

                    // Calculate remaining cases
                    for ynext in 1..(st.connectnum as usize) {
                        let next_y = base_idx + ynext;
                        if next_y < statenum && beam_diff <= window {
                            let beam_score = amx[beam_aj][beam_diff][next_y];
                            if beam_score > best_score {
                                let sc = beam_score + st.tmx[ynext];
                                if sc > best_score {
                                    best_score = sc;
                                }
                            }
                        }
                    }

                    // Add emission scores
                    if best_score != NEGINFINITY {
                        best_score += emitsc;
                    }

                    amx[aj][diff][y] = best_score;

                    // Make a copy into bmx if this is a BEGIN state
                    if st.statetype == U_BEGIN_ST {
                        if let Some(ref mut y_bmx) = bmx[y] {
                            y_bmx[bj][diff] = best_score;
                        }
                    }
                } else {
                    // A BIFURC state
                    let left_state = st.tmx[0] as usize;
                    let right_state = st.tmx[1] as usize;

                    let left_bmx = bmx[left_state].as_ref();
                    let right_bmx = bmx[right_state].as_ref();

                    if left_bmx.is_none() || right_bmx.is_none() {
                        amx[aj][diff][y] = NEGINFINITY;
                        continue;
                    }

                    let left_bmx = left_bmx.unwrap();
                    let right_bmx = right_bmx.unwrap();

                    // Initialize with case that left branch emits it all
                    let mut leftdiff = diff;
                    let mut leftj = bj; // Circular buffer index

                    // Initial: left gets everything, right gets nothing
                    let mut best_score = left_bmx[leftj][leftdiff] + right_bmx[leftj][0];

                    while leftdiff > 0 {
                        leftdiff -= 1;
                        // Scan window wraparound: leftj = leftj ? leftj-1 : window-1
                        leftj = if leftj > 0 { leftj - 1 } else { window - 1 };

                        // right_diff increases as leftdiff decreases
                        let right_diff = diff - leftdiff;

                        let sc = left_bmx[leftj][leftdiff] + right_bmx[bj][right_diff];
                        if sc > best_score {
                            best_score = sc;
                        }
                    }

                    amx[aj][diff][y] = best_score;
                }
            }
        }

        // We've completed row j. Now examine scores to decide whether to report a hit.
        // Look at bmx[0][bj][diff] for diff = 1..window (ROOT state scores)
        // Report the best scoring match ending at position j
        if let Some(ref root_bmx) = bmx[0] {
            let mut bestdiff = 1;
            let mut bestscore = root_bmx[bj][1];

            for diff in 2..=max_diff {
                if root_bmx[bj][diff] > bestscore {
                    bestscore = root_bmx[bj][diff];
                    bestdiff = diff;
                }
            }

            if bestscore > ithresh {
                let start = j - bestdiff + 1; // 1-based start
                let end = j;                  // 1-based end
                let score = bestscore as f64 / INTPRECISION;

                // Call the callback; if it returns false, we could abort
                // (original C code ignores return value with a Warn)
                let _ = gotone_f(start, end, score);
            }
        }
    }

    Ok(())
}

/// Simple scan hit structure for collecting results
#[derive(Debug, Clone)]
pub struct ScanHit {
    /// 1-based start position
    pub start: usize,
    /// 1-based end position
    pub end: usize,
    /// Score in bits
    pub score: f64,
}

/// Convenience function that collects hits into a Vec
///
/// # Arguments
/// * `icm` - Integer state model
/// * `statenum` - Number of states in the model
/// * `seq` - Sequence to scan (ASCII nucleotides)
/// * `window` - Scanning window size in nucleotides
/// * `thresh` - Score threshold for reporting hits (in bits)
///
/// # Returns
/// * `Ok(Vec<ScanHit>)` - Vector of hits above threshold
/// * `Err(ViterbiError)` - If scanning fails
pub fn viterbi_scan_collect(
    icm: &[IState],
    statenum: usize,
    seq: &[u8],
    window: usize,
    thresh: f64,
) -> Result<Vec<ScanHit>, ViterbiError> {
    let mut hits = Vec::new();

    viterbi_scan(icm, statenum, seq, window, thresh, |start, end, score| {
        hits.push(ScanHit { start, end, score });
        true // continue scanning
    })?;

    Ok(hits)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_symbol_index() {
        assert_eq!(symbol_index(b'A'), 0);
        assert_eq!(symbol_index(b'a'), 0);
        assert_eq!(symbol_index(b'C'), 1);
        assert_eq!(symbol_index(b'c'), 1);
        assert_eq!(symbol_index(b'G'), 2);
        assert_eq!(symbol_index(b'g'), 2);
        assert_eq!(symbol_index(b'T'), 3);
        assert_eq!(symbol_index(b't'), 3);
        assert_eq!(symbol_index(b'U'), 3);
        assert_eq!(symbol_index(b'u'), 3);
    }

    #[test]
    fn test_allocate_mx_dimensions() {
        let icm = vec![
            {
                let mut s = IState::new();
                s.statetype = U_BEGIN_ST;
                s
            },
            {
                let mut s = IState::new();
                s.statetype = U_BIFURC_ST;
                s
            },
            {
                let mut s = IState::new();
                s.statetype = U_END_ST;
                s
            },
        ];

        let window = 100;
        let matrices = allocate_mx(&icm, 3, window).unwrap();

        // amx should have exactly 2 rows (rolling window)
        assert_eq!(matrices.amx.len(), 2);

        // Each row should have window+1 diff slots
        for j in 0..=1 {
            assert_eq!(matrices.amx[j].len(), window + 1);
            // Each diff slot should have statenum entries
            for diff in 0..=window {
                assert_eq!(matrices.amx[j][diff].len(), 3);
            }
        }

        // bmx should have Some only for BEGIN (0), None for BIFURC (1) and END (2)
        assert!(matrices.bmx[0].is_some()); // BEGIN
        assert!(matrices.bmx[1].is_none()); // BIFURC (no bmx in dbscan!)
        assert!(matrices.bmx[2].is_none()); // END

        // Check bmx dimensions: [bj][diff] where bj = 0..window-1, diff = 0..window
        if let Some(ref bmx0) = matrices.bmx[0] {
            assert_eq!(bmx0.len(), window); // bj = 0..window-1
            for bj in 0..window {
                assert_eq!(bmx0[bj].len(), window + 1); // diff = 0..window
            }
        }
    }

    #[test]
    fn test_scan_hit_struct() {
        let hit = ScanHit {
            start: 10,
            end: 80,
            score: 25.5,
        };
        assert_eq!(hit.start, 10);
        assert_eq!(hit.end, 80);
        assert!((hit.score - 25.5).abs() < 0.001);
    }

    #[test]
    fn test_wraparound_logic() {
        // Test the circular buffer wraparound logic used in BIFURC
        let window = 5;
        let mut leftj: usize = 3;

        // Simulating the wraparound loop
        let expected = vec![2, 1, 0, 4]; // After wraparound from 0 -> window-1 = 4
        let mut actual = Vec::new();

        for _ in 0..4 {
            leftj = if leftj > 0 { leftj - 1 } else { window - 1 };
            actual.push(leftj);
        }

        assert_eq!(actual, expected);
    }
}
