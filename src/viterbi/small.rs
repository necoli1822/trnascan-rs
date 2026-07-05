//! Small-memory Viterbi Algorithm
//!
//! This module implements SmallViterbiAlign from smallviterbi.c
//! It is a memory-efficient variant of the Viterbi algorithm using a 2-row rolling window.
//!
//! # Memory Optimization
//!
//! In the two-matrix version of the alignment algorithm, we keep:
//! - A matrix for BEGIN states (B matrix)
//! - A full matrix for all scores (A matrix) with only 2 rows (rolling window)
//!
//! The key insight is that for scoring, we only need the current row and the last row
//! in matrix A. We only need full information for BEGIN state scores in matrix B.
//!
//! # Traceback Information
//!
//! Matrix cells carry traceback information (tback) indicating the i,j coords that
//! a segment's BIFURC/END aligns to. The tback is a packed unsigned integer containing
//! two 16-bit values (i in high bits, j in low bits).

use crate::types::constants::*;
use crate::types::state::IState;
use crate::types::trace::Trace;

use super::ViterbiError;

/// Packed traceback type - stores i,j coordinates in a single u32
/// i in high 16 bits, j in low 16 bits
type TBack = u32;

/// Pack two indices into a single traceback value
/// i goes in high 16 bits, j goes in low 16 bits
#[inline]
fn pack_tb(i: usize, j: usize) -> TBack {
    ((i as u32) << 16) | (j as u32)
}

/// Unpack traceback value into (i, j) indices
/// Returns (i, j) where i was in high 16 bits, j was in low 16 bits
#[inline]
fn unpack_tb(tback: TBack) -> (usize, usize) {
    let j = (tback & 0xFFFF) as usize;
    let i = (tback >> 16) as usize;
    (i, j)
}

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

/// Small-memory Viterbi alignment matrices
struct SmallViterbiMatrices {
    /// Main score matrix: amx[aj][diff][y] where aj is 0 or 1 (rolling window)
    /// Dimensions: [2][seqlen+1][statenum]
    amx: Vec<Vec<Vec<i32>>>,

    /// Traceback pointers for amx: same dimensions as amx
    atr: Vec<Vec<Vec<TBack>>>,

    /// BEGIN/BIFURC score matrix: bmx[y] is Some([j][diff]) for BEGIN/BIFURC states
    /// Only BEGIN and BIFURC states have allocated matrices
    bmx: Vec<Option<Vec<Vec<i32>>>>,

    /// Traceback pointers for bmx: only for BEGIN states
    btr: Vec<Option<Vec<Vec<TBack>>>>,
}

/// Small-memory Viterbi alignment function
///
/// Implements SmallViterbiAlign() from smallviterbi.c lines 127-172
///
/// # Arguments
/// * `icm` - Integer state model
/// * `statenum` - Number of states in the model
/// * `seq` - Sequence to align (ASCII nucleotides)
///
/// # Returns
/// * `Ok((score, trace))` - Alignment score in bits and traceback tree
/// * `Err(ViterbiError)` - If alignment fails
pub fn small_viterbi_align(
    icm: &[IState],
    statenum: usize,
    seq: &[u8],
) -> Result<(f64, Trace), ViterbiError> {
    let n = seq.len();

    // Allocate matrices
    let mut matrices = allocate_mx(icm, statenum, n)?;

    // Initialize matrices
    init_mx(icm, statenum, n, &mut matrices)?;

    // Fill matrices via recursion
    recurse_mx(icm, statenum, seq, n, &mut matrices)?;

    // Get final score from bmx[0][N][N]
    // State 0 is always a BEGIN state (ROOT), so bmx[0] is always Some
    let score = matrices.bmx[0].as_ref().unwrap()[n][n] as f64 / INTPRECISION;

    // Perform traceback
    let trace = trace_mx(icm, seq, n, &matrices)?;

    Ok((score, trace))
}

/// Allocate the scoring matrices
///
/// Implements allocate_mx() from smallviterbi.c lines 198-281
///
/// Key differences from full Viterbi:
/// - amx/atr only have 2 rows (j = 0..1) for rolling window
/// - bmx allocated for BEGIN and BIFURC states
/// - btr allocated only for BEGIN states
fn allocate_mx(
    icm: &[IState],
    statenum: usize,
    seqlen: usize,
) -> Result<SmallViterbiMatrices, ViterbiError> {
    // Main matrix amx: only 2 rows for rolling window
    // amx[aj][diff][y] where aj is 0 or 1
    let mut amx: Vec<Vec<Vec<i32>>> = Vec::with_capacity(2);
    let mut atr: Vec<Vec<Vec<TBack>>> = Vec::with_capacity(2);

    for _j in 0..=1 {
        // diff ranges from 0 to seqlen
        let mut amx_row: Vec<Vec<i32>> = Vec::with_capacity(seqlen + 1);
        let mut atr_row: Vec<Vec<TBack>> = Vec::with_capacity(seqlen + 1);

        for _diff in 0..=seqlen {
            amx_row.push(vec![NEGINFINITY; statenum]);
            atr_row.push(vec![0; statenum]);
        }
        amx.push(amx_row);
        atr.push(atr_row);
    }

    // B auxiliary matrices: bmx[y][j][diff]
    // bmx keeps scores for BEGIN and BIFURC states
    // btr keeps traceback info only for BEGIN states
    let mut bmx: Vec<Option<Vec<Vec<i32>>>> = vec![None; statenum];
    let mut btr: Vec<Option<Vec<Vec<TBack>>>> = vec![None; statenum];

    for y in 0..statenum {
        // We keep score info for BEGIN and BIFURC states
        if icm[y].statetype == U_BEGIN_ST || icm[y].statetype == U_BIFURC_ST {
            let mut y_matrix: Vec<Vec<i32>> = Vec::with_capacity(seqlen + 1);
            for j in 0..=seqlen {
                y_matrix.push(vec![NEGINFINITY; j + 1]);
            }
            bmx[y] = Some(y_matrix);
        }

        // We keep traceback info only for BEGIN states
        if icm[y].statetype == U_BEGIN_ST {
            let mut y_matrix: Vec<Vec<TBack>> = Vec::with_capacity(seqlen + 1);
            for j in 0..=seqlen {
                y_matrix.push(vec![0; j + 1]);
            }
            btr[y] = Some(y_matrix);
        }
    }

    Ok(SmallViterbiMatrices { amx, atr, bmx, btr })
}

/// Initialize the scoring matrices
///
/// Implements init_mx() from smallviterbi.c lines 349-439
fn init_mx(
    icm: &[IState],
    statenum: usize,
    n: usize,
    matrices: &mut SmallViterbiMatrices,
) -> Result<(), ViterbiError> {
    let SmallViterbiMatrices {
        amx,
        atr,
        bmx,
        btr,
    } = matrices;

    // Init the whole amx to -Infinity (done in allocation)
    // Re-initialize to be safe
    for j in 0..=1 {
        for diff in 0..=n {
            for y in 0..statenum {
                amx[j][diff][y] = NEGINFINITY;
            }
        }
    }

    // atr END and BIFURC traceback pointers point to themselves
    // Just set everything to point at itself: pack_tb(diff, j)
    for j in 0..=1 {
        for diff in 0..=n {
            for y in 0..statenum {
                atr[j][diff][y] = pack_tb(diff, j);
            }
        }
    }

    // Init the whole bmx to -Infinity
    // State 0 is always a BEGIN (ROOT), start there
    for y in 0..statenum {
        if let Some(ref mut y_bmx) = bmx[y] {
            for j in 0..=n {
                for diff in 0..=j {
                    y_bmx[j][diff] = NEGINFINITY;
                }
            }
        }
    }

    // Set all btr traceback ptrs to point at themselves
    for y in 0..statenum {
        if let Some(ref mut y_btr) = btr[y] {
            for j in 0..=n {
                for diff in 0..=j {
                    y_btr[j][diff] = pack_tb(diff, j);
                }
            }
        }
    }

    // Init the off-diagonal (j = 0..N; diff == 0) with -log P scores
    // End state = 0; del, bifurc states are calc'd; begin states same as del's
    for j in 0..=n {
        let aj = j % 2;

        // Process states in reverse order
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

            // Make a copy into bmx if y is a BEGIN or BIFURC
            if icm[y].statetype == U_BEGIN_ST || icm[y].statetype == U_BIFURC_ST {
                if let Some(ref mut y_bmx) = bmx[y] {
                    y_bmx[j][0] = amx[aj][0][y];
                }
            }
        }
    }

    Ok(())
}

/// Fill the scoring matrices via recursion
///
/// Implements recurse_mx() from smallviterbi.c lines 450-607
///
/// CRITICAL: Uses rolling window with aj = j % 2 and !aj for previous row
fn recurse_mx(
    icm: &[IState],
    statenum: usize,
    seq: &[u8],
    n: usize,
    matrices: &mut SmallViterbiMatrices,
) -> Result<(), ViterbiError> {
    let SmallViterbiMatrices {
        amx,
        atr,
        bmx,
        btr,
    } = matrices;

    for j in 1..=n {
        let aj = j % 2;
        let prev_aj = 1 - aj; // !aj toggles 0<->1
        let symj = symbol_index(seq[j - 1]); // 1-based to 0-based

        // Initialize END and BIFURC states to point at themselves in this row
        for diff in 0..=j {
            for y in 0..statenum {
                if icm[y].statetype == U_BIFURC_ST || icm[y].statetype == U_END_ST {
                    atr[aj][diff][y] = pack_tb(diff, j);
                }
            }
        }

        for diff in 1..=j {
            let i = j - diff + 1; // 1-based i
            let symi = symbol_index(seq[i - 1]); // 1-based to 0-based

            // Process states in reverse order
            for y in (0..statenum).rev() {
                let st = &icm[y];

                if st.statetype != U_BIFURC_ST {
                    // A normal (non-BIFURC) state
                    // Connect beam pointer to appropriate starting place

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
                    let mut best_score = if base_idx < statenum && beam_diff <= n {
                        amx[beam_aj][beam_diff][base_idx] + st.tmx[0]
                    } else {
                        NEGINFINITY
                    };

                    let mut best_tback = if base_idx < statenum && beam_diff <= n {
                        atr[beam_aj][beam_diff][base_idx]
                    } else {
                        pack_tb(diff, j)
                    };

                    // Calculate remaining cases
                    for ynext in 1..(st.connectnum as usize) {
                        let next_y = base_idx + ynext;
                        if next_y < statenum && beam_diff <= n {
                            let beam_score = amx[beam_aj][beam_diff][next_y];
                            if beam_score > best_score {
                                let sc = beam_score + st.tmx[ynext];
                                if sc > best_score {
                                    best_score = sc;
                                    best_tback = atr[beam_aj][beam_diff][next_y];
                                }
                            }
                        }
                    }

                    // Add emission scores
                    if best_score != NEGINFINITY {
                        best_score += emitsc;
                    }

                    amx[aj][diff][y] = best_score;
                    atr[aj][diff][y] = best_tback;

                    // Make a copy into bmx, btr if this is a BEGIN state
                    if st.statetype == U_BEGIN_ST {
                        if let Some(ref mut y_bmx) = bmx[y] {
                            y_bmx[j][diff] = best_score;
                        }
                        if let Some(ref mut y_btr) = btr[y] {
                            y_btr[j][diff] = best_tback;
                        }
                    }
                } else {
                    // A BIFURC state
                    let left = st.tmx[0] as usize;
                    let right = st.tmx[1] as usize;

                    let left_bmx = bmx[left].as_ref();
                    let right_bmx = bmx[right].as_ref();

                    if left_bmx.is_none() || right_bmx.is_none() {
                        amx[aj][diff][y] = NEGINFINITY;
                        continue;
                    }

                    let left_bmx = left_bmx.unwrap();
                    let right_bmx = right_bmx.unwrap();

                    // Initialize with case that left branch emits it all
                    let mut leftdiff = diff;
                    let mut leftj = j;

                    let mut best_score = left_bmx[leftj][leftdiff] + right_bmx[j][0];

                    while leftdiff > 0 {
                        leftdiff -= 1;
                        leftj -= 1;

                        let right_diff = diff - leftdiff;
                        if right_diff <= j {
                            let sc = left_bmx[leftj][leftdiff] + right_bmx[j][right_diff];
                            if sc > best_score {
                                best_score = sc;
                            }
                        }
                    }

                    amx[aj][diff][y] = best_score;

                    // Keep copy of score in bmx for tracing
                    if let Some(ref mut y_bmx) = bmx[y] {
                        y_bmx[j][diff] = best_score;
                    }
                }
            }
        }
    }

    Ok(())
}

/// Trace back through the matrices to reconstruct the alignment
///
/// Implements trace_mx() from smallviterbi.c lines 617-720
///
/// Uses bmx as framework and recalculates bits of alignment between
/// BEGINs and BIFURC/ENDs.
fn trace_mx(
    icm: &[IState],
    seq: &[u8],
    n: usize,
    matrices: &SmallViterbiMatrices,
) -> Result<Trace, ViterbiError> {
    let SmallViterbiMatrices { bmx, btr, .. } = matrices;

    // Start at i = 1, j = N (0-indexed: i = 0, j = N-1)
    // State 0 is always the ROOT BEGIN state
    let root = trace_segment(icm, seq, n, bmx, btr, 0, n.saturating_sub(1), 0)?;

    Ok(root)
}

/// Trace a single segment from BEGIN to BIFURC/END
fn trace_segment(
    icm: &[IState],
    seq: &[u8],
    n: usize,
    bmx: &[Option<Vec<Vec<i32>>>],
    btr: &[Option<Vec<Vec<TBack>>>],
    i: usize,     // 0-indexed left position
    j: usize,     // 0-indexed right position
    y: usize,     // state index (BEGIN state)
) -> Result<Trace, ViterbiError> {
    // Convert to 1-based for internal calculations
    let i_1based = i + 1;
    let j_1based = j + 1;
    let diff = if j_1based >= i_1based {
        j_1based - i_1based + 1
    } else {
        0
    };

    let st = &icm[y];

    // Create trace node
    let mut trace = Trace::new(i as i32, j as i32, st.nodeidx, st.statetype);

    // Find which BIFURC/END terminates this segment
    let mut end_y = y + 1;
    while end_y < icm.len()
        && icm[end_y].statetype != U_BIFURC_ST
        && icm[end_y].statetype != U_END_ST
    {
        end_y += 1;
    }

    if end_y >= icm.len() {
        return Err(ViterbiError::TracebackFailed(format!(
            "Could not find END/BIFURC after state {}",
            y
        )));
    }

    // Get i2, j2 from the traceback pointer
    let (diff2, j2) = if let Some(ref y_btr) = btr[y] {
        if j_1based <= n && diff <= j_1based {
            unpack_tb(y_btr[j_1based][diff])
        } else {
            (diff, j_1based)
        }
    } else {
        (diff, j_1based)
    };

    let i2 = if j2 >= diff2 { j2 - diff2 + 1 } else { j2 + 1 };

    // If this segment ends at a BIFURC, push children
    if icm[end_y].statetype == U_BIFURC_ST {
        let left_state = icm[end_y].tmx[0] as usize;
        let right_state = icm[end_y].tmx[1] as usize;

        if i2 > j2 {
            // Off diagonal - both branches get same empty range
            let left_child = trace_segment(
                icm,
                seq,
                n,
                bmx,
                btr,
                i2.saturating_sub(1),
                j2.saturating_sub(1),
                left_state,
            )?;
            let right_child = trace_segment(
                icm,
                seq,
                n,
                bmx,
                btr,
                i2.saturating_sub(1),
                j2.saturating_sub(1),
                right_state,
            )?;
            trace.nxtl = Some(Box::new(left_child));
            trace.nxtr = Some(Box::new(right_child));
        } else {
            // Find the split point
            let left_bmx = bmx[left_state].as_ref();
            let right_bmx = bmx[right_state].as_ref();

            if left_bmx.is_none() || right_bmx.is_none() {
                return Err(ViterbiError::TracebackFailed(
                    "Missing bmx for bifurcation".to_string(),
                ));
            }

            let left_bmx = left_bmx.unwrap();
            let right_bmx = right_bmx.unwrap();

            let bifurc_bmx = bmx[end_y].as_ref();
            if bifurc_bmx.is_none() {
                return Err(ViterbiError::TracebackFailed(
                    "Missing bmx for bifurc state".to_string(),
                ));
            }
            let bifurc_bmx = bifurc_bmx.unwrap();

            let mut found = false;
            let mut leftdiff = diff2;
            let mut leftj = j2;

            while leftdiff > 0 && !found {
                let right_diff = diff2 - leftdiff;
                if leftj <= n && leftdiff <= leftj && j2 <= n && right_diff <= j2 {
                    if bifurc_bmx[j2][diff2]
                        == left_bmx[leftj][leftdiff] + right_bmx[j2][right_diff]
                    {
                        // Found the split point
                        let left_i = i2.saturating_sub(1);
                        let left_j = if leftdiff > 0 {
                            i2 + leftdiff - 2
                        } else {
                            left_i
                        };
                        let right_i = i2 + leftdiff - 1;
                        let right_j = j2.saturating_sub(1);

                        let left_child =
                            trace_segment(icm, seq, n, bmx, btr, left_i, left_j, left_state)?;
                        let right_child =
                            trace_segment(icm, seq, n, bmx, btr, right_i, right_j, right_state)?;

                        trace.nxtl = Some(Box::new(left_child));
                        trace.nxtr = Some(Box::new(right_child));
                        found = true;
                    }
                }

                leftdiff -= 1;
                if leftj > 0 {
                    leftj -= 1;
                }
            }

            // Check if left branch is empty (leftdiff == 0)
            if !found {
                // Left gets empty, right gets everything
                let left_i = i2.saturating_sub(1);
                let left_j = left_i.saturating_sub(1);
                let right_i = i2.saturating_sub(1);
                let right_j = j2.saturating_sub(1);

                let left_child =
                    trace_segment(icm, seq, n, bmx, btr, left_i, left_j, left_state)?;
                let right_child =
                    trace_segment(icm, seq, n, bmx, btr, right_i, right_j, right_state)?;

                trace.nxtl = Some(Box::new(left_child));
                trace.nxtr = Some(Box::new(right_child));
            }
        }
    }
    // END state - no children needed

    Ok(trace)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pack_unpack_tb() {
        // Test basic packing
        assert_eq!(pack_tb(0, 0), 0);
        assert_eq!(pack_tb(1, 0), 0x10000);
        assert_eq!(pack_tb(0, 1), 1);
        assert_eq!(pack_tb(1, 1), 0x10001);

        // Test round-trip
        for i in 0..100 {
            for j in 0..100 {
                let packed = pack_tb(i, j);
                let (unpacked_i, unpacked_j) = unpack_tb(packed);
                assert_eq!(unpacked_i, i);
                assert_eq!(unpacked_j, j);
            }
        }

        // Test max values (16-bit)
        let packed = pack_tb(0xFFFF, 0xFFFF);
        let (i, j) = unpack_tb(packed);
        assert_eq!(i, 0xFFFF);
        assert_eq!(j, 0xFFFF);
    }

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

        let matrices = allocate_mx(&icm, 3, 5).unwrap();

        // amx should have exactly 2 rows (rolling window)
        assert_eq!(matrices.amx.len(), 2);
        assert_eq!(matrices.atr.len(), 2);

        // Each row should have seqlen+1 = 6 diff slots
        for j in 0..=1 {
            assert_eq!(matrices.amx[j].len(), 6);
            assert_eq!(matrices.atr[j].len(), 6);
        }

        // bmx should have Some for BEGIN (0) and BIFURC (1), None for END (2)
        assert!(matrices.bmx[0].is_some()); // BEGIN
        assert!(matrices.bmx[1].is_some()); // BIFURC
        assert!(matrices.bmx[2].is_none()); // END

        // btr should have Some only for BEGIN
        assert!(matrices.btr[0].is_some()); // BEGIN
        assert!(matrices.btr[1].is_none()); // BIFURC (no traceback)
        assert!(matrices.btr[2].is_none()); // END

        // Check bmx dimensions: [j][diff] where diff <= j
        if let Some(ref bmx0) = matrices.bmx[0] {
            assert_eq!(bmx0.len(), 6); // j = 0..5
            for j in 0..=5 {
                assert_eq!(bmx0[j].len(), j + 1); // diff = 0..j
            }
        }
    }
}
