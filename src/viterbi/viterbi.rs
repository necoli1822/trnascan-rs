//! Viterbi Dynamic Programming Algorithm
//!
//! This module implements ViterbiAlign from viterbi.c lines 56-122.
//! It performs optimal alignment of a sequence to a covariance model
//! using a 3D dynamic programming algorithm.
//!
//! # Matrix Layout
//!
//! The algorithm uses two matrices:
//! - `amx[j][diff][y]`: Main score matrix
//! - `bmx[y][j][diff]`: BEGIN state matrix for bifurcation handling
//!
//! Where:
//! - j = right sequence position (0..N)
//! - diff = j - i + 1 (difference, 0 = off-diagonal, 1 = diagonal)
//! - y = state index (0..statenum-1)

use crate::types::constants::*;
use crate::types::state::IState;
use crate::types::trace::Trace;

/// Error type for Viterbi operations
#[derive(Debug, Clone)]
pub enum ViterbiError {
    AllocationFailed(String),
    InitializationFailed(String),
    RecursionFailed(String),
    TracebackFailed(String),
}

impl std::fmt::Display for ViterbiError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ViterbiError::AllocationFailed(msg) => write!(f, "Allocation failed: {}", msg),
            ViterbiError::InitializationFailed(msg) => write!(f, "Initialization failed: {}", msg),
            ViterbiError::RecursionFailed(msg) => write!(f, "Recursion failed: {}", msg),
            ViterbiError::TracebackFailed(msg) => write!(f, "Traceback failed: {}", msg),
        }
    }
}

impl std::error::Error for ViterbiError {}

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

/// Main Viterbi alignment function
///
/// Implements ViterbiAlign() from viterbi.c lines 70-122
///
/// # Arguments
/// * `icm` - Integer state model from rearrange_cm()
/// * `statenum` - Number of states in the model
/// * `seq` - Sequence to align (ASCII nucleotides)
///
/// # Returns
/// * `Ok((score, trace))` - Alignment score in bits and traceback tree
/// * `Err(ViterbiError)` - If alignment fails
pub fn viterbi_align(
    icm: &[IState],
    statenum: usize,
    seq: &[u8],
) -> Result<(f64, Trace), ViterbiError> {
    let n = seq.len();

    // Allocate matrices
    let (mut amx, mut bmx) = allocate_mx(icm, statenum, n)?;

    // Initialize matrices
    init_mx(icm, statenum, n, &mut amx, &mut bmx)?;

    // Fill matrices via recursion
    recurse_mx(icm, statenum, seq, n, &mut amx, &mut bmx)?;

    // Get final score from bmx[0][N][N]
    // State 0 is always a BEGIN state (ROOT), so bmx[0] is always Some
    let score = bmx[0].as_ref().unwrap()[n][n] as f64 / INTPRECISION;

    // Perform traceback
    let trace = trace_mx(icm, seq, n, &amx, &bmx)?;

    Ok((score, trace))
}

/// Allocate the scoring matrices
///
/// Implements allocate_mx() from viterbi.c lines 144-199
///
/// amx is [j = 0..N] [diff = 0..j] [y = 0..statenum-1]
/// bmx is [y = 0..statenum-1] [j = 0..N] [diff = 0..j] (only for BEGIN states)
fn allocate_mx(
    icm: &[IState],
    statenum: usize,
    seqlen: usize,
) -> Result<(Vec<Vec<Vec<i32>>>, Vec<Option<Vec<Vec<i32>>>>), ViterbiError> {
    // Main matrix amx[j][diff][y]
    let mut amx: Vec<Vec<Vec<i32>>> = Vec::with_capacity(seqlen + 1);

    for j in 0..=seqlen {
        let mut row: Vec<Vec<i32>> = Vec::with_capacity(j + 1);
        for _diff in 0..=j {
            row.push(vec![NEGINFINITY; statenum]);
        }
        amx.push(row);
    }

    // BEGIN auxiliary matrix bmx[y][j][diff]
    let mut bmx: Vec<Option<Vec<Vec<i32>>>> = vec![None; statenum];

    for y in 0..statenum {
        if icm[y].statetype == U_BEGIN_ST {
            let mut y_matrix: Vec<Vec<i32>> = Vec::with_capacity(seqlen + 1);
            for j in 0..=seqlen {
                y_matrix.push(vec![NEGINFINITY; j + 1]);
            }
            bmx[y] = Some(y_matrix);
        }
    }

    Ok((amx, bmx))
}

/// Initialize the scoring matrices
///
/// Implements init_mx() from viterbi.c lines 249-318
fn init_mx(
    icm: &[IState],
    statenum: usize,
    n: usize,
    amx: &mut [Vec<Vec<i32>>],
    bmx: &mut [Option<Vec<Vec<i32>>>],
) -> Result<(), ViterbiError> {
    // Init the whole amx to -Infinity (already done in allocation)

    // Init the whole bmx to -Infinity (already done in allocation)

    // Init the off-diagonal (j = 0..N; diff == 0) with -log P scores
    // End state = 0; del, bifurc states are calc'd; begin states same as del's
    for j in 0..=n {
        // Process states in reverse order (statenum-1 down to 0)
        for y in (0..statenum).rev() {
            if icm[y].statetype == U_END_ST {
                amx[j][0][y] = 0;
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
                    amx[j][0][y] = left_score + right_score;
                }
            } else if icm[y].statetype == U_DEL_ST || icm[y].statetype == U_BEGIN_ST {
                // Only calc DEL-DEL and BEGIN-DEL transitions
                // Find the connection to a non-infinite score
                let offset = icm[y].offset as usize;
                let base_idx = y + offset;

                for ynext in 0..icm[y].connectnum as usize {
                    let next_y = base_idx + ynext;
                    if next_y < statenum && amx[j][0][next_y] != NEGINFINITY {
                        amx[j][0][y] = amx[j][0][next_y] + icm[y].tmx[ynext];
                        break;
                    }
                }

                // Make a copy into bmx if y is a BEGIN
                if icm[y].statetype == U_BEGIN_ST {
                    if let Some(ref mut y_bmx) = bmx[y] {
                        y_bmx[j][0] = amx[j][0][y];
                    }
                }
            }
        }
    }

    Ok(())
}

/// Fill the scoring matrices via recursion
///
/// Implements recurse_mx() from viterbi.c lines 329-460
fn recurse_mx(
    icm: &[IState],
    statenum: usize,
    seq: &[u8],
    n: usize,
    amx: &mut [Vec<Vec<i32>>],
    bmx: &mut [Option<Vec<Vec<i32>>>],
) -> Result<(), ViterbiError> {
    for j in 1..=n {
        let symj = symbol_index(seq[j - 1]); // 0-indexed sequence

        for diff in 1..=j {
            let i = j - diff + 1;
            if i < 1 {
                break;
            }

            let symi = symbol_index(seq[i - 1]); // 0-indexed sequence

            // Process states in reverse order
            for y in (0..statenum).rev() {
                let st = &icm[y];

                if st.statetype != U_BIFURC_ST {
                    // Normal (non-BIFURC) state
                    let (beam_j, beam_diff, emitsc): (usize, usize, i32) = match st.statetype {
                        U_BEGIN_ST | U_DEL_ST => (j, diff, 0),
                        U_MATP_ST => {
                            if diff == 1 {
                                continue;
                            }
                            (j - 1, diff - 2, st.emit[symi * ALPHASIZE + symj])
                        }
                        U_MATR_ST | U_INSR_ST => (j - 1, diff - 1, st.emit[symj]),
                        U_MATL_ST | U_INSL_ST => (j, diff - 1, st.emit[symi]),
                        U_END_ST => continue,
                        _ => continue,
                    };

                    let offset = st.offset as usize;
                    let base_idx = y + offset;

                    // Initialize with first connection
                    let mut best_score = if base_idx < statenum && beam_diff <= beam_j {
                        amx[beam_j][beam_diff][base_idx] + st.tmx[0]
                    } else {
                        NEGINFINITY
                    };

                    // Check remaining connections
                    for ynext in 1..(st.connectnum as usize) {
                        let next_y = base_idx + ynext;
                        if next_y < statenum && beam_diff <= beam_j {
                            let beam_score = amx[beam_j][beam_diff][next_y];
                            if beam_score > best_score {
                                let sc = beam_score + st.tmx[ynext];
                                if sc > best_score {
                                    best_score = sc;
                                }
                            }
                        }
                    }

                    // Add emission score
                    if best_score != NEGINFINITY {
                        best_score += emitsc;
                    }

                    amx[j][diff][y] = best_score;

                    // Copy to bmx if BEGIN state
                    if st.statetype == U_BEGIN_ST {
                        if let Some(ref mut y_bmx) = bmx[y] {
                            y_bmx[j][diff] = best_score;
                        }
                    }
                } else {
                    // BIFURC state
                    let left = st.tmx[0] as usize;
                    let right = st.tmx[1] as usize;

                    let left_bmx = bmx[left].as_ref();
                    let right_bmx = bmx[right].as_ref();

                    if left_bmx.is_none() || right_bmx.is_none() {
                        amx[j][diff][y] = NEGINFINITY;
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

                    amx[j][diff][y] = best_score;
                }
            }
        }
    }

    Ok(())
}

/// Trace back through the matrices to reconstruct the alignment
///
/// Implements trace_mx() from viterbi.c lines 470-619
fn trace_mx(
    icm: &[IState],
    seq: &[u8],
    n: usize,
    amx: &[Vec<Vec<i32>>],
    bmx: &[Option<Vec<Vec<i32>>>],
) -> Result<Trace, ViterbiError> {
    // Start trace from BEGIN state 0 covering the whole sequence
    let root = trace_recursive(icm, seq, n, amx, bmx, 0, n - 1, 0)?;
    Ok(root)
}

/// Recursive traceback helper
fn trace_recursive(
    icm: &[IState],
    seq: &[u8],
    n: usize,
    amx: &[Vec<Vec<i32>>],
    bmx: &[Option<Vec<Vec<i32>>>],
    i: usize,     // 0-indexed left position
    j: usize,     // 0-indexed right position
    y: usize,     // state index
) -> Result<Trace, ViterbiError> {
    // Convert to 1-indexed for internal calculations (matching C code)
    let i_1based = i + 1;
    let j_1based = j + 1;
    let diff = if j_1based >= i_1based {
        j_1based - i_1based + 1
    } else {
        0
    };

    let st = &icm[y];

    // Create trace node for current state
    let mut trace = Trace::new(i as i32, j as i32, st.nodeidx, st.statetype);

    // END state - we're done
    if st.statetype == U_END_ST {
        return Ok(trace);
    }

    // BIFURC state - find the split point
    if st.statetype == U_BIFURC_ST {
        let left_state = st.tmx[0] as usize;
        let right_state = st.tmx[1] as usize;

        if i_1based > j_1based {
            // Off diagonal - both branches get same empty range
            let left_child = trace_recursive(icm, seq, n, amx, bmx, i, j, left_state)?;
            let right_child = trace_recursive(icm, seq, n, amx, bmx, i, j, right_state)?;
            trace.nxtl = Some(Box::new(left_child));
            trace.nxtr = Some(Box::new(right_child));
        } else {
            let left_bmx = bmx[left_state].as_ref().unwrap();
            let right_bmx = bmx[right_state].as_ref().unwrap();

            // Find split point
            let mut found = false;
            let mut leftdiff = diff;
            let mut leftj = j_1based;

            while leftdiff > 0 && !found {
                let right_diff = diff - leftdiff;
                if amx[j_1based][diff][y]
                    == left_bmx[leftj][leftdiff] + right_bmx[j_1based][right_diff]
                {
                    // Found the split
                    let left_i = i;
                    let left_j = i + leftdiff - 1;
                    let right_i = left_j + 1;
                    let right_j = j;

                    let left_child =
                        trace_recursive(icm, seq, n, amx, bmx, left_i, left_j, left_state)?;
                    let right_child =
                        trace_recursive(icm, seq, n, amx, bmx, right_i, right_j, right_state)?;

                    trace.nxtl = Some(Box::new(left_child));
                    trace.nxtr = Some(Box::new(right_child));
                    found = true;
                }

                leftdiff -= 1;
                leftj -= 1;
            }

            if !found && diff > 0 {
                // Check if left branch is empty
                if amx[j_1based][diff][y] == left_bmx[i_1based - 1][0] + right_bmx[j_1based][diff] {
                    // Left is empty
                    let left_child =
                        trace_recursive(icm, seq, n, amx, bmx, i, i.saturating_sub(1), left_state)?;
                    let right_child =
                        trace_recursive(icm, seq, n, amx, bmx, i, j, right_state)?;

                    trace.nxtl = Some(Box::new(left_child));
                    trace.nxtr = Some(Box::new(right_child));
                    found = true;
                }
            }

            if !found {
                return Err(ViterbiError::TracebackFailed(format!(
                    "Bifurc reconstruction failed at i={}, j={}, y={}",
                    i, j, y
                )));
            }
        }
    } else {
        // Normal state - find the next state
        let symi = if i < n { symbol_index(seq[i]) } else { 0 };
        let symj = if j < n { symbol_index(seq[j]) } else { 0 };

        let (next_i, next_j, beam_j, beam_diff): (usize, usize, usize, usize) = match st.statetype {
            U_BEGIN_ST | U_DEL_ST => (i, j, j_1based, diff),
            U_MATP_ST => (i + 1, j.saturating_sub(1), j_1based - 1, diff.saturating_sub(2)),
            U_MATR_ST | U_INSR_ST => (i, j.saturating_sub(1), j_1based - 1, diff.saturating_sub(1)),
            U_MATL_ST | U_INSL_ST => (i + 1, j, j_1based, diff.saturating_sub(1)),
            _ => return Err(ViterbiError::TracebackFailed("Invalid state type".to_string())),
        };

        // Calculate expected score (subtract emission)
        let mut sc = amx[j_1based][diff][y];
        match st.statetype {
            U_MATP_ST => sc -= st.emit[symi * ALPHASIZE + symj],
            U_MATR_ST | U_INSR_ST => sc -= st.emit[symj],
            U_MATL_ST | U_INSL_ST => sc -= st.emit[symi],
            _ => {}
        }

        // Find matching connection
        let offset = st.offset as usize;
        let base_idx = y + offset;
        let mut found_next = false;

        for ynext in 0..(st.connectnum as usize) {
            let next_y = base_idx + ynext;
            if next_y < icm.len() && beam_diff <= beam_j {
                if sc == amx[beam_j][beam_diff][next_y] + st.tmx[ynext] {
                    let child = trace_recursive(icm, seq, n, amx, bmx, next_i, next_j, next_y)?;
                    trace.nxtl = Some(Box::new(child));
                    found_next = true;
                    break;
                }
            }
        }

        if !found_next {
            return Err(ViterbiError::TracebackFailed(format!(
                "Can't continue traceback at i={}, j={}, y={}",
                i, j, y
            )));
        }
    }

    Ok(trace)
}

/// Prepare a sequence for alignment
///
/// Converts to uppercase and handles special characters
pub fn prepare_sequence(seq: &str) -> Vec<u8> {
    seq.bytes()
        .filter(|c| c.is_ascii_alphabetic())
        .map(|c| c.to_ascii_uppercase())
        .collect()
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
    fn test_prepare_sequence() {
        let seq = "ACGTacgt\nUu";
        let prepared = prepare_sequence(seq);
        assert_eq!(prepared, b"ACGTACGTUU");
    }

    #[test]
    fn test_allocate_mx() {
        let icm = vec![
            {
                let mut s = IState::new();
                s.statetype = U_BEGIN_ST;
                s
            },
            {
                let mut s = IState::new();
                s.statetype = U_END_ST;
                s
            },
        ];

        let (amx, bmx) = allocate_mx(&icm, 2, 5).unwrap();

        // Check amx dimensions
        assert_eq!(amx.len(), 6); // j = 0..5
        for j in 0..=5 {
            assert_eq!(amx[j].len(), j + 1); // diff = 0..j
            for diff in 0..=j {
                assert_eq!(amx[j][diff].len(), 2); // y = 0..1
            }
        }

        // Check bmx - only BEGIN states have matrices
        assert!(bmx[0].is_some()); // BEGIN
        assert!(bmx[1].is_none()); // END
    }
}
