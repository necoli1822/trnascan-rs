//! CM to Integer State Model Conversion
//!
//! This module implements RearrangeCM from model.c lines 147-380.
//! It converts a probability-based CM into an integer log-odds state model
//! optimized for the Viterbi alignment algorithm.

use crate::types::cm::CM;
use crate::types::constants::*;
use crate::types::state::IState;

/// Convert a CM to an integer state-based model.
///
/// This is the Rust implementation of RearrangeCM() from model.c lines 175-380.
///
/// The integer CM is an array of IState structures in state-oriented form
/// (as opposed to node-oriented). The transition tables are rearranged
/// to optimize the recursion in recurse_mx():
/// - INSL, INSR transitions come first
/// - Order: INSL, INSR, DEL, MATP, MATL, MATR
///
/// # Arguments
/// * `cm` - The covariance model (probability form)
/// * `rfreq` - Background frequencies for log-odds calculation (typically [0.25; 4])
///
/// # Returns
/// * `(Vec<IState>, usize)` - The integer state model and number of states
pub fn rearrange_cm(cm: &CM, rfreq: &[f64; ALPHASIZE]) -> (Vec<IState>, usize) {
    // Allocate max possible states (nodes * STATETYPES)
    let max_states = cm.nodes * STATETYPES;
    let mut icm: Vec<IState> = Vec::with_capacity(max_states);

    // Stack for deferring bifurcation right connection assignment
    let mut bifstack: Vec<usize> = Vec::new();

    let mut y: usize = 0; // State counter

    for k in 0..cm.nodes {
        // Figure out what we're connected to (tflags)
        let mut tflags: u32 = if cm.nd[k].nxt == -1 {
            U_END_ST
        } else {
            let next_type = cm.nd[cm.nd[k].nxt as usize].node_type;
            match next_type as usize {
                BIFURC_NODE => U_BIFURC_ST,
                MATP_NODE => U_DEL_ST | U_MATP_ST | U_MATR_ST | U_MATL_ST,
                MATL_NODE => U_DEL_ST | U_MATL_ST,
                MATR_NODE => U_DEL_ST | U_MATR_ST,
                BEGINL_NODE | BEGINR_NODE => U_BEGIN_ST,
                _ => panic!("No such node type {}", next_type),
            }
        };

        // Figure out what we're coming from (fflags) and offset
        let (fflags, mut offset): (u32, i32) = match cm.nd[k].node_type as usize {
            BIFURC_NODE => (U_BIFURC_ST, 1),
            MATP_NODE => {
                tflags |= U_INSL_ST | U_INSR_ST;
                (
                    U_DEL_ST | U_MATP_ST | U_MATL_ST | U_MATR_ST | U_INSL_ST | U_INSR_ST,
                    4,
                )
            }
            MATL_NODE => {
                tflags |= U_INSL_ST;
                (U_DEL_ST | U_MATL_ST | U_INSL_ST, 2)
            }
            MATR_NODE => {
                tflags |= U_INSR_ST;
                (U_DEL_ST | U_MATR_ST | U_INSR_ST, 2)
            }
            BEGINL_NODE => (U_BEGIN_ST, 1),
            BEGINR_NODE => {
                tflags |= U_INSL_ST;
                (U_BEGIN_ST | U_INSL_ST, 1)
            }
            ROOT_NODE => {
                tflags |= U_INSL_ST | U_INSR_ST;
                (U_BEGIN_ST | U_INSL_ST | U_INSR_ST, 1)
            }
            _ => panic!("No such node type {}", cm.nd[k].node_type),
        };

        // Create states based on fflags

        // DEL state (or the first non-emit state)
        if fflags & U_DEL_ST != 0 {
            let mut st = IState::new();
            fill_state(&mut st, k as i32, U_DEL_ST, offset);
            copy_state_transitions(&mut st, &cm.nd[k].tmx[DEL_ST], tflags);
            icm.push(st);
            offset -= 1;
            y += 1;
        } else if fflags & U_BIFURC_ST != 0 {
            let mut st = IState::new();
            fill_state(&mut st, k as i32, U_BIFURC_ST, offset);
            // tmx[0] = left child (next state), tmx[1] = right child (deferred)
            st.tmx[0] = (y + 1) as i32;
            bifstack.push(y);
            icm.push(st);
            y += 1;
        } else if fflags & U_BEGIN_ST != 0 {
            let mut st = IState::new();
            fill_state(&mut st, k as i32, U_BEGIN_ST, offset);
            copy_state_transitions(&mut st, &cm.nd[k].tmx[DEL_ST], tflags);

            // If we're a right BEGIN, pop parent bifurc and set its right child
            if cm.nd[k].node_type as usize == BEGINR_NODE {
                if let Some(bifidx) = bifstack.pop() {
                    icm[bifidx].tmx[1] = y as i32;
                }
            }

            icm.push(st);
            offset -= 1;
            y += 1;
        }

        // MATP state
        if fflags & U_MATP_ST != 0 {
            let mut st = IState::new();
            fill_state(&mut st, k as i32, U_MATP_ST, offset);
            copy_pairwise_emissions(&mut st, &cm.nd[k].mp_emit, rfreq);
            copy_state_transitions(&mut st, &cm.nd[k].tmx[MATP_ST], tflags);
            icm.push(st);
            offset -= 1;
            y += 1;
        }

        // MATL state
        if fflags & U_MATL_ST != 0 {
            let mut st = IState::new();
            fill_state(&mut st, k as i32, U_MATL_ST, offset);
            copy_singlet_emissions(&mut st, &cm.nd[k].ml_emit, rfreq);
            copy_state_transitions(&mut st, &cm.nd[k].tmx[MATL_ST], tflags);
            icm.push(st);
            offset -= 1;
            y += 1;
        }

        // MATR state
        if fflags & U_MATR_ST != 0 {
            let mut st = IState::new();
            fill_state(&mut st, k as i32, U_MATR_ST, offset);
            copy_singlet_emissions(&mut st, &cm.nd[k].mr_emit, rfreq);
            copy_state_transitions(&mut st, &cm.nd[k].tmx[MATR_ST], tflags);
            icm.push(st);
            // offset is not used after this point, but we maintain the pattern
            // from the C code for potential future extensions
            let _ = offset - 1;
            y += 1;
        }

        // INSL state
        if fflags & U_INSL_ST != 0 {
            let mut st = IState::new();
            fill_state(&mut st, k as i32, U_INSL_ST, 0);
            copy_singlet_emissions(&mut st, &cm.nd[k].il_emit, rfreq);
            copy_state_transitions(&mut st, &cm.nd[k].tmx[INSL_ST], tflags);
            icm.push(st);
            y += 1;
        }

        // INSR state
        if fflags & U_INSR_ST != 0 {
            let mut st = IState::new();
            fill_state(&mut st, k as i32, U_INSR_ST, 0);
            copy_singlet_emissions(&mut st, &cm.nd[k].ir_emit, rfreq);
            // Note asymmetry: INSR->INSL transitions are disallowed
            copy_state_transitions(&mut st, &cm.nd[k].tmx[INSR_ST], tflags & !U_INSL_ST);
            icm.push(st);
            y += 1;
        }

        // End states must be added explicitly
        if cm.nd[k].nxt == -1 {
            let mut st = IState::new();
            fill_state(&mut st, -1, U_END_ST, 0);
            icm.push(st);
            y += 1;
        }
    }

    (icm, y)
}

/// Fill basic state information (model.c fill_state lines 394-403)
fn fill_state(st: &mut IState, nodeidx: i32, statetype: u32, offset: i32) {
    st.nodeidx = nodeidx;
    st.statetype = statetype;
    st.offset = offset;
}

/// Copy singlet emissions to state, converting to integer log-odds
/// (model.c copy_singlet_emissions lines 407-419)
fn copy_singlet_emissions(st: &mut IState, emvec: &[f64; ALPHASIZE], rfreq: &[f64; ALPHASIZE]) {
    for x in 0..ALPHASIZE {
        st.emit[x] = ilog2(emvec[x] / rfreq[x]);
    }
}

/// Copy pairwise emissions to state, converting to integer log-odds
/// (model.c copy_pairwise_emissions lines 424-439)
fn copy_pairwise_emissions(
    st: &mut IState,
    em: &[[f64; ALPHASIZE]; ALPHASIZE],
    rfreq: &[f64; ALPHASIZE],
) {
    for x in 0..(ALPHASIZE * ALPHASIZE) {
        let row = x / ALPHASIZE;
        let col = x % ALPHASIZE;
        st.emit[x] = ilog2(em[row][col] / (rfreq[col] * rfreq[row]));
    }
}

/// Copy state transitions, rearranging order and converting to integer log-odds
/// (model.c copy_state_transitions lines 444-480)
///
/// The transition vector is rearranged for optimization:
/// INSL, INSR are placed first, then DEL/BIFURC/BEGIN/END, then MATP, MATL, MATR
fn copy_state_transitions(st: &mut IState, tvec: &[f64; STATETYPES], tflags: u32) {
    let mut stx: usize = 0;

    if tflags & U_INSL_ST != 0 {
        st.tmx[stx] = ilog2(tvec[INSL_ST]);
        stx += 1;
    }

    if tflags & U_INSR_ST != 0 {
        st.tmx[stx] = ilog2(tvec[INSR_ST]);
        stx += 1;
    }

    if tflags & U_DEL_ST != 0
        || tflags & U_BIFURC_ST != 0
        || tflags & U_BEGIN_ST != 0
        || tflags & U_END_ST != 0
    {
        st.tmx[stx] = ilog2(tvec[DEL_ST]);
        stx += 1;
    }

    if tflags & U_MATP_ST != 0 {
        st.tmx[stx] = ilog2(tvec[MATP_ST]);
        stx += 1;
    }

    if tflags & U_MATL_ST != 0 {
        st.tmx[stx] = ilog2(tvec[MATL_ST]);
        stx += 1;
    }

    if tflags & U_MATR_ST != 0 {
        st.tmx[stx] = ilog2(tvec[MATR_ST]);
        stx += 1;
    }

    st.connectnum = stx as i32;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fill_state() {
        let mut st = IState::new();
        fill_state(&mut st, 5, U_MATP_ST, 3);
        assert_eq!(st.nodeidx, 5);
        assert_eq!(st.statetype, U_MATP_ST);
        assert_eq!(st.offset, 3);
    }

    #[test]
    fn test_copy_singlet_emissions() {
        let mut st = IState::new();
        let emvec = [0.25, 0.25, 0.25, 0.25];
        let rfreq = [0.25, 0.25, 0.25, 0.25];

        copy_singlet_emissions(&mut st, &emvec, &rfreq);

        // log2(1.0) = 0
        for x in 0..ALPHASIZE {
            assert_eq!(st.emit[x], 0);
        }
    }

    #[test]
    fn test_copy_singlet_emissions_biased() {
        let mut st = IState::new();
        let emvec = [0.5, 0.25, 0.125, 0.125];
        let rfreq = [0.25, 0.25, 0.25, 0.25];

        copy_singlet_emissions(&mut st, &emvec, &rfreq);

        // log2(0.5/0.25) = log2(2) = 1 * 1000 = 1000
        assert_eq!(st.emit[0], 1000);
        // log2(0.25/0.25) = log2(1) = 0
        assert_eq!(st.emit[1], 0);
        // log2(0.125/0.25) = log2(0.5) = -1 * 1000 = -1000
        assert_eq!(st.emit[2], -1000);
        assert_eq!(st.emit[3], -1000);
    }

    #[test]
    fn test_copy_state_transitions() {
        let mut st = IState::new();
        let tvec = [0.1, 0.2, 0.3, 0.2, 0.1, 0.1]; // DEL, MATP, MATL, MATR, INSL, INSR

        // Test with INSL and DEL flags
        copy_state_transitions(&mut st, &tvec, U_INSL_ST | U_DEL_ST);
        assert_eq!(st.connectnum, 2);

        // Test full flags
        let mut st2 = IState::new();
        copy_state_transitions(
            &mut st2,
            &tvec,
            U_INSL_ST | U_INSR_ST | U_DEL_ST | U_MATP_ST | U_MATL_ST | U_MATR_ST,
        );
        assert_eq!(st2.connectnum, 6);
    }

    #[test]
    fn test_rearrange_simple_cm() {
        // Create a minimal CM with just a root node
        let mut cm = CM::new(2);

        // Root node
        cm.nd[0].node_type = ROOT_NODE as i32;
        cm.nd[0].nxt = 1;
        cm.nd[0].nxt2 = -1;
        for i in 0..ALPHASIZE {
            cm.nd[0].il_emit[i] = 0.25;
            cm.nd[0].ir_emit[i] = 0.25;
        }
        for i in 0..STATETYPES {
            cm.nd[0].tmx[DEL_ST][i] = if i == DEL_ST { 1.0 } else { 0.0 };
            cm.nd[0].tmx[INSL_ST][i] = if i == DEL_ST { 1.0 } else { 0.0 };
            cm.nd[0].tmx[INSR_ST][i] = if i == DEL_ST { 1.0 } else { 0.0 };
        }

        // MATL node that ends
        cm.nd[1].node_type = MATL_NODE as i32;
        cm.nd[1].nxt = -1;
        cm.nd[1].nxt2 = -1;
        for i in 0..ALPHASIZE {
            cm.nd[1].ml_emit[i] = 0.25;
            cm.nd[1].il_emit[i] = 0.25;
        }
        for i in 0..STATETYPES {
            cm.nd[1].tmx[DEL_ST][i] = if i == DEL_ST { 1.0 } else { 0.0 };
            cm.nd[1].tmx[MATL_ST][i] = if i == DEL_ST { 1.0 } else { 0.0 };
            cm.nd[1].tmx[INSL_ST][i] = if i == DEL_ST { 1.0 } else { 0.0 };
        }

        let rfreq = [0.25, 0.25, 0.25, 0.25];
        let (icm, statenum) = rearrange_cm(&cm, &rfreq);

        // Should have: BEGIN, INSL, INSR (from root) + DEL, MATL, INSL, END (from MATL node)
        assert!(statenum > 0);
        assert_eq!(icm.len(), statenum);

        // First state should be BEGIN
        assert_eq!(icm[0].statetype, U_BEGIN_ST);
    }
}
