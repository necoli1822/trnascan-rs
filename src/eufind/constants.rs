// EuFindtRNA constants and weight matrices
// Based on original eufind_const.h and pavesi.c

// tRNA scanning cutoffs
pub const BBOX_CUTOFF: f32 = -14.14;
pub const BBOX_START_IDX: usize = 45;

pub const SEC_LOBOUND: f32 = -4.9;
pub const SEC_HIBOUND: f32 = -2.1;
pub const MAX_PENALTY: f32 = -5.442; // log(1/231)

pub const INT_SCORE_THRESH: f32 = -31.25;
pub const TOT_SCORE_THRESH: f32 = -31.8;

pub const MIN_AB_BOX_DIST: usize = 24;
pub const AB_BOX_DIST_RANGE: usize = 116;
pub const SEC_AB_BOX_DIST: usize = 26;
pub const SEC_BBOX_DIST_CORR: usize = 12;

pub const MIN_BTERM_DIST: usize = 11;
pub const MAX_TERM_SEARCH: usize = 133;

pub const ABOX_LEN: usize = 21;
pub const BBOX_LEN: usize = 11;

pub const MAX_OVLAP: usize = 10;

// Row indices for weight matrices
pub const GAP_ROW: usize = 4;
pub const AMBIG_ROW: usize = 5;

// A-box weight matrix [6 rows x 21 columns]
// Rows: 0=A, 1=C, 2=G, 3=T, 4=Gap, 5=Ambiguous
pub const ABOX_MAT: [[f32; ABOX_LEN]; 6] = [
    // A row
    [
        -1.268, -3.651, -0.899, -4.749, -5.442, -2.351, -3.363, -0.009, -1.977, -3.497, -5.442,
        -5.442, -5.442, -2.498, -4.749, -5.442, -0.031, -1.417, -1.180, -1.048, -4.344,
    ],
    // C row
    [
        -3.651, -5.442, -4.056, -2.958, -0.480, -1.073, -0.857, -5.442, -5.442, -1.887, -2.498,
        -5.442, -5.442, -2.958, -2.224, -5.442, -5.442, -3.363, -1.417, -3.651, -0.393,
    ],
    // G row
    [
        -0.779, -5.442, -0.598, -0.076, -3.651, -1.435, -1.614, -4.749, -0.154, -2.803, -5.442,
        0.000, 0.000, -3.363, -3.651, -5.442, -3.497, -0.672, -1.012, -0.473, -3.651,
    ],
    // T row
    [
        -1.453, -0.026, -3.651, -4.344, -1.036, -1.125, -1.073, -5.442, -5.442, -0.278, -1.399,
        -5.442, -5.442, -0.185, -0.827, -2.041, -5.442, -1.551, -2.447, -5.442, -1.253,
    ],
    // Gap row
    [
        -5.442, -5.442, -5.442, -5.442, -5.442, -5.442, -5.442, -5.442, -5.442, -5.442, -0.412,
        -5.442, -5.442, -5.442, -0.868, -0.144, -5.442, -5.442, -5.442, -5.442, -5.442,
    ],
    // Ambiguous row (min of ACGT for each position)
    [
        -0.779, -0.026, -0.598, -0.076, -0.480, -1.073, -0.857, -0.009, -0.154, -0.278, -1.399,
        0.000, 0.000, -0.185, -0.827, -2.041, -0.031, -0.672, -1.012, -0.473, -0.393,
    ],
];

// B-box weight matrix [6 rows x 11 columns]
// Rows: 0=A, 1=C, 2=G, 3=T, 4=Gap, 5=Ambiguous
pub const BBOX_MAT: [[f32; BBOX_LEN]; 6] = [
    // A row
    [
        -2.351, -5.442, -2.670, -5.442, -5.442, -1.472, 0.000, -0.798, -2.498, -5.442, -3.497,
    ],
    // C row
    [
        -3.245, -5.442, -5.442, -5.442, -0.004, -5.442, -5.442, -2.498, -1.435, -0.009, -0.190,
    ],
    // G row
    [
        -0.175, -0.004, -5.442, -5.442, -5.442, -0.272, -5.442, -2.147, -5.442, -5.442, -3.651,
    ],
    // T row
    [
        -3.651, -5.442, -0.072, 0.000, -5.442, -4.749, -5.442, -1.048, -0.393, -5.442, -2.147,
    ],
    // Gap row
    [
        -5.442, -5.442, -5.442, -5.442, -5.442, -5.442, -5.442, -5.442, -5.442, -5.442, -5.442,
    ],
    // Ambiguous row (min of ACGT for each position)
    [
        -0.175, -0.004, -0.072, 0.000, -0.004, -0.272, 0.000, -0.798, -0.393, -0.009, -0.190,
    ],
];

// A-B box distance weight matrices
pub const ABDIST_MAT_SIZE: usize = 7;
pub const AB_DIST_IDX_MAT: [i32; ABDIST_MAT_SIZE] = [30, 36, 42, 48, 54, 60, 66];
pub const AB_DIST_SC_MAT: [f32; ABDIST_MAT_SIZE] = [-0.46, -1.83, -2.35, -3.24, -4.06, -3.83, -4.75];

// B-box to termination distance weight matrices
pub const BTERM_MAT_SIZE: usize = 9;
pub const BTERM_DIST_IDX_MAT: [i32; BTERM_MAT_SIZE] = [17, 23, 29, 35, 41, 47, 53, 59, 100];
pub const BTERM_DIST_SC_MAT: [f32; BTERM_MAT_SIZE] =
    [-0.54, -1.40, -2.80, -3.36, -3.24, -5.44, -5.44, -4.06, -5.44];
