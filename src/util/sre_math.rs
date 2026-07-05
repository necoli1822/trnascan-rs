//! sre_math.rs - Math utilities matching sre_math.c
//!
//! This module implements the math functions from the original SQUID library's
//! sre_math.c, maintaining exact compatibility with the C implementation for
//! reproducible random number generation.

use std::cell::Cell;
use std::f64::consts::PI;

// Random number generator constants (sre_math.c lines 352-354)
const RANGE: i64 = 268435456; // 2^28
const DIV: i64 = 16384;        // sqrt(RANGE)
const MULT: i64 = 72530821;    // LCG multiplier

// Thread-local state for random number generator
thread_local! {
    static SRE_RANDSEED: Cell<i64> = const { Cell::new(666) };
    static SRE_RESEED: Cell<bool> = const { Cell::new(false) };
    static FIRST_TIME: Cell<bool> = const { Cell::new(true) };
    static RND: Cell<i64> = const { Cell::new(0) };
}

/// Set random seed - sre_srandom()
///
/// Corresponds to sre_srandom() in sre_math.c (lines 387-393)
pub fn sre_srandom(mut seed: i64) {
    if seed < 0 {
        seed = -seed;
    }
    SRE_RESEED.with(|r| r.set(true));
    SRE_RANDSEED.with(|s| s.set(seed));
}

/// Linear congruential random number generator - sre_random()
///
/// Returns uniform random in [0.0, 1.0)
/// Based on sre_random() in sre_math.c (lines 355-376)
///
/// Uses a simple linear congruential generator with period 2^28.
/// Based on discussion in Robert Sedgewick's _Algorithms in C_.
///
/// Returns f64 but computed via f32 to match C's float precision
pub fn sre_random() -> f64 {
    let need_init = SRE_RESEED.with(|r| r.get()) || FIRST_TIME.with(|f| f.get());

    if need_init {
        // Initialize rnd from seed (lines 363-370)
        SRE_RESEED.with(|r| r.set(false));
        FIRST_TIME.with(|f| f.set(false));

        let mut seed = SRE_RANDSEED.with(|s| s.get());
        if seed <= 0 {
            seed = 666; // seeds of zero break me
            SRE_RANDSEED.with(|s| s.set(seed));
        }

        let high1 = seed / DIV;
        let low1 = seed % DIV;
        let high2 = MULT / DIV;
        let low2 = MULT % DIV;

        let rnd_val = (((high2 * low1 + high1 * low2) % DIV) * DIV + low1 * low2) % RANGE;
        RND.with(|r| r.set(rnd_val));
    }

    // Generate next random number (lines 371-373)
    let rnd_val = RND.with(|r| r.get());
    let high1 = rnd_val / DIV;
    let low1 = rnd_val % DIV;
    let high2 = MULT / DIV;
    let low2 = MULT % DIV;

    let new_rnd = (((high2 * low1 + high1 * low2) % DIV) * DIV + low1 * low2) % RANGE;
    RND.with(|r| r.set(new_rnd));

    // Return as float in [0, 1) (line 375)
    // C returns float, so we compute via f32 then convert to f64
    (new_rnd as f32 / RANGE as f32) as f64
}

/// Gaussian random via Box-Muller transform - Gaussrandom()
///
/// Returns N(mean, stddev) distributed random variable.
/// This is a simplified Box-Muller implementation for basic testing.
///
/// The C version (lines 28-144) uses a more complex Ahrens-Dieter method.
/// For initial compatibility testing, we use the simpler Box-Muller transform.
pub fn gaussrandom(mean: f64, stddev: f64) -> f64 {
    let r1 = sre_random();
    let r2 = sre_random();

    let snorm = (-2.0 * r1.ln()).sqrt() * (2.0 * PI * r2).cos();
    stddev * snorm + mean
}

/// Normalize a probability distribution - DNorm()
///
/// Corresponds to DNorm() in sre_math.c (lines 252-265)
/// Returns true if successful, false if sum is zero
pub fn dnorm(vec: &mut [f64]) -> bool {
    let sum: f64 = vec.iter().sum();
    if sum != 0.0 {
        for v in vec.iter_mut() {
            *v /= sum;
        }
        true
    } else {
        false
    }
}

/// Normalize f32 distribution - FNorm()
///
/// Corresponds to FNorm() in sre_math.c (lines 267-280)
/// Returns true if successful, false if sum is zero
pub fn fnorm(vec: &mut [f32]) -> bool {
    let sum: f32 = vec.iter().sum();
    if sum != 0.0 {
        for v in vec.iter_mut() {
            *v /= sum;
        }
        true
    } else {
        false
    }
}

/// Scale a double vector - DScale()
///
/// Corresponds to DScale() in sre_math.c (lines 282-288)
pub fn dscale(vec: &mut [f64], scale: f64) {
    for v in vec.iter_mut() {
        *v *= scale;
    }
}

/// Scale a float vector - FScale()
///
/// Corresponds to FScale() in sre_math.c (lines 289-295)
pub fn fscale(vec: &mut [f32], scale: f32) {
    for v in vec.iter_mut() {
        *v *= scale;
    }
}

/// Set all values in a double vector - DSet()
///
/// Note: Original C code has a bug - it multiplies instead of sets.
/// This implementation matches the C code exactly for compatibility.
/// Corresponds to DSet() in sre_math.c (lines 297-303)
pub fn dset(vec: &mut [f64], value: f64) {
    for v in vec.iter_mut() {
        *v *= value; // Match C code bug
    }
}

/// Set all values in a float vector - FSet()
///
/// Note: Original C code has a bug - it multiplies instead of sets.
/// This implementation matches the C code exactly for compatibility.
/// Corresponds to FSet() in sre_math.c (lines 304-310)
pub fn fset(vec: &mut [f32], value: f32) {
    for v in vec.iter_mut() {
        *v *= value; // Match C code bug
    }
}

/// Corrected version that actually sets values
pub fn dset_corrected(vec: &mut [f64], value: f64) {
    for v in vec.iter_mut() {
        *v = value;
    }
}

/// Corrected version that actually sets values
pub fn fset_corrected(vec: &mut [f32], value: f32) {
    for v in vec.iter_mut() {
        *v = value;
    }
}

/// Sum of a double vector - DSum()
///
/// Corresponds to DSum() in sre_math.c (lines 312-320)
pub fn dsum(vec: &[f64]) -> f64 {
    vec.iter().sum()
}

/// Sum of a float vector - FSum()
///
/// Corresponds to FSum() in sre_math.c (lines 321-329)
pub fn fsum(vec: &[f32]) -> f32 {
    vec.iter().sum()
}

/// Choose from a normalized distribution - DChoose()
///
/// Given a normalized probability distribution, randomly sample an index.
/// Corresponds to DChoose() in sre_math.c (lines 403-418)
pub fn dchoose(p: &[f64]) -> usize {
    let roll = sre_random();
    let mut sum = 0.0;

    for (i, &prob) in p.iter().enumerate() {
        sum += prob;
        if roll < sum {
            return i;
        }
    }

    // Bulletproof fallback
    (sre_random() * p.len() as f64) as usize
}

/// Choose from a normalized distribution - FChoose()
///
/// Given a normalized probability distribution, randomly sample an index.
/// Corresponds to FChoose() in sre_math.c (lines 419-434)
pub fn fchoose(p: &[f32]) -> usize {
    let roll = sre_random() as f32;
    let mut sum = 0.0f32;

    for (i, &prob) in p.iter().enumerate() {
        sum += prob;
        if roll < sum {
            return i;
        }
    }

    // Bulletproof fallback
    (sre_random() * p.len() as f64) as usize
}

/// Log-sum of a log vector (normal space) - DLogSum()
///
/// Calculate sum of values in normal space from log values,
/// returning the log of the sum. Uses the log-sum-exp trick
/// to avoid overflow.
/// Corresponds to DLogSum() in sre_math.c (lines 441-456)
pub fn dlog_sum(logp: &[f64]) -> f64 {
    let max = logp.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mut sum = 0.0;

    for &lp in logp {
        if lp > max - 50.0 {
            sum += (lp - max).exp();
        }
    }

    sum.ln() + max
}

/// Log-sum of a log vector (normal space) - FLogSum()
///
/// Calculate sum of values in normal space from log values,
/// returning the log of the sum.
/// Corresponds to FLogSum() in sre_math.c (lines 457-472)
pub fn flog_sum(logp: &[f32]) -> f32 {
    let max = logp.iter().copied().fold(f32::NEG_INFINITY, f32::max);
    let mut sum = 0.0f32;

    for &lp in logp {
        if lp > max - 50.0 {
            sum += (lp - max).exp();
        }
    }

    sum.ln() + max
}

/// Log-sum of two values in log space
///
/// Computes log(exp(a) + exp(b)) without overflow.
/// More efficient than dlog_sum for just two values.
pub fn log_sum(a: f64, b: f64) -> f64 {
    if a > b {
        a + (1.0 + (b - a).exp()).ln()
    } else {
        b + (1.0 + (a - b).exp()).ln()
    }
}

/// Natural log of gamma function - Gammln()
///
/// Returns ln(Gamma(x)) for x > 0.
/// Adapted from NCBI core math library.
/// Corresponds to Gammln() in sre_math.c (lines 205-243)
pub fn gammln(x: f64) -> f64 {
    // Coefficients from NCBI implementation
    const COF: [f64; 11] = [
        4.694580336184385e+04,
        -1.560605207784446e+05,
        2.065049568014106e+05,
        -1.388934775095388e+05,
        5.031796415085709e+04,
        -9.601592329182778e+03,
        8.785855930895250e+02,
        -3.155153906098611e+01,
        2.908143421162229e-01,
        -2.319827630494973e-04,
        1.251639670050933e-10,
    ];

    // Handle x <= 0 case (severe hack but effective)
    if x <= 0.0 {
        return 999999.0;
    }

    let xx = x - 1.0;
    let tx = xx + 11.0;
    let mut tmp = tx;
    let mut value = 1.0;

    // Sum least significant terms first
    for coef in COF.iter().rev() {
        value += coef / tmp;
        tmp -= 1.0;
    }

    value = value.ln();
    let tx2 = tx + 0.5;
    value += 0.918938533 + (xx + 0.5) * tx2.ln() - tx2;

    value
}

/// Log gamma function (alias for gammln)
pub fn log_gamma(x: f64) -> f64 {
    gammln(x)
}

/// Linear regression fit - Linefit()
///
/// Given points (x, y), fits to y = a + bx.
/// Returns (intercept, slope, correlation_coefficient).
/// Corresponds to Linefit() in sre_math.c (lines 163-192)
pub fn linefit(x: &[f64], y: &[f64]) -> Option<(f64, f64, f64)> {
    let n = x.len();
    if n == 0 || n != y.len() {
        return None;
    }

    let n_f64 = n as f64;

    // Calculate averages
    let xbar: f64 = x.iter().sum::<f64>() / n_f64;
    let ybar: f64 = y.iter().sum::<f64>() / n_f64;

    // Calculate sums
    let mut sxx = 0.0;
    let mut syy = 0.0;
    let mut sxy = 0.0;

    for i in 0..n {
        let dx = x[i] - xbar;
        let dy = y[i] - ybar;
        sxx += dx * dx;
        // Note: Original C code has a bug here: syy += (y[i] - ybar) * (x[i] - xbar)
        // We implement the corrected version
        syy += dy * dy;
        sxy += dx * dy;
    }

    if sxx == 0.0 {
        return None;
    }

    let b = sxy / sxx;
    let a = ybar - xbar * b;
    let r = sxy / (sxx.sqrt() * syy.sqrt());

    Some((a, b, r))
}

/// Gaussian probability density function
///
/// Returns the probability density at x for a Gaussian
/// with given mean and standard deviation.
pub fn gaussian_pdf(x: f64, mean: f64, stddev: f64) -> f64 {
    let diff = x - mean;
    let exponent = -(diff * diff) / (2.0 * stddev * stddev);
    (exponent.exp()) / (stddev * (2.0 * PI).sqrt())
}

/// Incomplete gamma function P(a, x) using series expansion
///
/// Used for chi-squared p-value calculation.
/// Based on Numerical Recipes algorithm.
pub fn incomplete_gamma(a: f64, x: f64) -> f64 {
    if x < 0.0 || a <= 0.0 {
        return 0.0;
    }

    if x == 0.0 {
        return 0.0;
    }

    let gln = gammln(a);

    if x < a + 1.0 {
        // Series representation
        let mut sum = 1.0 / a;
        let mut del = sum;
        let mut ap = a;

        for _ in 0..100 {
            ap += 1.0;
            del *= x / ap;
            sum += del;
            if del.abs() < sum.abs() * 1e-10 {
                break;
            }
        }

        sum * (-x + a * x.ln() - gln).exp()
    } else {
        // Continued fraction representation
        let mut b = x + 1.0 - a;
        let mut c = 1.0 / 1e-30;
        let mut d = 1.0 / b;
        let mut h = d;

        for i in 1..=100 {
            let an = -(i as f64) * (i as f64 - a);
            b += 2.0;
            d = an * d + b;
            if d.abs() < 1e-30 {
                d = 1e-30;
            }
            c = b + an / c;
            if c.abs() < 1e-30 {
                c = 1e-30;
            }
            d = 1.0 / d;
            let del = d * c;
            h *= del;
            if (del - 1.0).abs() < 1e-10 {
                break;
            }
        }

        1.0 - h * (-x + a * x.ln() - gln).exp()
    }
}

/// Chi-squared p-value
///
/// Returns the probability of observing a chi-squared value
/// at least as extreme as x2 with df degrees of freedom.
pub fn chi_squared_pvalue(x2: f64, df: usize) -> f64 {
    if x2 <= 0.0 || df == 0 {
        return 1.0;
    }

    let a = df as f64 / 2.0;
    let x = x2 / 2.0;

    1.0 - incomplete_gamma(a, x)
}

/// Sample from a Dirichlet distribution
///
/// Given alpha parameters, returns a sample from the Dirichlet.
/// Requires a sum function for the gamma samples.
pub fn sample_dirichlet(alpha: &[f64]) -> Vec<f64> {
    let mut samples: Vec<f64> = alpha
        .iter()
        .map(|&a| {
            // Generate gamma(a, 1) sample using Marsaglia method
            let d = a - 1.0 / 3.0;
            let c = 1.0 / (9.0 * d).sqrt();

            loop {
                let x = gaussrandom(0.0, 1.0);
                let v = 1.0 + c * x;
                if v > 0.0 {
                    let v3 = v * v * v;
                    let u = sre_random();
                    if u < 1.0 - 0.0331 * (x * x) * (x * x)
                        || u.ln() < 0.5 * x * x + d * (1.0 - v3 + v3.ln())
                    {
                        return d * v3;
                    }
                }
            }
        })
        .collect();

    let sum: f64 = samples.iter().sum();
    for s in &mut samples {
        *s /= sum;
    }

    samples
}

/// Dirichlet log probability
///
/// Returns the log probability of observing counts given alpha parameters.
pub fn dirichlet_lnp(alpha: &[f64], counts: &[f64]) -> f64 {
    if alpha.len() != counts.len() {
        return f64::NEG_INFINITY;
    }

    let alpha_sum: f64 = alpha.iter().sum();
    let count_sum: f64 = counts.iter().sum();

    let mut lnp = gammln(alpha_sum) - gammln(alpha_sum + count_sum);

    for i in 0..alpha.len() {
        lnp += gammln(alpha[i] + counts[i]) - gammln(alpha[i]);
    }

    lnp
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sre_random_basic() {
        // Test with seed 666 (default)
        sre_srandom(666);

        let r1 = sre_random();
        let r2 = sre_random();

        // Values should be in [0, 1)
        assert!(r1 >= 0.0 && r1 < 1.0);
        assert!(r2 >= 0.0 && r2 < 1.0);

        // Should be different
        assert_ne!(r1, r2);
    }

    #[test]
    fn test_sre_random_reproducible() {
        // Same seed should give same sequence
        sre_srandom(12345);
        let seq1: Vec<f64> = (0..10).map(|_| sre_random()).collect();

        sre_srandom(12345);
        let seq2: Vec<f64> = (0..10).map(|_| sre_random()).collect();

        assert_eq!(seq1, seq2);
    }

    #[test]
    fn test_dnorm() {
        let mut vec = [1.0, 2.0, 3.0, 4.0];
        assert!(dnorm(&mut vec));

        let sum: f64 = vec.iter().sum();
        assert!((sum - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_fnorm() {
        let mut vec = [1.0f32, 2.0, 3.0, 4.0];
        assert!(fnorm(&mut vec));

        let sum: f32 = vec.iter().sum();
        assert!((sum - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_dscale() {
        let mut vec = [1.0, 2.0, 3.0];
        dscale(&mut vec, 2.0);

        assert_eq!(vec, [2.0, 4.0, 6.0]);
    }
}
