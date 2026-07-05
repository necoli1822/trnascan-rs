use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Compare two f64 values with epsilon tolerance
pub fn float_eq(a: f64, b: f64, epsilon: f64) -> bool {
    (a - b).abs() < epsilon
}

/// Compare two f32 values with epsilon tolerance
pub fn float32_eq(a: f32, b: f32, epsilon: f32) -> bool {
    (a - b).abs() < epsilon
}

/// Compare two f64 slices element-wise with epsilon tolerance
pub fn vec_float_eq(a: &[f64], b: &[f64], epsilon: f64) -> bool {
    if a.len() != b.len() {
        return false;
    }
    a.iter().zip(b.iter()).all(|(x, y)| float_eq(*x, *y, epsilon))
}

/// Compare two i32 slices for exact equality
pub fn i32_vec_eq(a: &[i32], b: &[i32]) -> bool {
    a == b
}

/// Load golden values from a file
/// Returns a vector of trimmed, non-empty lines
pub fn load_golden_values(path: &str) -> Vec<String> {
    let file = File::open(path).expect(&format!("Failed to open golden file: {}", path));
    let reader = BufReader::new(file);

    reader
        .lines()
        .filter_map(|line| {
            line.ok().and_then(|l| {
                let trimmed = l.trim().to_string();
                if trimmed.is_empty() || trimmed.starts_with('#') {
                    None
                } else {
                    Some(trimmed)
                }
            })
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_float_eq() {
        assert!(float_eq(1.0, 1.0, 1e-10));
        assert!(float_eq(1.0, 1.0000000001, 1e-8));
        assert!(!float_eq(1.0, 1.1, 1e-10));
    }

    #[test]
    fn test_float32_eq() {
        assert!(float32_eq(1.0f32, 1.0f32, 1e-6));
        assert!(float32_eq(1.0f32, 1.000001f32, 1e-5));
        assert!(!float32_eq(1.0f32, 1.1f32, 1e-6));
    }

    #[test]
    fn test_vec_float_eq() {
        let a = vec![1.0, 2.0, 3.0];
        let b = vec![1.0, 2.0, 3.0];
        let c = vec![1.0, 2.0, 3.1];
        let d = vec![1.0, 2.0];

        assert!(vec_float_eq(&a, &b, 1e-10));
        assert!(!vec_float_eq(&a, &c, 1e-10));
        assert!(!vec_float_eq(&a, &d, 1e-10));
    }

    #[test]
    fn test_i32_vec_eq() {
        let a = vec![1, 2, 3];
        let b = vec![1, 2, 3];
        let c = vec![1, 2, 4];

        assert!(i32_vec_eq(&a, &b));
        assert!(!i32_vec_eq(&a, &c));
    }
}
