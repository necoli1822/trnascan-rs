//! CM model file I/O (save.c)
//!
//! This module implements reading and writing covariance models to disk.
//! Supports both flat text format (version 2.0) and binary format.

use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::Path;

use crate::types::cm::CM;
use crate::types::constants::*;

/// Read a CM model from file (text or binary format)
///
/// Implements ReadCM() from save.c lines 138-197
///
/// # Format Detection
/// - Binary: First 4 bytes = 0xe3edb2b0 (V20_MAGIC)
/// - Text: First line = "### cove V2"
///
/// # Arguments
/// * `filename` - Path to CM model file
///
/// # Returns
/// * `Ok(CM)` - Successfully loaded model
/// * `Err(String)` - Error message
pub fn read_cm<P: AsRef<Path>>(filename: P) -> Result<CM, String> {
    let path = filename.as_ref();
    let mut file = File::open(path)
        .map_err(|e| format!("Cannot open model file {:?}: {}", path, e))?;

    // Read first 4 bytes to detect format
    let mut magic_bytes = [0u8; 4];
    file.read_exact(&mut magic_bytes)
        .map_err(|e| format!("Failed to read magic number: {}", e))?;

    let magic_number = u32::from_le_bytes(magic_bytes);

    if magic_number == V20_MAGIC {
        // Binary format
        read_binary_cm(file)
    } else {
        // Text format - rewind and check header
        let file = File::open(path)
            .map_err(|e| format!("Cannot reopen file: {}", e))?;
        let reader = BufReader::new(file);
        let mut lines = reader.lines();

        // Read first line
        let first_line = lines
            .next()
            .ok_or_else(|| "Empty model file".to_string())?
            .map_err(|e| format!("Failed to read header: {}", e))?;

        if first_line.starts_with("### cove V2") {
            read_text_cm(lines)
        } else {
            Err(format!("Unrecognized model format: {}", first_line))
        }
    }
}

/// Read text format CM (read_cm20 from save.c lines 200-266)
///
/// Format:
/// ```text
/// ### cove V2
/// N     nodes
/// ### node 0 type T
/// nxt nxt2
/// [6x6 transition matrix]
/// [4 INSL emissions] # INSL
/// [4 INSR emissions] # INSR
/// [4x4 MATP emissions] # MATP (4 lines)
/// [4 MATL emissions] # MATL
/// [4 MATR emissions] # MATR
/// ```
fn read_text_cm<I>(mut lines: I) -> Result<CM, String>
where
    I: Iterator<Item = io::Result<String>>,
{
    // Read number of nodes
    let nodes_line = lines
        .next()
        .ok_or_else(|| "Missing nodes line".to_string())?
        .map_err(|e| format!("Read error: {}", e))?;

    let nodes: usize = nodes_line
        .split_whitespace()
        .next()
        .ok_or_else(|| "Invalid nodes line".to_string())?
        .parse()
        .map_err(|e| format!("Invalid node count: {}", e))?;

    let mut cm = CM::new(nodes);

    // Read each node
    for k in 0..nodes {
        // Read node header: "### node N type T"
        let node_header = lines
            .next()
            .ok_or_else(|| format!("Missing node {} header", k))?
            .map_err(|e| format!("Read error: {}", e))?;

        let node_type = node_header
            .split_whitespace()
            .last()
            .ok_or_else(|| format!("Invalid node {} header", k))?
            .parse::<i32>()
            .map_err(|e| format!("Invalid node type: {}", e))?;

        cm.nd[k].node_type = node_type;

        // Read nxt and nxt2
        let nxt_line = lines
            .next()
            .ok_or_else(|| format!("Missing nxt line for node {}", k))?
            .map_err(|e| format!("Read error: {}", e))?;

        let nxt_parts: Vec<&str> = nxt_line.split_whitespace().collect();
        if nxt_parts.len() < 2 {
            return Err(format!("Invalid nxt line for node {}", k));
        }

        cm.nd[k].nxt = nxt_parts[0]
            .parse()
            .map_err(|e| format!("Invalid nxt: {}", e))?;
        cm.nd[k].nxt2 = nxt_parts[1]
            .parse()
            .map_err(|e| format!("Invalid nxt2: {}", e))?;

        // Read transition matrix (6x6)
        for i in 0..STATETYPES {
            let tmx_line = lines
                .next()
                .ok_or_else(|| format!("Missing tmx line {} for node {}", i, k))?
                .map_err(|e| format!("Read error: {}", e))?;

            let values: Result<Vec<f64>, _> = tmx_line
                .split_whitespace()
                .take(STATETYPES)
                .map(|s| s.parse::<f64>())
                .collect();

            let values = values.map_err(|e| format!("Invalid transition value: {}", e))?;
            if values.len() != STATETYPES {
                return Err(format!("Invalid transition row length for node {}", k));
            }

            for j in 0..STATETYPES {
                cm.nd[k].tmx[i][j] = values[j];
            }
        }

        // Read INSL emissions
        let insl_line = lines
            .next()
            .ok_or_else(|| format!("Missing INSL line for node {}", k))?
            .map_err(|e| format!("Read error: {}", e))?;

        let insl_values: Result<Vec<f64>, _> = insl_line
            .split_whitespace()
            .take(ALPHASIZE)
            .map(|s| s.parse::<f64>())
            .collect();

        let insl_values = insl_values.map_err(|e| format!("Invalid INSL value: {}", e))?;
        for i in 0..ALPHASIZE {
            cm.nd[k].il_emit[i] = insl_values[i];
        }

        // Read INSR emissions
        let insr_line = lines
            .next()
            .ok_or_else(|| format!("Missing INSR line for node {}", k))?
            .map_err(|e| format!("Read error: {}", e))?;

        let insr_values: Result<Vec<f64>, _> = insr_line
            .split_whitespace()
            .take(ALPHASIZE)
            .map(|s| s.parse::<f64>())
            .collect();

        let insr_values = insr_values.map_err(|e| format!("Invalid INSR value: {}", e))?;
        for i in 0..ALPHASIZE {
            cm.nd[k].ir_emit[i] = insr_values[i];
        }

        // Read MATP emissions (4 lines)
        for i in 0..ALPHASIZE {
            let matp_line = lines
                .next()
                .ok_or_else(|| format!("Missing MATP line {} for node {}", i, k))?
                .map_err(|e| format!("Read error: {}", e))?;

            let matp_values: Result<Vec<f64>, _> = matp_line
                .split_whitespace()
                .take(ALPHASIZE)
                .map(|s| s.parse::<f64>())
                .collect();

            let matp_values = matp_values.map_err(|e| format!("Invalid MATP value: {}", e))?;
            for j in 0..ALPHASIZE {
                cm.nd[k].mp_emit[i][j] = matp_values[j];
            }
        }

        // Read MATL emissions
        let matl_line = lines
            .next()
            .ok_or_else(|| format!("Missing MATL line for node {}", k))?
            .map_err(|e| format!("Read error: {}", e))?;

        let matl_values: Result<Vec<f64>, _> = matl_line
            .split_whitespace()
            .take(ALPHASIZE)
            .map(|s| s.parse::<f64>())
            .collect();

        let matl_values = matl_values.map_err(|e| format!("Invalid MATL value: {}", e))?;
        for i in 0..ALPHASIZE {
            cm.nd[k].ml_emit[i] = matl_values[i];
        }

        // Read MATR emissions
        let matr_line = lines
            .next()
            .ok_or_else(|| format!("Missing MATR line for node {}", k))?
            .map_err(|e| format!("Read error: {}", e))?;

        let matr_values: Result<Vec<f64>, _> = matr_line
            .split_whitespace()
            .take(ALPHASIZE)
            .map(|s| s.parse::<f64>())
            .collect();

        let matr_values = matr_values.map_err(|e| format!("Invalid MATR value: {}", e))?;
        for i in 0..ALPHASIZE {
            cm.nd[k].mr_emit[i] = matr_values[i];
        }
    }

    Ok(cm)
}

/// Read binary format CM (read_bincm20 from save.c lines 270-315)
///
/// Binary format (little-endian):
/// - 4 bytes: magic number (0xe3edb2b0) [already read]
/// - 4 bytes: number of nodes
/// - For each node:
///   - 4 bytes: node type
///   - 4 bytes: nxt
///   - 4 bytes: nxt2
///   - 36 * 8 bytes: transition matrix (6x6 doubles)
///   - 4 * 8 bytes: INSL emissions
///   - 4 * 8 bytes: INSR emissions
///   - 16 * 8 bytes: MATP emissions (4x4 doubles)
///   - 4 * 8 bytes: MATL emissions
///   - 4 * 8 bytes: MATR emissions
fn read_binary_cm(mut file: File) -> Result<CM, String> {
    // Read number of nodes
    let mut nodes_bytes = [0u8; 4];
    file.read_exact(&mut nodes_bytes)
        .map_err(|e| format!("Failed to read node count: {}", e))?;
    let nodes = i32::from_le_bytes(nodes_bytes) as usize;

    let mut cm = CM::new(nodes);

    for k in 0..nodes {
        // Read node type
        let mut type_bytes = [0u8; 4];
        file.read_exact(&mut type_bytes)
            .map_err(|e| format!("Failed to read node {} type: {}", k, e))?;
        cm.nd[k].node_type = i32::from_le_bytes(type_bytes);

        // Read nxt
        let mut nxt_bytes = [0u8; 4];
        file.read_exact(&mut nxt_bytes)
            .map_err(|e| format!("Failed to read node {} nxt: {}", k, e))?;
        cm.nd[k].nxt = i32::from_le_bytes(nxt_bytes);

        // Read nxt2
        let mut nxt2_bytes = [0u8; 4];
        file.read_exact(&mut nxt2_bytes)
            .map_err(|e| format!("Failed to read node {} nxt2: {}", k, e))?;
        cm.nd[k].nxt2 = i32::from_le_bytes(nxt2_bytes);

        // Read transition matrix (6x6 = 36 doubles)
        for i in 0..STATETYPES {
            for j in 0..STATETYPES {
                let mut tmx_bytes = [0u8; 8];
                file.read_exact(&mut tmx_bytes)
                    .map_err(|e| format!("Failed to read tmx[{}][{}]: {}", i, j, e))?;
                cm.nd[k].tmx[i][j] = f64::from_le_bytes(tmx_bytes);
            }
        }

        // Read INSL emissions
        for i in 0..ALPHASIZE {
            let mut il_bytes = [0u8; 8];
            file.read_exact(&mut il_bytes)
                .map_err(|e| format!("Failed to read il_emit[{}]: {}", i, e))?;
            cm.nd[k].il_emit[i] = f64::from_le_bytes(il_bytes);
        }

        // Read INSR emissions
        for i in 0..ALPHASIZE {
            let mut ir_bytes = [0u8; 8];
            file.read_exact(&mut ir_bytes)
                .map_err(|e| format!("Failed to read ir_emit[{}]: {}", i, e))?;
            cm.nd[k].ir_emit[i] = f64::from_le_bytes(ir_bytes);
        }

        // Read MATP emissions (4x4 = 16 doubles)
        for i in 0..ALPHASIZE {
            for j in 0..ALPHASIZE {
                let mut mp_bytes = [0u8; 8];
                file.read_exact(&mut mp_bytes)
                    .map_err(|e| format!("Failed to read mp_emit[{}][{}]: {}", i, j, e))?;
                cm.nd[k].mp_emit[i][j] = f64::from_le_bytes(mp_bytes);
            }
        }

        // Read MATL emissions
        for i in 0..ALPHASIZE {
            let mut ml_bytes = [0u8; 8];
            file.read_exact(&mut ml_bytes)
                .map_err(|e| format!("Failed to read ml_emit[{}]: {}", i, e))?;
            cm.nd[k].ml_emit[i] = f64::from_le_bytes(ml_bytes);
        }

        // Read MATR emissions
        for i in 0..ALPHASIZE {
            let mut mr_bytes = [0u8; 8];
            file.read_exact(&mut mr_bytes)
                .map_err(|e| format!("Failed to read mr_emit[{}]: {}", i, e))?;
            cm.nd[k].mr_emit[i] = f64::from_le_bytes(mr_bytes);
        }
    }

    Ok(cm)
}

/// Write CM to text format (WriteCM from save.c lines 25-80)
///
/// # Arguments
/// * `filename` - Path to output file
/// * `cm` - CM model to write
///
/// # Returns
/// * `Ok(())` - Success
/// * `Err(String)` - Error message
pub fn write_cm<P: AsRef<Path>>(filename: P, cm: &CM) -> Result<(), String> {
    let mut file = File::create(filename.as_ref())
        .map_err(|e| format!("Cannot create file: {}", e))?;

    // Header
    writeln!(file, "### cove V2").map_err(|e| format!("Write error: {}", e))?;
    writeln!(file, "{} \tnodes", cm.nodes).map_err(|e| format!("Write error: {}", e))?;

    // Write each node
    for k in 0..cm.nodes {
        writeln!(
            file,
            "### node {} type {}",
            k, cm.nd[k].node_type
        )
        .map_err(|e| format!("Write error: {}", e))?;
        writeln!(file, "{}  {}", cm.nd[k].nxt, cm.nd[k].nxt2)
            .map_err(|e| format!("Write error: {}", e))?;

        // Transition matrix
        for i in 0..STATETYPES {
            for j in 0..STATETYPES {
                write!(file, "{:.5} ", cm.nd[k].tmx[i][j])
                    .map_err(|e| format!("Write error: {}", e))?;
            }
            writeln!(file).map_err(|e| format!("Write error: {}", e))?;
        }

        // INSL emissions
        for i in 0..ALPHASIZE {
            write!(file, "{:.5} ", cm.nd[k].il_emit[i])
                .map_err(|e| format!("Write error: {}", e))?;
        }
        writeln!(file, "# INSL").map_err(|e| format!("Write error: {}", e))?;

        // INSR emissions
        for i in 0..ALPHASIZE {
            write!(file, "{:.5} ", cm.nd[k].ir_emit[i])
                .map_err(|e| format!("Write error: {}", e))?;
        }
        writeln!(file, "# INSR").map_err(|e| format!("Write error: {}", e))?;

        // MATP emissions
        for i in 0..ALPHASIZE {
            for j in 0..ALPHASIZE {
                write!(file, "{:.5} ", cm.nd[k].mp_emit[i][j])
                    .map_err(|e| format!("Write error: {}", e))?;
            }
            writeln!(file, "# MATP").map_err(|e| format!("Write error: {}", e))?;
        }

        // MATL emissions
        for i in 0..ALPHASIZE {
            write!(file, "{:.5} ", cm.nd[k].ml_emit[i])
                .map_err(|e| format!("Write error: {}", e))?;
        }
        writeln!(file, "# MATL").map_err(|e| format!("Write error: {}", e))?;

        // MATR emissions
        for i in 0..ALPHASIZE {
            write!(file, "{:.5} ", cm.nd[k].mr_emit[i])
                .map_err(|e| format!("Write error: {}", e))?;
        }
        writeln!(file, "# MATR").map_err(|e| format!("Write error: {}", e))?;
    }

    Ok(())
}

/// Write CM to binary format (WriteBinaryCM from save.c lines 84-133)
pub fn write_binary_cm<P: AsRef<Path>>(filename: P, cm: &CM) -> Result<(), String> {
    let mut file = File::create(filename.as_ref())
        .map_err(|e| format!("Cannot create file: {}", e))?;

    // Write magic number
    file.write_all(&V20_MAGIC.to_le_bytes())
        .map_err(|e| format!("Write error: {}", e))?;

    // Write number of nodes
    file.write_all(&(cm.nodes as i32).to_le_bytes())
        .map_err(|e| format!("Write error: {}", e))?;

    // Write each node
    for k in 0..cm.nodes {
        file.write_all(&cm.nd[k].node_type.to_le_bytes())
            .map_err(|e| format!("Write error: {}", e))?;
        file.write_all(&cm.nd[k].nxt.to_le_bytes())
            .map_err(|e| format!("Write error: {}", e))?;
        file.write_all(&cm.nd[k].nxt2.to_le_bytes())
            .map_err(|e| format!("Write error: {}", e))?;

        // Transition matrix
        for i in 0..STATETYPES {
            for j in 0..STATETYPES {
                file.write_all(&cm.nd[k].tmx[i][j].to_le_bytes())
                    .map_err(|e| format!("Write error: {}", e))?;
            }
        }

        // Emissions
        for i in 0..ALPHASIZE {
            file.write_all(&cm.nd[k].il_emit[i].to_le_bytes())
                .map_err(|e| format!("Write error: {}", e))?;
        }
        for i in 0..ALPHASIZE {
            file.write_all(&cm.nd[k].ir_emit[i].to_le_bytes())
                .map_err(|e| format!("Write error: {}", e))?;
        }
        for i in 0..ALPHASIZE {
            for j in 0..ALPHASIZE {
                file.write_all(&cm.nd[k].mp_emit[i][j].to_le_bytes())
                    .map_err(|e| format!("Write error: {}", e))?;
            }
        }
        for i in 0..ALPHASIZE {
            file.write_all(&cm.nd[k].ml_emit[i].to_le_bytes())
                .map_err(|e| format!("Write error: {}", e))?;
        }
        for i in 0..ALPHASIZE {
            file.write_all(&cm.nd[k].mr_emit[i].to_le_bytes())
                .map_err(|e| format!("Write error: {}", e))?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;

    #[test]
    fn test_write_read_text_roundtrip() {
        let mut cm = CM::new(2);
        cm.nd[0].node_type = 1;
        cm.nd[0].nxt = 1;
        cm.nd[0].nxt2 = -1;
        cm.nd[0].tmx[0][0] = 0.5;
        cm.nd[0].il_emit[0] = 0.25;

        cm.nd[1].node_type = 2;
        cm.nd[1].nxt = -1;
        cm.nd[1].nxt2 = -1;

        let tmpfile = NamedTempFile::new().unwrap();
        write_cm(tmpfile.path(), &cm).unwrap();

        let cm2 = read_cm(tmpfile.path()).unwrap();
        assert_eq!(cm2.nodes, 2);
        assert_eq!(cm2.nd[0].node_type, 1);
        assert_eq!(cm2.nd[0].nxt, 1);
        assert_eq!(cm2.nd[0].tmx[0][0], 0.5);
    }

    #[test]
    fn test_write_read_binary_roundtrip() {
        let mut cm = CM::new(1);
        cm.nd[0].node_type = 3;
        cm.nd[0].nxt = 0;
        cm.nd[0].nxt2 = 0;
        cm.nd[0].tmx[2][3] = 0.125;
        cm.nd[0].mp_emit[1][2] = 0.0625;

        let tmpfile = NamedTempFile::new().unwrap();
        write_binary_cm(tmpfile.path(), &cm).unwrap();

        let cm2 = read_cm(tmpfile.path()).unwrap();
        assert_eq!(cm2.nodes, 1);
        assert_eq!(cm2.nd[0].node_type, 3);
        assert_eq!(cm2.nd[0].tmx[2][3], 0.125);
        assert_eq!(cm2.nd[0].mp_emit[1][2], 0.0625);
    }
}
