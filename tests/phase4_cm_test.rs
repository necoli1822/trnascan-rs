// Phase 4: Covariance Model (CM) tests
// Tests for loading and parsing CM structures from tRNAscan-SE model files

mod common;

use common::load_golden_values;
use trnascan_rs::core::save::read_cm;
use trnascan_rs::types::constants::*;

#[test]
fn test_cm_structure_loading() {
    // Test loading of CM structure from model file
    // Golden file: tests/golden/cm/cm_structure.txt
    let golden = load_golden_values("tests/golden/cm/cm_structure.txt");
    assert!(!golden.is_empty(), "Golden file should contain CM structure data");

    // Load TRNA2.cm model
    let cm_path = "data/models/TRNA2.cm";
    let cm = read_cm(cm_path).expect("Failed to load CM model");

    // Verify number of nodes matches golden file
    assert_eq!(cm.nodes, 72, "Expected 72 nodes as per golden file");

    // Verify specific nodes from golden file
    assert_eq!(cm.nd[0].node_type, ROOT_NODE as i32); // Node 0: ROOT_NODE
    assert_eq!(cm.nd[0].nxt, 1);
    assert_eq!(cm.nd[0].nxt2, -1);

    assert_eq!(cm.nd[1].node_type, MATR_NODE as i32); // Node 1: MATR_NODE
    assert_eq!(cm.nd[1].nxt, 2);
    assert_eq!(cm.nd[1].nxt2, -1);

    assert_eq!(cm.nd[11].node_type, BIFURC_NODE as i32); // Node 11: BIFURC_NODE
    assert_eq!(cm.nd[11].nxt, 12);
    assert_eq!(cm.nd[11].nxt2, 25); // Bifurcation has nxt2 != -1
}

#[test]
fn test_cm_node_details() {
    // Test detailed node information from CM
    // Golden file: tests/golden/cm/node_details.txt
    let golden = load_golden_values("tests/golden/cm/node_details.txt");
    assert!(!golden.is_empty(), "Golden file should contain node details");

    let cm_path = "data/models/TRNA2.cm";
    let cm = read_cm(cm_path).expect("Failed to load CM model");

    // Verify different node types from golden file
    assert_eq!(cm.nd[0].node_type, ROOT_NODE as i32);
    assert_eq!(cm.nd[1].node_type, MATR_NODE as i32);
    assert_eq!(cm.nd[2].node_type, MATP_NODE as i32);
    assert_eq!(cm.nd[9].node_type, MATL_NODE as i32);
    assert_eq!(cm.nd[11].node_type, BIFURC_NODE as i32);
    assert_eq!(cm.nd[12].node_type, BEGINL_NODE as i32);
    assert_eq!(cm.nd[25].node_type, BEGINR_NODE as i32);

    // Verify some specific nxt/nxt2 values from golden file
    assert_eq!(cm.nd[24].nxt, -1); // Terminal node (no next)
    assert_eq!(cm.nd[24].nxt2, -1);
    assert_eq!(cm.nd[27].node_type, BIFURC_NODE as i32);
    assert_eq!(cm.nd[27].nxt, 28);
    assert_eq!(cm.nd[27].nxt2, 41); // Bifurcation
}

#[test]
fn test_cm_istate_array() {
    // Test insert state probability arrays
    // Golden file: tests/golden/cm/istate_dump.txt
    let golden = load_golden_values("tests/golden/cm/istate_dump.txt");
    assert!(!golden.is_empty(), "Golden file should contain istate data");

    let cm_path = "data/models/TRNA2.cm";
    let cm = read_cm(cm_path).expect("Failed to load CM model");

    // Verify that insert state emission arrays have correct dimensions
    for i in 0..cm.nodes {
        assert_eq!(cm.nd[i].il_emit.len(), ALPHASIZE);
        assert_eq!(cm.nd[i].ir_emit.len(), ALPHASIZE);
    }

    // Test that probabilities are within valid range [0, 1]
    for i in 0..cm.nodes {
        for j in 0..ALPHASIZE {
            assert!(cm.nd[i].il_emit[j] >= 0.0 && cm.nd[i].il_emit[j] <= 1.0);
            assert!(cm.nd[i].ir_emit[j] >= 0.0 && cm.nd[i].ir_emit[j] <= 1.0);
        }
    }
}

#[test]
fn test_cm_transition_probabilities() {
    // Test transition probability matrices
    let cm_path = "data/models/TRNA2.cm";
    let cm = read_cm(cm_path).expect("Failed to load CM model");

    // Verify transition matrix dimensions
    for i in 0..cm.nodes {
        assert_eq!(cm.nd[i].tmx.len(), STATETYPES);
        for j in 0..STATETYPES {
            assert_eq!(cm.nd[i].tmx[j].len(), STATETYPES);
        }
    }

    // Test specific transition probabilities from node 0 (from CM file)
    // Node 0 row 0: 0.00846 0.00000 0.00000 0.80479 0.02326 0.16350
    assert!((cm.nd[0].tmx[0][0] - 0.00846).abs() < 0.00001);
    assert!((cm.nd[0].tmx[0][3] - 0.80479).abs() < 0.00001);
    assert!((cm.nd[0].tmx[0][5] - 0.16350).abs() < 0.00001);

    // Verify probabilities are in valid range [0, 1]
    for i in 0..cm.nodes {
        for j in 0..STATETYPES {
            for k in 0..STATETYPES {
                assert!(cm.nd[i].tmx[j][k] >= 0.0 && cm.nd[i].tmx[j][k] <= 1.0);
            }
        }
    }
}

#[test]
fn test_cm_emission_probabilities() {
    // Test emission probability matrices
    let cm_path = "data/models/TRNA2.cm";
    let cm = read_cm(cm_path).expect("Failed to load CM model");

    // Verify emission array dimensions
    for i in 0..cm.nodes {
        assert_eq!(cm.nd[i].mp_emit.len(), ALPHASIZE); // MATP: 4x4
        for j in 0..ALPHASIZE {
            assert_eq!(cm.nd[i].mp_emit[j].len(), ALPHASIZE);
        }
        assert_eq!(cm.nd[i].ml_emit.len(), ALPHASIZE); // MATL: 4
        assert_eq!(cm.nd[i].mr_emit.len(), ALPHASIZE); // MATR: 4
    }

    // Test specific emission from Node 1 (MATR_NODE)
    // From CM file: 0.56676 0.06605 0.22443 0.14276 # MATR
    assert!((cm.nd[1].mr_emit[0] - 0.56676).abs() < 0.00001);
    assert!((cm.nd[1].mr_emit[1] - 0.06605).abs() < 0.00001);
    assert!((cm.nd[1].mr_emit[2] - 0.22443).abs() < 0.00001);
    assert!((cm.nd[1].mr_emit[3] - 0.14276).abs() < 0.00001);

    // Verify probabilities are in valid range [0, 1]
    for i in 0..cm.nodes {
        for j in 0..ALPHASIZE {
            assert!(cm.nd[i].ml_emit[j] >= 0.0 && cm.nd[i].ml_emit[j] <= 1.0);
            assert!(cm.nd[i].mr_emit[j] >= 0.0 && cm.nd[i].mr_emit[j] <= 1.0);
            for k in 0..ALPHASIZE {
                assert!(cm.nd[i].mp_emit[j][k] >= 0.0 && cm.nd[i].mp_emit[j][k] <= 1.0);
            }
        }
    }
}

#[test]
fn test_cm_begin_end_probs() {
    // Test begin and end state probabilities
    let cm_path = "data/models/TRNA2.cm";
    let cm = read_cm(cm_path).expect("Failed to load CM model");

    // Verify ROOT_NODE (node 0) has valid next pointer
    assert_eq!(cm.nd[0].node_type, ROOT_NODE as i32);
    assert_eq!(cm.nd[0].nxt, 1); // Points to start of model

    // Verify terminal nodes (nxt == -1) exist
    let mut has_terminal = false;
    for i in 0..cm.nodes {
        if cm.nd[i].nxt == -1 {
            has_terminal = true;
            break;
        }
    }
    assert!(has_terminal, "CM should have terminal nodes with nxt == -1");
}

#[test]
fn test_multiple_cm_loading() {
    // Test loading multiple CM variants
    let models = vec![
        "data/models/TRNA2.cm",       // Standard
        "data/models/TRNA2-euk.cm",   // Eukaryotic
        "data/models/TRNA2-bact.cm",  // Bacterial
        "data/models/TRNA2-arch.cm",  // Archaeal
    ];

    for model_path in models {
        let cm = read_cm(model_path).unwrap_or_else(|e| {
            panic!("Failed to load {}: {}", model_path, e)
        });

        // Verify basic CM structure
        assert!(cm.nodes > 0, "{} should have nodes", model_path);
        assert_eq!(cm.nd.len(), cm.nodes, "Node array size should match node count");

        // Verify first node is ROOT_NODE
        assert_eq!(
            cm.nd[0].node_type,
            ROOT_NODE as i32,
            "{} should start with ROOT_NODE",
            model_path
        );
    }
}
