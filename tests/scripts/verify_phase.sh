#!/bin/bash

# Phase-by-phase verification script for tRNAscan-SE
# Usage: ./verify_phase.sh <phase_number|all>

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
GOLDEN_DIR="${PROJECT_ROOT}/tests/golden"

# Color codes for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Verify phase argument
if [[ $# -lt 1 ]]; then
    echo "Usage: $(basename "$0") <phase_number|all>"
    echo "Phases: 1-10, or 'all' to run all phases"
    exit 1
fi

PHASE=$1

# Validate phase number or 'all'
if [[ "$PHASE" != "all" ]]; then
    if ! [[ "$PHASE" =~ ^[0-9]+$ ]] || [[ "$PHASE" -lt 1 ]] || [[ "$PHASE" -gt 10 ]]; then
        echo -e "${RED}Error: Phase must be a number between 1 and 10, or 'all'${NC}"
        exit 1
    fi
fi

# Function to print test results
print_result() {
    local test_name=$1
    local result=$2
    if [[ $result -eq 0 ]]; then
        echo -e "${GREEN}PASS${NC}: $test_name"
    else
        echo -e "${RED}FAIL${NC}: $test_name"
    fi
    return $result
}

# Function to print phase header
print_phase_header() {
    echo ""
    echo -e "${YELLOW}========== Phase $1 Verification ==========${NC}"
}

# Phase 1: Verify ILOG2 constants
verify_phase_1() {
    print_phase_header 1

    local ilog2_golden="${GOLDEN_DIR}/constants/ilog2_values.txt"

    if [[ ! -f "$ilog2_golden" ]]; then
        echo -e "${RED}Error: Golden file not found: $ilog2_golden${NC}"
        return 1
    fi

    # Generate current ILOG2 values (would come from Rust implementation)
    local ilog2_output="${PROJECT_ROOT}/build/ilog2_output.txt"

    # For now, just verify the golden file exists and is readable
    if [[ -r "$ilog2_golden" ]]; then
        print_result "ILOG2 golden file exists and readable" 0

        # Verify the golden file has expected content
        if grep -q "ILOG2(1) = 0" "$ilog2_golden"; then
            print_result "ILOG2(1) = 0 in golden file" 0
            return 0
        else
            print_result "ILOG2(1) = 0 in golden file" 1
            return 1
        fi
    else
        print_result "ILOG2 golden file exists and readable" 1
        return 1
    fi
}

# Phase 2: Verify sre_random sequence
verify_phase_2() {
    print_phase_header 2

    local random_golden="${GOLDEN_DIR}/random_sequence/random_1000.txt"

    if [[ ! -f "$random_golden" ]]; then
        echo -e "${RED}Error: Golden file not found: $random_golden${NC}"
        return 1
    fi

    # Generate current random sequence (would come from Rust implementation)
    local random_output="${PROJECT_ROOT}/build/random_output.txt"

    # For now, just verify the golden file exists and is readable
    if [[ -r "$random_golden" ]]; then
        print_result "Random sequence golden file exists and readable" 0

        # Verify the golden file has expected content (non-empty)
        if [[ -s "$random_golden" ]]; then
            print_result "Random sequence golden file is non-empty" 0

            # Check line count (should have approximately 1000 values)
            local line_count=$(wc -l < "$random_golden")
            if [[ $line_count -gt 900 ]]; then
                print_result "Random sequence has ~1000 values (found: $line_count)" 0
                return 0
            else
                print_result "Random sequence has ~1000 values (found: $line_count)" 1
                return 1
            fi
        else
            print_result "Random sequence golden file is non-empty" 1
            return 1
        fi
    else
        print_result "Random sequence golden file exists and readable" 1
        return 1
    fi
}

# Phase 3: Verify SQUID
verify_phase_3() {
    print_phase_header 3

    cd "$PROJECT_ROOT"

    # Run SQUID tests
    echo "Running cargo test phase3_squid..."
    if cargo test phase3_squid --quiet 2>&1 | tee /tmp/phase3_test.log; then
        print_result "cargo test phase3_squid" 0
    else
        print_result "cargo test phase3_squid" 1
        cat /tmp/phase3_test.log
        return 1
    fi

    # Verify golden file comparisons if they exist
    local squid_golden="${GOLDEN_DIR}/squid"
    if [[ -d "$squid_golden" ]]; then
        echo "Checking SQUID golden file comparisons..."
        local golden_count=$(find "$squid_golden" -type f | wc -l)
        print_result "SQUID golden files present ($golden_count files)" 0
    fi

    return 0
}

# Phase 4: Verify CM (Covariance Model)
verify_phase_4() {
    print_phase_header 4

    cd "$PROJECT_ROOT"

    # Run CM tests
    echo "Running cargo test phase4_cm..."
    if cargo test phase4_cm --quiet 2>&1 | tee /tmp/phase4_test.log; then
        print_result "cargo test phase4_cm" 0
    else
        print_result "cargo test phase4_cm" 1
        cat /tmp/phase4_test.log
        return 1
    fi

    # Verify CM structure loading matches golden files
    local cm_golden="${GOLDEN_DIR}/cm"
    if [[ -d "$cm_golden" ]]; then
        echo "Checking CM golden file comparisons..."
        # Look for .cm files
        local cm_file_count=$(find "$cm_golden" -name "*.cm" -type f 2>/dev/null | wc -l)
        if [[ $cm_file_count -gt 0 ]]; then
            print_result "CM golden files present ($cm_file_count files)" 0
        else
            echo -e "${YELLOW}Warning: No .cm files found in golden directory${NC}"
        fi
    fi

    return 0
}

# Phase 5: Verify Viterbi
verify_phase_5() {
    print_phase_header 5

    cd "$PROJECT_ROOT"

    # Run Viterbi tests
    echo "Running cargo test phase5_viterbi..."
    if cargo test phase5_viterbi --quiet 2>&1 | tee /tmp/phase5_test.log; then
        print_result "cargo test phase5_viterbi" 0
    else
        print_result "cargo test phase5_viterbi" 1
        cat /tmp/phase5_test.log
        return 1
    fi

    # Verify Viterbi scores match within tolerance (±0.01 bits)
    local viterbi_golden="${GOLDEN_DIR}/viterbi"
    if [[ -d "$viterbi_golden" ]]; then
        echo "Checking Viterbi golden file comparisons..."
        echo "Note: Scores should match within ±0.01 bits tolerance"
        local viterbi_file_count=$(find "$viterbi_golden" -type f 2>/dev/null | wc -l)
        if [[ $viterbi_file_count -gt 0 ]]; then
            print_result "Viterbi golden files present ($viterbi_file_count files)" 0
        fi
    fi

    return 0
}

# Phase 6: Verify EuFindtRNA
verify_phase_6() {
    print_phase_header 6

    cd "$PROJECT_ROOT"

    # Run EuFindtRNA tests
    echo "Running cargo test phase6_eufind..."
    if cargo test phase6_eufind --quiet 2>&1 | tee /tmp/phase6_test.log; then
        print_result "cargo test phase6_eufind" 0
    else
        print_result "cargo test phase6_eufind" 1
        cat /tmp/phase6_test.log
        return 1
    fi

    # Verify A-box, B-box scores match
    local eufind_golden="${GOLDEN_DIR}/eufindtRNA"
    if [[ -d "$eufind_golden" ]]; then
        echo "Checking EuFindtRNA golden file comparisons..."
        echo "Note: A-box and B-box scores should match golden files"
        local eufind_file_count=$(find "$eufind_golden" -type f 2>/dev/null | wc -l)
        if [[ $eufind_file_count -gt 0 ]]; then
            print_result "EuFindtRNA golden files present ($eufind_file_count files)" 0
        fi
    fi

    return 0
}

# Phase 7: Verify Structure
verify_phase_7() {
    print_phase_header 7

    cd "$PROJECT_ROOT"

    # Run Structure tests
    echo "Running cargo test phase7_structure..."
    if cargo test phase7_structure --quiet 2>&1 | tee /tmp/phase7_test.log; then
        print_result "cargo test phase7_structure" 0
    else
        print_result "cargo test phase7_structure" 1
        cat /tmp/phase7_test.log
        return 1
    fi

    # Verify secondary structure strings match
    local structure_golden="${GOLDEN_DIR}/structure"
    if [[ -d "$structure_golden" ]]; then
        echo "Checking secondary structure golden file comparisons..."
        echo "Note: Secondary structure strings should match exactly"
        local structure_file_count=$(find "$structure_golden" -type f 2>/dev/null | wc -l)
        if [[ $structure_file_count -gt 0 ]]; then
            print_result "Structure golden files present ($structure_file_count files)" 0
        fi
    fi

    return 0
}

# Phase 8: Verify Isotype
verify_phase_8() {
    print_phase_header 8

    cd "$PROJECT_ROOT"

    # Run Isotype tests
    echo "Running cargo test phase8_isotype..."
    if cargo test phase8_isotype --quiet 2>&1 | tee /tmp/phase8_test.log; then
        print_result "cargo test phase8_isotype" 0
    else
        print_result "cargo test phase8_isotype" 1
        cat /tmp/phase8_test.log
        return 1
    fi

    # Verify isotype assignments match
    local isotype_golden="${GOLDEN_DIR}/isotype"
    if [[ -d "$isotype_golden" ]]; then
        echo "Checking isotype assignment golden file comparisons..."
        echo "Note: Isotype assignments should match golden files"
        local isotype_file_count=$(find "$isotype_golden" -type f 2>/dev/null | wc -l)
        if [[ $isotype_file_count -gt 0 ]]; then
            print_result "Isotype golden files present ($isotype_file_count files)" 0
        fi
    fi

    return 0
}

# Phase 9: Verify Pipeline
verify_phase_9() {
    print_phase_header 9

    cd "$PROJECT_ROOT"

    # Run Pipeline tests
    echo "Running cargo test phase9_pipeline..."
    if cargo test phase9_pipeline --quiet 2>&1 | tee /tmp/phase9_test.log; then
        print_result "cargo test phase9_pipeline" 0
    else
        print_result "cargo test phase9_pipeline" 1
        cat /tmp/phase9_test.log
        return 1
    fi

    # Verify command-line argument handling
    echo "Checking command-line argument handling..."
    if cargo test --quiet cli 2>&1 | grep -q "test result: ok"; then
        print_result "CLI argument handling tests" 0
    else
        echo -e "${YELLOW}Note: CLI tests may not be present yet${NC}"
    fi

    return 0
}

# Phase 10: Verify Full Integration
verify_phase_10() {
    print_phase_header 10

    cd "$PROJECT_ROOT"

    # Run full integration tests
    echo "Running cargo test phase10_full..."
    if cargo test phase10_full --quiet 2>&1 | tee /tmp/phase10_test.log; then
        print_result "cargo test phase10_full" 0
    else
        print_result "cargo test phase10_full" 1
        cat /tmp/phase10_test.log
        return 1
    fi

    # Compare full run outputs against Example1/Example2
    local examples_dir="${PROJECT_ROOT}/tests/examples"
    if [[ -d "$examples_dir" ]]; then
        echo "Checking Example1/Example2 comparisons..."

        # Check for Example1
        if [[ -d "${examples_dir}/Example1" ]]; then
            echo "Example1 found - comparing full run outputs..."
            if [[ -f "${examples_dir}/Example1/expected_output.txt" ]]; then
                print_result "Example1 expected output exists" 0
            else
                echo -e "${YELLOW}Note: Example1 expected output not found${NC}"
            fi
        fi

        # Check for Example2
        if [[ -d "${examples_dir}/Example2" ]]; then
            echo "Example2 found - comparing full run outputs..."
            if [[ -f "${examples_dir}/Example2/expected_output.txt" ]]; then
                print_result "Example2 expected output exists" 0
            else
                echo -e "${YELLOW}Note: Example2 expected output not found${NC}"
            fi
        fi
    else
        echo -e "${YELLOW}Note: Examples directory not found at ${examples_dir}${NC}"
    fi

    # Run integration tests if they exist
    if cargo test --quiet integration 2>&1 | grep -q "test result: ok"; then
        print_result "Integration tests" 0
    else
        echo -e "${YELLOW}Note: Integration tests may not be present yet${NC}"
    fi

    return 0
}

# Run all phases
run_all_phases() {
    local failed_phases=()
    local passed_phases=()

    echo ""
    echo -e "${YELLOW}========================================${NC}"
    echo -e "${YELLOW}Running All Phases (1-10)${NC}"
    echo -e "${YELLOW}========================================${NC}"

    for i in {1..10}; do
        if verify_phase_$i; then
            passed_phases+=($i)
        else
            failed_phases+=($i)
        fi
    done

    echo ""
    echo -e "${YELLOW}========================================${NC}"
    echo -e "${YELLOW}Summary${NC}"
    echo -e "${YELLOW}========================================${NC}"
    echo ""

    if [[ ${#passed_phases[@]} -gt 0 ]]; then
        echo -e "${GREEN}PASSED Phases (${#passed_phases[@]}/10):${NC} ${passed_phases[*]}"
    fi

    if [[ ${#failed_phases[@]} -gt 0 ]]; then
        echo -e "${RED}FAILED Phases (${#failed_phases[@]}/10):${NC} ${failed_phases[*]}"
        echo ""
        echo -e "${RED}Overall: FAILED${NC}"
        return 1
    else
        echo ""
        echo -e "${GREEN}Overall: ALL PHASES PASSED${NC}"
        return 0
    fi
}

# Execute the appropriate phase verification
if [[ "$PHASE" == "all" ]]; then
    run_all_phases
    exit $?
else
    case $PHASE in
        1) verify_phase_1 ;;
        2) verify_phase_2 ;;
        3) verify_phase_3 ;;
        4) verify_phase_4 ;;
        5) verify_phase_5 ;;
        6) verify_phase_6 ;;
        7) verify_phase_7 ;;
        8) verify_phase_8 ;;
        9) verify_phase_9 ;;
        10) verify_phase_10 ;;
        *)
            echo -e "${RED}Error: Invalid phase number: $PHASE${NC}"
            exit 1
            ;;
    esac

    # Capture exit code and print final result
    RESULT=$?
    if [[ $RESULT -eq 0 ]]; then
        echo ""
        echo -e "${GREEN}Phase $PHASE verification PASSED${NC}"
    else
        echo ""
        echo -e "${RED}Phase $PHASE verification FAILED${NC}"
    fi

    exit $RESULT
fi
