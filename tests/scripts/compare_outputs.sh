#!/bin/bash

# Script to compare tRNA detection outputs
# Usage: compare_outputs.sh expected_file actual_file
# Exits with 0 if outputs are equivalent, 1 if different

set -e

# Check arguments
if [[ $# -ne 2 ]]; then
    echo "Usage: $0 expected_file actual_file" >&2
    exit 1
fi

EXPECTED_FILE="$1"
ACTUAL_FILE="$2"

# Check if files exist
if [[ ! -f "$EXPECTED_FILE" ]]; then
    echo "Error: Expected file not found: $EXPECTED_FILE" >&2
    exit 1
fi

if [[ ! -f "$ACTUAL_FILE" ]]; then
    echo "Error: Actual file not found: $ACTUAL_FILE" >&2
    exit 1
fi

# Create temporary files with headers stripped
EXPECTED_STRIPPED=$(mktemp)
ACTUAL_STRIPPED=$(mktemp)

# Clean up temp files on exit
trap "rm -f $EXPECTED_STRIPPED $ACTUAL_STRIPPED" EXIT

# Strip header lines (lines starting with "Sequence", "Name", or "---")
grep -v "^Sequence\|^Name\|^---" "$EXPECTED_FILE" > "$EXPECTED_STRIPPED" || true
grep -v "^Sequence\|^Name\|^---" "$ACTUAL_FILE" > "$ACTUAL_STRIPPED" || true

# Function to normalize and compare lines with floating-point tolerance
compare_lines() {
    local line_num=0
    local has_diff=0

    # Compare using awk for numeric field comparison with tolerance
    awk -v tolerance=0.1 '
    BEGIN {
        diff_count = 0
    }
    NR == FNR {
        expected[NR] = $0
        next
    }
    {
        actual_line = $0
        expected_line = expected[NR]

        if (expected_line == "") {
            print "Line " NR ": Extra line in actual output: " actual_line
            diff_count++
            next
        }

        if (expected_line != actual_line) {
            # Try numeric comparison for scores
            # Count fields
            n_expected = split(expected_line, expected_fields)
            n_actual = split(actual_line, actual_fields)

            if (n_expected == n_actual) {
                numeric_diff = 0
                for (i = 1; i <= n_expected; i++) {
                    if (expected_fields[i] != actual_fields[i]) {
                        # Check if both are numeric
                        if (expected_fields[i] ~ /^[0-9.-]+$/ && actual_fields[i] ~ /^[0-9.-]+$/) {
                            diff = (expected_fields[i] - actual_fields[i])
                            if (diff < 0) diff = -diff
                            if (diff > tolerance) {
                                numeric_diff = 1
                                break
                            }
                        } else {
                            numeric_diff = 1
                            break
                        }
                    }
                }
                if (numeric_diff) {
                    print "Line " NR ": Mismatch"
                    print "  Expected: " expected_line
                    print "  Actual:   " actual_line
                    diff_count++
                }
            } else {
                print "Line " NR ": Different number of fields"
                print "  Expected: " expected_line
                print "  Actual:   " actual_line
                diff_count++
            }
        }
    }
    END {
        # Check for missing lines in actual output
        for (i = NR + 1; i <= length(expected); i++) {
            if (expected[i] != "") {
                print "Line " i ": Missing from actual output: " expected[i]
                diff_count++
            }
        }

        if (diff_count > 0) {
            print "Total differences: " diff_count
            exit 1
        }
        exit 0
    }
    ' "$EXPECTED_STRIPPED" "$ACTUAL_STRIPPED"

    return $?
}

# Run comparison
if compare_lines; then
    echo "OK: Outputs are equivalent"
    exit 0
else
    echo "FAIL: Outputs differ"
    exit 1
fi
