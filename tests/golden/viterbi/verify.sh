#!/bin/bash
# verify.sh - Quick verification script for Viterbi golden files

set -e

echo "=== Viterbi Golden Files Verification ==="
echo

# Check all required files exist
echo "Checking file existence..."
files=(
    "gen_viterbi.c"
    "Makefile"
    "viterbi_scores.txt"
    "traceback.txt"
    "test_sequences.fa"
    "README.md"
    "SUMMARY.md"
)

for file in "${files[@]}"; do
    if [ -f "$file" ]; then
        echo "✓ $file"
    else
        echo "✗ $file - MISSING!"
        exit 1
    fi
done

echo
echo "Checking golden file contents..."

# Verify viterbi_scores.txt has correct number of test cases
score_count=$(grep -c "^[A-Za-z]" viterbi_scores.txt)
if [ "$score_count" -eq 3 ]; then
    echo "✓ viterbi_scores.txt contains 3 test cases"
else
    echo "✗ viterbi_scores.txt should have 3 test cases, found $score_count"
    exit 1
fi

# Verify expected scores are present
if grep -q "73.8820" viterbi_scores.txt; then
    echo "✓ Phe-73bp score present (73.8820)"
else
    echo "✗ Phe-73bp score missing"
    exit 1
fi

if grep -q "66.0850" viterbi_scores.txt; then
    echo "✓ Ser-82bp score present (66.0850)"
else
    echo "✗ Ser-82bp score missing"
    exit 1
fi

if grep -q -- "-12.8010" viterbi_scores.txt; then
    echo "✓ Fragment-50bp score present (-12.8010)"
else
    echo "✗ Fragment-50bp score missing"
    exit 1
fi

# Verify traceback.txt has all three test cases
echo
echo "Checking traceback file..."
traceback_count=$(grep -c "^===" traceback.txt)
if [ "$traceback_count" -eq 3 ]; then
    echo "✓ traceback.txt contains 3 traceback trees"
else
    echo "✗ traceback.txt should have 3 trees, found $traceback_count"
    exit 1
fi

# Check that traceback has proper structure
if grep -q "node=0 type=64.*uBEGIN_ST" traceback.txt; then
    echo "✓ Traceback trees start with uBEGIN_ST"
else
    echo "✗ Traceback trees missing uBEGIN_ST root"
    exit 1
fi

if grep -q "type=128.*uEND_ST" traceback.txt; then
    echo "✓ Traceback trees contain uEND_ST terminals"
else
    echo "✗ Traceback trees missing uEND_ST"
    exit 1
fi

# Verify test_sequences.fa
echo
echo "Checking test sequences..."
seq_count=$(grep -c "^>" test_sequences.fa)
if [ "$seq_count" -eq 3 ]; then
    echo "✓ test_sequences.fa contains 3 sequences"
else
    echo "✗ test_sequences.fa should have 3 sequences, found $seq_count"
    exit 1
fi

# Check executable exists
echo
echo "Checking compiled binary..."
if [ -x "gen_viterbi" ]; then
    echo "✓ gen_viterbi executable present"
else
    echo "⚠ gen_viterbi not compiled (run 'make gen_viterbi')"
fi

echo
echo "=== All checks passed! ==="
echo
echo "Summary:"
echo "  - 7 required files present"
echo "  - 3 test cases with valid scores"
echo "  - 3 complete traceback trees"
echo "  - All structural checks passed"
echo
echo "Golden files are ready for Rust implementation testing."
