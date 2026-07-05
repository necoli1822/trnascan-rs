#!/bin/bash
# Verification script for EuFindtRNA golden files

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "=== EuFindtRNA Golden Files Verification ==="
echo

# Check all required files exist
echo "Checking file existence..."
required_files=(
    "bbox_scores.txt"
    "abox_scores.txt"
    "trna_detection.txt"
    "weight_matrices.txt"
    "gen_eufind.c"
    "gen_weight_matrices.c"
    "Makefile"
    "README.md"
    "SUMMARY.md"
)

missing=0
for file in "${required_files[@]}"; do
    if [ -f "$file" ]; then
        size=$(stat -f%z "$file" 2>/dev/null || stat -c%s "$file" 2>/dev/null)
        echo "  ✓ $file ($size bytes)"
    else
        echo "  ✗ $file (MISSING)"
        missing=$((missing + 1))
    fi
done

echo

if [ $missing -gt 0 ]; then
    echo "ERROR: $missing file(s) missing!"
    exit 1
fi

# Verify key values in golden files
echo "Verifying key values..."

# Check B-box detection
if grep -q "116.*-1.9170.*FOUND" bbox_scores.txt; then
    echo "  ✓ B-box position 116 with score -1.9170"
else
    echo "  ✗ B-box detection values incorrect"
    exit 1
fi

# Check A-box detection
if grep -q "A-box score:.*-13.7640" abox_scores.txt; then
    echo "  ✓ A-box score -13.7640"
else
    echo "  ✗ A-box score incorrect"
    exit 1
fi

if grep -q "AB dist score:.*-5.4420" abox_scores.txt; then
    echo "  ✓ AB distance score -5.4420"
else
    echo "  ✗ AB distance score incorrect"
    exit 1
fi

# Check total detection
if grep -q "Total score:.*-21.67" trna_detection.txt; then
    echo "  ✓ Total score -21.67"
else
    echo "  ✗ Total score incorrect"
    exit 1
fi

# Check weight matrices
if grep -q "Abox_Mat\[0\]\[0\] = -1.268" weight_matrices.txt; then
    echo "  ✓ Abox_Mat[0][0] = -1.268"
else
    echo "  ✗ Weight matrix values incorrect"
    exit 1
fi

# Count matrix entries
abox_entries=$(grep "^Abox_Mat" weight_matrices.txt | wc -l)
bbox_entries=$(grep "^Bbox_Mat" weight_matrices.txt | wc -l)

if [ "$abox_entries" -eq 126 ]; then  # 6 rows × 21 positions
    echo "  ✓ Abox_Mat has 126 entries (6×21)"
else
    echo "  ✗ Abox_Mat has $abox_entries entries (expected 126)"
    exit 1
fi

if [ "$bbox_entries" -eq 66 ]; then  # 6 rows × 11 positions
    echo "  ✓ Bbox_Mat has 66 entries (6×11)"
else
    echo "  ✗ Bbox_Mat has $bbox_entries entries (expected 66)"
    exit 1
fi

echo
echo "=== All Verifications Passed ==="
echo
echo "Golden files are ready for use in Rust implementation testing."
