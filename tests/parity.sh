#!/bin/bash
# tRNAscan-rs faithful-port parity harness.
#
# Parity target = the C/Perl reference tRNAscan-SE run LOCALLY (infernal 1.1.5),
# NOT the shipped Demo/*.out (those were made with infernal 1.1.2 and are not
# reproducible here — see docs/faithful_port_spec.md "PARITY TARGET CORRECTION").
#
# Usage: bash tests/parity.sh [regen]
#   regen : (re)generate the reference goldens from the Perl driver first.
set -uo pipefail
ROOT=/mnt/DAS/sunju/programme/bactars/tRNAscan-SE
ORIG=$ROOT/original
GOLD=$ROOT/tests/reference_golden
RS=$ROOT/target/release/trnascan-rs
export PERL5LIB=$ORIG/lib

regen() {
  mkdir -p "$GOLD"
  for ex in Example1 Example2; do
    perl "$ORIG/tRNAscan-SE" -B -H --detail -c "$ORIG/tRNAscan-SE.conf" \
      -o "$GOLD/$ex-Bhd.out" "$ORIG/Demo/$ex.fa" 2>/dev/null
  done
  echo "regenerated reference goldens in $GOLD"
}

[ "${1:-}" = "regen" ] && regen

# Compare rust output vs reference golden, column-aware.
run_one() {
  local ex=$1
  local rout; rout=$(mktemp)
  # NOTE: the rust CLI flag surface for -H/--detail may still be WIP; adjust here
  # as the port matures. For now request the same 15-column detailed -B output.
  "$RS" -B -H --detail --models-dir "$ORIG/lib/models" -o "$rout" "$ORIG/Demo/$ex.fa" 2>/dev/null
  echo "===== $ex ====="
  if diff -u "$GOLD/$ex-Bhd.out" "$rout" > /dev/null; then
    echo "BYTE-IDENTICAL ✓"
  else
    echo "DIFF (rust vs reference):"
    diff -u "$GOLD/$ex-Bhd.out" "$rout" | head -40
  fi
  rm -f "$rout"
}

[ -x "$RS" ] || { echo "build first: (cd $ROOT && cargo build --release)"; exit 1; }
run_one Example1
run_one Example2
