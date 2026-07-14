# trnascan-rs

A faithful, byte-parity Rust reimplementation of **tRNAscan-SE 2.0** (Chan, Lin &
Lowe) for bacterial (`-B`) and archaeal (`-A`) genomes. The covariance-model
search stage is provided in-process by the
[`infernox`](https://crates.io/crates/infernox) library (a byte-parity Rust
port of Infernal 1.1.x `cmsearch`/`cmscan`), pulled from crates.io (`>= 0.1.4`)
— there is **no subprocess** and no dependency on an external Infernal install.

The port reproduces the C/Perl reference pipeline mechanism-for-mechanism:
first-pass prescan → per-candidate Phase-II verification → decode → NS rescore →
isotype cmscan → truncated-tRNA search → archaeal BHB noncanonical-intron scan →
output formatting.

## Parity status

Output is **byte-identical** to the local reference C/Perl tRNAscan-SE
(Infernal 1.1.5) across:

- **GTDB 38-genome sweep** (one representative per bacterial/archaeal phylum):
  `OK = 38 / 38`, key-level `Conly = 0, Ronly = 0`.
- `.out`, `--detail` struct, and `-a` fasta formats — verified with direct
  full-contig `diff` on *E. coli*, Micrarchaeota and Thermoproteota.
- Unit + integration suite: `260` lib tests + all phase tests pass.

**Sole known numerical divergence:** one Thermoproteota Cys score cell renders
`76.5` vs C's `76.6`. This is an infernox Inside-DP floating-point rounding
difference (the true value `76.5492` straddles the `.55` display-rounding
boundary), not a port-logic defect; it appears identically in every output
format.

**Self-identification (by design):** trnascan-rs reports *its own* name and
version in the fields that name the producing program — the `.stats` (`-m`)
banner and header, the GFF3 (`--gff`) source column, and the ACeDB (`--acedb`)
program ID now read `trnascan-rs` / `v0.2.0` instead of C's `tRNAscan-SE` /
`2.0.12`. These are the only deliberate string differences; every biological
result field (coordinates, scores, structures, isotypes, anticodons) — in
`.out`, `--detail`, `-a`, GFF, and ACeDB — is unchanged.

The parity target is the **locally-run** reference (Infernal 1.1.5), *not* the
shipped `Demo/*.out` goldens (produced with Infernal 1.1.2 and not reproducible
here). See `docs/faithful_port_spec.md`.

## Layout

This package is a **crates.io source crate** (`cargo publish` / `cargo install`
basis). The covariance models are *not* bundled — supply them at runtime with
`--models-dir` (they ship with a tRNAscan-SE 2.0 install, under `lib/models`).

```
trnascan-rs/
├── src/                Rust source
├── Cargo.toml          (infernox from crates.io, >= 0.1.4)
├── Cargo.lock          pinned dependency graph
└── LICENSE, COPYING    GPLv3 + upstream attribution
```

## Installation

```sh
cargo install trnascan-rs        # from crates.io
# or, from this source tree:
cargo build --release            # binary at target/release/trnascan-rs
```

## Usage

The model directory must be supplied with `--models-dir` (the built-in default
is the cwd-relative `models/`); point it at a tRNAscan-SE 2.0 `lib/models`:

```sh
MODELS=/path/to/tRNAscan-SE/lib/models

# bacterial
trnascan-rs -B --models-dir "$MODELS" -o out.txt genome.fna

# archaeal
trnascan-rs -A --models-dir "$MODELS" -o out.txt genome.fna

# detailed (-H = HMM/2'str score cols, --detail = isotype cols)
trnascan-rs -B -H --detail --models-dir "$MODELS" -o out.txt genome.fna

# fasta of predicted tRNA precursors
trnascan-rs -B --models-dir "$MODELS" -a trnas.fa -o out.txt genome.fna
```

## Building from source

Requires a Rust toolchain (1.70+); `infernox` resolves from crates.io, so no
external Infernal/Easel install and no C toolchain are needed. From the package
root:

```sh
cargo build --release
# binary at target/release/trnascan-rs
```

To build against a local `infernox` checkout instead of the published crate, add
`[patch.crates-io]` `infernox = { path = "../infernox" }` to `Cargo.toml`.

## Changelog

### 0.2.0
- Faithful, byte-parity Rust port of tRNAscan-SE 2.0 for `-B`/`-A` genomes,
  with the covariance-model stage provided in-process by `infernox` (Infernal
  1.1.5 parity). Byte-identical `.out`, `--detail`, and `-a` output verified on
  a GTDB 38-genome bacterial/archaeal sweep.
- **Self-identification.** The tool reports its own name/version — `trnascan-rs`
  `v0.2.0` — in the `--version`/credits, the `.stats` (`-m`) banner and header,
  the GFF3 (`--gff`) source column, and the ACeDB (`--acedb`) program ID, rather
  than impersonating C's `tRNAscan-SE 2.0.12`. Biological result fields are
  unchanged; the scan banner keeps C's layout.
- **Argv-faithful drop-in CLI.** Every flag name matches C tRNAscan-SE exactly
  — no Rust-only aliases: mode selectors (`-B`/`-A`/`-E`/`-O`/`-G`/`-M`/`-L`/
  `-I`/`-S`/`-U` with their long forms), `-D` (`--nopseudo`), `--nomerge`,
  `-b`/`-a`/`-f`/`-s`/`-j`/`-o`, `--brief`, `-X` (score cutoff, compared on the
  rounded display value as C does), `-z` (`--pad`), `--len`, `--codons`,
  `--isocm`, `--acedb`, `-w` (odd-struct file), `-Y` (PID file).
- `--acedb` output matches C `save_Acedb_from_secpass` (pre-promotion isotype
  in the `tRNAscan_id`, raw ascending coordinates, `SelCys→Z`/`fMet→∅` one-letter
  map, pseudo `Remark`). `-X` boundary hits at exactly the cutoff are kept
  (C compares the 1-decimal-rounded score).
- Requires `infernox >= 0.1.4` (the `--mid`-pass `--Fmid` threshold fix that
  keeps the archaeal isotype scan — e.g. Hsalinarum tRNA-Gln — byte-identical
  to C).
- Removed the unwired legacy pre-parity code that the faithful `core::scanner`
  path superseded: the tRNAscan-1.4 Fichant-Burks first-pass scaffolding, and
  the old `pipeline`/`scan_result`/`isotype::scorer` modules (the latter held a
  placeholder `score_cm_model` returning a `-999.0` dummy — never on the live
  path; real per-isotype CM scoring runs in-process through infernox cmscan).
  No behavioral change; `-B`/`-A` output stays byte-identical to C.

## License

GPLv3 (see `LICENSE`). tRNAscan-SE is Copyright © Patricia P. Chan, Brian Lin,
and Todd M. Lowe (see `COPYING`); this is a derivative Rust port and inherits
the GPLv3 terms.
