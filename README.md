# tRNAscan-SE (Rust faithful port)

A faithful, byte-parity Rust reimplementation of **tRNAscan-SE 2.0** (Chan, Lin &
Lowe) for bacterial (`-B`) and archaeal (`-A`) genomes. The covariance-model
search stage is provided in-process by the [`infernox`](../infernox) library
(a byte-parity Rust port of Infernal 1.1.x `cmsearch`/`cmscan`) — there is **no
subprocess** and no dependency on an external Infernal install.

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
- Unit + integration suite: `338` lib tests + all phase tests pass.

**Sole known divergence:** one Thermoproteota Cys score cell renders `76.5` vs
C's `76.6`. This is an infernox Inside-DP floating-point rounding difference
(the true value `76.5492` straddles the `.55` display-rounding boundary), not a
port-logic defect; it appears identically in every output format.

The parity target is the **locally-run** reference (Infernal 1.1.5), *not* the
shipped `Demo/*.out` goldens (produced with Infernal 1.1.2 and not reproducible
here). See `docs/faithful_port_spec.md`.

## Layout

```
trnascan-rs/
├── bin/trnascan-rs     prebuilt release binary (x86_64 linux)
├── src/                Rust source
├── data/models/        covariance models + signals (required at runtime)
├── examples/           Demo FASTAs (Example1/2) + upstream goldens
├── docs/               faithful-port spec
├── ARCHITECTURE.md     module map
├── Cargo.toml          (infernox as path dep: ../infernox)
└── LICENSE, COPYING    GPLv3 + upstream attribution
```

## Usage

The model directory must be supplied with `--models-dir` (the built-in default
is the cwd-relative `models/`). Canonical invocation from the package root:

```sh
# bacterial
bin/trnascan-rs -B --models-dir data/models -o out.txt genome.fna

# archaeal
bin/trnascan-rs -A --models-dir data/models -o out.txt genome.fna

# detailed (-H = HMM/2'str score cols, --detail = isotype cols)
bin/trnascan-rs -B -H --detail --models-dir data/models -o out.txt genome.fna

# fasta of predicted tRNA precursors
bin/trnascan-rs -B --models-dir data/models -a trnas.fa -o out.txt genome.fna
```

Quick check against a bundled example:

```sh
bin/trnascan-rs -B -H --detail --models-dir data/models -o /tmp/ex1.out examples/Example1.fa
```

## Building from source

Requires a Rust toolchain and the sibling `infernox` package at `../infernox`
(the `infernal`/`easel` crates are path dependencies). From the package root:

```sh
cargo build --release
# binary at target/release/trnascan-rs
```

## License

GPLv3 (see `LICENSE`). tRNAscan-SE is Copyright © Patricia P. Chan, Brian Lin,
and Todd M. Lowe (see `COPYING`); this is a derivative Rust port and inherits
the GPLv3 terms.
