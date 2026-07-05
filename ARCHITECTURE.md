# tRNAscan-SE 2.0 — Working Architecture (for the Rust re-development)

Source of truth: the Perl implementation under `original/` (driver `tRNAscan-SE.src`
1966 lines + `lib/tRNAscanSE/*.pm` ~14.6k lines). This maps the ACTUAL control
flow / algorithm, with `file:line` anchors, so the Rust port can reproduce it.

Analyzed 2026-07-02.

---

## 0. TL;DR — two-phase pipeline

```
main (tRNAscan-SE.src)
  set_options()                         # CLI → mode preset + engine flags
  initialize_process()
  ── PHASE I  first_pass_prescan()  (line 113 / def 185)   # fast candidate finding
  ── PHASE II run_cm_scan()         (line 139 / def 341)   # CM confirmation + post-proc + output
  cleanup()
```
- **Phase I** finds candidate tRNA regions cheaply (tRNAscan1.4 / EufindtRNA / Infernal-HMM).
- **Phase II** confirms & scores each candidate with a covariance model (Cove *or* Infernal
  cmsearch), then does intron / truncation / isotype / pseudogene analysis and writes output.
- The **mode flag (-E/-B/-A/-G/-O/-M/…) is a PRESET** that selects which scanners, which CM
  models, and which post-processors run.

Our byte-identical Rust cmsearch (`infernal/rust` `faithful_search`) is the engine behind
Phase II `analyze_with_cmsearch` (and behind the isotype / intron / truncation cmsearch calls).

---

## 1. Modes = presets (tRNAscan-SE.src ~960-1330, Options.pm)

| Flag | Mode | First-pass | Second-pass | Isotype | Notes |
|---|---|---|---|---|---|
| -E (default) | Eukaryotic | tRNAscan + EufindtRNA | Infernal cmsearch | yes | full support |
| -B | Bacterial | tRNAscan + EufindtRNA | Infernal cmsearch | yes | tuned eufind_intscore |
| -A | Archaeal | tRNAscan + EufindtRNA | Infernal cmsearch | yes | +noncanonical introns, +fragments |
| -G | General (3-domain) | tRNAscan + EufindtRNA | Infernal cmsearch | no | all domains, no isotype |
| -O | Organellar | none | Infernal/Cove | no | custom cutoff, no first-pass |
| -M | Mitochondrial | none (Infernal FP) | `analyze_mito` | — | mito models, alt genetic code |
| -T metagenome / -N NUMT / -U alternate | | Infernal-HMM FP | Infernal cmsearch | no | Infernal both passes |
| -L legacy | | tRNAscan(+Eufind) | **Cove** | | Cove-only 2nd pass |

Engine selection flags: `-C` → Cove (`CM_mode="cove"`), `-I` → Infernal FP,
default → Infernal cmsearch. `second_pass_label` = "Cove"/"Infernal".
Padding around first-pass hits before 2nd pass: default **10 bp** (Options `default_Padding`).

---

## 2. PHASE I — first_pass_prescan (tRNAscan-SE.src 185-338)

Reads FASTA in buffered chunks; per sequence, dispatch:

| Cond | Call | Engine / binary |
|---|---|---|
| `infernal_fp()` | `$cm->first_pass_scan()` (231) → `run_first_pass_cmsearch` (CM.pm 3105) | cmsearch scan_flag=3 `-g --mid --notrunc -T fp_cutoff` |
| `tscan_mode()` | `$tscan->run_tRNAscan()` (238) + `process_tRNAscan_hits()` (243) | `trnascan-1.4 -i <idx> -c <params>` (Tscan.pm 148) |
| `eufind_mode()` | `$eufind->run_eufind()` (249) + `process_Eufind_hits()` (252) | `eufindtRNA -i -F -I <intscore> -l <maxintron>` (Eufind.pm 110) |

- Parsing: Tscan regex on raw file (`start/end position=`, `tRNA predict as a tRNA-<T>: anticodon <AC>`, `potential intron between positions` → CI intron); Eufind tab-split 9 fields.
- Storage: `ArraytRNA $fp_tRNAs`; `merge_repeat_hit()` dedups overlaps (union bounds, OR `hit_source` mask: Ts=1, Eu=2, Both=3).
- Persist: `FpScanResultFile::save_firstpass_output()` (tab file). `get_next_tRNA_candidate()` applies ±padding (plus/minus-strand aware) when handing to Phase II.

---

## 3. PHASE II — run_cm_scan (tRNAscan-SE.src 341-499)

Per first-pass-indexed sequence:
```
cove_mode:  loop each candidate → prepare_tRNA_to_scan → analyze_with_cove (394)
else:       prepare_multi_tRNAs_to_scan
            mito_mode      → analyze_mito           (413)
            alternate_mode → analyze_alternate      (417)
            infernal_mode  → analyze_with_cmsearch  (429)   # DEFAULT
── if introns enabled:            scan_noncanonical_introns (441)
── if euk/bact/arch:              truncated_tRNA_search (449)
                                  isotype_cmsearch (453)  unless --no-isotype
── output_tRNA (457)              # ScanResult.pm
── if split halves:               scan_split_tRNAs (473)
── BED / GFF output (477-491)
```

### 3a. Infernal engine calls (CM.pm) — `exec_cmsearch` (2530) / `exec_cmscan` (2490)
cmd: `cmsearch <cm_options> <cm_file> <seq> > out`; cmscan: `... --fmt 2 --tblout <tab> -o <out>`.
`cm_options` by **scan_flag**:
| flag | options | use |
|---|---|---|
|0|`-g --nohmm --toponly --notrunc`|standard|
|1|`-g --mid --notrunc`|HMM filter|
|2|`-g --mid --toponly --notrunc`|isotype / tscan-eufind|
|3|`-g --mid --notrunc -T fp_cutoff`|first-pass|
|5|`-g --max --toponly --notrunc --notextw -T BHB`|intron (split halves)|
|6|`-g --toponly --notextw`|truncation check|
|7|`-g --max --toponly --notrunc --notextw -T 0`|rescore|
(+ conditional `-T score_cutoff` when ≤10 for flags ≠3,5,7.)

- `analyze_with_cmsearch` (3584): run each domain model → merge (`ArrayCMscanResults`) → filter by `cm_cutoff` → adjust coords/strand → **CCA-suffix adjustment** → `rescore_tRNA` if modified → `decode_tRNA_properties`.
- Output parsing: `CMscanResultFile` / `ArrayCMscanResults` parse `--fmt 2 --tblout` into tRNA records (seqname, score, start, end, strand, ss, seq, model).

### 3b. Isotype / anticodon (CM.pm)
- `isotype_cmsearch` (2761) → `scan_isotype_cm` (2797): run isotype-specific CM DB (cyto + mito) via `exec_cmscan` scan_flag=2 → per-amino-acid scores → best model = isotype.
- `find_anticodon` (731): parse the anticodon stem-loop from the secondary-structure string; extract middle 3 bases → anticodon. `decode_tRNA_properties` (1258) maps anticodon→type via `GeneticCode::get_tRNA_type`, handles Met/iMet/fMet/Ile2 disambiguation.

### 3c. Introns
- Canonical (CI): from tRNAscan first-pass output + `find_intron` (960, lowercase run in anticodon loop, min length).
- Non-canonical (NCI): `scan_noncanonical_introns` (1598) — 2 rounds of `run_cmsearch_intron` (scan_flag=5) + `check_intron_validity` (1885, SS regex splits pre/intron/post).

### 3d. Truncation & split
- `truncated_tRNA_search` (scan_flag=6) — 5'/3' truncated tRNAs (arch/euk/bact).
- `scan_split_tRNAs` — split half tRNAs (archaea).

---

## 4. Scoring model (tRNA.pm)
Per-tRNA scores decomposed:
- `score` (`_score`): total covariance bit score (Cove or Infernal).
- `hmm_score`: primary-structure-only score, from a **no-secondary-structure model** (`mainNS_cm`) run (`cmsearch_scoring` / coves).
- `ss_score` = `score − hmm_score` (secondary-structure contribution).
- Stored per engine in `_h_domain_models{cove|infernal}`; `set_default_scores` (742) picks final.

### Pseudogene filter (`is_pseudo_gene`, CM.pm 994)
Flag pseudo if `(ss_score < min_ss_score  OR  hmm_score < min_hmm_score)  AND  score < min_pseudo_filter_score`
(cove≈40 / infernal≈50; min_hmm≈20, min_ss≈5). Euk/Mito HighConfidenceFilter binaries add stricter post-filters.

---

## 5. Output (ScanResult.pm)
- Default `.out` (`construct_tab_output` 639): SeqName, tRNA#, Begin, End, Type, Codon/AC, Intron Begin/End, Inf Score, HMM, 2'Str, (Inf), Hit-origin, Isotype-CM, Isotype-score, Note.
- `-f .ss` secondary structure (`save_allStruct_output`), `-b` BED12 (introns as blocks), GFF3, `-m` stats, `.iso` per-isotype score table.
- `set_mature_tRNA` (tRNA.pm 797): splice out introns → mature seq/ss.

---

## 6. Implications for the Rust re-development
1. **Engine**: replace the current Rust port's shell-out to external `cmsearch` (and the legacy
   `infernal::cmsearch` stub) with our byte-identical `faithful_search` pipeline, exposed as a
   **library entry** callable per (model, seq, options=scan_flag). Needs to support the scan_flag
   option sets above (esp. `-g --nohmm --toponly --notrunc` default; `--mid`; `--max`; `-T`).
2. **Multi-model orchestration**: Phase II runs the engine once per domain model, merges, filters
   by cutoff; isotype runs it again over the isotype CM DB; introns/truncation are more engine runs.
3. **Post-processing is substantial** and currently broken/incomplete in the Rust port: CCA
   adjustment, anticodon extraction from SS, isotype assignment (incl. Met disambiguation), CI/NCI
   introns, pseudogene filter, score decomposition (needs the no-SS model run for hmm_score),
   mature-seq splicing, and all output formats.
4. **Models**: use the Infernal-format (`INFERNAL1/a`) models in `original/lib/models/TRNAinf*.cm`
   (domain + isotype + mito) — directly readable by our engine. `TRNA2*.cm` are legacy Cove format
   (different engine; only needed for `-C`/legacy mode).
5. **Verification**: golden outputs in `original/Demo/Example{1,2}-tRNAs.{out,ss,bed,iso,stats}`.
