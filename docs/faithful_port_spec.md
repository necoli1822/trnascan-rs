# tRNAscan-SE `-B` Faithful Rust Port — Implementation Spec

**Goal:** byte-parity with the C/Perl reference `-B` output for bacterial mode.

> ## ⚠️ PARITY TARGET CORRECTION (empirically verified 2026-07-04)
> The shipped `Demo/Example1-tRNAs.out` (Inf 74.2/81.6, Type=Leu, Isotype Score 119.9)
> is **NOT reproducible with the local `original/` model set** — it was generated with a
> different/newer model. Running the reference driver locally with the SAME columns
> (`-B -H --detail`) yields DIFFERENT values (Inf 61.3, Type=Undet, Isotype 61.1 for
> tRNA#1). Therefore the **real parity oracle is the reference implementation's LOCAL
> output**, regenerated on demand:
> ```
> PERL5LIB=<orig>/lib perl <orig>/tRNAscan-SE -B -H --detail -c <orig>/tRNAscan-SE.conf -o OUT input.fa
> ```
> (C cmsearch/cmscan = `bactars/infernal/original/src/`.) Saved canonical targets:
> `bactars/tRNAscan-SE/tests/reference_golden/Example1-Bhd.out`, `Example2-Bhd.out`
> (deterministic across reruns). All algorithm/milestone content below is valid; only
> substitute these local targets for the shipped `.out` files as the diff oracle.
> **Verified:** Inf Score column == C cmsearch `-g --nohmm --toponly --notrunc` bit score
> (e.g. 19480→68.3), so infernox needs a `-g --nohmm` search mode (not just alidisplay).

Original goal statement (target now the local-reference outputs above):
`bactars/tRNAscan-SE/original/Demo/Example1-tRNAs.out` and `Example2-tRNAs.out`.

**Source-of-truth trees**
- C/Perl reference: `/mnt/DAS/sunju/programme/bactars/tRNAscan-SE/original/` (driver `tRNAscan-SE`, `lib/tRNAscanSE/*.pm`, `tRNAscan-SE.conf`, `lib/tRNAscan-SE/models/`)
- Rust port under construction: `/mnt/DAS/sunju/programme/bactars/tRNAscan-SE/src/` (`core/`, `pipeline/`, `cm_scan/`, `isotype/`, `structure/`, `sprinzl/`, `trna/`)
- In-process Infernal engine: `/mnt/DAS/sunju/programme/infernox/infernal/src/` (`faithful_search.rs`, `cp9_faithful.rs`, `lib.rs`)

**Golden reference rows (Example1, the parity target):**
```
CELF22B7 	1	12619	12738	Leu	CAA	12657	12692	74.2	51.20	23.00	Inf	Leu	119.9	
CELF22B7 	2	19480	19561	Ser	AGA	0	0	81.6	47.50	34.10	Inf	Ser	125.0	
CELF22B7 	3	26367	26439	Phe	GAA	0	0	82.5	56.60	25.90	Inf	Phe	112.1	
CELF22B7 	4	26992	26920	Phe	GAA	0	0	82.5	56.60	25.90	Inf	Phe	112.1	
CELF22B7 	5	23765	23694	Pro	CGG	0	0	71.5	48.20	23.30	Inf	Pro	113.0	
```
Example2 exercises: minus/plus mix, an intron (MySeq3 `51 69`), a SeC hit with `HMM=2'Str=0.00` (MySeq5), and an **ISM note** (MySeq6 `Type=Lys` but `Isotype CM=Leu 75.7`, Note `ISM (-73.90)`). These two files together cover every column branch we must reproduce.

---

## 1. The exact ordered `-B` pipeline

All model files live in `cm_dir = {lib_dir}/models = original/lib/tRNAscan-SE/models/`. For `-B` (`set_options`, driver `:998-1029`; finalize `:1514-1568`): `search_mode="bacteria"`, `CM_mode="infernal"`, `opt_inf=1`, `opt_bact=1` ⇒ `infernal_fp=1`, `tscan_mode=0`, `eufind_mode=0`, `hmm_filter=0`, `no_isotype=0`, `CM_check_for_introns=0`, `CM_check_for_split_halves=0`. **Consequence: no tScan/EufindtRNA, no non-canonical-intron scan, no split-half scan.**

**Model set for `-B`** (loaded `CM.pm:425-440`, `sort keys` ⇒ deterministic order `Domain` then `SeC`):
| role | key | file |
|---|---|---|
| main domain CM | `main_cm_file_path{Domain}` | `TRNAinf-bact.cm` |
| main SeC CM | `main_cm_file_path{SeC}` | `TRNAinf-bact-SeC.cm` |
| no-structure CM (pseudogene) | `mainNS_cm_file_path{Domain}` | `TRNAinf-bact-ns.cm` |
| isotype CM db (cmpress'd) | `isotype_cm_db_file_path` | `TRNAinf-bact-iso` (+ `.i1f/.i1i/.i1m/.i1p`) |

**Resolved constants (`tRNAscan-SE.conf`, verified):** `infernal_fp_cutoff=10` (`:24`), `cm_cutoff=20` (`:26`), `isotype_cm_cutoff.bact=20` (`:31`), `default_Padding=10` (`Options.pm:123`), `upstream_len=downstream_len=70` (`:18-19`), `min_cmsearch_pseudo_filter_score=55` (`:50`), `min_ss_score=5`, `min_hmm_score=10` (`:52-54`), `BHB_cm_cutoff=6.5` (unused in `-B`).

**`exec_cmsearch` flag → options table** (`CM.pm:2530-2577`); base = `-g --nohmm --toponly --notrunc`:
| flag | options | phase |
|---|---|---|
| 3 | `-g --mid --notrunc -T 10` | **Phase I first-pass** |
| 0 | `-g --nohmm --toponly --notrunc` | **Phase II verify** (no `-T`; cutoff enforced in post) |
| 7 | `-g --max --toponly --notrunc --notextw -T 0` | pseudogene NS rescore (bact/arch) |
| 6 | `-g --toponly --notextw` | truncation labeling |
| 2 | `-g --mid --toponly --notrunc` (cmscan `--fmt 2 --tblout`) | isotype |

`-T score_cutoff` is appended for flags ∉{3,5,7} **only if `score_cutoff ≤ 10`** (`CM.pm:2565-2571`). In Phase II `score_cutoff=cm_cutoff=20` ⇒ **no `-T`**. `--cpu N` appended only if `--thread` given (`:2572`). Final command form (`:2577`): `cmsearch <opts> <cm_file> <tmp_seq> > <out>`.

### Ordered execution (driver `tRNAscan-SE` main loop `:111-315` then `:341-494`)

1. **Phase I — first-pass prescan** (`first_pass_prescan` `:185` → `CM.pm:3105 run_first_pass_cmsearch`). Log `"Phase I: Searching for tRNAs with HMM-enabled Infernal"`. For each model (Domain, then SeC), run **flag 3**: `cmsearch -g --mid --notrunc -T 10 <model> <full_input_seq> > out`. Input is the **whole input sequence**, not a window. Merge both models' hit tables → `*_fp_cm_merge.out`. `process_fp_cmsearch_hits` (`CM.pm:3137`) makes each hit a candidate; overlapping same-strand hits merged via `merge_repeat_hit` (`:3190-3215`), keeping extended span + latest score. Output: candidate regions (seqname, start, end, strand, score≥10, model). **Padding NOT applied here.**

2. **Phase II — candidate extraction** (`prepare_multi_tRNAs_to_scan` `:507`). Each candidate re-read (`FpScanResultFile.pm:297`), padded **±10 bp** (clamped `[1,seqlen]`, strand-aware, `:349-363`), then `get_tRNA_sequence` (`Sequence.pm:619`) extracts candidate **+70 bp upstream +70 bp downstream flank** (clamped at contig ends). Written one FASTA record per candidate to `tmp_trnaseq_file`, named `<seqname>.t<zeropad6(id)>`. Flanks stored separately.

3. **Phase II — verify** (`run_cm_scan` `:341` → `analyze_with_cmsearch` `CM.pm:3584` → `run_cmsearch` `:2591`). Log `"Phase II: Infernal verification..."`. For each model (Domain, SeC), run **flag 0**: `cmsearch -g --nohmm --toponly --notrunc <model> <tmp_trnaseq_file> > out`. Merge → `*_cm_merge.out` (written with `format=1` ⇒ SS normalization applied, §2). Per merged hit (matched to parent candidate by `seqname eq prescan.seqname().".t".pad_num(id,6)`, `:3626`):
   - Remap hit coords to genomic (`:3634-3662`).
   - **Reporting cutoff**: `score < cm_cutoff(20)` ⇒ reject `"Low cmsearch score"` (`:3664-3669`); else `curseq_trnact++`.
   - `set_domain_model("infernal", score)` (`:3677`) — seeds Inf score.
   - Bacterial CCA/3′ fixups (`:3685-3752`): append trailing `CCA` if present downstream, or trim spurious 3′ bases; if changed, `rescore_tRNA` (`:3754-3757`).
   - Flank re-trim block (`:3769-3776`) **skipped** (infernal_fp true) — flanks kept.
   - `decode_tRNA_properties` (`CM.pm:1258`, called `:3777`): anticodon, isotype, intron, pseudogene filter (§2).
   - Write accepted tRNA to `sp_int_results` (`:3789`).

4. **Phase III — truncation labeling** (only if `curseq_trnact>0`; `bact_mode` true ⇒ `truncated_tRNA_search` `:449`/`CM.pm:2653`). For each model, **flag 6**: `cmsearch -g --toponly --notextw <model> <tmp_trnaseq_file>`. Merged hits → `check_truncation` labels each tRNA's 5′/3′ trunc status; replaces `sp_int_results` (`:2694-2714`).

5. **Phase III — isotype scan** (`!no_isotype && bact_mode` ⇒ `isotype_cmsearch` `:453`/`CM.pm:2761` → `scan_isotype_cm` `:2797`). Write confirmed mature tRNAs to `tmp_trnaseq_file`; **flag 2 cmscan**: `cmscan -g --mid --toponly --notrunc --fmt 2 --tblout <tab> -o <out> TRNAinf-bact-iso <tmp_trnaseq_file>`. Consumes the **`--fmt 2 --tblout` table**, not alignment. Per-isotype bit scores appended as columns (§4). Bacterial isotype cutoff 20.

6. **Output** (`output_tRNA` `:457`): merge primary + hmm + isotype/trunc → final `.out` (§3).

7. **Reset** `curseq_trnact = 0` (driver `:494`) — **the tRNA-# counter is reset at the end of each source-sequence iteration.** (See §3.4 correction.)

---

## 2. Per-column derivation from a cmsearch hit + alignment

### 2.0 Alignment parsing (main search) — `CMscanResultFile.pm`
Main search parses the **human-readable alignment stdout via regex**, NOT tblout. Two passes:

**Pass A — hit summary** (`:113 sort_cmsearch_records`; only lines after `/Hit alignments:/`). Verbatim summary regex (`CMscanResultFile.pm:146`, dup `:259`; `CM.pm:3086`):
```
/^\s+\(\d+\)\s+\S+\s+([e0-9.\-]+)\s+([0-9.\-]+)\s+\S+\s+\S+\s+(\d+)\s+(\d+)\s+\S+\s+(\d+)\s+(\d+)\s+([+-])\s+\S+\s+\S+\s+(\S+)\s+([0-9.]+)/
```
`$1`=E-value, **`$2`=bit score (→ Inf Score)**, `$3/$4`=model coords, **`$5/$6`=seq from/to**, **`$7`=strand**, `$8`=trunc, `$9`=GC%. Target from `/^>>\s*(\S+)/` (`:141`). If strand `-`, swap start/end (`:154-159`); `trunc="no"`→`""`.

**Pass B — alidisplay reconstruction** (`:224 get_cmsearch_record`). Per alignment chunk, extract and concatenate:
- **SS_cons (CS) line** `:278`: `/^\s{5,}([(),<>._\-,\[\]\{\}\:\~]{1,250}) CS$/` → `$ss`.
- **NC line** `:274`: `/^(.+) NC$/` → marks low-confidence `v` positions.
- **consensus/model line** `:283`: `/^\s+\S+\s+\d+\s+([a-zA-Z\.0-9\>\<\[\]\*]{1,250})\s+\d+/` → `$model`.
- **aligned target seq line** (skip match line, `:289-291` replace `*[0]*`→`-----`) `:292-293`: `/^\s+\S+\s+\d+\s+([a-zA-Z\-]{1,250})\s+\d+/` → `$seq`.
- PP and RF lines matched and **discarded** (`:305-310`).

### 2.1 SS normalization (order-dependent, LOSSY) — `format_cmsearch_output` `:321`
Applied on merge (format=1). Exact order:
1. `seq`: `U→T`, `u→t` (`:327`).
2. `ss = fix_mismatch_ss(ss, seq, nc)` (`:360`): NC `v`→`.`; stack-pair `<>`/`()`; for each pair, if bases not Watson-Crick (G:U/T allowed) or either is `-`, demote both to `.` (`:405-415`).
3. For every `-` gap in `seq`, mark seq+ss with `*`, then strip all `*` from both — **deletes gap columns from both strings** (`:332-344`).
4. SS char remap (`:346-349`): `[,_\-:]→.`; then `[>)]→@`, `[(<]→>`, `@→<`. Net: pair-open→`>`, pair-close→`<`, unpaired→`.` (the `>>>...<<<` convention).
5. Pad ss with `.` to seq length (`:351-355`).

**All downstream anticodon/intron regexes operate on this post-normalization string.**

Merge-file tab columns (`ArrayCMscanResults.pm:260-278`): `seqname  start  end  strand  score  trunc  ss  seq  model  nc  type`. Cross-model overlapping hits deduped keeping higher score (`merge_indexes` `:179-200`).

### 2.2 Column: Anticodon + Type (isotype)
**`find_anticodon`** (`CM.pm:731`), input normalized `seq`,`ss`. Core regex (`:748`):
```perl
if ($ss =~ /^([>.]+<[<.]+>[>.]*)>([.]{4,})<+.+[>.]+<[<.]+/o) {
    $antiloop_index = length($1) + 1;   # 0-based start of AC loop
    $antiloop_len   = length($2);       # loop length (run of >=4 '.')
}
```
Then (`:755-789`):
```perl
$antiloop_end = $antiloop_index + $antiloop_len - 1;
$antiloop = substr($seq, $antiloop_index, $antiloop_len);
$antiloop =~ s/[\-]//g;  $antiloop =~ s/[a-z]//g;      # strip gaps + introns/lowercase
if (length($antiloop) < 5 || length($antiloop)%2 == 0) { return undef; }
$ac_index  = (length($antiloop) - 3) / 2;              # center 3 nt
$anticodon = substr($antiloop, $ac_index, 3);
$verify_ac = substr($seq, $ac_index + $antiloop_index, 3);
if ($verify_ac ne $anticodon) { category("undetermined_ac"); return undef; }
return ($anticodon, $antiloop_index, $antiloop_end, $ac_index+$antiloop_index+1);
```
`undef` ⇒ anticodon `"NNN"` (`GeneticCode.pm:35`). `decode_tRNA_properties` stores `anticodon()` and `add_ac_pos` (`CM.pm:1274-1276`).

**Anticodon → Type** — `GeneticCode::get_tRNA_type` (`:267`): `NNN`→`"Undet"`; SeC CM path→`"SeC"`; else expand ambiguity (`expand_ambig`) and look up **`trans_map`** (built from **reverse-complement of the `__DATA__` codon table**, `:190-217`) per expansion; disagreement→`"Undet"`. If `type=="SeC"` && model≠SeC && !cove → downgrade to `"Sup"` (`:299-302`). Post-check (`CM.pm:1332`): if `anticodon ne "TCA"` but isotype `SeC`, reset both to undef. `tRNAscan_id = seqname.tRNA<N>-<isotype><anticodon>` (`:3779`).

### 2.3 Column: Inf / HMM / 2'Str scores
- **Inf Score** = cmsearch bit score (`$2`) → merge col[4] → `set_domain_model("infernal")` → `tRNA->score()`. Output raw (no printf) `ScanResult.pm:723`.
- **HMM Score & 2'Str**: only when pseudogene filter runs — `is_pseudo_gene` (`CM.pm:999`), engaged iff `score < 55` OR `-H` (`:1017-1025`). Rescore vs NS model `TRNAinf-bact-ns.cm` via `cmsearch_scoring` (`:3061`, **flag 7** for bact: `-g --max --toponly --notrunc --notextw -T 0`), parse same summary regex, take **max bit score** across hits (`besthit_score` `:3094-3099`) → `hmm_score`. Then `ss_score = Inf_score − hmm_score` (`:1053`). `update_domain_model` → `tRNA->hmm_score()/ss_score()`. Output `sprintf "\t%.2f\t%.2f"` (`ScanResult.pm:727`). Pseudo flag if `(ss_score<5 || hmm_score<10) && score<55` (`:1061-1064`). *(SeC hits like MySeq5 skip the filter ⇒ `0.00 0.00`.)*

### 2.4 Column: Intron Bounds
**`find_intron`** (`CM.pm:960`, called `:1301` with `(seq, antiloop_index, antiloop_end)`). If `antiloop_index==-1`, no intron. Else (`:977`):
```perl
$antiloop_seq = substr($trna_seq, $antiloop_index, $antiloop_end-$antiloop_index+1);
if ($antiloop_seq =~ /^(.*[^a-z]+)([a-z]{$min_intron_length,})[^a-z]+/o) {
    $intron = $2;                                # lowercase run inside AC loop
    $istart = index(substr($trna_seq,0,$antiloop_end+1), $intron) + 1;  # 1-based
    $iend   = length($intron) + $istart - 1;
}
```
Genomic coords (`:1307-1317`): `+` → `intron_start=istart+trna.start-1`, `intron_end=iend+trna.start-1`; `-` → `intron_start=trna.end-iend+1`, `intron_end=trna.end-istart+1`. `add_intron(istart,iend,gstart,gend,"CI",intron)`. **No NCI/BHB scan in `-B`.**

### 2.5 Columns: Isotype CM / Isotype Score (last two)
From the **cmscan `--fmt 2 --tblout`** run (§1.5). tblout parse `get_next_tab_seq_hits` (`CMscanResultFile.pm:421`): per hit record `[0]=col[1]` target/model, `[1]=col[3]` query, `[2]=col[9]/[3]=col[10]` coords, `[4]=col[11]` strand, **`[5]=col[16]` bit score** (0-indexed 16 ⇒ field 17). Assemble one column per isotype model (`scan_isotype_cm:2836-2932`, mito prefixed `mito_`). At output `get_highest_score_model` (`tRNA.pm:1130`) returns top (model,score) by desc score → **Isotype CM = model name**, **Isotype Score = its bit score**. Met special-case (`ScanResult.pm:311-327`): iMet/fMet/Ile2 promotion (Ile2 needs `(score−ile2)≤5 && (ile2−met)≥5 && tRNA.score>50`). **ISM Note** emitted when `detail && model set && isotype≠Undet && model ne isotype` (excluding alias cases): `sprintf("IPD:%0.2f",$iso_score−$score)` — in Example2 rendered `ISM (-73.90)`; port must reproduce the sign/format exactly.

---

## 3. `.out` output format (`ScanResult.pm construct_tab_output :639-868`, write `:370`)

Field separators are literal TAB; only spaces are `%-Ns`/`%-Nd` padding. Demo mode = infernal + isotype-detail + `get_hmm_score=1` + `save_source=1`, `infernal_score()=false`.

### 3.1 Column widths (computed once from FIRST tRNA; frozen file-globals, `:344-348`)
```
$max_seq_name_width = max(length(src_seqid)+1, 8);   # CELF22B7=8 → 9
$max_seq_len_width  = length(src_seqlen());          # digits in seq length
```

### 3.2 Header (printed once, unless `brief_output`; `:81-197`). Exact byte layout for demo mode:
```
Line1: "Sequence"(%-9s)"\t\t" "tRNA"(%-Ls)"\t" "Bounds"(%-Ls)"\t" "tRNA\tAnti\tIntron Bounds" "\tInf" "\tHMM\t2'Str" "\tHit" "\tIsotype\tIsotype" "\t      " "\n"
Line2: "Name"(%-9s)"\t" "tRNA #\t" "Begin"(%-Ls)"\t" "End"(%-Ls)"\t" "Type\tCodon\tBegin\tEnd\tScore" "\tScore\tScore" "\tOrigin" "\tCM\tScore" "\tNote" "\n"
Line3: "--------"(%-9s)"\t" "------\t" "-----"(%-Ls)"\t" "------"(%-Ls)"\t" "----\t-----\t-----\t----\t------" "\t-----\t-----" "\t------" "\t-------\t-------" "\t------" "\n"
```
(`L = max_seq_len_width`.) Conditional blocks must mirror data-row conditionals exactly: `\tInf`/`\tCove`/`\tEufind` by mode; `\tHMM\t2'Str` iff `get_hmm_score`; a *second* `\tInf` iff `infernal_score()` (off in demo); `\tHit` iff `save_source`; `\tIsotype\tIsotype` iff `(euk|bact|arch)&&!no_isotype&&detail`. With `output_codon()` the `Anti` label → `"   "`.

### 3.3 Data row (`:644-865`)
Pre-computed: `($type,$model,$score,$ss)=get_highest_score_model()` (isotype CM pair); `($iso_model,$iso_score,$iso_ss)=get_model_hit("cyto",isotype())` (for IPD note).
| # | Column | code (line) | format |
|---|---|---|---|
|1|Seq Name|`sprintf "%-${w}s\t", seqname()` (648)|left-just min 8|
|2|tRNA #|`id()."\t"` (649)|int|
|3|Begin|(653/658)|`%-Ld` left-just|
|4|End|(654/659)|`%-Ld`|
|5|Type|`isotype()."\t"` (662)|string|
|6|Codon|anticodon (or revcomp if `output_codon`) (664-671)|string|
|7|Intron Begin|(674-722)|int(s) or `0`|
|8|Intron End|(674-722)|int(s) or `0`|
|9|Score|`"\t".score()` (723)|**raw stored value, no printf**|
|10|HMM Score|`sprintf "\t%.2f"` (727)|`%.2f` (iff get_hmm_score)|
|11|2'Str Score|`sprintf "\t%.2f"` (727)|`%.2f` (iff get_hmm_score)|
|(Inf)|Inf Score|`"\t".$inf->{score}` (734)|raw (iff infernal_score; off)|
|12|Hit Origin|`"\t".hit_source()` (739)|`Inf` (iff save_source)|
|13|Isotype CM|`"\t".$model` (743)|model name|
|14|Isotype Score|`"\t".$score` (744)|raw|
|15|Note|`$note` (838)|leading `\t` + note or empty|
Row ends `.= "\n"` (865).

**Begin/End per strand** (`:651-660`): `+` ⇒ Begin=start, End=end; `-` ⇒ Begin=end, End=start (so minus strand has Begin>End — Example1 rows 4,5). **Intron Begin/End** (`:673-722`): 0 introns ⇒ `"0\t0"`; else comma-joined lists, `+`: Begin=`{start}`,End=`{end}`; `-`: swapped. **Note** (`:758-864`): `pseudo` if pseudo; then ISM/IPD (§2.5); then trunc text; (archaea-only intron tags — N/A for `-B`). Demo rows have empty Note (row still ends `\t\n`).

### 3.4 Numbering & sort order — **CORRECTION to a source report**
The output-format report claimed `curseq_trnact` is "never reset between sequences." **This is wrong.** Driver line `494` executes `$curseq_trnact = 0;` at the end of every source-sequence iteration, so **tRNA # restarts at 1 for each FASTA record** — confirmed by Example2 (MySeq1..6 each show `tRNA # = 1`). The Rust port MUST reset the counter per source sequence.

**Sort order within a sequence** (`ArrayCMscanResults.pm:203-219 sort_by_tRNAscanSE_output`; re-index `IntResultFile:336-389`): key = ordered_seqname, then `+` strand ascending `start` / `-` strand descending `end`, then descending score. `output_tRNA` iterates `sp_int_results->get_indexes()` in that order and assigns id.

### 3.5 Secondary outputs (out of primary scope but same data model)
`.ss` (`save_allStruct_output :408-608`), `.stats` (`Stats.pm :282+`), `.bed` (BED12 `write_bed :997`, `convert_bed_score=min(1000,max(0,score*10))`), `.gff` (GFF3 `write_gff :1118`). The faithful Rust versions already exist in `pipeline/output.rs` (§5).

---

## 4. infernox `FaithfulSearcher` alidisplay requirement

**Why:** every §2 column is derived from the cmsearch **alignment** (SS_cons + aligned target seq + lowercase introns), not just tblout coords/score. The Rust port therefore needs `FaithfulSearcher` to expose, per hit, a faithful `cm_alidisplay` equivalent. Today the F7 stage (`faithful_search.rs:512-545`) runs banded **CYK** (`cyk_align_hb_cmbounds`, `cp9_faithful.rs:4370`) and **discards the traceback**, returning only `(cfrom_emit, cto_emit)` — no aligned strings, no SS, no PP. The C display path uses **HMM-banded OptAcc + posterior decoding** (`pli_align_hit`, `cm_AlignHB`), not CYK, so current bounds are only an approximation.

### 4.1 Fields `FaithfulSearcher` MUST expose (add `alignment: Option<CmAliDisplay>` to `FaithfulHit`, `faithful_search.rs:59`)
| field | consumer in §2 | C source |
|---|---|---|
|`aseq` (aligned target residues, lowercase=insert/intron)|§2.1 seq, §2.4 intron lowercase run|`cm_alidisplay.c:408/429/443/454/468`|
|`csline` (SS_cons)|§2.1 ss normalization → §2.2/2.4 regexes|`cm->cmcons->cstr` `:451/465`|
|`model` (aligned consensus query)|§2.0 model line (used by fix_mismatch_ss / display)|`cm->cmcons->cseq` `:452/466`|
|`ppline` (posterior-prob string)|discarded by main parse (PP), but needed for `.ss`/faithful display parity|`cm_PostCodeHB` `cm_dpalign.c:784`|
|`ncline` (noncanonical bp `v` marks)|§2.1 `fix_mismatch_ss` NC→`.`|display traversal|
|`cfrom_emit`/`cto_emit` (model from/to)|alignment coords|`ParsetreeToCMBounds` (already inlined `cp9_faithful.rs:4715`)|
|`sqfrom/sqto,sc,gc,avgpp,clen`|coords/score/GC%|already on `FaithfulHit` / trivially derivable|

The port must **reproduce the exact `cmsearch` stdout text** these regexes consume, OR expose structured equivalents and re-implement §2.0 parsing against them (preferred — avoids brittle text round-trip).

### 4.2 Implementation sketch (new `infernox/.../cm_alidisplay.rs`)
Data types: `FaithfulParsetree{state,emitl,emitr,mode,nxtl,nxtr,prv}` (faithful twin of `legacy/parsetree.rs:12` + `mode`); `CmConsensus{cseq,cstr,lpos,rpos}` (port `CreateCMConsensus`, `display.c:653` — argmax `cm.esc[v]` per MATP/MATL/MATR, lowercase below `pthresh=3.0/sthresh=1.0`, structure chars from multifurcation-order chart; `lpos/rpos` reuse existing `CMEmitMap`); `CmAliDisplay{aseq,model,csline,mline,ppline,cfrom_emit,cto_emit,sqfrom,sqto,sc,avgpp,gc,clen}`.

Missing DP in dependency order (the bulk of the work; legacy non-parity refs at `cm_dp_hb.rs:812/1381/1497` are a *reference*, must be re-ported onto the **faithful** `CP9Bands` `cp9_faithful.rs:1373` after `shift_cm_bands` `:4331`):
1. Extend traceback to build a real `FaithfulParsetree` (record `emitl=i, emitr=j, state=v, mode=TRMODE_J` — the existing traceback `cp9_faithful.rs:4653-4713` already computes i/j/v, it just drops them).
2. Port OptAcc+Post quartet: `cm_inside_align_hb` → `cm_outside_align_hb` → `cm_posterior_hb` → `cm_emitter_posterior_hb` (l_pp/r_pp) → `cm_optacc_align_hb` + OptAcc traceback → `cm_postcode_hb` (`Fscore2postcode`: `(p+0.05>=1.0)?'*':(char)((p+0.05)*10)+'0'`, exact f32). The banded loop skeleton / `hdmin/hdmax` indexing / EL handling in `cyk_align_hb_cmbounds` (`:4400-4645`) is the template — same band arithmetic, `FLogsum` for Inside/Outside, max-of-posterior for OptAcc.
3. Port `cm_alidisplay_Create` J-mode traversal (`cm_alidisplay.c:56`, main loop `:339-479`, EL `:365-401`, length `:148-175`): PDA walk emitting left/right columns per node from `cmcons.lpos/rpos`, `dsq[emitl/emitr]`, `ppstr`; `mline` from emission-score sign.

Gate behind `Option`/config flag (roughly doubles per-hit DP but only on F7 survivors ⇒ small total overhead). f32 exact arithmetic mandatory for byte-parity (`Fscore2postcode`, EL PP-averaging `cm_alidisplay.c:382-398`).

---

## 5. Gap analysis vs existing Rust (`bactars/tRNAscan-SE/src/`)

**Headline:** two disconnected stacks. The binary (`bin/trnascan.rs`) wires the **heuristic** `core::TrnaScanner` (`core/scanner.rs`); a **more faithful but dead** stack (`pipeline/`, `cm_scan/` decode fns, `trna/`, `structure/konings`, `sprinzl/`) has no caller. Fastest route to parity = wire the faithful modules, delete the heuristics.

| module | state | recommendation |
|---|---|---|
| `core/scanner.rs` orchestration (FP→CM 2-pass) | active | **keep skeleton, rewire** decode/isotype/intron/output steps |
| `core/scanner.rs::extract_anticodon_from_seq` (T[NNN]A motif, ignores SS) | active | **delete → replace** with SS-based `find_anticodon` |
| `core/scanner.rs::determine_isotypes` (shells out to external `cmscan`, best-score only, no cutoffs/ISM) | active | **rewrite** on in-process `CMScan::scan_isotype` (`cm_scan/mod.rs:1043`) + faithful `get_tRNA_type` |
| `core/scanner.rs::TrnaResult` + its formatters (crude, no header, wrong `.out`/`.ss`/`.bed`) | active | **delete → replace** with `TRna` + `pipeline/output.rs` |
| `cm_scan::find_anticodon`/`find_anticodon_loop` | dormant | **reuse `find_anticodon`; REWRITE `find_anticodon_loop`** to exact Perl regex (`CM.pm:748`); fix off-by-one (Perl `antiloop_index=len($1)+1`) |
| `cm_scan::find_intron` (lowercase-run concept) | dormant | **reuse once alignment exposed**; NCI/BHB not needed for `-B` |
| `cm_scan::is_pseudogene` (NS rescore, ss=total−hmm) | dormant | **reuse**; register `TRNAinf-bact-ns.cm`, call per hit |
| `cm_scan::batch_search_external` (in-process FaithfulSearcher) | active | **extend** to surface aseq+SS+PP (§4) so decode can run |
| `isotype/anticodon.rs` (static table; TCA→SeC unconditional, no Sup/undet/iMet/Ile2) | dormant | **rewrite** to mirror `GeneticCode::get_tRNA_type` |
| `isotype/scorer.rs` (invented heuristics, placeholder `-999`) | dead | **discard** |
| `structure/konings.rs` (faithful port of `konings.c`) | dead | keep (low priority for `-B` MVP) |
| `sprinzl/mod.rs` (Sprinzl 1-76, unverified vs `SprinzlPos.pm`) | dead | keep (low priority; Note-field mismatch annotation not needed for demo) |
| `pipeline/output.rs` (**faithful** port of `ScanResult.pm`/`Stats.pm`) | dead | **adopt as output layer**; fix `chrono_lite_now` timestamp |
| `pipeline/scanner.rs` (alt orchestration, `TrnaHit`) | dead | fold useful bits into `core`; avoid 3rd result type |
| `trna/mod.rs` (`TRna`, faithful port of `tRNA.pm`) | orphan | **promote to single result type** |

**Cross-cutting gaps the active path cannot currently reproduce:** (1) anticodon from CM SS; (2) intron coords (always 0 — decode never runs); (3) in-process isotype CM scoring + per-isotype cutoffs + ISM note; (4) pseudogene HMM/2'Str breakdown (cols 10-11 currently 0/absent); (5) faithful `.out` header + `.stats` + intron-block `.bed`; (6) Sprinzl Note annotation (not needed for demo). Net fix = expose alignment from FaithfulSearcher → run existing `decode_trna_properties` + `is_pseudogene` → route isotype through `scan_isotype` → swap output layer to `pipeline/output.rs`, consolidating on `TRna`.

---

## 6. Milestone plan with byte-parity checkpoints

Parity harness (build first): `diff <(rust -B Demo/Example1.fa) Demo/Example1-tRNAs.out` and same for Example2; assert **byte-identical**. Diff column-by-column during development.

**Critical path (sequential):** M1 → M2 → M3 → M6 → M8. M4, M5, M7 parallelizable.

| M | Deliverable | Parity checkpoint | Parallel? |
|---|---|---|---|
| **M1** | **FaithfulSearcher alidisplay** (§4): OptAcc+Post DP quartet, `FaithfulParsetree`, `CmConsensus`, `cm_alidisplay_Create`; expose `aseq/csline/model/ppline/cfrom/cto` on `FaithfulHit`. | Dump one hit's alidisplay; byte-match `cmsearch` stdout alignment block for a known CELF22B7 tRNA. | No — gates everything. Biggest/riskiest. |
| **M2** | **Pipeline wiring** (§1): `core::scanner` runs Phase I (flag 3) → extract ±10/±70 → Phase II (flag 0) → cutoff 20 → CCA fixups → truncation (flag 6). Consolidate on `TRna`; **reset counter per sequence** (§3.4). | Correct **count & genomic coords** (cols 1-4) for all 5 Example1 rows; sort/order/id match. | No — depends M1 for hits. |
| **M3** | **Decode: anticodon + Type + intron** (§2.1/2.2/2.4): port SS normalization (`format_cmsearch_output`) exactly; fix `find_anticodon_loop` to Perl regex; wire `find_intron`; port `get_tRNA_type` + `trans_map` (revcomp codon table). | Cols 5,6,7,8 match Example1 (incl. Leu intron `12657 12692`) + Example2 MySeq3 `51 69`. | Partly — decode logic (SS/regex/genetic-code) can be built against fixtures in parallel with M2; integration needs M1+M2. |
| **M4** | **Pseudogene / HMM / 2'Str** (§2.3): register `TRNAinf-bact-ns.cm`; run `is_pseudogene` (flag 7, `-T 0`, max bit score) when score<55; `%.2f` format. | Cols 10,11 match all Example1 rows + SeC `0.00 0.00` (MySeq5). | **Yes** — independent once M1 alidisplay/score exists. |
| **M5** | **Isotype cmscan** (§2.5): in-process `scan_isotype` vs `TRNAinf-bact-iso` (flag 2, `--fmt 2 --tblout`); parse field 17; `get_highest_score_model`; ISM/IPD note + Met/Ile2 special-case. | Cols 13,14 + Note match Example1 + Example2 MySeq6 `ISM (-73.90)`. | **Yes** — independent module; needs mature-seq input from M2. |
| **M6** | **Output layer** (§3): adopt `pipeline/output.rs`; exact header bytes, widths from first tRNA, per-strand Begin/End, raw vs `%.2f` fields, Note. | Full `.out` header + all Example1 & Example2 rows byte-identical. | No — integrates M2-M5. |
| **M7** | **Secondary outputs** `.ss/.stats/.bed/.gff` (§3.5): wire faithful formatters; fix `chrono_lite_now`. | Byte-match Demo `.ss`/`.bed` if golden provided (non-`.out` deterministic parts). | **Yes** — after M6 data model stable. |
| **M8** | **Full parity + hardening**: delete dead heuristics (`extract_anticodon_from_seq`, `determine_isotypes`, `IsotypeScorer`, `TrnaResult` formatters); CI parity gate on both demos. | `diff` == empty for Example1 **and** Example2. | No — final gate. |

**Parallelization:** after M1 lands, three agents can run concurrently — Agent A on M2+M3 (critical path), Agent B on M4 (pseudogene), Agent C on M5 (isotype). M6 joins them; M7 follows; M8 closes.
