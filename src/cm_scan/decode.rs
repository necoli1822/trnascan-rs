//! Faithful port of the tRNAscan-SE `-B` decode logic (anticodon / isotype /
//! intron) from a single cmsearch alignment.
//!
//! This is a self-contained, byte-parity-intent reimplementation of the Perl
//! decode path, taking the fields of one `cmsearch` alidisplay block:
//!   - `aseq`    : aligned target-sequence line
//!   - `ss_cons` : CS (SS_cons) line
//!   - `nc`      : NC line (low-confidence / non-canonical `v` marks)
//!   - `model`   : consensus/model line (an alidisplay field; not consumed by
//!                 the faithful decode itself, kept for interface completeness)
//!
//! Verbatim Perl sources reproduced here:
//!   - `CMscanResultFile.pm::format_cmsearch_output` (~321) + `fix_mismatch_ss` (~360)
//!   - `CM.pm::find_anticodon` (~731), `find_intron` (~960)
//!   - `GeneticCode.pm::get_tRNA_type` (~267), `trans_map` build (~190-217),
//!     `expand_ambig` / `rev_comp_seq`
//!
//! The order of operations in SS normalization is LOSSY and order-dependent;
//! it is preserved exactly.

use once_cell::sync::Lazy;
use regex::Regex;
use std::collections::BTreeMap;

use crate::trna::{Intron, Strand};

/// Minimum canonical intron length (`CM.pm` `min_intron_length`, conf = 3).
/// `find_intron` (CM.pm:977) requires a lowercase (insert) run of at least this
/// many bases in the anticodon loop, so a 3 bp intron must be detected.
pub const MIN_INTRON_LENGTH: usize = 3;

/// Undetermined anticodon sentinel (`GeneticCode.pm:35`).
pub const UNDEF_ANTICODON: &str = "NNN";
/// Undetermined isotype sentinel (`GeneticCode.pm:36`).
pub const UNDEF_ISOTYPE: &str = "Undet";

// ============================================================================
// Input / output interface
// ============================================================================

/// The four alidisplay fields consumed from one cmsearch hit alignment.
///
/// These correspond to the strings parsed in `CMscanResultFile.pm`
/// `get_cmsearch_record` (:224). The integration layer (infernox alidisplay)
/// must provide `aseq`, `ss_cons`, and `nc`; `model` is accepted for interface
/// parity but is not read by the faithful decode.
#[derive(Debug, Clone)]
pub struct AliDisplay {
    /// Aligned target residues (lowercase = insert/intron), raw from cmsearch.
    pub aseq: String,
    /// SS_cons / CS line (`(`,`)`,`<`,`>`,`.`,`,`,`_`,`-`,`:`,`~`,`[`,`]`,`{`,`}`).
    pub ss_cons: String,
    /// NC line; `v` marks flag columns to demote in `fix_mismatch_ss`.
    pub nc: String,
    /// Consensus/model line (alidisplay field; unused by decode).
    pub model: String,
}

/// Decoded tRNA properties.
#[derive(Debug, Clone, PartialEq)]
pub struct DecodedTrna {
    /// Anticodon (uppercase, T-form), or `"NNN"` if undetermined.
    pub anticodon: String,
    /// Isotype / Type string (e.g. `"Ser"`, `"Sup"`, `"Undet"`).
    pub isotype: String,
    /// 0-based start index of the anticodon loop in the normalized seq, or -1.
    pub antiloop_index: i32,
    /// 0-based end index of the anticodon loop, or -1.
    pub antiloop_end: i32,
    /// 1-based position (in normalized seq) of the anticodon's first base, or -1.
    pub ac_pos: i32,
    /// Intron, if a canonical lowercase-run intron was found in the AC loop.
    pub intron: Option<Intron>,
    /// The normalized (post-`format_cmsearch_output`) sequence.
    pub norm_seq: String,
    /// The normalized (post-`format_cmsearch_output`) secondary structure.
    pub norm_ss: String,
}

// ============================================================================
// §2.1  SS normalization — format_cmsearch_output / fix_mismatch_ss
// ============================================================================

/// Port of `CMscanResultFile.pm::fix_mismatch_ss` (:360).
///
/// Operates on the RAW (pre-gap-removal) alignment columns. `ss`, `seq`, and
/// `nc` must be column-aligned. `seq` is expected to have already had U->T /
/// u->t applied (as in `format_cmsearch_output`).
fn fix_mismatch_ss(ss: &str, seq: &str, nc: &str) -> String {
    let mut ss: Vec<u8> = ss.bytes().collect();
    let seq: &[u8] = seq.as_bytes();
    let nc: &[u8] = nc.as_bytes();

    // NC 'v' columns -> demote to '.'
    for (i, &c) in nc.iter().enumerate() {
        if c == b'v' && i < ss.len() {
            ss[i] = b'.';
        }
    }

    // Pair up '<'/'(' with '>'/')' via the exact Perl stack walk.
    let mut left: Vec<usize> = Vec::new();
    let mut right: Vec<usize> = Vec::new();
    // pairs[left_stack_index] = right_stack_index, or -1 if unpaired.
    let mut pairs: Vec<i64> = Vec::new();

    for (pos, &c) in ss.iter().enumerate() {
        if c == b'<' || c == b'(' {
            left.push(pos);
            pairs.push(-1);
        } else if c == b'>' || c == b')' {
            right.push(pos);
            let mut li: i64 = left.len() as i64 - 1;
            while li > -1 && pairs[li as usize] > -1 {
                li -= 1;
            }
            if li > -1 && pairs[li as usize] == -1 {
                pairs[li as usize] = right.len() as i64 - 1;
            }
        }
    }

    // For each paired column, demote both ends if not Watson-Crick (G:U/T wobble
    // allowed) or if either base is a gap.
    for li in 0..pairs.len() {
        let ri = pairs[li];
        if ri < 0 {
            continue;
        }
        let lpos = left[li];
        let rpos = right[ri as usize];
        let lb = seq.get(lpos).copied().unwrap_or(b'-').to_ascii_uppercase();
        let rb = seq.get(rpos).copied().unwrap_or(b'-').to_ascii_uppercase();

        let bad = (lb == b'A' && rb != b'U' && rb != b'T')
            || (lb == b'T' && rb != b'A' && rb != b'G')
            || (lb == b'U' && rb != b'A' && rb != b'G')
            || (lb == b'G' && rb != b'C' && rb != b'U' && rb != b'T')
            || (lb == b'C' && rb != b'G')
            || (lb == b'-')
            || (rb == b'-');
        if bad {
            ss[lpos] = b'.';
            ss[rpos] = b'.';
        }
    }

    String::from_utf8(ss).unwrap()
}

/// Port of `CMscanResultFile.pm::format_cmsearch_output` (:321).
///
/// Returns the normalized `(ss, seq)` pair used by all downstream anticodon /
/// intron regexes. Lossy and order-dependent — order preserved exactly.
pub fn format_cmsearch_output(ss_cons: &str, aseq: &str, nc: &str) -> (String, String) {
    // (a) seq: U->T, u->t
    let mut seq: Vec<u8> = aseq
        .bytes()
        .map(|b| match b {
            b'U' => b'T',
            b'u' => b't',
            other => other,
        })
        .collect();
    let seq_str = String::from_utf8(seq.clone()).unwrap();

    // (b) fix_mismatch_ss
    let mut ss: Vec<u8> = fix_mismatch_ss(ss_cons, &seq_str, nc).into_bytes();

    // (c) delete every '-' gap column from BOTH seq and ss.
    for i in 0..seq.len() {
        if seq[i] == b'-' {
            seq[i] = b'*';
            if ss.len() > i {
                ss[i] = b'*';
            }
        }
    }
    seq.retain(|&b| b != b'*');
    ss.retain(|&b| b != b'*');

    // (d) SS char remap: [,_-:] -> '.' ; then [>)] -> '@', [(<] -> '>', '@' -> '<'.
    for b in ss.iter_mut() {
        match *b {
            b',' | b'_' | b'-' | b':' => *b = b'.',
            _ => {}
        }
    }
    for b in ss.iter_mut() {
        match *b {
            b'>' | b')' => *b = b'@',
            _ => {}
        }
    }
    for b in ss.iter_mut() {
        match *b {
            b'(' | b'<' => *b = b'>',
            _ => {}
        }
    }
    for b in ss.iter_mut() {
        if *b == b'@' {
            *b = b'<';
        }
    }

    // (e) pad ss with '.' to seq length.
    if seq.len() > ss.len() {
        ss.extend(std::iter::repeat(b'.').take(seq.len() - ss.len()));
    }

    (
        String::from_utf8(ss).unwrap(),
        String::from_utf8(seq).unwrap(),
    )
}

// ============================================================================
// §2.2  find_anticodon
// ============================================================================

/// The anticodon-stem-loop regex from `CM.pm:748` (operates on normalized ss).
static AC_RE: Lazy<Regex> =
    Lazy::new(|| Regex::new(r"^([>.]+<[<.]+>[>.]*)>([.]{4,})<+.+[>.]+<[<.]+").unwrap());

/// Port of `CM.pm::find_anticodon` (:731).
///
/// Returns `(anticodon, antiloop_index, antiloop_end, ac_pos)` where indices are
/// into `norm_seq`. On failure returns `("NNN", -1, -1, -1)`.
pub fn find_anticodon(norm_seq: &str, norm_ss: &str) -> (String, i32, i32, i32) {
    let caps = match AC_RE.captures(norm_ss) {
        Some(c) => c,
        None => return (UNDEF_ANTICODON.to_string(), -1, -1, -1),
    };
    // antiloop_index = length($1) + 1 ; antiloop_len = length($2)
    let antiloop_index = caps.get(1).unwrap().as_str().len() + 1;
    let antiloop_len = caps.get(2).unwrap().as_str().len();

    if antiloop_index == 0 || antiloop_len == 0 {
        return (UNDEF_ANTICODON.to_string(), -1, -1, -1);
    }

    let seq = norm_seq.as_bytes();
    let antiloop_end = antiloop_index + antiloop_len - 1;

    // antiloop = substr(seq, antiloop_index, antiloop_len)
    let raw_loop: &[u8] = if antiloop_index < seq.len() {
        let end = (antiloop_index + antiloop_len).min(seq.len());
        &seq[antiloop_index..end]
    } else {
        &[]
    };
    // strip '-' gaps, then strip lowercase (introns / non-canonical)
    let stripped: Vec<u8> = raw_loop
        .iter()
        .copied()
        .filter(|&b| b != b'-' && !b.is_ascii_lowercase())
        .collect();

    if stripped.len() < 5 || stripped.len() % 2 == 0 {
        return (UNDEF_ANTICODON.to_string(), -1, -1, -1);
    }

    let ac_index = (stripped.len() - 3) / 2;
    let anticodon = String::from_utf8(stripped[ac_index..ac_index + 3].to_vec()).unwrap();

    // verify_ac = substr(seq, ac_index + antiloop_index, 3) -- RAW slice, no strip.
    let vstart = ac_index + antiloop_index;
    let verify_ac = if vstart + 3 <= seq.len() {
        String::from_utf8(seq[vstart..vstart + 3].to_vec()).unwrap()
    } else {
        String::new()
    };
    if verify_ac != anticodon {
        return (UNDEF_ANTICODON.to_string(), -1, -1, -1);
    }

    (
        anticodon,
        antiloop_index as i32,
        antiloop_end as i32,
        (ac_index + antiloop_index + 1) as i32,
    )
}

// ============================================================================
// §2.4  find_intron
// ============================================================================

static INTRON_RE: Lazy<Regex> =
    Lazy::new(|| Regex::new(&format!(r"^(.*[^a-z]+)([a-z]{{{},}})[^a-z]+", MIN_INTRON_LENGTH)).unwrap());

/// Port of `CM.pm::find_intron` (:960).
///
/// Returns `(intron_seq, istart, iend)` in 1-based coordinates within
/// `norm_seq`, or `("", 0, 0)` if none. `antiloop_index` == -1 => no intron.
pub fn find_intron(norm_seq: &str, antiloop_index: i32, antiloop_end: i32) -> (String, i32, i32) {
    if antiloop_index == -1 {
        return (String::new(), 0, 0);
    }
    let seq = norm_seq.as_bytes();
    let ai = antiloop_index.max(0) as usize;
    let ae = antiloop_end.max(0) as usize;
    if ai >= seq.len() {
        return (String::new(), 0, 0);
    }
    let end = (ae + 1).min(seq.len());
    if end <= ai {
        return (String::new(), 0, 0);
    }
    let antiloop_seq = &norm_seq[ai..end];

    let caps = match INTRON_RE.captures(antiloop_seq) {
        Some(c) => c,
        None => return (String::new(), 0, 0),
    };
    let intron = caps.get(2).unwrap().as_str().to_string();

    // istart = index(substr(seq, 0, antiloop_end+1), intron) + 1  (1-based)
    let search_region = &norm_seq[0..(ae + 1).min(seq.len())];
    let istart = match search_region.find(&intron) {
        Some(p) => p + 1,
        None => return (String::new(), 0, 0),
    };
    let iend = intron.len() + istart - 1;
    (intron, istart as i32, iend as i32)
}

// ============================================================================
// §2.2  Genetic code — trans_map (revcomp of codon table), get_tRNA_type
// ============================================================================

/// Default genetic-code codon table (`GeneticCode.pm` `__DATA__`), bacterial /
/// standard. Each entry: (codon-pattern-with-ambig, isotype).
const CODON_TABLE: &[(&str, &str)] = &[
    ("GCN", "Ala"),
    ("TGY", "Cys"),
    ("GAY", "Asp"),
    ("GAR", "Glu"),
    ("TTY", "Phe"),
    ("GGN", "Gly"),
    ("CAY", "His"),
    ("ATH", "Ile"),
    ("AAR", "Lys"),
    ("TTR", "Leu"),
    ("CTN", "Leu"),
    ("ATG", "Met"),
    ("AAY", "Asn"),
    ("CCN", "Pro"),
    ("CAR", "Gln"),
    ("AGR", "Arg"),
    ("CGN", "Arg"),
    ("AGY", "Ser"),
    ("TCN", "Ser"),
    ("ACN", "Thr"),
    ("GTN", "Val"),
    ("TGG", "Trp"),
    ("TAY", "Tyr"),
    ("TAR", "Sup"),
    ("TGA", "SeC"),
];

/// Reverse-complement including IUPAC ambiguity, port of `rev_comp_seq`.
fn rev_comp_ambig(s: &str) -> String {
    s.to_uppercase()
        .chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' | 'U' => 'A',
            'G' => 'C',
            'C' => 'G',
            'N' => 'N',
            'R' => 'Y',
            'Y' => 'R',
            'W' => 'W',
            'S' => 'S',
            'M' => 'K',
            'K' => 'M',
            'B' => 'V',
            'V' => 'B',
            'D' => 'H',
            'H' => 'D',
            other => other,
        })
        .collect()
}

/// IUPAC ambiguity expansion, port of `expand_ambig` (+ `expand2`/`expand3`).
fn expand_ambig(ac: &str) -> Vec<String> {
    fn subs(c: char) -> Option<&'static [char]> {
        match c {
            'N' => Some(&['A', 'C', 'G', 'T']),
            'Y' => Some(&['C', 'T']),
            'R' => Some(&['A', 'G']),
            'W' => Some(&['A', 'T']),
            'S' => Some(&['C', 'G']),
            'M' => Some(&['A', 'C']),
            'K' => Some(&['G', 'T']),
            'V' => Some(&['A', 'C', 'G']),
            'B' => Some(&['C', 'G', 'T']),
            'H' => Some(&['A', 'C', 'T']),
            'D' => Some(&['A', 'G', 'T']),
            _ => None,
        }
    }
    let mut res: Vec<String> = vec![String::new()];
    for ch in ac.chars() {
        match subs(ch) {
            Some(list) => {
                let mut next = Vec::with_capacity(res.len() * list.len());
                for prefix in &res {
                    for &b in list {
                        let mut s = prefix.clone();
                        s.push(b);
                        next.push(s);
                    }
                }
                res = next;
            }
            None => {
                for s in res.iter_mut() {
                    s.push(ch);
                }
            }
        }
    }
    res
}

/// The anticodon->isotype translation map, built from the reverse-complement of
/// the codon table then ambiguity-expanded (`GeneticCode.pm:190-217`).
/// Uses `BTreeMap` (sorted keys) so that, on any expansion collision, the
/// later-in-sorted-order entry wins — matching Perl's `sort keys` iteration.
static TRANS_MAP: Lazy<BTreeMap<String, String>> = Lazy::new(|| {
    // ambig_trans_map keyed by revcomp(codon); iterated in sorted key order.
    let mut ambig: BTreeMap<String, String> = BTreeMap::new();
    for &(codon, aa) in CODON_TABLE {
        ambig.insert(rev_comp_ambig(codon), aa.to_string());
    }
    let mut trans: BTreeMap<String, String> = BTreeMap::new();
    for (key, aa) in ambig.iter() {
        for expanded in expand_ambig(key) {
            trans.insert(expanded, aa.clone());
        }
    }
    trans
});

/// Port of `GeneticCode.pm::get_tRNA_type` (:267).
///
/// `cm_model_name` is the CM's model name (`"Domain"`, `"SeC"`, ...). `is_sec_cm`
/// short-circuits to `"SeC"` (the Perl `Pselc`/`Eselc` file-path check).
pub fn get_trna_type(anticodon: &str, cm_model_name: &str, is_sec_cm: bool, cove_mode: bool) -> String {
    get_trna_type_with_map(&TRANS_MAP, anticodon, cm_model_name, is_sec_cm, cove_mode)
}

/// Shared body of `GeneticCode.pm::get_tRNA_type` (:267), parameterized on the
/// active anticodon->isotype map (standard vs. an alt-gcode-overridden map).
fn get_trna_type_with_map(
    map: &BTreeMap<String, String>,
    anticodon: &str,
    cm_model_name: &str,
    is_sec_cm: bool,
    cove_mode: bool,
) -> String {
    if anticodon == UNDEF_ANTICODON {
        return UNDEF_ISOTYPE.to_string();
    }
    if is_sec_cm {
        return "SeC".to_string();
    }
    let mut prev: Option<String> = None; // None models Perl 'INIT'
    let mut typ = UNDEF_ISOTYPE.to_string();
    for exp in expand_ambig(&anticodon.to_uppercase()) {
        typ = map
            .get(&exp)
            .cloned()
            .unwrap_or_else(|| UNDEF_ISOTYPE.to_string());
        if typ == "SeC" && cm_model_name != "SeC" && !cove_mode {
            typ = "Sup".to_string();
        }
        if let Some(p) = &prev {
            if &typ != p {
                return UNDEF_ISOTYPE.to_string();
            }
        }
        prev = Some(typ.clone());
    }
    typ
}

/// The vertebrate-mitochondrial anticodon->isotype map: the standard [`TRANS_MAP`]
/// overridden by `gcode.vertmito` (GeneticCode.pm::read_transl_table:219-264, driver
/// :1189 loads gc_vert_mito with alt_gcode for BOTH `-M vert` and `-M mammal`). Each
/// override key is `expand_ambig(rev_comp_seq(codon))`:
///   TGA->Trp => anticodon TCA->Trp   (standard TCA->SeC is replaced)
///   ATA->Met => anticodon TAT->Met
///   AGR->Stp => anticodons TCT,CCT->Stp
static MITO_TRANS_MAP: Lazy<BTreeMap<String, String>> = Lazy::new(|| {
    let mut m = TRANS_MAP.clone();
    m.insert("TCA".to_string(), "Trp".to_string());
    m.insert("TAT".to_string(), "Met".to_string());
    m.insert("TCT".to_string(), "Stp".to_string());
    m.insert("CCT".to_string(), "Stp".to_string());
    m
});

/// `get_tRNA_type` under the vertebrate-mitochondrial genetic code (`-M` mode).
pub fn get_mito_trna_type(anticodon: &str, cm_model_name: &str) -> String {
    get_trna_type_with_map(&MITO_TRANS_MAP, anticodon, cm_model_name, false, false)
}

// ============================================================================
// Top-level decode
// ============================================================================

/// Decode anticodon, isotype and intron for one cmsearch tRNA hit.
///
/// `strand`, `trna_start`, `trna_end` are the genomic tRNA bounds (1-based,
/// `trna_start`/`trna_end` as stored on the tRNA record: for `+` start<end,
/// for `-` start<end too — see spec §2.4 which uses `trna.start()`/`trna.end()`
/// as the ascending genomic bounds). `cm_model_name` selects the SeC->Sup path.
pub fn decode_trna_properties(
    ali: &AliDisplay,
    cm_model_name: &str,
    is_sec_cm: bool,
    cove_mode: bool,
    strand: Strand,
    trna_start: i64,
    trna_end: i64,
) -> DecodedTrna {
    let (norm_ss, norm_seq) = format_cmsearch_output(&ali.ss_cons, &ali.aseq, &ali.nc);
    let (anticodon, antiloop_index, antiloop_end, ac_pos) = find_anticodon(&norm_seq, &norm_ss);

    if anticodon == UNDEF_ANTICODON {
        return DecodedTrna {
            anticodon,
            isotype: UNDEF_ISOTYPE.to_string(),
            antiloop_index,
            antiloop_end,
            ac_pos,
            intron: None,
            norm_seq,
            norm_ss,
        };
    }

    let isotype = get_trna_type(&anticodon, cm_model_name, is_sec_cm, cove_mode);

    // Post-check (CM.pm:1332): SeC isotype requires TCA anticodon.
    let (anticodon, isotype) = if isotype == "SeC" && anticodon != "TCA" {
        (UNDEF_ANTICODON.to_string(), UNDEF_ISOTYPE.to_string())
    } else {
        (anticodon, isotype)
    };

    let intron = {
        let (iseq, istart, iend) = find_intron(&norm_seq, antiloop_index, antiloop_end);
        if iseq.is_empty() {
            None
        } else {
            // Genomic coords (spec §2.4).
            let (gstart, gend) = match strand {
                Strand::Minus => (
                    trna_end - iend as i64 + 1,
                    trna_end - istart as i64 + 1,
                ),
                _ => (
                    istart as i64 + trna_start - 1,
                    iend as i64 + trna_start - 1,
                ),
            };
            Some(Intron {
                rel_start: istart,
                rel_end: iend,
                start: gstart,
                end: gend,
                intron_type: "CI".to_string(),
                seq: iseq,
            })
        }
    };

    DecodedTrna {
        anticodon,
        isotype,
        antiloop_index,
        antiloop_end,
        ac_pos,
        intron,
        norm_seq,
        norm_ss,
    }
}

// ============================================================================
// Mitochondrial decode — CM.pm::find_mito_anticodon (:799) +
// decode_mito_tRNA_properties (:1482)
// ============================================================================

/// The four mito anticodon-stem-loop regexes from `CM.pm::find_mito_anticodon`
/// (:817/824/831/838), operating on the normalized ss. These differ from the
/// nuclear `AC_RE`: the mito set adds "No D-arm" / "No T-arm" tolerant patterns
/// and a more permissive final branch (no required `>` after the D-stem).
static MITO_AC_RE_NODARM_SER: Lazy<Regex> =
    Lazy::new(|| Regex::new(r"^([>.]+)>([.]{4,})<+.+[>.]+<[<.]+").unwrap());
static MITO_AC_RE_NOTARM: Lazy<Regex> =
    Lazy::new(|| Regex::new(r"^([>.]+<[<.]+[>.]*)>([.]{4,})<[<.]+[.]{4,}<[<.]+$").unwrap());
static MITO_AC_RE_NODARM: Lazy<Regex> =
    Lazy::new(|| Regex::new(r"^([>.]+[.]{4,}[>.]+)>([.]{4,})<[<.]+\.+[>.]+<[<.]+$").unwrap());
static MITO_AC_RE_STD: Lazy<Regex> =
    Lazy::new(|| Regex::new(r"^([>.]+<[<.]+[>.]*)>([.]{4,})<+.+[>.]+<[<.]+").unwrap());

/// Result of [`find_mito_anticodon`], mirroring the Perl side-effects.
#[derive(Debug, Clone)]
pub struct MitoAnticodon {
    pub anticodon: String,
    pub antiloop_index: i32,
    pub antiloop_end: i32,
    /// 1-based position of the anticodon's first base (`ac_index+antiloop_index+1`).
    pub ac_pos: i32,
    /// `trna->note()` side-effect ("No D-arm" / "No T-arm" / "").
    pub note: String,
    /// `trna->category()` side-effect ("undetermined_ac" / "mito_ac_mislocation" / "").
    pub category: String,
}

/// True if the CM model key is consistent with a decoded isotype (CM.pm:872/920).
fn mito_model_iso_match(model: &str, isotype: &str) -> bool {
    model == isotype
        || (model == "SerGCT" && isotype == "Ser")
        || (model == "SerTGA" && isotype == "Ser")
        || (model == "LeuTAG" && isotype == "Leu")
        || (model == "LeuTAA" && isotype == "Leu")
        || (model == "Cys_NoDarm" && isotype == "Cys")
}

/// Perl `substr($s, $off, $len)` on a byte string, clamped like Perl (offset past
/// end => empty; len past end => truncated). `off` must be >= 0.
fn perl_substr(s: &[u8], off: usize, len: usize) -> String {
    if off >= s.len() {
        return String::new();
    }
    let end = (off + len).min(s.len());
    String::from_utf8_lossy(&s[off..end]).to_string()
}

/// Port of `CM.pm::find_mito_anticodon` (:799).
///
/// `norm_seq`/`norm_ss` are the normalized (post-`format_cmsearch_output`) strings;
/// `model` is the CM key (e.g. `"Pro"`, `"SerGCT"`, `"Cys_NoDarm"`). Returns the
/// decoded anticodon plus the note/category side-effects the Perl sets on the tRNA.
pub fn find_mito_anticodon(norm_seq: &str, norm_ss: &str, model: &str) -> MitoAnticodon {
    let seq = norm_seq.as_bytes();
    let mut note = String::new();
    let mut antiloop_index: usize = 0;
    let mut antiloop_len: usize = 0;

    // Match one of the four stem-loop patterns, in priority order (CM.pm:817-843).
    if model == "SerGCT" || model == "Cys_NoDarm" {
        if let Some(c) = MITO_AC_RE_NODARM_SER.captures(norm_ss) {
            antiloop_index = c.get(1).unwrap().as_str().len() + 1;
            antiloop_len = c.get(2).unwrap().as_str().len();
            note = "No D-arm".to_string();
        }
    }
    if antiloop_index == 0 {
        if let Some(c) = MITO_AC_RE_NOTARM.captures(norm_ss) {
            antiloop_index = c.get(1).unwrap().as_str().len() + 1;
            antiloop_len = c.get(2).unwrap().as_str().len();
            note = "No T-arm".to_string();
        } else if let Some(c) = MITO_AC_RE_NODARM.captures(norm_ss) {
            antiloop_index = c.get(1).unwrap().as_str().len() + 1;
            antiloop_len = c.get(2).unwrap().as_str().len();
            note = "No D-arm".to_string();
        } else if let Some(c) = MITO_AC_RE_STD.captures(norm_ss) {
            antiloop_index = c.get(1).unwrap().as_str().len() + 1;
            antiloop_len = c.get(2).unwrap().as_str().len();
        }
    }

    let undet = || MitoAnticodon {
        anticodon: UNDEF_ANTICODON.to_string(),
        antiloop_index: -1,
        antiloop_end: -1,
        ac_pos: -1,
        note: note.clone(),
        category: "undetermined_ac".to_string(),
    };

    if antiloop_index == 0 || antiloop_len == 0 {
        return undet();
    }

    let antiloop_end = antiloop_index + antiloop_len - 1;
    // antiloop = substr(seq, antiloop_index, antiloop_len); strip '-' gaps; uc.
    let raw_loop = perl_substr(seq, antiloop_index, antiloop_len);
    let antiloop: Vec<u8> = raw_loop
        .bytes()
        .filter(|&b| b != b'-')
        .map(|b| b.to_ascii_uppercase())
        .collect();

    let mut category = String::new();
    let anticodon: String;
    let ac_index: usize;
    let mut verify_ac: String;

    if antiloop.len() < 5 {
        return undet();
    } else if antiloop.len() % 2 == 0 {
        // Even-length loop: search for a model-consistent triplet (CM.pm:863-897).
        let n = antiloop.len();
        let mut found = false;
        let mut fi: i64 = ((n as i64 - 3) / 2) as i64;
        let mut j: i64 = 0;
        let mut cur_ai: usize = 0;
        let mut cur_ac = String::new();
        while fi <= (n as i64 - 3) && fi >= 0 {
            cur_ai = fi as usize;
            cur_ac = String::from_utf8_lossy(&antiloop[cur_ai..cur_ai + 3]).to_string();
            let iso = get_mito_trna_type(&cur_ac, model);
            if mito_model_iso_match(model, &iso) {
                category = "mito_ac_mislocation".to_string();
                found = true;
                break;
            }
            j = j.abs();
            j += 1;
            if j % 2 == 0 {
                j = -j;
            }
            fi += j;
        }
        if !found {
            return undet();
        }
        ac_index = cur_ai;
        anticodon = cur_ac;
        verify_ac = perl_substr(seq, ac_index + antiloop_index, 3).to_ascii_uppercase();
    } else {
        // Odd-length loop: centered triplet, with a consistency re-search
        // (CM.pm:898-940).
        let n = antiloop.len();
        let mut cur_ai = (n - 3) / 2;
        let mut cur_ac = String::from_utf8_lossy(&antiloop[cur_ai..cur_ai + 3]).to_string();
        verify_ac = perl_substr(seq, cur_ai + antiloop_index, 3).to_ascii_uppercase();
        let iso = get_mito_trna_type(&cur_ac, model);
        let model_iso: &str = if model.len() > 3 { &model[0..3] } else { model };
        if iso != model_iso {
            verify_ac = String::new();
            let mut found = false;
            let mut fi: i64 = ((n as i64 - 3) / 2) - 1;
            let mut j: i64 = 1;
            while fi <= (n as i64 - 3) && fi >= 0 {
                cur_ai = fi as usize;
                cur_ac = String::from_utf8_lossy(&antiloop[cur_ai..cur_ai + 3]).to_string();
                let iso2 = get_mito_trna_type(&cur_ac, model);
                if mito_model_iso_match(model, &iso2) {
                    category = "mito_ac_mislocation".to_string();
                    found = true;
                    break;
                }
                j = j.abs();
                j += 1;
                if j % 2 == 1 {
                    j = -j;
                }
                fi += j;
            }
            if found {
                verify_ac = perl_substr(seq, cur_ai + antiloop_index, 3).to_ascii_uppercase();
            }
        }
        ac_index = cur_ai;
        anticodon = cur_ac;
    }

    // verify_ac must equal the loop-derived anticodon (CM.pm:946).
    if verify_ac != anticodon {
        return MitoAnticodon {
            anticodon: UNDEF_ANTICODON.to_string(),
            antiloop_index: -1,
            antiloop_end: -1,
            ac_pos: -1,
            note,
            category: "undetermined_ac".to_string(),
        };
    }

    MitoAnticodon {
        anticodon,
        antiloop_index: antiloop_index as i32,
        antiloop_end: antiloop_end as i32,
        ac_pos: (ac_index + antiloop_index + 1) as i32,
        note,
        category,
    }
}

/// The result of the mito property decode (CM.pm:decode_mito_tRNA_properties:1482).
///
/// Note: `decode_mito_tRNA_properties` computes intron bounds into LOCAL variables
/// (CM.pm:1517-1531) that are never stored on the tRNA, so mito tRNAs carry NO
/// intron and the `.out` Intron Bounds are always `0  0`. We therefore do not
/// expose an intron here.
#[derive(Debug, Clone)]
pub struct DecodedMitoTrna {
    /// Anticodon (uppercase T-form) or `"NNN"`.
    pub anticodon: String,
    /// Output Type = the model-key isotype (first 3 chars). CM.pm:1595.
    pub isotype: String,
    /// `trna->note()` value after the decode's "(...)" prepends (CM.pm:1560-1584).
    pub note: String,
    /// `trna->category()` value (used only by the `--detail` mito Note column,
    /// ScanResult.pm:843-857). Empty when consistent.
    pub category: String,
    pub norm_seq: String,
    pub norm_ss: String,
}

/// `vert_mito_aa_list` (GeneticCode.pm:114): anticodon -> isotype for the
/// vertebrate-mito expected-anticodon check (CM.pm:1534/1587).
fn vert_mito_type(ac: &str) -> &'static str {
    match ac {
        "TGC" => "Ala", "TCC" => "Gly", "TGG" => "Pro", "TGT" => "Thr", "TAC" => "Val",
        "TGA" => "Ser", "GCT" => "Ser", "TCG" => "Arg", "TAG" => "Leu", "TAA" => "Leu",
        "GAA" => "Phe", "GTT" => "Asn", "TTT" => "Lys", "GTC" => "Asp", "TTC" => "Glu",
        "GTG" => "His", "TTG" => "Gln", "GTA" => "Tyr", "GAT" => "Ile", "TAT" => "Met",
        "CAT" => "Met", "GCA" => "Cys", "TCA" => "Trp", "GCC" => "Asp",
        _ => "",
    }
}

/// Port of `CM.pm::decode_mito_tRNA_properties` (:1482).
///
/// `ali` is the normalized-input alidisplay for the winning model; `model` is the
/// CM key. Intron bounds are computed-but-discarded in the Perl, so no intron is
/// returned (see [`DecodedMitoTrna`]).
#[allow(non_snake_case)] // mirrors the C `decode_mito_tRNA_properties` sub name
pub fn decode_mito_tRNA_properties(ali: &AliDisplay, model: &str) -> DecodedMitoTrna {
    let (norm_ss, norm_seq) = format_cmsearch_output(&ali.ss_cons, &ali.aseq, &ali.nc);
    let ac = find_mito_anticodon(&norm_seq, &norm_ss, model);
    let mut note = ac.note.clone();
    let mut category = ac.category.clone();

    // check for problem parsing anticodon loop (CM.pm:1501): undef => intron off.
    let anticodon = if ac.anticodon == UNDEF_ANTICODON {
        UNDEF_ANTICODON.to_string()
    } else {
        ac.anticodon.clone()
    };

    // isotype = get_tRNA_type(anticodon, model) (CM.pm:1533).
    let isotype_decoded = get_mito_trna_type(&anticodon, model);
    let vert_iso = vert_mito_type(&anticodon);

    // model_iso / model_ac extraction (CM.pm:1536-1554).
    let mut model_iso = model.to_string();
    let mut model_ac = String::new();
    if model.len() > 3 {
        model_iso = model[0..3].to_string();
        if let Some(us) = model.find('_') {
            let temp = &model[0..us];
            if temp.len() > 3 {
                model_ac = temp[3..].to_string();
            }
        } else {
            model_ac = model[3..].to_string();
        }
    }

    // Consistency notes (CM.pm:1555-1593). Perl prepends "(...)" only when
    // category is still empty at that check.
    if model_iso != isotype_decoded {
        if category.is_empty() {
            category = "mito_inconsistent_isotype".to_string();
            note = if !note.is_empty() {
                format!("({}); {}", isotype_decoded, note)
            } else {
                format!("({})", isotype_decoded)
            };
        }
    }
    if !model_ac.is_empty() && model_ac != anticodon {
        if category.is_empty() {
            category = "mito_inconsistent_ac".to_string();
            note = if !note.is_empty() {
                format!("({}); {}", model_ac, note)
            } else {
                format!("({})", model_ac)
            };
        }
    }
    if vert_iso.is_empty() && category.is_empty() {
        // category = "mito_unexpected_ac" (no note side-effect).
    }

    DecodedMitoTrna {
        anticodon,
        isotype: model_iso,
        note,
        category,
        norm_seq,
        norm_ss,
    }
}

// ============================================================================
// Tests — real C-cmsearch (INFERNAL 1.1.5) alignment fixtures
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    /// One fixture: (ss_cons, aseq, nc) exactly as parsed from a
    /// `cmsearch -g --nohmm --toponly --notrunc TRNAinf-bact.cm` alidisplay.
    struct Fx {
        ss: &'static str,
        seq: &'static str,
        nc: &'static str,
        strand: Strand,
        start: i64,
        end: i64,
        exp_ac: &'static str,
        exp_type: &'static str,
        exp_intron: Option<(i32, i32, i64, i64)>, // (istart,iend,gstart,gend)
    }

    fn run(fx: &Fx) -> DecodedTrna {
        let ali = AliDisplay {
            aseq: fx.seq.to_string(),
            ss_cons: fx.ss.to_string(),
            nc: fx.nc.to_string(),
            model: String::new(),
        };
        decode_trna_properties(&ali, "Domain", false, false, fx.strand, fx.start, fx.end)
    }

    fn fixtures() -> Vec<(&'static str, Fx)> {
        vec![
            (
                "E1_1",
                Fx {
                    ss: "(((((((,,<<<<_______.._>>>>,<<<<<__...................................._____>>>>>,,<<<<<<<____>>>>>>>,,<<<<<_______>>>>>)))))))::::",
                    seq: "GCACGGAUGGCCGAGU-GGUcuAAGGCGCCAGACUcaagcgaaaugcuugccucaugcucgaggucgacugGGUGUUCUGGU-ACUCGU------AUGGGUGCGUGGGUUCGAAUCCCACUUCGUGCA---",
                    nc:  "                                                                                                                                       ",
                    strand: Strand::Plus, start: 12619, end: 12738,
                    exp_ac: "NNN", exp_type: "Undet", exp_intron: None,
                },
            ),
            (
                "E1_2",
                Fx {
                    ss: "(((((((,,<<<<_______._>>>>,<<<<<_______>>>>>,,<<<<<<<____>>>>>>>,,<<<<<_______>>>>>)))))))::::",
                    seq: "GCAGUCAUGUCCGAGU-GGUuAAGGAGAUUGACUAGAAAUCAAUU-GGGC---UCU----GCCCGCGUAGGUUCGAAUCCUGCUGACUGCG---",
                    nc:  "                                                                                                               ",
                    strand: Strand::Plus, start: 19480, end: 19561,
                    exp_ac: "AGA", exp_type: "Ser", exp_intron: None,
                },
            ),
            (
                "E1_3",
                Fx {
                    ss: "(((((((,,<<<<________>>>>,<<<<<_______>>>>>,,<<<<<<<____>>>>>>>,.,<<<<<_______>>>>>)))))))::::",
                    seq: "GCCUCGAUAGCUCAGUUGGGAGAGCGUACGACUGAAGAUCGUAAG------------------GuCACCAGUUCGAUCCUGGUUCGGGGCA---",
                    nc:  "                                                                                                               ",
                    strand: Strand::Plus, start: 26367, end: 26439,
                    exp_ac: "GAA", exp_type: "Phe", exp_intron: None,
                },
            ),
            (
                "E1_4",
                Fx {
                    ss: "(((((((,,<<<<________>>>>,<<<<<_______>>>>>,,<<<<<<<____>>>>>>>,.,<<<<<_______>>>>>)))))))::::",
                    seq: "GCCUCGAUAGCUCAGUUGGGAGAGCGUACGACUGAAGAUCGUAAG------------------GuCACCAGUUCGAUCCUGGUUCGGGGCA---",
                    nc:  "                                                                                                               ",
                    strand: Strand::Minus, start: 26920, end: 26992,
                    exp_ac: "GAA", exp_type: "Phe", exp_intron: None,
                },
            ),
            (
                "E1_5",
                Fx {
                    ss: "(((((((,,<<<<________>>>>,<<<<<_______>>>>>,,<<<<<<<____>>>>>>>,.,<<<<<_______>>>>>)))))))::::",
                    seq: "GGCCGGAUGGUCUAGA-GGUAUGAUUCUCGCUUCGGGUGCGAGAG------------------GuCCCGGGUUCGAUUCCCGGUUCGGCCC---",
                    nc:  "                             v        v                                                                        ",
                    strand: Strand::Minus, start: 23694, end: 23765,
                    exp_ac: "CGG", exp_type: "Pro", exp_intron: None,
                },
            ),
            (
                "MySeq1",
                Fx {
                    ss: "(((((((,,<<<<_______._>>>>,<<<<<_______>>>>>,,<<<<<<<____>>>>>>>,.,<<<<<_______>>>>>)))))))::::",
                    seq: "GGCCCUAUAGCUCAGG-GGUuAGAGCACUGGUCUUGUAAACCAGGG------------------GuCGCGAGUUCAAAUCUCGCUGGGGCCU---",
                    nc:  "                                                                                                                ",
                    strand: Strand::Plus, start: 13, end: 85,
                    exp_ac: "TGT", exp_type: "Thr", exp_intron: None,
                },
            ),
            (
                "MySeq2",
                Fx {
                    ss: "(((((((,,<<<<_______.._>>>>,<<<<<_______>>>>>,,<<<<<<<____>>>>>>>,.,<<<<<_______>>>>>)))))))::::",
                    seq: "GUCUCUGUGGCGCAAU-GGAcgAGCGCGCUGGACUUCUAAUCCAGAG------------------GuUCUGGGUUCGAGUCCCGGCAGAGAUG---",
                    nc:  "                                                                                                                 ",
                    strand: Strand::Plus, start: 6, end: 79,
                    exp_ac: "TCT", exp_type: "Arg", exp_intron: None,
                },
            ),
            (
                "MySeq3",
                Fx {
                    ss: "(((((((,,<<<<_______._>>>>,<<<<<______..................._>>>>>,,<<<<<<<____>>>>>>>,,<<<<<_______>>>>>)))))))::::",
                    seq: "GGCACUAUGGCCGAGU-GGUuAAGGCGAGAGACUCGAAuggaauaaaaaguucggcuAUCUCUU-GGGC---UCU----GCCCGCGCUGGUUCAAAUCCUGCUGGUGUCG---",
                    nc:  "                                                                                                         v                             v             ",
                    strand: Strand::Plus, start: 14, end: 114,
                    exp_ac: "CGA", exp_type: "Ser", exp_intron: Some((38, 56, 51, 69)),
                },
            ),
            (
                "MySeq4",
                Fx {
                    ss: "(((((((,,<<<<_______.._>>>>,<<<<<______._>>>>>,,<<<<<<<____>>>>>>>,,<<<<<_______>>>>>)))))))::::",
                    seq: "GGAGAGAUGGCCGAGC-GGUccAAGGCGCUGGUUUAAGGcAACCAGUAGCUUC--------GGGGG-CGUGGGUUCGAAUCCCACUCUCUUCA---",
                    nc:  "                                                                                                                 ",
                    strand: Strand::Plus, start: 6, end: 88,
                    exp_ac: "AAG", exp_type: "Leu", exp_intron: None,
                },
            ),
            (
                "MySeq5",
                Fx {
                    ss: "(((((((,,..<<<<_______._>>>>,<<<<<_______>>>>>,,<<<<<..<<____>>>>>>>,,<<<<<_______>>>.>>)))))))::::",
                    seq: "GCCCGGAUGauCCUCAGU-GGUcUGGGGUGCAGGCUUCAAACCUGUAGCUGUCuaG------CGACAGAG--UGGUUCAAUUCCAcCUUUCGGGCG---",
                    nc:  "                                                                                       vv              vv           ",
                    strand: Strand::Plus, start: 3, end: 89,
                    exp_ac: "TCA", exp_type: "Sup", exp_intron: None,
                },
            ),
            (
                "MySeq6",
                Fx {
                    ss: "(((((((,,<<<<_______.._>>>>,<<<<<_______>>>>>,,<<<<<<<____>>>>>>>,.,<<<<<_______>>>>>)))))))::::",
                    seq: "GACACGGUGGCCGAGU-GGUuuAAGGCAUGAGACACUUGAUCUCAAACGGU---UCUA---ACCGAaCGCAGGUUCGAAUCCUGCCCGUGUCA---",
                    nc:  "                                                                                                                 ",
                    strand: Strand::Plus, start: 7, end: 92,
                    exp_ac: "CTT", exp_type: "Lys", exp_intron: None,
                },
            ),
        ]
    }

    #[test]
    fn test_all_fixtures_decode() {
        let mut failures = Vec::new();
        for (name, fx) in fixtures() {
            let d = run(&fx);
            let mut errs = Vec::new();
            if d.anticodon != fx.exp_ac {
                errs.push(format!("anticodon {} != {}", d.anticodon, fx.exp_ac));
            }
            if d.isotype != fx.exp_type {
                errs.push(format!("type {} != {}", d.isotype, fx.exp_type));
            }
            match (&d.intron, fx.exp_intron) {
                (None, None) => {}
                (Some(i), Some((is, ie, gs, ge))) => {
                    if i.rel_start != is || i.rel_end != ie || i.start != gs || i.end != ge {
                        errs.push(format!(
                            "intron ({},{},{},{}) != ({},{},{},{})",
                            i.rel_start, i.rel_end, i.start, i.end, is, ie, gs, ge
                        ));
                    }
                }
                (got, exp) => errs.push(format!("intron {:?} != {:?}", got.is_some(), exp.is_some())),
            }
            if !errs.is_empty() {
                failures.push(format!("{}: {}", name, errs.join("; ")));
            } else {
                eprintln!(
                    "{:8} ac={:4} type={:6} intron={:?}  OK",
                    name,
                    d.anticodon,
                    d.isotype,
                    d.intron.as_ref().map(|i| (i.rel_start, i.rel_end, i.start, i.end))
                );
            }
        }
        assert!(failures.is_empty(), "decode failures:\n{}", failures.join("\n"));
    }

    #[test]
    fn test_ser_fixture_normalization() {
        // The concrete verified Ser fixture (Example1 #2).
        let (ss, seq) = format_cmsearch_output(
            "(((((((,,<<<<_______._>>>>,<<<<<_______>>>>>,,<<<<<<<____>>>>>>>,,<<<<<_______>>>>>)))))))::::",
            "GCAGUCAUGUCCGAGU-GGUuAAGGAGAUUGACUAGAAAUCAAUU-GGGC---UCU----GCCCGCGUAGGUUCGAAUCCUGCUGACUGCG---",
            "                                                                                                               ",
        );
        assert_eq!(ss.len(), seq.len());
        let (ac, ai, _ae, _acp) = find_anticodon(&seq, &ss);
        assert_eq!(ac, "AGA");
        assert!(ai > 0);
    }

    #[test]
    fn test_trans_map_key_examples() {
        // Sanity: a few revcomp-built anticodon -> isotype entries.
        assert_eq!(TRANS_MAP.get("AGA").map(String::as_str), Some("Ser"));
        assert_eq!(TRANS_MAP.get("GAA").map(String::as_str), Some("Phe"));
        assert_eq!(TRANS_MAP.get("CGG").map(String::as_str), Some("Pro"));
        assert_eq!(TRANS_MAP.get("CTT").map(String::as_str), Some("Lys"));
        assert_eq!(TRANS_MAP.get("CAT").map(String::as_str), Some("Met"));
        assert_eq!(TRANS_MAP.get("TCA").map(String::as_str), Some("SeC"));
    }

    #[test]
    fn test_sec_to_sup_downgrade() {
        // TCA under a non-SeC domain model -> Sup.
        assert_eq!(get_trna_type("TCA", "Domain", false, false), "Sup");
        // TCA under the SeC CM path -> SeC.
        assert_eq!(get_trna_type("TCA", "SeC", true, false), "SeC");
    }
}
