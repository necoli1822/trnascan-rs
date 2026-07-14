//! Statistics module for tRNAscan-SE
//!
//! This module tracks and reports statistics for tRNA scanning runs,
//! including counts, timings, scores, and detailed breakdowns by isotype/anticodon.
//!
//! Ported from tRNAscanSE::Stats.pm

use std::collections::HashMap;
use std::io::{self, Write};
use std::time::Instant;

// ============================================================================
// Main Statistics Struct
// ============================================================================

/// Comprehensive statistics for a tRNAscan-SE run
#[derive(Debug, Clone)]
pub struct ScanStats {
    // === File handling ===
    file_name: String,

    // === Timing (first pass) ===
    fp_start_time: Option<Instant>,
    fp_end_time: Option<Instant>,
    fp_script_cpu: f64,
    fp_scan_cpu: f64,

    // === Timing (second pass) ===
    sp_start_time: Option<Instant>,
    sp_end_time: Option<Instant>,
    sp_script_cpu: f64,
    sp_scan_cpu: f64,

    // === Sequence counts ===
    pub numscanned: usize,           // Total sequences scanned
    pub seqs_hit: usize,             // Sequences with at least one hit
    pub trnatotal: usize,            // Total tRNAs found by first pass

    // === Base counts ===
    pub first_pass_base_ct: usize,   // Total bases in all sequences
    pub fpass_trna_base_ct: usize,   // Bases in tRNAs (first pass)
    pub fpos_base_ct: usize,         // Bases in false positives
    pub secpass_base_ct: usize,      // Bases scanned in second pass
    pub coves_base_ct: usize,        // Bases scanned by coves/infernal
    pub total_secpass_ct: usize,     // Total confirmed by second pass
}

impl ScanStats {
    /// Create new statistics object
    pub fn new() -> Self {
        ScanStats {
            file_name: String::new(),
            fp_start_time: None,
            fp_end_time: None,
            fp_script_cpu: 0.0,
            fp_scan_cpu: 0.0,
            sp_start_time: None,
            sp_end_time: None,
            sp_script_cpu: 0.0,
            sp_scan_cpu: 0.0,
            numscanned: 0,
            seqs_hit: 0,
            trnatotal: 0,
            first_pass_base_ct: 0,
            fpass_trna_base_ct: 0,
            fpos_base_ct: 0,
            secpass_base_ct: 0,
            coves_base_ct: 0,
            total_secpass_ct: 0,
        }
    }

    /// Clear all statistics to initial state
    pub fn clear(&mut self) {
        *self = ScanStats::new();
    }

    // ========================================================================
    // File Management
    // ========================================================================

    /// Set output file name
    pub fn set_file_name(&mut self, name: &str) {
        self.file_name = name.to_string();
    }

    /// Get output file name
    pub fn file_name(&self) -> &str {
        &self.file_name
    }

    // ========================================================================
    // Timer Management
    // ========================================================================

    /// Start first-pass timer
    pub fn start_fp_timer(&mut self) {
        self.fp_start_time = Some(Instant::now());
        self.fp_end_time = None;
        self.sp_start_time = None;
        self.sp_end_time = None;
    }

    /// End first-pass timer
    pub fn end_fp_timer(&mut self) {
        if let Some(start) = self.fp_start_time {
            self.fp_end_time = Some(Instant::now());
            // In production, would get CPU times from OS
            // For now, use wall clock time as approximation
            if let Some(end) = self.fp_end_time {
                let elapsed = end.duration_since(start);
                self.fp_script_cpu = elapsed.as_secs_f64();
                self.fp_scan_cpu = elapsed.as_secs_f64();
            }
        }
    }

    /// Start second-pass timer
    pub fn start_sp_timer(&mut self) {
        self.sp_start_time = Some(Instant::now());
        // If first pass never ended, end it now
        if self.fp_end_time.is_none() {
            self.fp_end_time = self.fp_start_time;
        }
    }

    /// End second-pass timer
    pub fn end_sp_timer(&mut self) {
        if let Some(start) = self.sp_start_time {
            self.sp_end_time = Some(Instant::now());
            // In production, would get CPU times from OS
            if let Some(end) = self.sp_end_time {
                let elapsed = end.duration_since(start);
                self.sp_script_cpu = elapsed.as_secs_f64();
                self.sp_scan_cpu = elapsed.as_secs_f64();
            }
        }
    }

    // ========================================================================
    // Sequence Counters
    // ========================================================================

    /// Get sequences hit count
    pub fn seqs_hit(&self) -> usize {
        self.seqs_hit
    }

    /// Set sequences hit count
    pub fn set_seqs_hit(&mut self, count: usize) {
        self.seqs_hit = count;
    }

    /// Increment sequences hit
    pub fn increment_seqs_hit(&mut self, count: usize) {
        self.seqs_hit += count;
    }

    /// Increment sequences hit by 1
    pub fn inc_seqs_hit(&mut self) {
        self.seqs_hit += 1;
    }

    /// Get total sequences scanned
    pub fn numscanned(&self) -> usize {
        self.numscanned
    }

    /// Set total sequences scanned
    pub fn set_numscanned(&mut self, count: usize) {
        self.numscanned = count;
    }

    /// Increment sequences scanned
    pub fn increment_numscanned(&mut self, count: usize) {
        self.numscanned += count;
    }

    /// Increment sequences scanned by 1
    pub fn inc_numscanned(&mut self) {
        self.numscanned += 1;
    }

    /// Get total tRNAs found
    pub fn trnatotal(&self) -> usize {
        self.trnatotal
    }

    /// Set total tRNAs found
    pub fn set_trnatotal(&mut self, count: usize) {
        self.trnatotal = count;
    }

    /// Increment tRNA total
    pub fn increment_trnatotal(&mut self, count: usize) {
        self.trnatotal += count;
    }

    /// Increment tRNA total by 1
    pub fn inc_trnatotal(&mut self) {
        self.trnatotal += 1;
    }

    /// Decrement tRNA total
    pub fn decrement_trnatotal(&mut self, count: usize) {
        self.trnatotal = self.trnatotal.saturating_sub(count);
    }

    /// Decrement tRNA total by 1
    pub fn dec_trnatotal(&mut self) {
        self.trnatotal = self.trnatotal.saturating_sub(1);
    }

    // ========================================================================
    // Base Count Getters/Setters
    // ========================================================================

    /// Get first pass base count
    pub fn first_pass_base_ct(&self) -> usize {
        self.first_pass_base_ct
    }

    /// Set first pass base count
    pub fn set_first_pass_base_ct(&mut self, count: usize) {
        self.first_pass_base_ct = count;
    }

    /// Increment first pass base count
    pub fn increment_first_pass_base_ct(&mut self, count: usize) {
        self.first_pass_base_ct += count;
    }

    /// Increment first pass base count by 1
    pub fn inc_first_pass_base_ct(&mut self) {
        self.first_pass_base_ct += 1;
    }

    /// Get first pass tRNA base count
    pub fn fpass_trna_base_ct(&self) -> usize {
        self.fpass_trna_base_ct
    }

    /// Set first pass tRNA base count
    pub fn set_fpass_trna_base_ct(&mut self, count: usize) {
        self.fpass_trna_base_ct = count;
    }

    /// Increment first pass tRNA base count
    pub fn increment_fpass_trna_base_ct(&mut self, count: usize) {
        self.fpass_trna_base_ct += count;
    }

    /// Increment first pass tRNA base count by 1
    pub fn inc_fpass_trna_base_ct(&mut self) {
        self.fpass_trna_base_ct += 1;
    }

    /// Get false positive base count
    pub fn fpos_base_ct(&self) -> usize {
        self.fpos_base_ct
    }

    /// Set false positive base count
    pub fn set_fpos_base_ct(&mut self, count: usize) {
        self.fpos_base_ct = count;
    }

    /// Increment false positive base count
    pub fn increment_fpos_base_ct(&mut self, count: usize) {
        self.fpos_base_ct += count;
    }

    /// Increment false positive base count by 1
    pub fn inc_fpos_base_ct(&mut self) {
        self.fpos_base_ct += 1;
    }

    /// Get second pass base count
    pub fn secpass_base_ct(&self) -> usize {
        self.secpass_base_ct
    }

    /// Set second pass base count
    pub fn set_secpass_base_ct(&mut self, count: usize) {
        self.secpass_base_ct = count;
    }

    /// Increment second pass base count
    pub fn increment_secpass_base_ct(&mut self, count: usize) {
        self.secpass_base_ct += count;
    }

    /// Increment second pass base count by 1
    pub fn inc_secpass_base_ct(&mut self) {
        self.secpass_base_ct += 1;
    }

    /// Get coves/infernal base count
    pub fn coves_base_ct(&self) -> usize {
        self.coves_base_ct
    }

    /// Set coves/infernal base count
    pub fn set_coves_base_ct(&mut self, count: usize) {
        self.coves_base_ct = count;
    }

    /// Increment coves/infernal base count
    pub fn increment_coves_base_ct(&mut self, count: usize) {
        self.coves_base_ct += count;
    }

    /// Increment coves/infernal base count by 1
    pub fn inc_coves_base_ct(&mut self) {
        self.coves_base_ct += 1;
    }

    /// Get total second pass count
    pub fn total_secpass_ct(&self) -> usize {
        self.total_secpass_ct
    }

    /// Set total second pass count
    pub fn set_total_secpass_ct(&mut self, count: usize) {
        self.total_secpass_ct = count;
    }

    /// Increment total second pass count
    pub fn increment_total_secpass_ct(&mut self, count: usize) {
        self.total_secpass_ct += count;
    }

    /// Increment total second pass count by 1
    pub fn inc_total_secpass_ct(&mut self) {
        self.total_secpass_ct += 1;
    }

    // ========================================================================
    // File I/O
    // ========================================================================

    /// Write a line to the stats file
    pub fn write_line<W: Write>(&self, w: &mut W, line: &str) -> io::Result<()> {
        writeln!(w, "{}", line)
    }

    /// Save first-pass statistics
    pub fn save_firstpass_stats<W: Write>(&self, w: &mut W) -> io::Result<()> {
        writeln!(w, "First-pass Stats:")?;
        writeln!(w, "---------------")?;
        writeln!(w, "Sequences read:         {}", self.numscanned)?;
        writeln!(w, "Seqs w/at least 1 hit:  {}", self.seqs_hit)?;
        writeln!(w, "Bases read:             {} (x2 for both strands)", self.first_pass_base_ct)?;
        writeln!(w, "Bases in tRNAs:         {}", self.fpass_trna_base_ct)?;
        writeln!(w, "tRNAs predicted:        {}", self.trnatotal)?;

        let avg_len = if self.trnatotal > 0 {
            self.fpass_trna_base_ct / self.trnatotal
        } else {
            0
        };
        writeln!(w, "Av. tRNA length:        {}", avg_len)?;

        writeln!(w, "Script CPU time:        {:.2} s", self.fp_script_cpu)?;
        writeln!(w, "Scan CPU time:          {:.2} s", self.fp_scan_cpu)?;

        let scan_speed = if self.fp_scan_cpu > 0.001 {
            (self.first_pass_base_ct * 2) as f64 / self.fp_scan_cpu / 1000.0
        } else {
            0.0
        };
        writeln!(w, "Scan speed:             {:.1} Kbp/sec", scan_speed)?;

        // Get current time for end message. C: `print "\nFirst pass search(es)
        // ended: ",`date`,"\n"` — `date` supplies one trailing newline and the extra
        // "\n" a blank line, so a blank line follows this line (before "Infernal Stats:").
        let now = chrono::Local::now();
        writeln!(w, "\nFirst pass search(es) ended: {}\n", now.format("%a %b %d %H:%M:%S %Y"))?;

        Ok(())
    }

    /// Save the second-pass ("Infernal Stats") block + "Overall scan speed" line, for
    /// the faithful `.stats` writer. Faithful port of the CM-mode / infernal_fp branch
    /// of C Stats.pm `save_final_stats` (label = "Infernal"), EXCLUDING output_summary
    /// (the caller emits the isotype-count summary via `write_faithful_stats`). The
    /// date / CPU-time / scan-speed lines carry runtime values (normalized in parity
    /// checks); the labels, counts, and "% seq scanned" are deterministic.
    pub fn save_secondpass_stats<W: Write>(&self, w: &mut W) -> io::Result<()> {
        let label = "Infernal";
        writeln!(w, "{} Stats:", label)?;
        writeln!(w, "-----------")?;
        // C passes the first-pass hit count as prescan_count; `trnatotal` holds it.
        writeln!(w, "Candidate tRNAs read:     {}", self.trnatotal)?;
        writeln!(w, "{}-confirmed tRNAs:     {}", label, self.total_secpass_ct)?;
        writeln!(w, "Bases scanned by {}:  {}", label, self.secpass_base_ct)?;

        let pct_scanned = ((self.secpass_base_ct as f64
            / (self.first_pass_base_ct * 2).max(1) as f64)
            * 100.0)
            .min(100.0);
        writeln!(w, "% seq scanned by {}:  {:.1} %", label, pct_scanned)?;

        writeln!(w, "Script CPU time:          {:.2} s", self.sp_script_cpu)?;
        writeln!(w, "{} CPU time:            {:.2} s", label, self.sp_scan_cpu)?;

        let scan_speed = self.secpass_base_ct as f64 / self.sp_scan_cpu.max(0.001);
        writeln!(w, "Scan speed:               {:.1} bp/sec", scan_speed)?;

        let now = chrono::Local::now();
        writeln!(
            w,
            "\n{} analysis of tRNAs ended: {}\n",
            label,
            now.format("%a %b %d %H:%M:%S %Y")
        )?;

        let total_time = self.fp_script_cpu + self.fp_scan_cpu + self.sp_script_cpu + self.sp_scan_cpu;
        let total_bases = (self.first_pass_base_ct * 2).max(self.secpass_base_ct);
        let overall_speed = total_bases as f64 / total_time.max(0.001);
        writeln!(w, "Overall scan speed: {:.1} bp/sec", overall_speed)?;

        Ok(())
    }

    /// Save final statistics with summary
    pub fn save_final_stats<W: Write>(
        &self,
        w: &mut W,
        second_pass_label: &str,
        is_cm_mode: bool,
        is_prescan_mode: bool,
        prescan_count: usize,
        tab_results: &[String],
        gc_isotypes: &[&str],
        gc_ac_list: &HashMap<String, Vec<Vec<String>>>,
        show_mismatch: bool,
    ) -> io::Result<()> {
        if is_cm_mode {
            writeln!(w, "{} Stats:", second_pass_label)?;
            writeln!(w, "-----------")?;

            if is_prescan_mode {
                writeln!(w, "Candidate tRNAs read:     {}", prescan_count)?;
            } else {
                writeln!(w, "Sequences read:           {}", self.numscanned)?;
            }

            writeln!(w, "{}-confirmed tRNAs:     {}", second_pass_label, self.total_secpass_ct)?;
            writeln!(w, "Bases scanned by {}:  {}", second_pass_label, self.secpass_base_ct)?;

            let pct_scanned = if self.first_pass_base_ct > 0 {
                ((self.secpass_base_ct as f64 / (self.first_pass_base_ct * 2) as f64) * 100.0).min(100.0)
            } else {
                0.0
            };
            writeln!(w, "% seq scanned by {}:  {:.1} %", second_pass_label, pct_scanned)?;

            writeln!(w, "Script CPU time:          {:.2} s", self.sp_script_cpu)?;
            writeln!(w, "{} CPU time:            {:.2} s", second_pass_label, self.sp_scan_cpu)?;

            let scan_speed = if self.sp_scan_cpu > 0.001 {
                self.secpass_base_ct as f64 / self.sp_scan_cpu
            } else {
                0.0
            };
            writeln!(w, "Scan speed:               {:.1} bp/sec", scan_speed)?;

            let now = chrono::Local::now();
            writeln!(w, "\n{} analysis of tRNAs ended: {}", second_pass_label, now.format("%a %b %d %H:%M:%S %Y"))?;

            if is_prescan_mode {
                writeln!(w, "Summary")?;
                writeln!(w, "--------")?;
            }
        }

        // Overall scan speed
        let total_time = self.fp_script_cpu + self.fp_scan_cpu + self.sp_script_cpu + self.sp_scan_cpu;
        let total_bases = (self.first_pass_base_ct * 2).max(self.secpass_base_ct);
        let overall_speed = if total_time > 0.001 {
            total_bases as f64 / total_time
        } else {
            0.0
        };
        writeln!(w, "Overall scan speed: {:.1} bp/sec", overall_speed)?;

        // Output summary
        self.output_summary(w, tab_results, gc_isotypes, gc_ac_list, show_mismatch)?;

        Ok(())
    }

    /// Output detailed summary of tRNA counts by isotype/anticodon
    fn output_summary<W: Write>(
        &self,
        w: &mut W,
        tab_results: &[String],
        gc_isotypes: &[&str],
        gc_ac_list: &HashMap<String, Vec<Vec<String>>>,
        show_mismatch: bool,
    ) -> io::Result<()> {
        // Parse tab results to count isotypes and anticodons
        let mut iso_counts: HashMap<String, usize> = HashMap::new();
        let iso_cm_counts: HashMap<String, usize> = HashMap::new();
        let mut ac_counts: HashMap<String, HashMap<String, usize>> = HashMap::new();
        let mut intron_ac_counts: HashMap<String, HashMap<String, usize>> = HashMap::new();

        let mut trna_ct = 0;
        let mut selcys_ct = 0;
        let mut pseudo_ct = 0;
        let mut mismatch_ct = 0;
        let mut undet_ct = 0;
        let mut stop_sup_ct = 0;
        let mut intron_ct = 0;

        for line in tab_results {
            if line.is_empty() {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 9 {
                continue;
            }

            let iso = fields[4];
            let ac = fields[5];
            let istart = fields[6].trim();

            // Get note field (last column)
            let note_col = fields.len() - 1;
            let note = fields.get(note_col).unwrap_or(&"");

            // Classify tRNA
            if note.contains("pseudo") {
                pseudo_ct += 1;
                *iso_counts.entry("Pseudo".to_string()).or_insert(0) += 1;
            } else if note.contains("IPD") {
                mismatch_ct += 1;
            } else if iso == "???" {  // Undetermined isotype
                undet_ct += 1;
            } else if iso.contains("SeC") {
                selcys_ct += 1;
            } else if iso == "Sup" {
                stop_sup_ct += 1;
            } else {
                trna_ct += 1;
            }

            // Count isotype
            let iso_key = if iso.contains("SeC") {
                "SelCys".to_string()
            } else if iso == "Sup" {
                "Supres".to_string()
            } else {
                iso.to_string()
            };

            *iso_counts.entry(iso_key.clone()).or_insert(0) += 1;

            // Count anticodon
            ac_counts.entry(iso_key.clone()).or_insert_with(HashMap::new)
                .entry(ac.to_string()).or_insert(0);
            *ac_counts.get_mut(&iso_key).unwrap().get_mut(ac).unwrap() += 1;

            // Count introns
            if istart != "0" && !istart.is_empty() {
                let introns: Vec<&str> = istart.split(',').collect();
                intron_ct += introns.len();

                *intron_ac_counts.entry(iso_key.clone()).or_insert_with(HashMap::new)
                    .entry(ac.to_string()).or_insert(0) += introns.len();
            }
        }

        let total = trna_ct + selcys_ct + pseudo_ct + mismatch_ct + undet_ct + stop_sup_ct;

        // Output summary counts
        writeln!(w)?;
        writeln!(w, "tRNAs decoding Standard 20 AA:              {}", trna_ct)?;
        writeln!(w, "Selenocysteine tRNAs (TCA):                 {}", selcys_ct)?;
        writeln!(w, "Possible suppressor tRNAs (CTA,TTA,TCA):    {}", stop_sup_ct)?;
        writeln!(w, "tRNAs with undetermined/unknown isotypes:   {}", undet_ct)?;

        if show_mismatch {
            writeln!(w, "tRNAs with mismatch isotypes:               {}", mismatch_ct)?;
        }

        writeln!(w, "Predicted pseudogenes:                      {}", pseudo_ct)?;
        writeln!(w, "                                            -------")?;
        writeln!(w, "Total tRNAs:                                {}", total)?;
        writeln!(w)?;
        writeln!(w, "tRNAs with introns:     \t{}", intron_ct)?;
        writeln!(w)?;

        // Output intron breakdown
        for aa in gc_isotypes {
            if let Some(acsets) = gc_ac_list.get(*aa) {
                for acset in acsets {
                    for ac in acset {
                        if let Some(intron_acs) = intron_ac_counts.get(*aa) {
                            if let Some(&count) = intron_acs.get(ac) {
                                write!(w, "| {}-{}: {} ", aa, ac, count)?;
                            }
                        }
                    }
                }
            }
        }

        if !intron_ac_counts.is_empty() {
            writeln!(w, "|")?;
            writeln!(w)?;
        }

        // Output isotype/anticodon table
        writeln!(w, "Isotype / Anticodon Counts:")?;
        writeln!(w)?;

        for aa in gc_isotypes {
            let mut label = aa.to_string();
            let mut iso_count = *iso_counts.get(*aa).unwrap_or(&0);
            let mut iso_cm_count = *iso_cm_counts.get(*aa).unwrap_or(&0);

            // Handle special cases: Met/iMet, Met/fMet
            if *aa == "Met" {
                if iso_counts.contains_key("iMet") {
                    label = "Met/iMet".to_string();
                    iso_count += iso_counts.get("iMet").unwrap_or(&0);
                    iso_cm_count += iso_cm_counts.get("iMet").unwrap_or(&0);
                } else if iso_counts.contains_key("fMet") {
                    label = "Met/fMet".to_string();
                    iso_count += iso_counts.get("fMet").unwrap_or(&0);
                    iso_cm_count += iso_cm_counts.get("fMet").unwrap_or(&0);
                }
            } else if *aa == "Ile" {
                if iso_counts.contains_key("Ile2") {
                    iso_count += iso_counts.get("Ile2").unwrap_or(&0);
                    iso_cm_count += iso_cm_counts.get("Ile2").unwrap_or(&0);
                }
            }

            if show_mismatch {
                write!(w, "{:<8}: {} ({})\t", label, iso_count, iso_cm_count)?;
            } else {
                write!(w, "{:<8}: {}\t", label, iso_count)?;
            }

            // Output anticodon counts
            if let Some(acsets) = gc_ac_list.get(*aa) {
                for acset in acsets {
                    for ac in acset {
                        if ac == "&nbsp" {
                            write!(w, "             ")?;
                        } else {
                            let mut count = 0;
                            if let Some(ac_map) = ac_counts.get(*aa) {
                                count = *ac_map.get(ac.as_str()).unwrap_or(&0);
                            }

                            // Handle Met/iMet and Ile2 special cases
                            if *aa == "Met" {
                                if let Some(ac_map) = ac_counts.get("iMet") {
                                    count += ac_map.get(ac.as_str()).unwrap_or(&0);
                                } else if let Some(ac_map) = ac_counts.get("fMet") {
                                    count += ac_map.get(ac.as_str()).unwrap_or(&0);
                                }
                            } else if *aa == "Ile" {
                                if let Some(ac_map) = ac_counts.get("Ile2") {
                                    count += ac_map.get(ac.as_str()).unwrap_or(&0);
                                }
                            }

                            write!(w, "{:>5}: {:<6}", ac, count)?;
                        }
                    }
                }
            }

            writeln!(w)?;
        }

        writeln!(w)?;

        Ok(())
    }
}

impl Default for ScanStats {
    fn default() -> Self {
        Self::new()
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_stats_new() {
        let stats = ScanStats::new();
        assert_eq!(stats.numscanned(), 0);
        assert_eq!(stats.seqs_hit(), 0);
        assert_eq!(stats.trnatotal(), 0);
    }

    #[test]
    fn test_stats_counters() {
        let mut stats = ScanStats::new();

        stats.inc_numscanned();
        assert_eq!(stats.numscanned(), 1);

        stats.increment_numscanned(5);
        assert_eq!(stats.numscanned(), 6);

        stats.inc_trnatotal();
        assert_eq!(stats.trnatotal(), 1);

        stats.inc_trnatotal();
        stats.dec_trnatotal();
        assert_eq!(stats.trnatotal(), 1);
    }

    #[test]
    fn test_stats_base_counts() {
        let mut stats = ScanStats::new();

        stats.increment_first_pass_base_ct(10000);
        assert_eq!(stats.first_pass_base_ct(), 10000);

        stats.increment_fpass_trna_base_ct(500);
        assert_eq!(stats.fpass_trna_base_ct(), 500);

        stats.increment_secpass_base_ct(300);
        assert_eq!(stats.secpass_base_ct(), 300);
    }

    #[test]
    fn test_stats_file_name() {
        let mut stats = ScanStats::new();
        stats.set_file_name("test.log");
        assert_eq!(stats.file_name(), "test.log");
    }

    #[test]
    fn test_stats_timers() {
        let mut stats = ScanStats::new();

        stats.start_fp_timer();
        std::thread::sleep(std::time::Duration::from_millis(10));
        stats.end_fp_timer();

        // Should have non-zero CPU time
        assert!(stats.fp_script_cpu > 0.0);

        stats.start_sp_timer();
        std::thread::sleep(std::time::Duration::from_millis(10));
        stats.end_sp_timer();

        assert!(stats.sp_script_cpu > 0.0);
    }

    #[test]
    fn test_stats_clear() {
        let mut stats = ScanStats::new();
        stats.inc_numscanned();
        stats.inc_trnatotal();
        stats.increment_first_pass_base_ct(1000);

        stats.clear();

        assert_eq!(stats.numscanned(), 0);
        assert_eq!(stats.trnatotal(), 0);
        assert_eq!(stats.first_pass_base_ct(), 0);
    }

    #[test]
    fn test_write_firstpass_stats() {
        let mut stats = ScanStats::new();
        stats.set_numscanned(10);
        stats.set_seqs_hit(8);
        stats.increment_first_pass_base_ct(100000);
        stats.set_fpass_trna_base_ct(5000);
        stats.set_trnatotal(65);
        stats.start_fp_timer();
        std::thread::sleep(std::time::Duration::from_millis(10));
        stats.end_fp_timer();

        let mut output = Vec::new();
        stats.save_firstpass_stats(&mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(result.contains("First-pass Stats:"));
        assert!(result.contains("Sequences read:         10"));
        assert!(result.contains("Seqs w/at least 1 hit:  8"));
        assert!(result.contains("tRNAs predicted:        65"));
    }
}
