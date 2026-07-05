//! Output formatting for tRNAscan-SE results.
//!
//! This module provides formatters for all output formats:
//! - Tabular output (.out): Main results in tab-delimited format
//! - Secondary structure (.ss): Detailed sequence and structure
//! - Statistics (.stats): Summary statistics
//! - BED format (.bed): UCSC Genome Browser compatible
//! - Isotype output (.iso): Isotype-specific results

use std::io::{self, Write};

use super::config::TrnaScanConfig;
use super::scanner::{ScanResults, TrnaHit};

// ============================================================================
// Tabular Output (.out)
// ============================================================================

/// Format the header line for tabular output.
pub fn format_tabular_header(config: &TrnaScanConfig) -> String {
    let mut header = String::new();

    // Main header line
    header.push_str("Sequence \t\ttRNA \tBounds\ttRNA\tAnti\tIntron Bounds\t");

    if config.show_breakdown {
        header.push_str("Inf\tHMM\t2'Str\t");
    } else {
        header.push_str("Inf\t");
    }

    if config.show_hit_source {
        header.push_str("Hit\t");
    }

    header.push_str("Isotype\tIsotype\t      \n");

    // Field names line
    header.push_str("Name     \ttRNA #\tBegin\tEnd  \tType\tCodon\tBegin\tEnd\t");

    if config.show_breakdown {
        header.push_str("Score\tScore\tScore\t");
    } else {
        header.push_str("Score\t");
    }

    if config.show_hit_source {
        header.push_str("Origin\t");
    }

    header.push_str("CM\tScore\tNote\n");

    // Separator line
    header.push_str("-------- \t------\t-----\t------\t----\t-----\t-----\t----\t");

    if config.show_breakdown {
        header.push_str("------\t-----\t-----\t");
    } else {
        header.push_str("------\t");
    }

    if config.show_hit_source {
        header.push_str("------\t");
    }

    header.push_str("-------\t-------\t------\n");

    header
}

/// Format a single hit as a tabular line.
pub fn format_tabular_line(hit: &TrnaHit, config: &TrnaScanConfig) -> String {
    let mut line = String::new();

    // Sequence name and tRNA number
    line.push_str(&format!("{} \t{}\t", hit.seq_name, hit.trna_num));

    // Bounds
    line.push_str(&format!("{}\t{}\t", hit.start, hit.end));

    // Type and anticodon
    line.push_str(&format!("{}\t{}\t", hit.isotype, hit.anticodon));

    // Intron bounds
    line.push_str(&format!("{}\t{}\t", hit.intron_start, hit.intron_end));

    // Scores
    if config.show_breakdown {
        line.push_str(&format!(
            "{:.1}\t{:.2}\t{:.2}\t",
            hit.inf_score, hit.hmm_score, hit.ss_score
        ));
    } else {
        line.push_str(&format!("{:.1}\t", hit.inf_score));
    }

    // Hit origin
    if config.show_hit_source {
        line.push_str(&format!("{}\t", hit.origin.short_name()));
    }

    // Isotype CM and score
    if hit.isotype_cm.is_empty() {
        line.push_str(&format!("{}\t", hit.isotype));
    } else {
        line.push_str(&format!("{}\t", hit.isotype_cm));
    }
    line.push_str(&format!("{:.1}\t", hit.isotype_score));

    // Note
    line.push_str(&hit.note);
    line.push('\n');

    line
}

/// Tabular output formatter.
pub struct TabularFormatter<W: Write> {
    writer: W,
    config: TrnaScanConfig,
    header_written: bool,
}

impl<W: Write> TabularFormatter<W> {
    /// Create a new tabular formatter.
    pub fn new(writer: W, config: TrnaScanConfig) -> Self {
        Self {
            writer,
            config,
            header_written: false,
        }
    }

    /// Write the header if not already written.
    pub fn write_header(&mut self) -> io::Result<()> {
        if !self.header_written {
            let header = format_tabular_header(&self.config);
            self.writer.write_all(header.as_bytes())?;
            self.header_written = true;
        }
        Ok(())
    }

    /// Write a single hit.
    pub fn write_hit(&mut self, hit: &TrnaHit) -> io::Result<()> {
        self.write_header()?;
        let line = format_tabular_line(hit, &self.config);
        self.writer.write_all(line.as_bytes())
    }

    /// Write all hits from results.
    pub fn write_results(&mut self, results: &ScanResults) -> io::Result<()> {
        self.write_header()?;
        for hit in &results.hits {
            self.write_hit(hit)?;
        }
        Ok(())
    }

    /// Finish writing and return the inner writer.
    pub fn finish(self) -> W {
        self.writer
    }
}

// ============================================================================
// Secondary Structure Output (.ss)
// ============================================================================

/// Format a secondary structure entry.
pub fn format_ss_entry(hit: &TrnaHit) -> String {
    let mut entry = String::new();

    // Header line: name (coords) Length: N bp
    entry.push_str(&format!(
        "{}.trna{} ({}-{})\tLength: {} bp\n",
        hit.seq_name,
        hit.trna_num,
        hit.start,
        hit.end,
        hit.full_length()
    ));

    // Type line
    entry.push_str(&format!(
        "Type: {}\tAnticodon: {} at {}-{} ({}-{})\tScore: {:.1}\n",
        hit.isotype,
        hit.anticodon,
        hit.anticodon_pos_start,
        hit.anticodon_pos_end,
        hit.start + hit.anticodon_pos_start - 1,
        hit.start + hit.anticodon_pos_end - 1,
        hit.inf_score
    ));

    // Intron line (if present)
    if hit.has_intron() {
        entry.push_str(&format!(
            "Possible intron: {}-{} ({}-{})\n",
            hit.intron_start - hit.start + 1,
            hit.intron_end - hit.start + 1,
            hit.intron_start,
            hit.intron_end
        ));
    }

    // HMM/SS scores
    entry.push_str(&format!(
        "HMM Sc={:.2}\tSec struct Sc={:.2}\n",
        hit.hmm_score, hit.ss_score
    ));

    // Ruler line
    let seq_len = hit.sequence.len();
    let mut ruler = String::with_capacity(seq_len + 10);
    ruler.push_str("         ");
    for i in 0..seq_len {
        if (i + 1) % 10 == 5 {
            ruler.push('*');
        } else if (i + 1) % 10 == 0 {
            ruler.push('|');
        } else {
            ruler.push(' ');
        }
    }
    entry.push_str(&ruler);
    entry.push('\n');

    // Sequence line
    entry.push_str(&format!("Seq: {}\n", hit.sequence));

    // Structure line
    entry.push_str(&format!("Str: {}\n", hit.secondary_structure));

    // Blank line separator
    entry.push('\n');

    entry
}

/// Secondary structure output formatter.
pub struct SsFormatter<W: Write> {
    writer: W,
}

impl<W: Write> SsFormatter<W> {
    /// Create a new secondary structure formatter.
    pub fn new(writer: W) -> Self {
        Self { writer }
    }

    /// Write a single hit's secondary structure.
    pub fn write_hit(&mut self, hit: &TrnaHit) -> io::Result<()> {
        let entry = format_ss_entry(hit);
        self.writer.write_all(entry.as_bytes())
    }

    /// Write all hits from results.
    pub fn write_results(&mut self, results: &ScanResults) -> io::Result<()> {
        for hit in &results.hits {
            self.write_hit(hit)?;
        }
        Ok(())
    }

    /// Finish writing and return the inner writer.
    pub fn finish(self) -> W {
        self.writer
    }
}

// ============================================================================
// Statistics Output (.stats)
// ============================================================================

/// Format statistics output.
pub fn format_stats(
    results: &ScanResults,
    config: &TrnaScanConfig,
    input_file: &str,
) -> String {
    let mut stats = String::new();

    // Header
    stats.push('\n');
    stats.push_str("tRNAscan-SE v.2.0 (Rust) scan results\n");
    stats.push_str(&format!(
        "Started: {}\n",
        chrono_lite_now()
    ));
    stats.push('\n');
    stats.push_str("------------------------------------------------------------\n");

    // Configuration
    stats.push_str(&format!(
        "Sequence file(s) to search:        {}\n",
        input_file
    ));
    stats.push_str(&format!(
        "Search Mode:                       {}\n",
        config.search_mode
    ));

    if let Some(ref output_file) = config.output.output_file {
        stats.push_str(&format!(
            "Results written to:                {}\n",
            output_file.display()
        ));
    }

    stats.push_str("Output format:                     Tabular\n");
    stats.push_str(&format!(
        "Searching with:                    {} First Pass->Infernal\n",
        config.first_pass
    ));
    stats.push_str(&format!(
        "Isotype-specific model scan:       {}\n",
        if config.isotype_specific { "Yes" } else { "No" }
    ));
    stats.push_str(&format!(
        "Covariance model:                  {}\n",
        config.cm_model.display()
    ));
    stats.push_str(&format!(
        "Infernal first pass cutoff score:  {}\n",
        config.infernal_first_pass_cutoff
    ));
    stats.push('\n');
    stats.push_str(&format!(
        "Temporary directory:               {}\n",
        config.temp_dir.display()
    ));

    if let Some(ref ss_file) = config.output.ss_file {
        stats.push_str(&format!(
            "tRNA secondary structure\n    predictions saved to:          {}\n",
            ss_file.display()
        ));
    }

    if let Some(ref bed_file) = config.output.bed_file {
        stats.push_str(&format!(
            "tRNA predictions saved to:         {}\n",
            bed_file.display()
        ));
    }

    if let Some(ref iso_file) = config.output.iso_file {
        stats.push_str(&format!(
            "Isotype specific\n    predictions saved to:          {}\n",
            iso_file.display()
        ));
    }

    if config.show_breakdown {
        stats.push_str("\nReporting HMM/2' structure score breakdown\n");
    }

    stats.push_str("------------------------------------------------------------\n");
    stats.push('\n');

    // First-pass stats
    stats.push_str("First-pass Stats:\n");
    stats.push_str("---------------\n");
    stats.push_str(&format!(
        "Sequences read:         {}\n",
        results.sequences_scanned
    ));
    stats.push_str(&format!(
        "Seqs w/at least 1 hit:  {}\n",
        results.sequences_with_hits
    ));
    stats.push_str(&format!(
        "Bases read:             {} (x2 for both strands)\n",
        results.bases_scanned
    ));

    let total_trna_bases: i32 = results.hits.iter().map(|h| h.full_length()).sum();
    stats.push_str(&format!("Bases in tRNAs:         {}\n", total_trna_bases));
    stats.push_str(&format!(
        "tRNAs predicted:        {}\n",
        results.hits.len()
    ));
    stats.push_str(&format!(
        "Av. tRNA length:        {}\n",
        results.avg_trna_length.round() as i32
    ));
    stats.push_str(&format!(
        "Script CPU time:        {:.2} s\n",
        results.first_pass_time * 0.1
    ));
    stats.push_str(&format!(
        "Scan CPU time:          {:.2} s\n",
        results.first_pass_time
    ));

    let scan_speed = if results.first_pass_time > 0.0 {
        results.bases_scanned as f64 / results.first_pass_time / 1000.0
    } else {
        0.0
    };
    stats.push_str(&format!("Scan speed:             {:.1} Kbp/sec\n", scan_speed));
    stats.push('\n');

    // Infernal stats
    stats.push_str("Infernal Stats:\n");
    stats.push_str("-----------\n");
    stats.push_str(&format!(
        "Candidate tRNAs read:     {}\n",
        results.hits.len()
    ));
    stats.push_str(&format!(
        "Infernal-confirmed tRNAs:     {}\n",
        results.total_trnas()
    ));

    let bases_scanned_inf = total_trna_bases + (results.hits.len() as i32 * 100);
    stats.push_str(&format!(
        "Bases scanned by Infernal:  {}\n",
        bases_scanned_inf
    ));

    let pct_scanned = if results.bases_scanned > 0 {
        bases_scanned_inf as f64 / results.bases_scanned as f64 * 100.0
    } else {
        0.0
    };
    stats.push_str(&format!("% seq scanned by Infernal:  {:.1} %\n", pct_scanned));
    stats.push_str(&format!(
        "Script CPU time:          {:.2} s\n",
        results.infernal_time * 0.1
    ));
    stats.push_str(&format!(
        "Infernal CPU time:            {:.2} s\n",
        results.infernal_time
    ));

    let inf_speed = if results.infernal_time > 0.0 {
        bases_scanned_inf as f64 / results.infernal_time
    } else {
        0.0
    };
    stats.push_str(&format!("Scan speed:               {:.1} bp/sec\n", inf_speed));
    stats.push('\n');

    // Overall speed
    let overall_speed = results.scan_speed();
    stats.push_str(&format!(
        "Overall scan speed: {:.1} bp/sec\n",
        overall_speed
    ));
    stats.push('\n');

    // Summary counts
    stats.push_str(&format!(
        "tRNAs decoding Standard 20 AA:              {}\n",
        results.standard_aa_count()
    ));
    stats.push_str(&format!(
        "Selenocysteine tRNAs (TCA):                 {}\n",
        results.sec_count
    ));
    stats.push_str(&format!(
        "Possible suppressor tRNAs (CTA,TTA):        {}\n",
        results.sup_count
    ));

    let unk_count = results.isotype_counts.get("Unk").copied().unwrap_or(0);
    stats.push_str(&format!(
        "tRNAs with undetermined/unknown isotypes:   {}\n",
        unk_count
    ));
    stats.push_str(&format!(
        "tRNAs with mismatch isotypes:               {}\n",
        results.mismatch_count
    ));
    stats.push_str(&format!(
        "Predicted pseudogenes:                      {}\n",
        results.pseudogene_count
    ));
    stats.push_str("                                            -------\n");
    stats.push_str(&format!(
        "Total tRNAs:                                {}\n",
        results.hits.len()
    ));
    stats.push('\n');

    // Intron count
    stats.push_str(&format!(
        "tRNAs with introns:     \t{}\n",
        results.trnas_with_introns
    ));
    stats.push('\n');

    // List intron-containing tRNAs
    let intron_hits: Vec<_> = results.hits.iter().filter(|h| h.has_intron()).collect();
    if !intron_hits.is_empty() {
        stats.push_str("| ");
        for hit in intron_hits {
            stats.push_str(&format!("{}-{}: 1 ", hit.isotype, hit.anticodon));
        }
        stats.push_str("|\n");
        stats.push('\n');
    }

    // Isotype/Anticodon table
    stats.push_str("Isotype / Anticodon Counts:\n\n");
    stats.push_str(&format_isotype_anticodon_table(results));
    stats.push('\n');

    stats
}

/// Format the isotype/anticodon count table.
fn format_isotype_anticodon_table(results: &ScanResults) -> String {
    let mut table = String::new();

    // Define amino acids and their anticodons in display order
    let aa_anticodon_order = [
        ("Ala", vec!["AGC", "GGC", "CGC", "TGC"]),
        ("Gly", vec!["ACC", "GCC", "CCC", "TCC"]),
        ("Pro", vec!["AGG", "GGG", "CGG", "TGG"]),
        ("Thr", vec!["AGT", "GGT", "CGT", "TGT"]),
        ("Val", vec!["AAC", "GAC", "CAC", "TAC"]),
        (
            "Ser",
            vec!["AGA", "GGA", "CGA", "TGA", "ACT", "GCT"],
        ),
        (
            "Arg",
            vec!["ACG", "GCG", "CCG", "TCG", "CCT", "TCT"],
        ),
        (
            "Leu",
            vec!["AAG", "GAG", "CAG", "TAG", "CAA", "TAA"],
        ),
        ("Phe", vec!["AAA", "GAA"]),
        ("Asn", vec!["ATT", "GTT"]),
        ("Lys", vec!["CTT", "TTT"]),
        ("Asp", vec!["ATC", "GTC"]),
        ("Glu", vec!["CTC", "TTC"]),
        ("His", vec!["ATG", "GTG"]),
        ("Gln", vec!["CTG", "TTG"]),
        ("Ile", vec!["AAT", "GAT", "TAT"]),
        ("Met", vec!["CAT"]),
        ("Tyr", vec!["ATA", "GTA"]),
        ("Supres", vec!["CTA", "TTA"]),
        ("Cys", vec!["ACA", "GCA"]),
        ("Trp", vec!["CCA"]),
        ("SelCys", vec!["TCA"]),
    ];

    let iso_table = results.get_isotype_anticodon_table();

    for (aa, anticodons) in &aa_anticodon_order {
        let total: usize = if let Some(aa_map) = iso_table.get(*aa) {
            aa_map.values().sum()
        } else {
            0
        };

        // Left column: AA name and count
        let aa_display = if *aa == "Supres" {
            "Supres"
        } else if *aa == "SelCys" {
            "SelCys"
        } else {
            aa
        };

        // Count in parentheses is non-pseudo count
        table.push_str(&format!("{:6}: {} ({})\t", aa_display, total, total));

        // Anticodon columns
        for anticodon in anticodons {
            let count = iso_table
                .get(*aa)
                .and_then(|m| m.get(*anticodon))
                .copied()
                .unwrap_or(0);

            if count > 0 {
                table.push_str(&format!("  {}: {}", anticodon, count));
            } else {
                table.push_str(&format!("  {}:  ", anticodon));
            }
        }
        table.push('\n');
    }

    table
}

/// Statistics formatter.
pub struct StatsFormatter<W: Write> {
    writer: W,
    config: TrnaScanConfig,
    input_file: String,
}

impl<W: Write> StatsFormatter<W> {
    /// Create a new statistics formatter.
    pub fn new(writer: W, config: TrnaScanConfig, input_file: impl Into<String>) -> Self {
        Self {
            writer,
            config,
            input_file: input_file.into(),
        }
    }

    /// Write statistics for results.
    pub fn write_stats(&mut self, results: &ScanResults) -> io::Result<()> {
        let stats = format_stats(results, &self.config, &self.input_file);
        self.writer.write_all(stats.as_bytes())
    }

    /// Finish writing and return the inner writer.
    pub fn finish(self) -> W {
        self.writer
    }
}

// ============================================================================
// BED Format Output (.bed)
// ============================================================================

/// Format a single hit as a BED line.
///
/// BED format: chrom start end name score strand thickStart thickEnd itemRgb blockCount blockSizes blockStarts
pub fn format_bed_line(hit: &TrnaHit) -> String {
    let (bed_start, bed_end) = hit.bed_coords();
    let strand = if hit.forward_strand { '+' } else { '-' };
    let name = hit.identifier();
    let score = (hit.inf_score * 10.0).min(1000.0) as i32;

    // For tRNAs with introns, use blocks
    if hit.has_intron() {
        let intron_start = if hit.forward_strand {
            hit.intron_start - hit.start
        } else {
            hit.end - hit.intron_end
        };
        let intron_end = if hit.forward_strand {
            hit.intron_end - hit.start + 1
        } else {
            hit.end - hit.intron_start + 1
        };

        let block1_size = intron_start;
        let block2_size = hit.full_length() - intron_end;

        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t0\t2\t{},{},\t0,{},\n",
            hit.seq_name,
            bed_start,
            bed_end,
            name,
            score,
            strand,
            bed_start,
            bed_end,
            block1_size,
            block2_size,
            intron_end
        )
    } else {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t0\t1\t{},\t0,\n",
            hit.seq_name,
            bed_start,
            bed_end,
            name,
            score,
            strand,
            bed_start,
            bed_end,
            hit.full_length()
        )
    }
}

/// BED output formatter.
pub struct BedFormatter<W: Write> {
    writer: W,
}

impl<W: Write> BedFormatter<W> {
    /// Create a new BED formatter.
    pub fn new(writer: W) -> Self {
        Self { writer }
    }

    /// Write a single hit.
    pub fn write_hit(&mut self, hit: &TrnaHit) -> io::Result<()> {
        let line = format_bed_line(hit);
        self.writer.write_all(line.as_bytes())
    }

    /// Write all hits from results.
    pub fn write_results(&mut self, results: &ScanResults) -> io::Result<()> {
        // Sort hits by chromosome and position for BED output
        let mut sorted_hits: Vec<_> = results.hits.iter().collect();
        sorted_hits.sort_by(|a, b| {
            a.seq_name
                .cmp(&b.seq_name)
                .then_with(|| a.start.cmp(&b.start))
        });

        for hit in sorted_hits {
            self.write_hit(hit)?;
        }
        Ok(())
    }

    /// Finish writing and return the inner writer.
    pub fn finish(self) -> W {
        self.writer
    }
}

// ============================================================================
// Isotype Output (.iso)
// ============================================================================

/// Format a single hit as an isotype line.
pub fn format_iso_line(hit: &TrnaHit) -> String {
    format!(
        "{}\t{}\t{}\t{}-{}\t{:.1}\t{}\t{}\t{}\t{:.1}\n",
        hit.seq_name,
        hit.trna_num,
        hit.isotype,
        hit.isotype,
        hit.anticodon,
        hit.inf_score,
        hit.start,
        hit.end,
        hit.isotype_cm,
        hit.isotype_score
    )
}

/// Isotype output formatter.
pub struct IsoFormatter<W: Write> {
    writer: W,
}

impl<W: Write> IsoFormatter<W> {
    /// Create a new isotype formatter.
    pub fn new(writer: W) -> Self {
        Self { writer }
    }

    /// Write header.
    pub fn write_header(&mut self) -> io::Result<()> {
        self.writer.write_all(b"SeqName\ttRNA#\tIsotype\tType\tScore\tStart\tEnd\tIsoCM\tIsoScore\n")
    }

    /// Write a single hit.
    pub fn write_hit(&mut self, hit: &TrnaHit) -> io::Result<()> {
        let line = format_iso_line(hit);
        self.writer.write_all(line.as_bytes())
    }

    /// Write all hits from results.
    pub fn write_results(&mut self, results: &ScanResults) -> io::Result<()> {
        for hit in &results.hits {
            self.write_hit(hit)?;
        }
        Ok(())
    }

    /// Finish writing and return the inner writer.
    pub fn finish(self) -> W {
        self.writer
    }
}

// ============================================================================
// Helper Functions
// ============================================================================

/// Get a simple timestamp string (without external dependencies).
fn chrono_lite_now() -> String {
    use std::time::SystemTime;

    // Get current time as duration since UNIX_EPOCH
    match SystemTime::now().duration_since(std::time::UNIX_EPOCH) {
        Ok(duration) => {
            // Just return a simple timestamp
            let secs = duration.as_secs();
            format!("Timestamp: {}", secs)
        }
        Err(_) => "Unknown time".to_string(),
    }
}

/// Write all output files based on configuration.
pub fn write_all_outputs(
    results: &ScanResults,
    config: &TrnaScanConfig,
    input_file: &str,
) -> io::Result<()> {
    // Write tabular output (.out)
    if let Some(ref path) = config.output.output_file {
        let file = std::fs::File::create(path)?;
        let mut formatter = TabularFormatter::new(file, config.clone());
        formatter.write_results(results)?;
    }

    // Write secondary structure (.ss)
    if let Some(ref path) = config.output.ss_file {
        let file = std::fs::File::create(path)?;
        let mut formatter = SsFormatter::new(file);
        formatter.write_results(results)?;
    }

    // Write statistics (.stats)
    if let Some(ref path) = config.output.stats_file {
        let file = std::fs::File::create(path)?;
        let mut formatter = StatsFormatter::new(file, config.clone(), input_file);
        formatter.write_stats(results)?;
    }

    // Write BED format (.bed)
    if let Some(ref path) = config.output.bed_file {
        let file = std::fs::File::create(path)?;
        let mut formatter = BedFormatter::new(file);
        formatter.write_results(results)?;
    }

    // Write isotype output (.iso)
    if let Some(ref path) = config.output.iso_file {
        let file = std::fs::File::create(path)?;
        let mut formatter = IsoFormatter::new(file);
        formatter.write_results(results)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_test_hit() -> TrnaHit {
        let mut hit = TrnaHit::new("CELF22B7", 1, 12619, 12738, "Leu", "CAA");
        hit.inf_score = 74.2;
        hit.hmm_score = 51.20;
        hit.ss_score = 23.00;
        hit.intron_start = 12657;
        hit.intron_end = 12692;
        hit.isotype_cm = "Leu".to_string();
        hit.isotype_score = 119.9;
        hit.origin = super::super::scanner::HitOrigin::Infernal;
        hit.sequence = "GCACGGATGGCCGAGTGGTCTAAGGCGCCAGACTCAAG".to_string();
        hit.secondary_structure = ">>>>>>>..>>>...........<<<.>>>>>".to_string();
        hit
    }

    #[test]
    fn test_tabular_header() {
        let config = TrnaScanConfig::default();
        let header = format_tabular_header(&config);
        assert!(header.contains("Sequence"));
        assert!(header.contains("tRNA"));
        assert!(header.contains("Bounds"));
        assert!(header.contains("Isotype"));
    }

    #[test]
    fn test_tabular_line() {
        let hit = make_test_hit();
        let config = TrnaScanConfig::default();
        let line = format_tabular_line(&hit, &config);
        assert!(line.contains("CELF22B7"));
        assert!(line.contains("12619"));
        assert!(line.contains("Leu"));
        assert!(line.contains("CAA"));
    }

    #[test]
    fn test_tabular_line_with_breakdown() {
        let hit = make_test_hit();
        let config = TrnaScanConfig::default().with_breakdown(true);
        let line = format_tabular_line(&hit, &config);
        assert!(line.contains("74.2"));
        assert!(line.contains("51.20"));
        assert!(line.contains("23.00"));
    }

    #[test]
    fn test_ss_entry() {
        let hit = make_test_hit();
        let entry = format_ss_entry(&hit);
        assert!(entry.contains("CELF22B7.trna1"));
        assert!(entry.contains("12619-12738"));
        assert!(entry.contains("Type: Leu"));
        assert!(entry.contains("Anticodon: CAA"));
    }

    #[test]
    fn test_bed_line() {
        let hit = make_test_hit();
        let line = format_bed_line(&hit);
        assert!(line.contains("CELF22B7"));
        assert!(line.contains("12618")); // 0-based start
        assert!(line.contains("12738")); // End
        assert!(line.contains("+")); // Forward strand
    }

    #[test]
    fn test_bed_line_with_intron() {
        let hit = make_test_hit();
        let line = format_bed_line(&hit);
        // Should have 2 blocks due to intron
        assert!(line.contains("\t2\t"));
    }

    #[test]
    fn test_iso_line() {
        let hit = make_test_hit();
        let line = format_iso_line(&hit);
        assert!(line.contains("CELF22B7"));
        assert!(line.contains("Leu"));
        assert!(line.contains("74.2"));
    }

    #[test]
    fn test_stats_format() {
        let results = ScanResults::new();
        let config = TrnaScanConfig::default();
        let stats = format_stats(&results, &config, "test.fa");
        assert!(stats.contains("tRNAscan-SE"));
        assert!(stats.contains("test.fa"));
        assert!(stats.contains("Eukaryotic"));
    }

    #[test]
    fn test_tabular_formatter() {
        let config = TrnaScanConfig::default();
        let mut output = Vec::new();
        {
            let mut formatter = TabularFormatter::new(&mut output, config);
            let hit = make_test_hit();
            formatter.write_hit(&hit).unwrap();
        }
        let output_str = String::from_utf8(output).unwrap();
        assert!(output_str.contains("CELF22B7"));
    }

    #[test]
    fn test_bed_formatter() {
        let mut output = Vec::new();
        {
            let mut formatter = BedFormatter::new(&mut output);
            let hit = make_test_hit();
            formatter.write_hit(&hit).unwrap();
        }
        let output_str = String::from_utf8(output).unwrap();
        assert!(output_str.contains("CELF22B7"));
        assert!(output_str.contains("tRNA1-LeuCAA"));
    }
}
