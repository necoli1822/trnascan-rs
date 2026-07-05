// SQUID sequence information types and constants
// Ported from original squid.h

use std::fmt;

// Maximum lengths for SQINFO fields
pub const SQINFO_NAMELEN: usize = 64;
pub const SQINFO_DESCLEN: usize = 128;

// SQINFO flags - indicate which fields are set
pub const SQINFO_NAME: u32 = 1 << 0;
pub const SQINFO_ID: u32 = 1 << 1;
pub const SQINFO_ACC: u32 = 1 << 2;
pub const SQINFO_DESC: u32 = 1 << 3;
pub const SQINFO_START: u32 = 1 << 4;
pub const SQINFO_STOP: u32 = 1 << 5;
pub const SQINFO_LEN: u32 = 1 << 6;
pub const SQINFO_TYPE: u32 = 1 << 7;
pub const SQINFO_WGT: u32 = 1 << 8;
pub const SQINFO_OLEN: u32 = 1 << 9;
pub const SQINFO_SS: u32 = 1 << 10;
pub const SQINFO_SA: u32 = 1 << 12;

// Sequence types
pub const K_OTHER_SEQ: i32 = 0;
pub const K_DNA: i32 = 1;
pub const K_RNA: i32 = 2;
pub const K_AMINO: i32 = 3;

// Sequence file formats
pub const K_UNKNOWN: i32 = 0;
pub const K_IG: i32 = 1;
pub const K_GENBANK: i32 = 2;
pub const K_NBRF: i32 = 3;
pub const K_EMBL: i32 = 4;
pub const K_GCG: i32 = 5;
pub const K_STRIDER: i32 = 6;
pub const K_PEARSON: i32 = 7;
pub const K_ZUKER: i32 = 8;
pub const K_IDRAW: i32 = 9;
pub const K_SELEX: i32 = 10;
pub const K_MSF: i32 = 11;
pub const K_PIR: i32 = 12;
pub const K_RAW: i32 = 13;
pub const K_SQUID: i32 = 14;
pub const K_XPEARSON: i32 = 15;
pub const K_GCGDATA: i32 = 16;
pub const K_CLUSTAL: i32 = 17;

pub const K_MIN_FORMAT: i32 = 1;
pub const K_MAX_FORMAT: i32 = 17;
pub const K_NUM_FORMATS: i32 = K_MAX_FORMAT + 1;
pub const K_NOFORMAT: i32 = -1;

/// Sequence information structure
/// Corresponds to SQINFO in original squid.h
#[derive(Clone, Default, Debug)]
pub struct SqInfo {
    pub flags: u32,
    pub name: String,
    pub id: String,
    pub acc: String,
    pub desc: String,
    pub len: usize,
    pub start: i32,
    pub stop: i32,
    pub olen: i32,
    pub seq_type: i32,
    pub weight: f32,
    pub ss: Option<String>,
    pub sa: Option<String>,
}

impl SqInfo {
    pub fn new() -> Self {
        Self::default()
    }

    /// Set a string field with proper truncation and flag setting
    pub fn set_string_field(&mut self, value: &str, flag: u32) {
        let trimmed = value.trim();
        if trimmed == "-" || trimmed.is_empty() {
            return;
        }

        match flag {
            SQINFO_NAME => {
                self.name = Self::truncate_string(trimmed, SQINFO_NAMELEN);
                self.flags |= SQINFO_NAME;
            }
            SQINFO_ID => {
                self.id = Self::truncate_string(trimmed, SQINFO_NAMELEN);
                self.flags |= SQINFO_ID;
            }
            SQINFO_ACC => {
                self.acc = Self::truncate_string(trimmed, SQINFO_NAMELEN);
                self.flags |= SQINFO_ACC;
            }
            SQINFO_DESC => {
                if self.flags & SQINFO_DESC != 0 {
                    // Append to existing description
                    if self.desc.len() < SQINFO_DESCLEN - 2 {
                        self.desc.push(' ');
                        let remaining = SQINFO_DESCLEN - 1 - self.desc.len();
                        let to_add = &trimmed[..trimmed.len().min(remaining)];
                        self.desc.push_str(to_add);
                    }
                } else {
                    self.desc = Self::truncate_string(trimmed, SQINFO_DESCLEN);
                }
                self.flags |= SQINFO_DESC;
            }
            _ => {}
        }
    }

    /// Set an integer field
    pub fn set_int_field(&mut self, value: i32, flag: u32) {
        match flag {
            SQINFO_START => {
                self.start = value;
                if value != 0 {
                    self.flags |= SQINFO_START;
                }
            }
            SQINFO_STOP => {
                self.stop = value;
                if value != 0 {
                    self.flags |= SQINFO_STOP;
                }
            }
            SQINFO_OLEN => {
                self.olen = value;
                if value != 0 {
                    self.flags |= SQINFO_OLEN;
                }
            }
            _ => {}
        }
    }

    /// Truncate string to maximum length
    fn truncate_string(s: &str, max_len: usize) -> String {
        if s.len() >= max_len {
            s[..max_len - 1].to_string()
        } else {
            s.to_string()
        }
    }

    /// Check if a flag is set
    pub fn has_flag(&self, flag: u32) -> bool {
        self.flags & flag != 0
    }
}

impl fmt::Display for SqInfo {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "SqInfo {{ ")?;
        write!(f, "flags: 0x{:x}", self.flags)?;
        if self.has_flag(SQINFO_NAME) {
            write!(f, ", name: \"{}\"", self.name)?;
        }
        if self.has_flag(SQINFO_ID) {
            write!(f, ", id: \"{}\"", self.id)?;
        }
        if self.has_flag(SQINFO_ACC) {
            write!(f, ", acc: \"{}\"", self.acc)?;
        }
        if self.has_flag(SQINFO_DESC) {
            write!(f, ", desc: \"{}\"", self.desc)?;
        }
        if self.has_flag(SQINFO_LEN) {
            write!(f, ", len: {}", self.len)?;
        }
        if self.has_flag(SQINFO_START) {
            write!(f, ", start: {}", self.start)?;
        }
        if self.has_flag(SQINFO_STOP) {
            write!(f, ", stop: {}", self.stop)?;
        }
        if self.has_flag(SQINFO_OLEN) {
            write!(f, ", olen: {}", self.olen)?;
        }
        write!(f, " }}")
    }
}

/// SQFILE - file handle for reading sequences
/// Corresponds to ReadSeqVars in original squid.h
pub struct SqFile {
    pub sbuffer: String,
    pub seqlen: usize,
    pub maxseq: usize,
    pub dash_equals_n: bool,
    pub seq: Vec<u8>,
    pub sqinfo: SqInfo,
}

impl SqFile {
    pub const START_LENGTH: usize = 500;

    pub fn new() -> Self {
        Self {
            sbuffer: String::with_capacity(4096),
            seqlen: 0,
            maxseq: Self::START_LENGTH,
            dash_equals_n: false,
            seq: Vec::with_capacity(Self::START_LENGTH),
            sqinfo: SqInfo::new(),
        }
    }
}

impl Default for SqFile {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sqinfo_flags() {
        let mut info = SqInfo::new();
        assert_eq!(info.flags, 0);
        assert!(!info.has_flag(SQINFO_NAME));

        info.set_string_field("test", SQINFO_NAME);
        assert!(info.has_flag(SQINFO_NAME));
        assert_eq!(info.name, "test");
    }

    #[test]
    fn test_sqinfo_truncation() {
        let mut info = SqInfo::new();
        let long_name = "a".repeat(100);
        info.set_string_field(&long_name, SQINFO_NAME);
        assert_eq!(info.name.len(), SQINFO_NAMELEN - 1);
    }

    #[test]
    fn test_sqinfo_desc_append() {
        let mut info = SqInfo::new();
        info.set_string_field("First", SQINFO_DESC);
        info.set_string_field("Second", SQINFO_DESC);
        assert_eq!(info.desc, "First Second");
    }
}
