/// Logging utilities for tRNAscan-SE
///
/// This module provides structured logging with file output and console display.

use chrono::Local;
use std::fs::{File, OpenOptions};
use std::io::{self, Write};
use std::path::{Path, PathBuf};

/// Log file handler with support for file and console output
pub struct LogFile {
    file: Option<File>,
    file_name: Option<PathBuf>,
    quiet_mode: bool,
}

impl LogFile {
    /// Create a new log file handler
    ///
    /// # Arguments
    /// * `log_path` - Optional path to log file. If None, logging is disabled.
    /// * `quiet` - If true, suppress console output
    pub fn new(log_path: Option<&Path>, quiet: bool) -> io::Result<Self> {
        let file = if let Some(path) = log_path {
            Some(
                OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open(path)?,
            )
        } else {
            None
        };

        Ok(Self {
            file,
            file_name: log_path.map(|p| p.to_path_buf()),
            quiet_mode: quiet,
        })
    }

    /// Create log file with default naming (app-timestamp-pid.log)
    pub fn with_default_name(app: &str, log_dir: &Path, quiet: bool) -> io::Result<Self> {
        let timestamp = Local::now().format("%Y%m%d-%H:%M:%S");
        let pid = std::process::id();
        let filename = format!("{}-{}-{}.log", app, timestamp, pid);
        let path = log_dir.join(filename);

        Self::new(Some(&path), quiet)
    }

    /// Set quiet mode (suppress console output)
    pub fn set_quiet(&mut self, quiet: bool) {
        self.quiet_mode = quiet;
    }

    /// Get the log file path
    pub fn file_name(&self) -> Option<&Path> {
        self.file_name.as_deref()
    }

    /// Write a line to the log file only
    fn write_to_file(&mut self, line: &str) -> io::Result<()> {
        if let Some(file) = &mut self.file {
            writeln!(file, "{}", line)?;
            file.flush()?;
        }
        Ok(())
    }

    /// Write a line to both file and console
    fn write_broadcast(&mut self, line: &str) -> io::Result<()> {
        self.write_to_file(line)?;
        if !self.quiet_mode {
            println!("{}", line);
        }
        Ok(())
    }

    /// Initialize log file with application header
    pub fn initialize(&mut self, app: &str, cmd: Option<&str>) -> io::Result<()> {
        if self.file.is_none() {
            return Ok(());
        }

        self.write_to_file(&format!("Application: {}", app))?;

        // Get username
        let user = std::env::var("USER")
            .or_else(|_| std::env::var("USERNAME"))
            .unwrap_or_else(|_| "unknown".to_string());
        self.write_to_file(&format!("User: {}", user))?;

        // Get hostname
        let hostname = hostname::get()
            .ok()
            .and_then(|h| h.into_string().ok())
            .unwrap_or_else(|| "unknown".to_string());
        self.write_to_file(&format!("Host: {}", hostname))?;

        self.write_to_file(&format!("Start Time: {}", Local::now()))?;
        self.write_to_file("")?;

        if let Some(command) = cmd {
            self.write_to_file(&format!("Command: {}", command))?;
            self.write_to_file("")?;
        }

        Ok(())
    }

    /// Finish logging and write end time
    pub fn finish(&mut self) -> io::Result<()> {
        if self.file.is_some() {
            self.write_to_file("")?;
            self.write_to_file(&format!("End Time: {}", Local::now()))?;
        }
        Ok(())
    }

    /// Log a command
    pub fn command(&mut self, cmd: &str) -> io::Result<()> {
        self.write_to_file(&format!("Command: {}", cmd))
    }

    /// Broadcast a message (to both file and console)
    pub fn broadcast(&mut self, msg: &str) -> io::Result<()> {
        self.write_broadcast(msg)
    }

    /// Log a status message
    pub fn status(&mut self, msg: &str) -> io::Result<()> {
        self.write_broadcast(&format!("Status: {}", msg))
    }

    /// Log an error message
    pub fn error(&mut self, msg: &str) -> io::Result<()> {
        self.write_broadcast(&format!("Error: {}", msg))
    }

    /// Log a warning message
    pub fn warning(&mut self, msg: &str) -> io::Result<()> {
        self.write_broadcast(&format!("Warning: {}", msg))
    }

    /// Log a debug message (file only)
    pub fn debug(&mut self, msg: &str) -> io::Result<()> {
        self.write_to_file(&format!("Debug: {}", msg))
    }

    /// Log an info message
    pub fn info(&mut self, msg: &str) -> io::Result<()> {
        self.write_to_file(&format!("Info: {}", msg))
    }

    /// Close the log file
    pub fn close(&mut self) -> io::Result<()> {
        if let Some(mut file) = self.file.take() {
            file.flush()?;
        }
        Ok(())
    }
}

impl Drop for LogFile {
    fn drop(&mut self) {
        let _ = self.close();
    }
}

/// Simple logger that writes to stderr
pub struct SimpleLogger {
    quiet: bool,
}

impl SimpleLogger {
    pub fn new(quiet: bool) -> Self {
        Self { quiet }
    }

    pub fn log(&self, msg: &str) {
        if !self.quiet {
            eprintln!("{}", msg);
        }
    }

    pub fn status(&self, msg: &str) {
        self.log(&format!("Status: {}", msg));
    }

    pub fn error(&self, msg: &str) {
        self.log(&format!("Error: {}", msg));
    }

    pub fn warning(&self, msg: &str) {
        self.log(&format!("Warning: {}", msg));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use tempfile::tempdir;

    #[test]
    fn test_logfile_creation() {
        let dir = tempdir().unwrap();
        let log_path = dir.path().join("test.log");

        let mut log = LogFile::new(Some(&log_path), true).unwrap();
        log.initialize("test-app", Some("test command")).unwrap();
        log.status("Test message").unwrap();
        log.finish().unwrap();
        drop(log);

        let content = fs::read_to_string(&log_path).unwrap();
        assert!(content.contains("Application: test-app"));
        assert!(content.contains("Command: test command"));
        assert!(content.contains("Status: Test message"));
    }

    #[test]
    fn test_logfile_no_file() {
        let mut log = LogFile::new(None, false).unwrap();
        assert!(log.status("Test").is_ok());
        assert!(log.file_name().is_none());
    }

    #[test]
    fn test_simple_logger() {
        let logger = SimpleLogger::new(true);
        logger.status("This should not panic");
        logger.error("Nor should this");
    }
}
