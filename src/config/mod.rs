/// Configuration management for tRNAscan-SE
///
/// This module handles reading and parsing configuration files that specify
/// paths to models, binaries, and other resources.

use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

/// Configuration error type
#[derive(Debug)]
pub enum ConfigError {
    FileNotFound(PathBuf),
    ParseError(String),
    IoError(std::io::Error),
}

impl std::fmt::Display for ConfigError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ConfigError::FileNotFound(path) => write!(f, "Config file not found: {:?}", path),
            ConfigError::ParseError(msg) => write!(f, "Parse error: {}", msg),
            ConfigError::IoError(e) => write!(f, "IO error: {}", e),
        }
    }
}

impl std::error::Error for ConfigError {}

impl From<std::io::Error> for ConfigError {
    fn from(err: std::io::Error) -> Self {
        ConfigError::IoError(err)
    }
}

/// Main configuration structure
pub struct Config {
    /// Installation directory
    pub install_dir: PathBuf,
    /// Models directory (covariance models)
    pub models_dir: PathBuf,
    /// Temporary directory
    pub temp_dir: PathBuf,
    /// Binary directory
    pub bin_dir: PathBuf,
    /// Library directory
    pub lib_dir: PathBuf,
    /// Additional key-value pairs
    values: HashMap<String, String>,
    /// Nested key-value pairs (key.subkey format)
    nested_values: HashMap<String, HashMap<String, String>>,
}

impl Config {
    /// Load configuration with sensible defaults
    pub fn load() -> Self {
        Self::default()
    }

    /// Load configuration from a file
    pub fn from_file(path: &Path) -> Result<Self, ConfigError> {
        if !path.exists() {
            return Err(ConfigError::FileNotFound(path.to_path_buf()));
        }

        let content = fs::read_to_string(path)?;
        Self::parse(&content)
    }

    /// Parse configuration from string content
    fn parse(content: &str) -> Result<Self, ConfigError> {
        let mut config = Self::default();

        for line in content.lines() {
            let line = line.trim();

            // Skip empty lines and comments
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            // Parse "key: value" format
            if let Some((key, value)) = line.split_once(':') {
                let key = key.trim();
                let mut value = value.trim().to_string();

                // Variable substitution {var}
                value = config.substitute_vars(&value);

                // Process $$ for PID
                if value.contains("$$") {
                    let pid = std::process::id();
                    value = value.replace("$$", &pid.to_string());
                }

                // Check for nested key (key.subkey format)
                if let Some((parent, subkey)) = key.split_once('.') {
                    config
                        .nested_values
                        .entry(parent.to_string())
                        .or_insert_with(HashMap::new)
                        .insert(subkey.to_string(), value);
                } else {
                    // Update standard fields
                    match key {
                        "install_dir" => config.install_dir = PathBuf::from(&value),
                        "models_dir" => config.models_dir = PathBuf::from(&value),
                        "temp_dir" => {
                            // Don't override if already set
                            if config.temp_dir == std::env::temp_dir() {
                                config.temp_dir = PathBuf::from(&value);
                            }
                        }
                        "bin_dir" => config.bin_dir = PathBuf::from(&value),
                        "lib_dir" => config.lib_dir = PathBuf::from(&value),
                        _ => {
                            config.values.insert(key.to_string(), value);
                        }
                    }
                }
            }
        }

        Ok(config)
    }

    /// Substitute {var} references with actual values
    fn substitute_vars(&self, value: &str) -> String {
        let mut result = value.to_string();

        // Find all {var} patterns and replace them
        while let Some(start) = result.find('{') {
            if let Some(end) = result[start..].find('}') {
                let var_name = &result[start + 1..start + end];
                if let Some(var_value) = self.values.get(var_name) {
                    result = result.replace(&format!("{{{}}}", var_name), var_value);
                } else {
                    break; // Avoid infinite loop on undefined vars
                }
            } else {
                break;
            }
        }

        result
    }

    /// Get a configuration value
    pub fn get(&self, key: &str) -> Option<&String> {
        self.values.get(key)
    }

    /// Get a nested configuration value (key.subkey)
    pub fn get_subvalue(&self, key: &str, subkey: &str) -> Option<&String> {
        self.nested_values.get(key)?.get(subkey)
    }

    /// Set temporary directory
    pub fn set_temp_dir(&mut self, dir: PathBuf) {
        self.temp_dir = dir;
    }

    /// Get path to a covariance model file
    pub fn get_cm_path(&self, name: &str) -> PathBuf {
        self.models_dir.join(format!("{}.cm", name))
    }

    /// Get path to a binary
    pub fn get_bin_path(&self, name: &str) -> PathBuf {
        self.bin_dir.join(name)
    }

    /// Get home directory
    pub fn home_dir() -> Option<PathBuf> {
        std::env::var("HOME").ok().map(PathBuf::from)
    }
}

impl Default for Config {
    fn default() -> Self {
        // Try to determine install directory from environment or defaults
        let install_dir = std::env::var("TRNASCAN_INSTALL_DIR")
            .ok()
            .map(PathBuf::from)
            .unwrap_or_else(|| PathBuf::from("/usr/local"));

        let models_dir = install_dir.join("share/tRNAscan-SE/models");
        let bin_dir = install_dir.join("bin");
        let lib_dir = install_dir.join("lib/tRNAscan-SE");
        let temp_dir = std::env::temp_dir();

        Self {
            install_dir,
            models_dir,
            temp_dir,
            bin_dir,
            lib_dir,
            values: HashMap::new(),
            nested_values: HashMap::new(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_default_config() {
        let config = Config::default();
        assert!(config.install_dir.to_str().is_some());
        assert!(config.temp_dir.exists());
    }

    #[test]
    fn test_parse_config() {
        let content = r#"
# Comment line
install_dir: /opt/trnascan

models_dir: {install_dir}/models
temp_dir: /tmp

# Nested values
search.evalue: 0.01
search.threads: 4
"#;

        let config = Config::parse(content).unwrap();
        assert_eq!(config.install_dir, PathBuf::from("/opt/trnascan"));
        assert_eq!(config.get_subvalue("search", "evalue"), Some(&"0.01".to_string()));
        assert_eq!(config.get_subvalue("search", "threads"), Some(&"4".to_string()));
    }

    #[test]
    fn test_load_from_file() {
        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "install_dir: /test/path").unwrap();
        writeln!(temp_file, "models_dir: /test/models").unwrap();

        let config = Config::from_file(temp_file.path()).unwrap();
        assert_eq!(config.install_dir, PathBuf::from("/test/path"));
        assert_eq!(config.models_dir, PathBuf::from("/test/models"));
    }

    #[test]
    fn test_get_cm_path() {
        let config = Config::default();
        let path = config.get_cm_path("trna");
        assert!(path.to_str().unwrap().ends_with("trna.cm"));
    }

    #[test]
    fn test_variable_substitution() {
        let content = r#"
base_dir: /opt/app
models_dir: {base_dir}/models
config_file: {base_dir}/config.txt
"#;

        let config = Config::parse(content).unwrap();
        assert_eq!(config.models_dir, PathBuf::from("/opt/app/models"));
        assert_eq!(config.get("config_file"), Some(&"/opt/app/config.txt".to_string()));
    }
}
