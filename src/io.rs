//! IO for `mdax`
//!
//! Used to be bigger, but now we just open the fasta reader
//! which is essentially a needletail wrapper.
use anyhow::Result;
use needletail::FastxReader;
use std::path::Path;

/// Open a FASTA file that may be plain or gzipped.
pub fn open_fasta_reader<P: AsRef<Path>>(path: P) -> Result<Box<dyn FastxReader>> {
    let path = path.as_ref();
    needletail::parse_fastx_file(path).map_err(|e| anyhow::anyhow!(e))
}
