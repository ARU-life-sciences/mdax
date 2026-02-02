use anyhow::Result;
use needletail::FastxReader;
use std::{
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

// Open a FASTA file that may be plain or gzipped.
pub fn open_fasta_reader<P: AsRef<Path>>(path: P) -> Result<Box<dyn FastxReader>> {
    let path = path.as_ref();
    needletail::parse_fastx_file(path).map_err(|e| anyhow::anyhow!(e))
}

pub fn open_tsv_writer<P: AsRef<Path>>(path: P) -> Result<Box<dyn Write>> {
    let path = path.as_ref();

    // stdout
    // TODO: maybe buffer output
    if path.to_string_lossy() == "-" {
        return Ok(Box::new(std::io::stdout()));
    }

    let file = File::create(path)?;
    Ok(Box::new(BufWriter::new(file)))
}

// Construct a split record ID like "read123|p1".
pub fn split_id(base: &str, part_index_1based: usize) -> String {
    format!("{base}|p{part_index_1based}")
}
