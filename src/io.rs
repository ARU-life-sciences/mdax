use anyhow::Result;
use needletail::FastxReader;
use noodles::fasta::{
    self as fasta,
    record::{Definition, Sequence},
};
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

// Open an output FASTA writer (plain, not gzipped yet).
pub fn open_fasta_writer<P: AsRef<Path>>(path: P) -> Result<fasta::io::Writer<Box<dyn Write>>> {
    let path = path.as_ref();

    // fail if asked for gz
    if path.extension().and_then(|e| e.to_str()) == Some("gz") {
        anyhow::bail!("Output .gz not supported yet (got {})", path.display());
    }

    let file = File::create(path)?;
    let w = BufWriter::new(file);
    Ok(fasta::io::Writer::new(Box::new(w)))
}

// Write one FASTA record.
pub fn write_fasta_record<W: Write>(
    w: &mut fasta::io::Writer<W>,
    id: &str,
    seq: &[u8],
) -> Result<()> {
    let definition = Definition::new(id, None);
    let sequence = Sequence::from(seq.to_vec());
    let record = fasta::Record::new(definition, sequence);
    w.write_record(&record)?;
    Ok(())
}

// Construct a split record ID like "read123|p1".
pub fn split_id(base: &str, part_index_1based: usize) -> String {
    format!("{base}|p{part_index_1based}")
}
