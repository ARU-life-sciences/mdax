// reading and writing fasta files

use anyhow::Result;
use bio::io::fasta;
use flate2::read::MultiGzDecoder;
use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::Path,
};

/// Open a FASTA file that may be plain or gzipped.
pub fn open_fasta_reader<P: AsRef<Path>>(
    path: P,
) -> Result<fasta::Reader<Box<dyn BufRead + Send>>> {
    let path = path.as_ref();
    let file = File::open(path)?;

    // Decide whether to decompress based on extension.
    // You could also sniff the magic bytes 0x1f, 0x8b if you want.
    let boxed_reader: Box<dyn BufRead + Send> = match path.extension().and_then(|e| e.to_str()) {
        Some("gz") => {
            let gz = MultiGzDecoder::new(file);
            Box::new(BufReader::new(gz))
        }
        _ => Box::new(BufReader::new(file)),
    };

    Ok(fasta::Reader::from_bufread(boxed_reader))
}

pub fn open_tsv_writer<P: AsRef<Path>>(path: P) -> Result<Box<dyn Write>> {
    let path = path.as_ref();

    // stdout
    // TODO: maybe buffer output
    if path.to_string_lossy() == "-" {
        return Ok(Box::new(std::io::stdout())); // fine for now (unbuffered), or wrap in BufWriter if you prefer
    }

    let file = File::create(path)?;
    Ok(Box::new(BufWriter::new(file)))
}

// Open an output FASTA writer (plain, not gzipped yet).
pub fn open_fasta_writer<P: AsRef<Path>>(path: P) -> Result<fasta::Writer<Box<dyn Write>>> {
    let path = path.as_ref();

    // You said .gz not supported yet for output, so fail loudly if asked.
    if path.extension().and_then(|e| e.to_str()) == Some("gz") {
        anyhow::bail!("Output .gz not supported yet (got {})", path.display());
    }

    let file = File::create(path)?;
    let w = BufWriter::new(file);
    Ok(fasta::Writer::new(Box::new(w)))
}

// Write one FASTA record.
pub fn write_fasta_record<W: Write>(w: &mut fasta::Writer<W>, id: &str, seq: &[u8]) -> Result<()> {
    // bio::io::fasta::Record wants owned Strings/Vec<u8>.
    // This is fine for now; we can optimize later if needed.
    let rec = fasta::Record::with_attrs(id, None, seq);
    w.write_record(&rec)?;
    Ok(())
}

// Construct a split record ID like "read123|p1".
pub fn split_id(base: &str, part_index_1based: usize) -> String {
    format!("{base}|p{part_index_1based}")
}
