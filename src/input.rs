use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::Path;

#[derive(Clone, Copy)]
/// Compression detected from an input stream's leading bytes.
pub enum Compression {
    Plain,
    Gzip,
}

impl Compression {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Plain => "plain",
            Self::Gzip => "gzip",
        }
    }
}

/// Opens a plain or gzip-compressed input as a buffered stream.
pub fn open(path: &Path) -> Result<(Box<dyn BufRead>, Compression), String> {
    let mut file = File::open(path).map_err(|error| error.to_string())?;
    let mut magic = [0; 2];
    let count = file.read(&mut magic).map_err(|error| error.to_string())?;
    file.seek(SeekFrom::Start(0))
        .map_err(|error| error.to_string())?;
    let compression = if count == 2 && magic == [0x1f, 0x8b] {
        Compression::Gzip
    } else {
        Compression::Plain
    };
    let reader: Box<dyn BufRead> = match compression {
        Compression::Gzip => Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(file))),
        Compression::Plain => Box::new(BufReader::new(file)),
    };
    Ok((reader, compression))
}
