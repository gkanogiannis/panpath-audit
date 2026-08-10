#[derive(Clone, Copy, Debug, PartialEq, Eq)]
/// Direction in which a graph segment is traversed.
pub enum Orientation {
    Forward,
    Reverse,
}

impl Orientation {
    /// Parses a GFA P-line orientation suffix.
    pub fn from_path(symbol: char) -> Result<Self, String> {
        match symbol {
            '+' => Ok(Self::Forward),
            '-' => Ok(Self::Reverse),
            _ => Err("unsupported orientation".into()),
        }
    }

    /// Parses a GFA W-line orientation marker.
    pub fn from_walk(symbol: u8) -> Result<Self, String> {
        match symbol {
            b'>' => Ok(Self::Forward),
            b'<' => Ok(Self::Reverse),
            _ => Err("missing orientation".into()),
        }
    }

    /// Converts a traversal offset to a one-based segment coordinate.
    pub fn segment_position(self, index: usize, length: usize) -> usize {
        match self {
            Self::Forward => index + 1,
            Self::Reverse => length - index,
        }
    }
}

impl std::fmt::Display for Orientation {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter.write_str(match self {
            Self::Forward => "+",
            Self::Reverse => "-",
        })
    }
}

/// Returns whether a byte is a valid DNA IUPAC symbol.
pub fn is_iupac(symbol: u8) -> bool {
    matches!(
        symbol.to_ascii_uppercase(),
        b'A' | b'C'
            | b'G'
            | b'T'
            | b'R'
            | b'Y'
            | b'S'
            | b'W'
            | b'K'
            | b'M'
            | b'B'
            | b'D'
            | b'H'
            | b'V'
            | b'N'
    )
}

/// Returns whether a valid nucleotide byte is non-canonical.
pub fn is_ambiguous(symbol: u8) -> bool {
    !matches!(symbol.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')
}

/// Complements one DNA IUPAC symbol while preserving case.
pub fn complement(base: u8) -> Result<u8, String> {
    let lower = base.is_ascii_lowercase();
    let complemented = match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'R' => b'Y',
        b'Y' => b'R',
        b'S' => b'S',
        b'W' => b'W',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'D' => b'H',
        b'H' => b'D',
        b'V' => b'B',
        b'N' => b'N',
        _ => return Err("invalid nucleotide".into()),
    };
    Ok(if lower {
        complemented.to_ascii_lowercase()
    } else {
        complemented
    })
}

/// Parses the oriented segment steps encoded by a GFA W-line walk.
pub fn walk_steps(walk: &str) -> Result<Vec<(&str, Orientation)>, String> {
    let mut steps = Vec::new();
    let mut start = 0;
    let bytes = walk.as_bytes();
    while start < bytes.len() {
        let orientation = Orientation::from_walk(bytes[start])?;
        let name_start = start + 1;
        let end = bytes[name_start..]
            .iter()
            .position(|byte| matches!(byte, b'>' | b'<'))
            .map_or(bytes.len(), |offset| name_start + offset);
        if end == name_start {
            return Err("empty walk step".into());
        }
        steps.push((&walk[name_start..end], orientation));
        start = end;
    }
    if steps.is_empty() {
        Err("empty walk".into())
    } else {
        Ok(steps)
    }
}
