use sha2::Digest;

#[derive(Clone, Debug, PartialEq, Eq)]
/// SHA-256 and BLAKE3 digests of one normalized sequence.
pub struct SequenceDigest {
    pub sha256: String,
    pub blake3: String,
}

/// Incrementally hashes sequence bytes after ASCII case normalization.
pub struct SequenceHasher {
    sha256: sha2::Sha256,
    blake3: blake3::Hasher,
}

impl SequenceHasher {
    /// Creates hashers for both supported digest algorithms.
    pub fn new() -> Self {
        Self {
            sha256: sha2::Sha256::new(),
            blake3: blake3::Hasher::new(),
        }
    }

    /// Adds sequence bytes, normalized to uppercase, to both digests.
    pub fn update(&mut self, sequence: &[u8]) {
        for base in sequence {
            let normalized = base.to_ascii_uppercase();
            self.sha256.update([normalized]);
            self.blake3.update(&[normalized]);
        }
    }

    /// Finalizes both digests as lowercase hexadecimal strings.
    pub fn finish(self) -> SequenceDigest {
        SequenceDigest {
            sha256: format!("{:x}", self.sha256.finalize()),
            blake3: self.blake3.finalize().to_hex().to_string(),
        }
    }
}

/// Computes both supported digests for a complete sequence.
pub fn sequence_digest(sequence: &[u8]) -> SequenceDigest {
    let mut hasher = SequenceHasher::new();
    hasher.update(sequence);
    hasher.finish()
}

#[cfg(test)]
mod tests {
    use super::sequence_digest;

    #[test]
    fn normalizes_case() {
        let digest = sequence_digest(b"acgtn");
        assert_eq!(digest, sequence_digest(b"ACGTN"));
        assert_eq!(
            digest.sha256,
            "d254552eaf2579aa2ecb2a56439c41472f8e7de08ab3f15e898705eada76fc2d"
        );
        assert_eq!(
            digest.blake3,
            "118ce2eb6ea55868dc0f8b0fb1c5a83dcc6b0967bb5407523e72795bebed3478"
        );
    }
}
