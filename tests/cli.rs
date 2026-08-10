use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};

static CASE_ID: AtomicU64 = AtomicU64::new(0);

struct Case {
    dir: PathBuf,
}

impl Case {
    fn new(fasta: &str, gfa: &str) -> Self {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let case_id = CASE_ID.fetch_add(1, Ordering::Relaxed);
        let dir = std::env::temp_dir().join(format!(
            "panpath-audit-{}-{nonce}-{case_id}",
            std::process::id()
        ));
        fs::create_dir(&dir).unwrap();
        fs::write(dir.join("source.fa"), fasta).unwrap();
        fs::write(dir.join("graph.gfa"), gfa).unwrap();
        Self { dir }
    }

    fn run(&self) -> Output {
        self.run_with(&[])
    }

    fn run_with(&self, arguments: &[&str]) -> Output {
        Command::new(env!("CARGO_BIN_EXE_panpath-audit"))
            .args(arguments)
            .args([self.path("source.fa"), self.path("graph.gfa")])
            .output()
            .unwrap()
    }

    fn run_with_sources(&self, sources: &[&str], arguments: &[&str]) -> Output {
        let paths: Vec<_> = sources
            .iter()
            .enumerate()
            .map(|(index, source)| {
                let path = self.dir.join(format!("source-{index}.fa"));
                fs::write(&path, source).unwrap();
                path
            })
            .collect();
        Command::new(env!("CARGO_BIN_EXE_panpath-audit"))
            .args(arguments)
            .args(paths)
            .arg(self.path("graph.gfa"))
            .output()
            .unwrap()
    }

    fn write_source(&self, name: &str, content: &str) -> PathBuf {
        let path = self.dir.join(name);
        fs::write(&path, content).unwrap();
        path
    }

    fn gzip_source(&self, fasta: &str) {
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        encoder.write_all(fasta.as_bytes()).unwrap();
        fs::write(self.path("source.fa"), encoder.finish().unwrap()).unwrap();
    }

    fn gzip_graph(&self, gfa: &str) {
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        encoder.write_all(gfa.as_bytes()).unwrap();
        fs::write(self.path("graph.gfa"), encoder.finish().unwrap()).unwrap();
    }

    fn path(&self, name: &str) -> &Path {
        self.dir.join(name).leak()
    }
}

impl Drop for Case {
    fn drop(&mut self) {
        fs::remove_dir_all(&self.dir).unwrap();
    }
}

#[test]
fn help_exits_without_inputs() {
    for option in ["--help", "-h"] {
        let output = Command::new(env!("CARGO_BIN_EXE_panpath-audit"))
            .arg(option)
            .output()
            .unwrap();
        let stdout = String::from_utf8(output.stdout).unwrap();

        assert_eq!(output.status.code(), Some(0));
        assert!(stdout.contains("Usage:\n"));
        assert!(stdout.contains("--fasta-pansn SAMPLE HAP FILE"));
        assert!(output.stderr.is_empty());
    }
}

#[test]
fn version_exits_without_inputs() {
    for option in ["--version", "-V"] {
        let output = Command::new(env!("CARGO_BIN_EXE_panpath-audit"))
            .arg(option)
            .output()
            .unwrap();

        assert_eq!(output.status.code(), Some(0));
        assert_eq!(
            String::from_utf8(output.stdout).unwrap(),
            format!("panpath-audit {}\n", env!("CARGO_PKG_VERSION"))
        );
        assert!(output.stderr.is_empty());
    }
}

#[test]
fn double_dash_allows_paths_starting_with_hyphen() {
    let case = Case::new("", "");
    fs::write(case.path("-source.fa"), ">sample\nAC\n").unwrap();
    fs::write(case.path("-graph.gfa"), "S\t1\tAC\nP\tsample\t1+\t*\n").unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_panpath-audit"))
        .current_dir(&case.dir)
        .args(["--", "-source.fa", "-graph.gfa"])
        .output()
        .unwrap();

    assert_eq!(output.status.code(), Some(0));
}

#[test]
fn identical_path_exits_zero() {
    let case = Case::new(
        ">sampleA#1#chr1\nacgtttgacc\n",
        "H\tVN:Z:1.0\nS\t1\tACGT\nS\t2\tTTGA\nS\t3\tCC\nP\tsampleA#1#chr1\t1+,2+,3+\t*\n",
    );

    let output = case.run();

    assert_eq!(output.status.code(), Some(0));
    assert!(
        String::from_utf8(output.stdout)
            .unwrap()
            .ends_with("IDENTICAL 1\n")
    );
    assert!(output.stderr.is_empty());
}

#[test]
fn explicit_thread_count_is_reported() {
    let case = Case::new(">sample\nAC\n", "S\t1\tAC\nP\tsample\t1+\t*\n");

    let output = case.run_with(&["--format", "json", "--threads", "3"]);
    let document: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();

    assert_eq!(output.status.code(), Some(0));
    assert_eq!(document["provenance"]["threads"], 3);
}

#[test]
fn automatic_thread_count_is_capped_at_eight() {
    let case = Case::new(">sample\nAC\n", "S\t1\tAC\nP\tsample\t1+\t*\n");

    let output = case.run_with(&["--format", "json"]);
    let document: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    let threads = document["provenance"]["threads"].as_u64().unwrap();

    assert_eq!(output.status.code(), Some(0));
    assert!((1..=8).contains(&threads));
}

#[test]
fn zero_threads_is_rejected() {
    let case = Case::new(">sample\nAC\n", "S\t1\tAC\nP\tsample\t1+\t*\n");

    let output = case.run_with(&["--threads", "0"]);

    assert_eq!(output.status.code(), Some(4));
    assert!(
        String::from_utf8(output.stderr)
            .unwrap()
            .contains("threads must be a positive integer")
    );
}

#[test]
fn parallel_and_single_thread_outcomes_are_identical() {
    let case = Case::new(
        ">a\nAC\n>b\nGT\n",
        "S\t1\tAC\nS\t2\tGT\nP\ta\t1+\t*\nP\tb\t2+\t*\n",
    );

    let single = case.run_with(&["--format", "json", "--threads", "1"]);
    let parallel = case.run_with(&["--format", "json", "--threads", "4"]);
    let single: serde_json::Value = serde_json::from_slice(&single.stdout).unwrap();
    let parallel: serde_json::Value = serde_json::from_slice(&parallel.stdout).unwrap();

    assert_eq!(single["summary"], parallel["summary"]);
    assert_eq!(single["outcomes"], parallel["outcomes"]);
}

#[test]
fn multiple_fasta_inputs_form_one_source_set() {
    let case = Case::new(
        "",
        "S\t1\tAC\nS\t2\tGT\nP\tfirst\t1+\t*\nP\tsecond\t2+\t*\n",
    );

    let output = case.run_with_sources(&[">first\nAC\n", ">second\nGT\n"], &["--format", "json"]);
    let document: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();

    assert_eq!(output.status.code(), Some(0));
    assert_eq!(document["summary"]["identical"], 2);
    assert_eq!(
        document["provenance"]["source_paths"]
            .as_array()
            .unwrap()
            .len(),
        2
    );
    assert_eq!(
        document["provenance"]["fasta_compressions"],
        serde_json::json!(["plain", "plain"])
    );
}

#[test]
fn duplicate_identifiers_across_fasta_inputs_are_rejected() {
    let case = Case::new("", "S\t1\tAC\nP\tduplicate\t1+\t*\n");

    let output = case.run_with_sources(&[">duplicate\nAC\n", ">duplicate\nAC\n"], &[]);

    assert_eq!(output.status.code(), Some(4));
    assert!(
        String::from_utf8(output.stderr)
            .unwrap()
            .contains("DUPLICATE_SEQUENCE_IDENTIFIER identifier=duplicate")
    );
}

#[test]
fn pansn_and_prefix_mappings_disambiguate_reference_headers() {
    let case = Case::new(
        "",
        "H\tVN:Z:1.0\nS\t1\tAC\nS\t2\tGT\nW\tCHM13\t0\tchr1\t0\t2\t>1\nP\tGRCh38.chr1\t2+\t*\n",
    );
    let chm13 = case.write_source("chm13.fa", ">chr1\nAC\n");
    let grch38 = case.write_source("grch38.fa", ">chr1 description\nGT\n");

    let output = Command::new(env!("CARGO_BIN_EXE_panpath-audit"))
        .args(["--format", "json", "--fasta-pansn", "CHM13", "0"])
        .arg(chm13)
        .args(["--fasta-prefix", "GRCh38."])
        .arg(grch38)
        .arg(case.path("graph.gfa"))
        .output()
        .unwrap();
    let document: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();

    assert_eq!(output.status.code(), Some(0));
    assert_eq!(document["summary"]["identical"], 2);
    assert_eq!(document["sources"][0], serde_json::Value::Null);
    assert_eq!(document["provenance"]["sources"][0]["mapping"], "pansn");
    assert_eq!(document["provenance"]["sources"][1]["mapping"], "prefix");
}

#[test]
fn divergence_in_reverse_segment_reports_both_coordinates() {
    let case = Case::new(
        ">sampleB#1#chr1\nACGTGCAACC\n",
        "H\tVN:Z:1.0\nS\t1\tACGT\nS\t2\tTTGA\nS\t3\tCC\nP\tsampleB#1#chr1\t1+,2-,3+\t*\n",
    );

    let output = case.run();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(output.status.code(), Some(1));
    assert!(stdout.contains("DIVERGENT 1\n"));
    assert!(stdout.contains("sampleB#1#chr1"));
    assert!(stdout.contains("source_position=5"));
    assert!(stdout.contains("path_position=5"));
    assert!(stdout.contains("segment=2-"));
    assert!(stdout.contains("traversal_position=1"));
    assert!(stdout.contains("segment_position=4"));
    assert!(output.stderr.is_empty());
}

#[test]
fn reverse_path_supports_the_complete_iupac_alphabet() {
    let case = Case::new(
        ">iupac\nnbdhvkmwsryacgt\n",
        "S\t1\tACGTRYSWKMBDHVN\nP\tiupac\t1-\t*\n",
    );

    let output = case.run();

    assert_eq!(output.status.code(), Some(0));
    assert!(
        String::from_utf8(output.stdout)
            .unwrap()
            .ends_with("IDENTICAL 1\n")
    );
    assert!(output.stderr.is_empty());
}

#[test]
fn audit_reports_complete_correspondence_outcomes() {
    let case = Case::new(
        ">matched\nAC\n>source-only\nGT\n",
        "S\t1\tAC\nP\tmatched\t1+\t*\nP\tpath-only\t1+\t*\n",
    );

    let output = case.run();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(output.status.code(), Some(2));
    assert!(stdout.contains("IDENTICAL 1\n"));
    assert!(stdout.contains("MISSING_PATH 1\n"));
    assert!(stdout.contains("MISSING_SOURCE 1\n"));
    assert!(stdout.contains("source-only MISSING_PATH\n"));
    assert!(stdout.contains("path-only MISSING_SOURCE\n"));
    assert!(output.stderr.is_empty());
}

#[test]
fn truncated_path_reports_exhaustion_without_segment_coordinates() {
    let case = Case::new(">truncated\nACGT\n", "S\t1\tACG\nP\ttruncated\t1+\t*\n");

    let output = case.run();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(output.status.code(), Some(1));
    assert!(stdout.contains("truncated source_position=4 path_position=4"));
    assert!(stdout.contains("source_length=4 path_length=3"));
    assert!(stdout.contains("graph=<END>"));
    assert!(!stdout.contains("segment="));
    assert!(output.stderr.is_empty());
}

#[test]
fn invalid_unused_segment_blocks_the_entire_audit() {
    let case = Case::new(">valid\nAC\n", "S\t1\tAC\nP\tvalid\t1+\t*\nS\tunused\tAX\n");

    let output = case.run();
    let stderr = String::from_utf8(output.stderr).unwrap();

    assert_eq!(output.status.code(), Some(3));
    assert!(output.stdout.is_empty());
    assert!(stderr.contains("INVALID_NUCLEOTIDE"));
    assert!(stderr.contains("line=3"));
    assert!(stderr.contains("symbol=X"));
}

#[test]
fn duplicate_fasta_identifier_is_a_preflight_error() {
    let case = Case::new(
        ">dup\nAC\n>dup description\nGT\n",
        "S\t1\tAC\nP\tdup\t1+\t*\n",
    );

    let output = case.run();
    let stderr = String::from_utf8(output.stderr).unwrap();

    assert_eq!(output.status.code(), Some(4));
    assert!(output.stdout.is_empty());
    assert!(stderr.contains("DUPLICATE_SEQUENCE_IDENTIFIER identifier=dup"));
}

#[test]
fn gfa_preflight_reports_all_detectable_errors() {
    let case = Case::new(
        ">valid\nAC\n",
        "S\t1\tAC\nS\t1\tGT\nP\tvalid\t1+\t*\nP\tvalid\t1+\t*\nW\tsample\tx\tchr1\t0\t2\t>1\nS\tunused\tAX\n",
    );

    let output = case.run();
    let stderr = String::from_utf8(output.stderr).unwrap();

    assert_eq!(output.status.code(), Some(3));
    assert!(output.stdout.is_empty());
    assert!(stderr.contains("DUPLICATE_SEGMENT identifier=1"));
    assert!(stderr.contains("DUPLICATE_PATH identifier=valid"));
    assert!(stderr.contains("INVALID_HAPLOTYPE line=5 value=x"));
    assert!(stderr.contains("INVALID_NUCLEOTIDE line=6 symbol=X"));
}

#[test]
fn gfa_preflight_rejects_unreconstructable_paths() {
    let case = Case::new(
        ">valid\nAC\n",
        "S\t1\t*\nS\t2\tAC\nP\tvalid\t2+,missing+\t1M\n",
    );

    let output = case.run();
    let stderr = String::from_utf8(output.stderr).unwrap();

    assert_eq!(output.status.code(), Some(3));
    assert!(output.stdout.is_empty());
    assert!(stderr.contains("SEQUENCELESS_SEGMENT line=1 identifier=1"));
    assert!(stderr.contains("UNSUPPORTED_OVERLAP line=3 overlap=1M"));
    assert!(stderr.contains("UNRESOLVED_SEGMENT path=valid segment=missing"));
}

#[test]
fn json_mode_emits_one_versioned_document() {
    let case = Case::new(">json-path\nAC\n", "S\t1\tAC\nP\tjson-path\t1+\t*\n");

    let output = case.run_with(&["--format", "json"]);
    let document: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();

    assert_eq!(output.status.code(), Some(0));
    assert_eq!(document["schema"], "panpath-audit-report");
    assert_eq!(document["schema_version"], "0.1.0");
    assert_eq!(document["state"], "completed");
    assert_eq!(document["summary"]["identical"], 1);
    assert_eq!(document["outcomes"][0]["identifier"], "json-path");
    assert_eq!(document["outcomes"][0]["status"], "IDENTICAL");
    assert_eq!(document["provenance"]["topology_validated"], false);
    assert!(output.stderr.is_empty());
}

#[test]
fn json_mode_remains_valid_when_preflight_fails() {
    let case = Case::new(">invalid\nAC\n", "S\t1\tAX\nP\tinvalid\t1+\t*\n");

    let output = case.run_with(&["--format", "json"]);
    let document: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();

    assert_eq!(output.status.code(), Some(3));
    assert_eq!(document["schema"], "panpath-audit-report");
    assert_eq!(document["schema_version"], "0.1.0");
    assert_eq!(document["state"], "invalid");
    assert_eq!(document["outcomes"], serde_json::json!([]));
    assert_eq!(document["errors"][0]["code"], "INVALID_NUCLEOTIDE");
    assert!(output.stderr.is_empty());
}

#[test]
fn gzip_fasta_is_detected_by_content() {
    let case = Case::new(">compressed\nACGT\n", "S\t1\tACGT\nP\tcompressed\t1+\t*\n");
    case.gzip_source(">compressed\nACGT\n");

    let output = case.run();

    assert_eq!(output.status.code(), Some(0));
    assert!(
        String::from_utf8(output.stdout)
            .unwrap()
            .ends_with("IDENTICAL 1\n")
    );
    assert!(output.stderr.is_empty());
}

#[test]
fn gzip_gfa_is_streamed_and_detected_by_content() {
    let case = Case::new(">compressed\nACGT\n", "");
    case.gzip_graph("H\tVN:Z:1.0\nS\t1\tACGT\nP\tcompressed\t1+\t*\n");

    let output = case.run_with(&["--format", "json"]);
    let document: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();

    assert_eq!(output.status.code(), Some(0));
    assert_eq!(document["summary"]["identical"], 1);
    assert_eq!(document["provenance"]["gfa_compression"], "gzip");
    assert!(output.stderr.is_empty());
}

#[test]
fn full_length_walk_reconstructs_a_pansn_source() {
    let case = Case::new(
        ">sample#1#chr1\nACCTTGAGATT\n",
        "H\tVN:Z:1.1\nS\ts11\tACCTT\nS\ts12\tTC\nS\ts13\tGATT\nW\tsample\t1\tchr1\t0\t11\t>s11<s12>s13\n",
    );

    let output = case.run();

    assert_eq!(output.status.code(), Some(0));
    assert!(
        String::from_utf8(output.stdout)
            .unwrap()
            .ends_with("IDENTICAL 1\n")
    );
    assert!(output.stderr.is_empty());
}

#[test]
fn walk_ranges_are_stitched_in_coordinate_order() {
    let case = Case::new(
        ">sample#1#chr1\nACGTTCAACC\n",
        "H\tVN:Z:1.1\nS\t1\tACGT\nS\t2\tTTGA\nS\t3\tCC\nW\tsample\t1\tchr1\t4\t10\t<2>3\nW\tsample\t1\tchr1\t0\t4\t>1\n",
    );

    let output = case.run();

    assert_eq!(output.status.code(), Some(0));
    assert!(
        String::from_utf8(output.stdout)
            .unwrap()
            .ends_with("IDENTICAL 1\n")
    );
    assert!(output.stderr.is_empty());
}

#[test]
fn clipped_gaps_between_walk_ranges_are_ignored() {
    let case = Case::new(
        ">sample#1#chr1\nACGTGT\n",
        "H\tVN:Z:1.1\nS\t1\tAC\nS\t2\tGT\nW\tsample\t1\tchr1\t0\t2\t>1\nW\tsample\t1\tchr1\t4\t6\t>2\n",
    );

    let output = case.run();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(output.status.code(), Some(0));
    assert!(stdout.ends_with("IDENTICAL 1\n"));
    assert!(output.stderr.is_empty());
}

#[test]
fn overlapping_walk_ranges_are_a_preflight_error() {
    let case = Case::new(
        ">sample#1#chr1\nACGTA\n",
        "H\tVN:Z:1.1\nS\t1\tACGT\nS\t2\tTA\nW\tsample\t1\tchr1\t0\t4\t>1\nW\tsample\t1\tchr1\t3\t5\t>2\n",
    );

    let output = case.run();
    let stderr = String::from_utf8(output.stderr).unwrap();

    assert_eq!(output.status.code(), Some(2));
    assert!(output.stdout.is_empty());
    assert!(stderr.contains("W_RANGE_OVERLAP identifier=sample#1#chr1"));
}

#[test]
fn path_and_walk_cannot_claim_the_same_identifier() {
    let case = Case::new(
        ">sample#1#chr1\nAC\n",
        "H\tVN:Z:1.1\nS\t1\tAC\nP\tsample#1#chr1\t1+\t*\nW\tsample\t1\tchr1\t0\t2\t>1\n",
    );

    let output = case.run();
    let stderr = String::from_utf8(output.stderr).unwrap();

    assert_eq!(output.status.code(), Some(2));
    assert!(output.stdout.is_empty());
    assert!(stderr.contains("DUPLICATE_EMBEDDED_PATH identifier=sample#1#chr1"));
}

#[test]
fn walk_is_accepted_with_legacy_gfa_1_0_header() {
    let case = Case::new(
        ">sample#1#chr1\nAC\n",
        "H\tVN:Z:1.0\nS\t1\tAC\nW\tsample\t1\tchr1\t0\t2\t>1\n",
    );

    let output = case.run();
    assert_eq!(output.status.code(), Some(0));
    assert!(
        String::from_utf8(output.stdout)
            .unwrap()
            .ends_with("IDENTICAL 1\n")
    );
    assert!(output.stderr.is_empty());
}

#[test]
fn walk_range_cannot_exceed_the_source_sequence() {
    let case = Case::new(
        ">sample#1#chr1\nACGT\n",
        "H\tVN:Z:1.1\nS\t1\tACGT\nW\tsample\t1\tchr1\t0\t5\t>1\n",
    );

    let output = case.run();
    let stderr = String::from_utf8(output.stderr).unwrap();

    assert_eq!(output.status.code(), Some(3));
    assert!(output.stdout.is_empty());
    assert!(
        stderr.contains("W_RANGE_OUT_OF_BOUNDS identifier=sample#1#chr1 end=5 source_length=4")
    );
}

#[test]
fn walk_length_must_match_its_declared_range() {
    let case = Case::new(
        ">sample#1#chr1\nACGTCC\n",
        "H\tVN:Z:1.1\nS\t1\tACGTC\nS\t2\tCC\nW\tsample\t1\tchr1\t0\t4\t>1\nW\tsample\t1\tchr1\t4\t6\t>2\n",
    );

    let output = case.run();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(output.status.code(), Some(1));
    assert!(stdout.contains("source_position=5"));
    assert!(stdout.contains("range_length=4 walk_length=5"));
    assert!(output.stderr.is_empty());
}

#[test]
fn unspecified_walk_range_cannot_be_combined_with_another_record() {
    let case = Case::new(
        ">sample#1#chr1\nACGT\n",
        "H\tVN:Z:1.1\nS\t1\tACGT\nS\t2\tGT\nW\tsample\t1\tchr1\t*\t*\t>1\nW\tsample\t1\tchr1\t2\t4\t>2\n",
    );

    let output = case.run();
    let stderr = String::from_utf8(output.stderr).unwrap();

    assert_eq!(output.status.code(), Some(2));
    assert!(output.stdout.is_empty());
    assert!(stderr.contains("AMBIGUOUS_W_RANGES identifier=sample#1#chr1"));
}

#[test]
fn walk_haplotype_index_must_be_numeric() {
    let case = Case::new(
        ">sample#x#chr1\nAC\n",
        "H\tVN:Z:1.1\nS\t1\tAC\nW\tsample\tx\tchr1\t0\t2\t>1\n",
    );

    let output = case.run();
    let stderr = String::from_utf8(output.stderr).unwrap();

    assert_eq!(output.status.code(), Some(3));
    assert!(output.stdout.is_empty());
    assert!(stderr.contains("INVALID_HAPLOTYPE line=3 value=x"));
}

#[test]
fn json_walk_divergence_reports_range_and_segment_coordinates() {
    let case = Case::new(
        ">sample#1#chr1\nACGTGCAA\n",
        "H\tVN:Z:1.1\nS\t1\tACGT\nS\t2\tTTGA\nW\tsample\t1\tchr1\t0\t8\t>1<2\n",
    );

    let output = case.run_with(&["--format", "json"]);
    let document: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    let outcome = &document["outcomes"][0];

    assert_eq!(output.status.code(), Some(1));
    assert_eq!(outcome["record_type"], "W");
    assert_eq!(outcome["source_range"]["start_0based"], 0);
    assert_eq!(outcome["source_range"]["end_0based_exclusive"], 8);
    assert_eq!(outcome["divergence"]["segment"], "2");
    assert_eq!(outcome["divergence"]["orientation"], "-");
    assert_eq!(outcome["divergence"]["traversal_position_1based"], 1);
    assert_eq!(outcome["divergence"]["segment_position_1based"], 4);
    assert!(output.stderr.is_empty());
}

#[test]
fn json_provenance_counts_path_and_walk_records() {
    let case = Case::new(
        ">sample#1#chr1\nAC\n>named-path\nGT\n",
        "H\tVN:Z:1.1\nS\t1\tAC\nS\t2\tGT\nW\tsample\t1\tchr1\t*\t*\t>1\nP\tnamed-path\t2+\t*\n",
    );

    let output = case.run_with(&["--format", "json"]);
    let document: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    let provenance = &document["provenance"];

    assert_eq!(output.status.code(), Some(0));
    assert_eq!(provenance["path_records"], 1);
    assert_eq!(provenance["walk_records"], 1);
    assert_eq!(provenance["record_types"], serde_json::json!(["P", "W"]));
    assert!(output.stderr.is_empty());
}

#[test]
fn json_identical_outcomes_report_their_record_type() {
    let case = Case::new(
        ">sample#1#chr1\nAC\n>named-path\nGT\n",
        "H\tVN:Z:1.1\nS\t1\tAC\nS\t2\tGT\nW\tsample\t1\tchr1\t0\t2\t>1\nP\tnamed-path\t2+\t*\n",
    );

    let output = case.run_with(&["--format", "json"]);
    let document: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    let outcomes = document["outcomes"].as_array().unwrap();
    let path = outcomes
        .iter()
        .find(|outcome| outcome["identifier"] == "named-path")
        .unwrap();
    let walk = outcomes
        .iter()
        .find(|outcome| outcome["identifier"] == "sample#1#chr1")
        .unwrap();

    assert_eq!(output.status.code(), Some(0));
    assert_eq!(path["record_type"], "P");
    assert_eq!(walk["record_type"], "W");
    assert!(output.stderr.is_empty());
}

#[test]
fn missing_source_walk_retains_its_record_type_in_json() {
    let case = Case::new(
        ">present\nGT\n",
        "H\tVN:Z:1.1\nS\t1\tAC\nS\t2\tGT\nW\tsample\t1\tchr1\t0\t2\t>1\nP\tpresent\t2+\t*\n",
    );

    let output = case.run_with(&["--format", "json"]);
    let document: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    let outcome = document["outcomes"]
        .as_array()
        .unwrap()
        .iter()
        .find(|outcome| outcome["identifier"] == "sample#1#chr1")
        .unwrap();

    assert_eq!(output.status.code(), Some(2));
    assert_eq!(outcome["status"], "MISSING_SOURCE");
    assert_eq!(outcome["record_type"], "W");
    assert!(output.stderr.is_empty());
}

#[test]
fn json_clipped_walk_reports_only_embedded_length() {
    let case = Case::new(
        ">sample#1#chr1\nACGTGT\n",
        "H\tVN:Z:1.1\nS\t1\tAC\nS\t2\tGT\nW\tsample\t1\tchr1\t0\t2\t>1\nW\tsample\t1\tchr1\t4\t6\t>2\n",
    );

    let output = case.run_with(&["--format", "json"]);
    let document: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(output.status.code(), Some(0));
    assert_eq!(document["outcomes"][0]["status"], "IDENTICAL");
    assert_eq!(document["outcomes"][0]["source_length"], 6);
    assert_eq!(document["outcomes"][0]["path_length"], 4);
    assert!(output.stderr.is_empty());
}

#[test]
fn json_walk_range_length_mismatch_is_structured() {
    let case = Case::new(
        ">sample#1#chr1\nACGTCC\n",
        "H\tVN:Z:1.1\nS\t1\tACGTC\nS\t2\tCC\nW\tsample\t1\tchr1\t0\t4\t>1\nW\tsample\t1\tchr1\t4\t6\t>2\n",
    );

    let output = case.run_with(&["--format", "json"]);
    let document: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    let divergence = &document["outcomes"][0]["divergence"];

    assert_eq!(output.status.code(), Some(1));
    assert_eq!(divergence["kind"], "RANGE_LENGTH");
    assert_eq!(divergence["range_length"], 4);
    assert_eq!(divergence["walk_length"], 5);
    assert!(output.stderr.is_empty());
}

#[test]
fn truncated_walk_record_is_a_preflight_error() {
    let case = Case::new(
        ">sample#1#chr1\nAC\n",
        "H\tVN:Z:1.1\nS\t1\tAC\nW\tsample\t1\tchr1\t0\t2\n",
    );

    let output = case.run();
    let stderr = String::from_utf8(output.stderr).unwrap();

    assert_eq!(output.status.code(), Some(3));
    assert!(output.stdout.is_empty());
    assert!(stderr.contains("MALFORMED_W_RECORD line=3"));
}

#[test]
fn comprehensive_statistics_partition_source_and_graph_bases() {
    let case = Case::new(
        ">same\nACGT\n>changed\nAAAA\n>absent\nNN\n",
        "S\t1\tACGT\nS\t2\tAATA\nS\t3\tGG\nP\tsame\t1+\t*\nP\tchanged\t2+\t*\nP\textra\t3+\t*\n",
    );
    let output = case.run_with(&["--format", "json", "--stats", "comprehensive"]);
    let document: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();

    assert_eq!(output.status.code(), Some(2));
    assert_eq!(document["schema_version"], "0.1.0");
    assert_eq!(document["statistics"]["source"]["total_bases"], 10);
    assert_eq!(document["statistics"]["source"]["bases"]["matched"], 7);
    assert_eq!(document["statistics"]["source"]["bases"]["substituted"], 1);
    assert_eq!(document["statistics"]["source"]["bases"]["missing_path"], 2);
    assert_eq!(
        document["statistics"]["source"]["exact_recovery_percent"],
        70.0
    );
    assert_eq!(document["statistics"]["graph"]["total_bases"], 10);
    assert_eq!(
        document["statistics"]["graph"]["bases"]["missing_source"],
        2
    );
    assert_eq!(document["statistics"]["alignment"]["edit_distance"], 1);
    assert_eq!(document["statistics"]["ambiguity"]["source_bases"], 2);
}

#[test]
fn comprehensive_statistics_count_indels() {
    let case = Case::new(">changed\nACGT\n", "S\t1\tAGT\nP\tchanged\t1+\t*\n");
    let output = case.run_with(&["--format", "json", "--stats", "comprehensive"]);
    let document: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    let outcome = &document["outcomes"][0];

    assert_eq!(outcome["base_statistics"]["source"]["matched"], 3);
    assert_eq!(outcome["base_statistics"]["source"]["source_only"], 1);
    assert_eq!(outcome["base_statistics"]["graph"]["graph_only"], 0);
    assert_eq!(outcome["length_delta"], -1);
}

#[test]
fn comprehensive_statistics_count_clipped_walk_bases() {
    let case = Case::new(
        ">sample#1#chr1\nACGTGT\n",
        "H\tVN:Z:1.1\nS\t1\tAC\nS\t2\tGT\nW\tsample\t1\tchr1\t0\t2\t>1\nW\tsample\t1\tchr1\t4\t6\t>2\n",
    );
    let output = case.run_with(&["--format", "json", "--stats", "comprehensive"]);
    let document: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();

    assert_eq!(output.status.code(), Some(0));
    assert_eq!(document["statistics"]["source"]["bases"]["matched"], 4);
    assert_eq!(document["statistics"]["source"]["bases"]["not_embedded"], 2);
    assert_eq!(
        document["statistics"]["source"]["embedded_coverage_percent"],
        66.666667
    );
}

#[test]
fn alignment_limit_uses_unavailable_buckets() {
    let case = Case::new(">changed\nAAAA\n", "S\t1\tTTTT\nP\tchanged\t1+\t*\n");
    let output = case.run_with(&[
        "--format",
        "json",
        "--stats",
        "comprehensive",
        "--alignment-max-cells",
        "1",
    ]);
    let document: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    let outcome = &document["outcomes"][0];

    assert_eq!(outcome["alignment_status"], "limit_exceeded");
    assert_eq!(outcome["edit_distance"], serde_json::Value::Null);
    assert_eq!(
        outcome["base_statistics"]["source"]["unaligned_divergent"],
        4
    );
    assert_eq!(
        outcome["base_statistics"]["graph"]["unaligned_divergent"],
        4
    );
}

#[test]
fn comprehensive_invalid_json_has_null_statistics() {
    let case = Case::new(">invalid\nAC\n", "S\t1\tAX\nP\tinvalid\t1+\t*\n");
    let output = case.run_with(&["--format", "json", "--stats", "comprehensive"]);
    let document: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();

    assert_eq!(output.status.code(), Some(3));
    assert_eq!(document["schema_version"], "0.1.0");
    assert_eq!(document["statistics"], serde_json::Value::Null);
}

#[test]
fn missing_source_lengths_support_forward_and_repeated_segments() {
    let case = Case::new(
        ">present\nAC\n",
        "P\textra\t1+,1+\t*\nS\t1\tAN\nP\tpresent\t1+\t*\n",
    );
    let output = case.run_with(&["--format", "json", "--stats", "comprehensive"]);
    let document: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    let extra = document["outcomes"]
        .as_array()
        .unwrap()
        .iter()
        .find(|outcome| outcome["identifier"] == "extra")
        .unwrap();

    assert_eq!(extra["base_statistics"]["graph"]["missing_source"], 4);
    assert_eq!(extra["ambiguity"]["graph_bases"], 2);
}

#[test]
fn comprehensive_statistics_are_parallel_deterministic() {
    let case = Case::new(
        ">a\nACGT\n>b\nAAAA\n",
        "S\t1\tACGT\nS\t2\tAATA\nP\ta\t1+\t*\nP\tb\t2+\t*\n",
    );
    let single = case.run_with(&[
        "--format",
        "json",
        "--stats",
        "comprehensive",
        "--threads",
        "1",
    ]);
    let parallel = case.run_with(&[
        "--format",
        "json",
        "--stats",
        "comprehensive",
        "--threads",
        "4",
    ]);
    let mut single: serde_json::Value = serde_json::from_slice(&single.stdout).unwrap();
    let mut parallel: serde_json::Value = serde_json::from_slice(&parallel.stdout).unwrap();
    single["provenance"]["threads"] = 0.into();
    parallel["provenance"]["threads"] = 0.into();

    assert_eq!(single, parallel);
}

#[test]
fn alignment_limit_requires_comprehensive_mode() {
    let case = Case::new(">a\nAC\n", "S\t1\tAC\nP\ta\t1+\t*\n");
    let output = case.run_with(&["--alignment-max-cells", "10"]);

    assert_eq!(output.status.code(), Some(4));
    assert!(
        String::from_utf8(output.stderr)
            .unwrap()
            .contains("alignment max cells requires --stats comprehensive")
    );
}

#[test]
fn memory_budget_is_reported_and_limits_concurrency() {
    let sequence = "A".repeat(700_000);
    let source = format!(">a\n{sequence}\n>b\n{sequence}\n");
    let graph = format!("S\t1\t{sequence}\nP\ta\t1+\t*\nP\tb\t1+\t*\n");
    let case = Case::new(&source, &graph);
    let output = case.run_with(&["--format", "json", "--memory-mib", "1", "--threads", "4"]);
    let document: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();

    assert_eq!(output.status.code(), Some(0));
    assert_eq!(document["provenance"]["memory_mib"], 1);
    assert_eq!(document["provenance"]["peak_tracked_bytes"], 700_000);
    assert_eq!(document["provenance"]["oversized_contigs"], 0);
}

#[test]
fn clipped_walk_reports_full_and_addressed_digests() {
    let case = Case::new(
        ">sample#1#chr1\nACGTGT\n",
        "S\t1\tAC\nS\t2\tGT\nW\tsample\t1\tchr1\t0\t2\t>1\nW\tsample\t1\tchr1\t4\t6\t>2\n",
    );
    let output = case.run_with(&["--format", "json"]);
    let document: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    let digests = &document["outcomes"][0]["digests"];

    assert_eq!(output.status.code(), Some(0));
    assert_ne!(
        digests["sha256"]["source"],
        digests["sha256"]["addressed_source"]
    );
    assert_eq!(
        digests["sha256"]["addressed_source"],
        digests["sha256"]["graph"]
    );
    assert_ne!(
        digests["blake3"]["source"],
        digests["blake3"]["addressed_source"]
    );
    assert_eq!(
        digests["blake3"]["addressed_source"],
        digests["blake3"]["graph"]
    );
}

#[test]
fn tsv_reports_outcomes_and_digests() {
    let case = Case::new(">sample\nac\n", "S\t1\tAC\nP\tsample\t1+\t*\n");
    let output = case.run_with(&["--format", "tsv"]);
    let stdout = String::from_utf8(output.stdout).unwrap();
    let lines: Vec<_> = stdout.lines().collect();

    assert_eq!(output.status.code(), Some(0));
    assert!(lines[0].contains("addressed_source_sha256"));
    assert!(lines[0].contains("graph_blake3"));
    assert!(lines[1].starts_with("outcome\tsample\tIDENTICAL\tP\t"));

    let fields: std::collections::HashMap<_, _> =
        lines[0].split('\t').zip(lines[1].split('\t')).collect();
    assert_eq!(fields["identifier"], "sample");
    assert_eq!(fields["source_length"], "2");
    assert_eq!(fields["graph_length"], "2");
    assert_eq!(fields["source_sha256"], fields["graph_sha256"]);
    assert_eq!(fields["source_blake3"], fields["graph_blake3"]);
    assert_eq!(fields["topology_validated"], "false");
    assert_eq!(fields["unspecified_overlaps_as_blunt"], "true");
}

#[test]
fn invalid_tsv_uses_the_outcome_header() {
    let case = Case::new(">sample\nAC\n", "S\t1\tAX\nP\tsample\t1+\t*\n");
    let output = case.run_with(&["--format", "tsv"]);
    let stdout = String::from_utf8(output.stdout).unwrap();
    let lines: Vec<_> = stdout.lines().collect();

    assert_eq!(output.status.code(), Some(3));
    assert!(lines[0].contains("addressed_source_sha256"));
    assert!(lines[1].starts_with("error\t\t\t\tINVALID_NUCLEOTIDE\t"));
    assert_eq!(lines[0].split('\t').count(), lines[1].split('\t').count());
}

#[test]
fn human_output_discloses_provenance_limits() {
    let case = Case::new(">sample\nAC\n", "S\t1\tAC\nP\tsample\t1+\t*\n");
    let output = case.run();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert!(
        stdout.contains("PROVENANCE topology_validated=false unspecified_overlaps_as_blunt=true\n")
    );
}

#[test]
fn preflight_combines_fasta_and_gfa_errors() {
    let case = Case::new(
        ">dup\nAC\n>dup\nGT\n",
        "S\t1\tAX\nP\tdup\tbroken\t*\nS\t2\tAZ\n",
    );
    let output = case.run();
    let stderr = String::from_utf8(output.stderr).unwrap();

    assert_eq!(output.status.code(), Some(4));
    assert!(stderr.contains("DUPLICATE_SEQUENCE_IDENTIFIER"));
    assert!(stderr.contains("MALFORMED_PATH"));
    assert!(stderr.contains("INVALID_NUCLEOTIDE line=3"));
}

#[test]
fn invalid_source_symbols_are_preflight_errors() {
    let case = Case::new(
        ">bad-one\nAC\nGX\n>valid\nryswkmbdhvn\n>bad-two\nA?\n",
        "S\t1\tACGX\nP\tbad-one\t1+\t*\n",
    );
    let output = case.run_with(&["--format", "json"]);
    let document: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    let errors = document["errors"].as_array().unwrap();

    assert_eq!(output.status.code(), Some(4));
    assert!(errors.iter().any(|error| {
        error["code"] == "INVALID_SOURCE_NUCLEOTIDE"
            && error["message"]
                .as_str()
                .unwrap()
                .contains("identifier=bad-one position=4 symbol=X")
    }));
    assert!(errors.iter().any(|error| {
        error["code"] == "INVALID_SOURCE_NUCLEOTIDE"
            && error["message"]
                .as_str()
                .unwrap()
                .contains("identifier=bad-two position=2 symbol=?")
    }));
    assert_eq!(
        errors
            .iter()
            .filter(|error| error["code"] == "INVALID_SOURCE_NUCLEOTIDE")
            .count(),
        2
    );
    assert_eq!(document["outcomes"], serde_json::json!([]));
}

#[test]
fn invalid_gzip_source_symbol_is_a_source_error() {
    let case = Case::new(">bad\nACX\n", "S\t1\tACG\nP\tbad\t1+\t*\n");
    case.gzip_source(">bad\nACX\n");
    let output = case.run();
    let stderr = String::from_utf8(output.stderr).unwrap();

    assert_eq!(output.status.code(), Some(4));
    assert!(stderr.contains("INVALID_SOURCE_NUCLEOTIDE"));
    assert!(stderr.contains("identifier=bad position=3 symbol=X"));
}

#[test]
fn anonymous_scratch_files_leave_no_directory_entries() {
    let case = Case::new(">sample\nAC\n", "S\t1\tAC\nP\tsample\t1+\t*\n");
    let before = fs::read_dir(&case.dir).unwrap().count();
    let temp_dir = case.dir.to_string_lossy().into_owned();
    let output = case.run_with(&["--temp-dir", &temp_dir]);
    let after = fs::read_dir(&case.dir).unwrap().count();

    assert_eq!(output.status.code(), Some(0));
    assert_eq!(before, after);
}

#[test]
fn oversized_contig_runs_alone_and_is_reported() {
    let sequence = "A".repeat(1_100_000);
    let source = format!(">large\n{sequence}\n");
    let graph = format!("S\t1\t{sequence}\nP\tlarge\t1+\t*\n");
    let case = Case::new(&source, &graph);
    let output = case.run_with(&["--format", "json", "--memory-mib", "1"]);
    let document: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();

    assert_eq!(output.status.code(), Some(0));
    assert_eq!(document["provenance"]["peak_tracked_bytes"], 1_100_000);
    assert_eq!(document["provenance"]["oversized_contigs"], 1);
}

#[test]
fn alignment_respects_memory_budget() {
    let source_sequence = "A".repeat(600_000);
    let graph_sequence = "T".repeat(600_000);
    let source = format!(">changed\n{source_sequence}\n");
    let graph = format!("S\t1\t{graph_sequence}\nP\tchanged\t1+\t*\n");
    let case = Case::new(&source, &graph);
    let output = case.run_with(&[
        "--format",
        "json",
        "--stats",
        "comprehensive",
        "--memory-mib",
        "1",
    ]);
    let document: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();

    assert_eq!(output.status.code(), Some(1));
    assert_eq!(
        document["outcomes"][0]["alignment_status"],
        "memory_limit_exceeded"
    );
    assert_eq!(
        document["statistics"]["alignment"]["unavailable_traversals"],
        1
    );
    assert_eq!(
        document["outcomes"][0]["base_statistics"]["source"]["unaligned_divergent"],
        600_000
    );
}

#[test]
fn zero_memory_budget_is_rejected() {
    let case = Case::new(">a\nAC\n", "S\t1\tAC\nP\ta\t1+\t*\n");
    let output = case.run_with(&["--memory-mib", "0"]);

    assert_eq!(output.status.code(), Some(4));
    assert!(
        String::from_utf8(output.stderr)
            .unwrap()
            .contains("memory must be a positive integer in MiB")
    );
}
