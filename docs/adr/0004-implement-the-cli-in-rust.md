# Implement the CLI in Rust

Panpath Audit is implemented in Rust and distributed as a self-contained executable, while allowing audited and pinned Rust crate dependencies. C++ was a credible alternative with comparable performance and stronger integration into some bioinformatics libraries, but v0.1 requires no such integration; Rust better supports safe parsing, checked coordinate bookkeeping, explicit validated states, and standalone deployment for a correctness verifier handling large or malformed inputs.
