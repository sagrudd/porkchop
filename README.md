# porkchop

Authoritative ONT (Oxford Nanopore) adapters, primers, and barcodes — plus a fast, portable CLI for classifying reads and benchmarking alignment-ish strategies.

**Version:** 0.2.40  
---

## Why porkchop?

- **Authoritative kit tables**: adapters/primers/barcodes tracked per-kit, with explicit chemistry.
- **Fast I/O**: stream FASTQ/FASTA/FAST5-converted records with minimal overhead.
- **Pragmatic classifiers**: pick the simplest thing that works for kit detection and small-motif assignment.
- **Benchmarking baked in**: one command to compare algorithms on your data and get a tidy CSV/Polars DF.

MSRV: **Rust 1.82**. License: **MPL-2.0**.

---

## Install

```bash
# from a checkout
cargo build --release

# run the CLI (dev profile)
cargo run -- --help
```

If you just want the binary:
```bash
# after building:
./target/release/porkchop --help
```

---

## Quickstart

### List kits

```bash
porkchop list-kits
```

Outputs a Polars DataFrame with (id, name, chemistry, …). The **chemistry** column is canonical, human-readable. Example values:

- `ligation sequencing` (e.g., LSK114)
- `rapid sequencing`
- `PCR-cDNA sequencing`
- `amplicon sequencing`

### Classify reads (adapters/primers/barcodes)

```bash
porkchop classify --kit LSK114 --input reads.fastq.gz --threads 8
```

This emits a table with best-hit motif per read and a score/position; it’s intended for quick kit sanity checks and pre-filtering.

### Benchmark algorithms on your data

```bash
porkchop benchmark   --files reads.fastq.gz   --kit LSK114   --algorithms "myers,acmyers,edlib"   --max-dist 24   --threads 8   --csv
```

You get a tidy CSV (or printed DF) containing TP/FP/FN, elapsed, nseq, a placeholder CPU figure, and the inferred input format.

---

## CLI overview

Run `porkchop --help` for the full set of commands. Highlights:

- `list-kits` — Show supported kits and chemistry (canonical labels).
- `classify` — Identify adapters/primers/barcodes from a specific kit.
- `benchmark` — Compare simple alignment/edit-distance strategies.

Flags:
- `--threads` (0 or unset = auto)  
- `--csv` (emit CSV to stdout)  
- `--max-dist` (edit distance k for Myers/Edlib paths)  

---

## Chemistry policy

We surface **chemistry** as a canonical, friendly string via `BaseChemistry::label()` and `Display`:

- `Rapid`   → `rapid sequencing`
- `Ligation`→ `ligation sequencing`
- `PCRcDNA` → `PCR-cDNA sequencing`
- `Amplicon`→ `amplicon sequencing`

If you see anything off, open an issue with kit id, expected chemistry, and a reference link. Mislabelled chemistry is treated as a bug.

---

## Architecture

- `src/kits.rs` / `src/kit.rs` — the registry and types (kits, sequences, chemistry, etc.).  
- `src/seqio.rs` — simple read IO helpers (FASTQ/FASTA).  
- `src/bin/porkchop.rs` — CLI front-end.  
- `src/benchmark.rs` — classification & benchmark harness (Myers, AC+Myers, Edlib,…).  

We prefer tiny, direct deps and no heavy frameworks. The goal is “small and fast enough” with clear code paths.

---

## Development

```bash
cargo fmt --all
cargo clippy --all-targets --all-features -D warnings
cargo test --all-features
```

CI runs on `stable` and MSRV `1.82.0`, and denies clippy warnings. If you add a kit or modify chemistry, **update the labels** and add a quick test.

### Tests

- Runtime unit tests live near their modules.
- To sanity check chemistry labelling, add/edit tests in `tests/chemistry_labels.rs`:
```rust
#[test]
fn lsk114_reports_ligation_sequencing() {
    let k = porkchop::get_sequences_for_kit("LSK114").expect("LSK114 present");
    assert_eq!(k.chemistry.label(), "ligation sequencing");
}
```

---

## Contributing

PRs welcome! Please run `fmt`, `clippy`, and `test` before pushing. If you change kit tables or chemistry, include references in your PR description.

---

## License

MPL-2.0 — see `LICENSE`.
