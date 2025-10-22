# porkchop

**Version:** 0.2.5  
**License:** MPL-2.0

This package combines the authoritative Oxford Nanopore **kit registry** with a high-performance IO layer (FASTQ/FASTQ.GZ/SAM/BAM) and a **benchmark** framework to evaluate adapter/primer/barcode classifiers.

> This archive builds with **placeholder kit data** so you can compile and run the CLI immediately.
> In your environment we'll merge the full registry from the upstream authoritative repo:
> `https://github.com/sagrudd/porkchop`.

## CLI

```bash
porkchop list-kits
porkchop describe-kit --kit NBD114.24
porkchop benchmark --kit NBD114.24 --algorithms aho,myers,edlib,parasail,ac+myers,ac+parasail reads.fastq.gz
```

The `describe-kit` output is **not truncated** and includes provenance columns (`provenance_source`, `provenance_appendix`).

## Truth set (benchmark)
CSV/TSV with headers: `read_id, expected_labels, [kind]` where `expected_labels` are `;`-separated label symbols.

## Notes
- For a patch release from 0.2.0, we preserved all benchmark code and CLI. Only the crate version and docs were bumped.


---

## Minimum Supported Rust Version (MSRV)

- **Rust 1.90** (declared via `rust-version = "1.90"`).

## Build

```bash
cargo build --release
```

## CLI usage (verbose)

### List kits
```bash
porkchop list-kits
```

### Describe a kit (no truncation; provenance columns included)
```bash
porkchop describe-kit --kit NBD114.24
# Add --csv to emit machine‑readable output
```

### Benchmark classifiers
```bash
# With a truth set
porkchop benchmark --kit NBD114.24 --truth truth.csv \
  --algorithms aho,myers,edlib,parasail,ac+myers,ac+parasail \
  --threads 16 reads1.fastq.gz reads2.bam > results.tsv

# Without a truth set (timing/throughput + CPU utilisation only)
porkchop benchmark --kit PCS114 reads.fastq.gz --threads 32
```

#### Truth set schema
A CSV/TSV with headers:
- `read_id` — the read identifier
- `expected_labels` — semicolon‑separated label symbols (e.g., `NB01;SSPII`)
- `kind` *(optional)* — one of `adapter`, `primer`, `barcode`; restricts evaluation

#### Performance tuning
- `--threads` (default: all logical cores) controls Rayon and BAM reader threads.
- Classifiers using Myers/Parasail **prebuild** automata/profiles once per run to reduce per‑read overhead.
- Keep the callback in `seqio::for_each_parallel` lightweight; aggregate metrics outside the closure when possible.

## Provenance

Each sequence record carries a `provenance_source` (URL/tag) and optional `provenance_appendix` to track origin.
