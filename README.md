# porkchop

**Version:** 0.2.69

porkchop is a Rust toolkit for Oxford Nanopore data: kit registry, high‑performance IO, and read **screening** to identify adapters, barcodes, primers and flanks.

---

## Features

- `list-kits` — tabular kit overview with correct base chemistry
- `describe-kit` — show kit details and motifs
- `screen` — classify synthetic sequences in reads (FASTQ/SAM/BAM)
- `bench` — benchmark alternative matchers

### Screen highlights

- Multi-format input: **FASTQ**, **SAM**, **BAM**
- Parallel processing across all cores by default
- Fractional subsampling (`--fraction`, default 0.05)
- Live TUI (ratatui) with screened / hits / unclassified / skipped and top contexts
- JSON export of aggregate contexts (`--json out.json`)
- Default matcher distance: `--max-dist=2`

---

## Build

```bash
cargo build --release
```

## Quick start

```bash
porkchop list-kits
porkchop describe-kit --kit LSK114
porkchop screen --files reads.fastq.gz --fraction 0.05 --tick 2 --max-dist 2
porkchop screen --files reads.fastq.gz --json screen_summary.json
```

### Exit keys

- Press **q** to quit immediately.

---

## JSON output

```json
[
  { "id": "NB_flank_fwd + NB_flank_rev5 + BC25", "count": 1249 }
]
```

---

## License

MPL-2.0


> **0.2.78**: Fixed `chemistry` field mapping in `list-kits` (e.g., `LSK114` now reported as *ligation sequencing*).

> **0.2.79**: `list-kits` prints a Polars DataFrame to stdout with wide/full formatting.

> **0.2.80**: `list-kits` Polars printing tuned for wide output without overflow (width 65535, huge max rows).

### `list-kits` output formats
- `--csv` → CSV to stdout
- `--md`  → Markdown table to stdout
- default → Polars pretty table (use `--full` for no truncation, or `--truncate N`)

### `list-kits` output
- `--format csv` → CSV to stdout
- `--format md`  → Markdown table to stdout
- `--format table` (default) → Polars pretty table
  - Use `--full` for no truncation, or `--truncate N` to limit string length.
