# porkchop

**Authoritative reference** of Oxford Nanopore adapters, primers and barcodes, mapped to kit identifiers, with a fuzzy matcher to infer kit-of-origin.

- License: MPL-2.0
- Version: 0.1.3
- Stability: public API is additive within 0.1.x
- No optional features: all functionality is enabled by default.

## Install
```toml
porkchop = "0.1"
```

## Kit coverage

- PCS111 (PCR‑cDNA Sequencing)
- PCS114 (PCR‑cDNA Sequencing V14)

**Current (Kit 14)**

- LSK114 / LSK114‑XL — LA adapter.
- NBD114.24 / NBD114.96 — NA adapter; NB01–NB96 (+ NB flanks).
- RBK114.24 / RBK114.96 — RA adapter; RB01–RB96 core barcodes (+ RB flanks).
- PCB114.24 — RA + cDNA adapters/primers (RTP, SSPII, CRTA, cPRM); BP01–BP24 (+ PCB flanks).
- RPB114.24 — RA + RLB01–RLB24 core barcodes (+ RPB flank).
- 16S114.24 — RA + 16S01–16S24 core barcodes (+ 16S primer flanks).

**Expansions**

- EXP‑PBC001 (PBC001) — BC01–BC12 (+ PCB flanks).
- EXP‑PBC096 (PBC096) — BC01–BC96 (+ PCB flanks).

**Legacy**

- RBK004 — RB01–RB12 (+ RB flank).
- RBK110.96 — RB01–RB96 (+ RB flank).
- LSK108 / LSK109 — Y‑adapter trunk fragments (from Porechop forks).
- LSK308 (1D^2) — 1D^2 adapter fragments (from Porechop forks).

## API (high level)
```rust
// Retrieve a kit definition (adapters/primers + barcodes).
let k = porkchop::get_sequences_for_kit("NBD114.96").unwrap();

// Scan a read for known motifs (adapters/barcodes/flanks).
let hits = porkchop::infer_kit_from_sequence("...TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCT...", 2, None);

// Rank kits by evidence.
let ranking = porkchop::infer_kits_from_sequence("...AAGGTTAA" , 2);
```

## Provenance of sequences
- **Chemistry Technical Document** (CHTD_500_v1_revAR_25Nov2024), Appendix 14 (barcodes & flanks), Appendix 15 (adapters and primers). We transcribed sequences verbatim.
- **Porechop** forks (e.g. Sn0flingan/Poresnip `porechop/adapters.py`) for legacy adapter fragments (1D^2 / Y‑adapter trunks). We include only raw strings (no GPL code).
- **qcat** (archived) used to cross‑check legacy kit naming/coverage.


## MSRV
Rust 1.70.

## License
MPL-2.0


## CLI (binary)

This crate ships a small CLI, `porkchop`, built with **clap**.

### Install & run
```bash
# From the repo root:
cargo run -- list-kits

# After installing the binary:
cargo install --path .
porkchop list-kits
porkchop list-kits --csv > kits.csv
porkchop list-kits --markdown
```

### `list-kits`
Builds a tiny **Polars** table of all supported kits with:
- `kit_id`
- `legacy` (bool)
- `base_chemistry` ("rapid" or "ligation")
- `description`

The CLI is intentionally thin: it calls library functions
[`list_supported_kits()`], [`kit_is_legacy()`] and [`base_chemistry_of()`] to assemble
the table with minimal duplication.


> **Note:** `porkchop list-kits` prints a **non‑truncated** full-width table by default
> (all rows, full cell contents). Use `--csv` for machine-readable output.


### `describe-kit`
List all **primers**, **adapters**, and **barcodes** used by a specific kit, without truncation.

```bash
porkchop describe-kit LSK114
porkchop describe-kit NBD114.24
porkchop describe-kit RBK114.96 --csv > kit.csv
```

Columns:
- `name`
- `kind` (`adapter_top`, `adapter_bottom`, `primer`, `barcode`)
- `sequence`
- `provenance_url`
- `provenance_ref`

### PCS111 & PCS114 support

> **Note:** PCS114 (Kit 14) is a current product and is **not** considered legacy.

- **PCS111 (SQK‑PCS111)** — classified as *rapid* (uses Rapid adapters). Primers included:
  - `SSP` and `VNP` (legacy PCR‑cDNA), sourced from **pychopper** default primer file (`cDNA_SSP_VNP.fas`).
  - Also includes `CRTA` and `RTP`.
- **PCS114 (SQK‑PCS114, Kit 14)** — uses `SSPII`, `RTP`, `CRTA` and Rapid adapters.

`porkchop describe-kit PCS111` and `porkchop describe-kit PCS114` will list **all primers/adapters** with
full **provenance** columns (`provenance_url`, `provenance_ref`). SSP/VNP rows cite the pychopper file;
Kit 14 primers cite the Chemistry Technical Document Appendix 15.


### PCB111 & PCB114 (PCR‑cDNA Barcoding) — pychopper primers

- **PCB111.24 (Kit 11)** — includes *pychopper*-style primer tokens **SSP** and **VNP**, alongside `CRTA`/`RTP` and Rapid adapters; barcodes **BP01–24**.
- **PCB114.24 (Kit 14)** — uses `SSPII`/`RTP`/`CRTA` (no VNP) with Rapid adapters; barcodes **BP01–24**.

Use:
```bash
porkchop describe-kit PCB111.24
porkchop describe-kit PCB114.24
```
You’ll see **complete, non‑truncated** tables with sequences and provenance (`provenance_url`, `provenance_ref`). SSP/VNP rows cite the pychopper primer FASTA; SSPII/RTP rows cite the Chemistry Technical Document (Appendix 15).


### MAB114 (Microbial Amplicon Barcoding 24 V14)

- **MAB114.24 (SQK‑MAB114.24)** — Rapid‑based workflow for **16S** and **ITS** amplicons with **24 barcodes**.
- `porkchop list-kits` includes MAB114; `porkchop describe-kit MAB114.24` prints Rapid adapter and **BP01–24** barcode sequences.
- Primer oligo sequences for MAB114 are shipped with the kit; where not published in the CHTD appendices, they are omitted here and will be added when formally published by ONT.

## High‑performance sequence I/O (FASTQ/SAM/BAM)

- FASTQ/FASTQ.GZ via `needletail`
- SAM/BAM via `rust-htslib`
- Multithreading: `rayon` threadpool (defaults to all logical cores); BAM enables htslib BGZF threads when possible.
- CLI integrated via `clap`.

### Examples
```bash
# Count reads across multiple inputs using all cores by default
porkchop reads --count foo.fastq.gz bar.bam

# Limit threads (e.g., 8)
porkchop reads --count --threads 8 foo.fastq.gz bar.bam
```

### Library
```rust
use porkchop::seqio::{for_each_parallel, NARead};
let (_, n) = for_each_parallel("reads.fastq.gz", None, |r: NARead| {
    // Your per-read work here
}).unwrap();
println!("processed {n} reads");
```
