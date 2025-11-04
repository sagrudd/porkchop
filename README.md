# porkchop


Oxford Nanopore™ library/kit inspection, screening and reporting.

**Highlights**
- List & describe kits
- Live TUI screening with strand‑specific hit counts and throughput
- Polars dataframes for kit predictions
- HTML reporting (`--html`)
- Benchmarking of algorithms (edlib, Myers, aCMYers)



## Build
```bash
cargo build --release
```


## Subcommands


### `listkits`
List all supported kits and output in CSV/Markdown/plain table.


**Options**  
- `--format <csv|md|table>` (default: table)  
- `--full` (no truncation for table)  
- `--truncate <N>` (truncate cell width when table & not `--full`)

**Examples**
```bash
porkchop list-kits --format table --full
porkchop list-kits --format csv
porkchop list-kits --format md --truncate 32
```

### `describe`
Describe a specific kit by `--id` with adapters/primers/barcodes.


**Args**
- `--id <KIT_ID>` (e.g., LSK114)

**Example**
```bash
porkchop describe --id LSK114
```

### `benchmark`
Benchmark edit distance algorithms on dataset(s) with an optional truth set.


**Options**
- `--algorithms edlib,myers,acmyers` (default in code)
- `--max-dist <N>` (default: 24)
- `--threads <N>` (0/None = all)

**Args**
- `--files <...>` (FASTQ/FASTA/FASTQ.GZ/SAM/BAM)  
- `--kit <KIT_ID>`  
- `--truth <CSV>` (optional)

**Example**
```bash
porkchop benchmark --files reads.fastq.gz --kit LSK114   --algorithms edlib,myers,acmyers --max-dist 24 --threads 8
```

### `screen`
Screen reads to infer the most likely sequencing kit and display a live TUI.


**Options**
- `--algorithm <edlib|myers|acmyers>` (default: edlib)
- `--max-dist <N>` (default: 24)
- `--fraction <0.0-1.0>` (default: 0.05)
- `--tick <seconds>` (default: 2)
- `--threads <N>` (0/None = all)
- `--json <PATH>` (write contexts as JSON)
- `--kit-prob-min <P>` (default: 0.1) — hide kits with P ≤ threshold
- `--html <PATH>` (write HTML report at end of run)

**Example**
```bash
porkchop screen --files reads.fastq.gz --algorithm edlib --max-dist 24   --fraction 0.05 --tick 2 --kit-prob-min 0.1 --html screen_report.html
```


## Supported Sequencing Kits


- **LSK114** — Ligation Sequencing Kit V14 (LSK114). Uses LA adapter; pairs with Native Barcoding kits.

  
  Line‑art:
  ```
5' ── LA_TOP ── INSERT ── LA_BOTTOM ── 3'
      ↑ Adapter ligation at both ends; barcodes via NBD kits
  ```

- **PCS111** — PCR‑cDNA Sequencing Kit (SQK‑PCS111). Uses legacy SSP/VNP primers and RA; CRTA+RTP included.

  
  Line‑art:
  ```
5' ── RA (transposase‑attached) ── INSERT ── RA ── 3'
      ↑ Rapid attachment; barcodes via RB kits
  ```

- **PCS114** — PCR‑cDNA Sequencing Kit V14 (SQK‑PCS114). Uses SSPII/RTP/CRTA and RA.

  
  Line‑art:
  ```
5' ── VNP/SSPII primer → RT → CRTA/RTP tails ── INSERT ── RA ── 3'
      ↑ cDNA workflow; PCB/NB barcodes optional
  ```

- **LSK114-XL** — Ligation Sequencing Kit V14 XL (LSK114-XL). LA adapter; typical with NBD114.x sets.

  
  Line‑art:
  ```
5' ── LA_TOP ── INSERT ── LA_BOTTOM ── 3'
      ↑ Adapter ligation at both ends; barcodes via NBD kits
  ```

- **NBD114.24** — Native Barcoding Kit 24 V14. Uses NA adapter + NB01–24; NB flanks.

  
  Line‑art:
  ```
5' ── LA_TOP ── INSERT ── LA_BOTTOM ── 3'
      ↑ Adapter ligation at both ends; barcodes via NBD kits
  ```

- **NBD114.96** — Native Barcoding Kit 96 V14. Uses NA adapter + NB01–96; NB flanks.

  
  Line‑art:
  ```
5' ── LA_TOP ── INSERT ── LA_BOTTOM ── 3'
      ↑ Adapter ligation at both ends; barcodes via NBD kits
  ```

- **RBK114.24** — Rapid Barcoding Kit 24 V14. Uses RA adapter + RB01–24 core barcodes; RB flanks.

  
  Line‑art:
  ```
5' ── RA (transposase‑attached) ── INSERT ── RA ── 3'
      ↑ Rapid attachment; barcodes via RB kits
  ```

- **RBK114.96** — Rapid Barcoding Kit 96 V14. Uses RA adapter + RB01–96 core barcodes; RB flanks.

  
  Line‑art:
  ```
5' ── RA (transposase‑attached) ── INSERT ── RA ── 3'
      ↑ Rapid attachment; barcodes via RB kits
  ```

- **PCB111.24** — PCR‑cDNA Barcoding Kit 24 (Kit 11). Uses SSP/VNP (pychopper) with CRTA/RTP and PCB flanks; BP01–24.

  
  Line‑art:
  ```
5' ── RA (transposase‑attached) ── INSERT ── RA ── 3'
      ↑ Rapid attachment; barcodes via RB kits
  ```

- **PCB114.24** — PCR‑cDNA Barcoding Kit 24 V14. Uses cDNA adapters/primers + BP01–24; PCB flanks.

  
  Line‑art:
  ```
5' ── VNP/SSPII primer → RT → CRTA/RTP tails ── INSERT ── RA ── 3'
      ↑ cDNA workflow; PCB/NB barcodes optional
  ```

- **RPB114.24** — Rapid PCR Barcoding Kit 24 V14. Uses RA + RLB01–24 core barcodes; RPB flank.

  
  Line‑art:
  ```
5' ── RA (transposase‑attached) ── INSERT ── RA ── 3'
      ↑ Rapid attachment; barcodes via RB kits
  ```

- **PBC001** — PCR Barcoding Expansion 1–12 (EXP‑PBC001): BC01–BC12; PBC flanks.

  
  Line‑art:
  ```
5' ── RA (transposase‑attached) ── INSERT ── RA ── 3'
      ↑ Rapid attachment; barcodes via RB kits
  ```

- **PBC096** — PCR Barcoding Expansion 1–96 (EXP‑PBC096): BC01–BC96; PBC flanks.

  
  Line‑art:
  ```
5' ── RA (transposase‑attached) ── INSERT ── RA ── 3'
      ↑ Rapid attachment; barcodes via RB kits
  ```

- **RBK004** — Rapid Barcoding Kit (RBK004): RB01–RB12; RB flank.

  
  Line‑art:
  ```
5' ── RA (transposase‑attached) ── INSERT ── RA ── 3'
      ↑ Rapid attachment; barcodes via RB kits
  ```

- **RBK110.96** — Rapid Barcoding Kit 96 (RBK110.96): RB01–RB96; RB flank.

  
  Line‑art:
  ```
5' ── RA (transposase‑attached) ── INSERT ── RA ── 3'
      ↑ Rapid attachment; barcodes via RB kits
  ```

- **LSK109** — Ligation Sequencing Kit LSK109 (legacy). Y‑adapter trunk per Porechop forks.

  
  Line‑art:
  ```
5' ── LA_TOP ── INSERT ── LA_BOTTOM ── 3'
      ↑ Adapter ligation at both ends; barcodes via NBD kits
  ```

- **LSK108** — Ligation Sequencing Kit LSK108 (legacy). Y‑adapter trunk per Porechop forks.

  
  Line‑art:
  ```
5' ── LA_TOP ── INSERT ── LA_BOTTOM ── 3'
      ↑ Adapter ligation at both ends; barcodes via NBD kits
  ```

- **LSK308** — 1D^2 kit LSK308 (legacy). 1D^2 adapter fragments per Porechop forks.

  
  Line‑art:
  ```
5' ── LA_TOP ── INSERT ── LA_BOTTOM ── 3'
      ↑ Adapter ligation at both ends; barcodes via NBD kits
  ```

- **MAB114.24** — Microbial Amplicon Barcoding 24 V14 (SQK‑MAB114.24). Rapid‑based; 16S and ITS targets; 24 barcodes.

  
  Line‑art:
  ```
5' ── PCR primers (with adapters) ── INSERT ── PCR primers ── 3'
      ↑ Amplicon library; barcodes via EXP‑NBD or native
  ```


## HTML Report
When `--html <path>` is provided, a self‑contained report is written at the end of screening (after TUI teardown). It includes:
- Run parameters
- Top synthetic sequences (complete; `(+ )` forward and `(−)` reverse hit counts; deduped reads)
- Top co‑occurrence contexts (complete)
- Sequencing kit predictions (filtered by `--kit-prob-min`)
