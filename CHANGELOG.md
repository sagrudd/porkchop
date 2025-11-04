
## 0.2.87 - 2025-11-04
### Fixed
- Removed stray duplicated block after `tui_loop` that re-opened the draw loop (`} if done.load(...)`), fixing brace mismatch.

## 0.2.86 - 2025-11-04
### Screen TUI
- Rewrote `tui_loop` with a clean, bracket‑balanced implementation.
- Proper `terminal.draw(|f| { ... })` closure, explicit layout (header + 2 columns), and safe event loop.
- Header shows: screened, total hits, hits/read, reads with ≥1 hit, unclassified %, skipped %.
- Top sequences table now: name, kind, (+), (−), reads (per‑read deduped).
- Top co‑occurrence table unchanged, sorted by count.

## 0.2.85 - 2025-11-04
### Fixed
- Repaired TUI draw loop braces and restored valid stats line (no ellipses).
- Removed accidental Polars pretty‑print block from inside the draw loop; final kit table still printed **after** the TUI exits.

## 0.2.84 - 2025-11-04
### Screen improvements
- Counts in **Top synthetic sequences** are now per **read** (deduped by `(name, kind)`), avoiding inflated counts when multiple hits occur in one read.
- Added columns **(+)**/**(-)** with the motif sequence and its **reverse complement**.
- Header now shows **reads with ≥1 hit** and **hits/read** for clarity.
- After the TUI closes, the **kit-likelihood analysis** is printed as a wide **Polars DataFrame** at the console (as well as JSON on disk).

## 0.2.83 - 2025-11-04
### Fixed
- Attached `#[derive(Subcommand)]` to `enum Commands` to restore clap parsing.
- Corrected Markdown pipe escaping in `print_df_markdown` to a valid Rust string.

## 0.2.82 - 2025-11-04
### Added
- `list-kits` now supports `--format <csv|md|table>` (default: `table`).
  - `--full` and `--truncate <N>` apply **only** when `--format table` is selected.

## 0.2.81 - 2025-11-04
### Added
- `list-kits` flags:
  - `--csv`: write CSV to stdout
  - `--md`: write GitHub‑flavoured Markdown table to stdout
  - `--full`: show all columns/rows with no truncation (wide table)
  - `--truncate <N>`: truncate string cells to N characters (ignored if `--full`)

## 0.2.80 - 2025-11-04
### Fixed
- Prevented Polars pretty-printer overflow by capping `POLARS_TABLE_WIDTH=65535` and using `POLARS_FMT_MAX_ROWS=1000000`.

## 0.2.79 - 2025-11-04
### Changed
- `list-kits` now prints a **Polars DataFrame** (not CSV) with **all columns** and **full cell width**.
  Formatting is controlled via env vars set in the command (`POLARS_FMT_MAX_COLS`, `POLARS_FMT_STR_LEN`, etc.).

## 0.2.78 - 2025-11-04
### Fixed
- Corrected `chemistry` values in `list-kits` output (Polars table):
  - `LSK114`, `LSK114-XL`, `LSK109`, `LSK108`, `LSK308`, `NBD114.24`, `NBD114.96` → **ligation sequencing**
  - `PCS114`, `PCB114.24` → **pcr-cdna sequencing**
  - `MAB114.24` → **amplicon sequencing**
  - Confirmed rapid kits remain **rapid sequencing** (`RBK114.24`, `RBK114.96`, `RPB114.24`).

## [0.2.38] - 2025-11-03
- CLI: “chemistry” column now uses canonical labels (`BaseChemistry::label()` / `Display`).
- Core: add `BaseChemistry::label()` + `Display` (stable human-readable names).
- Example: LSK114 → "ligation sequencing".


## [0.2.37] - 2025-10-23
- Remove unnecessary `mut` on DataFrame bindings where CSV writing isn't used (silences `unused_mut` warning).


## [0.2.36] - 2025-10-23
- CLI: avoid moving `truth_map` by using `truth_map.clone()` per iteration.
- CLI: prefix `input_format` with underscore to silence unused warning.
- CLI: make the DataFrame binding mutable for `CsvWriter::finish(&mut df)`.


## [0.2.35] - 2025-10-23
- Add `BenchmarkAlgo::as_str()` and `Display` impl to support `algo.as_str()` usage in CLI output.


## [0.2.34] - 2025-10-23
- Remove line 54 from `src/benchmark.rs` (previous content: `}`).


## [0.2.33] - 2025-10-23
- Restore `enum BenchmarkAlgo` block in `src/benchmark.rs` (previously corrupted line removed).
- Keep `impl BenchmarkAlgo::from_list(&str)` and `benchmark_file(..., max_dist)` adjustments.


## [0.2.32] - 2025-10-23
- Fix CLI benchmark invocation: add `BenchmarkAlgo::from_list`, pass `&Kit` not `&&str`.
- Thread `max_dist` through `benchmark::benchmark_file` and remove hard-coded threshold.
- Minor: drop unnecessary `mut` on CSV writer.

## [0.2.31] - 2025-10-23
- Hotfix: repair corrupted `version` line in `Cargo.toml`; ensure README sync.

## [0.2.29] - 2025-10-23
- Add CI for MSRV + stable
- Add rust-toolchain.toml (stable + rustfmt + clippy)
- Add Makefile helper targets
- Sync README **Version**

# Changelog

## 0.2.1 - 2025-10-22
- Patch release: preserves 0.2.0 benchmark framework and CLI.
- Merge-ready hooks for upstream kit registry (this archive includes placeholder kits so it compiles).

## 0.2.0 - 2025-10-22
- Added benchmark command and framework (Aho-Corasick, Myers, Edlib, Parasail, AC two-stage).
- High-performance IO for FASTQ/FASTQ.GZ/SAM/BAM.


## 0.2.2 - 2025-10-22
### Improved
- **Benchmark performance:** prebuilt Myers automata and Parasail profiles (per‑run) to cut per‑read setup costs.
- **Docs:** expanded crate/module/CLI documentation; added MSRV section and tuning tips.
- **MSRV:** declare `rust-version = "1.90"` in `Cargo.toml`.

### Maintained
- Preserved v0.2.1 APIs and CLI behavior.


## 0.2.3 - 2025-10-22
### Fixed
- Set `edlib_rs = "0.1.2"` and adapted `benchmark` to use `edlibrs::edlibAlignRs` API (HW mode, distance task).


## 0.2.4 - 2025-10-22
### Fixed
- Pin `parasailors = "0.3.1"` to match crates.io availability; no source changes required (API compatible with our usage).


## 0.2.5 - 2025-10-22
### Fixed
- Removed invalid `hts-sys` feature from `rust-htslib` (using only valid features: `bzip2`, `lzma`, `libdeflate`).

## 0.2.8 - 2025-10-22
### Changed
- Fully refactored `src/benchmark.rs`: removed residual code, simplified algorithm dispatch, added exhaustive rustdoc.
- AC+Myers, Myers, and Smith–Waterman paths are now cleanly implemented and portable.
- edlib/parasail entrypoints remain but route to Smith–Waterman for portability.

## 0.2.10 - 2025-10-22
### Fixed
- `kit.rs`: corrected `Provenance` type to `&'static str` for `source`.
- `benchmark.rs`: now uses `(end, dist)` from `find_all_end`; removed incorrect `distance(end)` calls; unwrapped AC builder result; consistent `pos: Some(end)`.
- `seqio.rs`: removed unused BAM record allocation.

## 0.2.88
### Fixed
- Removed remaining orphan `} if done.load(..)` fragment following `tui_loop`.
