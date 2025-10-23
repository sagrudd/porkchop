
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
