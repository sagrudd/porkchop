# Changelog


## 0.1.20 - 2025-10-21
### Fixed
- Eliminated `threads` unused-variable warning by explicitly deriving pool size and passing to `with_threadpool` in FASTQ path.


## 0.1.19 - 2025-10-21
### Fixed
- Respect `--threads` in `seqio` by running scoped tasks inside a sized rayon pool (no unused warning).
- CLI: import `rayon::prelude::*` and use the record count returned from `for_each_parallel` per file.
- Ensured archive version and `Cargo.toml` match: 0.1.19.


## 0.1.15 - 2025-10-21
### Fixed
- Corrected `KITS` registry formatting for **MAB114.24** (malformed struct caused parse errors).


## 0.1.14 - 2025-10-21
### Added
- **MAB114.24** kit (Rapid-based Microbial Amplicon Barcoding for 16S & ITS) using shared barcodes BP01–24; Rapid adapter listed.
- CLI and docs updated: `list-kits` shows MAB114; `describe-kit MAB114.24` prints barcodes and adapter.
### Note
- Primer sequences for MAB114 (16S/ITS) will be added when available in the Chemistry Technical Document.


## 0.1.13 - 2025-10-21
### Added
- **PCB111.24** support (pychopper SSP/VNP + CRTA/RTP; PCB flanks; BP01–24).
- Docs refreshed for PCB111/PCB114 and CLI examples.


## 0.1.11 - 2025-10-21
### Fixed
- Restored `pub mod data { ... }` export in `lib.rs` so `crate::data::...` paths resolve.
- Removed accidental inline primer definitions from `lib.rs` (now properly in `data::cdna_legacy`).
- Rewrote `data/cdna_legacy.rs` to a complete, compile-ready module (no placeholders), with full `Provenance` fields.


## 0.1.10 - 2025-10-21
### Fixed
- Resolved duplicate `pub mod data` declaration causing E0428.
- Completed `Provenance` fields for SSP/VNP (added `source`) to fix E0063.


## 0.1.9 - 2025-10-21
### Fixed
- Clarified that **PCS114** is **not legacy**; added a unit test to assert the legacy flag is false.

## 0.1.8 - 2025-10-21
### Added
- New kits: **PCS111** and **PCS114** (PCR‑cDNA Sequencing, Kit 11 & Kit 14).
- Legacy cDNA primers: **SSP** and **VNP** (from epi2me‑labs/pychopper), exposed with full provenance.
- `describe-kit` now reports SSP/VNP for PCS111; PCS114 lists SSPII/RTP/CRTA per CHTD.

### Notes
- Base chemistry classification reports these kits as *rapid* (Rapid adapters used post‑PCR).

## 0.1.7 - 2025-10-21
- README affiliation fix.
