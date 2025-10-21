# Changelog


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
