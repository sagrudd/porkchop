# porkchop

**Authoritative reference** of Oxford Nanopore adapters, primers and barcodes, mapped to kit identifiers, with a fuzzy matcher to infer kit-of-origin.

- License: MPL-2.0
- Version: 0.1.2
- Stability: public API is additive within 0.1.x
- No optional features: all functionality is enabled by default.

## Install
```toml
porkchop = "0.1"
```

## Kit coverage
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

This project is **not affiliated** with Oxford Nanopore Technologies.

## MSRV
Rust 1.70.

## License
MPL-2.0
