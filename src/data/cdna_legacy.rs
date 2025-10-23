
//! Legacy cDNA primer sequences used by historical PCR‑cDNA workflows (PCS109/PCS111).
//!
//! These provide the canonical tokens **SSP** (strand‑switching primer) and **VNP**
//! (VN primer / oligo‑dT) as used in the community tool **pychopper**.
//!
//! # Provenance
//! - Repository: https://github.com/epi2me-labs/pychopper
//! - File: `pychopper/primer_data/cDNA_SSP_VNP.fas`
//! - Purpose: default primers distributed with pychopper (token names `SSP` and `VNP`).
//!
//! For Kit 14 (PCS114), use the **SSPII** / **RTP** primers from the Chemistry
//! Technical Document (see `crate::data::adapters`).

use crate::kit::{SequenceRecord, SeqKind, Provenance};

const PYCHOPPER_DEFAULTS: Provenance = Provenance {
    source: "https://github.com/epi2me-labs/pychopper/blob/master/pychopper/primer_data/cDNA_SSP_VNP.fas",
    appendix: Some("primer_data/cDNA_SSP_VNP.fas"),
    notes: Some("Community primer file providing SSP/VNP tokens."),
};

/// SSP — classic strand‑switching primer (historical PCR‑cDNA).
pub const SSP: SequenceRecord = SequenceRecord {
    name: "SSP",
    kind: SeqKind::Primer,
    sequence: "TTTCTGTTGGTGCTGATATTGCTGGG",
    provenance: PYCHOPPER_DEFAULTS,
};

/// VNP — classic oligo‑dT VN primer (historical PCR‑cDNA).
pub const VNP: SequenceRecord = SequenceRecord {
    name: "VNP",
    kind: SeqKind::Primer,
    sequence: "ACTTGCCTGTCGCTCTATCTTCTTTTTTTTT",
    provenance: PYCHOPPER_DEFAULTS,
};

/// Convenience slice for legacy PCR‑cDNA primers.
pub const LEGACY_CDNA_PRIMERS: &[SequenceRecord] = &[SSP, VNP];