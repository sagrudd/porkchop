//! Adapter & primer sequences for **current kits** (CHTD Appendix 15).
//!
//! Source: Chemistry Technical Document (CHTD_500_v1_revAR_25Nov2024) →
//! Appendix 15: *Adapter sequences*.
//!
//! - LA/NA/RA top/bottom sequences
//! - Rapid Adapter T (RAP T)
//! - RNA/cDNA: RTP, SSPII (with wobble codes), CRTA, cPRM forward/reverse
//!
//! Notes:
//! - Sequences are uppercase as published (wobble codes retained: V = A/C/G; mG = riboguanosine).
//! - These are sufficient for adapter detection and kit mapping.

use crate::kit::{SequenceRecord, SeqKind, Provenance};

const CHTD_A15: Provenance = Provenance {
    source: "https://nanoporetech.com/document/chemistry-technical-document",
    appendix: Some("Appendix 15: Adapter sequences"),
    notes: Some("Sequences transcribed verbatim from ONT documentation."),
};

/// Ligation Adapter (LA) top strand. 5'-TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCT-3'
/// const `LA_TOP` — auto‑generated rustdoc.
pub const LA_TOP: SequenceRecord = SequenceRecord {
    name: "LA_top",
    kind: SeqKind::AdapterTop,
    sequence: "TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCT",
    provenance: CHTD_A15,
};

/// Ligation Adapter (LA) bottom strand. 5’-GCAATACGTAACTGAACGAAGTACAGG-3’
/// const `LA_BOTTOM` — auto‑generated rustdoc.
pub const LA_BOTTOM: SequenceRecord = SequenceRecord {
    name: "LA_bottom",
    kind: SeqKind::AdapterBottom,
    sequence: "GCAATACGTAACTGAACGAAGTACAGG",
    provenance: CHTD_A15,
};

/// Native Adapter (NA) top strand (same top as LA).
/// const `NA_TOP` — auto‑generated rustdoc.
pub const NA_TOP: SequenceRecord = SequenceRecord {
    name: "NA_top",
    kind: SeqKind::AdapterTop,
    sequence: "TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCT",
    provenance: CHTD_A15,
};

/// Native Adapter (NA) bottom strand. 5'-ACGTAACTGAACGAAGTACAGG-3'
/// const `NA_BOTTOM` — auto‑generated rustdoc.
pub const NA_BOTTOM: SequenceRecord = SequenceRecord {
    name: "NA_bottom",
    kind: SeqKind::AdapterBottom,
    sequence: "ACGTAACTGAACGAAGTACAGG",
    provenance: CHTD_A15,
};

/// Rapid Adapter (RA) top strand.
/// const `RA_TOP` — auto‑generated rustdoc.
pub const RA_TOP: SequenceRecord = SequenceRecord {
    name: "RA_top",
    kind: SeqKind::AdapterTop,
    sequence: "TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCT",
    provenance: CHTD_A15,
};

/// Rapid Adapter T (RAP T).
/// const `RAP_T` — auto‑generated rustdoc.
pub const RAP_T: SequenceRecord = SequenceRecord {
    name: "RAP_T",
    kind: SeqKind::AdapterTop,
    sequence: "TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCT",
    provenance: CHTD_A15,
};

/// RNA/cDNA RT Primer (RTP).
/// const `RTP` — auto‑generated rustdoc.
pub const RTP: SequenceRecord = SequenceRecord {
    name: "RTP",
    kind: SeqKind::Primer,
    sequence: "CTTGCCTGTCGCTCTATCTTCAGAGGAG",
    provenance: CHTD_A15,
};

/// Strand Switching Primer II (SSPII). V=A/C/G; mG=riboguanosine (ONT notation).
/// const `SSPII` — auto‑generated rustdoc.
pub const SSPII: SequenceRecord = SequenceRecord {
    name: "SSPII",
    kind: SeqKind::Primer,
    sequence: "TTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTmGmGmG",
    provenance: Provenance { notes: Some("V = A/C/G; mG = riboguanosine"), ..CHTD_A15 },
};

/// cDNA RT Adapter (CRTA).
/// const `CRTA` — auto‑generated rustdoc.
pub const CRTA: SequenceRecord = SequenceRecord {
    name: "CRTA",
    kind: SeqKind::Primer,
    sequence: "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAG",
    provenance: CHTD_A15,
};

/// cDNA Primer (cPRM) forward sequence.
/// const `CPRM_FWD` — auto‑generated rustdoc.
pub const CPRM_FWD: SequenceRecord = SequenceRecord {
    name: "cPRM_forward",
    kind: SeqKind::Primer,
    sequence: "ATCGCCTACCGTGACAAGAAAGTTGTCGGTGTCTTTGTGACTTGCCTGTCGCTCTATCTTC",
    provenance: CHTD_A15,
};

/// cDNA Primer (cPRM) reverse sequence.
/// const `CPRM_REV` — auto‑generated rustdoc.
pub const CPRM_REV: SequenceRecord = SequenceRecord {
    name: "cPRM_reverse",
    kind: SeqKind::Primer,
    sequence: "ATCGCCTACCGTGACAAGAAAGTTGTCGGTGTCTTTGTGTTTCTGTTGGTGCTGATATTGC",
    provenance: CHTD_A15,
};

/// Convenience: all adapter/primer records for current kits.
/// const `CURRENT_ADAPTERS_AND_PRIMERS` — auto‑generated rustdoc.
pub const CURRENT_ADAPTERS_AND_PRIMERS: &[SequenceRecord] = &[
    LA_TOP, LA_BOTTOM, NA_TOP, NA_BOTTOM, RA_TOP, RAP_T, RTP, SSPII, CRTA, CPRM_FWD, CPRM_REV,
];