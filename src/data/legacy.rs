//! Legacy adapter fragments derived from public Porechop forks.
//!
//! We **only** include sequence strings (facts) and *no GPL code*. Sources:
//! - Sn0flingan/Poresnip `porechop/adapters.py`
//!   (Y-adapter trunk and 1D^2 sequences)
use crate::kit::{SequenceRecord, SeqKind, Provenance};

const PORECHOP_FORK: Provenance = Provenance {
    source: "https://github.com/Sn0flingan/Poresnip/blob/master/porechop/adapters.py",
    appendix: Some("adapters.py (kit_adapters)"),
    notes: Some("Motor-binding regions are omitted in the fork; we reference ‘trunk’ fragments."),
};

/// Legacy Y‑adapter trunk (SQK‑NSK007/LSK108/LSK109) — start/top.
pub const NSK007_Y_TOP_TRUNK: SequenceRecord = SequenceRecord {
    name: "SQK-NSK007_Y_Top_trunk",
    kind: SeqKind::AdapterTop,
    sequence: "AATGTACTTCGTTCAGTTACGTATTGCT",
    provenance: PORECHOP_FORK,
};

/// Legacy Y‑adapter bottom.
pub const NSK007_Y_BOTTOM: SequenceRecord = SequenceRecord {
    name: "SQK-NSK007_Y_Bottom",
    kind: SeqKind::AdapterBottom,
    sequence: "GCAATACGTAACTGAACGAAGT",
    provenance: PORECHOP_FORK,
};

/// 1D^2 top fragment (LSK308).
pub const LSK308_1D2_TOP: SequenceRecord = SequenceRecord {
    name: "SQK-LSK308_1D2_Top",
    kind: SeqKind::AdapterTop,
    sequence: "GTCAGAGAGGTTCCAAGTCAGAGAGGTTCCT",
    provenance: PORECHOP_FORK,
};

/// 1D^2 bottom fragment (LSK308).
pub const LSK308_1D2_BOTTOM: SequenceRecord = SequenceRecord {
    name: "SQK-LSK308_1D2_Bottom",
    kind: SeqKind::AdapterBottom,
    sequence: "GGCGTCTGCTTGGGTGTTTAACCTTTTTGTCAGAGAGGTTCCAAGTCAGAGAGGTTCCT",
    provenance: PORECHOP_FORK,
};

/// Group of legacy fragments.
pub const LEGACY: &[SequenceRecord] = &[
    NSK007_Y_TOP_TRUNK, NSK007_Y_BOTTOM, LSK308_1D2_TOP, LSK308_1D2_BOTTOM,
];