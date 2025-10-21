//! Core types for **kits**, **sequence records** and **provenance**.
//!
//! This module holds the data model used across the crate. It is intentionally
//! simple and `no_std`-friendly (we stay on owned `&'static str` to allow all
//! sequences to live in the binary as constants).
//!
//! # Provenance
//! Every `SequenceRecord` contains a [`Provenance`] entry that records the
//! original source for the sequence (e.g. the *Chemistry Technical Document* or
//! a historical adapter list in a Porechop fork).
//!
//! See the crate-level docs and README for references.
use core::fmt;

/// Canonical identifier for an ONT kit (e.g. `"LSK114"`, `"NBD114.24"`, `"RBK114.96"`).
#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub struct KitId(pub &'static str);

impl fmt::Display for KitId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result { write!(f, "{}", self.0) }
}

/// Type of motif recorded in this crate.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum SeqKind {
    /// Top strand of a double-stranded adapter.
    AdapterTop,
    /// Bottom strand of a double-stranded adapter.
    AdapterBottom,
    /// Primer sequence (e.g. RTP/SSPII/cPRM).
    Primer,
    /// Barcode core sequence (as published by ONT).
    Barcode,
    /// Flanking motif used around barcodes/primers.
    Flank,
}

/// Where a sequence string came from.
#[derive(Clone, Debug)]
pub struct Provenance {
    /// Human-readable source (e.g. document name).
    pub source: &'static str,
    /// Public URL for the source.
    pub url: &'static str,
    /// Section, appendix or file reference within the source.
    pub reference: &'static str,
    /// Any helpful notes (ambiguity, wobble bases, legacy status, etc.).
    pub notes: &'static str,
}

/// A single sequence record with name, type, nucleotide string and provenance.
#[derive(Clone, Debug)]
pub struct SequenceRecord {
    /// Short stable name (e.g. `"LA_top"`, `"NB01"`, `"RB_flank"`).
    pub name: &'static str,
    /// Sequence category.
    pub kind: SeqKind,
    /// Uppercase IUPAC string (DNA alphabet, may include wobble codes, see docs).
    pub sequence: &'static str,
    /// Source information for auditability.
    pub provenance: Provenance,
}

/// A kit definition: ties a kit ID to its adapters/primers and (optionally) barcodes.
#[derive(Clone, Debug)]
pub struct Kit {
    /// Identifier such as `"LSK114"` or `"NBD114.96"`.
    pub id: KitId,
    /// One-line description (chemistry family and scope).
    pub description: &'static str,
    /// Sequences integral to the kit chemistry (adapters, primers, flanks).
    pub adapters_and_primers: &'static [SequenceRecord],
    /// Barcode cores available for the kit (empty for non-barcoding kits).
    pub barcodes: &'static [SequenceRecord],
}

/// A match against any known sequence in this crate.
#[derive(Clone, Debug)]
pub struct Match {
    /// Known kit if the search was constrained to a kit; otherwise `None`.
    pub kit: Option<KitId>,
    /// The named element that matched (e.g. `"NB07"`, `"LA_top"`).
    pub element: &'static str,
    /// Category of the element.
    pub kind: SeqKind,
    /// Start index (0-based) of the best window within the query.
    pub start: usize,
    /// End index (exclusive) of the best window within the query.
    pub end: usize,
    /// Number of mismatches (Levenshtein edits) across the window.
    pub mismatches: usize,
}

/// A likelihood summary used by [`crate::infer_kits_from_sequence`].
#[derive(Clone, Debug)]
pub struct KitLikelihood {
    /// Kit identifier.
    pub kit: KitId,
    /// Aggregate score (lower is better) derived from edit distances for matched elements.
    pub score: usize,
    /// How many elements from this kit contributed to the score.
    pub n_hits: usize,
}
