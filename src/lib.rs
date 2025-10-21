#![forbid(unsafe_code)]
//! # porkchop
//!
//! Authoritative, well-documented reference of Oxford Nanopore **adapters**, **primers**
//! and **barcodes**, with a strongly-typed mapping from **kit identifiers** to
//! sequences and a fuzzy matcher to infer kit-of-origin from a nucleotide string.
//!
//! ## Highlights
//! - ❗ **No feature flags**: all capabilities are always enabled.
//! - 📚 **Excessive rustdoc**: every item has provenance and examples.
//! - 🧭 **Deterministic data**: sequences are embedded as `&'static str` constants.
//!
//! ## Primary sources
//! - *Chemistry Technical Document* (CHTD_500_v1_revAR_25Nov2024), Appendix 14 (barcodes)
//!   and Appendix 15 (adapter/primer sequences). Public URL in each record's provenance.
//! - Historical adapter fragments from public Porechop forks (we include **only** strings).
//!
//! ## Examples
//! ```rust
//! // Discover kits:
//! for k in porkchop::list_supported_kits() { println!("{} — {}", k.id, k.description); }
//! // Retrieve NB barcodes for a kit:
//! let nb = porkchop::get_sequences_for_kit("NBD114.24").unwrap();
//! assert_eq!(nb.barcodes.len(), 24);
//! // Infer likely kit(s) from a read end containing an LA top adapter:
//! let hits = porkchop::infer_kit_from_sequence("ACGT...TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCT...TGCA", 2, None);
//! assert!(hits.iter().any(|h| h.element == "LA_top"));
//! ```
//!
//! ## Version
//! This build is "0.1.2" ("2025-10-21").

pub mod kit;
pub mod detect;
pub mod data { pub mod adapters; pub mod barcodes; pub mod legacy; }
pub mod kits;

use kit::*;

/// High-level classification for a kit's **base sequencing chemistry**.
///
/// This is intentionally coarse and only differentiates *ligation* vs *rapid* families.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum BaseChemistry { Ligation, Rapid }

/// Return `true` if the kit is considered **legacy** in this crate.
///
/// This is determined by well-known, superseded kit identifiers from historical workflows.
/// The list is small and explicit to avoid false positives.
///
/// # Examples
/// ```
/// let l = porkchop::get_sequences_for_kit("LSK109").map(porkchop::kit_is_legacy).unwrap();
/// assert!(l);
/// ```
pub fn kit_is_legacy(k: &crate::kit::Kit) -> bool {
    matches!(k.id.0, "LSK108" | "LSK109" | "LSK308" | "RBK004" | "RBK110.96")
}

/// Infer the **base chemistry** for a kit by inspecting its adapters/primers.
/// If the Rapid Adapter (RA) is present, the kit is treated as *Rapid*; otherwise *Ligation*.
///
/// # Examples
/// ```
/// let rapid = porkchop::get_sequences_for_kit("RBK114.24").map(porkchop::base_chemistry_of).unwrap();
/// assert_eq!(rapid, porkchop::BaseChemistry::Rapid);
/// let lig = porkchop::get_sequences_for_kit("LSK114").map(porkchop::base_chemistry_of).unwrap();
/// assert_eq!(lig, porkchop::BaseChemistry::Ligation);
/// ```
pub fn base_chemistry_of(k: &crate::kit::Kit) -> BaseChemistry {
    let is_rapid = k.adapters_and_primers.iter().any(|r| r.name.eq_ignore_ascii_case("RA_top") || r.name.eq_ignore_ascii_case("RAP_T"));
    if is_rapid { BaseChemistry::Rapid } else { BaseChemistry::Ligation }
}

/// Convenience: return a vector of rows describing each kit (for CLI/UX).
/// Each row is `(kit_id, description, legacy, base_chemistry)`.
pub fn list_kits_rows() -> Vec<(String, String, bool, BaseChemistry)> {
    list_supported_kits()
        .iter()
        .map(|k| (k.id.0.to_string(), k.description.to_string(), kit_is_legacy(k), base_chemistry_of(k)))
        .collect()
}


/// Return the static registry of supported kits.
pub fn list_supported_kits() -> &'static [crate::kit::Kit] { kits::KITS }

/// Retrieve the full kit definition (adapters/primers + barcodes) for a kit id.
///
/// Kit identifiers are case-insensitive and may be shorthand like `LSK114` or
/// include the pack size like `NBD114.96`.
///
/// # Examples
/// ```rust
/// let k = porkchop::get_sequences_for_kit("RBK114.24").unwrap();
/// assert!(k.description.contains("Rapid Barcoding Kit"));
/// ```
pub fn get_sequences_for_kit(id: &str) -> Option<&'static crate::kit::Kit> {
    kits::KITS.iter().find(|k| k.id.0.eq_ignore_ascii_case(id))
}

/// Fuzzy-match **all known** ONT motifs (adapters, barcodes, flanks) within a query
/// string, returning per-element matches up to `max_edits` edits.
///
/// If `kit_hint` is supplied, the reported [`Match::kit`] is set accordingly.
pub fn infer_kit_from_sequence(seq: &str, max_edits: usize, kit_hint: Option<&str>) -> Vec<Match> {
    let hint = kit_hint.and_then(|s| get_sequences_for_kit(s).map(|k| k.id.clone()));
    let mut records: Vec<SequenceRecord> = Vec::new();
    for k in kits::KITS.iter() {
        records.extend_from_slice(k.adapters_and_primers);
        records.extend_from_slice(k.barcodes);
    }
    detect::find_matches(seq, &records, max_edits, hint)
}

/// Aggregate matches from [`infer_kit_from_sequence`] into **kit likelihoods**.
pub fn infer_kits_from_sequence(seq: &str, max_edits: usize) -> Vec<kit::KitLikelihood> {
    let mut scores: Vec<(KitId, usize, usize)> = kits::KITS.iter().map(|k| (k.id.clone(), 0usize, 0usize)).collect();
    let hits = infer_kit_from_sequence(seq, max_edits, None);
    // Count a hit towards any kit that contains the element.
    for h in hits.iter() {
        for (kid, score, n) in scores.iter_mut() {
            let k = kits::KITS.iter().find(|k| &k.id == kid).unwrap();
            if k.adapters_and_primers.iter().any(|r| r.name == h.element) || k.barcodes.iter().any(|r| r.name == h.element) {
                *score += h.mismatches;
                *n += 1;
            }
        }
    }
    let mut out: Vec<kit::KitLikelihood> = scores.into_iter()
        .filter(|(_,_,n)| *n > 0)
        .map(|(kit, score, n_hits)| kit::KitLikelihood { kit, score, n_hits })
        .collect();
    out.sort_by_key(|kl| (kl.score, core::cmp::Reverse(kl.n_hits), kl.kit.0));
    out
}

/// Crate version string (from `CARGO_PKG_VERSION`).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
