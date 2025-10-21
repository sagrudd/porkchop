//! Sequence detection utilities (exact and fuzzy).
//!
//! The implementation is intentionally dependency-free to keep builds lean and
//! predictable. We expose a fast Levenshtein edit-distance over short motifs
//! (adapters/barcodes/flanks are <= ~80 nt), and a sliding-window scan.
//!
//! # Examples
//! ```
//! use porkchop::detect::{edit_distance, best_window_edit};
//! assert_eq!(edit_distance("ACGT", "ACGT"), 0);
//! assert!(best_window_edit("NNNNACGTNN", "ACGT").unwrap().2 <= 0);
//! ```
use crate::kit::{Match, SequenceRecord, KitId};

/// Compute Levenshtein edit distance between two ASCII strings (DNA alphabet).
#[inline]
pub fn edit_distance(a: &str, b: &str) -> usize {
    let (a, b) = (a.as_bytes(), b.as_bytes());
    let mut prev: Vec<usize> = (0..=b.len()).collect();
    let mut curr = vec![0usize; b.len() + 1];
    for (i, &ca) in a.iter().enumerate() {
        curr[0] = i + 1;
        for (j, &cb) in b.iter().enumerate() {
            let cost = if ca == cb { 0 } else { 1 };
            let ins = curr[j] + 1;
            let del = prev[j + 1] + 1;
            let sub = prev[j] + cost;
            curr[j + 1] = ins.min(del).min(sub);
        }
        prev.clone_from(&curr);
    }
    prev[b.len()]
}

/// Slide `needle` across `haystack`, returning best (lowest) edit distance and span.
pub fn best_window_edit(haystack: &str, needle: &str) -> Option<(usize, usize, usize)> {
    if needle.is_empty() || haystack.len() < needle.len() { return None; }
    let h = haystack.as_bytes();
    let nlen = needle.len();
    let mut best: Option<(usize, usize, usize)> = None;
    for i in 0..=h.len() - nlen {
        let window = &haystack[i..i + nlen];
        let d = edit_distance(window, needle);
        if best.map_or(true, |(_, _, bd)| d < bd) {
            best = Some((i, i + nlen, d));
        }
    }
    best
}

/// Find matches to any of the provided records in `query`, allowing up to `max_edits`.
pub fn find_matches<'a>(query: &str, records: &'a [SequenceRecord], max_edits: usize, kit_hint: Option<KitId>) -> Vec<Match> {
    let q = query.to_ascii_uppercase();
    let mut hits = Vec::new();
    for r in records {
        if let Some((s, e, d)) = best_window_edit(&q, r.sequence) {
            if d <= max_edits {
                hits.push(Match { kit: kit_hint.clone(), element: r.name, kind: r.kind, start: s, end: e, mismatches: d });
            }
        }
    }
    hits
}
