//! Benchmark framework for porkchop.
//!
//! This module benchmarks adapter/primer/barcode assignment performance
//! using several alternative classifiers. It is intentionally simple,
//! dependency-light and portable, so it builds on common developer hosts.
//!
//! Algorithms implemented:
//! - `Myers` (edit distance with k threshold, bio crate)
//! - `ACMyers` (Aho–Corasick prefilter + per-candidate Myers)
//! - `Edlib` (C FFI via edlib_rs bindings, distance-only, semiglobal)
//! - `Parasail` (placeholder; returns None to keep build portable)
//!
//! The benchmarking entrypoint is [`benchmark_file`].

use std::time::{Duration, Instant};
use std::sync::{Arc};
use std::sync::atomic::{AtomicU64, Ordering};
use std::path::Path;
use std::collections::HashMap;

use aho_corasick::{AhoCorasick, AhoCorasickBuilder, AhoCorasickKind};
use bio::pattern_matching::myers::{Myers, MyersBuilder};
use edlib_rs::edlibrs::{
    edlibAlign, edlibDefaultAlignConfig, edlibFreeAlignResult,
    EdlibAlignConfig, EdlibAlignMode_EDLIB_MODE_HW, EdlibAlignTask_EDLIB_TASK_DISTANCE,
};

use crate::kit::{SequenceRecord, SeqKind, Kit};
use crate::seqio;

/// Algorithms available to the benchmark.
#[derive(Clone, Copy, Debug)]
pub enum BenchmarkAlgo {
    Myers,
    ACMyers,
    Edlib,
    Parasail,
}

impl std::str::FromStr for BenchmarkAlgo {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_ascii_lowercase().as_str() {
            "myers" => Ok(Self::Myers),
            "acmyers" | "ac-myers" | "aho-myers" => Ok(Self::ACMyers),
            "edlib" => Ok(Self::Edlib),
            "parasail" => Ok(Self::Parasail),
            other => Err(format!("Unknown algorithm: {}", other)),
        }
    }
}

/// A single top hit label from a classifier.
#[derive(Clone, Debug)]
pub struct LabelHit {
    /// Name/identifier of the matched sequence (e.g., "NB01", "RA-top").
    pub name: String,
    /// Kind/category of the sequence (adapter / primer / barcode).
    pub kind: SeqKind,
    /// A score where **larger is better**. We use negative edit distance for
    /// distance-based classifiers (e.g. `-dist`).
    pub score: i32,
    /// Optional end position of the match in the read.
    pub pos: Option<usize>,
}

/// Prebuilt state shared across many classifications (optional).
///
/// We keep this lightweight; only AC is immutable and free to share
/// across threads. We build Myers per-candidate to avoid interior mutability.
pub struct Prebuilt<'a> {
    pub records: &'a [SequenceRecord],
    pub ac: AhoCorasick,
}

/// Build an `AhoCorasick` automaton across all kit motifs.
pub fn prebuild_for<'a>(records: &'a [SequenceRecord]) -> Prebuilt<'a> {
    let patterns = records.iter().map(|r| r.sequence.as_bytes());
    let ac = AhoCorasickBuilder::new()
        .kind(Some(AhoCorasickKind::DFA)) // prefer DFA for short motifs
        .build(patterns)
        .expect("failed to build Aho-Corasick automaton");
    Prebuilt { records, ac }
}

/// Return the best label according to the requested algorithm.
pub fn classify_best(
    algo: BenchmarkAlgo,
    seq: &[u8],
    records: &[SequenceRecord],
    max_dist: usize,
) -> Option<LabelHit> {
    match algo {
        BenchmarkAlgo::Myers => myers_best(seq, records, max_dist),
        BenchmarkAlgo::ACMyers => {
            let pre = prebuild_for(records);
            ac_myers_best(seq, &pre, max_dist)
        }
        BenchmarkAlgo::Edlib => edlib_best(seq, records, max_dist),
        BenchmarkAlgo::Parasail => parasail_best(seq, records),
    }
}

/// Pure Myers (build per-record).
fn myers_best(seq: &[u8], records: &[SequenceRecord], max_dist: usize) -> Option<LabelHit> {
    let mut best: Option<LabelHit> = None;
    for r in records {
        let mut m: Myers<u64> = MyersBuilder::new().build_64(r.sequence.as_bytes().iter().copied());
        if let Some((_, end, dist)) = m.find_all(seq, max_dist as u8).next() {
            let score = -(dist as i32);
            let hit = LabelHit { name: r.name.to_string(), kind: r.kind, score, pos: Some(end) };
            if best.as_ref().map(|b| hit.score > b.score).unwrap_or(true) {
                best = Some(hit);
            }
        }
    }
    best
}

/// Aho–Corasick prefilter then Myers per-candidate.
fn ac_myers_best(seq: &[u8], pre: &Prebuilt, max_dist: usize) -> Option<LabelHit> {
    let mut best: Option<LabelHit> = None;
    // De-duplicate pattern IDs we’ll verify with Myers.
    let mut seen = std::collections::HashSet::new();
    for m in pre.ac.find_iter(seq) {
        let pid = m.pattern();
        if !seen.insert(pid) { continue; }
        let r = &pre.records[pid];
        let mut my: Myers<u64> = MyersBuilder::new().build_64(r.sequence.as_bytes().iter().copied());
        if let Some((_, end, dist)) = my.find_all(seq, max_dist as u8).next() {
            let score = -(dist as i32);
            let hit = LabelHit { name: r.name.to_string(), kind: r.kind, score, pos: Some(end) };
            if best.as_ref().map(|b| hit.score > b.score).unwrap_or(true) {
                best = Some(hit);
            }
        }
    }
    best
}

/// Edlib distance (C FFI; distance-only, semiglobal).
fn edlib_best(seq: &[u8], records: &[SequenceRecord], max_dist: usize) -> Option<LabelHit> {
    let mut best: Option<LabelHit> = None;
    for r in records {
        let mut cfg: EdlibAlignConfig = unsafe { edlibDefaultAlignConfig() };
        cfg.mode = EdlibAlignMode_EDLIB_MODE_HW; // semiglobal (end-free) is close to adapter matching
        cfg.task = EdlibAlignTask_EDLIB_TASK_DISTANCE;
        cfg.k = max_dist as i32;

        let q = r.sequence.as_bytes();
        let res = unsafe {
            edlibAlign(
                q.as_ptr() as *const i8, q.len() as i32,
                seq.as_ptr() as *const i8, seq.len() as i32,
                cfg,
            )
        };
        let dist = res.editDistance;
        // Free internal allocations, if any (safe even when distance-only).
        unsafe { edlibFreeAlignResult(res) };

        if dist >= 0 {
            let score = -(dist as i32);
            let hit = LabelHit { name: r.name.to_string(), kind: r.kind, score, pos: None };
            if best.as_ref().map(|b| hit.score > b.score).unwrap_or(true) {
                best = Some(hit);
            }
        }
    }
    best
}

/// Placeholder to keep the crate portable. Swap in a real parasail-rs
/// implementation when the native library is available on the host.
fn parasail_best(_seq: &[u8], _records: &[SequenceRecord]) -> Option<LabelHit> {
    None
}

/// Load a simple truth map (read_id -> expected_label). Supports CSV or TSV.
pub fn load_truth<P: AsRef<Path>>(path: P) -> anyhow::Result<HashMap<String, String>> {
    let p = path.as_ref();
    let delim = if p.extension().map(|e| e == "tsv").unwrap_or(false) { b'\t' } else { b',' };
    let mut rdr = csv::ReaderBuilder::new().has_headers(true).delimiter(delim).from_path(p)?;
    let mut map = HashMap::new();
    for rec in rdr.records() {
        let r = rec?;
        if r.len() >= 2 {
            map.insert(r[0].to_string(), r[1].to_string());
        }
    }
    Ok(map)
}

/// Benchmark a single file.
///
/// Returns:
/// `(tp, fp, fn, elapsed, nseq, cpu_util (placeholder), input_format)`
pub fn benchmark_file<P: AsRef<Path>>(
    path: P,
    kit: &Kit,
    algo: BenchmarkAlgo,
    truth: Option<HashMap<String, String>>,
    threads: Option<usize>,
) -> anyhow::Result<(u64, u64, u64, Duration, usize, f32, seqio::InputFormat)> {
    let start = Instant::now();

    // Atomic counters to be shared by worker threads.
    let tp = Arc::new(AtomicU64::new(0));
    let fp = Arc::new(AtomicU64::new(0));
    let fn_ = Arc::new(AtomicU64::new(0));
    let nseq = Arc::new(AtomicU64::new(0));

    // Owned truth map moved into closure (if any).
    let truth_owned = truth;

    // Prebuild AC for ACMyers (immutable, thread-safe).
    let pre: Option<Prebuilt> = match algo {
        BenchmarkAlgo::ACMyers => Some(prebuild_for(&kit.adapters_and_primers)),
        _ => None,
    };

    let tp_c = tp.clone();
    let fp_c = fp.clone();
    let fn_c = fn_.clone();
    let nseq_c = nseq.clone();

        // Own a copy of the static records so the closure can capture without borrowing `kit`.
    let records_arc: Arc<Vec<SequenceRecord>> = Arc::new(kit.adapters_and_primers.to_vec());
let fmt_n = seqio::for_each_parallel(path.as_ref(), threads, move |rec: seqio::NARead| {
        nseq_c.fetch_add(1, Ordering::Relaxed);

        let records = records_arc.as_slice();
        let label = match algo {
            BenchmarkAlgo::Myers => myers_best(&rec.seq, records, 24),
            BenchmarkAlgo::ACMyers => {
                // Rebuild a minimal pre each call (safe if `pre` is None),
                // otherwise use the computed AC.
                let local_pre = if let Some(ref pr) = pre { Some(pr) } else { None };
                if let Some(pr) = local_pre { ac_myers_best(&rec.seq, pr, 24) } else { myers_best(&rec.seq, records, 24) }
            }
            BenchmarkAlgo::Edlib => edlib_best(&rec.seq, records, 24),
            BenchmarkAlgo::Parasail => parasail_best(&rec.seq, records),
        };

        if let Some(ref tmap) = truth_owned {
            let id = rec.id.as_str();
            let expected = tmap.get(id);
            match (label.as_ref().map(|l| l.name.as_str()), expected) {
                (Some(found), Some(true_label)) => {
                    if found == *true_label { tp_c.fetch_add(1, Ordering::Relaxed); }
                    else { fp_c.fetch_add(1, Ordering::Relaxed); }
                }
                (Some(_), None) => { fp_c.fetch_add(1, Ordering::Relaxed); },
                (None, Some(_)) => { fn_c.fetch_add(1, Ordering::Relaxed); },
                (None, None) => {},
            }
        }
    })?;

    let elapsed = start.elapsed();

    // Snapshot results
    let tp_v = tp.load(Ordering::Relaxed);
    let fp_v = fp.load(Ordering::Relaxed);
    let fn_v = fn_.load(Ordering::Relaxed);
    let nseq_v = nseq.load(Ordering::Relaxed) as usize;

    // Portable placeholder for CPU util (can wire sysinfo back if desired)
    let cpu_util = 0.0_f32;

    Ok((tp_v, fp_v, fn_v, elapsed, nseq_v, cpu_util, fmt_n.0))
}    // Own a copy of the static records so the closure can capture without borrowing `kit`.