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
impl BenchmarkAlgo {
    /// Return a stable lowercase name for display/CSV.
    pub fn as_str(&self) -> &'static str {
        match self {
            BenchmarkAlgo::Myers => "myers",
            BenchmarkAlgo::ACMyers => "acmyers",
            BenchmarkAlgo::Edlib => "edlib",
            BenchmarkAlgo::Parasail => "parasail",
        }
    }

    /// Parse a comma-separated list of algorithms into a Vec<BenchmarkAlgo>.
    /// Unknown names are ignored with a warning.
    pub fn from_list(s: &str) -> Vec<BenchmarkAlgo> {
        let mut v = Vec::new();
        for t in s.split(',').map(|t| t.trim()).filter(|t| !t.is_empty()) {
            match t.parse::<BenchmarkAlgo>() {
                Ok(a) => v.push(a),
                Err(_) => eprintln!("Warning: unknown algorithm '{}'; skipping", t),
            }
        }
        v
    }
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
fn revcomp_bytes(seq: &[u8]) -> Vec<u8> {
    fn comp(b: u8) -> u8 {
        match b {
            b'A'|b'a' => b'T',
            b'C'|b'c' => b'G',
            b'G'|b'g' => b'C',
            b'T'|b't' => b'A',
            b'U'|b'u' => b'A',
            b'N'|b'n' => b'N',
            _ => b'N',
        }
    }
    let mut out = Vec::with_capacity(seq.len());
    for &b in seq.iter().rev() { out.push(comp(b)); }
    out
}

pub struct Prebuilt {
    pub records: Arc<Vec<SequenceRecord>>,
    pub ac: AhoCorasick,
    pub pat2rec: Vec<usize>,
    pub pat_is_rc: Vec<bool>,
}

/// Build an `AhoCorasick` automaton across all kit motifs.
pub fn prebuild_for(records: &[SequenceRecord]) -> Prebuilt {
    let mut pats: Vec<Vec<u8>> = Vec::new();
    let mut pat2rec: Vec<usize> = Vec::new();
    let mut pat_is_rc: Vec<bool> = Vec::new();
    let owned: Arc<Vec<SequenceRecord>> = Arc::new(records.to_vec());
    for (i, r) in owned.iter().enumerate() {
        let fwd = r.sequence.as_bytes().to_vec();
        pats.push(fwd); pat2rec.push(i); pat_is_rc.push(false);
        let rc = revcomp_bytes(r.sequence.as_bytes());
        pats.push(rc); pat2rec.push(i); pat_is_rc.push(true);
    }
    let pat_refs: Vec<&[u8]> = pats.iter().map(|v| v.as_slice()).collect();
    let ac = AhoCorasickBuilder::new()
        .kind(Some(AhoCorasickKind::DFA))
        .build(pat_refs)
        .expect("failed to build Aho-Corasick automaton");
    Prebuilt { records: owned, ac, pat2rec, pat_is_rc }
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
        // forward
        let mut m: Myers<u64> = MyersBuilder::new().build_64(r.sequence.as_bytes().iter().copied());
        if let Some((_, end, dist)) = m.find_all(seq, max_dist as u8).next() {
            let score = -(dist as i32);
            let hit = LabelHit { name: r.name.to_string(), kind: r.kind, score, pos: Some(end) };
            if best.as_ref().map(|b| hit.score > b.score).unwrap_or(true) { best = Some(hit); }
        } else {
            // reverse-complement of reference
            let rc = revcomp_bytes(r.sequence.as_bytes());
            let mut mrc: Myers<u64> = MyersBuilder::new().build_64(rc.into_iter());
            if let Some((_, end, dist)) = mrc.find_all(seq, max_dist as u8).next() {
                let score = -(dist as i32);
                let hit = LabelHit { name: r.name.to_string(), kind: r.kind, score, pos: Some(end) };
                if best.as_ref().map(|b| hit.score > b.score).unwrap_or(true) { best = Some(hit); }
            }
        }
    }
    best
}

/// Aho–Corasick prefilter then Myers per-candidate.
fn ac_myers_best(seq: &[u8], pre: &Prebuilt, max_dist: usize) -> Option<LabelHit> {
    let mut best: Option<LabelHit> = None;
    let mut seen = std::collections::HashSet::new();
    for m in pre.ac.find_iter(seq) {
        let pid = m.pattern();
        if !seen.insert(pid) { continue; }
        let ridx = pre.pat2rec[pid];
        let r = &&pre.records[ridx];
        let pat_bytes: Vec<u8> = if pre.pat_is_rc[pid] {
            revcomp_bytes(r.sequence.as_bytes())
        } else { r.sequence.as_bytes().to_vec() };
        let mut my: Myers<u64> = MyersBuilder::new().build_64(pat_bytes.into_iter());
        if let Some((_, end, dist)) = my.find_all(seq, max_dist as u8).next() {
            let score = -(dist as i32);
            let hit = LabelHit { name: r.name.to_string(), kind: r.kind, score, pos: Some(end) };
            if best.as_ref().map(|b| hit.score > b.score).unwrap_or(true) { best = Some(hit); }
        }
    }
    best
}


/// Edlib distance (C FFI; distance-only, semiglobal).
fn edlib_best(seq: &[u8], records: &[SequenceRecord], max_dist: usize) -> Option<LabelHit> {
    let mut best: Option<LabelHit> = None;
    for r in records {
        let mut cfg: EdlibAlignConfig = unsafe { edlibDefaultAlignConfig() };
        cfg.mode = EdlibAlignMode_EDLIB_MODE_HW; // semiglobal (end-free)
        cfg.task = EdlibAlignTask_EDLIB_TASK_DISTANCE;
        cfg.k = max_dist as i32;

        // forward
        let q = r.sequence.as_bytes();
        let res = unsafe { edlibAlign(q.as_ptr() as *const i8, q.len() as i32, seq.as_ptr() as *const i8, seq.len() as i32, cfg) };
        let mut best_local: Option<(i32, Option<i32>)> = None;
        if res.editDistance >= 0 { best_local = Some((res.editDistance, None)); }
        unsafe { edlibFreeAlignResult(res) };

        // reverse-complement of reference
        let rc = revcomp_bytes(r.sequence.as_bytes());
        let res2 = unsafe { edlibAlign(rc.as_ptr() as *const i8, rc.len() as i32, seq.as_ptr() as *const i8, seq.len() as i32, cfg) };
        if res2.editDistance >= 0 {
            if let Some((d,_)) = best_local {
                if res2.editDistance < d { best_local = Some((res2.editDistance, None)); }
            } else {
                best_local = Some((res2.editDistance, None));
            }
        }
        unsafe { edlibFreeAlignResult(res2) };

        if let Some((d,_)) = best_local {
            let score = -(d as i32);
            let hit = LabelHit { name: r.name.to_string(), kind: r.kind, score, pos: None };
            if best.as_ref().map(|b| hit.score > b.score).unwrap_or(true) { best = Some(hit); }
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
    max_dist: usize,
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
            BenchmarkAlgo::Myers => myers_best(&rec.seq, records, max_dist),
            BenchmarkAlgo::ACMyers => {
                // Rebuild a minimal pre each call (safe if `pre` is None),
                // otherwise use the computed AC.
                let local_pre = if let Some(ref pr) = pre { Some(pr) } else { None };
                if let Some(pr) = local_pre { ac_myers_best(&rec.seq, pr, max_dist) } else { myers_best(&rec.seq, records, max_dist) }
            }
            BenchmarkAlgo::Edlib => edlib_best(&rec.seq, records, max_dist),
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

impl std::fmt::Display for BenchmarkAlgo {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}


pub fn classify_all(
    algo: BenchmarkAlgo,
    seq: &[u8],
    records: &[SequenceRecord],
    prebuilt: Option<&Prebuilt>,
    max_dist: usize,
) -> Vec<(String, SeqKind, bool)> {
    let mut out: Vec<(String, SeqKind, bool)> = Vec::new();
    match algo {
        BenchmarkAlgo::ACMyers => {
            if let Some(pre) = prebuilt {
                let mut seen = std::collections::HashSet::new();
                for m in pre.ac.find_iter(seq) {
                    let pid = m.pattern();
                    if !seen.insert(pid) { continue; }
                    let ridx = pre.pat2rec[pid];
                    let is_rc = pre.pat_is_rc[pid];
                    let r = &pre.records[ridx];
                    let pat_bytes: Vec<u8> = if is_rc {
                        revcomp_bytes(r.sequence.as_bytes())
                    } else {
                        r.sequence.as_bytes().to_vec()
                    };
                    let mut my: Myers<u64> = MyersBuilder::new().build_64(pat_bytes.into_iter());
                    if let Some((_s, _e, dist)) = my.find_all(seq, max_dist as u8).next() {
                        let _ = dist;
                        out.push((r.name.to_string(), r.kind, is_rc));
                    }
                }
            }
        }
        BenchmarkAlgo::Myers => {
            for r in records {
                // forward
                let mut m: Myers<u64> = MyersBuilder::new().build_64(r.sequence.as_bytes().iter().copied());
                if let Some((_s,_e,dist)) = m.find_all(seq, max_dist as u8).next() {
                    let _ = dist;
                    out.push((r.name.to_string(), r.kind, false));
                    continue;
                }
                // reverse-complement motif
                let rc = revcomp_bytes(r.sequence.as_bytes());
                let mut mrc: Myers<u64> = MyersBuilder::new().build_64(rc.into_iter());
                if let Some((_s,_e,dist)) = mrc.find_all(seq, max_dist as u8).next() {
                    let _ = dist;
                    out.push((r.name.to_string(), r.kind, true));
                }
            }
        }
        BenchmarkAlgo::Edlib => {
            unsafe {
                for r in records {
                    let mut cfg: EdlibAlignConfig = edlibDefaultAlignConfig();
                    cfg.mode = EdlibAlignMode_EDLIB_MODE_HW; // semiglobal (end-free)
                    cfg.task = EdlibAlignTask_EDLIB_TASK_DISTANCE;
                    cfg.k = max_dist as i32;

                    // forward
                    let q = r.sequence.as_bytes();
                    let res = edlibAlign(q.as_ptr() as *const i8, q.len() as i32, seq.as_ptr() as *const i8, seq.len() as i32, cfg);
                    let mut matched = false;
                    if res.editDistance >= 0 {
                        out.push((r.name.to_string(), r.kind, false));
                        matched = true;
                    }
                    edlibFreeAlignResult(res);

                    if !matched {
                        let rc = revcomp_bytes(r.sequence.as_bytes());
                        let res2 = edlibAlign(rc.as_ptr() as *const i8, rc.len() as i32, seq.as_ptr() as *const i8, seq.len() as i32, cfg);
                        if res2.editDistance >= 0 {
                            out.push((r.name.to_string(), r.kind, true));
                        }
                        edlibFreeAlignResult(res2);
                    }
                }
            }
        }
        BenchmarkAlgo::Parasail => {
            // Fallback: behave like Myers (forward + RC motifs), do not RC the read
            for r in records {
                let mut m: Myers<u64> = MyersBuilder::new().build_64(r.sequence.as_bytes().iter().copied());
                if let Some((_s,_e,dist)) = m.find_all(seq, max_dist as u8).next() {
                    let _ = dist;
                    out.push((r.name.to_string(), r.kind, false));
                    continue;
                }
                let rc = revcomp_bytes(r.sequence.as_bytes());
                let mut mrc: Myers<u64> = MyersBuilder::new().build_64(rc.into_iter());
                if let Some((_s,_e,dist)) = mrc.find_all(seq, max_dist as u8).next() {
                    let _ = dist;
                    out.push((r.name.to_string(), r.kind, true));
                }
            }
        }
    }
    out
}
