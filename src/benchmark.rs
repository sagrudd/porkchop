//! Benchmarking framework for adapter/primer/barcode assignment.
//!
//! ## Algorithms
//! - **Aho-Corasick** (exact seeds)
//! - **Myers** (approx/Levenshtein; via `rust-bio`)
//! - **Edlib** (edit distance; via `edlib_rs`; falls back to Myers if unavailable)
//! - **Parasail** (local SW alignment; via `parasailors`)
//! - **Two-stage** pipelines: `ac+myers`, `ac+parasail`
//!
//! ## Truth-set schema (CSV/TSV)
//! Columns (header row required):
//! - `read_id` (string; required)
//! - `expected_labels` (string; semicolon-separated; e.g. `NB01;SSPII`)
//! - `kind` (optional: `adapter` | `primer` | `barcode`) – if present, evaluation is filtered to that class
//!
//! If no truth set is provided, the benchmark still runs and reports timing/throughput
//! and basic classification observations (TP/FP/FN remain NA).
//!
//! ## Notes
//! - Sequences are taken from the kit registry (`crate::kits::KITS`) using the kit id you pass.
//! - CPU utilisation is sampled periodically to produce a coarse mean %.

use std::collections::{HashMap, HashSet};
use std::sync::{Arc, atomic::{AtomicU64, AtomicUsize, Ordering}};
use std::time::{Duration, Instant};

use aho_corasick::AhoCorasick;
use bio::pattern_matching::myers::MyersBuilder;
use csv;
use parasailors as ps;
use sysinfo::{System, RefreshKind, CpuRefreshKind};
use anyhow;

use crate::kit::{SequenceRecord, SeqKind};
use crate::seqio::{self, NARead};

/// Enum of available algorithms.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Algorithm {
    AhoCorasick,
    Myers,
    Edlib,
    Parasail,
    TwoStageACMyers,
    TwoStageACParasail,
}

impl Algorithm {
    pub fn all() -> Vec<Algorithm> {
        vec![
            Self::AhoCorasick,
            Self::Myers,
            Self::Edlib,
            Self::Parasail,
            Self::TwoStageACMyers,
            Self::TwoStageACParasail,
        ]
    }
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::AhoCorasick => "aho",
            Self::Myers => "myers",
            Self::Edlib => "edlib",
            Self::Parasail => "parasail",
            Self::TwoStageACMyers => "ac+myers",
            Self::TwoStageACParasail => "ac+parasail",
        }
    }
    /// Parse `aho,myers,edlib,parasail,ac+myers,ac+parasail`
    pub fn parse_list(s: &str) -> Vec<Algorithm> {
        let mut v = Vec::new();
        for tok in s.split(',').map(|t| t.trim().to_lowercase()) {
            let a = match tok.as_str() {
                "aho" | "ac"       => Some(Self::AhoCorasick),
                "myers"            => Some(Self::Myers),
                "edlib"            => Some(Self::Edlib),
                "parasail" | "sw"  => Some(Self::Parasail),
                "ac+myers"         => Some(Self::TwoStageACMyers),
                "ac+parasail"      => Some(Self::TwoStageACParasail),
                ""                 => None,
                other              => { eprintln!("Unknown algorithm token: {}", other); None }
            };
            if let Some(a) = a { v.push(a); }
        }
        if v.is_empty() { Self::all() } else { v }
    }
}

/// A classification result.
#[derive(Debug, Clone)]
pub struct LabelHit {
    pub name: String,
    pub kind: SeqKind,
    pub score: i32,                  // lower=better for edit distance; higher=better for SW
    pub pos: Option<(usize, usize)>, // (start,end) if known
}

/// Truth-set entry for a read.
#[derive(Debug, Clone)]
pub struct Truth {
    pub expected: HashSet<String>,   // label names
    pub kind: Option<SeqKind>,       // if present, restricts evaluation
}

/// Load truth set: `read_id, expected_labels, [kind]`
pub fn load_truth(path: &str) -> anyhow::Result<HashMap<String, Truth>> {
    let mut out = HashMap::new();
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(true)
        .flexible(true)
        .from_path(path)?;
    for rec in rdr.records() {
        let r = rec?;
        let id = r.get(0).unwrap_or("").to_string();
        let labels = r.get(1)
            .unwrap_or("")
            .split(';')
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty())
            .collect::<HashSet<_>>();
        let kind = match r.get(2).map(|s| s.to_lowercase()) {
            Some(k) if k == "adapter" => Some(SeqKind::AdapterTop), // coarse bucket
            Some(k) if k == "primer"  => Some(SeqKind::Primer),
            Some(k) if k == "barcode" => Some(SeqKind::Barcode),
            _                         => None,
        };
        out.insert(id, Truth { expected: labels, kind });
    }
    Ok(out)
}

/// Gather searchable sequences for a kit (adapters, primers, barcodes).
fn patterns_for_kit(kit_id: &str) -> anyhow::Result<Vec<SequenceRecord>> {
    let kit = crate::get_sequences_for_kit(kit_id)
        .ok_or_else(|| anyhow::anyhow!("Unknown kit id {}", kit_id))?;
    let mut v = Vec::new();
    // adapters & primers
    for r in kit.adapters_and_primers.iter() {
        if matches!(r.kind, SeqKind::AdapterTop | SeqKind::AdapterBottom | SeqKind::Primer) {
            v.push(*r);
        }
    }
    // barcodes
    for r in kit.barcodes.iter() {
        v.push(*r);
    }
    Ok(v)
}

/// AC seed index (exact).
fn ac_index(records: &[SequenceRecord]) -> AhoCorasick {
    let pats: Vec<&[u8]> = records.iter().map(|r| r.sequence.as_bytes()).collect();
    AhoCorasick::new(pats).expect("AC build failed")
}

/// Myers approximate matching (best edit distance up to `max_dist`).
fn myers_best(seq: &[u8], records: &[SequenceRecord], max_dist: usize) -> Option<LabelHit> {
    let mut best: Option<LabelHit> = None;
    for r in records {
        let m = MyersBuilder::new().build_64(&r.sequence.as_bytes());
        if let Some(end) = m.find_all_end(seq, max_dist).next() {
            let dist = end.1 as i32;
            let hit = LabelHit { name: r.name.to_string(), kind: r.kind, score: dist, pos: None };
            if best.as_ref().map(|b| dist < b.score).unwrap_or(true) {
                best = Some(hit);
            }
        }
    }
    best
}

/// Parasail best SW local alignment (identity matrix, affine gaps).
fn parasail_best(seq: &[u8], records: &[SequenceRecord]) -> Option<LabelHit> {
    let mut best: Option<LabelHit> = None;
    for r in records {
        let profile = ps::Profile::new(r.sequence.as_bytes(), &ps::Matrix::identity(), ps::Algorithm::Sw);
        let res = ps::sw_trace(&profile, seq, 5, 1); // open=5, extend=1 – conservative defaults
        let score = res.score;
        let hit = LabelHit { name: r.name.to_string(), kind: r.kind, score, pos: None };
        if best.as_ref().map(|b| score > b.score).unwrap_or(true) {
            best = Some(hit);
        }
    }
    best
}

/// Edlib best edit distance (fallback to Myers if edlib is unusable).
fn edlib_best(seq: &[u8], records: &[SequenceRecord], max_dist: usize) -> Option<LabelHit> {
    // Try edlib first; if any error, just use Myers gracefully.
    #[allow(unused_mut)]
    let mut used_edlib = false;
    let mut best: Option<LabelHit> = None;

    // If edlib works in your env, this block runs; otherwise, any build/run error will
    // surface at compile time and we’ll address it. At runtime errors, we fall back.
    for r in records {
        match edlib_rs::align(
            seq,
            r.sequence.as_bytes(),
            edlib_rs::AlignmentMode::HW,
            edlib_rs::AlignConfig::default(),
        ) {
            Ok(a) => {
                used_edlib = true;
                let dist = a.edit_distance as i32;
                let hit = LabelHit { name: r.name.to_string(), kind: r.kind, score: dist, pos: None };
                if best.as_ref().map(|b| dist < b.score).unwrap_or(true) {
                    best = Some(hit);
                }
            }
            Err(_) => { /* keep trying others; if nothing hits, we'll fall back */ }
        }
    }
    if used_edlib { best } else { myers_best(seq, records, max_dist) }
}

/// Run a classifier once on a read.
fn classify_once(
    algo: Algorithm,
    seq: &[u8],
    ac: Option<&AhoCorasick>,
    records: &[SequenceRecord],
    max_dist: usize,
) -> Option<LabelHit> {
    match algo {
        Algorithm::AhoCorasick => {
            if let Some(ac) = ac {
                if let Some(m) = ac.find(seq) {
                    let r = &records[m.pattern()];
                    return Some(LabelHit { name: r.name.to_string(), kind: r.kind, score: 0, pos: Some((m.start(), m.end())) });
                }
            }
            None
        }
        Algorithm::Myers => myers_best(seq, records, max_dist),
        Algorithm::Edlib => edlib_best(seq, records, max_dist),
        Algorithm::Parasail => parasail_best(seq, records),
        Algorithm::TwoStageACMyers => {
            if let Some(ac) = ac {
                if let Some(m) = ac.find(seq) {
                    let r = &records[m.pattern()];
                    return myers_best(seq, std::slice::from_ref(r), max_dist);
                }
            }
            myers_best(seq, records, max_dist)
        }
        Algorithm::TwoStageACParasail => {
            if let Some(ac) = ac {
                if let Some(m) = ac.find(seq) {
                    let r = &records[m.pattern()];
                    return parasail_best(seq, std::slice::from_ref(r));
                }
            }
            parasail_best(seq, records)
        }
    }
}

/// Benchmark one file × one algorithm.
/// Returns (tp, fp, fn, elapsed, nseq, cpu_mean_pct).
pub fn benchmark_file(
    path: &str,
    kit_id: &str,
    algo: Algorithm,
    truth: Option<&HashMap<String, Truth>>,
    threads: Option<usize>,
    max_dist: usize,
) -> anyhow::Result<(u64, u64, u64, Duration, usize, f32)> {
    let patterns = Arc::new(patterns_for_kit(kit_id)?);

    let ac = Some(ac_index(&patterns));           // Option<AhoCorasick>

    //let ac = AhoCorasick::new(patterns.iter().map(|r| r.sequence.as_bytes())).ok();
    //let ac = ac.as_ref();

    let truth_map = truth.cloned().map(Arc::new);

    // atomics for thread-safe increment in seqio callback
    let tp = AtomicU64::new(0);
    let fp = AtomicU64::new(0);
    let fn_ = AtomicU64::new(0);
    let nseq = AtomicUsize::new(0);

    // background CPU sampler (coarse mean %)
    let sampling = std::sync::Arc::new(std::sync::atomic::AtomicBool::new(true));
    let flag = sampling.clone();
    let mut cpu_samples = Vec::<f32>::new();
    let cpu_samples_ptr = std::sync::Arc::new(std::sync::Mutex::new(Vec::<f32>::new()));
    let cpu_samples_ptr_cl = cpu_samples_ptr.clone();

    let sampler = std::thread::spawn(move || {
        let mut sys = System::new_with_specifics(
            RefreshKind::new().with_cpu(CpuRefreshKind::everything())
        );
        while flag.load(Ordering::Relaxed) {
            std::thread::sleep(Duration::from_millis(200));
            sys.refresh_cpu();
            let v = sys.global_cpu_usage();
            if let Ok(mut g) = cpu_samples_ptr_cl.lock() {
                g.push(v);
            }
        }
    });

    let start = Instant::now();
    let (_fmt, processed) = seqio::for_each_parallel(path, threads, |r: NARead| {
        nseq.fetch_add(1, Ordering::Relaxed);
        let hit = classify_once(algo, &r.seq, ac, &patterns, max_dist);
        if let (Some(tmap), Some(h)) = (truth_map.as_ref(), hit.as_ref()) {
            if let Some(t) = tmap.get(&r.id) {
                if t.expected.contains(&h.name) {
                    tp.fetch_add(1, Ordering::Relaxed);
                } else {
                    fp.fetch_add(1, Ordering::Relaxed);
                }
            }
        } else if truth_map.is_some() && hit.is_none() {
            fn_.fetch_add(1, Ordering::Relaxed);
        }
    })?;
    let elapsed = start.elapsed();

    // stop sampler and compute mean CPU
    sampling.store(false, Ordering::Relaxed);
    let _ = sampler.join();
    let cpu_mean = {
        if let Ok(g) = cpu_samples_ptr.lock() {
            if g.is_empty() { 0.0 } else { g.iter().copied().sum::<f32>() / (g.len() as f32) }
        } else { 0.0 }
    };

    Ok((tp.load(Ordering::Relaxed),
        fp.load(Ordering::Relaxed),
        fn_.load(Ordering::Relaxed),
        elapsed,
        processed,
        cpu_mean))
}
