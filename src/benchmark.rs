
//! Benchmarking framework for adapter/primer/barcode assignment.
//!
//! Algorithms: Aho-Corasick, Myers, Edlib, Parasail, plus two-stage AC+{Myers,Parasail}.
//! Truth set CSV: read_id, expected_labels (;-separated), [kind].

use std::collections::{HashMap, HashSet};
use std::sync::{Arc, atomic::{AtomicU64, AtomicUsize, Ordering}};
use std::time::{Duration, Instant};

use anyhow::Result;
use aho_corasick::AhoCorasick;
use bio::pattern_matching::myers::MyersBuilder;
use bio::alignment::pairwise::Aligner;
use csv;
use edlib_rs::edlibrs::{edlibAlignRs, edlibDefaultAlignConfig, EdlibAlignMode_EDLIB_MODE_HW, EdlibAlignTask_EDLIB_TASK_DISTANCE};
use sysinfo::{System, RefreshKind, CpuRefreshKind};

use crate::kit::{SequenceRecord, SeqKind};
use crate::seqio::{self, NARead};

         // Prebuilt search objects for faster per-read classification
         struct Prebuilt<'a> {
             records: &'a [SequenceRecord],
             myers: Vec<bio::pattern_matching::myers::Myers<u64>>,
         }
         impl<'a> Prebuilt<'a> {
                fn build(records: &'a [SequenceRecord]) -> Self {
                         let myers = records.iter()
                             .map(|r| bio::pattern_matching::myers::MyersBuilder::new().build_64(&r.sequence.as_bytes()))
                             .collect();
                         Self { records, myers }
                     }
         }


#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Algorithm { AhoCorasick, Myers, Edlib, Parasail, TwoStageACMyers, TwoStageACParasail }
impl Algorithm {
    pub fn all() -> Vec<Algorithm> { vec![Self::AhoCorasick, Self::Myers, Self::Edlib, Self::Parasail, Self::TwoStageACMyers, Self::TwoStageACParasail] }
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
    pub fn parse_list(s: &str) -> Vec<Algorithm> {
        let mut v = Vec::new();
        for tok in s.split(',').map(|t| t.trim().to_lowercase()) {
            let a = match tok.as_str() {
                "aho" | "ac" => Some(Self::AhoCorasick),
                "myers" => Some(Self::Myers),
                "edlib" => Some(Self::Edlib),
                "parasail" | "sw" => Some(Self::Parasail),
                "ac+myers" => Some(Self::TwoStageACMyers),
                "ac+parasail" => Some(Self::TwoStageACParasail),
                "" => None,
                _ => None,
            };
            if let Some(a) = a { v.push(a); }
        }
        if v.is_empty() { Self::all() } else { v }
    }
}

#[derive(Debug, Clone)]
pub struct LabelHit { pub name: String, pub kind: SeqKind, pub score: i32, pub pos: Option<(usize,usize)> }

#[derive(Debug, Clone)]
pub struct Truth { pub expected: HashSet<String>, pub kind: Option<SeqKind> }

pub fn load_truth(path: &str) -> Result<HashMap<String, Truth>> {
    let mut out = HashMap::new();
    let mut rdr = csv::ReaderBuilder::new().has_headers(true).flexible(true).from_path(path)?;
    for rec in rdr.records() {
        let r = rec?;
        let id = r.get(0).unwrap_or("").to_string();
        let labels = r.get(1).unwrap_or("").split(';').map(|s| s.trim().to_string()).filter(|s| !s.is_empty()).collect::<HashSet<_>>();
        let kind = match r.get(2).map(|s| s.to_lowercase()) {
            Some(k) if k == "adapter" => Some(SeqKind::AdapterTop),
            Some(k) if k == "primer"  => Some(SeqKind::Primer),
            Some(k) if k == "barcode" => Some(SeqKind::Barcode),
            _ => None,
        };
        out.insert(id, Truth { expected: labels, kind });
    }
    Ok(out)
}

fn patterns_for_kit(kit_id: &str) -> Result<Vec<SequenceRecord>> {
    let kit = crate::get_sequences_for_kit(kit_id).ok_or_else(|| anyhow::anyhow!("Unknown kit id {}", kit_id))?;
    let mut v = Vec::new();
    for r in kit.adapters_and_primers.iter() {
        if matches!(r.kind, SeqKind::AdapterTop | SeqKind::AdapterBottom | SeqKind::Primer) { v.push(*r); }
    }
    for r in kit.barcodes.iter() { v.push(*r); }
    Ok(v)
}

fn ac_index(records: &[SequenceRecord]) -> AhoCorasick {
    let pats: Vec<&[u8]> = records.iter().map(|r| r.sequence.as_bytes()).collect();
    AhoCorasick::new(pats).expect("AC build")
}

fn myers_best(seq: &[u8], records: &[SequenceRecord], max_dist: usize) -> Option<LabelHit> {
    let mut best: Option<LabelHit> = None;
    for r in records {
        let m = MyersBuilder::new().build_64(&r.sequence.as_bytes());
        if let Some(end) = m.find_all_end(seq, max_dist).next() {
            let dist = end.1 as i32;
            let hit = LabelHit { name: r.name.to_string(), kind: r.kind, score: dist, pos: None };
            if best.as_ref().map(|b| dist < b.score).unwrap_or(true) { best = Some(hit); }
        }
    }
    best
}

fn parasail_best(seq: &[u8], records: &[SequenceRecord]) -> Option<LabelHit> { sw_best(seq, records) }

fn sw_best(seq: &[u8], records: &[SequenceRecord]) -> Option<LabelHit> {
    let mut best: Option<LabelHit> = None;
    let mut aligner = Aligner::new(-5, -1, |a: u8, b: u8| if a == b { 2 } else { -1 });
    for r in records {
        let score = aligner.local(r.sequence.as_bytes(), seq).score;
        let hit = LabelHit { name: r.name.to_string(), kind: r.kind, score, pos: None };
        if best.as_ref().map(|b| score > b.score).unwrap_or(true) { best = Some(hit); }
    }
    best
}

fn sw_best_prebuilt(seq: &[u8], pre: &Prebuilt) -> Option<LabelHit> {
    let mut best: Option<LabelHit> = None;
    let mut aligner = Aligner::new(-5, -1, |a: u8, b: u8| if a == b { 2 } else { -1 });
    for r in pre.records {
        let score = aligner.local(r.sequence.as_bytes(), seq).score;
        let hit = LabelHit { name: r.name.to_string(), kind: r.kind, score, pos: None };
        if best.as_ref().map(|b| score > b.score).unwrap_or(true) { best = Some(hit); }
    }
    best
}



fn edlib_best(seq: &[u8], records: &[SequenceRecord], max_dist: usize) -> Option<LabelHit> {
    let mut best: Option<LabelHit> = None;
    let mut any = false;
    for r in records {
        let mut cfg = unsafe { edlibDefaultAlignConfig() };
        cfg.k = max_dist as i32;
        cfg.mode = EdlibAlignMode_EDLIB_MODE_HW; // "HW" / infix
        cfg.task = EdlibAlignTask_EDLIB_TASK_DISTANCE; // distance only
        // Align pattern against read; distance is symmetric for our purpose.
        let res = unsafe { edlibAlignRs(r.sequence.as_bytes(), seq, cfg) };
        if res.editDistance >= 0 {
            any = true;
            let dist = res.editDistance as i32;
            let hit = LabelHit { name: r.name.to_string(), kind: r.kind, score: dist, pos: None };
            if best.as_ref().map(|b| dist < b.score).unwrap_or(true) { best = Some(hit); }
        }
    }
    if any { best } else { None }
}


fn classify_once(algo: Algorithm, seq: &[u8], ac: Option<&AhoCorasick>, records: &[SequenceRecord], pre: Option<&Prebuilt>, max_dist: usize) -> Option<LabelHit> {
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
        Algorithm::Myers => pre.map(|p| myers_best_prebuilt(seq, p, max_dist)).unwrap_or_else(|| myers_best(seq, records, max_dist)),
        Algorithm::Edlib => edlib_best(seq, records, max_dist),
        
        Algorithm::Parasail => pre.map(|p| parasail_best_prebuilt(seq, p)).unwrap_or_else(|| parasail_best(seq, records)),
Algorithm::TwoStageACParasail => {
    if let Some(ac) = ac {
        if let Some(m) = ac.find(seq) {
            let r = &records[m.pattern()];
            return parasail_best(seq, std::slice::from_ref(r));
        }
    }
    parasail_best(seq, records)
}
        
        Algorithm::TwoStageACMyers => {
            if let Some(ac) = ac {
                if let Some(m) = ac.find(seq) {
                    let r = &records[m.pattern()];
                    return pre.map(|p| { let sub = Prebuilt{ records: std::slice::from_ref(r), myers: vec![p.myers[ m.pattern() ].clone()], parasail: vec![p.parasail[ m.pattern() ].clone()] }; myers_best_prebuilt(seq, &sub, max_dist) }).unwrap_or_else(|| myers_best(seq, std::slice::from_ref(r), max_dist));
                }
            }
            myers_best(seq, records, max_dist)
        }
    
    }
}

pub fn benchmark_file(path: &str, kit_id: &str, algo: Algorithm, truth: Option<&HashMap<String, Truth>>, threads: Option<usize>, max_dist: usize)
-> Result<(u64, u64, u64, Duration, usize, f32)>
{
    let patterns = Arc::new(patterns_for_kit(kit_id)?);
    let ac_idx = ac_index(&patterns);
    let ac = Some(&ac_idx);
    let prebuilt = Prebuilt::build(&patterns);

    let truth_map = truth.cloned().map(Arc::new);

    let tp = AtomicU64::new(0);
    let fp = AtomicU64::new(0);
    let fn_ = AtomicU64::new(0);
    let nseq = AtomicUsize::new(0);

    let sampling = std::sync::Arc::new(std::sync::atomic::AtomicBool::new(true));
    let flag = sampling.clone();
    let cpu_samples = std::sync::Arc::new(std::sync::Mutex::new(Vec::<f32>::new()));
    let cpu_samples_cl = cpu_samples.clone();

    let sampler = std::thread::spawn(move || {
        let mut sys = System::new_with_specifics(RefreshKind::new().with_cpu(CpuRefreshKind::everything()));
        while flag.load(Ordering::Relaxed) {
            std::thread::sleep(Duration::from_millis(200));
            sys.refresh_cpu();
            let v = sys.global_cpu_usage();
            if let Ok(mut g) = cpu_samples_cl.lock() { g.push(v); }
        }
    });

    let start = Instant::now();
    let (_fmt, _processed) = seqio::for_each_parallel(path, threads, |r: NARead| {
        nseq.fetch_add(1, Ordering::Relaxed);
        let hit = classify_once(algo, &r.seq, ac, &patterns, Some(&prebuilt), max_dist);
        if let (Some(tmap), Some(h)) = (truth_map.as_ref(), hit.as_ref()) {
            if let Some(t) = tmap.get(&r.id) {
                if t.expected.contains(&h.name) { tp.fetch_add(1, Ordering::Relaxed); }
                else { fp.fetch_add(1, Ordering::Relaxed); }
            }
        } else if truth_map.is_some() && hit.is_none() {
            fn_.fetch_add(1, Ordering::Relaxed);
        }
    })?;
    let elapsed = start.elapsed();

    sampling.store(false, Ordering::Relaxed);
    let _ = sampler.join();
    let cpu_mean = {
        if let Ok(g) = cpu_samples.lock() {
            if g.is_empty() { 0.0 } else { g.iter().copied().sum::<f32>() / (g.len() as f32) }
        } else { 0.0 }
    };

    Ok((tp.load(Ordering::Relaxed), fp.load(Ordering::Relaxed), fn_.load(Ordering::Relaxed), elapsed, nseq.load(Ordering::Relaxed), cpu_mean))
}


fn myers_best_prebuilt(seq: &[u8], pre: &Prebuilt, max_dist: usize) -> Option<LabelHit> {
    let mut best: Option<LabelHit> = None;
    for (i, m) in pre.myers.iter().enumerate() {
        if let Some(end) = m.find_all_end(seq, max_dist).next() {
            let dist = end.1 as i32;
            let r = &pre.records[i];
            let hit = LabelHit { name: r.name.to_string(), kind: r.kind, score: dist, pos: None };
            if best.as_ref().map(|b| dist < b.score).unwrap_or(true) { best = Some(hit); }
        }
    }
    best
}

fn parasail_best_prebuilt(seq: &[u8], pre: &Prebuilt) -> Option<LabelHit> { sw_best_prebuilt(seq, pre) }
