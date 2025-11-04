use std::collections::{HashMap, HashSet};
use std::sync::{Arc, Mutex, atomic::{AtomicBool, AtomicUsize, Ordering}};
use std::time::{Duration, Instant};

use crate::benchmark::{self, BenchmarkAlgo};
use crate::kit::SeqKind;
use crate::list_supported_kits;
use crate::seqio::{for_each_parallel, NARead};
use polars::prelude::*;
use rayon::ThreadPoolBuilder;
use std::sync::mpsc;
fn canonical_barcode(name: &str) -> Option<String> {
    // Map tokens like "BP05/BC05/RB05/16S05/RLB05" -> "BC05"
    for tok in name.split('/') {
        if tok.len() == 4 {
            let (pre, num) = tok.split_at(2);
            if (pre == "BP" || pre == "BC" || pre == "RB" || pre == "16" || pre == "RL") && num.chars().all(|c| c.is_ascii_digit()) {
                return Some(format!("BC{}", num));
            }
        } else if tok.starts_with("BC") && tok[2..].chars().all(|c| c.is_ascii_digit()) {
            return Some(tok.to_string());
        }
    }
    None
}

fn weight_of(kind: SeqKind) -> f64 {
    match kind {
        SeqKind::AdapterTop | SeqKind::AdapterBottom => 3.0,
        SeqKind::Primer => 2.0,
        SeqKind::Barcode => 1.0,
        SeqKind::Flank => 0.5,
    }
}

fn infer_kits_df(tally: &std::collections::HashMap<(String, SeqKind), usize>) -> polars::prelude::PolarsResult<DataFrame> {
    use std::collections::{HashMap, HashSet};
    use std::cmp::Ordering;

    let kits = crate::list_supported_kits();
    let mut rows: Vec<(String, String, String, f64, f64, usize, usize)> = Vec::new();
    let mut scores: Vec<f64> = Vec::new();

    let total_hits: usize = tally.values().copied().sum();

    for k in kits {
        let mut sig: HashSet<(String, SeqKind)> = HashSet::new();
        for r in k.adapters_and_primers {
            sig.insert((r.name.to_string(), r.kind));
        }
        for r in k.barcodes {
            let nm = canonical_barcode(r.name).unwrap_or_else(|| r.name.to_string());
            sig.insert((nm, SeqKind::Barcode));
            if matches!(r.kind, SeqKind::Flank) {
                sig.insert((r.name.to_string(), SeqKind::Flank));
            }
        }

        let mut score = 0.0f64;
        let mut matched = 0usize;
        for ((nm, kind), cnt) in tally.iter() {
            let key = if *kind == SeqKind::Barcode {
                (canonical_barcode(nm).unwrap_or_else(|| nm.clone()), *kind)
            } else {
                (nm.clone(), *kind)
            };
            if sig.contains(&key) {
                matched += 1;
                score += weight_of(*kind) * (*cnt as f64);
            }
        }
        scores.push(score);

        rows.push((
            k.id.0.to_string(),
            k.description.to_string(),
            k.chemistry.to_string(),
            score,
            0.0, // prob, filled below
            matched,
            total_hits,
        ));
    }

    // Softmax probabilities
    let max_s = scores.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let exps: Vec<f64> = scores.iter().map(|s| (s - max_s).exp()).collect();
    let z: f64 = exps.iter().sum::<f64>().max(1e-12);
    let probs: Vec<f64> = exps.iter().map(|e| e / z).collect();

    for (i, row) in rows.iter_mut().enumerate() {
        row.4 = probs[i];
    }

    // Sort rows by probability desc
    rows.sort_by(|a, b| b.4.partial_cmp(&a.4).unwrap_or(Ordering::Equal));

    // Build DataFrame via df! macro
    let kit_v: Vec<String> = rows.iter().map(|r| r.0.clone()).collect();
    let desc_v: Vec<String> = rows.iter().map(|r| r.1.clone()).collect();
    let chem_v: Vec<String> = rows.iter().map(|r| r.2.clone()).collect();
    let score_v: Vec<f64> = rows.iter().map(|r| r.3).collect();
    let prob_v: Vec<f64> = rows.iter().map(|r| r.4).collect();
    let matched_v: Vec<u64> = rows.iter().map(|r| r.5 as u64).collect();
    let total_v: Vec<u64> = rows.iter().map(|r| r.6 as u64).collect();

    let df = df!(
        "kit"             => kit_v,
        "description"     => desc_v,
        "chemistry"       => chem_v,
        "score"           => score_v,
        "probability"     => prob_v,
        "matched_motifs"  => matched_v,
        "total_hits"      => total_v,
    )?;

    Ok(df)
}


// TUI
use crossterm::event::{self, Event, KeyCode};
use crossterm::terminal::{disable_raw_mode, enable_raw_mode};
use ratatui::prelude::*;
use ratatui::widgets::{Block, Borders, Row, Table};
use serde_json;

#[derive(Debug, Clone)]
pub struct ScreenOpts {
    pub files: Vec<String>,
    pub threads: Option<usize>,
    pub fraction: f64,
    pub tick_secs: u64,
    pub algo: BenchmarkAlgo,
    pub max_dist: usize,
    pub json: Option<String>,
    pub kit_prob_min: f64,
}

fn collect_all_sequences() -> Vec<crate::kit::SequenceRecord> {
    let mut v = Vec::new();
    for k in list_supported_kits() {
        v.extend_from_slice(k.adapters_and_primers);
        v.extend_from_slice(k.barcodes);
    }
    v
}

fn kind_suffix(k: SeqKind) -> &'static str {
    match k {
        SeqKind::Primer => "target",
        SeqKind::Flank => "flank",
        SeqKind::Barcode => "barcode",
        SeqKind::AdapterTop | SeqKind::AdapterBottom => "adapter",
    }
}

fn revcomp(s: &str) -> String {
    fn comp(c: u8) -> u8 {
        match c {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'U' | b'u' => b'A',
            b'R' | b'r' => b'N',
            b'Y' | b'y' => b'N',
            b'S' | b's' => b'N',
            b'W' | b'w' => b'N',
            b'K' | b'k' => b'N',
            b'M' | b'm' => b'N',
            b'B' | b'b' => b'N',
            b'D' | b'd' => b'N',
            b'H' | b'h' => b'N',
            b'V' | b'v' => b'N',
            _ => b'N',
        }
    }
    let bytes = s.as_bytes();
    let mut out = Vec::with_capacity(bytes.len());
    for &b in bytes.iter().rev() {
        out.push(comp(b));
    }
    String::from_utf8(out).unwrap_or_default()
}


pub fn run_screen(opts: ScreenOpts) -> anyhow::Result<()> {
    let records = Arc::new(collect_all_sequences());
    let rec_map: std::collections::HashMap<String, String> = records.iter()
        .map(|r| (r.name.to_string(), r.sequence.to_string()))
        .collect();
    let rec_map = Arc::new(rec_map);

    // Tallies
    let unit_tally: Arc<Mutex<HashMap<(String, SeqKind), usize>>> = Arc::new(Mutex::new(HashMap::new()));
    let combo_tally: Arc<Mutex<HashMap<String, usize>>> = Arc::new(Mutex::new(HashMap::new()));
    let fwd_tally: Arc<Mutex<HashMap<(String, SeqKind), usize>>> = Arc::new(Mutex::new(HashMap::new()));
    let rev_tally: Arc<Mutex<HashMap<(String, SeqKind), usize>>> = Arc::new(Mutex::new(HashMap::new()));
    let done = Arc::new(AtomicBool::new(false));
    let screened = Arc::new(AtomicUsize::new(0));
    let unclassified = Arc::new(AtomicUsize::new(0));
    let skipped = Arc::new(AtomicUsize::new(0));
    let reads_with_hits = Arc::new(AtomicUsize::new(0));

    // Optional prebuilt for ACMyers
    let prebuilt = if let BenchmarkAlgo::ACMyers = opts.algo {
        Some(Arc::new(benchmark::prebuild_for(records.as_slice())))
    } else {
        None
    };

    // UI thread
    let unit_ui = unit_tally.clone();
    let fwd_ui = fwd_tally.clone();
    let rev_ui = rev_tally.clone();
    let combo_ui = combo_tally.clone();
    let done_ui = done.clone();
    let screened_ui = screened.clone();
    let unclassified_ui = unclassified.clone();
    let skipped_ui = skipped.clone();
    let tick = Duration::from_secs(opts.tick_secs.max(1));
let rwh_ui = reads_with_hits.clone();
    let mut ui_handle_opt: Option<std::thread::JoinHandle<()>> = Some(std::thread::spawn(move || {
        let _ = tui_loop(unit_ui, fwd_ui, rev_ui, combo_ui, done_ui, screened_ui, unclassified_ui, skipped_ui, rwh_ui, rec_map.clone(), tick);
}));
// Sampling params
    let p = opts.fraction.clamp(0.0, 1.0);
    let threads = opts.threads;


    // Build a dedicated Rayon pool for classification
    let threads_n = opts.threads.unwrap_or_else(num_cpus::get).max(1);
    let pool = ThreadPoolBuilder::new().num_threads(threads_n).build()?;

    // Bounded work queue to decouple IO from CPU
    let (tx, rx) = mpsc::sync_channel::<NARead>(threads_n * 1024);
    let rx = Arc::new(Mutex::new(rx));
    // Spawn a producer per file (reading thread); keep IO single-threaded inside each reader,
    // classification happens in the Rayon pool below.
    let mut producers = Vec::new();
    for file in &opts.files {
        let file = file.clone();
        let tx_p = tx.clone();
        let done_p = done.clone();
        let skipped_p = skipped.clone();
        let p_sample = p;
        producers.push(std::thread::spawn(move || {
            let _ = for_each_parallel(file, Some(1), move |read| {
                if done_p.load(Ordering::SeqCst) { return; }
                // Bernoulli(p) sampling via deterministic hash of read id
                let take = if p_sample >= 1.0 {
                    true
                } else {
                    let mut h: u64 = 0xcbf29ce484222325;
                    for b in read.id.as_bytes() { h = h.wrapping_mul(0x100000001b3) ^ (*b as u64); }
                    (h as f64 / std::u64::MAX as f64) < p_sample
                };
                if !take { skipped_p.fetch_add(1, Ordering::Relaxed); return; }
                if tx_p.send(read).is_err() {
                    // channel closed; stop producing
                }
            });
        }));
    }
    drop(tx); // close when producers join

    // Start a minimal Tokio runtime for ticking in TUI (interval sleeps)

    // Spawn classification workers in Rayon pool
    let unit_w = unit_tally.clone();
    let fwd_w = fwd_tally.clone();
    let rev_w = rev_tally.clone();
    let combo_w = combo_tally.clone();
    let rwh = reads_with_hits.clone();
    let records_arc = records.clone();
    let prebuilt_c = prebuilt.clone();
    let algo = opts.algo;
    let max_dist = opts.max_dist;
    let screened_c = screened.clone();
    let unclassified_c = unclassified.clone();

    pool.install(|| {
        rayon::scope(|s| {
            for _ in 0..threads_n {
                let rx_c = rx.clone();  // Arc<Mutex<Receiver<NARead>>>
                let unit_wc = unit_w.clone();
                let fwd_wc = fwd_w.clone();
                let rev_wc = rev_w.clone();
                let combo_wc = combo_w.clone();
                let records_c = records_arc.clone();
                let pre_c = prebuilt_c.clone();
                let done_c = done.clone();
                let screened_wc = screened_c.clone();
                let unclassified_wc = unclassified_c.clone();
let rwh = reads_with_hits.clone();
                s.spawn(move |_| {
                    loop {
                        let read = { let guard = rx_c.lock().unwrap(); guard.recv() };
                        let read = match read { Ok(r) => r, Err(_) => break };
                        if done_c.load(Ordering::SeqCst) { break; }
                        // Enumerate all motif hits for this read using requested algorithm
                        let hits = benchmark::classify_all(algo, &read.seq, records_c.as_slice(), pre_c.as_deref(), max_dist);

                        if hits.is_empty() {
                            screened_wc.fetch_add(1, Ordering::Relaxed);
                            unclassified_wc.fetch_add(1, Ordering::Relaxed);
                            continue;
                        }

                        // Tally individual hits (dedupe per read by (name, kind))
                        rwh.fetch_add(1, Ordering::Relaxed);
                        {
                            let mut uniq = std::collections::HashSet::new();
                            for (name, kind, _is_rc, _pos) in &hits {
                                uniq.insert((name.clone(), *kind));
                            }
                            let mut g = unit_wc.lock().unwrap();
                            for key in uniq {
                                *g.entry(key).or_insert(0) += 1;
                            }
                        }

                        // Strand-specific counters
                        {
                            let mut gf = fwd_wc.lock().unwrap();
                            let mut gr = rev_wc.lock().unwrap();
                            for (name, kind, is_rc, _pos) in &hits {
                                if *is_rc { *gr.entry((name.clone(), *kind)).or_insert(0) += 1; }
                                else { *gf.entry((name.clone(), *kind)).or_insert(0) += 1; }
                            }
                        }

                        // Compose aggregate identifier for this read ordered by position
                        let mut labels_pos: Vec<(usize, String)> = Vec::new();
                        let mut seen = HashSet::new();
                        for (name, _kind, _is_rc, pos) in hits {
                            if seen.insert(name.clone()) {
                                labels_pos.push((pos, name));
                            }
                        }
                        labels_pos.sort_by_key(|(pos, _)| *pos);
                        let id = labels_pos.into_iter().map(|(_, nm)| nm).collect::<Vec<_>>().join(" + ");
                        {
                            let mut g = combo_wc.lock().unwrap();
                            *g.entry(id).or_insert(0) += 1;
                        }

                        screened_wc.fetch_add(1, Ordering::Relaxed);
                    }
                });
            }
        });
    });

    // Ensure producers finished
    for jh in producers { let _ = jh.join(); }



    done.store(true, Ordering::SeqCst);
    std::thread::sleep(Duration::from_millis(150));
    
    if let Some(path) = &opts.json {
        // contexts (aggregate identifiers)
        let entries: Vec<(String, usize)> = {
            let g = combo_tally.lock().unwrap();
            let mut v: Vec<(String, usize)> = g.iter().map(|(k,v)| (k.clone(), *v)).collect();
            v.sort_by(|a,b| b.1.cmp(&a.1));
            v
        };
        let contexts: Vec<serde_json::Value> = entries.into_iter()
            .map(|(id, count)| serde_json::json!({"id": id, "count": count}))
            .collect();

        // predicted kits table -> JSON
        // predicted kits table -> JSON
        let kits_json: Vec<serde_json::Value> = if let Ok(unit_map) = unit_tally.lock() {
            if let Ok(df) = infer_kits_df(&*unit_map) {
                use polars::prelude::AnyValue;
                let c_kit = df.column("kit").ok();
                let c_desc = df.column("description").ok();
                let c_chem = df.column("chemistry").ok();
                let c_score = df.column("score").ok();
                let c_prob = df.column("probability").ok();
                let c_matched = df.column("matched_motifs").ok();
                let c_total = df.column("total_hits").ok();
                let mut out = Vec::new();
                for i in 0..df.height() {
                    let kit_av = c_kit.as_ref().and_then(|s| s.get(i).ok()).unwrap_or(AnyValue::Null);
                    let desc_av = c_desc.as_ref().and_then(|s| s.get(i).ok()).unwrap_or(AnyValue::Null);
                    let chem_av = c_chem.as_ref().and_then(|s| s.get(i).ok()).unwrap_or(AnyValue::Null);
                    let score_av = c_score.as_ref().and_then(|s| s.get(i).ok()).unwrap_or(AnyValue::Null);
                    let prob_av = c_prob.as_ref().and_then(|s| s.get(i).ok()).unwrap_or(AnyValue::Null);
                    let matched_av = c_matched.as_ref().and_then(|s| s.get(i).ok()).unwrap_or(AnyValue::Null);
                    let total_av = c_total.as_ref().and_then(|s| s.get(i).ok()).unwrap_or(AnyValue::Null);
                    let kit = kit_av.to_string();
                    let desc = desc_av.to_string();
                    let chem = chem_av.to_string();
                    let score_f = match score_av {
                        AnyValue::Float64(v) => v,
                        AnyValue::Float32(v) => v as f64,
                        AnyValue::Int64(v) => v as f64,
                        AnyValue::Int32(v) => v as f64,
                        AnyValue::UInt64(v) => v as f64,
                        AnyValue::UInt32(v) => v as f64,
                        _ => 0.0,
                    };
                    let prob_f = match prob_av {
                        AnyValue::Float64(v) => v,
                        AnyValue::Float32(v) => v as f64,
                        AnyValue::Int64(v) => v as f64,
                        AnyValue::Int32(v) => v as f64,
                        AnyValue::UInt64(v) => v as f64,
                        AnyValue::UInt32(v) => v as f64,
                        _ => 0.0,
                    };
                    let matched_u = match matched_av {
                        AnyValue::UInt64(v) => v,
                        AnyValue::UInt32(v) => v as u64,
                        AnyValue::Int64(v) => v as u64,
                        AnyValue::Int32(v) => v as u64,
                        _ => 0,
                    };
                    let total_u = match total_av {
                        AnyValue::UInt64(v) => v,
                        AnyValue::UInt32(v) => v as u64,
                        AnyValue::Int64(v) => v as u64,
                        AnyValue::Int32(v) => v as u64,
                        _ => 0,
                    };
                    out.push(serde_json::json!({
                        "kit": kit,
                        "description": desc,
                        "chemistry": chem,
                        "score": score_f,
                        "probability": prob_f,
                        "matched_motifs": matched_u,
                        "total_hits": total_u
                    }));
                }
                out
            } else { Vec::new() }
        } else { Vec::new() };

        // write a single object combining both sections
        let combined = serde_json::json!({
            "contexts": contexts,
            "kits": kits_json
        });
        let mut f = std::fs::File::create(path)?;

    // Ensure the TUI is fully torn down before printing tables (idempotent)
    done.store(true, Ordering::SeqCst);
    if let Some(h) = ui_handle_opt.take() {
        let _ = h.join();
    }
    let _ = crossterm::terminal::disable_raw_mode();
    let _ = crossterm::execute!(std::io::stdout(), crossterm::cursor::Show, crossterm::terminal::LeaveAlternateScreen);
serde_json::to_writer_pretty(&mut f, &combined)?;
    }

    // Also print the kit-likelihood table to stdout as a wide Polars DataFrame
    if let Ok(unit_map) = unit_tally.lock() {
        if let Ok(df) = infer_kits_df(&*unit_map) {
            // Ensure full-width display and no truncation for Polars 0.42
            std::env::set_var("POLARS_FMT_TABLE_FORMATTING", "UTF8_FULL");
            std::env::set_var("POLARS_FMT_MAX_COLS", "100000");
            std::env::set_var("POLARS_FMT_MAX_ROWS", "1000000");
            std::env::set_var("POLARS_FMT_STR_LEN", "1000000");
            std::env::set_var("POLARS_TABLE_WIDTH", "65535");

    // Ensure the TUI is fully torn down before printing tables (idempotent)
    done.store(true, Ordering::SeqCst);
    if let Some(h) = ui_handle_opt.take() {
        let _ = h.join();
    }
    let _ = crossterm::terminal::disable_raw_mode();
    let _ = crossterm::execute!(std::io::stdout(), crossterm::cursor::Show, crossterm::terminal::LeaveAlternateScreen);
println!("\n=== Sequencing kit prediction (p > {:.3}) ===", opts.kit_prob_min);
// Filter to only show kits with probability above user threshold
let df = match df.column("probability") {
    Ok(prob) => {
        if let Ok(prob) = prob.f64() {
            let mask = prob.gt(opts.kit_prob_min);
            df.filter(&mask).unwrap_or(df.clone())
        } else {
            df.clone()
        }
    }
    Err(_) => df.clone(),
};
// Drop the verbose kit description column from the console view
let df = match df.drop("description") {
    Ok(d) => d,
    Err(_) => df,
};
println!("{}", df);
        }
    }

    Ok(())
}


fn tui_loop(
    unit: Arc<Mutex<HashMap<(String, SeqKind), usize>>>,
    fwd: Arc<Mutex<HashMap<(String, SeqKind), usize>>>,
    rev: Arc<Mutex<HashMap<(String, SeqKind), usize>>>,
    combos: Arc<Mutex<HashMap<String, usize>>>,
    done: Arc<AtomicBool>,
    screened: Arc<AtomicUsize>,
    unclassified: Arc<AtomicUsize>,
    skipped: Arc<AtomicUsize>,
    reads_with_hits: Arc<AtomicUsize>,
    rec_map: Arc<HashMap<String, String>>,
    tick: Duration,
) -> anyhow::Result<()> {
    enable_raw_mode()?;
    let mut stdout = std::io::stdout();
    crossterm::execute!(stdout, crossterm::terminal::EnterAlternateScreen)?;
    let backend = ratatui::backend::CrosstermBackend::new(stdout);
    let mut terminal = ratatui::Terminal::new(backend)?;
    let started = Instant::now();

    loop {
        terminal.draw(|f| {
            let size = f.size();
            let block = Block::default()
                .title("porkchop::screen — observed synthetic sequences")
                .borders(Borders::ALL);
            f.render_widget(block, size);

            let layout = ratatui::layout::Layout::default()
                .direction(ratatui::layout::Direction::Vertical)
                .constraints([
                    ratatui::layout::Constraint::Length(3),
                    ratatui::layout::Constraint::Min(3),
                    ratatui::layout::Constraint::Length(1),
                ])
                .margin(1)
                .split(size);

            // Header stats
            let hits_sum: usize = {
                let g = unit.lock().unwrap();
                g.values().sum()
            };
            let scr = screened.load(Ordering::Relaxed) as f64;
            let hits = hits_sum as f64;
            let uncls = unclassified.load(Ordering::Relaxed) as f64;
            let skip = skipped.load(Ordering::Relaxed) as f64;
            let tot_seen = scr + skip;
            let rwh = reads_with_hits.load(Ordering::Relaxed) as f64;
            let hp = if scr > 0.0 { 100.0 * hits / scr } else { 0.0 };
            let up = if scr > 0.0 { 100.0 * uncls / scr } else { 0.0 };
            let sp = if tot_seen > 0.0 { 100.0 * skip / tot_seen } else { 0.0 };
            let rwp = if scr > 0.0 { 100.0 * rwh / scr } else { 0.0 };

            let stats = format!(
                "screened: {}  total hits: {} ({:.1} hits/read)  reads with ≥1 hit: {} ({:.1}%)    unclassified: {} ({:.1}%)    skipped (not sampled): {} ({:.1}%)",
                scr as u64, hits_sum as u64, if scr > 0.0 { hits / scr } else { 0.0 },
                reads_with_hits.load(Ordering::Relaxed), rwp,
                unclassified.load(Ordering::Relaxed), up,
                skipped.load(Ordering::Relaxed), sp
            );
            let stats_para = ratatui::widgets::Paragraph::new(stats);
            f.render_widget(stats_para, layout[0]);

            let cols = ratatui::layout::Layout::default()
                .direction(ratatui::layout::Direction::Horizontal)
                .constraints([
                    ratatui::layout::Constraint::Percentage(50),
                    ratatui::layout::Constraint::Percentage(50),
                ])
                .split(layout[1]);

            // Top synthetic sequences
            let mut unit_items: Vec<(String, SeqKind, usize)> = {
                let g = unit.lock().unwrap();
                g.iter().map(|((name, kind), c)| (name.clone(), *kind, *c)).collect()
            };
            unit_items.sort_by(|a, b| b.2.cmp(&a.2));

            let mut unit_rows: Vec<Row> = Vec::new();
            for (name, kind, c) in unit_items.into_iter().take(12) {
                let f = { let g = fwd.lock().unwrap(); *g.get(&(name.clone(), kind)).unwrap_or(&0) };
                let r = { let g = rev.lock().unwrap(); *g.get(&(name.clone(), kind)).unwrap_or(&0) };
                unit_rows.push(Row::new(vec![
                    name,
                    match kind {
                        SeqKind::AdapterTop | SeqKind::AdapterBottom => "Adapter".to_string(),
                        SeqKind::Primer => "Primer".to_string(),
                        SeqKind::Barcode => "Barcode".to_string(),
                        SeqKind::Flank => "Flank".to_string(),
                    },
                    format!("{}", f),
                    format!("{}", r),
                    format!("{}", c),
                ]));
            }

            let unit_table = Table::new(
                unit_rows,
                [
                    ratatui::layout::Constraint::Percentage(22),
                    ratatui::layout::Constraint::Percentage(14),
                    ratatui::layout::Constraint::Percentage(14),
                    ratatui::layout::Constraint::Percentage(14),
                    ratatui::layout::Constraint::Percentage(10),
                ],
            )
            .header(Row::new(vec!["name", "kind", "(+)", "(-)", "reads"]).bold())
            .block(Block::default().borders(Borders::ALL).title("Top synthetic sequences"));
            f.render_widget(unit_table, cols[0]);

            // Top co-occurrence
            let mut combo_items: Vec<(String, usize)> = {
                let g = combos.lock().unwrap();
                g.iter().map(|(k, v)| (k.clone(), *v)).collect()
            };
            combo_items.sort_by(|a, b| b.1.cmp(&a.1));

            let mut combo_rows: Vec<Row> = Vec::new();
            for (id, c) in combo_items.into_iter().take(12) {
                combo_rows.push(Row::new(vec![id, format!("{}", c)]));
            }
            let combo_table = Table::new(
                combo_rows,
                [
                    ratatui::layout::Constraint::Percentage(86),
                    ratatui::layout::Constraint::Percentage(14),
                ],
            )
            .header(Row::new(vec!["aggregate identifier", "count"]).bold())
            .block(Block::default().borders(Borders::ALL).title("Top co-occurrence contexts"));
            f.render_widget(combo_table, cols[1]);
        
            // Footer: performance indicator
            let elapsed = started.elapsed().as_secs_f64();
            let rate = if elapsed > 0.0 { (screened.load(Ordering::Relaxed) as f64) / elapsed } else { 0.0 };
            let footer = format!("rate: {:.1} seq/s   elapsed: {:.1}s   screened: {}",
                                 rate, elapsed, screened.load(Ordering::Relaxed));
            let foot_para = ratatui::widgets::Paragraph::new(footer);
            f.render_widget(foot_para, layout[2]);
})?;

        // Keys
        if crossterm::event::poll(tick)? {
            if let Event::Key(k) = event::read()? {
                match k.code {
                    KeyCode::Char('q') | KeyCode::Esc => {
                        done.store(true, Ordering::SeqCst);
                        break;
                    }
                    _ => {}
                }
            }
        }

        if done.load(Ordering::SeqCst) {
            break;
        }
    }

    terminal.show_cursor().ok();
    disable_raw_mode().ok();
    let _ = crossterm::execute!(std::io::stdout(), crossterm::cursor::Show, crossterm::terminal::LeaveAlternateScreen);
    Ok(())
}
