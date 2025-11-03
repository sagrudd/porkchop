use std::collections::{HashMap, HashSet};
use std::sync::{Arc, Mutex, atomic::{AtomicBool, AtomicUsize, Ordering}};
use std::time::{Duration, Instant};

use crate::benchmark::{self, BenchmarkAlgo};
use crate::kit::SeqKind;
use crate::list_supported_kits;
use crate::seqio::for_each_parallel;

// TUI
use crossterm::event::{self, Event, KeyCode};
use crossterm::terminal::{disable_raw_mode, enable_raw_mode};
use ratatui::prelude::*;
use ratatui::widgets::{Block, Borders, Row, Table};

#[derive(Debug, Clone)]
pub struct ScreenOpts {
    pub files: Vec<String>,
    pub threads: Option<usize>,
    pub fraction: f64,
    pub tick_secs: u64,
    pub algo: BenchmarkAlgo,
    pub max_dist: usize,
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

pub fn run_screen(opts: ScreenOpts) -> anyhow::Result<()> {
    let records = Arc::new(collect_all_sequences());

    // Tallies
    let unit_tally: Arc<Mutex<HashMap<(String, SeqKind), usize>>> = Arc::new(Mutex::new(HashMap::new()));
    let combo_tally: Arc<Mutex<HashMap<String, usize>>> = Arc::new(Mutex::new(HashMap::new()));
    let done = Arc::new(AtomicBool::new(false));
    let screened = Arc::new(AtomicUsize::new(0));
    let unclassified = Arc::new(AtomicUsize::new(0));
    let skipped = Arc::new(AtomicUsize::new(0));

    // Optional prebuilt for ACMyers
    let prebuilt = if let BenchmarkAlgo::ACMyers = opts.algo {
        Some(Arc::new(benchmark::prebuild_for(records.as_slice())))
    } else {
        None
    };

    // UI thread
    let unit_ui = unit_tally.clone();
    let combo_ui = combo_tally.clone();
    let done_ui = done.clone();
    let screened_ui = screened.clone();
    let unclassified_ui = unclassified.clone();
    let skipped_ui = skipped.clone();
    let tick = Duration::from_secs(opts.tick_secs.max(1));
    std::thread::spawn(move || {
        let _ = tui_loop(unit_ui, combo_ui, done_ui, screened_ui, unclassified_ui, skipped_ui, tick);
    });

    // Sampling params
    let p = opts.fraction.clamp(0.0, 1.0);
    let threads = opts.threads;

    for file in &opts.files {
        let file = file.clone();
        let unit_w = unit_tally.clone();
        let combo_w = combo_tally.clone();
        let algo = opts.algo;
        let records_arc = records.clone();
        let prebuilt_c = prebuilt.clone();
        let md = opts.max_dist;
        let p_sample = p;
        let done_c = done.clone();
        let screened_c = screened.clone();
        let unclassified_c = unclassified.clone();
        let skipped_c = skipped.clone();

        let _ = for_each_parallel(file, threads, move |read| {
            if done_c.load(Ordering::SeqCst) { return; }

            // Bernoulli(p) sampling via deterministic hash of read id
            let take = if p_sample >= 1.0 {
                true
            } else {
                let mut h: u64 = 0xcbf29ce484222325;
                for b in read.id.as_bytes() { h = h.wrapping_mul(0x100000001b3) ^ (*b as u64); }
                (h as f64 / std::u64::MAX as f64) < p_sample
            };
            if !take { skipped_c.fetch_add(1, Ordering::Relaxed); return; }

            // Enumerate all motif hits for this read using requested algorithm
            let hits = benchmark::classify_all(algo, &read.seq, records_arc.as_slice(), prebuilt_c.as_deref(), md);

            if hits.is_empty() {
                screened_c.fetch_add(1, Ordering::Relaxed);
                unclassified_c.fetch_add(1, Ordering::Relaxed);
                return;
            }

            // Tally individual hits
            {
                let mut g = unit_w.lock().unwrap();
                for (name, kind, _is_rc) in &hits {
                    *g.entry((name.clone(), *kind)).or_insert(0) += 1;
                }
            }

            // Compose aggregate identifier for this read: e.g., "NB_flank_fwd + 16S_rev_target"
            let mut labels: Vec<String> = Vec::new();
            let mut seen = HashSet::new();
            for (name, kind, is_rc) in hits {
                let orient = if is_rc { "rev" } else { "fwd" };
                let leaf = format!("{}_{}_{}", name, orient, kind_suffix(kind));
                if seen.insert(leaf.clone()) {
                    labels.push(leaf);
                }
            }
            labels.sort();
            let combo = labels.join(" + ");

            {
                let mut g = combo_w.lock().unwrap();
                *g.entry(combo).or_insert(0) += 1;
            }

            screened_c.fetch_add(1, Ordering::Relaxed);
        });
    }

    done.store(true, Ordering::SeqCst);
    std::thread::sleep(Duration::from_millis(150));
    Ok(())
}

fn tui_loop(
    unit: Arc<Mutex<HashMap<(String, SeqKind), usize>>>,
    combos: Arc<Mutex<HashMap<String, usize>>>,
    done: Arc<AtomicBool>,
    screened: Arc<AtomicUsize>,
    unclassified: Arc<AtomicUsize>,
    skipped: Arc<AtomicUsize>,
    tick: Duration
) -> anyhow::Result<()> {
    enable_raw_mode()?;
    let mut stdout = std::io::stdout();
    crossterm::execute!(stdout, crossterm::terminal::EnterAlternateScreen)?;
    let backend = ratatui::backend::CrosstermBackend::new(stdout);
    let mut terminal = ratatui::Terminal::new(backend)?;

    let start = Instant::now();

    loop {
        terminal.draw(|f| {
            let size = f.size();
            let block = Block::default().title("porkchop::screen — observed synthetic sequences").borders(Borders::ALL);
            f.render_widget(block, size);

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
            let hp = if scr > 0.0 { 100.0 * hits / scr } else { 0.0 };
            let up = if scr > 0.0 { 100.0 * uncls / scr } else { 0.0 };
            let sp = if tot_seen > 0.0 { 100.0 * skip / tot_seen } else { 0.0 };
            let stats = format!(
                "screened: {}  hits: {} ({:.1}%)  unclassified: {} ({:.1}%)    skipped (not sampled): {} ({:.1}%)",
                scr as u64, hits_sum as u64, hp, unclassified.load(Ordering::Relaxed), up, skipped.load(Ordering::Relaxed), sp
            );
            let stats_para = ratatui::widgets::Paragraph::new(stats);
            let stats_area = Rect::new(size.x + 2, size.y + 1, size.width.saturating_sub(4), 1);
            f.render_widget(stats_para, stats_area);

            // Top individual sequences (name, kind, count)
            let mut unit_rows: Vec<Row> = Vec::new();
            let mut unit_items: Vec<(String, SeqKind, usize)> = {
                let g = unit.lock().unwrap();
                g.iter().map(|((name, kind), c)| (name.clone(), *kind, *c)).collect()
            };
            unit_items.sort_by(|a, b| b.2.cmp(&a.2));
            for (name, kind, c) in unit_items.into_iter().take(12) {
                unit_rows.push(Row::new(vec![
                    name,
                    match kind {
                        SeqKind::AdapterTop | SeqKind::AdapterBottom => "Adapter".to_string(),
                        SeqKind::Primer => "Primer".to_string(),
                        SeqKind::Barcode => "Barcode".to_string(),
                        SeqKind::Flank => "Flank".to_string(),
                    },
                    format!("{}", c),
                ]));
            }
            let unit_table = Table::new(
                unit_rows,
                [Constraint::Percentage(50), Constraint::Percentage(25), Constraint::Percentage(25)]
            )
                .header(Row::new(vec!["name", "kind", "count"]).bold())
                .block(Block::default().borders(Borders::ALL).title("Top synthetic sequences"));

            // Top combos table (aggregate identifiers)
            let mut combo_rows: Vec<Row> = Vec::new();
            let mut combo_items: Vec<(String, usize)> = {
                let g = combos.lock().unwrap();
                g.iter().map(|(k, v)| (k.clone(), *v)).collect()
            };
            combo_items.sort_by(|a, b| b.1.cmp(&a.1));
            for (id, c) in combo_items.into_iter().take(12) {
                combo_rows.push(Row::new(vec![id, format!("{}", c)]));
            }
            let combo_table = Table::new(
                combo_rows,
                [Constraint::Percentage(75), Constraint::Percentage(25)]
            )
                .header(Row::new(vec!["aggregate identifier", "count"]).bold())
                .block(Block::default().borders(Borders::ALL).title("Top co-occurrence contexts"));

            // Layout: header line + two stacked tables + footer
            let unit_area = Rect::new(size.x + 2, size.y + 3, size.width.saturating_sub(4), (size.height.saturating_sub(7)) / 2);
            let combo_area = Rect::new(size.x + 2, unit_area.y + unit_area.height, size.width.saturating_sub(4), (size.height.saturating_sub(7)) - unit_area.height);
            f.render_widget(unit_table, unit_area);
            f.render_widget(combo_table, combo_area);

            // Footer rate
            let elapsed = start.elapsed().as_secs_f64().max(1e-6);
            let rate = (screened.load(Ordering::Relaxed) as f64) / elapsed;
            let footer = ratatui::widgets::Paragraph::new(format!("rate: {:.1} screened/s   (q/Esc to quit)", rate));
            let farea = Rect::new(size.x + 2, size.y + size.height.saturating_sub(2), size.width.saturating_sub(4), 1);
            f.render_widget(footer, farea);
        })?;

        // Quit on 'q' or Esc
        if event::poll(tick)? {
            if let Event::Key(k) = event::read()? {
                if k.code == KeyCode::Char('q') || k.code == KeyCode::Esc {
                    done.store(true, Ordering::SeqCst);
                    disable_raw_mode()?;
                    crossterm::execute!(std::io::stdout(), crossterm::terminal::LeaveAlternateScreen)?;
                    std::process::exit(0);
                }
            }
        }

        if done.load(Ordering::SeqCst) {
            break;
        }
    }

    disable_raw_mode()?;
    crossterm::execute!(std::io::stdout(), crossterm::terminal::LeaveAlternateScreen)?;
    Ok(())
}
