use std::collections::HashMap;
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


fn revcomp(seq: &str) -> String {
    fn comp(b: u8) -> u8 {
        match b {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            b'U' | b'u' => b'A',
            b'N' | b'n' => b'N',
            _ => b'N',
        }
    }
    let mut out = Vec::with_capacity(seq.len());
    for &b in seq.as_bytes().iter().rev() {
        out.push(comp(b));
    }
    String::from_utf8(out).unwrap_or_default()
}

fn collect_all_sequences() -> Vec<crate::kit::SequenceRecord> {
    let mut v = Vec::new();
    for k in list_supported_kits() {
        v.extend_from_slice(k.adapters_and_primers);
        v.extend_from_slice(k.barcodes);
    }
    v
}

pub fn run_screen(opts: ScreenOpts) -> anyhow::Result<()> {
    let records = Arc::new(collect_all_sequences());

    let tally: Arc<Mutex<HashMap<(String, SeqKind), usize>>> = Arc::new(Mutex::new(HashMap::new()));
    let done = Arc::new(AtomicBool::new(false));
    let screened = Arc::new(AtomicUsize::new(0));
    let unclassified = Arc::new(AtomicUsize::new(0));
    let skipped = Arc::new(AtomicUsize::new(0));
// UI thread
    let tally_ui = tally.clone();
    let done_ui = done.clone();
    let screened_ui = screened.clone();
    let unclassified_ui = unclassified.clone();
    let skipped_ui = skipped.clone();
let tick = Duration::from_secs(opts.tick_secs.max(1));
    std::thread::spawn(move || {
        let _ = tui_loop(tally_ui, done_ui, screened_ui, unclassified_ui, skipped_ui, tick);
    });

    // Sampling params
    let p = opts.fraction.clamp(0.0, 1.0);
    let threads = opts.threads;

    // Process each file with seqio::for_each_parallel
    for file in &opts.files {
        let file = file.clone();
        let tally_w = tally.clone();
        let algo = opts.algo;
        let records_arc = records.clone();
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

            // Count as 'screened' only when sampled and assessed
            screened_c.fetch_add(1, Ordering::Relaxed);

            if let Some(hit) = benchmark::classify_best(algo, &read.seq, records_arc.as_slice(), md) {
            let mut g = tally_w.lock().unwrap();
            *g.entry((hit.name, hit.kind)).or_insert(0) += 1;
        } else {
            unclassified_c.fetch_add(1, Ordering::Relaxed);
            }
        });
    }

    done.store(true, Ordering::SeqCst);
    std::thread::sleep(Duration::from_millis(150));
    Ok(())
}

fn tui_loop(
    tally: Arc<Mutex<HashMap<(String, SeqKind), usize>>>,
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
        // Draw dashboard
        terminal.draw(|f| {
            let size = f.size();
            let block = Block::default().title("porkchop::screen — observed synthetic sequences").borders(Borders::ALL);
            f.render_widget(block, size);

            // Stats line
            let hits_sum: usize = {
                let g = tally.lock().unwrap();
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

            // Top-k table
            let mut rows: Vec<Row> = Vec::new();
            let mut items: Vec<(String, SeqKind, usize)> = {
                let g = tally.lock().unwrap();
                g.iter().map(|((name, kind), c)| (name.clone(), *kind, *c)).collect()
            };
            items.sort_by(|a, b| b.2.cmp(&a.2));
            for (name, kind, c) in items.into_iter().take(20) {
                rows.push(Row::new(vec![
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

            let table = Table::new(
                rows,
                [Constraint::Percentage(50), Constraint::Percentage(25), Constraint::Percentage(25)]
            )
                .header(Row::new(vec!["name", "kind", "count"]).bold())
                .block(Block::default().borders(Borders::ALL).title("Top synthetic sequences"));

            let area = Rect::new(size.x + 2, size.y + 3, size.width.saturating_sub(4), size.height.saturating_sub(5));
            f.render_widget(table, area);

        // Footer: rate per second (screened)
        let elapsed = start.elapsed().as_secs_f64().max(1e-6);
        let rate = (screened.load(Ordering::Relaxed) as f64) / elapsed;
        let footer = ratatui::widgets::Paragraph::new(format!("rate: {:.1} screened/s", rate));
        let farea = Rect::new(size.x + 2, size.y + size.height.saturating_sub(2), size.width.saturating_sub(4), 1);
        f.render_widget(footer, farea);
        })?;

        // Quit on 'q' or Esc: set cancel flag, restore terminal, exit process
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
