
use std::collections::{BTreeSet, HashMap};
use std::path::{Path, PathBuf};
use std::sync::mpsc;
use std::time::{Duration, Instant};

use rayon::prelude::*;

#[derive(Clone)]
struct OwnedRecord {
    id: String,
    seq: Vec<u8>,  // uppercase A/C/G/T/N
    qual: Vec<u8>, // phred+33, or 'I' if missing
}

fn write_fastq_record<W: std::io::Write>(w: &mut W, id: &str, seq: &[u8], qual: &[u8]) -> std::io::Result<()> {
    w.write_all(b"@")?;
    w.write_all(id.as_bytes())?;
    w.write_all(b"\n")?;
    w.write_all(seq)?;
    w.write_all(b"\n+\n")?;
    w.write_all(qual)?;
    w.write_all(b"\n")?;
    Ok(())
}

// ---------- edlib wrapper ----------
mod edwrap {
    use edlib_rs::edlibrs::{edlibAlignRs, EdlibAlignConfigRs, EdlibAlignModeRs, EdlibAlignTaskRs, EdlibEqualityPairRs};
    pub struct Hit { pub start: i32, pub end: i32, pub edits: i32 }
    pub fn locate(pattern: &[u8], text: &[u8], max_edits: i32) -> Option<Hit> {
        // edlib_rs 0.1.2 expects a slice for `additionalequalities`
        let empty: &[EdlibEqualityPairRs] = &[];
        let cfg = EdlibAlignConfigRs { k: max_edits, mode: EdlibAlignModeRs::EDLIB_MODE_HW, task: EdlibAlignTaskRs::EDLIB_TASK_LOC, additionalequalities: empty };
        let res = edlibAlignRs(pattern, text, &cfg);
        if res.editDistance < 0 { return None; }
        let start = res.startLocations.as_ref()?.get(0).copied()?;
        let end   = res.endLocations.as_ref()?.get(0).copied()?;
        Some(Hit { start, end, edits: res.editDistance })
    }
}

#[derive(Clone)]
struct Motif<'a> { name: &'a str, kind: &'a str, seq: &'a [u8] }

fn motifs_for_kit<'a>(kit: &'static crate::kit::Kit) -> Vec<Motif<'a>> {
    let mut m = Vec::new();
    for s in kit.adapters_and_primers {
        m.push(Motif { name: s.name, kind: "adapter_or_primer", seq: s.sequence.as_bytes() });
    }
    for s in kit.barcodes {
        m.push(Motif { name: s.name, kind: "barcode_or_flank", seq: s.sequence.as_bytes() });
    }
    m
}

fn normalize_seq(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|&b| match b { b'a'..=b'z' => b.to_ascii_uppercase(), _ => b }).collect()
}
#[allow(dead_code)]
#[allow(dead_code)]
fn max_edits_for(len: usize) -> i32 { ((len as f64 * 0.15).ceil() as i32).max(1) }

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct ModalityKey { left: String, right: String, barcode: String }

#[derive(Clone)]
struct CleanResult { rec: OwnedRecord, modality: ModalityKey, clipped: bool, structure: String }

fn annotate_and_trim_one(seq: &[u8], qual: &[u8], _kit_id: &str, motifs: &[Motif], edits: i32) -> CleanResult {
    let s = normalize_seq(seq);
    let n = s.len() as i32;

    let mut left_best: Option<(i32, i32, i32, &str)> = None;
    let mut right_best: Option<(i32, i32, i32, &str)> = None;
    let mut barcode_left: Option<(i32, i32, i32, &str)> = None;
    let mut barcode_right: Option<(i32, i32, i32, &str)> = None;
let s = normalize_seq(seq);
    let n = s.len() as i32;

    let mut left_best: Option<(i32, i32, i32, &str)> = None;
    let mut right_best: Option<(i32, i32, i32, &str)> = None;
    let mut barcode: Option<String> = None;

    for m in motifs {
        if let Some(hit) = edwrap::locate(m.seq, &s, edits) {
            let center = (hit.start + hit.end) / 2;
            match m.kind {
                "adapter_or_primer" => {
                    if center < 300 {
                        if left_best.map_or(true, |lb| hit.edits < lb.2) {
                            left_best = Some((hit.start, hit.end, hit.edits, m.name));
                        }
                    }
                    if center > n - 300 {
                        if right_best.map_or(true, |rb| hit.edits < rb.2) {
                            right_best = Some((hit.start, hit.end, hit.edits, m.name));
                        }
                    }
                }
                "barcode_or_flank" => {
                    if center <= n / 2 {
                        if barcode_left.map_or(true, |b| hit.edits < b.2) {
                            barcode_left = Some((hit.start, hit.end, hit.edits, m.name));
                        }
                    } else {
                        if barcode_right.map_or(true, |b| hit.edits < b.2) {
                            barcode_right = Some((hit.start, hit.end, hit.edits, m.name));
                        }
                    }
                }
                _ => {}
            }
        }
    }

    let mut left_cut:  i32 = 0;
    let mut right_cut: i32 = n;
    let mut notes: Vec<String> = Vec::new();

    if let Some((st, en, ed, nm)) = left_best  { left_cut = en + 1; notes.push(format!("L:{}:{}-{}:ed={}", nm, st, en, ed)); }
    if let Some((st, en, ed, nm)) = right_best { right_cut = st;   notes.push(format!("R:{}:{}-{}:ed={}", nm, st, en, ed)); }

    // Also clip barcodes at ends if detected
    if let Some((st, en, _ed, nm)) = barcode_left { if en + 1 > left_cut { left_cut = en + 1; notes.push(format!("BL:{}:{}-{}", nm, st, en)); } }
    if let Some((st, en, _ed, nm)) = barcode_right { if st < right_cut { right_cut = st; notes.push(format!("BR:{}:{}-{}", nm, st, en)); } }


    if left_cut < 0 { left_cut = 0; }
    if right_cut > n { right_cut = n; }
    if left_cut >= right_cut { left_cut = 0; right_cut = n; } // unclippable: pass-through

    let start = left_cut as usize;
    let end   = right_cut as usize;
    let new_seq  = s[start..end].to_vec();
    let new_qual = if !qual.is_empty() { qual[start..end].to_vec() } else { vec![b'I'; new_seq.len()] };

    let id = format!("trim={}..{};len={};{}", left_cut, right_cut, n, notes.join(";"));
    let modality = ModalityKey {
        left:    left_best.map(|t| t.3.to_string()).unwrap_or_else(|| "—".into()),
        right:   right_best.map(|t| t.3.to_string()).unwrap_or_else(|| "—".into()),
        barcode: barcode.unwrap_or_else(|| "—".into()),
    };
    let clipped = left_best.is_some() || right_best.is_some();
    let mut structure: Vec<&str> = Vec::new();
    if left_best.is_some() { structure.push("sequencing adapter"); }
    if barcode_left.is_some() { structure.push("barcode"); }
    structure.push("insert");
    if barcode_right.is_some() { structure.push("reverse barcode"); }
    if right_best.is_some() { structure.push("reverse adapter"); }
    CleanResult { rec: OwnedRecord { id, seq: new_seq, qual: new_qual }, modality, clipped, structure: structure.join(" > ") }
}

struct Tallies { total: u64, clipped: u64, unclippable: u64, by_structure: HashMap<String, u64>, clip5_hist: HashMap<usize,u64>, clip3_hist: HashMap<usize,u64> }

impl Default for Tallies {
    fn default() -> Self {
        Tallies {
            total: 0,
            clipped: 0,
            unclippable: 0,
            by_structure: HashMap::new(),
            clip5_hist: HashMap::new(),
            clip3_hist: HashMap::new(),
        }
    }
}
enum StatEvent { Seen(String, bool), Clip(usize, usize), Done }

fn expected_modalities(kit: &'static crate::kit::Kit) -> BTreeSet<(String,String)> {
    let mut names: Vec<String> = Vec::new();
    for s in kit.adapters_and_primers { names.push(s.name.to_string()); }
    for s in kit.barcodes            { names.push(s.name.to_string()); }
    let mut set = BTreeSet::new();
    for l in &names { for r in &names { set.insert((l.clone(), r.clone())); } }
    set
}

fn draw_dashboard<B: ratatui::backend::Backend>(terminal: &mut ratatui::Terminal<B>, tallies: &Tallies) -> std::io::Result<()> {
    use ratatui::layout::{Constraint, Direction, Layout};
    use ratatui::text::Text;
    use ratatui::widgets::{Block, Borders, Paragraph, Row, Table, BarChart};
    use std::collections::HashMap;

    terminal.draw(|f| {
        let size = f.size();
        let chunks = Layout::default()
            .direction(Direction::Vertical)
            .constraints([Constraint::Min(8), Constraint::Length(3), Constraint::Length(2), Constraint::Min(12)].as_ref())
            .split(size);

        // Observed contexts (top)
        let mut rows: Vec<(String, u64)> = tallies.by_structure.iter().map(|(k,v)| (k.clone(), *v)).collect();
        rows.sort_by(|a,b| b.1.cmp(&a.1));
        rows.truncate(20);
        let table_rows = rows.into_iter().map(|(k,c)| Row::new(vec![k, c.to_string()]));
        let table = Table::new(
            table_rows,
            [Constraint::Percentage(80), Constraint::Length(10)],
        )
            .header(Row::new(vec!["Structure", "Count"]))
            .block(Block::default().borders(Borders::ALL).title("Observed modalities (top 20)"));
        f.render_widget(table, chunks[0]);

        // Helper to build dynamic bins for a given histogram map
        fn build_binned<'a>(hm: &HashMap<usize,u64>, max_bars: usize) -> (Vec<(&'a str, u64)>, Vec<(String,u64)>, usize, usize, usize) {
            if hm.is_empty() {
                return (vec![("0", 0)], vec![("0".to_string(), 0)], 0, 0, 1);
            }
            let min_k = *hm.keys().min().unwrap();
            let max_k = *hm.keys().max().unwrap();
            let span = max_k - min_k + 1;
            let max_bars = if max_bars == 0 { 1 } else { max_bars };
            let bin_size = std::cmp::max(1usize, (span + max_bars - 1) / max_bars);
            let bin_count = (span + bin_size - 1) / bin_size;
            let mut bins: Vec<u64> = vec![0; bin_count];
            for (k, v) in hm.iter() {
                let idx = ((*k - min_k) / bin_size).min(bin_count - 1);
                bins[idx] += *v;
            }
            let mut data: Vec<(&'a str, u64)> = Vec::with_capacity(bin_count);
            let mut summary: Vec<(String,u64)> = Vec::with_capacity(bin_count);
            for i in 0..bin_count {
                let start = min_k + i*bin_size;
                let end = std::cmp::min(start + bin_size - 1, max_k);
                let lbl = if bin_size == 1 { format!("{}", start) } else { format!("{}-{}", start, end) };
                let leaked: &'a str = Box::leak(lbl.clone().into_boxed_str());
                let count = bins[i];
                data.push((leaked, count));
                summary.push((lbl, count));
            }
            (data, summary, min_k, max_k, bin_size)
        }

        // Summary line
        let summary = Paragraph::new(Text::from(format!(
            "total: {}   clipped: {}   unclippable: {}   modalities: {}",
            tallies.total, tallies.clipped, tallies.unclippable, tallies.by_structure.len()
        ))).block(Block::default().borders(Borders::ALL).title("Summary"));
        f.render_widget(summary, chunks[1]);

        // Bottom area split vertically into a small bin summary row and the charts row
        let bottom = Layout::default()
            .direction(Direction::Vertical)
            .constraints([Constraint::Length(6), Constraint::Min(6)].as_ref())
            .split(chunks[3]);

        // Compute dynamic bins using available width
        let chart_width = std::cmp::max(10usize, bottom[1].width as usize / 2); // approximate half width per chart
        let (left_data, left_bins, left_min, left_max, left_step) = build_binned(&tallies.clip5_hist, chart_width.saturating_sub(4));
        let (right_data, right_bins, right_min, right_max, right_step) = build_binned(&tallies.clip3_hist, chart_width.saturating_sub(4));

        // Legend with min/max and bin size for both ends
        let legend = Paragraph::new(Text::from(format!(
            "x: clipped nt | y: read count   |   5′: min={} max={} bin={}   |   3′: min={} max={} bin={}",
            left_min, left_max, left_step, right_min, right_max, right_step
        ))).block(Block::default().borders(Borders::ALL).title("Legend"));
        f.render_widget(legend, chunks[2]);

        // Bin summaries (top row of bottom area): show top bins with counts for each side
        fn top_rows<'a>(pairs: &[(String,u64)], n: usize) -> Vec<Row<'a>> {
            let mut v = pairs.to_vec();
            v.sort_by(|a,b| b.1.cmp(&a.1));
            v.truncate(n);
            v.into_iter().map(|(k,c)| Row::new(vec![k, c.to_string()])).collect()
        }

        let bin_tables = Layout::default().direction(Direction::Horizontal)
            .constraints([Constraint::Percentage(50), Constraint::Percentage(50)].as_ref())
            .split(bottom[0]);

        let left_table = Table::new(
                top_rows(&left_bins, 8),
                [Constraint::Percentage(70), Constraint::Length(10)],
            )
            .header(Row::new(vec!["5′ bin", "count"]))
            .block(Block::default().borders(Borders::ALL).title("5′ bin counts (top 8)"));
        let right_table = Table::new(
                top_rows(&right_bins, 8),
                [Constraint::Percentage(70), Constraint::Length(10)],
            )
            .header(Row::new(vec!["3′ bin", "count"]))
            .block(Block::default().borders(Borders::ALL).title("3′ bin counts (top 8)"));

        f.render_widget(left_table, bin_tables[0]);
        f.render_widget(right_table, bin_tables[1]);

        // Charts (bottom row of bottom area)
        let charts_row = Layout::default()
            .direction(Direction::Horizontal)
            .constraints([Constraint::Percentage(50), Constraint::Percentage(50)].as_ref())
            .split(bottom[1]);
        let max_all = std::cmp::max(
            left_data.iter().map(|(_,v)| *v).max().unwrap_or(0),
            right_data.iter().map(|(_,v)| *v).max().unwrap_or(0)
        );
        let left_chart = BarChart::default()
            .block(Block::default().borders(Borders::ALL).title("5′ clipped (nt) — count"))
            .data(&left_data)
            .bar_width(1)
            .bar_gap(0)
            .value_style(ratatui::style::Style::default().add_modifier(ratatui::style::Modifier::BOLD))
            .label_style(ratatui::style::Style::default())
            .max(max_all);
        let right_chart = BarChart::default()
            .block(Block::default().borders(Borders::ALL).title("3′ clipped (nt) — count"))
            .data(&right_data)
            .bar_width(1)
            .bar_gap(0)
            .value_style(ratatui::style::Style::default().add_modifier(ratatui::style::Modifier::BOLD))
            .label_style(ratatui::style::Style::default())
            .max(max_all);
        f.render_widget(left_chart, charts_row[0]);
        f.render_widget(right_chart, charts_row[1]);
    })?;
    Ok(())
}
fn stats_thread(rx: mpsc::Receiver<StatEvent>, _kit: &'static crate::kit::Kit) -> std::thread::JoinHandle<()> {
    use crossterm::{execute, terminal::{disable_raw_mode, enable_raw_mode, EnterAlternateScreen, LeaveAlternateScreen}};
    use ratatui::backend::CrosstermBackend;
    use std::io::stdout;

    std::thread::spawn(move || {
        let mut tallies = Tallies::default();
        let mut out = stdout();
        let _ = enable_raw_mode();
        let _ = execute!(out, EnterAlternateScreen);
        let backend = CrosstermBackend::new(out);
        let mut term = ratatui::Terminal::new(backend).expect("tui terminal");

        let tick = Duration::from_millis(200);
        let mut last = Instant::now();
        let mut done = false;

        while !done {
            while let Ok(ev) = rx.try_recv() {
                match ev {
                    StatEvent::Seen(modality, clipped) => {
                        tallies.total += 1;
                        if clipped { tallies.clipped += 1; } else { tallies.unclippable += 1; }
                        *tallies.by_structure.entry(modality).or_insert(0) += 1;
                    }
                    StatEvent::Clip(l5, l3) => { *tallies.clip5_hist.entry(l5).or_insert(0) += 1; *tallies.clip3_hist.entry(l3).or_insert(0) += 1; },
                    StatEvent::Done => { done = true; }
                }
            }
            if last.elapsed() >= tick {
                let _ = draw_dashboard(&mut term, &tallies);
                last = Instant::now();
            }
            std::thread::sleep(Duration::from_millis(25));
        }

        let _ = disable_raw_mode();
        let _ = term.show_cursor();
        let mut out2 = std::io::stdout();
        let _ = execute!(out2, LeaveAlternateScreen);
    })
}



fn parse_trim_from_id(id: &str) -> (usize, usize) {
    // expects id like "trim=a..b;len=n;..."
    let mut left = 0usize;
    let mut right = 0usize;
    let mut total = None::<usize>;

    if let Some(rest) = id.strip_prefix("trim=") {
        let parts: Vec<&str> = rest.split(';').collect();
        if let Some(range) = parts.get(0) {
            if let Some((a, b)) = range.split_once("..") {
                if let (Ok(l), Ok(r)) = (a.parse::<i32>(), b.parse::<i32>()) {
                    left = if l > 0 { l as usize } else { 0 };
                    right = if r >= 0 { r as usize } else { 0 };
                }
            }
        }
        for part in parts.iter().skip(1) {
            if let Some(val) = part.strip_prefix("len=") {
                if let Ok(n) = val.parse::<usize>() {
                    total = Some(n);
                    break;
                }
            }
        }
    }
    let total = total.unwrap_or_else(|| right.max(left));
    let clip5 = left;
    let clip3 = if right <= total { total - right } else { 0 };
    (clip5, clip3)
}
fn ensure_known_kit(kit: &str) -> anyhow::Result<()> {
    if crate::get_sequences_for_kit(kit).is_none() {
        anyhow::bail!("Unknown kit: {}. Use `porkchop list-kits --format table` to see valid kit ids.", kit);
    }
    Ok(())
}

fn split_supported_files(paths: Vec<PathBuf>) -> (Vec<PathBuf>, Vec<PathBuf>) {
    let mut ok = Vec::new();
    let mut bad = Vec::new();
    for p in paths {
        let name = p.file_name().and_then(|s| s.to_str()).unwrap_or("").to_ascii_lowercase();
        let ext  = p.extension().and_then(|s| s.to_str()).unwrap_or("").to_ascii_lowercase();
        let is_ok = name.ends_with(".fastq.gz") || name.ends_with(".fq.gz") || ext == "fastq" || ext == "fq" || ext == "sam" || ext == "bam";
        if is_ok { ok.push(p); } else { bad.push(p); }
    }
    (ok, bad)
}

fn process_fastx_to_gz(out_path: &Path, input_files: Vec<PathBuf>, kit_id: &str, edits: i32, kit_ref: &'static crate::kit::Kit, events: &mpsc::Sender<StatEvent>) -> anyhow::Result<()> {
    use std::fs::File;
    use std::io::BufWriter;
    use needletail::parser::parse_fastx_file;

    let motifs = motifs_for_kit(kit_ref);

    let ofh = File::create(out_path)?;
    let writer = BufWriter::new(ofh);
    let mut gz = flate2::write::GzEncoder::new(writer, flate2::Compression::default());

    const CHUNK: usize = 2000;

    for path in input_files {
        let lower = path.to_string_lossy().to_ascii_lowercase();

        if lower.ends_with(".sam") {
            use rust_htslib::bam::{self, Read};
            let mut reader = bam::Reader::from_path(&path)?;
            let mut buf: Vec<rust_htslib::bam::Record> = Vec::new();

            for r in reader.records() {
                if let Ok(rec) = r { buf.push(rec); }
                if buf.len() >= CHUNK {
                    let processed: Vec<(String, CleanResult)> = buf.par_iter().map(|r| {
                        let name = std::str::from_utf8(r.qname()).unwrap_or("SAM");
                        let seq = r.seq().as_bytes();
                        let qual = r.qual().iter().map(|q| (q + 33) as u8).collect::<Vec<u8>>();
                        let cr = annotate_and_trim_one(&seq, &qual, kit_id, &motifs, edits);
                        (name.to_string(), cr)
                    }).collect();
                    for (name, cr) in &processed { let _ = events.send(StatEvent::Seen(cr.structure.clone(), cr.clipped)); let (lc, rc) = parse_trim_from_id(&cr.rec.id); let _ = events.send(StatEvent::Clip(lc, rc)); let (lc, rc) = parse_trim_from_id(&cr.rec.id); let _ = events.send(StatEvent::Clip(lc, rc)); let rid = format!("{} {}", name, cr.rec.id); write_fastq_record(&mut gz, &rid, &cr.rec.seq, &cr.rec.qual)?; }
                    buf.clear();
                }
            }
            if !buf.is_empty() {
                let processed: Vec<(String, CleanResult)> = buf.par_iter().map(|r| {
                        let name = std::str::from_utf8(r.qname()).unwrap_or("SAM");
                        let seq = r.seq().as_bytes();
                        let qual = r.qual().iter().map(|q| (q + 33) as u8).collect::<Vec<u8>>();
                        let cr = annotate_and_trim_one(&seq, &qual, kit_id, &motifs, edits);
                        (name.to_string(), cr)
                    }).collect();
                for (name, cr) in &processed { let _ = events.send(StatEvent::Seen(cr.structure.clone(), cr.clipped)); let (lc, rc) = parse_trim_from_id(&cr.rec.id); let _ = events.send(StatEvent::Clip(lc, rc)); let (lc, rc) = parse_trim_from_id(&cr.rec.id); let _ = events.send(StatEvent::Clip(lc, rc)); let rid = format!("{} {}", name, cr.rec.id); write_fastq_record(&mut gz, &rid, &cr.rec.seq, &cr.rec.qual)?; }
            }

        } else if lower.ends_with(".bam") {
            use rust_htslib::bam::{self, Read};
            let mut reader = bam::Reader::from_path(&path)?;
            let mut buf: Vec<rust_htslib::bam::Record> = Vec::new();

            for r in reader.records() {
                if let Ok(rec) = r { buf.push(rec); }
                if buf.len() >= CHUNK {
                    let processed: Vec<(String, CleanResult)> = buf.par_iter().map(|r| {
                        let name = std::str::from_utf8(r.qname()).unwrap_or("BAM");
                        let seq = r.seq().as_bytes();
                        let qual = r.qual().iter().map(|q| (q + 33) as u8).collect::<Vec<u8>>();
                        let cr = annotate_and_trim_one(&seq, &qual, kit_id, &motifs, edits);
                        (name.to_string(), cr)
                    }).collect();
                    for (name, cr) in &processed { let _ = events.send(StatEvent::Seen(cr.structure.clone(), cr.clipped)); let (lc, rc) = parse_trim_from_id(&cr.rec.id); let _ = events.send(StatEvent::Clip(lc, rc)); let (lc, rc) = parse_trim_from_id(&cr.rec.id); let _ = events.send(StatEvent::Clip(lc, rc)); let rid = format!("{} {}", name, cr.rec.id); write_fastq_record(&mut gz, &rid, &cr.rec.seq, &cr.rec.qual)?; }
                    buf.clear();
                }
            }
            if !buf.is_empty() {
                let processed: Vec<(String, CleanResult)> = buf.par_iter().map(|r| {
                        let name = std::str::from_utf8(r.qname()).unwrap_or("BAM");
                        let seq = r.seq().as_bytes();
                        let qual = r.qual().iter().map(|q| (q + 33) as u8).collect::<Vec<u8>>();
                        let cr = annotate_and_trim_one(&seq, &qual, kit_id, &motifs, edits);
                        (name.to_string(), cr)
                    }).collect();
                for (name, cr) in &processed { let _ = events.send(StatEvent::Seen(cr.structure.clone(), cr.clipped)); let (lc, rc) = parse_trim_from_id(&cr.rec.id); let _ = events.send(StatEvent::Clip(lc, rc)); let (lc, rc) = parse_trim_from_id(&cr.rec.id); let _ = events.send(StatEvent::Clip(lc, rc)); let rid = format!("{} {}", name, cr.rec.id); write_fastq_record(&mut gz, &rid, &cr.rec.seq, &cr.rec.qual)?; }
            }

        } else {
            let mut reader = parse_fastx_file(&path)?;
            loop {
                let mut owned_chunk: Vec<OwnedRecord> = Vec::with_capacity(CHUNK);
                for _ in 0..CHUNK {
                    match reader.next() {
                        Some(Ok(rec)) => {
                            let id   = String::from_utf8_lossy(rec.id()).to_string();
                            let seq  = rec.seq().to_vec();
                            let qual = rec.qual().map(|q| q.to_vec()).unwrap_or_else(|| vec![b'I'; seq.len()]);
                            owned_chunk.push(OwnedRecord { id, seq, qual });
                        }
                        Some(Err(_)) => continue,
                        None => break,
                    }
                }
                if owned_chunk.is_empty() { break; }
                let processed: Vec<CleanResult> = owned_chunk.par_iter()
                    .map(|r| annotate_and_trim_one(&r.seq, &r.qual, kit_id, &motifs, edits))
                    .collect();
                for (src, cr) in owned_chunk.iter().zip(processed.iter()) {
                    let _ = events.send(StatEvent::Seen(cr.structure.clone(), cr.clipped)); let (lc, rc) = parse_trim_from_id(&cr.rec.id); let _ = events.send(StatEvent::Clip(lc, rc)); let (lc, rc) = parse_trim_from_id(&cr.rec.id); let _ = events.send(StatEvent::Clip(lc, rc));
                    let rid = format!("{} {}", src.id, cr.rec.id);
                    write_fastq_record(&mut gz, &rid, &cr.rec.seq, &cr.rec.qual)?;
                }
            }
        }
    }

    gz.finish()?;
    Ok(())
}

pub fn run(threads: usize, kit: &str, edits: i32, output: &Path, files: Vec<PathBuf>) -> anyhow::Result<()> {
        ensure_known_kit(kit)?;
if crate::get_sequences_for_kit(kit).is_none() {
        anyhow::bail!("Unknown kit: {}. Use `porkchop list-kits --format table` to see valid kit ids.", kit);
    }
    let (ok, bad) = split_supported_files(files);
    if !bad.is_empty() {
        let mut msg = String::from("Unsupported file type(s):\n");
        for p in &bad { msg.push_str(&format!("  - {}\n", p.display())); }
        msg.push_str("Allowed: SAM (.sam), BAM (.bam), FASTQ (.fastq/.fq), and gzipped FASTQ (.fastq.gz/.fq.gz).");
        anyhow::bail!(msg);
    }

    let threads_eff = if threads == 0 { std::cmp::max(1, num_cpus::get()) } else { threads };
    rayon::ThreadPoolBuilder::new().num_threads(threads_eff).build_global().ok();

    let kit_ref: &'static crate::kit::Kit = crate::get_sequences_for_kit(kit).expect("validated kit");
    let (tx, rx) = mpsc::channel::<StatEvent>();
    let ui_handle = stats_thread(rx, kit_ref);

    eprintln!("clean: kit={} | threads={} | inputs={} | output={}", kit, threads_eff, ok.len(), output.display());
    let ret = process_fastx_to_gz(output, ok, kit, edits, kit_ref, &tx);

    let _ = tx.send(StatEvent::Done);
    let _ = ui_handle.join();

    ret
}