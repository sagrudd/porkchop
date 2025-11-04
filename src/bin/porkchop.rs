//! Module `porkchop` — auto‑generated docs for porkchop.
use clap::{Parser, Subcommand, ValueEnum};
use polars::prelude::*;

/// Porkchop CLI
#[derive(Parser)]
#[command(name = "porkchop")]
#[command(version)]
#[command(about = "ONP kit registry, IO and benchmarking", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
enum OutputFormat { Csv, Md, Table }
#[derive(Subcommand)]
enum Commands {
    /// List all supported kits
    ListKits {
        /// Output format: csv | md | table
        #[arg(long, value_enum, default_value_t = OutputFormat::Table)]
        format: OutputFormat,
        /// Show full table width (no truncation) [only for --format table]
        #[arg(long)]
        full: bool,
        /// Truncate string cells to N characters (ignored if --format != table or if --full)
        #[arg(long)]
        truncate: Option<usize>,
    },

    /// Describe a kit by id (e.g., "LSK114")
    Describe {
        /// Kit id to describe
        id: String,
    },

    /// Benchmark classification algorithms against an optional truth set
    Benchmark {
        /// Input files (FASTQ/FASTA/FASTQ.GZ/SAM/BAM)
        #[arg(required = true)]
        files: Vec<String>,
        /// Kit id (e.g., "LSK114")
        kit: String,
        /// Truth CSV (optional)
        #[arg(long)]
        truth: Option<String>,
        /// Comma-separated algorithms (edlib,myers,acmyers). Default handled in code.
        #[arg(long, default_value = "edlib,myers,acmyers")]
        algorithms: String,
        /// Edit-distance threshold
        #[arg(long, default_value_t = 24)]
        max_dist: usize,
        /// Threads (0/None = all)
        #[arg(long)]
        threads: Option<usize>,
        /// Emit CSV to stdout
        #[arg(long)]
        csv: bool,
    },

    /// Screen a dataset to infer library chemistry by scoring adapters/primers/barcodes
    Screen {
        /// Input files (FASTQ/FASTQ.GZ/SAM/BAM)
        #[arg(required = true)]
        files: Vec<String>,
        /// Algorithm (default: edlib)
        #[arg(long, default_value = "edlib")]
        algorithm: String,
        /// Edit-distance threshold (default: 24)
        #[arg(long, default_value_t = 24)]
        #[arg(long, default_value_t = 2)]
        max_dist: usize,
        /// Fraction of reads to sample (0.0-1.0; default 0.05)
        #[arg(long, default_value_t = 0.05)]
        fraction: f64,
        /// UI refresh in seconds (default: 2)
        #[arg(long, default_value_t = 2)]
        tick: u64,
        /// Threads (0/None = all)
        #[arg(long)]
        threads: Option<usize>,
        /// Write aggregate identifiers and counts to JSON file
        #[arg(long)]
        json: Option<String>,
           /// Minimum probability to display a kit (default: 0.1)
        #[arg(long, default_value_t = 0.1)]
        kit_prob_min: f64,
            /// Write an HTML report to this path
        #[arg(long)]
        html: Option<String>,
    },
    /// Clean sequencing files (SAM/BAM/FASTQ/FASTQ.GZ) with kit-aware validation
    Clean {
        /// Thread count (0 = auto / all)
        #[arg(short = 't', long, default_value_t = 0)]
        threads: usize,
        /// Kit id (must match a known ONT kit, e.g. "LSK114")
        #[arg(short, long)]
        kit: String,
        /// Output FASTQ.GZ path
        #[arg(short, long, value_name = "OUT.fastq.gz")]
        output: std::path::PathBuf,
        /// One or more input files (SAM, BAM, FASTQ, or FASTQ.GZ)
        #[arg(value_name = "FILES", required = true)]
        files: Vec<std::path::PathBuf>,
    },

}

fn main() -> polars::prelude::PolarsResult<()> {
    let cli = Cli::parse();

    match cli.command {
        
        Commands::Clean { threads, kit, output, files } => {
            cmd_clean(threads, kit, output, files);
        }

        Commands::ListKits { format, full, truncate } => { 
            cmd_list_kits(format, full, truncate); 
        }

        Commands::Describe { id } => {
            cmd_describe(id);
        }

        Commands::Benchmark { files, kit, truth, algorithms, max_dist, threads, csv } => {
            use porkchop::benchmark::{self, BenchmarkAlgo};

            let algorithms = algorithms.to_lowercase();
            let algos = BenchmarkAlgo::from_list(algorithms.as_str());

            let mut rows = Vec::new();
let truth_map = match truth {
                Some(p) => benchmark::load_truth(p).ok(),
                None => None,
            };

            for file in files {
                for algo in &algos {
                    let kit_ref = match porkchop::get_sequences_for_kit(kit.as_str()) {
                        Some(k) => k,
                        None => {
                            return Err(polars::prelude::PolarsError::ComputeError(
                                format!("Unknown kit: {}", kit).into(),
                            ));
                        }
                    };

                    let (tp, fp, fn_, dur, nseq, cpu, _input_format) =
                        benchmark::benchmark_file(file.clone(), kit_ref, *algo, truth_map.clone(), threads, max_dist)
                        .map_err(|e| polars::prelude::PolarsError::ComputeError(e.to_string().into()))?;

                    rows.push((
                        file.clone(),
                        algo.as_str().to_string(),
                        tp,
                        fp,
                        fn_,
                        dur.as_millis(),
                        nseq,
                        cpu,
                        threads.unwrap_or(0),
                    ));
                }
            }

            if csv {
                let mut df = df!(
                    "file"        => rows.iter().map(|r| r.0.clone()).collect::<Vec<_>>(),
                    "algorithm"   => rows.iter().map(|r| r.1.clone()).collect::<Vec<_>>(),
                    "tp"          => rows.iter().map(|r| r.2 as u64).collect::<Vec<_>>(),
                    "fp"          => rows.iter().map(|r| r.3 as u64).collect::<Vec<_>>(),
                    "fn"          => rows.iter().map(|r| r.4 as u64).collect::<Vec<_>>(),
                    "elapsed_ms"  => rows.iter().map(|r| r.5 as u64).collect::<Vec<_>>(),
                    "nseq"        => rows.iter().map(|r| r.6 as u64).collect::<Vec<_>>(),
                    "cpu"         => rows.iter().map(|r| r.7 as f32).collect::<Vec<_>>(),
                    "threads"     => rows.iter().map(|r| r.8 as u64).collect::<Vec<_>>(),
                )?;
                let w = CsvWriter::new(std::io::stdout());
                w.include_header(true).finish(&mut df)?;
            } else {
                for (file, algo, tp, fp, fn_, dur, nseq, cpu, threads) in rows {
                    println!("{file}\t{algo}\tTP={tp}\tFP={fp}\tFN={fn_}\tms={dur}\tN={nseq}\tcpu={cpu}\tthreads={threads}");
                }
            }
        }

        Commands::Screen { files, algorithm, max_dist, fraction, tick, threads, json, kit_prob_min, html } => {
            let algo = match algorithm.parse::<porkchop::benchmark::BenchmarkAlgo>() {
                Ok(a) => a,
                Err(_) => porkchop::benchmark::BenchmarkAlgo::Edlib,
            };
            let opts = porkchop::screen::ScreenOpts {
                files,
                threads,
                fraction,
                tick_secs: tick,
                algo,
                max_dist,
                json,
                kit_prob_min,
                html,
            };
            if let Err(e) = porkchop::screen::run_screen(opts) {
                eprintln!("screen error: {e}");
            }
        }
    }

    Ok(())
}

fn cmd_list_kits(format: OutputFormat, full: bool, truncate: Option<usize>) {
    use porkchop::list_supported_kits;

    let kits = list_supported_kits();
    let ids: Vec<String> = kits.iter().map(|k| k.id.0.to_string()).collect();
    let descs: Vec<String> = kits.iter().map(|k| k.description.to_string()).collect();
    let legacy: Vec<bool> = kits.iter().map(|k| k.legacy).collect();
    let chems: Vec<String> = kits.iter().map(|k| k.chemistry.to_string()).collect();

    let df = df!( "kit" => ids,
                  "description" => descs,
                  "legacy" => legacy,
                  "chemistry" => chems, ).expect("dataframe");

    match format {
        OutputFormat::Csv => {
            let w = CsvWriter::new(std::io::stdout());
            w.include_header(true).finish(&mut df.clone()).expect("write csv");
        }
        OutputFormat::Md => {
            print_df_markdown(&df);
        }
        OutputFormat::Table => {
            // Pretty-print Polars DataFrame with configurable formatting.
            if full {
                std::env::set_var("POLARS_FMT_TABLE_FORMATTING", "UTF8_FULL");
                std::env::set_var("POLARS_FMT_MAX_COLS", "100000");
                std::env::set_var("POLARS_FMT_MAX_ROWS", "1000000");
                std::env::set_var("POLARS_FMT_STR_LEN", "1000000");
                std::env::set_var("POLARS_TABLE_WIDTH", "65535"); // safe cap for Polars 0.42
            } else {
                let trunc = truncate.unwrap_or(80).to_string();
                std::env::set_var("POLARS_FMT_TABLE_FORMATTING", "UTF8_FULL");
                std::env::set_var("POLARS_FMT_MAX_COLS", "200");
                std::env::set_var("POLARS_FMT_MAX_ROWS", "2000");
                std::env::set_var("POLARS_FMT_STR_LEN", &trunc);
                std::env::set_var("POLARS_TABLE_WIDTH", "200");
            }
            println!("{}", df);
        }
    }
}

fn print_df_markdown(df: &DataFrame) {

    let cols = df.get_columns();
    if cols.is_empty() {
        println!("(empty)");
        return;
    }

    print!("|");
    for s in cols {
        print!(" {} |", s.name());
    }
    println!();

    print!("|");
    for _ in cols.iter() {
        print!("---|");
    }
    println!();

    let height = df.height();
    for i in 0..height {
        print!("|");
        for s in cols {
            let v = s.get(i).map(|any| any.to_string()).unwrap_or_else(|_| "".to_string());
            let v = v.replace("|", r"\|");
            print!(" {} |", v);
        }
        println!();
    }
}


fn cmd_describe(id: String) {
    use porkchop::{get_sequences_for_kit, base_chemistry_of, kit_is_legacy};

    match get_sequences_for_kit(id.as_str()) {
        Some(kit) => {
            println!("id: {}", kit.id.0);
            println!("description: {}", kit.description);
            println!("legacy: {}", kit_is_legacy(kit));
            println!("chemistry: {}", base_chemistry_of(kit).to_string());
            println!("adapters/primers: {}", kit.adapters_and_primers.len());
            println!("barcodes: {}", kit.barcodes.len());
        }
        None => {
            eprintln!("Unknown kit: {}", id);
        }
    }
}






use rayon::prelude::*;

/// Validate that the kit exists in the registry; exit with code 2 if not.
fn ensure_known_kit(kit: &str) {
    if porkchop::get_sequences_for_kit(kit).is_none() {
        eprintln!("Unknown kit: {}. Use `porkchop list-kits --format table` to see valid kit ids.", kit);
        std::process::exit(2);
    }
}

/// Return (ok_files, bad_files) based on extension checks.
fn split_supported_files(paths: Vec<std::path::PathBuf>) -> (Vec<std::path::PathBuf>, Vec<std::path::PathBuf>) {
    let mut ok = Vec::new();
    let mut bad = Vec::new();
    for p in paths {
        let name = p.file_name().and_then(|s| s.to_str()).unwrap_or("").to_ascii_lowercase();
        let ext = p.extension().and_then(|s| s.to_str()).unwrap_or("").to_ascii_lowercase();
        let is_ok =
            name.ends_with(".fastq.gz") ||
            name.ends_with(".fq.gz") ||
            ext == "fastq" || ext == "fq" ||
            ext == "sam" || ext == "bam";
        if is_ok { ok.push(p); } else { bad.push(p); }
    }
    (ok, bad)
}

#[derive(Clone)]
struct OwnedRecord {
    id: String,
    seq: Vec<u8>,
    qual: Vec<u8>,
}

fn write_fastq_record<W: std::io::Write>(w: &mut W, id: &str, seq: &[u8], qual: &[u8]) -> std::io::Result<()> {
    w.write_all(b"@")?;
    w.write_all(id.as_bytes())?;
    w.write_all(b"
")?;
    w.write_all(seq)?;
    w.write_all(b"
+
")?;
    w.write_all(qual)?;
    w.write_all(b"
")?;
    Ok(())
}

// --- Edlib wrapper ---
mod edwrap {
    use edlib_rs::edlibrs::{EdlibAlignConfigRs, EdlibAlignModeRs, EdlibAlignTaskRs, EdlibEqualityPairRs, edlibAlignRs};
    pub struct Hit { pub start: i32, pub end: i32, pub edits: i32 }
    pub fn locate(pattern: &[u8], text: &[u8], max_edits: i32) -> Option<Hit> {
        let empty: &[EdlibEqualityPairRs] = &[];
        let cfg = EdlibAlignConfigRs {
            k: max_edits,
            mode: EdlibAlignModeRs::EDLIB_MODE_HW,
            task: EdlibAlignTaskRs::EDLIB_TASK_LOC,
            additionalequalities: empty,
        };
        let res = edlibAlignRs(pattern, text, &cfg);
        if res.editDistance < 0 { return None; }
        let start = res.startLocations.as_ref()?.get(0).copied()?;
        let end = res.endLocations.as_ref()?.get(0).copied()?;
        Some(Hit { start, end, edits: res.editDistance })
    }
}

#[derive(Clone)]
struct Motif<'a> {
    name: &'a str,
    kind: &'a str,
    seq: &'a [u8],
}

fn motifs_for_kit<'a>(kit: &'a porkchop::kit::Kit) -> Vec<Motif<'a>> {
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

fn max_edits_for(len: usize) -> i32 {
    let m = (len as f64 * 0.15).ceil() as i32;
    if m < 1 { 1 } else { m }
}

fn annotate_and_trim_one(seq: &[u8], qual: &[u8], kit_id: &str, motifs: &[Motif]) -> OwnedRecord {
    let s = normalize_seq(seq);
    let n = s.len() as i32;
    let mut left_best: Option<(i32, i32, i32, &str)> = None;  // (start,end,edits,name)
    let mut right_best: Option<(i32, i32, i32, &str)> = None;

    for m in motifs {
        let maxk = max_edits_for(m.seq.len()) as i32;
        if let Some(hit) = edwrap::locate(m.seq, &s, maxk) {
            // classify based on position
            let center = (hit.start + hit.end) / 2;
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
    }

    let mut left_cut: i32 = 0;
    let mut right_cut: i32 = n;

    let mut notes: Vec<String> = Vec::new();
    if let Some((st, en, ed, nm)) = left_best {
        left_cut = en + 1;
        notes.push(format!("L:{}:{}-{}:ed={}", nm, st, en, ed));
    }
    if let Some((st, en, ed, nm)) = right_best {
        right_cut = st;
        notes.push(format!("R:{}:{}-{}:ed={}", nm, st, en, ed));
    }
    if left_cut < 0 { left_cut = 0; }
    if right_cut > n { right_cut = n; }
    if left_cut >= right_cut { left_cut = 0; right_cut = n; }

    let start = left_cut as usize;
    let end = right_cut as usize;
    let new_seq = s[start..end].to_vec();
    let new_qual = if !qual.is_empty() {
        qual[start..end].to_vec()
    } else {
        vec![b'I'; new_seq.len()]
    };

    let id = format!("kit={};trim={}..{};{}", kit_id, left_cut, right_cut, notes.join(";"));
    OwnedRecord { id, seq: new_seq, qual: new_qual }
}

fn process_fastx_to_gz(out_path: &std::path::Path, input_files: Vec<std::path::PathBuf>, _threads_eff: usize, kit_id: &str) -> anyhow::Result<()> {
    use std::fs::File;
    use std::io::BufWriter;
    use needletail::parser::parse_fastx_file;

    let kit = porkchop::get_sequences_for_kit(kit_id).expect("validated kit");
    let motifs = motifs_for_kit(kit);

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
                    let processed: Vec<OwnedRecord> = buf.par_iter().map(|r| {
                        //let id = std::str::from_utf8(r.qname()).unwrap_or("SAM").to_string();
                        let seq = r.seq().as_bytes();
                        let qualv = r.qual().to_vec();
                        let qual = qualv.into_iter().map(|q| (q + 33) as u8).collect::<Vec<u8>>();
                        annotate_and_trim_one(&seq, &qual, kit_id, &motifs)
                    }).collect();
                    for pr in &processed { write_fastq_record(&mut gz, &pr.id, &pr.seq, &pr.qual)?; }
                    buf.clear();
                }
            }
            if !buf.is_empty() {
                let processed: Vec<OwnedRecord> = buf.par_iter().map(|r| {
                    //let id = std::str::from_utf8(r.qname()).unwrap_or("SAM").to_string();
                    let seq = r.seq().as_bytes();
                    let qualv = r.qual().to_vec();
                    let qual = qualv.into_iter().map(|q| (q + 33) as u8).collect::<Vec<u8>>();
                    annotate_and_trim_one(&seq, &qual, kit_id, &motifs)
                }).collect();
                for pr in &processed { write_fastq_record(&mut gz, &pr.id, &pr.seq, &pr.qual)?; }
            }
        } else if lower.ends_with(".bam") {
            use rust_htslib::bam::{self, Read};
            let mut reader = bam::Reader::from_path(&path)?;
            let mut buf: Vec<rust_htslib::bam::Record> = Vec::new();
            for r in reader.records() {
                if let Ok(rec) = r { buf.push(rec); }
                if buf.len() >= CHUNK {
                    let processed: Vec<OwnedRecord> = buf.par_iter().map(|r| {
                        //let id = std::str::from_utf8(r.qname()).unwrap_or("BAM").to_string();
                        let seq = r.seq().as_bytes();
                        let qualv = r.qual().to_vec();
                        let qual = qualv.into_iter().map(|q| (q + 33) as u8).collect::<Vec<u8>>();
                        annotate_and_trim_one(&seq, &qual, kit_id, &motifs)
                    }).collect();
                    for pr in &processed { write_fastq_record(&mut gz, &pr.id, &pr.seq, &pr.qual)?; }
                    buf.clear();
                }
            }
            if !buf.is_empty() {
                let processed: Vec<OwnedRecord> = buf.par_iter().map(|r| {
                    //let id = std::str::from_utf8(r.qname()).unwrap_or("BAM").to_string();
                    let seq = r.seq().as_bytes();
                    let qualv = r.qual().to_vec();
                    let qual = qualv.into_iter().map(|q| (q + 33) as u8).collect::<Vec<u8>>();
                    annotate_and_trim_one(&seq, &qual, kit_id, &motifs)
                }).collect();
                for pr in &processed { write_fastq_record(&mut gz, &pr.id, &pr.seq, &pr.qual)?; }
            }
        } else {
            // FASTA/FASTQ (optionally gz) via needletail
            let mut reader = parse_fastx_file(&path)?;
            loop {
                let mut owned_chunk: Vec<OwnedRecord> = Vec::with_capacity(CHUNK);
                for _ in 0..CHUNK {
                    match reader.next() {
                        Some(Ok(record)) => {
                            let id = String::from_utf8_lossy(record.id()).to_string();
                            let seq = record.seq().to_vec();
                            let qual = record.qual().map(|q| q.to_vec()).unwrap_or_else(|| vec![b'I'; seq.len()]);
                            owned_chunk.push(OwnedRecord { id, seq, qual });
                        }
                        Some(Err(_e)) => continue,
                        None => break,
                    }
                }
                if owned_chunk.is_empty() { break; }
                let processed: Vec<OwnedRecord> = owned_chunk.par_iter().map(|r| {
                    annotate_and_trim_one(&r.seq, &r.qual, kit_id, &motifs)
                }).collect();
                for pr in &processed { write_fastq_record(&mut gz, &pr.id, &pr.seq, &pr.qual)?; }
            }
        }
    }

    gz.finish()?;
    Ok(())
}


/// Implementation for `porkchop clean`
fn cmd_clean(threads: usize, kit: String, output: std::path::PathBuf, files: Vec<std::path::PathBuf>) {
    // Validate kit id
    ensure_known_kit(&kit);

    // Validate file extensions
    let (ok, bad) = split_supported_files(files);
    if !bad.is_empty() {
        eprintln!("Unsupported file type(s):");
        for p in &bad { eprintln!("  - {}", p.display()); }
        eprintln!("Allowed: SAM (.sam), BAM (.bam), FASTQ (.fastq/.fq), and gzipped FASTQ (.fastq.gz/.fq.gz).");
        std::process::exit(2);
    }

    // Determine effective thread count
    let threads_eff = if threads == 0 { std::cmp::max(1, num_cpus::get()) } else { threads };
    rayon::ThreadPoolBuilder::new().num_threads(threads_eff).build_global().ok();

    eprintln!("clean: kit={} | threads={} | inputs={} | output={}", kit, threads_eff, ok.len(), output.display());

    if let Err(e) = process_fastx_to_gz(&output, ok, threads_eff, &kit) {
        eprintln!("clean error: {:?}", e);
        std::process::exit(1);
    }
}

