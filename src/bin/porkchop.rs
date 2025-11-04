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

    /// Clean sequencing files and show a live dashboard of clip modalities
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
        /// One or more input files (SAM, BAM, FASTQ, FASTQ.GZ)
        #[arg(value_name = "FILES", required = true)]
        files: Vec<std::path::PathBuf>,
    },

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
}

fn main() -> polars::prelude::PolarsResult<()> {
    let cli = Cli::parse();

    match cli.command {
            Commands::Clean { threads, kit, output, files } => {
        if let Err(e) = porkchop::clean::run(threads, &kit, &output, files) {
            eprintln!("clean error: {:?}", e);
            std::process::exit(1);
        }
    },
Commands::ListKits { format, full, truncate } => { cmd_list_kits(format, full, truncate); }

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
