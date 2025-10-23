
use clap::{Parser, Subcommand};
use polars::prelude::*;

#[derive(Debug, Parser)]
#[command(name="porkchop")]
#[command(about="ONT adapters/primers/barcodes utilities")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Debug, Subcommand)]
enum Commands {
    /// List all kits known to porkchop
    ListKits,
    /// Describe all sequences and barcodes for a kit (no truncation; includes provenance)
    DescribeKit {
        /// Kit id (e.g., LSK114, PCS114, NBD114.24)
        #[arg(long, required=true)]
        kit: String,
        /// Emit CSV instead of table
        #[arg(long, action=clap::ArgAction::SetTrue)]
        csv: bool,
    },
    /// Benchmark classification algorithms against a truth set (optional).
    Benchmark {
        /// Input files (FASTQ/FASTQ.GZ/SAM/BAM)
        #[arg(required=true)]
        files: Vec<String>,
        /// Kit ID (e.g., LSK114, PCS114, NBD114.24)
        #[arg(long, required=true)]
        kit: String,
        /// Optional truth set CSV with columns: read_id, expected_labels, [kind]
        #[arg(long)]
        truth: Option<String>,
        /// BenchmarkAlgo list (comma-separated): aho,myers,edlib,parasail,ac+myers,ac+parasail
        #[arg(long, default_value = "aho,myers,edlib,parasail,ac+myers,ac+parasail")]
        algorithms: String,
        /// Max edit distance for approximate matchers
        #[arg(long, default_value_t = 5)]
        max_dist: usize,
        /// Threads to use (defaults to all cores)
        #[arg(long)]
        threads: Option<usize>,
        /// Emit CSV instead of pretty table
        #[arg(long, action=clap::ArgAction::SetTrue)]
        csv: bool,
    },
}

fn main() -> polars::prelude::PolarsResult<()> {
    let cli = Cli::parse();
    match cli.command {
        Commands::ListKits => cmd_list_kits(),
        Commands::DescribeKit { kit, csv } => cmd_describe(&kit, csv),
        Commands::Benchmark { files, kit, truth, algorithms, max_dist, threads, csv } => {
            cmd_benchmark(&files, &kit, truth.as_deref(), &algorithms, max_dist, threads, csv)
        }
    }
}

fn cmd_list_kits() -> polars::prelude::PolarsResult<()> {
    let kits = porkchop::list_supported_kits();
    let ids: Vec<&str> = kits.iter().map(|k| k.id.0).collect();
    let desc: Vec<&str> = kits.iter().map(|k| k.description).collect();
    let legacy: Vec<bool> = kits.iter().map(|k| porkchop::kit_is_legacy(k)).collect();
    let chemistry: Vec<&str> = kits.iter().map(|k| match porkchop::base_chemistry_of(k) {
        porkchop::BaseChemistry::Rapid => "Rapid",
        porkchop::BaseChemistry::Ligation => "Ligation",
        porkchop::BaseChemistry::PCRcDNA => "PCRcDNA",
        porkchop::BaseChemistry::Amplicon => "Amplicon",
    }).collect();

    let df = df!(
        "kit" => ids,
        "description" => desc,
        "legacy" => legacy,
        "chemistry" => chemistry
    )?;

    print_table(&df);
    Ok(())
}

fn cmd_describe(kit: &str, as_csv: bool) -> polars::prelude::PolarsResult<()> {
    let k = porkchop::get_sequences_for_kit(kit).ok_or_else(|| PolarsError::ComputeError(format!("Unknown kit: {}", kit).into()))?;

    let mut names = Vec::new();
    let mut kinds = Vec::new();
    let mut seqs = Vec::new();
    let mut prov_src = Vec::new();
    let mut prov_app = Vec::new();

    for r in k.adapters_and_primers.iter().chain(k.barcodes.iter()) {
        names.push(r.name);
        kinds.push(match r.kind {
            porkchop::SeqKind::AdapterTop => "AdapterTop",
            porkchop::SeqKind::AdapterBottom => "AdapterBottom",
            porkchop::SeqKind::Primer => "Primer",
            porkchop::SeqKind::Barcode => "Barcode",
            porkchop::SeqKind::Flank => todo!(),
        });
        seqs.push(r.sequence);
        prov_src.push(r.provenance.source);
        prov_app.push(r.provenance.appendix.unwrap_or(""));
    }

    let mut df = df!(
        "name" => names,
        "kind" => kinds,
        "sequence" => seqs,
        "provenance_source" => prov_src,
        "provenance_appendix" => prov_app
    )?;

    if as_csv {
        let mut w = CsvWriter::new(std::io::stdout());
        w.include_header(true).finish(&mut df)?;
    } else {
        print_table(&df);
    }
    Ok(())
}

fn cmd_benchmark(
    files: &Vec<String>,
    kit: &str,
    truth: Option<&str>,
    algorithms: &str,
    max_dist: usize,
    threads: Option<usize>,
    as_csv: bool
) -> polars::prelude::PolarsResult<()> {
    use porkchop::benchmark::{self, BenchmarkAlgo};

    let truth_map = match truth { Some(p) => benchmark::load_truth(p).ok(), None => None };
    let algos = BenchmarkAlgo::from_list(algorithms);

    let mut rows: Vec<(String,String,u64,u64,u64,u128,usize,f32,usize)> = Vec::new();
    for file in files {
        for algo in &algos {
            let (tp, fp, fn_, dur, nseq, cpu, input_format) = benchmark::benchmark_file(file, &kit, *algo, truth_map.map(|m| m.clone()), threads, max_dist)
                .map_err(|e| PolarsError::ComputeError(e.to_string().into()))?;
            rows.push((file.clone(), algo.as_str().to_string(), tp, fp, fn_, dur.as_millis(), nseq, cpu, threads.unwrap_or_else(num_cpus::get)));
        }
    }

    let df = df!(
        "file" => rows.iter().map(|r| r.0.as_str()).collect::<Vec<&str>>(),
        "algorithm" => rows.iter().map(|r| r.1.as_str()).collect::<Vec<&str>>(),
        "tp" => rows.iter().map(|r| r.2).collect::<Vec<u64>>(),
        "fp" => rows.iter().map(|r| r.3).collect::<Vec<u64>>(),
        "fn" => rows.iter().map(|r| r.4).collect::<Vec<u64>>(),
        "precision" => rows.iter().map(|r| { let tp=r.2 as f64; let fp=r.3 as f64; if tp+fp>0.0 { tp/(tp+fp) } else { f64::NAN } }).collect::<Vec<f64>>(),
        "recall" => rows.iter().map(|r| { let tp=r.2 as f64; let fn_=r.4 as f64; if tp+fn_>0.0 { tp/(tp+fn_) } else { f64::NAN } }).collect::<Vec<f64>>(),
        "f1" => {
            let p: Vec<f64> = rows.iter().map(|r| { let tp=r.2 as f64; let fp=r.3 as f64; if tp+fp>0.0 { tp/(tp+fp) } else { f64::NAN } }).collect();
            let q: Vec<f64> = rows.iter().map(|r| { let tp=r.2 as f64; let fn_=r.4 as f64; if tp+fn_>0.0 { tp/(tp+fn_) } else { f64::NAN } }).collect();
            p.iter().zip(q.iter()).map(|(pp,rr)| if pp.is_nan()||rr.is_nan()||(*pp+*rr)==0.0 { f64::NAN } else { 2.0*pp*rr/(pp+rr) }).collect::<Vec<f64>>()
        },
        "time_ms" => rows.iter().map(|r| r.5).map(|v| v as u64).collect::<Vec<u64>>(),
        "time_per_seq_us" => rows.iter().map(|r| if r.6>0 { (r.5 as f64 * 1000.0)/(r.6 as f64) } else { f64::NAN }).collect::<Vec<f64>>(),
        "nseq" => rows.iter().map(|r| r.6 as u64).collect::<Vec<u64>>(),
        "cpu_mean_pct" => rows.iter().map(|r| r.7).collect::<Vec<f32>>(),
        "threads" => rows.iter().map(|r| r.8 as u64).collect::<Vec<u64>>()
    )?;

    if as_csv {
        let mut w = CsvWriter::new(std::io::stdout());
        w.include_header(true).finish(&mut df)?;
    } else {
        print_table(&df);
    }
    Ok(())
}

/// Pretty table printing without truncation.
fn print_table(df: &DataFrame) {
    let headers = df.get_column_names_owned();
    let mut cols: Vec<Vec<String>> = Vec::new();
    for name in &headers {
        let s = df.column(name).unwrap();
        let v = s.len();
        let mut out = Vec::with_capacity(v);
        for i in 0..v { out.push(s.get(i).unwrap().to_string()); }
        cols.push(out);
    }
    let mut widths: Vec<usize> = headers.iter().map(|h| h.len()).collect();
    for (ci, col) in cols.iter().enumerate() {
        for cell in col { widths[ci] = widths[ci].max(cell.chars().count()); }
    }
    // header
    {
        let mut line = String::new();
        for (i, h) in headers.iter().enumerate() {
            if i>0 { line.push_str(" | "); }
            line.push_str(&format!("{:width$}", h, width=widths[i]));
        }
        println!("{}", line);
        let mut sep = String::new();
        for (i, _) in headers.iter().enumerate() {
            if i>0 { sep.push_str("-+-"); }
            sep.push_str(&"-".repeat(widths[i]));
        }
        println!("{}", sep);
    }
    for r in 0..df.height() {
        let mut line = String::new();
        for (i, _) in headers.iter().enumerate() {
            let cell = &cols[i][r];
            if i>0 { line.push_str(" | "); }
            line.push_str(&format!("{:width$}", cell, width=widths[i]));
        }
        println!("{}", line);
    }
}
