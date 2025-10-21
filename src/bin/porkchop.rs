
/*!
Command-line interface for the **porkchop** crate.

This binary is intentionally tiny and delegates logic to the library crate.
It uses **clap** for subcommand parsing and **polars** to build DataFrames,
but prints non-truncated tables by default for readability.
*/

use clap::{Parser, Subcommand, ArgAction};
use polars::prelude::*;
use rayon::prelude::*;
use porkchop::{self, BaseChemistry};

/// Main CLI entrypoint for `porkchop` utilities.
#[derive(Parser, Debug)]
#[command(name = "porkchop", version, about = "Utilities for the porkchop ONT reference")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

/// Available subcommands.
#[derive(Subcommand, Debug)]
enum Commands {
    /// Sequence IO utilities (FASTQ/SAM/BAM)
    Reads {
        /// Count records in one or more input files (FASTQ/FASTQ.GZ/SAM/BAM)
        #[arg(long, action=clap::ArgAction::SetTrue)]
        count: bool,
        /// Input files (one or more)
        #[arg(required=true)]
        files: Vec<String>,
        /// Number of threads (defaults to all logical cores)
        #[arg(long)]
        threads: Option<usize>,
    },
    /// List all supported kits in a table (no truncation).
    ListKits {
        /// Emit CSV instead of the full-width table.
        #[arg(long, action=ArgAction::SetTrue)]
        csv: bool,
    },
    /// Describe a specific kit: list all primers, adapters and barcodes for that kit.
    DescribeKit {
        /// Kit identifier (e.g., LSK114, NBD114.24, RBK114.96).
        kit: String,
        /// Emit CSV instead of the full-width table.
        #[arg(long, action=ArgAction::SetTrue)]
        csv: bool,
    },
}

fn main() -> polars::prelude::PolarsResult<()> {
    let cli = Cli::parse();
    match cli.command {
    Commands::Reads { count, files, threads } => {
        if count {
            return cmd_reads_count(&files, threads);
        } else {
            eprintln!("Specify a reads subcommand (e.g., --count).");
        }
    },
        Commands::ListKits { csv } => cmd_list_kits(csv)?,
        Commands::DescribeKit { kit, csv } => cmd_describe_kit(&kit, csv)?,
    }
    Ok(())
}

/// Build and print the kit table using Polars.
///
/// By default, we print a **non-truncated** table (no ellipses) by constructing
/// it ourselves from the DataFrame columns. Use `--csv` for machine-readable output.
fn cmd_list_kits(as_csv: bool) -> PolarsResult<()> {
    let rows = porkchop::list_kits_rows();

    // Columns backed by the library API
    let kit_id: Vec<String> = rows.iter().map(|r| r.0.clone()).collect();
    let legacy: Vec<String> = rows.iter().map(|r| r.2.to_string()).collect();
    let base_chem: Vec<String> = rows.iter().map(|r| match r.3 { BaseChemistry::Rapid => "rapid".to_string(), BaseChemistry::Ligation => "ligation".to_string() }).collect();
    let description: Vec<String> = rows.iter().map(|r| r.1.clone()).collect();

    // Also expose as a DataFrame (for potential future operations).
    let _df = df!(
        "kit_id" => kit_id.clone(),
        "legacy" => legacy.iter().map(|s| s.as_str()).collect::<Vec<&str>>(),
        "base_chemistry" => base_chem.iter().map(|s| s.as_str()).collect::<Vec<&str>>(),
        "description" => description.iter().map(|s| s.as_str()).collect::<Vec<&str>>()
    )?;

    if as_csv {
        // Machine-readable, non-truncated CSV
        let headers = ["kit_id", "legacy", "base_chemistry", "description"];
        println!("{}", headers.join(","));
        for i in 0..kit_id.len() {
            let row = [&kit_id[i], &legacy[i], &base_chem[i], &description[i]];
            let mut out: Vec<String> = Vec::with_capacity(row.len());
            for cell in row {
                let mut v = cell.clone();
                if v.contains(',') || v.contains('\"') {
                    v = format!("\"{}\"", v.replace('\"', "\"\""));
                }
                out.push(v);
            }
            println!("{}", out.join(","));
        }
        return Ok(());
    }

    // Full-width, non-truncated table printing.
    // Compute column widths (at least header width).
    let headers = ["kit_id", "legacy", "base_chemistry", "description"];
    let cols = [&kit_id, &legacy, &base_chem, &description];

    let mut widths: Vec<usize> = headers.iter().map(|h| h.len()).collect();
    for (ci, col) in cols.iter().enumerate() {
        for cell in col.iter() {
            let w = cell.chars().count();
            if w > widths[ci] { widths[ci] = w; }
        }
    }

    // Render header
    let mut line = String::new();
    for (i, h) in headers.iter().enumerate() {
        if i > 0 { line.push_str(" | "); }
        line.push_str(&format!("{:width$}", h, width = widths[i]));
    }
    println!("{}", line);
    // separator
    let mut sep = String::new();
    for (i, _) in headers.iter().enumerate() {
        if i > 0 { sep.push_str("-+-"); }
        sep.push_str(&"-".repeat(widths[i]));
    }
    println!("{}", sep);

    // rows
    for r in 0..kit_id.len() {
        let row = [&kit_id[r], &legacy[r], &base_chem[r], &description[r]];
        let mut line = String::new();
        for (i, cell) in row.iter().enumerate() {
            if i > 0 { line.push_str(" | "); }
            line.push_str(&format!("{:width$}", cell, width = widths[i]));
        }
        println!("{}", line);
    }

    Ok(())
}

/// Describe a kit by printing a table of all **primers**, **adapters**, and **barcodes**.
/// The output is **not truncated**; every cell is printed in full.

/// Describe a kit by printing a table of all **primers**, **adapters**, and **barcodes**,
/// including **provenance URL** and **reference** columns. The output is **not truncated**.
fn cmd_describe_kit(kit_id: &str, as_csv: bool) -> PolarsResult<()> {
    let rows = porkchop::kit_elements_rows_with_provenance(kit_id)
        .ok_or_else(|| PolarsError::NoData("unknown kit id".into()))?;

    // Split columns
    let name: Vec<&str> = rows.iter().map(|r| r.0.as_str()).collect();
    let kind: Vec<&str> = rows.iter().map(|r| r.1.as_str()).collect();
    let sequence: Vec<&str> = rows.iter().map(|r| r.2.as_str()).collect();
    let prov_url: Vec<&str> = rows.iter().map(|r| r.3.as_str()).collect();
    let prov_ref: Vec<&str> = rows.iter().map(|r| r.4.as_str()).collect();

    // DF (not required for printing, but kept for parity)
    let _df = df!(
        "name" => name.clone(),
        "kind" => kind.clone(),
        "sequence" => sequence.clone(),
        "provenance_url" => prov_url.clone(),
        "provenance_ref" => prov_ref.clone()
    )?;

    if as_csv {
        let headers = ["name", "kind", "sequence", "provenance_url", "provenance_ref"];
        println!("{}", headers.join(","));
        for i in 0..rows.len() {
            let row = [name[i], kind[i], sequence[i], prov_url[i], prov_ref[i]];
            let mut out: Vec<String> = Vec::with_capacity(row.len());
            for cell in row {
                let mut v = cell.to_string();
                if v.contains(',') || v.contains('\"') {
                    v = format!("\"{}\"", v.replace('\"', "\"\""));
                }
                out.push(v);
            }
            println!("{}", out.join(","));
        }
        return Ok(());
    }

    // Non-truncated pretty table with provenance columns.
    let headers = ["name", "kind", "sequence", "provenance_url", "provenance_ref"];
    let cols: [&Vec<&str>; 5] = [&name, &kind, &sequence, &prov_url, &prov_ref];
    let mut widths: Vec<usize> = headers.iter().map(|h| h.len()).collect();
    for (ci, col) in cols.iter().enumerate() {
        for cell in col.iter() {
            let w = cell.chars().count();
            if w > widths[ci] { widths[ci] = w; }
        }
    }
    // header
    let mut line = String::new();
    for (i, h) in headers.iter().enumerate() {
        if i > 0 { line.push_str(" | "); }
        line.push_str(&format!("{:width$}", h, width = widths[i]));
    }
    println!("{}", line);
    // sep
    let mut sep = String::new();
    for (i, _) in headers.iter().enumerate() {
        if i > 0 { sep.push_str("-+-"); }
        sep.push_str(&"-".repeat(widths[i]));
    }
    println!("{}", sep);
    // rows
    for r in 0..rows.len() {
        let row = [name[r], kind[r], sequence[r], prov_url[r], prov_ref[r]];
        let mut line = String::new();
        for (i, cell) in row.iter().enumerate() {
            if i > 0 { line.push_str(" | "); }
            line.push_str(&format!("{:width$}", cell, width = widths[i]));
        }
        println!("{}", line);
    }

    Ok(())
}


/// Count reads across input files in parallel; prints per-file and total counts.

fn cmd_reads_count(files: &Vec<String>, threads: Option<usize>) -> polars::prelude::PolarsResult<()> {
    use porkchop::seqio;
    let results: Vec<(String, usize)> = files.par_iter().map(|p| {
        match seqio::for_each_parallel(p, threads, |_r| {}) {
            Ok((_fmt, n)) => (p.clone(), n),
            Err(_) => (p.clone(), 0),
        }
    }).collect();

    let mut total = 0usize;
    for (p, n) in &results {
        println!("{}	{}", p, n);
        total += *n;
    }
    println!("TOTAL	{}", total);
    Ok(())
}

