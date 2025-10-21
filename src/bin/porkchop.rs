
/*!
Command-line interface for the **porkchop** crate.

This binary is intentionally tiny and delegates logic to the library crate.
It uses **clap** for subcommand parsing and **polars** to build a DataFrame, but
prints a non-truncated, full-width table by default for readability.
*/

use clap::{Parser, Subcommand, ArgAction};
use polars::prelude::*;
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
    /// List all supported kits in a table (no truncation).
    ListKits {
        /// Emit CSV instead of the full-width table.
        #[arg(long, action=ArgAction::SetTrue)]
        csv: bool,
    }
}

fn main() -> polars::prelude::PolarsResult<()> {
    let cli = Cli::parse();
    match cli.command {
        Commands::ListKits { csv } => cmd_list_kits(csv)?,
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

    // Also expose as a DataFrame (for potential future operations) as requested.
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
            let w = cell.chars().count(); // character count (not perfect for grapheme clusters, but fine)
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
