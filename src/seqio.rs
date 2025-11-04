
//! High‑performance IO for **FASTQ / FASTQ.GZ / SAM / BAM**.
//!
//! ### Design
//! - **FASTQ/FASTQ.GZ** parsed with `needletail`
//! - **SAM/BAM** parsed with `rust-htslib` (optionally multithreaded via `set_threads`)
//! - **Parallelism**: uses a local Rayon pool; `threads = None` uses all logical cores.
//!
//! ### Callback contract
//! The `on_record` callback must be `Fn(NARead) + Send + Sync + 'static` and should be
//! **fast**; heavy work should batch/accumulate in outer scope to avoid blocking IO.
//!
//! ### Errors
//! Parsing/IO errors are bubbled via `anyhow::Result` to the caller.
//!
//! ### Example
//! ```no_run
//! use porkchop::seqio;
//! let (_fmt, n) = seqio::for_each_parallel("reads.fastq.gz", Some(16), |r| {
//!     // r.id, r.seq, r.qual
//! }).unwrap();
//! println!("processed {n} records");
//! ```
//!
//!
//! `for_each_parallel` detects format and iterates records, invoking a user callback.
//! BAM/SAM via rust-htslib; FASTQ/FASTQ.GZ via needletail.
//!
//! The callback must be `Fn(NARead) + Send + Sync + 'static`.
//! Parallelism uses Rayon; `--threads` controls thread count.

use std::path::Path;
use anyhow::Result;
use rayon::ThreadPoolBuilder;
use needletail::parse_fastx_file;
use rust_htslib::bam;
use rust_htslib::bam::Read;

/// Input format detected from path.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputFormat { Fastq, Bam, Sam }

/// A normalized read passed to callbacks.
#[derive(Debug, Clone)]
pub struct NARead {
    pub id: String,
    pub seq: Vec<u8>,
    pub qual: Option<Vec<u8>>,
}

/// Core driver: parse and iterate records, potentially in parallel (rayon pool size).
/// fn `for_each_parallel` — auto‑generated rustdoc.
pub fn for_each_parallel<P, F>(path: P, threads: Option<usize>, on_record: F) -> Result<(InputFormat, usize)>
where
    P: AsRef<Path>,
    F: Fn(NARead) + Send + Sync + 'static,
{
    let p = path.as_ref();
    let fmt = if let Some(ext) = p.extension().and_then(|s| s.to_str()) {
        match ext.to_ascii_lowercase().as_str() {
            "fq" | "fastq" | "gz" => InputFormat::Fastq,
            "bam" => InputFormat::Bam,
            "sam" => InputFormat::Sam,
            _ => {
                if p.to_string_lossy().contains("fastq") || p.to_string_lossy().contains(".fq.") {
                    InputFormat::Fastq
                } else { InputFormat::Bam }
            }
        }
    } else { InputFormat::Fastq };

    let n = threads.unwrap_or_else(num_cpus::get).max(1);
    let pool = ThreadPoolBuilder::new().num_threads(n).build()?;

    let counter = std::sync::atomic::AtomicUsize::new(0);
    let cb = &on_record;

    pool.install(|| -> Result<()> {
        match fmt {
            InputFormat::Fastq => {
                let mut reader = parse_fastx_file(p)?;
                while let Some(record) = reader.next() {
                    let rec = record?;
                    let id = String::from_utf8_lossy(rec.id()).to_string();
                    let seq = rec.seq().to_vec();
                    let qual = rec.qual().map(|q| q.to_vec());
                    let naread = NARead { id, seq, qual };
                    cb(naread);
                    counter.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                }
            }
            InputFormat::Bam | InputFormat::Sam => {
                let mut reader = bam::Reader::from_path(p)?;
                if n > 1 { let _ = reader.set_threads(n); }
for result in reader.records() {
                let rec = result?;
                    let id = String::from_utf8_lossy(rec.qname()).to_string();
                    let seq = rec.seq().as_bytes();
                    let qual = {
                        let q = rec.qual();
                        if q.is_empty() { None } else { Some(q.to_vec()) }
                    };
                    let naread = NARead { id, seq, qual };
                    cb(naread);
                    counter.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                }
            }
        }
        Ok(())
    })?;

    Ok((fmt, counter.load(std::sync::atomic::Ordering::Relaxed)))
}
