
/*!
High‑performance Sequence I/O:
- FASTQ/FASTQ.GZ via **needletail**
- SAM/BAM via **rust-htslib**
- Multithreading via **rayon** (defaults to all logical CPUs)

You can override the number of threads via the CLI `--threads` flag or by
calling [`with_threadpool`] directly.
*/

use std::path::Path;

use rayon::{ThreadPoolBuilder, scope};
use needletail::parse_fastx_file;
use rust_htslib::bam;
use rust_htslib::bam::Read as BamRead;

/// Representation of a sequencing read (subset). Qualities are optional.
#[derive(Debug, Clone)]
pub struct NARead {
    pub id: String,
    pub seq: Vec<u8>,          // raw bases as bytes (ACGTN)
    pub qual: Option<Vec<u8>>, // ASCII Phred qualities if present
}

/// Input format, inferred from file extension.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputFormat { Fastq, Sam, Bam }

/// Detect the input format from the file extension (case-insensitive).
pub fn detect_format<P: AsRef<Path>>(p: P) -> Option<InputFormat> {
    let s = p.as_ref().to_string_lossy().to_lowercase();
    if s.ends_with(".fastq") || s.ends_with(".fq") || s.ends_with(".fastq.gz") || s.ends_with(".fq.gz") {
        Some(InputFormat::Fastq)
    } else if s.ends_with(".bam") {
        Some(InputFormat::Bam)
    } else if s.ends_with(".sam") {
        Some(InputFormat::Sam)
    } else {
        None
    }
}

/// Run a closure inside a rayon threadpool of the requested size (defaults to all CPUs).
pub fn with_threadpool<T: Send, F: FnOnce() -> T + Send>(threads: Option<usize>, f: F) -> T {
    let n = threads.unwrap_or_else(num_cpus::get);
    ThreadPoolBuilder::new().num_threads(n).build().unwrap().install(f)
}

/// Iterate FASTQ/FASTQ.GZ records and call `on_record` using a scoped spawn.
/// The scope guarantees all spawned tasks finish before returning.
pub fn fastq_for_each_parallel<P: AsRef<Path>, F>(path: P, threads: Option<usize>, on_record: F) -> anyhow::Result<usize>
where
    F: Fn(NARead) + Send + Sync,
{
    let mut reader = parse_fastx_file(path.as_ref()).map_err(|e| anyhow::anyhow!(e))?;
    let counter = std::sync::atomic::AtomicUsize::new(0);
    let mut first_err: Option<anyhow::Error> = None;

    // Drive parsing sequentially; fan out user work with scoped tasks
    scope(|s| {
        while let Some(record) = reader.next() {
            match record {
                Ok(rec) => {
                    let id = String::from_utf8_lossy(rec.id()).to_string();
                    let seq = rec.seq().to_vec();
                    let qual = rec.qual().map(|q| q.to_vec());
                    let naread = NARead { id, seq, qual };
                    s.spawn(|_| {
                        on_record(naread);
                    });
                    counter.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                }
                Err(e) => {
                    first_err = Some(anyhow::anyhow!(e));
                    break;
                }
            }
        }
    });

    if let Some(e) = first_err { return Err(e); }
    Ok(counter.into_inner())
}

/// Iterate SAM/BAM records and call `on_record` using a scoped spawn.
pub fn bam_for_each_parallel<P: AsRef<Path>, F>(path: P, threads: Option<usize>, on_record: F) -> anyhow::Result<usize>
where
    F: Fn(NARead) + Send + Sync,
{
    let mut reader = bam::Reader::from_path(path.as_ref())?;
    let n_threads = threads.unwrap_or_else(num_cpus::get);
    if n_threads > 1 {
        let _ = reader.set_threads(n_threads);
    }
    let counter = std::sync::atomic::AtomicUsize::new(0);
    let mut first_err: Option<anyhow::Error> = None;

    scope(|s| {
        for result in reader.records() {
            match result {
                Ok(rec) => {
                    let id = std::str::from_utf8(rec.qname()).unwrap_or("").to_string();
                    let seq = rec.seq().as_bytes();
                    let qual = if rec.qual().is_empty() { None } else { Some(rec.qual().to_vec()) };
                    let naread = NARead { id, seq, qual };
                    s.spawn(|_| {
                        on_record(naread);
                    });
                    counter.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                }
                Err(e) => {
                    first_err = Some(anyhow::anyhow!(e));
                    break;
                }
            }
        }
    });

    if let Some(e) = first_err { return Err(e); }
    Ok(counter.into_inner())
}

/// High-level entry point: dispatch by extension and parse with the appropriate reader.
pub fn for_each_parallel<P: AsRef<Path>, F>(path: P, threads: Option<usize>, on_record: F) -> anyhow::Result<(InputFormat, usize)>
where
    F: Fn(NARead) + Send + Sync,
{
    match detect_format(&path) {
        Some(InputFormat::Fastq) => Ok((InputFormat::Fastq, fastq_for_each_parallel(path, threads, on_record)?)),
        Some(InputFormat::Bam) | Some(InputFormat::Sam) => Ok((InputFormat::Bam, bam_for_each_parallel(path, threads, on_record)?)),
        None => Err(anyhow::anyhow!("Unsupported input format (expect .fastq(.gz)/.fq(.gz)/.sam/.bam)")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn fastq_counts_records() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, "@r1\nACGT\n+\nIIII").unwrap();
        writeln!(tmp, "@r2\nNNNN\n+\n!!!!").unwrap();
        let p = tmp.into_temp_path();
        let n = fastq_for_each_parallel(&p, Some(2), |_r| {}).unwrap();
        assert_eq!(n, 2);
    }
}
