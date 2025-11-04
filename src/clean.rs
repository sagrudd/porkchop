
use std::path::{Path, PathBuf};

use rayon::prelude::*;

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
    w.write_all(b"\n")?;
    w.write_all(seq)?;
    w.write_all(b"\n+\n")?;
    w.write_all(qual)?;
    w.write_all(b"\n")?;
    Ok(())
}

mod edwrap {
    use edlib_rs::edlibrs::{edlibAlignRs, EdlibAlignConfigRs, EdlibAlignModeRs, EdlibAlignTaskRs, EdlibEqualityPairRs};
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
    _kind: &'a str,
    seq: &'a [u8],
}

fn motifs_for_kit<'a>(kit: &'a crate::kit::Kit) -> Vec<Motif<'a>> {
    let mut m = Vec::new();
    for s in kit.adapters_and_primers {
        m.push(Motif { name: s.name, _kind: "adapter_or_primer", seq: s.sequence.as_bytes() });
    }
    for s in kit.barcodes {
        m.push(Motif { name: s.name, _kind: "barcode_or_flank", seq: s.sequence.as_bytes() });
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
    let mut left_best: Option<(i32, i32, i32, &str)> = None;
    let mut right_best: Option<(i32, i32, i32, &str)> = None;

    for m in motifs {
        let maxk = max_edits_for(m.seq.len()) as i32;
        if let Some(hit) = edwrap::locate(m.seq, &s, maxk) {
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

fn process_fastx_to_gz(out_path: &Path, input_files: Vec<PathBuf>, _threads_eff: usize, kit_id: &str) -> anyhow::Result<()> {
    use std::fs::File;
    use std::io::BufWriter;
    use needletail::parser::parse_fastx_file;

    let kit = crate::get_sequences_for_kit(kit_id).expect("validated kit");
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
                        let qname = std::str::from_utf8(r.qname()).unwrap_or("SAM");
                        let seq = r.seq().as_bytes();
                        let qualv = r.qual().to_vec();
                        let qual = qualv.into_iter().map(|q| (q + 33) as u8).collect::<Vec<u8>>();
                        let trimmed = annotate_and_trim_one(&seq, &qual, kit_id, &motifs);
                        OwnedRecord { id: format!("{} {}", qname, trimmed.id), seq: trimmed.seq, qual: trimmed.qual }
                    }).collect();
                    for pr in &processed { write_fastq_record(&mut gz, &pr.id, &pr.seq, &pr.qual)?; }
                    buf.clear();
                }
            }
            if !buf.is_empty() {
                let processed: Vec<OwnedRecord> = buf.par_iter().map(|r| {
                    let qname = std::str::from_utf8(r.qname()).unwrap_or("SAM");
                    let seq = r.seq().as_bytes();
                    let qualv = r.qual().to_vec();
                    let qual = qualv.into_iter().map(|q| (q + 33) as u8).collect::<Vec<u8>>();
                    let trimmed = annotate_and_trim_one(&seq, &qual, kit_id, &motifs);
                    OwnedRecord { id: format!("{} {}", qname, trimmed.id), seq: trimmed.seq, qual: trimmed.qual }
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
                        let qname = std::str::from_utf8(r.qname()).unwrap_or("BAM");
                        let seq = r.seq().as_bytes();
                        let qualv = r.qual().to_vec();
                        let qual = qualv.into_iter().map(|q| (q + 33) as u8).collect::<Vec<u8>>();
                        let trimmed = annotate_and_trim_one(&seq, &qual, kit_id, &motifs);
                        OwnedRecord { id: format!("{} {}", qname, trimmed.id), seq: trimmed.seq, qual: trimmed.qual }
                    }).collect();
                    for pr in &processed { write_fastq_record(&mut gz, &pr.id, &pr.seq, &pr.qual)?; }
                    buf.clear();
                }
            }
            if !buf.is_empty() {
                let processed: Vec<OwnedRecord> = buf.par_iter().map(|r| {
                    let qname = std::str::from_utf8(r.qname()).unwrap_or("BAM");
                    let seq = r.seq().as_bytes();
                    let qualv = r.qual().to_vec();
                    let qual = qualv.into_iter().map(|q| (q + 33) as u8).collect::<Vec<u8>>();
                    let trimmed = annotate_and_trim_one(&seq, &qual, kit_id, &motifs);
                    OwnedRecord { id: format!("{} {}", qname, trimmed.id), seq: trimmed.seq, qual: trimmed.qual }
                }).collect();
                for pr in &processed { write_fastq_record(&mut gz, &pr.id, &pr.seq, &pr.qual)?; }
            }

        } else {
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

                let processed: Vec<OwnedRecord> = owned_chunk.par_iter()
                    .map(|r| annotate_and_trim_one(&r.seq, &r.qual, kit_id, &motifs))
                    .collect();

                for pr in &processed { write_fastq_record(&mut gz, &pr.id, &pr.seq, &pr.qual)?; }
            }
        }
    }

    gz.finish()?;
    Ok(())
}

pub fn run(threads: usize, kit: &str, output: &Path, files: Vec<PathBuf>) -> anyhow::Result<()> {
    ensure_known_kit(kit)?;

    let (ok, bad) = split_supported_files(files);
    if !bad.is_empty() {
        let mut msg = String::from("Unsupported file type(s):\n");
        for p in &bad { msg.push_str(&format!("  - {}\n", p.display())); }
        msg.push_str("Allowed: SAM (.sam), BAM (.bam), FASTQ (.fastq/.fq), and gzipped FASTQ (.fastq.gz/.fq.gz).");
        anyhow::bail!(msg);
    }

    let threads_eff = if threads == 0 { std::cmp::max(1, num_cpus::get()) } else { threads };
    rayon::ThreadPoolBuilder::new().num_threads(threads_eff).build_global().ok();

    eprintln!("clean: kit={} | threads={} | inputs={} | output={}", kit, threads_eff, ok.len(), output.display());
    process_fastx_to_gz(output, ok, threads_eff, kit)
}
