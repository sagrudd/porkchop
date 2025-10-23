//! # porkchop: Oxford Nanopore kit registry, IO and benchmarking
//!
//! This crate provides:
//!
//! - **Authoritative kit registry** (adapters, primers, barcodes; with provenance)
//! - **High‑performance sequence IO** for FASTQ/FASTQ.GZ/SAM/BAM
//! - **Benchmarking** of label assignment across several algorithms
//!
//! ## Minimum Supported Rust Version (MSRV)
//!
//! This crate targets **Rust 1.90** or newer (`rust-version = "1.90"` in `Cargo.toml`).
//!
//! ## CLI tools
//! See the README for `list-kits`, `describe-kit`, and `benchmark` usage examples.


pub mod kit;
pub mod kits;
pub mod seqio;
/// Benchmarking framework.
pub mod benchmark;


pub mod data { pub mod adapters; pub mod barcodes; pub mod cdna_legacy; pub mod legacy; }
pub use kit::{Kit, KitId, SeqKind, SequenceRecord, BaseChemistry};

/// Return static registry of supported kits.
pub fn list_supported_kits() -> &'static [kit::Kit] { kits::KITS }

/// Lookup a kit by id (case-sensitive).
pub fn get_sequences_for_kit(id: &str) -> Option<&'static kit::Kit> {
    kits::KITS.iter().find(|k| k.id.0 == id)
}

/// Is a kit legacy?
pub fn kit_is_legacy(k: &kit::Kit) -> bool { k.legacy }

/// Base chemistry for a kit.
pub fn base_chemistry_of(k: &kit::Kit) -> kit::BaseChemistry { k.chemistry }
