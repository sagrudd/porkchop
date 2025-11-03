
//! Core types for kits, sequences and provenance.

/// Where a sequence definition came from.
#[derive(Debug, Clone, Copy)]
pub struct Provenance {
    pub source: &'static str,
    pub appendix: Option<&'static str>,
    pub notes: Option<&'static str>,
}

/// High-level category of sequence.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SeqKind {
    AdapterTop,
    AdapterBottom,
    Primer,
    Barcode,

    Flank,
}

/// A named nucleotide sequence with kind and provenance.
#[derive(Debug, Clone, Copy)]
pub struct SequenceRecord {
    pub name: &'static str,
    pub kind: SeqKind,
    pub sequence: &'static str,
    pub provenance: Provenance,
}

/// Newtype for kit identifiers (e.g., "LSK114", "PCS114", "NBD114.24").
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct KitId(pub &'static str);

/// Base sequencing chemistry.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BaseChemistry {
    Rapid,
    Ligation,
    PCRcDNA,
    Amplicon,
}

/// A kit bundles known adapters/primers and optional barcodes.
#[derive(Debug, Clone, Copy)]
pub struct Kit {
    pub id: KitId,
    pub description: &'static str,
    pub legacy: bool,
    pub chemistry: BaseChemistry,
    pub adapters_and_primers: &'static [SequenceRecord],
    pub barcodes: &'static [SequenceRecord],
}


impl std::fmt::Display for BaseChemistry {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            BaseChemistry::Rapid => "rapid sequencing",
            BaseChemistry::Ligation => "ligation sequencing",
            BaseChemistry::PCRcDNA => "pcr-cdna sequencing",
            BaseChemistry::Amplicon => "amplicon sequencing",
        };
        write!(f, "{}", s)
    }
}
