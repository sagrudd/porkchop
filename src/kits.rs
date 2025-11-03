//! Registry mapping kit identifiers to their sequences.
//!
//! This includes current **Kit 14** families and selected legacy kits to help
//! interpret older datasets.
use crate::kit::{Kit, KitId};
use crate::BaseChemistry;
use crate::data::adapters::{RA_TOP, RTP, SSPII, CRTA};
use crate::data::cdna_legacy::{SSP, VNP};
use crate::data::{adapters::*, barcodes::*, legacy::*};

pub const KITS: &[Kit] = &[
    // Current ligation chemistry (Kit 14)
    Kit{
        id: KitId("LSK114"),
        description: "Ligation Sequencing Kit V14 (LSK114). Uses LA adapter; pairs with Native Barcoding kits.",
        adapters_and_primers: &[LA_TOP, LA_BOTTOM],
        chemistry: BaseChemistry::Rapid,
        legacy: false,
        barcodes: &[],
    },


// PCR‑cDNA Sequencing Kit (Kit 11)
Kit{
    id: KitId("PCS111"),
    description: "PCR‑cDNA Sequencing Kit (SQK‑PCS111). Uses legacy SSP/VNP primers and RA; CRTA+RTP included.",
    adapters_and_primers: &[RA_TOP, CRTA, RTP, SSP, VNP],
    chemistry: BaseChemistry::Rapid,
        legacy: false,
        barcodes: &[],
},

// PCR‑cDNA Sequencing Kit V14
Kit{
    id: KitId("PCS114"),
    description: "PCR‑cDNA Sequencing Kit V14 (SQK‑PCS114). Uses SSPII/RTP/CRTA and RA.",
    adapters_and_primers: &[RA_TOP, CRTA, RTP, SSPII],
    chemistry: BaseChemistry::Rapid,
        legacy: false,
        barcodes: &[],
},

    Kit{
        id: KitId("LSK114-XL"),
        description: "Ligation Sequencing Kit V14 XL (LSK114-XL). LA adapter; typical with NBD114.x sets.",
        adapters_and_primers: &[LA_TOP, LA_BOTTOM],
        chemistry: BaseChemistry::Rapid,
        legacy: false,
        barcodes: &[],
    },

    // Native barcoding (Kit 14)
    Kit{
        id: KitId("NBD114.24"),
        description: "Native Barcoding Kit 24 V14. Uses NA adapter + NB01–24; NB flanks.",
        adapters_and_primers: &[NA_TOP, NA_BOTTOM, NB_FLANK_FWD, NB_FLANK_REV5, NB_FLANK_REV3],
        chemistry: BaseChemistry::Rapid,
        legacy: false,
        barcodes: &NB_BARCODES_24,
    },
    Kit{
        id: KitId("NBD114.96"),
        description: "Native Barcoding Kit 96 V14. Uses NA adapter + NB01–96; NB flanks.",
        adapters_and_primers: &[NA_TOP, NA_BOTTOM, NB_FLANK_FWD, NB_FLANK_REV5, NB_FLANK_REV3],
        chemistry: BaseChemistry::Rapid,
        legacy: false,
        barcodes: NB_BARCODES,
    },

    // Rapid barcoding (Kit 14)
    Kit{
        id: KitId("RBK114.24"),
        description: "Rapid Barcoding Kit 24 V14. Uses RA adapter + RB01–24 core barcodes; RB flanks.",
        adapters_and_primers: &[RA_TOP, RB_FLANK_LEFT, RB_FLANK_RIGHT],
        chemistry: BaseChemistry::Rapid,
        legacy: false,
        barcodes: &SHARED_1_TO_24,
    },
    Kit{
        id: KitId("RBK114.96"),
        description: "Rapid Barcoding Kit 96 V14. Uses RA adapter + RB01–96 core barcodes; RB flanks.",
        adapters_and_primers: &[RA_TOP, RB_FLANK_LEFT, RB_FLANK_RIGHT],
        chemistry: BaseChemistry::Rapid,
        legacy: false,
        barcodes: &SHARED_BARCODE_SET,
    },

    // PCR‑cDNA barcoding (Kit 14)
    

// PCR‑cDNA barcoding (Kit 11)
Kit{
    id: KitId("PCB111.24"),
    description: "PCR‑cDNA Barcoding Kit 24 (Kit 11). Uses SSP/VNP (pychopper) with CRTA/RTP and PCB flanks; BP01–24.",
    adapters_and_primers: &[RA_TOP, RTP, CRTA, SSP, VNP, PCB_FLANK_TOP, PCB_FLANK_BOT_A, PCB_FLANK_BOT_B],
    chemistry: BaseChemistry::Rapid,
        legacy: false,
        barcodes: &SHARED_1_TO_24,
},
Kit{
        id: KitId("PCB114.24"),
        description: "PCR‑cDNA Barcoding Kit 24 V14. Uses cDNA adapters/primers + BP01–24; PCB flanks.",
        adapters_and_primers: &[RA_TOP, RTP, SSPII, CRTA, CPRM_FWD, CPRM_REV, PCB_FLANK_TOP, PCB_FLANK_BOT_A, PCB_FLANK_BOT_B],
        chemistry: BaseChemistry::Rapid,
        legacy: false,
        barcodes: &SHARED_1_TO_24,
    },

    // Rapid PCR barcoding (Kit 14)
    Kit{
        id: KitId("RPB114.24"),
        description: "Rapid PCR Barcoding Kit 24 V14. Uses RA + RLB01–24 core barcodes; RPB flank.",
        adapters_and_primers: &[RA_TOP, RPB_FLANK],
        chemistry: BaseChemistry::Rapid,
        legacy: false,
        barcodes: &SHARED_1_TO_24,
    },

    // 16S barcoding (Kit 14)
    Kit{
        id: 

// Microbial Amplicon Barcoding (16S & ITS) — Rapid‑based, 24 barcodes
// Provenance: ONT protocol page (SQK‑MAB114.24), Rapid workflow, up to 24 barcodes.

KitId("16S114.24"),
        description: "16S Barcoding Kit 24 V14. Uses RA + 16S01–24 core barcodes; 16S primer flanks.",
        adapters_and_primers: &[RA_TOP, SIXTEENS_FLANK, SIXTEENS_FWD_TARGET, SIXTEENS_REV_TARGET],
        chemistry: BaseChemistry::Rapid,
        legacy: false,
        barcodes: &SHARED_1_TO_24,
    },

    // Expansions (BC01–96)
    Kit{
        id: KitId("PBC001"),
        description: "PCR Barcoding Expansion 1–12 (EXP‑PBC001): BC01–BC12; PBC flanks.",
        adapters_and_primers: &[PCB_FLANK_TOP, PCB_FLANK_BOT_A, PCB_FLANK_BOT_B],
        chemistry: BaseChemistry::Rapid,
        legacy: false,
        barcodes: &SHARED_1_TO_12,
    },
    Kit{
        id: KitId("PBC096"),
        description: "PCR Barcoding Expansion 1–96 (EXP‑PBC096): BC01–BC96; PBC flanks.",
        adapters_and_primers: &[PCB_FLANK_TOP, PCB_FLANK_BOT_A, PCB_FLANK_BOT_B],
        chemistry: BaseChemistry::Rapid,
        legacy: false,
        barcodes: &SHARED_BARCODE_SET,
    },

    // Legacy rapid kits (for historic data)
    Kit{
        id: KitId("RBK004"),
        description: "Rapid Barcoding Kit (RBK004): RB01–RB12; RB flank.",
        adapters_and_primers: &[RA_TOP, RB_FLANK_LEFT, RB_FLANK_RIGHT],
        chemistry: BaseChemistry::Rapid,
        legacy: false,
        barcodes: &SHARED_1_TO_12,
    },
    Kit{
        id: KitId("RBK110.96"),
        description: "Rapid Barcoding Kit 96 (RBK110.96): RB01–RB96; RB flank.",
        adapters_and_primers: &[RA_TOP, RB_FLANK_LEFT, RB_FLANK_RIGHT],
        chemistry: BaseChemistry::Rapid,
        legacy: false,
        barcodes: &SHARED_BARCODE_SET,
    },

    // Legacy ligation adapters (informational)
    Kit{
        id: KitId("LSK109"),
        description: "Ligation Sequencing Kit LSK109 (legacy). Y‑adapter trunk per Porechop forks.",
        adapters_and_primers: &[NSK007_Y_TOP_TRUNK, NSK007_Y_BOTTOM],
        chemistry: BaseChemistry::Rapid,
        legacy: false,
        barcodes: &[],
    },
    Kit{
        id: KitId("LSK108"),
        description: "Ligation Sequencing Kit LSK108 (legacy). Y‑adapter trunk per Porechop forks.",
        adapters_and_primers: &[NSK007_Y_TOP_TRUNK, NSK007_Y_BOTTOM],
        chemistry: BaseChemistry::Rapid,
        legacy: false,
        barcodes: &[],
    },
    Kit{
        id: KitId("LSK308"),
        description: "1D^2 kit LSK308 (legacy). 1D^2 adapter fragments per Porechop forks.",
        adapters_and_primers: &[LSK308_1D2_TOP, LSK308_1D2_BOTTOM],
        chemistry: BaseChemistry::Rapid,
        legacy: false,
        barcodes: &[],
    },


// Microbial Amplicon Barcoding (16S & ITS) — Rapid‑based, 24 barcodes
// Provenance: ONT protocol page (SQK‑MAB114.24), Rapid workflow, up to 24 barcodes.
Kit{
    id: KitId("MAB114.24"),
    description: "Microbial Amplicon Barcoding 24 V14 (SQK‑MAB114.24). Rapid‑based; 16S and ITS targets; 24 barcodes.",
    adapters_and_primers: &[RA_TOP],
    chemistry: BaseChemistry::Rapid,
        legacy: false,
        barcodes: &SHARED_1_TO_24,
},
];
