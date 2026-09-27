"""Titans Med genomics/protein identity spine.

This is the shared, research-planning-only identity catalog for genes, proteins,
processes, and tissues. Program overlays such as the Mito Evidence Atlas consume
a projection of this catalog; they do not own it. The catalog never records
claim status or eligibility.
"""
from __future__ import annotations

import copy
import hashlib
import html
import json
import re
from collections.abc import Iterable, Mapping
from typing import Any, NamedTuple

import networkx as nx
import plotly.graph_objects as go

ACCESSED_DATE = "2026-08-30"
BUNDLE_SCHEMA = "titans-med-genomics-protein-spine/v2"
HUMAN_TAXONOMY = "NCBITaxon:9606"
NON_SPECIES_SCOPE = "non-species"
ALLOWED_KINDS = frozenset({"gene", "transcript", "variant", "haplotype", "protein", "complex", "process", "tissue"})
ALLOWED_LOCI = frozenset({"nuclear", "mitochondrial"})
ALLOWED_PROGRAMS = frozenset({"mito", "neuro", "pgx", "cardio", "onco"})
ALLOWED_EPISTEMIC = frozenset({"identity", "hypothesis", "context", "assay_context"})
ALLOWED_VARIANT_CLASSES = frozenset({"snv", "haplotype"})
MITO_ATLAS_KINDS = frozenset({"gene", "protein", "process", "tissue"})
_PDB_RE = re.compile(r"^[0-9][A-Z0-9]{3}$")
_HASH_RE = re.compile(r"^[0-9a-f]{64}$")
_HGVS_RE = re.compile(r"^(NC_\d{6}\.\d+):g\.([0-9]+)([ACGT])>([ACGT])$")
_ACCESSION_RE = re.compile(r"^(NC_\d{6})\.(\d+)$")
_NM_RE = re.compile(r"^NM_\d+\.\d+$")
_ENST_RE = re.compile(r"^ENST\d{11}\.\d+$")
_DBSNP_RE = re.compile(r"^dbSNP:rs\d+$")
_CDNA_HGVS_RE = re.compile(r"^(NM_\d+\.\d+):c\.[0-9]+[ACGT]>[ACGT]$")
_PROTEIN_HGVS_RE = re.compile(r"^NP_\d+\.\d+:p\.[A-Za-z]+[0-9]+([A-Za-z]+|=)$")
_FORBIDDEN_KEYS = frozenset(
    {"raw_text", "full_text", "abstract", "claim_status", "eligibility", "automatic_status_change"}
)
_ALPHAFOLD_RE = re.compile(r"^AF-[A-Z0-9]{6,10}-F1$")

EDGE_STYLES = {
    "identity": {"color": "#8aa6a8", "dash": "solid", "opacity": 0.72},
    "context": {"color": "#6f9f9a", "dash": "dash", "opacity": 0.58},
    "assay_context": {"color": "#708792", "dash": "longdash", "opacity": 0.52},
    "hypothesis": {"color": "#526873", "dash": "dot", "opacity": 0.34},
}
KIND_COLORS = {
    "compound": "#55c5d9",
    "project": "#e8e6d9",
    "gene": "#e8b94a",
    "protein": "#9b82e8",
    "transcript": "#d7c07a",
    "variant": "#d98b6a",
    "haplotype": "#c47c8a",
    "complex": "#c4a35a",
    "process": "#86c98d",
    "tissue": "#e77d82",
}


class GeneSpec(NamedTuple):
    symbol: str
    hgnc: str
    ncbi_gene: str
    ensembl: str
    uniprot: str
    locus: str
    programs: tuple[str, ...]
    process_ids: tuple[str, ...]


class ProcessSpec(NamedTuple):
    go_id: str
    label: str
    programs: tuple[str, ...]


class TissueSpec(NamedTuple):
    uberon_id: str
    label: str
    programs: tuple[str, ...]


class ComplexSpec(NamedTuple):
    complexportal: str
    label: str
    go_id: str
    programs: tuple[str, ...]
    subunit_uniprots: tuple[str, ...]
    locus: str = "mitochondrial"


class VariantSpec(NamedTuple):
    variant_id: str
    refseq_accession: str
    hgvs_genomic: str
    source_identifiers: tuple[str, ...]
    review_status: str
    gene_hgnc: str
    transcript_nm: str = ""
    hgvs_cdna: str = ""
    hgvs_protein: str = ""


class HaplotypeSpec(NamedTuple):
    haplotype_id: str
    gene_hgnc: str
    label: str
    composition_of: tuple[str, ...]
    source_identifiers: tuple[str, ...]
    review_status: str


class TranscriptSpec(NamedTuple):
    gene_symbol: str
    mane_select: str
    ensembl_transcript: str


class FeatureSpec(NamedTuple):
    feature_type: str
    start: int
    end: int
    label: str


NUCLEAR_MITO_GENES: tuple[GeneSpec, ...] = (
    GeneSpec("PPARGC1A", "HGNC:9237", "10891", "ENSG00000109819", "Q9UBK2", "nuclear", ("mito",), ("GO:0007005", "GO:0000422", "GO:0022900")),
    GeneSpec("CREB1", "HGNC:2345", "1385", "ENSG00000118260", "P16220", "nuclear", ("mito",), ("GO:0006357",)),
    GeneSpec("PINK1", "HGNC:14581", "65018", "ENSG00000158828", "Q9BXM7", "nuclear", ("mito",), ("GO:0007005", "GO:0000422")),
    GeneSpec("PRKN", "HGNC:8607", "5071", "ENSG00000185345", "O60260", "nuclear", ("mito",), ("GO:0007005", "GO:0000422")),
    GeneSpec("MFSD2A", "HGNC:25897", "84879", "ENSG00000168389", "Q8NA29", "nuclear", ("mito",), ("GO:0015908",)),
    GeneSpec("SIRT1", "HGNC:14929", "23411", "ENSG00000096717", "Q96EB6", "nuclear", ("mito",), ("GO:0007005", "GO:0000422")),
    GeneSpec("NRF1", "HGNC:7996", "4899", "ENSG00000106459", "Q16656", "nuclear", ("mito",), ("GO:0007005", "GO:0022900")),
    GeneSpec("TFAM", "HGNC:11741", "7019", "ENSG00000108064", "Q00059", "nuclear", ("mito",), ("GO:0007005", "GO:0022900")),
)

# Official HGNC / NCBI Gene / Ensembl / UniProt identities for the 13 human
# mtDNA protein-coding genes (NC_012920.1 / GRCh38).
MTDNA_GENES: tuple[GeneSpec, ...] = (
    GeneSpec("MT-ND1", "HGNC:7455", "4535", "ENSG00000198888", "P03886", "mitochondrial", ("mito",), ("GO:0007005", "GO:0022900")),
    GeneSpec("MT-ND2", "HGNC:7456", "4536", "ENSG00000198763", "P03891", "mitochondrial", ("mito",), ("GO:0007005", "GO:0022900")),
    GeneSpec("MT-ND3", "HGNC:7458", "4537", "ENSG00000198840", "P03897", "mitochondrial", ("mito",), ("GO:0007005", "GO:0022900")),
    GeneSpec("MT-ND4", "HGNC:7459", "4538", "ENSG00000198886", "P03905", "mitochondrial", ("mito",), ("GO:0007005", "GO:0022900")),
    GeneSpec("MT-ND4L", "HGNC:7460", "4539", "ENSG00000212907", "P03901", "mitochondrial", ("mito",), ("GO:0007005", "GO:0022900")),
    GeneSpec("MT-ND5", "HGNC:7461", "4540", "ENSG00000198786", "P03915", "mitochondrial", ("mito",), ("GO:0007005", "GO:0022900")),
    GeneSpec("MT-ND6", "HGNC:7462", "4541", "ENSG00000198695", "P03923", "mitochondrial", ("mito",), ("GO:0007005", "GO:0022900")),
    GeneSpec("MT-CO1", "HGNC:7419", "4512", "ENSG00000198804", "P00395", "mitochondrial", ("mito",), ("GO:0007005", "GO:0022900")),
    GeneSpec("MT-CO2", "HGNC:7421", "4513", "ENSG00000198712", "P00403", "mitochondrial", ("mito",), ("GO:0007005", "GO:0022900")),
    GeneSpec("MT-CO3", "HGNC:7422", "4514", "ENSG00000198938", "P00414", "mitochondrial", ("mito",), ("GO:0007005", "GO:0022900")),
    GeneSpec("MT-ATP6", "HGNC:7414", "4508", "ENSG00000198899", "P00846", "mitochondrial", ("mito",), ("GO:0007005", "GO:0022900")),
    GeneSpec("MT-ATP8", "HGNC:7415", "4509", "ENSG00000228253", "P03928", "mitochondrial", ("mito",), ("GO:0007005", "GO:0022900")),
    GeneSpec("MT-CYB", "HGNC:7427", "4519", "ENSG00000198727", "P00156", "mitochondrial", ("mito",), ("GO:0007005", "GO:0022900")),
)

# Official HGNC / NCBI Gene / Ensembl / UniProt identities for neuro, PGx,
# cardio, and onco program genes. Representative PDBs are UniProt PDB
# cross-references fetched 2026-08-30.
NUCLEAR_PROGRAM_GENES: tuple[GeneSpec, ...] = (
    GeneSpec("APOE", "HGNC:613", "348", "ENSG00000130203", "P02649", "nuclear", ("neuro",), ("GO:0007611", "GO:0030301")),
    GeneSpec("APP", "HGNC:620", "351", "ENSG00000142192", "P05067", "nuclear", ("neuro",), ("GO:0007611",)),
    GeneSpec("MAPT", "HGNC:6893", "4137", "ENSG00000186868", "P10636", "nuclear", ("neuro",), ("GO:0007611",)),
    GeneSpec("SNCA", "HGNC:11138", "6622", "ENSG00000145335", "P37840", "nuclear", ("neuro",), ("GO:0007611",)),
    GeneSpec("TREM2", "HGNC:17761", "54209", "ENSG00000095970", "Q9NZC2", "nuclear", ("neuro",), ("GO:0007611",)),
    GeneSpec("GBA1", "HGNC:4177", "2629", "ENSG00000177628", "P04062", "nuclear", ("neuro",), ("GO:0006665", "GO:0007611")),
    GeneSpec("CYP2D6", "HGNC:2625", "1565", "ENSG00000100197", "P10635", "nuclear", ("pgx",), ("GO:0006805",)),
    GeneSpec("CYP2C19", "HGNC:2621", "1557", "ENSG00000165841", "P33261", "nuclear", ("pgx",), ("GO:0006805",)),
    GeneSpec("BRCA1", "HGNC:1100", "672", "ENSG00000012048", "P38398", "nuclear", ("onco",), ("GO:0006281", "GO:0006974")),
    GeneSpec("TP53", "HGNC:11998", "7157", "ENSG00000141510", "P04637", "nuclear", ("onco",), ("GO:0006974", "GO:0006915")),
    GeneSpec("PCSK9", "HGNC:20001", "255738", "ENSG00000169174", "Q8NBP7", "nuclear", ("cardio",), ("GO:0030301",)),
    GeneSpec("LDLR", "HGNC:6547", "3949", "ENSG00000130164", "P01130", "nuclear", ("cardio",), ("GO:0030301",)),
    GeneSpec("MYH7", "HGNC:7577", "4625", "ENSG00000092054", "P12883", "nuclear", ("cardio",), ("GO:0006936",)),
)

GENE_CATALOG: tuple[GeneSpec, ...] = NUCLEAR_MITO_GENES + MTDNA_GENES + NUCLEAR_PROGRAM_GENES

# GRCh38 chromosome accessions from NCBI Gene genomicinfo (2026-08-30).
NUCLEAR_CHROMOSOME_ACCESSION = {
    "PPARGC1A": "NC_000004.12",
    "CREB1": "NC_000002.12",
    "PINK1": "NC_000001.11",
    "PRKN": "NC_000006.12",
    "MFSD2A": "NC_000001.11",
    "SIRT1": "NC_000010.11",
    "NRF1": "NC_000007.14",
    "TFAM": "NC_000010.11",
    "APOE": "NC_000019.10",
    "APP": "NC_000021.9",
    "MAPT": "NC_000017.11",
    "SNCA": "NC_000004.12",
    "TREM2": "NC_000006.12",
    "GBA1": "NC_000001.11",
    "CYP2D6": "NC_000022.11",
    "CYP2C19": "NC_000010.11",
    "BRCA1": "NC_000017.11",
    "TP53": "NC_000017.11",
    "PCSK9": "NC_000001.11",
    "LDLR": "NC_000019.10",
    "MYH7": "NC_000014.9",
}

# HGNC MANE Select NM_ and ENST accessions fetched 2026-08-30. Mitochondrial
# protein-coding genes have no MANE Select transcript.
TRANSCRIPT_CATALOG: tuple[TranscriptSpec, ...] = (
    TranscriptSpec("PPARGC1A", "NM_013261.5", "ENST00000264867.7"),
    TranscriptSpec("CREB1", "NM_004379.5", "ENST00000353267.8"),
    TranscriptSpec("PINK1", "NM_032409.3", "ENST00000321556.5"),
    TranscriptSpec("PRKN", "NM_004562.3", "ENST00000366898.6"),
    TranscriptSpec("MFSD2A", "NM_032793.5", "ENST00000372811.10"),
    TranscriptSpec("SIRT1", "NM_012238.5", "ENST00000212015.11"),
    TranscriptSpec("NRF1", "NM_005011.5", "ENST00000393232.6"),
    TranscriptSpec("TFAM", "NM_003201.3", "ENST00000487519.6"),
    TranscriptSpec("APOE", "NM_000041.4", "ENST00000252486.9"),
    TranscriptSpec("APP", "NM_000484.4", "ENST00000346798.8"),
    TranscriptSpec("MAPT", "NM_001377265.1", "ENST00000262410.10"),
    TranscriptSpec("SNCA", "NM_000345.4", "ENST00000394991.8"),
    TranscriptSpec("TREM2", "NM_018965.4", "ENST00000373113.8"),
    TranscriptSpec("GBA1", "NM_000157.4", "ENST00000368373.8"),
    TranscriptSpec("CYP2D6", "NM_000106.6", "ENST00000645361.2"),
    TranscriptSpec("CYP2C19", "NM_000769.4", "ENST00000371321.9"),
    TranscriptSpec("BRCA1", "NM_007294.4", "ENST00000357654.9"),
    TranscriptSpec("TP53", "NM_000546.6", "ENST00000269305.9"),
    TranscriptSpec("PCSK9", "NM_174936.4", "ENST00000302118.5"),
    TranscriptSpec("LDLR", "NM_000527.5", "ENST00000558518.6"),
    TranscriptSpec("MYH7", "NM_000257.4", "ENST00000355349.4"),
)

# UniProt feature ranges fetched 2026-08-30. These are architecture context,
# not functional claims. Representative PDB coverage is residue-range identity
# of the experimental construct, not full-chain coverage.
PROTEIN_ARCHITECTURE = {
    "Q9UBK2": (FeatureSpec("Domain", 677, 753, "RRM"),),
    "P16220": (FeatureSpec("Domain", 87, 146, "KID"), FeatureSpec("Domain", 269, 327, "bZIP")),
    "Q9BXM7": (FeatureSpec("Domain", 156, 511, "Protein kinase"),),
    "O60260": (
        FeatureSpec("Domain", 1, 76, "Ubiquitin-like"),
        FeatureSpec("Zinc finger", 141, 225, "RING-type 0"),
        FeatureSpec("Zinc finger", 238, 293, "RING-type 1"),
        FeatureSpec("Zinc finger", 313, 377, "IBR-type"),
    ),
    "Q8NA29": (FeatureSpec("Transmembrane", 41, 70, "Helical"),),
    "Q96EB6": (FeatureSpec("Domain", 236, 496, "Deacetylase sirtuin-type"),),
    "Q16656": (FeatureSpec("DNA binding", 109, 305, "DNA binding"),),
    "Q00059": (
        FeatureSpec("DNA binding", 50, 118, "HMG box 1"),
        FeatureSpec("DNA binding", 155, 219, "HMG box 2"),
    ),
    "P02649": (
        FeatureSpec("Repeat", 80, 101, "1"),
        FeatureSpec("Repeat", 102, 123, "2"),
        FeatureSpec("Repeat", 124, 145, "3"),
        FeatureSpec("Repeat", 146, 167, "4"),
    ),
    "P05067": (FeatureSpec("Domain", 28, 189, "E1"), FeatureSpec("Transmembrane", 702, 722, "Helical")),
    "P10636": (
        FeatureSpec("Repeat", 561, 591, "Tau/MAP 1"),
        FeatureSpec("Repeat", 592, 622, "Tau/MAP 2"),
        FeatureSpec("Repeat", 623, 653, "Tau/MAP 3"),
        FeatureSpec("Repeat", 654, 685, "Tau/MAP 4"),
    ),
    "P37840": (FeatureSpec("Repeat", 20, 30, "1"), FeatureSpec("Repeat", 31, 41, "2")),
    "Q9NZC2": (FeatureSpec("Domain", 29, 112, "Ig-like V-type"), FeatureSpec("Transmembrane", 175, 195, "Helical")),
    "P38398": (
        FeatureSpec("Zinc finger", 24, 65, "RING-type"),
        FeatureSpec("Domain", 1642, 1736, "BRCT 1"),
        FeatureSpec("Domain", 1756, 1855, "BRCT 2"),
    ),
    "P04637": (FeatureSpec("DNA binding", 102, 292, "DNA binding"),),
    "Q8NBP7": (
        FeatureSpec("Domain", 77, 149, "Inhibitor I9"),
        FeatureSpec("Domain", 155, 461, "Peptidase S8"),
        FeatureSpec("Region", 450, 692, "C-terminal domain"),
    ),
    "P01130": (FeatureSpec("Domain", 25, 65, "LDL-receptor class A 1"), FeatureSpec("Transmembrane", 789, 810, "Helical")),
    "P12883": (FeatureSpec("Domain", 85, 778, "Myosin motor"), FeatureSpec("Domain", 781, 810, "IQ")),
    "P04062": (FeatureSpec("Chain", 40, 536, "Lysosomal acid glucosylceramidase"),),
    "P10635": (FeatureSpec("Chain", 1, 497, "Cytochrome P450 2D6"),),
    "P33261": (FeatureSpec("Chain", 1, 490, "Cytochrome P450 2C19"),),
    "P00846": (FeatureSpec("Transmembrane", 6, 26, "Helical"), FeatureSpec("Transmembrane", 68, 88, "Helical")),
    "P03905": (FeatureSpec("Transmembrane", 22, 42, "Helical"), FeatureSpec("Transmembrane", 60, 80, "Helical")),
}

PROCESS_CATALOG: tuple[ProcessSpec, ...] = (
    ProcessSpec("GO:0007005", "mitochondrion organization", ("mito",)),
    ProcessSpec("GO:0000422", "autophagy of mitochondrion", ("mito",)),
    ProcessSpec("GO:0022900", "electron transport chain", ("mito",)),
    ProcessSpec("GO:0015908", "fatty acid transport", ("mito",)),
    ProcessSpec("GO:0006357", "regulation of transcription by RNA polymerase II", ("mito",)),
    ProcessSpec("GO:0007611", "learning or memory", ("neuro",)),
    ProcessSpec("GO:0006665", "sphingolipid metabolic process", ("neuro",)),
    ProcessSpec("GO:0006805", "xenobiotic metabolic process", ("pgx",)),
    ProcessSpec("GO:0030301", "cholesterol transport", ("cardio", "neuro")),
    ProcessSpec("GO:0006936", "muscle contraction", ("cardio",)),
    ProcessSpec("GO:0006281", "DNA repair", ("onco",)),
    ProcessSpec("GO:0006974", "cellular response to DNA damage stimulus", ("onco",)),
    ProcessSpec("GO:0006915", "apoptotic process", ("onco",)),
)

TISSUE_CATALOG: tuple[TissueSpec, ...] = (
    TissueSpec("UBERON:0001134", "skeletal muscle tissue", ("mito", "cardio")),
    TissueSpec("UBERON:0002113", "kidney", ("mito",)),
    TissueSpec("UBERON:0000955", "brain", ("mito", "neuro")),
    TissueSpec("UBERON:0000120", "blood brain barrier", ("mito", "neuro")),
    TissueSpec("UBERON:0002107", "liver", ("pgx", "cardio")),
    TissueSpec("UBERON:0001133", "cardiac muscle tissue", ("cardio",)),
    TissueSpec("UBERON:0002108", "small intestine", ("pgx",)),
)

ASSAY_CONTEXT_MATRIX: dict[str, tuple[str, ...]] = {
    "GO:0007005": ("UBERON:0001134", "UBERON:0002113", "UBERON:0000955"),
    "GO:0000422": ("UBERON:0001134", "UBERON:0000955"),
    "GO:0022900": ("UBERON:0001134", "UBERON:0002113", "UBERON:0000955"),
    "GO:0015908": ("UBERON:0000955", "UBERON:0000120"),
    "GO:0007611": ("UBERON:0000955", "UBERON:0000120"),
    "GO:0006665": ("UBERON:0000955",),
    "GO:0006805": ("UBERON:0002107", "UBERON:0002108"),
    "GO:0030301": ("UBERON:0002107", "UBERON:0000955"),
    "GO:0006936": ("UBERON:0001133", "UBERON:0001134"),
    "GO:0006281": ("UBERON:0002107",),
    "GO:0006974": ("UBERON:0002107",),
    "GO:0006915": ("UBERON:0002107",),
}

MITO_ATLAS_GENE_SYMBOLS = frozenset(spec.symbol for spec in NUCLEAR_MITO_GENES)

# Representative experimental structures are UniProt PDB cross-references whose
# RCSB titles identify a human structure of that protein or complex. AlphaFold
# remains the complete-chain identity for every protein.
REPRESENTATIVE_PDB = {
    "Q9UBK2": "1XB7",  # PGC-1α peptide / ERRα LBD
    "P16220": "5ZK1",  # CRTC2-CREB-CRE
    "Q9BXM7": "9EIH",  # PINK1-TOM import complex
    "O60260": "5C1Z",  # Parkin UblR0RBR
    "Q8NA29": "7OIX",  # human MFSD2A
    "Q96EB6": "4I5I",  # SIRT1 catalytic domain
    "Q16656": "8K3D",  # NRF1 DNA-binding domain
    "Q00059": "3FGH",  # human TFAM
    "P03886": "5XTC",
    "P03891": "5XTC",
    "P03897": "5XTC",
    "P03905": "5XTC",
    "P03901": "5XTC",
    "P03915": "5XTC",
    "P03923": "5XTC",
    "P00395": "5Z62",
    "P00403": "5Z62",
    "P00414": "5Z62",
    "P00846": "8H9F",
    "P03928": "8H9F",
    "P00156": "5XTE",
    "P02649": "1LPE",  # human apolipoprotein E
    "P05067": "1AAP",  # Alzheimer amyloid precursor protease inhibitor domain
    "P10636": "2MZ7",  # tau repeat domain
    "P37840": "1XQ8",  # human micelle-bound alpha-synuclein
    "Q9NZC2": "5ELI",  # TREM2 Ig-like domain
    "P04062": "1OGS",  # acid-beta-glucosidase / GBA1
    "P10635": "2F9Q",  # human CYP2D6
    "P33261": "4GQS",  # human CYP2C19
    "P38398": "1JNX",  # BRCA1 BRCT
    "P04637": "1C26",  # p53 tetramerization domain
    "Q8NBP7": "2P4E",  # PCSK9
    "P01130": "1AJJ",  # LDL receptor ligand-binding module
    "P12883": "2FXM",  # cardiac myosin S2 fragment
}

COMPLEX_CATALOG: tuple[ComplexSpec, ...] = (
    ComplexSpec("CPX-577", "respiratory chain complex I", "GO:0045271", ("mito",), ("P03886", "P03891", "P03897", "P03905", "P03901", "P03915", "P03923")),
    ComplexSpec("CPX-560", "respiratory chain complex III", "GO:0045275", ("mito",), ("P00156",)),
    ComplexSpec("CPX-6123", "respiratory chain complex IV", "GO:0045277", ("mito",), ("P00395", "P00403", "P00414")),
    ComplexSpec("CPX-6151", "proton-transporting ATP synthase complex", "GO:0045259", ("mito",), ("P00846", "P03928")),
    ComplexSpec("CPX-128", "LDLR-PCSK9 complex", "", ("cardio",), ("P01130", "Q8NBP7"), "nuclear"),
)

COMPLEX_PROCESS_EDGES = {
    "CPX-577": ("process:GO:0022900", "identity"),
    "CPX-560": ("process:GO:0022900", "identity"),
    "CPX-6123": ("process:GO:0022900", "identity"),
    "CPX-6151": ("process:GO:0007005", "hypothesis"),
    "CPX-128": ("process:GO:0030301", "identity"),
}

# Curated GRCh38 / NC_012920.1 SNVs. Genomic HGVS accessions must match the
# versioned RefSeq accession. Source identifiers are dbSNP rsIDs from Ensembl
# variation mappings fetched 2026-08-30. These records are identity catalog
# entries, not pathogenicity or eligibility claims.
SEEDED_VARIANT_SPECS: tuple[VariantSpec, ...] = (
    VariantSpec("VAR-MT-ATP6-m8993TG", "NC_012920.1", "NC_012920.1:g.8993T>G", ("dbSNP:rs199476133",), "curated", "HGNC:7414"),
    VariantSpec("VAR-MT-ND4-m11778GA", "NC_012920.1", "NC_012920.1:g.11778G>A", ("dbSNP:rs199476112",), "curated", "HGNC:7459"),
    VariantSpec("VAR-APOE-rs429358", "NC_000019.10", "NC_000019.10:g.44908684T>C", ("dbSNP:rs429358",), "curated", "HGNC:613", "NM_000041.4", "NM_000041.4:c.388T>C", "NP_000032.1:p.Cys130Arg"),
    VariantSpec("VAR-APOE-rs7412", "NC_000019.10", "NC_000019.10:g.44908822C>T", ("dbSNP:rs7412",), "curated", "HGNC:613", "NM_000041.4", "NM_000041.4:c.526C>T", "NP_000032.1:p.Arg176Cys"),
    VariantSpec("VAR-CYP2C19-rs4244285", "NC_000010.11", "NC_000010.11:g.94781859G>A", ("dbSNP:rs4244285",), "curated", "HGNC:2621", "NM_000769.4", "NM_000769.4:c.681G>A", "NP_000760.1:p.Pro227="),
    VariantSpec("VAR-PCSK9-rs11591147", "NC_000001.11", "NC_000001.11:g.55039974G>T", ("dbSNP:rs11591147",), "curated", "HGNC:20001", "NM_174936.4"),
    VariantSpec("VAR-MYH7-rs121913624", "NC_000014.9", "NC_000014.9:g.23429278C>T", ("dbSNP:rs121913624",), "curated", "HGNC:7577", "NM_000257.4"),
    VariantSpec("VAR-TP53-rs28934578", "NC_000017.11", "NC_000017.11:g.7675088C>T", ("dbSNP:rs28934578",), "curated", "HGNC:11998", "NM_000546.6"),
    VariantSpec("VAR-BRCA1-rs28897672", "NC_000017.11", "NC_000017.11:g.43106487A>C", ("dbSNP:rs28897672",), "curated", "HGNC:1100", "NM_007294.4"),
)

# Haplotype / star-allele composition objects. These name genomic-SNV
# combinations; they do not encode function, metabolizer status, or risk.
SEEDED_HAPLOTYPE_SPECS: tuple[HaplotypeSpec, ...] = (
    HaplotypeSpec("HAP-APOE-e3", "HGNC:613", "APOE epsilon3 genomic composition", ("VAR-APOE-rs429358:T", "VAR-APOE-rs7412:C"), ("dbSNP:rs429358", "dbSNP:rs7412"), "curated"),
    HaplotypeSpec("HAP-APOE-e4", "HGNC:613", "APOE epsilon4 genomic composition", ("VAR-APOE-rs429358:C", "VAR-APOE-rs7412:C"), ("dbSNP:rs429358", "dbSNP:rs7412"), "curated"),
    HaplotypeSpec("HAP-APOE-e2", "HGNC:613", "APOE epsilon2 genomic composition", ("VAR-APOE-rs429358:T", "VAR-APOE-rs7412:T"), ("dbSNP:rs429358", "dbSNP:rs7412"), "curated"),
    HaplotypeSpec("HAP-CYP2C19-star2-defining", "HGNC:2621", "CYP2C19 star2 defining SNV composition", ("VAR-CYP2C19-rs4244285:A",), ("dbSNP:rs4244285",), "curated"),
)

COMPOUND_GENE_SYMBOLS = {
    "pqq": ("PPARGC1A", "CREB1"),
    "urolithin_a": ("PINK1", "PRKN"),
    "dha": ("MFSD2A",),
}
COMPOUND_PROCESS_IDS = {
    "coq10": ("GO:0022900",),
}
COMPOUND_COMPLEX_IDS = {
    "coq10": tuple(spec.complexportal for spec in COMPLEX_CATALOG if "mito" in spec.programs),
}


def official_url(kind: str, value: str) -> str:
    urls = {
        "hgnc": f"https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/{value}",
        "ncbi": f"https://www.ncbi.nlm.nih.gov/gene/{value}",
        "ensembl": f"https://www.ensembl.org/id/{value}",
        "uniprot": f"https://www.uniprot.org/uniprotkb/{value}/entry",
        "alphafold": f"https://alphafold.ebi.ac.uk/entry/{value}",
        "go": f"https://www.ebi.ac.uk/QuickGO/term/{value}",
        "uberon": (
            "https://www.ebi.ac.uk/ols4/ontologies/uberon/classes/"
            f"http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2F{value.replace(':', '_')}"
        ),
        "refseq": f"https://www.ncbi.nlm.nih.gov/nuccore/{value}",
        "pdb": f"https://www.rcsb.org/structure/{value}",
        "complexportal": f"https://www.ebi.ac.uk/complexportal/complex/{value}",
        "ensembl_transcript": f"https://www.ensembl.org/id/{value}",
        "dbsnp": f"https://www.ncbi.nlm.nih.gov/snp/{value}",
    }
    if kind not in urls:
        raise ValueError(f"unknown official URL kind: {kind}")
    return urls[kind]


def provenance_rows(refs: Iterable[str], *, accessed: str = ACCESSED_DATE) -> list[dict[str, str]]:
    return [{"ref": ref, "source_url": ref, "accessed": accessed} for ref in sorted(set(refs))]


def make_entity(
    entity_id: str,
    kind: str,
    label: str,
    identifiers: Mapping[str, str],
    urls: Iterable[str],
    **extra: Any,
) -> dict[str, Any]:
    return {
        "id": entity_id,
        "kind": kind,
        "label": label,
        "identifiers": dict(sorted(identifiers.items())),
        "source_urls": sorted(set(urls)),
        "accessed": ACCESSED_DATE,
        "provenance": provenance_rows(urls),
        **extra,
    }


def make_edge(
    source: str,
    target: str,
    relation: str,
    status: str,
    refs: Iterable[str],
    claims: Iterable[str] = (),
) -> dict[str, Any]:
    return {
        "source": source,
        "target": target,
        "relation": relation,
        "epistemic_status": status,
        "evidence_required": True,
        "provenance": provenance_rows(refs),
        "linked_claim_ids": sorted(set(claims)),
    }


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def digest_bundle(bundle: Mapping[str, Any]) -> str:
    payload = {key: value for key, value in bundle.items() if key != "content_hash"}
    return hashlib.sha256(canonical_json(payload).encode()).hexdigest()


def finalize_bundle(bundle: dict[str, Any]) -> dict[str, Any]:
    bundle["entities"] = sorted(bundle["entities"], key=lambda item: item["id"])
    bundle["edges"] = sorted(bundle["edges"], key=lambda item: (item["source"], item["target"], item["relation"]))
    bundle["content_hash"] = digest_bundle(bundle)
    return bundle


def _walk(value: Any):
    if isinstance(value, Mapping):
        for key, child in value.items():
            yield str(key), child
            yield from _walk(child)
    elif isinstance(value, list):
        for child in value:
            yield from _walk(child)


def _valid_provenance(rows: Any) -> bool:
    return (
        isinstance(rows, list)
        and bool(rows)
        and all(
            isinstance(row, Mapping)
            and isinstance(row.get("source_url"), str)
            and row["source_url"].startswith("https://")
            and isinstance(row.get("accessed"), str)
            and bool(row["accessed"])
            for row in rows
        )
    )


def parse_genomic_hgvs(hgvs: str, accession: str) -> tuple[int, str, str]:
    if not isinstance(hgvs, str) or not isinstance(accession, str):
        raise ValueError("RefSeq accession must exactly match HGVS accession")
    match = _HGVS_RE.fullmatch(hgvs)
    if not match or match.group(1) != accession:
        raise ValueError("RefSeq accession must exactly match HGVS accession")
    return int(match.group(2)), match.group(3), match.group(4)


def validate_variant_record(record: Mapping[str, Any]) -> dict[str, Any]:
    if not isinstance(record, Mapping):
        raise ValueError("variant record must be a mapping")
    if any(key in _FORBIDDEN_KEYS for key, _ in _walk(record)):
        raise ValueError("forbidden variant field")
    required = {"variant_id", "assembly", "refseq_accession", "hgvs_genomic", "source_identifiers", "review_status"}
    missing = required - set(record)
    if missing:
        raise ValueError(f"variant missing required fields: {sorted(missing)}")
    if record["assembly"] != "GRCh38":
        raise ValueError("variant assembly must be GRCh38")
    accession, hgvs = record["refseq_accession"], record["hgvs_genomic"]
    if not isinstance(accession, str) or not _ACCESSION_RE.fullmatch(accession):
        raise ValueError("variant requires a versioned RefSeq genomic accession")
    position, ref_allele, alt_allele = parse_genomic_hgvs(hgvs, accession)
    if not isinstance(record["source_identifiers"], list) or not record["source_identifiers"] or not all(
        isinstance(item, str) and item for item in record["source_identifiers"]
    ):
        raise ValueError("variant source_identifiers are required")
    if record["review_status"] not in {"unreviewed", "reviewed", "curated"}:
        raise ValueError("invalid variant review status")
    variant_class = record.get("variant_class", "snv")
    if variant_class not in ALLOWED_VARIANT_CLASSES or variant_class == "haplotype":
        raise ValueError("invalid variant class")
    payload = copy.deepcopy(dict(record))
    payload["variant_class"] = "snv"
    payload["position"] = position
    payload["ref_allele"] = ref_allele
    payload["alt_allele"] = alt_allele
    gene_id = payload.get("gene_id")
    if gene_id is not None and (not isinstance(gene_id, str) or not re.fullmatch(r"gene:HGNC:\d+", gene_id)):
        raise ValueError("variant gene_id must be a gene entity ID")
    transcript = payload.get("transcript_accession")
    if transcript is not None and (not isinstance(transcript, str) or not _NM_RE.fullmatch(transcript)):
        raise ValueError("variant transcript_accession must be a versioned NM_ accession")
    hgvs_cdna = payload.get("hgvs_cdna")
    if hgvs_cdna:
        if not isinstance(hgvs_cdna, str) or not _CDNA_HGVS_RE.fullmatch(hgvs_cdna):
            raise ValueError("variant hgvs_cdna must be versioned NM_ coding HGVS")
        if transcript and not hgvs_cdna.startswith(f"{transcript}:"):
            raise ValueError("variant hgvs_cdna accession must match transcript_accession")
    hgvs_protein = payload.get("hgvs_protein")
    if hgvs_protein and (not isinstance(hgvs_protein, str) or not _PROTEIN_HGVS_RE.fullmatch(hgvs_protein)):
        raise ValueError("variant hgvs_protein must be versioned NP_ protein HGVS")
    return payload


def construct_variant_record(**kwargs: Any) -> dict[str, Any]:
    return validate_variant_record(kwargs)


def validate_haplotype_record(record: Mapping[str, Any], variants: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    if not isinstance(record, Mapping):
        raise ValueError("haplotype record must be a mapping")
    if any(key in _FORBIDDEN_KEYS for key, _ in _walk(record)):
        raise ValueError("forbidden haplotype field")
    required = {"haplotype_id", "gene_id", "label", "composition_of", "source_identifiers", "review_status", "assembly"}
    missing = required - set(record)
    if missing:
        raise ValueError(f"haplotype missing required fields: {sorted(missing)}")
    if record["assembly"] != "GRCh38":
        raise ValueError("haplotype assembly must be GRCh38")
    if record["review_status"] not in {"unreviewed", "reviewed", "curated"}:
        raise ValueError("invalid haplotype review status")
    if not isinstance(record["gene_id"], str) or not re.fullmatch(r"gene:HGNC:\d+", record["gene_id"]):
        raise ValueError("haplotype gene_id must be a gene entity ID")
    if not isinstance(record["composition_of"], list) or not record["composition_of"]:
        raise ValueError("haplotype composition_of is required")
    if not isinstance(record["source_identifiers"], list) or not record["source_identifiers"]:
        raise ValueError("haplotype source_identifiers are required")
    for item in record["composition_of"]:
        if not isinstance(item, str) or ":" not in item:
            raise ValueError("haplotype composition alleles must be variant_id:base")
        variant_id, allele = item.rsplit(":", 1)
        if variant_id not in variants:
            raise ValueError(f"unknown haplotype member: {variant_id}")
        member = variants[variant_id]
        if member.get("gene_id") != record["gene_id"]:
            raise ValueError("haplotype members must belong to the same gene")
        if allele not in {member["ref_allele"], member["alt_allele"]}:
            raise ValueError("haplotype allele must match member ref or alt")
    return copy.deepcopy(dict(record))


def _validate_programs(value: Any) -> tuple[str, ...]:
    if not isinstance(value, list) or not value or not all(isinstance(item, str) and item in ALLOWED_PROGRAMS for item in value):
        raise ValueError("invalid programs")
    if len(value) != len(set(value)):
        raise ValueError("duplicate programs")
    return tuple(value)


def validate_genomics_protein_spine(bundle: Mapping[str, Any]) -> bool:
    if not isinstance(bundle, Mapping) or bundle.get("content_hash") != digest_bundle(bundle) or not _HASH_RE.fullmatch(str(bundle.get("content_hash", ""))):
        raise ValueError("invalid identity-spine content hash")
    if (
        bundle.get("schema") != BUNDLE_SCHEMA
        or bundle.get("research_planning_only") is not True
        or bundle.get("claim_eligibility_mutation") is not False
        or bundle.get("owner") != "titans-med"
    ):
        raise ValueError("invalid planning-only bundle contract")
    variants = bundle.get("seeded_variants")
    if not isinstance(variants, list):
        raise ValueError("seeded_variants must be a list")
    haplotypes = bundle.get("seeded_haplotypes", [])
    if not isinstance(haplotypes, list):
        raise ValueError("seeded_haplotypes must be a list")
    if any(key in _FORBIDDEN_KEYS for key, _ in _walk(bundle)):
        raise ValueError("forbidden field")
    entities = bundle.get("entities", [])
    ids = [entity.get("id") for entity in entities]
    if len(ids) != len(set(ids)):
        raise ValueError("duplicate entity IDs")
    known = set(ids)
    for entity in entities:
        kind = entity.get("kind")
        if kind not in ALLOWED_KINDS:
            raise ValueError("invalid entity kind")
        if kind in {"gene", "transcript", "variant", "haplotype", "protein", "complex"}:
            if entity.get("taxonomy") != HUMAN_TAXONOMY or entity.get("species") != "Homo sapiens":
                raise ValueError("molecular entities require human taxonomy and species")
            if entity.get("locus") not in ALLOWED_LOCI:
                raise ValueError("molecular entities require a nuclear or mitochondrial locus")
        elif entity.get("taxonomy") != NON_SPECIES_SCOPE or entity.get("scope") != NON_SPECIES_SCOPE or "species" in entity or "locus" in entity:
            raise ValueError("non-molecular entities require non-species scope")
        _validate_programs(entity.get("programs"))
        if not entity.get("source_urls") or not _valid_provenance(entity.get("provenance")) or not isinstance(entity.get("identifiers"), Mapping):
            raise ValueError("missing HTTPS provenance or identifiers")
        ident = entity["identifiers"]
        patterns = {
            "gene": r"gene:HGNC:\d+",
            "transcript": r"transcript:NM_\d+\.\d+",
            "variant": r"variant:VAR-[A-Za-z0-9-]+",
            "haplotype": r"haplotype:HAP-[A-Za-z0-9-]+",
            "protein": r"protein:[A-Z0-9]+",
            "complex": r"complex:CPX-\d+",
            "process": r"process:GO:\d{7}",
            "tissue": r"tissue:UBERON:\d{7}",
        }
        if not re.fullmatch(patterns[kind], entity["id"]):
            raise ValueError("malformed entity identifier")
        checks = {
            "gene": [("hgnc", r"HGNC:\d+"), ("ncbi_gene", r"Gene:\d+"), ("ensembl", r"ENSG\d{11}")],
            "transcript": [("mane_select", r"NM_\d+\.\d+"), ("ensembl_transcript", r"ENST\d{11}\.\d+"), ("gene_symbol", r"[A-Z0-9-]+")],
            "variant": [("dbsnp", r"rs\d+"), ("hgvs_genomic", r"NC_\d{6}\.\d+:g\.[0-9]+[ACGT]>[ACGT]")],
            "haplotype": [("gene_id", r"gene:HGNC:\d+")],
            "protein": [("uniprot", r"[A-Z0-9]{6,10}"), ("alphafold", r"AF-[A-Z0-9]{6,10}-F1"), ("gene_symbol", r"[A-Z0-9-]+")],
            "complex": [("complexportal", r"CPX-\d+")],
            "process": [("go", r"GO:\d{7}")],
            "tissue": [("uberon", r"UBERON:\d{7}")],
        }
        if any(not re.fullmatch(pattern, str(ident.get(key, ""))) for key, pattern in checks[kind]):
            raise ValueError("malformed identifiers")
        if kind == "protein" and ident.get("alphafold") != f"AF-{ident['uniprot']}-F1":
            raise ValueError("AlphaFold accession must match UniProt identity")
        if kind == "protein":
            pdb = ident.get("pdb")
            if not isinstance(pdb, str) or not _PDB_RE.fullmatch(pdb):
                raise ValueError("protein requires a representative PDB accession")
            if official_url("pdb", pdb) not in entity["source_urls"]:
                raise ValueError("protein PDB provenance URL is required")
            architecture = entity.get("architecture")
            if architecture is not None:
                if not isinstance(architecture, list):
                    raise ValueError("protein architecture must be a list")
                for feature in architecture:
                    if not isinstance(feature, Mapping) or feature.get("start", 0) > feature.get("end", -1):
                        raise ValueError("malformed protein architecture feature")
        if kind == "gene" and entity["locus"] == "mitochondrial" and ident.get("refseq_accession") != "NC_012920.1":
            raise ValueError("mitochondrial genes require NC_012920.1")
        if kind == "gene" and entity["locus"] == "nuclear":
            chrom = ident.get("refseq_accession")
            if not isinstance(chrom, str) or not _ACCESSION_RE.fullmatch(chrom):
                raise ValueError("nuclear genes require a GRCh38 chromosome accession")
        if kind == "complex":
            go_id = ident.get("go")
            if go_id and not re.fullmatch(r"GO:\d{7}", str(go_id)):
                raise ValueError("malformed complex GO identifier")
        if kind == "variant":
            if entity.get("gene_id") not in known:
                raise ValueError("variant entity requires a catalog gene_id")
            if entity["id"] != f"variant:{entity.get('variant_id', '')}":
                raise ValueError("variant entity ID must match variant_id")
        if kind == "haplotype":
            if entity.get("gene_id") not in known:
                raise ValueError("haplotype entity requires a catalog gene_id")
            if entity["id"] != f"haplotype:{entity.get('haplotype_id', '')}":
                raise ValueError("haplotype entity ID must match haplotype_id")
        if kind == "transcript" and ident.get("mane_select") != entity["id"].split(":", 1)[1]:
            raise ValueError("transcript ID must match MANE Select accession")
    seen_variant_ids: set[str] = set()
    variant_index: dict[str, Mapping[str, Any]] = {}
    for record in variants:
        validated = validate_variant_record(record)
        variant_id = validated["variant_id"]
        if not isinstance(variant_id, str) or not variant_id:
            raise ValueError("variant_id is required")
        if variant_id in seen_variant_ids:
            raise ValueError("duplicate variant_id")
        gene_id = validated.get("gene_id")
        if not isinstance(gene_id, str) or gene_id not in known:
            raise ValueError("variant requires a catalog gene_id")
        gene = next(entity for entity in entities if entity["id"] == gene_id)
        if gene["identifiers"].get("refseq_accession") != validated["refseq_accession"]:
            raise ValueError("variant accession must match gene chromosome accession")
        transcript = validated.get("transcript_accession")
        if transcript and f"transcript:{transcript}" not in known:
            raise ValueError("variant transcript is not in catalog")
        seen_variant_ids.add(variant_id)
        variant_index[variant_id] = validated
    seen_haplotypes: set[str] = set()
    for record in haplotypes:
        validated = validate_haplotype_record(record, variant_index)
        haplotype_id = validated["haplotype_id"]
        if haplotype_id in seen_haplotypes:
            raise ValueError("duplicate haplotype_id")
        if validated["gene_id"] not in known:
            raise ValueError("haplotype requires a catalog gene_id")
        seen_haplotypes.add(haplotype_id)
    for variant_id in seen_variant_ids:
        if f"variant:{variant_id}" not in known:
            raise ValueError("seeded variant missing graph node")
    for haplotype_id in seen_haplotypes:
        if f"haplotype:{haplotype_id}" not in known:
            raise ValueError("seeded haplotype missing graph node")
    seen = set()
    for edge in bundle.get("edges", []):
        key = (edge.get("source"), edge.get("target"), edge.get("relation"))
        if key in seen:
            raise ValueError("duplicate edges")
        seen.add(key)
        if (
            edge.get("source") not in known
            or edge.get("target") not in known
            or edge.get("evidence_required") is not True
            or edge.get("epistemic_status") not in ALLOWED_EPISTEMIC
            or not _valid_provenance(edge.get("provenance"))
            or edge.get("linked_claim_ids") != []
        ):
            raise ValueError("invalid edge contract")
    return True


def _gene_entity(spec: GeneSpec) -> dict[str, Any]:
    identifiers = {
        "hgnc": spec.hgnc,
        "ncbi_gene": f"Gene:{spec.ncbi_gene}",
        "ensembl": spec.ensembl,
    }
    urls = [
        official_url("hgnc", spec.hgnc),
        official_url("ncbi", spec.ncbi_gene),
        official_url("ensembl", spec.ensembl),
    ]
    if spec.locus == "mitochondrial":
        identifiers["refseq_accession"] = "NC_012920.1"
        urls.append(official_url("refseq", "NC_012920.1"))
    else:
        chrom = NUCLEAR_CHROMOSOME_ACCESSION.get(spec.symbol, "")
        if not _ACCESSION_RE.fullmatch(chrom):
            raise ValueError(f"missing GRCh38 chromosome accession for {spec.symbol}")
        identifiers["refseq_accession"] = chrom
        urls.append(official_url("refseq", chrom))
    return make_entity(
        f"gene:{spec.hgnc}",
        "gene",
        spec.symbol,
        identifiers,
        urls,
        taxonomy=HUMAN_TAXONOMY,
        species="Homo sapiens",
        locus=spec.locus,
        programs=list(spec.programs),
    )


def _transcript_entity(gene: GeneSpec, transcript: TranscriptSpec) -> dict[str, Any]:
    return make_entity(
        f"transcript:{transcript.mane_select}",
        "transcript",
        f"{gene.symbol} MANE Select",
        {
            "mane_select": transcript.mane_select,
            "ensembl_transcript": transcript.ensembl_transcript,
            "gene_symbol": gene.symbol,
        },
        [
            official_url("refseq", transcript.mane_select),
            official_url("ensembl_transcript", transcript.ensembl_transcript.split(".", 1)[0]),
        ],
        taxonomy=HUMAN_TAXONOMY,
        species="Homo sapiens",
        locus=gene.locus,
        programs=list(gene.programs),
        transcript_select="MANE Select",
    )


def _protein_entity(spec: GeneSpec) -> dict[str, Any]:
    alphafold = f"AF-{spec.uniprot}-F1"
    if not _ALPHAFOLD_RE.fullmatch(alphafold):
        raise ValueError("malformed AlphaFold accession")
    pdb = REPRESENTATIVE_PDB.get(spec.uniprot, "")
    if not _PDB_RE.fullmatch(pdb):
        raise ValueError(f"missing representative PDB for {spec.uniprot}")
    architecture = [
        {"type": feature.feature_type, "start": feature.start, "end": feature.end, "label": feature.label}
        for feature in PROTEIN_ARCHITECTURE.get(spec.uniprot, ())
    ]
    return make_entity(
        f"protein:{spec.uniprot}",
        "protein",
        f"{spec.symbol} protein",
        {"uniprot": spec.uniprot, "gene_symbol": spec.symbol, "alphafold": alphafold, "pdb": pdb},
        [official_url("uniprot", spec.uniprot), official_url("alphafold", spec.uniprot), official_url("pdb", pdb)],
        taxonomy=HUMAN_TAXONOMY,
        species="Homo sapiens",
        locus=spec.locus,
        programs=list(spec.programs),
        architecture=architecture,
    )


def _complex_entity(spec: ComplexSpec) -> dict[str, Any]:
    identifiers = {"complexportal": spec.complexportal}
    urls = [official_url("complexportal", spec.complexportal)]
    if spec.go_id:
        identifiers["go"] = spec.go_id
        urls.append(official_url("go", spec.go_id))
    return make_entity(
        f"complex:{spec.complexportal}",
        "complex",
        spec.label,
        identifiers,
        urls,
        taxonomy=HUMAN_TAXONOMY,
        species="Homo sapiens",
        locus=spec.locus,
        programs=list(spec.programs),
    )


def build_genomics_protein_spine() -> dict[str, Any]:
    entities: list[dict[str, Any]] = []
    edges: list[dict[str, Any]] = []
    transcripts_by_symbol = {spec.gene_symbol: spec for spec in TRANSCRIPT_CATALOG}
    for spec in GENE_CATALOG:
        entities.append(_gene_entity(spec))
        entities.append(_protein_entity(spec))
        edges.append(
            make_edge(
                f"gene:{spec.hgnc}",
                f"protein:{spec.uniprot}",
                "encodes",
                "identity",
                [official_url("hgnc", spec.hgnc), official_url("uniprot", spec.uniprot)],
            )
        )
        if spec.locus == "nuclear":
            transcript = transcripts_by_symbol.get(spec.symbol)
            if transcript is None:
                raise ValueError(f"missing MANE Select transcript for {spec.symbol}")
            entities.append(_transcript_entity(spec, transcript))
            edges.append(
                make_edge(
                    f"gene:{spec.hgnc}",
                    f"transcript:{transcript.mane_select}",
                    "transcribed_to",
                    "identity",
                    [official_url("hgnc", spec.hgnc), official_url("refseq", transcript.mane_select)],
                )
            )
            edges.append(
                make_edge(
                    f"transcript:{transcript.mane_select}",
                    f"protein:{spec.uniprot}",
                    "translates_to",
                    "identity",
                    [official_url("refseq", transcript.mane_select), official_url("uniprot", spec.uniprot)],
                )
            )
        for process_id in spec.process_ids:
            edges.append(
                make_edge(
                    f"protein:{spec.uniprot}",
                    f"process:{process_id}",
                    "program_process_context",
                    "hypothesis",
                    [official_url("uniprot", spec.uniprot), official_url("go", process_id)],
                )
            )
    process_ids = {spec.go_id for spec in PROCESS_CATALOG}
    for spec in PROCESS_CATALOG:
        entities.append(
            make_entity(
                f"process:{spec.go_id}",
                "process",
                spec.label,
                {"go": spec.go_id},
                [official_url("go", spec.go_id)],
                taxonomy=NON_SPECIES_SCOPE,
                scope=NON_SPECIES_SCOPE,
                programs=list(spec.programs),
            )
        )
    tissue_ids = {spec.uberon_id for spec in TISSUE_CATALOG}
    for spec in TISSUE_CATALOG:
        entities.append(
            make_entity(
                f"tissue:{spec.uberon_id}",
                "tissue",
                spec.label,
                {"uberon": spec.uberon_id},
                [official_url("uberon", spec.uberon_id)],
                taxonomy=NON_SPECIES_SCOPE,
                scope=NON_SPECIES_SCOPE,
                programs=list(spec.programs),
            )
        )
    for process_id, tissues in ASSAY_CONTEXT_MATRIX.items():
        if process_id not in process_ids:
            raise ValueError(f"unknown assay-context process: {process_id}")
        for tissue_id in tissues:
            if tissue_id not in tissue_ids:
                raise ValueError(f"unknown assay-context tissue: {tissue_id}")
            edges.append(
                make_edge(
                    f"process:{process_id}",
                    f"tissue:{tissue_id}",
                    "assay_context",
                    "hypothesis",
                    [official_url("go", process_id)],
                )
            )
    uniprot_to_spec = {spec.uniprot: spec for spec in GENE_CATALOG}
    for spec in COMPLEX_CATALOG:
        entities.append(_complex_entity(spec))
        process_id, status = COMPLEX_PROCESS_EDGES[spec.complexportal]
        process_go = process_id.split(":", 1)[1]
        refs = [official_url("complexportal", spec.complexportal), official_url("go", spec.go_id or process_go)]
        edges.append(
            make_edge(
                f"complex:{spec.complexportal}",
                process_id,
                "participates_in",
                status,
                refs,
            )
        )
        for uniprot in spec.subunit_uniprots:
            if uniprot not in uniprot_to_spec:
                raise ValueError(f"unknown complex subunit: {uniprot}")
            edges.append(
                make_edge(
                    f"protein:{uniprot}",
                    f"complex:{spec.complexportal}",
                    "subunit_of",
                    "identity",
                    [official_url("uniprot", uniprot), official_url("complexportal", spec.complexportal)],
                )
            )
    seeded_variants = []
    for spec in SEEDED_VARIANT_SPECS:
        record = {
            "variant_id": spec.variant_id,
            "assembly": "GRCh38",
            "refseq_accession": spec.refseq_accession,
            "hgvs_genomic": spec.hgvs_genomic,
            "source_identifiers": list(spec.source_identifiers),
            "review_status": spec.review_status,
            "gene_id": f"gene:{spec.gene_hgnc}",
        }
        if spec.transcript_nm:
            record["transcript_accession"] = spec.transcript_nm
        if spec.hgvs_cdna:
            record["hgvs_cdna"] = spec.hgvs_cdna
        if spec.hgvs_protein:
            record["hgvs_protein"] = spec.hgvs_protein
        seeded_variants.append(construct_variant_record(**record))
    seeded_haplotypes = [
        {
            "haplotype_id": spec.haplotype_id,
            "assembly": "GRCh38",
            "gene_id": f"gene:{spec.gene_hgnc}",
            "label": spec.label,
            "composition_of": list(spec.composition_of),
            "source_identifiers": list(spec.source_identifiers),
            "review_status": spec.review_status,
        }
        for spec in SEEDED_HAPLOTYPE_SPECS
    ]
    genes_by_id = {entity["id"]: entity for entity in entities if entity["kind"] == "gene"}
    for record in seeded_variants:
        gene = genes_by_id[record["gene_id"]]
        dbsnp = next(item.split(":", 1)[1] for item in record["source_identifiers"] if item.startswith("dbSNP:"))
        identifiers = {
            "dbsnp": dbsnp,
            "hgvs_genomic": record["hgvs_genomic"],
            "refseq_accession": record["refseq_accession"],
        }
        urls = [official_url("dbsnp", dbsnp), official_url("refseq", record["refseq_accession"])]
        if record.get("hgvs_cdna"):
            identifiers["hgvs_cdna"] = record["hgvs_cdna"]
        if record.get("hgvs_protein"):
            identifiers["hgvs_protein"] = record["hgvs_protein"]
        entities.append(
            make_entity(
                f"variant:{record['variant_id']}",
                "variant",
                record["variant_id"],
                identifiers,
                urls,
                taxonomy=HUMAN_TAXONOMY,
                species="Homo sapiens",
                locus=gene["locus"],
                programs=list(gene["programs"]),
                variant_id=record["variant_id"],
                gene_id=record["gene_id"],
            )
        )
        edges.append(
            make_edge(
                record["gene_id"],
                f"variant:{record['variant_id']}",
                "has_variant",
                "identity",
                [official_url("dbsnp", dbsnp), gene["source_urls"][0]],
            )
        )
        transcript = record.get("transcript_accession")
        if transcript:
            edges.append(
                make_edge(
                    f"transcript:{transcript}",
                    f"variant:{record['variant_id']}",
                    "has_coding_variant",
                    "identity",
                    [official_url("refseq", transcript), official_url("dbsnp", dbsnp)],
                )
            )
    for record in seeded_haplotypes:
        gene = genes_by_id[record["gene_id"]]
        entities.append(
            make_entity(
                f"haplotype:{record['haplotype_id']}",
                "haplotype",
                record["label"],
                {"gene_id": record["gene_id"], "haplotype_id": record["haplotype_id"]},
                [gene["source_urls"][0]],
                taxonomy=HUMAN_TAXONOMY,
                species="Homo sapiens",
                locus=gene["locus"],
                programs=list(gene["programs"]),
                haplotype_id=record["haplotype_id"],
                gene_id=record["gene_id"],
                composition_of=list(record["composition_of"]),
            )
        )
        edges.append(
            make_edge(
                record["gene_id"],
                f"haplotype:{record['haplotype_id']}",
                "has_haplotype",
                "identity",
                [gene["source_urls"][0]],
            )
        )
        seen_members = []
        for item in record["composition_of"]:
            variant_id = item.rsplit(":", 1)[0]
            if variant_id in seen_members:
                continue
            seen_members.append(variant_id)
            edges.append(
                make_edge(
                    f"haplotype:{record['haplotype_id']}",
                    f"variant:{variant_id}",
                    "composed_of",
                    "identity",
                    [gene["source_urls"][0]],
                )
            )
    bundle = {
        "schema": BUNDLE_SCHEMA,
        "owner": "titans-med",
        "accessed": ACCESSED_DATE,
        "research_planning_only": True,
        "claim_eligibility_mutation": False,
        "seeded_variants": seeded_variants,
        "seeded_haplotypes": seeded_haplotypes,
        "entities": entities,
        "edges": edges,
    }
    finalize_bundle(bundle)
    validate_genomics_protein_spine(bundle)
    return bundle


def _walk_related_ids(catalog: Mapping[str, Any], seeds: Iterable[str], allowed_kinds: Iterable[str] | None = None) -> set[str]:
    wanted = set(seeds)
    kinds = frozenset(allowed_kinds) if allowed_kinds is not None else None
    by_id = {entity["id"]: entity for entity in catalog["entities"]}
    changed = True
    while changed:
        changed = False
        for edge in catalog["edges"]:
            if edge["source"] not in wanted or edge["target"] in wanted:
                continue
            target = by_id.get(edge["target"])
            if target is None:
                continue
            if kinds is not None and target["kind"] not in kinds:
                continue
            wanted.add(edge["target"])
            changed = True
    return wanted


def _records_for_genes(records: Iterable[Mapping[str, Any]], gene_ids: set[str]) -> list[dict[str, Any]]:
    return [copy.deepcopy(dict(record)) for record in records if record.get("gene_id") in gene_ids]


def project_spine(
    bundle: Mapping[str, Any],
    entity_ids: Iterable[str],
    *,
    include_variants: bool = False,
) -> dict[str, Any]:
    validate_genomics_protein_spine(bundle)
    wanted = set(entity_ids)
    entities = [copy.deepcopy(entity) for entity in bundle["entities"] if entity["id"] in wanted]
    found = {entity["id"] for entity in entities}
    missing = wanted - found
    if missing:
        raise ValueError(f"unknown entity IDs: {sorted(missing)}")
    edges = [
        copy.deepcopy(edge)
        for edge in bundle["edges"]
        if edge["source"] in wanted and edge["target"] in wanted
    ]
    gene_ids = {entity["id"] for entity in entities if entity["kind"] == "gene"}
    seeded_variants = _records_for_genes(bundle.get("seeded_variants", []), gene_ids) if include_variants else []
    seeded_haplotypes = _records_for_genes(bundle.get("seeded_haplotypes", []), gene_ids) if include_variants else []
    projected = {
        "schema": BUNDLE_SCHEMA,
        "owner": "titans-med",
        "accessed": bundle["accessed"],
        "research_planning_only": True,
        "claim_eligibility_mutation": False,
        "seeded_variants": seeded_variants,
        "seeded_haplotypes": seeded_haplotypes,
        "entities": entities,
        "edges": edges,
    }
    finalize_bundle(projected)
    validate_genomics_protein_spine(projected)
    return projected


def mito_atlas_entity_ids(bundle: Mapping[str, Any] | None = None) -> set[str]:
    catalog = bundle if bundle is not None else build_genomics_protein_spine()
    seeds = set()
    for entity in catalog["entities"]:
        if entity["kind"] == "gene" and entity["label"] in MITO_ATLAS_GENE_SYMBOLS:
            seeds.add(entity["id"])
        if entity["kind"] == "protein" and entity["identifiers"].get("gene_symbol") in MITO_ATLAS_GENE_SYMBOLS:
            seeds.add(entity["id"])
    return _walk_related_ids(catalog, seeds, MITO_ATLAS_KINDS)


def program_entity_ids(bundle: Mapping[str, Any], program: str) -> set[str]:
    if program not in ALLOWED_PROGRAMS:
        raise ValueError("unknown program")
    seeds = {entity["id"] for entity in bundle["entities"] if program in entity.get("programs", [])}
    return _walk_related_ids(bundle, seeds)


def program_projection(bundle: Mapping[str, Any] | None = None, program: str = "mito") -> dict[str, Any]:
    catalog = bundle if bundle is not None else build_genomics_protein_spine()
    return project_spine(catalog, program_entity_ids(catalog, program), include_variants=True)


def identity_links_for_compound(compound: str) -> list[dict[str, Any]]:
    slug = str(compound or "").strip().lower()
    by_symbol = {spec.symbol: spec for spec in GENE_CATALOG}
    by_process = {spec.go_id: spec for spec in PROCESS_CATALOG}
    by_complex = {spec.complexportal: spec for spec in COMPLEX_CATALOG}
    links: list[dict[str, Any]] = []
    for symbol in COMPOUND_GENE_SYMBOLS.get(slug, ()):
        spec = by_symbol[symbol]
        links.append(
            {
                "entity_id": f"gene:{spec.hgnc}",
                "kind": "gene",
                "label": spec.symbol,
                "relation": "mechanistic_context",
                "epistemic_status": "hypothesis",
                "evidence_required": True,
            }
        )
        links.append(
            {
                "entity_id": f"protein:{spec.uniprot}",
                "kind": "protein",
                "label": f"{spec.symbol} protein",
                "relation": "mechanistic_context",
                "epistemic_status": "hypothesis",
                "evidence_required": True,
            }
        )
    for process_id in COMPOUND_PROCESS_IDS.get(slug, ()):
        spec = by_process[process_id]
        links.append(
            {
                "entity_id": f"process:{spec.go_id}",
                "kind": "process",
                "label": spec.label,
                "relation": "mechanistic_context",
                "epistemic_status": "hypothesis",
                "evidence_required": True,
            }
        )
    for complex_id in COMPOUND_COMPLEX_IDS.get(slug, ()):
        spec = by_complex[complex_id]
        links.append(
            {
                "entity_id": f"complex:{spec.complexportal}",
                "kind": "complex",
                "label": spec.label,
                "relation": "mechanistic_context",
                "epistemic_status": "hypothesis",
                "evidence_required": True,
            }
        )
    return sorted(links, key=lambda item: item["entity_id"])


def mito_atlas_projection(bundle: Mapping[str, Any] | None = None) -> dict[str, Any]:
    catalog = bundle if bundle is not None else build_genomics_protein_spine()
    return project_spine(catalog, mito_atlas_entity_ids(catalog))


def build_genomics_graph(bundle: Mapping[str, Any] | None = None) -> nx.DiGraph:
    data = build_genomics_protein_spine() if bundle is None else bundle
    validate_genomics_protein_spine(data)
    graph = nx.DiGraph(name="Titans Med genomics/protein spine")
    graph.add_nodes_from((entity["id"], {**entity, "type": entity["kind"], "display_label": entity["label"]}) for entity in data["entities"])
    graph.add_edges_from((edge["source"], edge["target"], edge) for edge in data["edges"])
    return graph


def filter_graph_by_program(graph: nx.Graph, program: str | None) -> nx.Graph:
    """Return a deterministic, detached graph view for one research program."""
    if not program:
        return graph.copy()
    nodes = sorted(node for node, data in graph.nodes(data=True) if program in data.get("programs", []))
    return graph.subgraph(nodes).copy()


def identity_neighborhood(graph: nx.Graph, stable_id: str, depth: int = 1) -> nx.Graph:
    """Return an undirected one/two-hop neighborhood without changing ``graph``."""
    if depth not in (1, 2):
        raise ValueError("neighborhood depth must be 1 or 2")
    if stable_id not in graph:
        raise ValueError(f"unknown stable ID: {stable_id}")
    nodes = {stable_id}
    frontier = {stable_id}
    undirected = graph.to_undirected() if graph.is_directed() else graph
    for _ in range(depth):
        frontier = set().union(*(set(undirected.neighbors(node)) for node in frontier)) - nodes
        nodes.update(frontier)
    return graph.subgraph(sorted(nodes)).copy()


build_identity_neighborhood = identity_neighborhood


def build_layered_identity_figure(
    graph: nx.DiGraph,
    *,
    layers: Mapping[str, int],
    layer_labels: tuple[str, ...],
    width: int = 1400,
    height: int = 820,
    show_labels: bool = False,
) -> go.Figure:
    nodes = sorted(graph.nodes, key=lambda node: (layers[graph.nodes[node]["type"]], node))
    positions: dict[str, tuple[float, float]] = {}
    layer_indexes = sorted(set(layers.values()))
    for layer in layer_indexes:
        group = [node for node in nodes if layers[graph.nodes[node]["type"]] == layer]
        for index, node in enumerate(group):
            positions[node] = (layer * 2.1, (len(group) - 1) / 2 - index)
    grouped_edges: dict[str, list[tuple[str, str, Mapping[str, Any]]]] = {
        status: [] for status in ("identity", "context", "assay_context", "hypothesis")
    }
    for source, target, edge in sorted(graph.edges(data=True)):
        grouped_edges.setdefault(edge["epistemic_status"], []).append((source, target, edge))
    traces = []
    for status, edges in grouped_edges.items():
        if not edges:
            continue
        x_values, y_values, hovertext = [], [], []
        for source, target, edge in edges:
            x_values.extend((positions[source][0], positions[target][0], None))
            y_values.extend((positions[source][1], positions[target][1], None))
            claims = ", ".join(edge["linked_claim_ids"]) or "None"
            provenance = ", ".join(item["source_url"] for item in edge["provenance"])
            detail = (
                f"<b>{html.escape(graph.nodes[source]['display_label'])}</b> → "
                f"<b>{html.escape(graph.nodes[target]['display_label'])}</b><br>"
                f"Relation: {html.escape(edge['relation'])}<br>Epistemic status: {status}<br>"
                f"Review boundary: {'evidence required; review-required' if edge['evidence_required'] else 'identity/context only'}<br>"
                f"Linked claim IDs: {claims}<br>Source: {provenance or 'none recorded'}"
            )
            hovertext.extend((detail, detail, None))
        style = EDGE_STYLES.get(status, EDGE_STYLES["context"])
        traces.append(
            go.Scatter(
                x=x_values,
                y=y_values,
                mode="lines",
                name=f"Edges · {status}",
                line={"color": style["color"], "width": 1.5, "dash": style["dash"]},
                opacity=style["opacity"],
                showlegend=True,
                hovertext=hovertext,
                hovertemplate="%{hovertext}<extra></extra>",
            )
        )
    kind_order = ("compound", "project", "gene", "transcript", "variant", "haplotype", "protein", "complex", "process", "tissue")
    kind_symbols = {"gene": "diamond", "transcript": "hexagon", "variant": "circle", "haplotype": "circle-open", "protein": "square", "complex": "hexagram", "process": "triangle-up", "tissue": "star", "compound": "cross", "project": "hourglass"}
    kind_sizes = {"gene": 18, "transcript": 15, "variant": 12, "haplotype": 13, "protein": 19, "complex": 22, "process": 17, "tissue": 20, "compound": 16, "project": 18}
    for kind in kind_order:
        selected = [node for node in nodes if graph.nodes[node]["type"] == kind]
        if not selected:
            continue
        hovertext = []
        for node in selected:
            data = graph.nodes[node]
            identifiers = ", ".join(f"{key}: {value}" for key, value in data["identifiers"].items())
            scope = data.get("taxonomy", data.get("scope", "not specified"))
            sources = "<br>".join(f'<a href="{url}">{url}</a>' for url in data["source_urls"])
            programs = ", ".join(data.get("programs") or []) or "none recorded"
            locus = data.get("locus", "not specified")
            hovertext.append(
                f"<b>{html.escape(data['display_label'])}</b><br>Kind: {kind}<br>Stable ID: {html.escape(node)}<br>"
                f"Official identifiers: {html.escape(identifiers) or 'none recorded'}<br>"
                f"Locus/scope: {html.escape(str(locus if locus != 'not specified' else scope))}<br>"
                f"Programs: {html.escape(programs)}<br>Taxonomy: {html.escape(str(data.get('taxonomy', 'not specified')))}<br>"
                f"Source count: {len(data.get('source_urls', []))}<br>Official sources:<br>{sources or 'none recorded'}"
            )
        traces.append(
            go.Scatter(
                x=[positions[node][0] for node in selected],
                y=[positions[node][1] for node in selected],
                mode="markers+text" if show_labels else "markers",
                name=kind.title(),
                text=[graph.nodes[node]["display_label"] for node in selected],
                textposition="top center",
                textfont={"color": "#F4F7F5", "size": 12},
                marker={"size": kind_sizes.get(kind, 14), "symbol": kind_symbols.get(kind, "circle"), "color": KIND_COLORS.get(kind, "#b9cfcd"), "line": {"color": "#0a1116", "width": 1.5}},
                hovertext=hovertext,
                hovertemplate="%{hovertext}<extra></extra>",
            )
        )
    y_values = [y for _, y in positions.values()]
    y_min = min(y_values) if y_values else -1
    y_max = max(y_values) if y_values else 1
    span = y_max - y_min if y_values else 2
    padding = max(1.0, span * 0.12)
    figure = go.Figure(traces)
    figure.update_layout(shapes=[
        {"type": "rect", "xref": "x", "yref": "paper", "x0": index * 2.1 - 0.72, "x1": index * 2.1 + 0.72, "y0": 0, "y1": 1, "fillcolor": "rgba(110,216,208,0.035)" if index % 2 == 0 else "rgba(232,185,74,0.035)", "line": {"width": 0}, "layer": "below"}
        for index in layer_indexes
    ])
    figure.update_layout(
        width=width,
        height=height,
        template="plotly_dark",
        paper_bgcolor="#101820",
        plot_bgcolor="#101820",
        showlegend=True,
        xaxis={"visible": False},
        yaxis={"visible": False, "range": [y_min - padding, y_max + padding], "fixedrange": True},
        margin={"l": 20, "r": 20, "t": 70, "b": 20},
        annotations=[
            {
                "x": index * 2.1,
                "y": 1.04,
                "xref": "x",
                "yref": "paper",
                "text": label,
                "showarrow": False,
                "font": {"color": "#F4F7F5", "size": 15},
            }
            for index, label in enumerate(layer_labels)
        ],
    )
    if not nodes:
        figure.add_annotation(
            x=0.5,
            y=0.5,
            xref="paper",
            yref="paper",
            text="No molecular identities match the current search or filters.",
            showarrow=False,
            font={"color": "#F4F7F5", "size": 16},
        )
    return figure


def build_genomics_figure(graph: nx.DiGraph | None = None, *, width: int = 1400, height: int = 820, show_labels: bool = False) -> go.Figure:
    graph = build_genomics_graph() if graph is None else graph
    return build_layered_identity_figure(
        graph,
        layers={"gene": 0, "transcript": 1, "variant": 2, "haplotype": 2, "protein": 3, "complex": 4, "process": 5, "tissue": 6},
        layer_labels=("Gene", "Transcript", "Allele", "Protein", "Complex", "Process", "Tissue"),
        width=width,
        height=height,
        show_labels=show_labels,
    )


def resolve_identity_query(query: str, bundle: Mapping[str, Any] | None = None, *, limit: int = 12) -> list[dict[str, Any]]:
    catalog = bundle if bundle is not None else build_genomics_protein_spine()
    validate_genomics_protein_spine(catalog)
    text = re.sub(r"\s+", " ", str(query or "")).strip()
    if not text:
        return []
    lowered = text.lower()
    tokens = set(re.findall(r"[a-z0-9:-]+", lowered))
    hits: list[dict[str, Any]] = []
    for entity in catalog["entities"]:
        ident = entity.get("identifiers") or {}
        haystacks = {
            entity["id"].lower(),
            entity["label"].lower(),
            str(ident.get("gene_symbol", "")).lower(),
            str(ident.get("dbsnp", "")).lower(),
            str(ident.get("uniprot", "")).lower(),
            str(ident.get("hgnc", "")).lower(),
            str(ident.get("haplotype_id", "")).lower(),
            str(entity.get("variant_id", "")).lower(),
        }
        score = 0
        if lowered in haystacks or any(value == lowered for value in haystacks if value):
            score = 100
        elif any(token and token in haystacks for token in tokens):
            score = 80
        elif any(token and token in value for token in tokens for value in haystacks if len(token) >= 4):
            score = 40
        if score:
            hits.append(
                {
                    "id": entity["id"],
                    "kind": entity["kind"],
                    "label": entity["label"],
                    "programs": list(entity.get("programs") or []),
                    "identifiers": dict(ident),
                    "score": score,
                    "research_planning_only": True,
                }
            )
    hits.sort(key=lambda item: (-item["score"], item["kind"], item["id"]))
    return hits[: max(1, int(limit))]


def catalog_snapshot(bundle: Mapping[str, Any] | None = None) -> dict[str, Any]:
    catalog = bundle if bundle is not None else build_genomics_protein_spine()
    validate_genomics_protein_spine(catalog)
    counts = {
        kind: sum(1 for entity in catalog["entities"] if entity["kind"] == kind)
        for kind in ("gene", "transcript", "variant", "haplotype", "protein", "complex", "process", "tissue")
    }
    return {
        "schema": catalog["schema"],
        "accessed": catalog["accessed"],
        "content_hash": catalog["content_hash"],
        "research_planning_only": True,
        "counts": counts,
        "seeded_variants": len(catalog.get("seeded_variants", [])),
        "seeded_haplotypes": len(catalog.get("seeded_haplotypes", [])),
    }


__all__ = [
    "ACCESSED_DATE",
    "BUNDLE_SCHEMA",
    "COMPLEX_CATALOG",
    "GENE_CATALOG",
    "MITO_ATLAS_GENE_SYMBOLS",
    "MTDNA_GENES",
    "NUCLEAR_MITO_GENES",
    "NUCLEAR_PROGRAM_GENES",
    "PROCESS_CATALOG",
    "SEEDED_HAPLOTYPE_SPECS",
    "SEEDED_VARIANT_SPECS",
    "REPRESENTATIVE_PDB",
    "TRANSCRIPT_CATALOG",
    "TISSUE_CATALOG",
    "identity_links_for_compound",
    "build_genomics_figure",
    "build_genomics_graph",
    "filter_graph_by_program",
    "identity_neighborhood",
    "build_identity_neighborhood",
    "build_genomics_protein_spine",
    "build_layered_identity_figure",
    "canonical_json",
    "catalog_snapshot",
    "construct_variant_record",
    "digest_bundle",
    "finalize_bundle",
    "make_edge",
    "make_entity",
    "mito_atlas_entity_ids",
    "mito_atlas_projection",
    "official_url",
    "parse_genomic_hgvs",
    "program_entity_ids",
    "program_projection",
    "resolve_identity_query",
    "project_spine",
    "provenance_rows",
    "validate_genomics_protein_spine",
    "validate_haplotype_record",
    "validate_variant_record",
]
