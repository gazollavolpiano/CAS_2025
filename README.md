## Study overview

The **Childhood Asthma Study (CAS)** is a longitudinal birth cohort. This repository contains analysis code supporting the manuscript's characterisation of the nasopharyngeal (NP) microbiome in **58 infants** sampled at **2, 6, and 12 months of age**.

- **1,036 bacterial isolates** were cultured from NP plate sweeps across five selective media (BLOO, CHOC, COBA, HAEM, MORA)
- **1,028 genomes** were assembled; **925** passed quality thresholds and were retained for analysis
- Isolates were assigned to **6 genera / 72 species** using GTDB-Tk r214
- **143 plate-sweep metagenomes** were processed with Kraken2/Bracken and StrainScan to compare strain detection with culture-based methods
- **11 dominant species** were profiled for strain-level dynamics using StrainScan

## Bioinformatics pipeline summary

```
Raw reads (Illumina paired-end, isolates)
        │
        ▼
  Unicycler (de novo assembly)  ──►  seqkit (filter contigs ≥500 bp)
        │
        ▼
  CheckM (completeness / contamination)
  Quality score = completeness − (5 × contamination)
        │
        ▼
  GTDB-Tk r214 (taxonomy)
        │
        ├──► mlst (MLST)                        [S. aureus, S. pneumoniae,
        │                                         S. agalactiae, M. catarrhalis]
        ├──► AgrVATE (agr group)                [S. aureus]
        ├──► spaTyper (spa type)                [S. aureus]
        ├──► Staphopia-SCCmec (SCCmec type)     [S. aureus]
        ├──► PneumoCaT (capsular serotype)      [S. pneumoniae]
        ├──► GBS-SBG (serotype)                 [S. agalactiae]
        ├──► MiST + moraxella_mist_to_partial_lincode.py
        │    (cgMLST → partial LIN codes)       [M. catarrhalis]
        └──► abricate — VFDB + ResFinder        [all CAS genomes]

Public genomes (NCBI Datasets)
        │
        ├──► StrainScan build (per-species databases, 11 species)
        └──► SKA2 + rapidNJ (core SNP phylogenies)

Raw reads (Illumina paired-end, plate-sweep metagenomes)
        │
        ├──► Kraken2 (HPRC db) → KrakenTools (human read removal)
        │         │
        │         ▼
        │    Kraken2 (GTDB r214 db, confidence=0.1)
        │         │
        │         ▼
        │    Bracken (species-level, read length 150, threshold 250)
        │         │
        │         ▼
        │    R: combine_bracken_outputs → phyloseq object
        │
        └──► StrainScan scan (11 species databases × 143 samples)
```

***

## Software and versions

### Command-line tools

| Tool | Version | Purpose |
|---|---|---|
| Unicycler | 0.5.0 | De novo genome assembly |
| seqkit | 2.6.0 | Contig filtering and assembly stats |
| CheckM | 1.2.3 | Genome completeness and contamination |
| GTDB-Tk | 2.3.0 (r214) | Taxonomic classification |
| mlst | 2.23.0 | Multi-locus sequence typing |
| AgrVATE | 1.0.2 | *S. aureus* agr group |
| spaTyper | 0.2.1 | *S. aureus* spa type |
| Staphopia-SCCmec | 1.0.0 | *S. aureus* SCCmec typing |
| PneumoCaT | 1.2.1 | *S. pneumoniae* capsular serotype |
| GBS-SBG | — | *S. agalactiae* serotype |
| MiST | 1.2.0 | Moraxella cgMLST (PubMLST scheme) |
| abricate | 1.0.1 | Virulence factors and ARG profiling |
| Kraken2 | 2.1.3 | Metagenomic taxonomic classification |
| Bracken | 2.9 | Species-level abundance re-estimation |
| StrainScan | 1.0.14 | Strain-level metagenomic profiling |
| SKA2 | 0.3.7 | Core SNP alignment |
| rapidNJ | 2.3.2 | Neighbour-joining phylogeny |
| NCBI Datasets | 16.31.0 | Public genome download |
