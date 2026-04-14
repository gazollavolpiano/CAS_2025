# Genome Assembly and Taxonomic Classification for Isolates

```bash
# Genome Assembly with Unicycler
unicycler -1 SAMPLE/cleaned_R1.fastq.gz -2 SAMPLE/cleaned_R2.fastq.gz -s SAMPLE/unpaired_R1.fastq.gz -o SAMPLE/unicycler_output 

# Filtering Small Contigs with seqkit
seqkit seq --min-len 500 SAMPLE/unicycler_output/assembly.fasta --remove-gaps > SAMPLE/filtered_assembly.fasta

# Completeness and Contamination with CheckM
checkm lineage_wf --extension fasta --file checkm_results.tsv SAMPLE/filtered_assembly.fasta SAMPLE/checkm_output

# N50 and Other Stats with seqkit
seqkit stats SAMPLE/filtered_assembly.fasta -T --all > SAMPLE/stats_report.tsv

# Taxonomic Classification with GTDB-Tk
gtdbtk classify_wf --genome_dir SAMPLE/filtered_assembly.fasta --extension fasta --out_dir SAMPLE/gtdbtk_output --mash_db gtdb-tk-r214.msh 
```

# Isolates Typing and Capsular Characterization

```bash
# MLST with mlst Tool
mlst --csv --legacy --scheme scheme SAMPLE/filtered_assembly.fasta

# Download Data on Clonal Complexes
wget "https://rest.pubmlst.org/db/pubmlst_scheme_seqdef/schemes/1/profiles_csv" # replace db

# SCCmec with staphopia-sccmec Tool
staphopia-sccmec SAMPLE/filtered_assembly.fasta

# spa Typing with spatyper Tool
spaTyper -f SAMPLE/filtered_assembly.fasta --output SAMPLE_spatyper.txt

# agr Groups with agrvate Tool
agrvate --typing-only -i SAMPLE/filtered_assembly.fasta

# Capsular Typing with PneumoCaT Tool
PneumoCaT.py  --fastq_1 SAMPLE/cleaned_R1.fastq.gz --fastq_2 SAMPLE/cleaned_R2.fastq.gz --output_dir SAMPLE/PneumoCaT_output

# Serotyping with GBS-SBG Tool
./GBS-SBG.pl SAMPLE/filtered_assembly.fasta -name SAMPLE -best > serotyping_results.txt

# cgMLST with MiST

## retrieve the scheme and its profiles:
mist download \
  --downloader bigsdb_auth \
  --url https://rest.pubmlst.org/db/pubmlst_moraxella_seqdef/schemes/2\
  --output moraxella_cgmlst \
  --include-profiles \
  --dir-tokens .bigsdb_tokens \
  --key-name PubMLST \
  --site PubMLST

## create the MiST index
mist index --fasta-list moraxella_cgmlst/fasta_list.txt --output db_moraxella_cgmlst 

## profile the samples with MiST
mist call --fasta SAMPLE/filtered_assembly.fasta --db db_moraxella_cgmlst --out-json SAMPLE.json

## determine the closest match and calculate a partial LINcode
python3 moraxella_mist_to_partial_lincode.py SAMPLE.json moraxella_cgmlst/profiles.tsv 
```

# Profiling of Virulence Factors and Antimicrobial Resistance Genes 

```bash
## Virulence Factors with ABRicate Tool
abricate SAMPLE/filtered_assembly.fasta --db vfdb --minid 90 --mincov 60 > SAMPLE_vfdb.tab

## Resistance Genes with ABRicate Tool
abricate SAMPLE/filtered_assembly.fasta --db resfinder --minid 90 --mincov 60 > SAMPLE_resfinder.tab
```

# Retrieval and Curation of Public Bacterial Genomes and Phylogenetic Analysis 

```bash
## Downloading Genomes with NCBI Datasets Genome Package 
datasets download genome accession GCA_000419095.2 # example accession

## Phylogenetic Analysis with SKA2 and rapidNJ
ska build -o SAMPLES_al_seqs -k 31 -f SAMPLES_path.list
ska align -o SAMPLES_al_seqs.aln SAMPLES_al_seqs.skf 
rapidnj SAMPLES_al_seqs.aln -i fa > SAMPLES_NJ_tree.nwk
```

# Taxonomic Classification of Species in Plate Sweep Metagenomes

```bash
## Kraken2 with HPRC indices
kraken2 --db /k2_HPRC_20230810/db/ --paired --gzip-compressed  SAMPLE_R1.fastq.gz SAMPLE_R2.fastq.gz --output SAMPLE.HPRC.kraken

## Read extraction with KrakenTools
python extract_kraken_reads.py -k SAMPLE.HPRC.kraken -s1 SAMPLE_R1.fastq.gz -s2 SAMPLE_R2.fastq.gz -t 0 -o nonhuman_SAMPLE.1.fastq -o2 nonhuman_SAMPLE.2.fastq --fastq-output

## Run Kraken2 with GTDB indices
kraken2 --db /gtdb_r214/128gb/ --output SAMPLE.kraken --report SAMPLE.kraken.report --paired tempnonhuman_SAMPLE.1.fastq nonhuman_SAMPLE.2.fastq --confidence 0.1 

## Run Bracken with Kraken2 report
bracken -d /gtdb_r214/128gb/ -i SAMPLE.kraken.report -o SAMPLE.bracken.S250.report -r 150 -l S -t 250

## Combine outputs in a single file
python combine_bracken_outputs_mod.py --files *bracken.S250.report -o combined_bracken.S250.report
```

Use R to create a phyloseq object from the combined bracken output file.

```R
# library
library(phyloseq)
library(tidyverse)

# read merged bracken output result and modify it
merged_bracken_output <- read.table("combined_bracken.S250.report", header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
                          select(-matches("frac|taxonomy_id|taxonomy_lvl")) %>% 
                          rename(Species = name) 

merged_bracken_output <- merged_bracken_output %>% rename_all(~ str_remove(., ".bracken.S250.report_num"))

# read the GTDB taxonomy file (available at https://data.gtdb.ecogenomic.org/releases/release214/214.0/)
taxonomy <- read.table("gtdb_v214_taxonomy.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# add the full taxonomy to the merged bracken output 
merged_bracken_output <- merge(taxonomy, merged_bracken_output, all.y = TRUE)

# create otu table and tax table
otu <- merged_bracken_output %>% 
      select(-Domain, -Phylum, -Class, -Order, -Family, -Genus, Species) %>%
      column_to_rownames("Species")

tax <- merged_bracken_output %>% select(Domain, Phylum, Class, Order, Family, Genus, Species) 
rownames(tax) <- tax$Species
 
# create phyloseq objectq
physeq <- phyloseq(otu_table(as.matrix(otu), taxa_are_rows = TRUE), tax_table(as.matrix(tax)))

# write phyloseq object
saveRDS(physeq, "phyloseq.rds")
```

# Temporal Analysis of Strains

```bash
## Database creation with StrainScan 
strainscan_build -i GENOME_DIRECTORY/ -o strainscan_db 

## Profile Strains with StrainScan
strainscan -i SAMPLE_R1.fastq.gz -j SAMPLE_R2.fastq.gz -d strainscan_db -o SAMPLE -o SAMPLE/
```
