# Create Table with Summary Data for 1,036 Isolates 

Only genomes from genera with more than four isolates were retained for analysis.

The available data to create the table can be found in the following sources:
  - `/CAS_genomes/genomes_fastq/`: 1,036 fastq files with reads for the isolates
  - `/CAS_genomes/genomes_fasta/assemblies/`: 1,028 genomes created from the fastq files
  - `/12_CAS/4_gtdbtk/maldi_results_per_isolate.csv`: has columns Sample_ID | Species for maldi results
  - `/12_CAS/4_gtdbtk/log_age_subjectid.csv`: has columns log.number | age_months | subject.id
  - `/12_CAS/CAS_425days_OTUmat_clean_culturedLog.csv`: has several columns with clinical/phenotype for each sample
  - `/12_CAS/2_checkm/check_res.tsv`: CheckM results for the genomes
  - `/12_CAS/2_checkm/fasta_stats.tsv`: has statistics on genome size, GC content, and N50 for the genomes etc.
  - `/12_CAS/4_gtdbtk/gtdbtk.bac120.summary.tsv`: taxonomical classification for the genomes
  - `/12_CAS/7_final_genome_db_and_metadata/genomes_db_925_CAS_2939_public/`: final collection of 925 genomes 

```R
# load libraries
library(tidyverse)

# list initial isolates (based on the fastq files)
isolates <- list.files("/CAS_genomes/genomes_fastq/", full.names = TRUE, pattern = "R1.fastq.gz") %>%
  str_extract("CAS_\\d+_\\d+_.*_R1.fastq.gz") %>%
  str_remove("_R1.fastq.gz") %>%
  tibble(Isolate = .)

isolates$Culture.Medium <- gsub("CAS_\\d+_\\d+_(.*)", "\\1", isolates$Isolate) # extract culture medium
isolates$Patch.Plate.Number <- gsub("CAS_\\d+_(\\d+)_.*", "\\1", isolates$Isolate) # extract patch plate number

dim(isolates) # 1036x3

# list genomes that were assembled
assembled <- list.files("/CAS_genomes/genomes_fasta/assemblies/", full.names = TRUE) %>% 
  str_extract("CAS_\\d{4}.*.fasta") %>% 
  str_remove(".fasta") %>%
  tibble(Isolate = .)

isolates$Assembly <- ifelse(isolates$Isolate %in% assembled$Isolate, "Yes", "No(low read count)")
table(isolates$Assembly) # Yes 1028, No(low read count) 8

# read MALDI-TOF results and join with isolates
maldi_results_per_isolate <- read.csv("/12_CAS/4_gtdbtk/maldi_results_per_isolate.csv") %>% filter(Species != "plate_sweep")
colnames(maldi_results_per_isolate) <- c("Isolate", "MALDI.TOF.Identification")
dim(maldi_results_per_isolate) # 1036x2
isolates <- left_join(isolates, maldi_results_per_isolate, by = "Isolate")

# read CheckM results and join with isolates
check_res <- read.delim("/12_CAS/2_checkm/check_res.tsv", skip = 3, header = FALSE, sep = "")
check_res <- check_res[-nrow(check_res), ]
colnames(check_res) <- c("Bin_Id", "Marker_lineage", "Marker_lineage_code",  "Genomes", "Markers", "Marker_sets", "0", "1", "2", "3", "4", "5_plus", "Completeness", "Contamination", "Strain_heterogeneity")
check_res <- check_res %>% select(Bin_Id, Completeness, Contamination) %>% mutate(Isolate = str_remove(Bin_Id, "filt_")) %>% select(-Bin_Id)
check_res$Quality.Score <- check_res$Completeness - 5 * check_res$Contamination
isolates <- left_join(isolates, check_res, by = "Isolate")

# read stats on genome size and N50 for the genomes, and join with isolates
fasta_stats <- read.delim("/2024/12_CAS/2_checkm/fasta_stats.tsv", sep = "\t") %>%
  mutate(Isolate = str_extract(file, "CAS_.*.fasta") %>% str_remove(".fasta")) %>%
  rename(Contigs.Number = num_seqs) %>%
  select(Isolate, Contigs.Number, N50) 
dim(fasta_stats) # 1028x3
isolates <- left_join(isolates, fasta_stats, by = "Isolate")

# read GTDB-Tk results and join with isolates
gtdbtk <- read.delim("/12_CAS/4_gtdbtk/gtdbtk.bac120.summary.tsv", sep = "\t") %>%
  select(user_genome = user_genome, classification = classification) %>%
  rename(Isolate = user_genome, Taxonomy.GTDB.r214 = classification) %>%
  mutate(Isolate = str_remove(Isolate, "filt_"))
dim(gtdbtk) # 944x2
isolates <- left_join(isolates, gtdbtk, by = "Isolate")

# add columns for Infant Log Number and Sampling Age (months)
log_age_subjectid <- read.csv("/12_CAS/4_gtdbtk/log_age_subjectid.csv") %>%
  rename(Subject.ID = subject.id, Sampling.Age = age_months, Sample.Log.Number = log.number) 
dim(log_age_subjectid) # 149x3

isolates$Sample.Log.Number <- as.numeric(gsub("CAS_(\\d+)_\\d+.*", "\\1", isolates$Isolate)) # extract log number

setdiff(isolates$Sample.Log.Number, log_age_subjectid$Sample.Log.Number) # 231 is missing from log_age_subjectid
read.csv("/12_CAS/CAS_425days_OTUmat_clean_culturedLog.csv") %>% filter(Log.No == 231) # 231 is missing from pheno too!

# join with log_age_subjectid with isolates
isolates <- left_join(isolates, log_age_subjectid, by = "Sample.Log.Number")

# add column for Retained for Analysis
final_collection <- list.files("/12_CAS/7_final_genome_db_and_metadata/genomes_db_925_CAS_2939_public/", full.names = TRUE, pattern = "filt_CAS") %>%
  str_extract("CAS_\\d+_\\d+_.*.fasta") %>%
  str_remove(".fasta") %>%
  tibble(Isolate = .)
length(final_collection$Isolate) # 925

isolates$`Retained for Analysis*` <- ifelse(isolates$Isolate %in% final_collection$Isolate, "Yes", "No")
table(isolates$`Retained for Analysis*`) # Yes 925, No 111

# tidy column order
isolates <- isolates %>% 
select(Isolate, MALDI.TOF.Identification, Culture.Medium, Patch.Plate.Number, Assembly, Completeness, 
        Contamination, Quality.Score, Contigs.Number, N50, Taxonomy.GTDB.r214, Subject.ID, Sample.Log.Number, Sampling.Age, `Retained for Analysis*`)

# recategorize Sampling.Age column
isolates <- isolates %>%
  mutate(Sampling.Age = case_when(
    Sampling.Age %in% c(1, 2, 3) ~ "2",
    Sampling.Age %in% c(6, 7) ~ "6",
    Sampling.Age %in% c(11, 12, 13) ~ "12",
    TRUE ~ as.character(Sampling.Age)
  )) 

# Replace cells with NA with ""
isolates <- isolates %>% mutate_all(~replace(., is.na(.), ""))

# save table
write_csv(isolates, "01_tableS1/01_tableS1.csv")

# quick check on number of infants
isolates %>% 
filter(`Retained for Analysis*` == "Yes") %>%
pull(Subject.ID) %>%
unique() %>%
length() # 58

# answer: how many infants sampled at 2 months, 6 months, and 12 months?
isolates %>%
filter(`Retained for Analysis*` == "Yes") %>%
select(Subject.ID, Sampling.Age) %>%
distinct() %>%
group_by(Sampling.Age) %>%
summarise(n = n()) %>%
as.data.frame()

#   Sampling.Age  n
# 1               1
# 2           12 40
# 3            2 50
# 4            6 43

# quick check on species distribution
isolates %>% 
filter(`Retained for Analysis*` == "Yes") %>%
mutate(Species = str_extract(Taxonomy.GTDB.r214, "g__(.*)")) %>%
group_by(Species) %>% 
summarise(n = n()) %>% 
arrange(desc(n)) %>% 
as.data.frame()

#                                                       Species   n
# 1                  g__Staphylococcus;s__Staphylococcus aureus 231
# 2                       g__Moraxella;s__Moraxella catarrhalis 165
# 3  g__Corynebacterium;s__Corynebacterium pseudodiphtheriticum  96
# 4                  g__Dolosigranulum;s__Dolosigranulum pigrum  58
# 5                g__Streptococcus;s__Streptococcus pneumoniae  46
# 6              g__Corynebacterium;s__Corynebacterium accolens  26
# 7            g__Corynebacterium;s__Corynebacterium propinquum  24
# 8                g__Streptococcus;s__Streptococcus agalactiae  23
# 9                             g__Rothia;s__Rothia sp902373285  22
# 10                  g__Moraxella;s__Moraxella nonliquefaciens  19
# 11                                       g__Streptococcus;s__  18
# 12            g__Staphylococcus;s__Staphylococcus lugdunensis  14
# 13                g__Streptococcus;s__Streptococcus lactarius  12
# 14                g__Staphylococcus;s__Staphylococcus hominis  10
# 15                 g__Streptococcus;s__Streptococcus mitis_AQ   9
# 16         g__Streptococcus;s__Streptococcus pseudopneumoniae   7
# 17          g__Corynebacterium;s__Corynebacterium sp943914355   6
# 18                g__Staphylococcus;s__Staphylococcus capitis   6
# 19            g__Staphylococcus;s__Staphylococcus ureilyticus   6
# 20                 g__Streptococcus;s__Streptococcus mitis_AD   6
# 21          g__Staphylococcus;s__Staphylococcus saprophyticus   5
# 22             g__Streptococcus;s__Streptococcus dysgalactiae   5
# 23                    g__Streptococcus;s__Streptococcus mitis   5
# 24                   g__Streptococcus;s__Streptococcus oralis   5
# 25              g__Streptococcus;s__Streptococcus sp001556435   5
# 26              g__Streptococcus;s__Streptococcus sp903645285   5
# 27                         g__Rothia;s__Rothia mucilaginosa_A   4
# 28               g__Staphylococcus;s__Staphylococcus pasteuri   4
# 29              g__Streptococcus;s__Streptococcus intermedius   4
# 30                  g__Streptococcus;s__Streptococcus mitis_I   4
# 31              g__Streptococcus;s__Streptococcus sp019448285   4
# 32                                              g__Rothia;s__   3
# 33                  g__Streptococcus;s__Streptococcus mitis_D   3
# 34                g__Streptococcus;s__Streptococcus oralis_BA   3
# 35          g__Streptococcus;s__Streptococcus parasanguinis_E   3
# 36               g__Streptococcus;s__Streptococcus salivarius   3
# 37              g__Streptococcus;s__Streptococcus sp000187445   3
# 38              g__Streptococcus;s__Streptococcus sp901875575   3
# 39                                     g__Corynebacterium;s__   2
# 40                           g__Rothia;s__Rothia mucilaginosa   2
# 41            g__Staphylococcus;s__Staphylococcus epidermidis   2
# 42       g__Staphylococcus;s__Staphylococcus pseudintermedius   2
# 43                 g__Streptococcus;s__Streptococcus infantis   2
# 44                  g__Streptococcus;s__Streptococcus mitis_A   2
# 45                 g__Streptococcus;s__Streptococcus mitis_AC   2
# 46                 g__Streptococcus;s__Streptococcus mitis_AF   2
# 47                 g__Streptococcus;s__Streptococcus mitis_BB   2
# 48                 g__Streptococcus;s__Streptococcus mitis_BM   2
# 49                  g__Streptococcus;s__Streptococcus mitis_K   2
# 50          g__Streptococcus;s__Streptococcus parasanguinis_L   2
# 51       g__Streptococcus;s__Streptococcus pseudopneumoniae_C   2
# 52              g__Corynebacterium;s__Corynebacterium coyleae   1
# 53      g__Corynebacterium;s__Corynebacterium kefirresidentii   1
# 54       g__Corynebacterium;s__Corynebacterium minutissimum_A   1
# 55             g__Corynebacterium;s__Corynebacterium otitidis   1
# 56             g__Corynebacterium;s__Corynebacterium simulans   1
# 57                           g__Rothia;s__Rothia dentocariosa   1
# 58                                 g__Rothia;s__Rothia terrae   1
# 59                 g__Staphylococcus;s__Staphylococcus cohnii   1
# 60           g__Staphylococcus;s__Staphylococcus haemolyticus   1
# 61               g__Staphylococcus;s__Staphylococcus simulans   1
# 62                g__Staphylococcus;s__Staphylococcus warneri   1
# 63                g__Streptococcus;s__Streptococcus anginosus   1
# 64                  g__Streptococcus;s__Streptococcus mitis_P   1
# 65                g__Streptococcus;s__Streptococcus oralis_BC   1
# 66                 g__Streptococcus;s__Streptococcus oralis_I   1
# 67                 g__Streptococcus;s__Streptococcus oralis_S   1
# 68             g__Streptococcus;s__Streptococcus pneumoniae_E   1
# 69       g__Streptococcus;s__Streptococcus pseudopneumoniae_E   1
# 70       g__Streptococcus;s__Streptococcus pseudopneumoniae_G   1
# 71              g__Streptococcus;s__Streptococcus shenyangsis   1
# 72              g__Streptococcus;s__Streptococcus sp000187745   1
# 73              g__Streptococcus;s__Streptococcus sp000314795   1
# 74              g__Streptococcus;s__Streptococcus sp902373455   1
# 75              g__Streptococcus;s__Streptococcus sp943912555   1

isolates %>% 
filter(`Retained for Analysis*` == "Yes") %>%
mutate(Genus = gsub(".*g__(.*);.*", "\\1", Taxonomy.GTDB.r214)) %>%
group_by(Genus) %>%
summarise(n_isolates = n(), n_unique_participants = n_distinct(Subject.ID)) %>%
arrange(desc(n_isolates)) %>%
as.data.frame()

#             Genus n_isolates n_unique_participants
# 1  Staphylococcus        284                    47
# 2   Streptococcus        207                    38
# 3       Moraxella        184                    28
# 4 Corynebacterium        159                    36
# 5  Dolosigranulum         58                    22
# 6          Rothia         33                    20
```
