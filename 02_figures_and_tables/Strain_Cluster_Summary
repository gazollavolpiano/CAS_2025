# Create Table S7: Summary Data of Strain Clusters Detected 

The columns will be:

  Species: The species being compared.
  Strain Cluster: The detected strain cluster (e.g., sau_C127), considering data from both isolates and metagenomes from plate sweeps.
  Existence Evidence: The evidence supporting the existence of the strain cluster (e.g., CAS Isolate, CAS Isolate & Metagenomic Detection, Metagenomic Detection Based on Public Genomes Only).
  Sequence Types in the Cluster (CAS only): The sequence types detected within the strain cluster, considering only the CAS genome profiles.
  Clonal Complexes in the Cluster (CAS only): The clonal complexes detected within the strain cluster, considering only the CAS genome profiles.
  Typing Results in the Cluster (CAS only): The typing results detected within the strain cluster, considering only the CAS genome profiles.
  Prevalence (Number of Colonized Infants):
      Adjusted for Multiple Cultivation Media.
      Adjusted for Multiple Cultivation Media and Longitudinal Sampling.

The available data to create the table can be found in the following sources:
  - `/CAS_2025/01_tableS1/01_tableS1.csv`: it is a summary data for isolates, also including infant age and log number
  - `/CAS_2025/04_figure1_tableS4_table2/04_tableS4.csv`: has data for the metagenomes available, will be use to get the Subject.IDs not present in the isolates data
  - `/2025/01_CAS/strain_scan/profile_plate_sweep/strain_scan_results_plate_sweep_11_databases.csv`: the strain clusters detected in the metagenomes of the infants
  - `/CAS_2025/06_figure2_tableS5_tableS6/06_tableS5.csv`: data for classification of genomes into strain clusters 
  - `/CAS_2025/02_tableS2/02_tableS2.csv`: data for classification of genomes into strain clusters 

```R
# library
library(tidyverse)

# increase the stdout limit
options(width = 200)

##############################
# Create a basic data frame with data from isolates and metagenomes
##############################

# read the data for the classification of genomes into strain clusters
table_s5 <- read.csv("/CAS_2025/06_figure2_tableS5_tableS6/06_tableS5.csv")

# list the species of interest
species_list <- table_s5$Species %>% unique()
species_list # correct list of species

# read summary data for isolates and filter for genomes of interest
isolates_data <- read.csv("/CAS_2025/01_tableS1/01_tableS1.csv") %>%
                    mutate(Species = gsub(".*;s__", "", Taxonomy.GTDB.r214)) %>% # remove the prefix
                    mutate(Species = gsub("Streptococcus mitis.*", "Streptococcus mitis", Species)) %>% # collapse all S. mitis species
                    filter(Species %in% species_list) %>% # filter for species of interest
                    select(Species, Isolate, Subject.ID, Sampling.Age) 
dim(isolates_data) # 750x4 --> correct number of CAS isolates for the 11 species

# add the Cluster ID to the isolates data
strains <- table_s5 %>%
                    separate_rows(Genomes, sep = ",") %>%
                    mutate(Genomes = gsub("filt_", "", Genomes)) %>%
                    select(Species, Genomes, Cluster.ID) 
setdiff(isolates_data$Isolate, strains$Genomes) # no missing genomes from CAS

isolates_data <- left_join(isolates_data, strains, by = c("Isolate" = "Genomes", "Species" = "Species")) %>% select(Species, Cluster.ID, Subject.ID, Sampling.Age)
isolates_data$Existence.Evidence <- "CAS Isolate"
head(isolates_data)
#                 Species Cluster.ID Subject.ID Sampling.Age Existence.Evidence
# 1 Staphylococcus aureus        C21   0003THCO            2        CAS Isolate
# 2 Staphylococcus aureus        C21   0003THCO            2        CAS Isolate
# 3 Staphylococcus aureus        C21   0003THCO            2        CAS Isolate
# 4 Staphylococcus aureus        C21   0003THCO            2        CAS Isolate
# 5 Staphylococcus aureus        C42   0007EVQU            2        CAS Isolate
# 6    Rothia sp902373285        C10   0007EVQU            2        CAS Isolate

# organize the data from plate sweeps metagenomes processed with strainscan and tidy it
strains_sweeps <- read.csv("/2025/01_CAS/strain_scan/profile_plate_sweep/strain_scan_results_plate_sweep_11_databases.csv", header = TRUE, stringsAsFactors = FALSE) %>%
                    filter(!database == "mitis_constrained") %>% # remove poorly performing database
                    mutate(Cluster.ID = as.character(cluster_ID)) %>%
                    mutate(metagenomic_sample = gsub("_plate_sweep", "\\1", metagenomic_sample)) %>%
                    mutate(Species = case_when(
                      database == "accolens" ~ "Corynebacterium accolens",
                      database == "agalactiae" ~ "Streptococcus agalactiae",
                      database == "aureus" ~ "Staphylococcus aureus",
                      database == "catarrhalis" ~ "Moraxella catarrhalis",
                      database == "mitis_relaxed" ~ "Streptococcus mitis",
                      database == "nonliquefaciens" ~ "Moraxella nonliquefaciens",
                      database == "pigrum" ~ "Dolosigranulum pigrum",
                      database == "pneumoniae" ~ "Streptococcus pneumoniae",
                      database == "propinquum" ~ "Corynebacterium propinquum",
                      database == "pseudodiphtheriticum" ~ "Corynebacterium pseudodiphtheriticum",
                      database == "rothia" ~ "Rothia sp902373285",
                      TRUE ~ database
                    )) %>%
                    select(metagenomic_sample, Species, Cluster.ID) %>%
                    filter(metagenomic_sample != "CAS_2152") # remove the plate that cannot be linked to age

## add columns for Infant Log Number and Sampling Age (months) to strains_sweeps
log_age_subjectid <- read.csv("/2024/12_CAS/4_gtdbtk/log_age_subjectid.csv") %>%
                      rename(Subject.ID = subject.id, Sampling.Age = age_months, Sample.Log.Number = log.number) %>%
                      mutate(Sampling.Age = case_when(Sampling.Age %in% c(1, 2, 3) ~ "2", Sampling.Age %in% c(6, 7) ~ "6", Sampling.Age %in% c(11, 12, 13) ~ "12", TRUE ~ as.character(Sampling.Age))) %>%
                      mutate(Sampling.Age = as.numeric(Sampling.Age)) 
strains_sweeps$Sample.Log.Number <- as.numeric(gsub("CAS_", "", strains_sweeps$metagenomic_sample))
strains_sweeps <- left_join(strains_sweeps, log_age_subjectid, by = "Sample.Log.Number") %>% select(Species, Cluster.ID, Subject.ID, Sampling.Age)
strains_sweeps$Existence.Evidence <- "Metagenome"

head(strains_sweeps)
#                    Species Cluster.ID Subject.ID Sampling.Age Existence.Evidence
# 1 Corynebacterium accolens         C9   0039SABI            2         Metagenome
# 2 Corynebacterium accolens        C17   0003THCO            6         Metagenome
# 3 Corynebacterium accolens        C14   0017CAST            2         Metagenome
# 4 Corynebacterium accolens        C17   0017CAST            2         Metagenome
# 5 Corynebacterium accolens         C7   0078SAWA            2         Metagenome
# 6 Corynebacterium accolens         C8   0137ELVA            2         Metagenome

# union the data
union_data <- rbind(isolates_data, strains_sweeps)

##############################
# Calculate Prevalence (Number of Colonized Infants):
#         Adjusted for Multiple Cultivation Media
#         Adjusted for Multiple Cultivation Media and Longitudinal Sampling
##############################

# Adjusted for Multiple Cultivation Media
tmp1 <- union_data %>% 
  select(Species, Cluster.ID, Subject.ID, Sampling.Age) %>%
  distinct() %>%
  group_by(Species, Cluster.ID) %>%
  summarise(Prevalence.CultMedia = n()) 

# Adjusted for Multiple Cultivation Media and Longitudinal Sampling
tmp2 <- union_data %>% 
  select(Species, Cluster.ID, Subject.ID) %>%
  distinct() %>%
  group_by(Species, Cluster.ID) %>%
  summarise(Prevalence.CultMedia.LongSamp = n()) 

# join the data
Prevalence <- left_join(tmp1, tmp2, by = c("Species", "Cluster.ID")) 

##############################
# Calculate Existence Evidence 
##############################

Existence <- union_data %>% 
  select(Species, Cluster.ID, Existence.Evidence) %>%
  distinct() %>%
  group_by(Species, Cluster.ID) %>%
  summarise(Existence.Evidence = paste(Existence.Evidence, collapse = ", "))

Existence$Existence.Evidence <- ifelse(Existence$Existence.Evidence == "Metagenome", "Metagenomic Detection Based on Public Genomes Only", Existence$Existence.Evidence)
Existence$Existence.Evidence <- ifelse(Existence$Existence.Evidence == "CAS Isolate, Metagenome", "CAS Isolate & Metagenomic Detection", Existence$Existence.Evidence)

table(Existence$Existence.Evidence)
#  CAS Isolate                CAS Isolate & Metagenomic Detection Metagenomic Detection Based on Public Genomes Only 
#           10                                                136                                                 88 

# add the Existence Evidence to the Prevalence data
Prevalence <- left_join(Prevalence, Existence, by = c("Species", "Cluster.ID"))

##############################
# Add Sequence Types, Clonal Complexes and Typing Results (considering only the CAS genome profiles)
##############################

# read 02_tableS2.csv and filter for CAS genomes
table_s2 <- read.csv("/CAS_2025/02_tableS2/02_tableS2.csv") %>%
            filter(Genome.Source == "CAS")          # filter for CAS genomes  

# add the Cluster.ID
table_s2 <- left_join(table_s2, strains, by = c("Genome" = "Genomes", "Species" = "Species")) 

# for each Cluster.ID, get the Sequence Types within the cluster
Sequence_Types <- table_s2 %>%
  select(Species, Cluster.ID, ST) %>%
  distinct() %>%
  group_by(Species, Cluster.ID) %>%
  summarise(Sequence.Types = paste(ST, collapse = ", "))
Sequence_Types$Sequence.Types <- gsub("Scheme not available", "", Sequence_Types$Sequence.Types)

# for each Cluster.ID, get the Clonal Complexes within the cluster
Clonal_Complexes <- table_s2 %>%
  select(Species, Cluster.ID, Clonal.Complex) %>%
  distinct() %>%
  group_by(Species, Cluster.ID) %>%
  summarise(Clonal.Complexes = paste(Clonal.Complex, collapse = ", "))

Clonal_Complexes$Clonal.Complexes <- gsub(", Clonal Complex not available", "", Clonal_Complexes$Clonal.Complexes)
Clonal_Complexes$Clonal.Complexes <- gsub("Clonal Complex not available", "", Clonal_Complexes$Clonal.Complexes)

# for each Cluster.ID, get the Typing Results within the cluster for agr.Group
agr.Group <- table_s2 %>%
  select(Species, Cluster.ID, agr.Group) %>%
  distinct() %>%
  group_by(Species, Cluster.ID) %>%
  summarise(agr.Group = paste(agr.Group, collapse = ", ")) %>%
  mutate(agr.Group = paste0("agr group = ", agr.Group))
  
agr.Group$agr.Group <- gsub("^agr group = $", "", agr.Group$agr.Group)

# for each Cluster.ID, get the Typing Results within the cluster for spa.Type
spa.Type <- table_s2 %>%
  select(Species, Cluster.ID, spa.Type) %>%
  distinct() %>%
  group_by(Species, Cluster.ID) %>%
  summarise(spa.Type = paste(spa.Type, collapse = ", ")) %>%
  mutate(spa.Type = paste0("spa type = ", spa.Type))

spa.Type$spa.Type <- gsub("^spa type = $", "", spa.Type$spa.Type)

# for each Cluster.ID, get the Typing Results for Capsular.Serotype
Capsular_Serotype <- table_s2 %>%
  select(Species, Cluster.ID, Capsular.Serotype) %>%
  distinct() %>%
  group_by(Species, Cluster.ID) %>%
  mutate(Capsular.Serotype = ifelse(Capsular.Serotype == "Non-typable (absence of a fully functional capsular locus)", "", Capsular.Serotype)) %>%
  summarise(Capsular.Serotype = paste(Capsular.Serotype, collapse = ", ")) %>%
  mutate(Capsular.Serotype = paste0("capsular serotype = ", Capsular.Serotype))

Capsular_Serotype$Capsular.Serotype <- gsub("^capsular serotype = $", "", Capsular_Serotype$Capsular.Serotype)
Capsular_Serotype$Capsular.Serotype <- gsub("= , ", "= ", Capsular_Serotype$Capsular.Serotype)
Capsular_Serotype$Capsular.Serotype <- gsub(", $", "", Capsular_Serotype$Capsular.Serotype)

# merge the data
Prevalence <- left_join(Prevalence, Sequence_Types, by = c("Species", "Cluster.ID"))
Prevalence <- left_join(Prevalence, Clonal_Complexes, by = c("Species", "Cluster.ID"))
Prevalence <- left_join(Prevalence, agr.Group, by = c("Species", "Cluster.ID"))
Prevalence <- left_join(Prevalence, spa.Type, by = c("Species", "Cluster.ID"))
Prevalence <- left_join(Prevalence, Capsular_Serotype, by = c("Species", "Cluster.ID"))
## replace NA with ""
Prevalence[is.na(Prevalence)] <- ""
Prevalence <- Prevalence %>% mutate(Typing.Results = paste(agr.Group, spa.Type, Capsular.Serotype, sep = "; ")) %>% select(-agr.Group, -spa.Type, -Capsular.Serotype)
Prevalence$Typing.Results  <- gsub("^; ; $", "", Prevalence$Typing.Results)
Prevalence$Typing.Results  <- gsub("; ; ", "", Prevalence$Typing.Results)
Prevalence$Typing.Results  <- gsub("; $", "", Prevalence$Typing.Results)

# tidy the data
Prevalence <- Prevalence %>% select(Species, Cluster.ID, Existence.Evidence, Sequence.Types, Clonal.Complexes, Typing.Results, Prevalence.CultMedia, Prevalence.CultMedia.LongSamp) %>% arrange(Species, Prevalence.CultMedia)

# write the table
write.csv(Prevalence, "08_tableS7.csv", row.names = FALSE)
```
