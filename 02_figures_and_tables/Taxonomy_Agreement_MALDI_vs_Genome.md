# Taxonomy Agreement MALDI vs Genome

Download https://data.gtdb.ecogenomic.org/releases/release214/214.0/bac120_metadata_r214.tar.gz to extract former NCBI names for GTDB species.

```R
# Load necessary libraries
library(tidyverse)

# Read TableS1.xlsx
# Ignore first row which contains variable descriptions
table1 <- readxl::read_xlsx("/CAS_V3/Supplementary_Tables/TableS1.xlsx", skip = 1)

# Remove genomes that were filtered out due to low quality or contamination
table1 <- table1 %>% filter(`Retained for Analysis` == "Yes")
dim(table1) # 925 assembled genomes

# How many distinct species?
table1 %>% 
mutate(Species = gsub(".*s__(.*)", "\\1", `Taxonomy.GTDB.r214`)) %>%
filter(Species != "") %>%
pull(Species) %>%
unique() %>%
length() # 72 distinct species

# How many distinct genera?
table1 %>% mutate(Genus = gsub(".*g__(.*)s__.*", "\\1", `Taxonomy.GTDB.r214`)) %>% distinct(Genus) %>% nrow() # 6 distinct genera

#############################################
# Table: comparison of MALDI-TOF vs GTDB species (considering former NCBI names)
#############################################

# Table for MALDI-TOF identification agreement at species level
spp_agreement <- table1 %>%
  mutate(GTDB_Species = gsub(".*s__(.*)", "\\1", `Taxonomy.GTDB.r214`)) %>%
  filter(GTDB_Species != "") %>%  # Remove entries with species not assigned in GTDB
  select(Isolate, GTDB_Species, MALDI.TOF.Identification)

dim(spp_agreement) # 902 genomes with species assigned in GTDB

# Correct small typos in MALDI-TOF identifications based on GTDB species names
spp_agreement %>%
distinct(GTDB_Species, MALDI.TOF.Identification) %>%
arrange(GTDB_Species) %>%
as.data.frame()
#                             GTDB_Species            MALDI.TOF.Identification
# 6   Corynebacterium pseudodiphtheriticum Corynebacterium pseudodiptheriticum
# 42       Staphylococcus pseudintermedius    Staphylococcus pseudointermedius

# Harmonize these names
spp_agreement <- spp_agreement %>%
  mutate(MALDI.TOF.Identification = case_when(
    MALDI.TOF.Identification == "Corynebacterium pseudodiptheriticum" ~ "Corynebacterium pseudodiphtheriticum",
    MALDI.TOF.Identification == "Staphylococcus pseudointermedius" ~ "Staphylococcus pseudintermedius",
    TRUE ~ MALDI.TOF.Identification
  )) 
  
# TEST AGREEMENT CONSIDERING ORIGINAL GTDB SPECIES NAMES (with suffixes)
spp_agreement <- spp_agreement %>% mutate(Agreement_GTDB_original_names = ifelse(MALDI.TOF.Identification == GTDB_Species, "Yes", "No"))
spp_agreement$Agreement_GTDB_original_names %>% table() # YES = 719

# TEST AGREEMENT CONSIDERING GTDB SPECIES NAMES WITHOUT SUFFIXES
spp_agreement$GTDB_Species_no_suffixes <- sub("_[A-Z]+$", "", spp_agreement$GTDB_Species) # Remove alphabetic suffixes (e.g. "_A", "_B" etc) 
spp_agreement$Agreement_GTDB_no_suffixes <- ifelse(spp_agreement$MALDI.TOF.Identification == spp_agreement$GTDB_Species_no_suffixes, "Yes", "No")
spp_agreement$Agreement_GTDB_no_suffixes %>% table() # YES = 755

# TEST AGREEMENT CONSIDERING FORMER NCBI NAMES FROM GTDB METADATA
bac120 <- read.delim("bac120_metadata_r214.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>% 
          mutate(GTDB_Species = gsub(".*s__(.*)", "\\1", gtdb_taxonomy)) %>%
          mutate(ncbi_taxonomy = gsub(".*s__(.*)", "\\1", ncbi_taxonomy)) %>%
          select(GTDB_Species, ncbi_taxonomy) %>%
          filter(GTDB_Species %in% spp_agreement$GTDB_Species) %>%
          distinct(GTDB_Species, ncbi_taxonomy)

bac120$GTDB_Species %>% unique() %>% length() # 72 distinct GTDB species
dim(bac120) # 144 distinct NCBI names for 72 GTDB species... I will need to concatenate them

bac120_collapsed <- bac120 %>%
  filter(ncbi_taxonomy != "") %>%
  group_by(GTDB_Species) %>%
  summarise(ncbi_taxonomy = paste(unique(ncbi_taxonomy), collapse = "; ")) %>%
  ungroup()

# Add former NCBI names to spp_agreement
spp_agreement <- spp_agreement %>% left_join(bac120_collapsed, by = "GTDB_Species")

# Save Table 
# I will check former NCBI names manually to make sure we are not considering misidentifications that are actually just nomenclature updates
write.csv(spp_agreement, file = "TableSB_MALDI_vs_GTDB_species.csv", row.names = FALSE)
```
