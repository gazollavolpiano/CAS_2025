# Create Table S3: Virulence Factors and Antimicrobial Resistance Genes 

We will add only CAS genomes that were used in the StrainScan analysis. 

The available data to create the table can be found in the following sources:
  - `/12_CAS/8_abricate/misc/CAS_genomes_resfinder.tab`: results of the resistance genes search
  - `/12_CAS/8_abricate/misc/CAS_genomes_vfdb.tab`: results of the virulence factors search
  - `/12_CAS/8_abricate/misc/VF_metadata.csv`: metadata for the virulence factors

```R
# call library
library(tidyverse)

##############################
# Create column Genome and Species
# (only for CAS genomes used in the StrainScan analysis)
##############################

# create a list with the public genomes that were used in the StrainScan analysis
databases <- c("accolens", "agalactiae", "aureus", "catarrhalis", "mitis_relaxed", "nonliquefaciens", "pigrum", "pneumoniae", "propinquum", "pseudodiphtheriticum", "rothia")

genomes_df <- data.frame()
for(d in databases) {
  tmp <- list.files(paste0("/01_CAS/strain_scan/database/",d,"/genomes"), pattern = ".fasta|.fna", full.names = TRUE) %>%
          map_chr(basename) %>%
          gsub("filt_(CAS_.*).fasta", "\\1", .) 
  tmp <- data.frame(Genome = tmp, Species = d)
  tmp <- tmp %>% filter(grepl("CAS", Genome)) # only CAS genomes
  genomes_df <- rbind(genomes_df, tmp)
}

# modify the Species column to give full names
replacements <- c("aureus"= "Staphylococcus aureus", "pneumoniae"= "Streptococcus pneumoniae", "agalactiae"= "Streptococcus agalactiae", "catarrhalis"= "Moraxella catarrhalis",
  "pigrum"= "Dolosigranulum pigrum", "accolens"= "Corynebacterium accolens", "propinquum"= "Corynebacterium propinquum", "pseudodiphtheriticum"= "Corynebacterium pseudodiphtheriticum",
  "rothia"= "Rothia sp902373285", "nonliquefaciens"= "Moraxella nonliquefaciens", "mitis_relaxed"= "Streptococcus mitis")
genomes_df$Species <- str_replace_all(genomes_df$Species, replacements)

table(genomes_df$Species)
#             Corynebacterium accolens           Corynebacterium propinquum 
#                                   26                                   24 
# Corynebacterium pseudodiphtheriticum                Dolosigranulum pigrum 
#                                   96                                   58 
#                Moraxella catarrhalis            Moraxella nonliquefaciens 
#                                  165                                   19 
#                   Rothia sp902373285                Staphylococcus aureus 
#                                   22                                  231 
#             Streptococcus agalactiae                  Streptococcus mitis 
#                                   23                                   40 
#             Streptococcus pneumoniae 
#                                   46 

##############################
# Process resfinder and vfdb results
##############################

# resfinder processing
resfinder <- read.delim("/12_CAS/8_abricate/misc/CAS_genomes_resfinder.tab", header=T) %>%
            mutate(Genome = gsub(".*_(CAS_.*).fasta", "\\1", X.FILE), Database = "ResFinder") %>%
            rename(Gene = PRODUCT, Coverage = X.COVERAGE, Identity = X.IDENTITY, Observation = RESISTANCE) %>%
            filter(Genome %in% genomes_df$Genome) %>%
            select(Genome, Database, Gene, Coverage, Identity, Observation) %>%
            mutate(Observation = paste0("Resistance to ", tolower(Observation)))

resfinder$Genome %>% unique() %>% length() # 509 genomes

# replace empty "Resistance to " observations with ""
resfinder <- resfinder %>% mutate(Observation = ifelse(Observation == "Resistance to ", "", Observation))

# vfdb processing
vfdb_metadata <- read.csv("/12_CAS/8_abricate/misc/VF_metadata.csv", header=T) %>% select(VFID, VF_Name, VFcategory)

vfdb <- read.delim("/12_CAS/8_abricate/misc/CAS_genomes_vfdb.tab", header=T) %>%
            # remove CANDIDATE VIRULENCE FACTORS (keep only virulence factors with a VFID)
            mutate(VFID = gsub(".*\\((VF\\d+).*", "\\1", PRODUCT)) %>%
            filter(grepl("^VF", VFID)) %>%
            mutate(Genome = gsub(".*_(CAS_.*).fasta", "\\1", X.FILE), Database = "VFDB") %>%
            rename(Gene = GENE, Coverage = X.COVERAGE, Identity = X.IDENTITY) %>%
            filter(Genome %in% genomes_df$Genome) %>%
            select(Genome, Database, Gene, Coverage, Identity, VFID) %>%
            left_join(vfdb_metadata) %>%
            mutate(Observation = paste0(VFID," - ",VF_Name, " (", VFcategory, ")")) %>%
            select(Genome, Database, Gene, Coverage, Identity, Observation)

vfdb$Genome %>% unique() %>% length() # 340 genomes

# combine resfinder and vfdb
resfinder_vfdb <- rbind(resfinder, vfdb)

# join with genomes_df
genomes_df <- left_join(genomes_df, resfinder_vfdb, by = "Genome")

# replace NA with ""
genomes_df[is.na(genomes_df)] <- ""

# write the table
write.csv(genomes_df, "03_tableS3.csv", row.names = FALSE)
```
