# Create Table S2: Typing Data for CAS and Public Genomes 

Use CAS genomes and public genomes that were used in the StrainScan analysis. 

The available data to create the table can be found in the following sources:
  - `/12_CAS/6_ST/saureus_231_CAS_14959_public_mlst_corrected.csv`: MLST for S. aureus
  - `/12_CAS/6_ST/spneumoniae_46_CAS_8895_public_mlst_corrected.csv`: MLST for S. pneumoniae
  - `/12_CAS/6_ST/sagalactiae_fasta_mlst_corrected.csv`: MLST for S. agalactiae 
  - `/12_CAS/6_ST/mcatarrhalis_165_CAS_208_public_mlst_corrected.csv`: MLST for M. catarrhalis 
  - `/12_CAS/6_ST/aug_12_2024_achtman_saureus_cc_profiles.csv`: CC for S. aureus
  - `/12_CAS/6_ST/aug_12_2024_achtman_mcatarrhalis_cc_profiles.csv`: CC for M. catarrhalis
  - `/12_CAS/6_ST/feb_28_2025_sagalactiae_cc_profiles.csv`: CC for S. agalactiae
  - `/12_CAS/6_ST/staphopia_CAS.csv`: SCCmec for S. aureus - CAS
  - `/12_CAS/6_ST/staphopia_public.csv`: SCCmec for S. aureus - public genomes
  - `/12_CAS/6_ST/spaTyper`: spa types for S. aureus, directory with multiple files
  - `/12_CAS/6_ST/agrvate`: agr groups for S. aureus, directory with multiple files
  - `/12_CAS/6_ST/pneumocat`: Capsular typing for S. pneumoniae, directory with multiple files --> CAS only because it needs the raw reads
  - `/01_CAS/strain_scan/database`: StrainScan genome databases (accolens, agalactiae, aureus, catarrhalis, mitis_relaxed, nonliquefaciens, pigrum, pneumoniae, propinquum, pseudodiphtheriticum, rothia)

```R
# call library
library(tidyverse)

##############################
# Create column Genome, Species and Genome Source
# (only for CAS and public genomes used in the StrainScan analysis)
##############################

# create a list with the public genomes that were used in the StrainScan analysis
databases <- c("accolens", "agalactiae", "aureus", "catarrhalis", "mitis_relaxed", "nonliquefaciens", "pigrum", "pneumoniae", "propinquum", "pseudodiphtheriticum", "rothia")

genomes_df <- data.frame()
for(d in databases) {
  tmp <- list.files(paste0("/01_CAS/strain_scan/database/",d,"/genomes"), pattern = ".fasta|.fna", full.names = TRUE) %>%
          map_chr(basename) %>%
          gsub("filt_(CAS_.*).fasta", "\\1", .) %>%
          gsub("(GCA_...........).*", "\\1", .)
  tmp <- data.frame(Genome = tmp, Species = d)
  genomes_df <- rbind(genomes_df, tmp)
}

# modify the Species column to give full names
replacements <- c("aureus"= "Staphylococcus aureus", "pneumoniae"= "Streptococcus pneumoniae", "agalactiae"= "Streptococcus agalactiae", "catarrhalis"= "Moraxella catarrhalis",
  "pigrum"= "Dolosigranulum pigrum", "accolens"= "Corynebacterium accolens", "propinquum"= "Corynebacterium propinquum", "pseudodiphtheriticum"= "Corynebacterium pseudodiphtheriticum",
  "rothia"= "Rothia sp902373285", "nonliquefaciens"= "Moraxella nonliquefaciens", "mitis_relaxed"= "Streptococcus mitis")
genomes_df$Species <- str_replace_all(genomes_df$Species, replacements)

table(genomes_df$Species)
#             Corynebacterium accolens           Corynebacterium propinquum 
#                                   42                                   40 
# Corynebacterium pseudodiphtheriticum                Dolosigranulum pigrum 
#                                  112                                   91 
#                Moraxella catarrhalis            Moraxella nonliquefaciens 
#                                  373                                   25 
#                   Rothia sp902373285                Staphylococcus aureus 
#                                   30                                 1479 
#             Streptococcus agalactiae                  Streptococcus mitis 
#                                 1658                                  230 
#             Streptococcus pneumoniae 
#                                 1480 

# add the source of the genome
genomes_df <- genomes_df %>% mutate(Genome.Source = ifelse(grepl("CAS", Genome), "CAS", "public"))

##############################
# Create column with ST
##############################

# list of 4 files with MLST data
files_to_read <- c("/12_CAS/6_ST/saureus_231_CAS_14959_public_mlst_corrected.csv",
                   "/12_CAS/6_ST/spneumoniae_46_CAS_8895_public_mlst_corrected.csv",
                   "/12_CAS/6_ST/sagalactiae_fasta_mlst_corrected.csv",
                   "/12_CAS/6_ST/mcatarrhalis_165_CAS_208_public_mlst_corrected.csv")

mlst_result <- data.frame()
for(file in files_to_read) {
  tmp <- read.csv(file, header = TRUE, sep = ",") %>%
    filter(FILE != "FILE") %>%
    mutate(Genome = gsub(".*filt_(CAS_.*).fasta", "\\1", FILE)) %>%
    mutate(Genome = gsub(".*(GCA_...........).*","\\1",Genome)) %>%
    filter(Genome %in% genomes_df$Genome) %>%
    select(Genome, ST)
  mlst_result <- rbind(mlst_result, tmp)
}

# replace "-" with "Alleles do not match any known STs"
mlst_result$ST[mlst_result$ST=="-"] <- "Alleles do not match any known STs"

# add the ST column to the genomes_df
genomes_df <- left_join(genomes_df, mlst_result, by = "Genome")

# replace NA with "Scheme not available"
genomes_df$ST[is.na(genomes_df$ST)] <- "Scheme not available"

##############################
# Create column with Clonal Complex
##############################

# read the clonal complex data
cc_files <- c("/12_CAS/6_ST/aug_12_2024_achtman_saureus_cc_profiles.csv",
              "/12_CAS/6_ST/aug_12_2024_achtman_mcatarrhalis_cc_profiles.csv",
              "/12_CAS/6_ST/feb_28_2025_sagalactiae_cc_profiles.csv")

cc_result <- data.frame()
for(file in cc_files) {
  tmp <- read.delim(file, header = TRUE) %>%
          select(ST, clonal_complex) %>%
          rename(Clonal.Complex = clonal_complex) %>%
          mutate(ST = as.character(ST), Species = gsub(".*/6_ST/(.*)","\\1",file))
  cc_result <- rbind(cc_result, tmp)
}

cc_result$Species <- gsub("aug_12_2024_achtman_saureus_cc_profiles.csv", "Staphylococcus aureus", cc_result$Species)
cc_result$Species <- gsub("aug_12_2024_achtman_mcatarrhalis_cc_profiles.csv", "Moraxella catarrhalis", cc_result$Species)
cc_result$Species <- gsub("feb_28_2025_sagalactiae_cc_profiles.csv", "Streptococcus agalactiae", cc_result$Species)

# check if all STs are present in the clonal complex data
setdiff(genomes_df %>% filter(Species == "Staphylococcus aureus") %>% pull(ST), cc_result %>% filter(Species == "Staphylococcus aureus") %>% pull(ST)) # ok, all STs are present
setdiff(genomes_df %>% filter(Species == "Moraxella catarrhalis") %>% pull(ST), cc_result %>% filter(Species == "Moraxella catarrhalis") %>% pull(ST)) # ok, all STs are present
setdiff(genomes_df %>% filter(Species == "Streptococcus agalactiae") %>% pull(ST), cc_result %>% filter(Species == "Streptococcus agalactiae") %>% pull(ST)) # ok, all STs are present

# standardize the Clonal.Complex column
cc_result$Clonal.Complex <- toupper(cc_result$Clonal.Complex)

# add the clonal complex to the genomes_df
genomes_df <- left_join(genomes_df, cc_result, by = c("ST" = "ST", "Species" = "Species"))

# replace NA with "Clonal Complex not available"
genomes_df$Clonal.Complex[is.na(genomes_df$Clonal.Complex)] <- "Clonal Complex not available"

##############################
# Create column with agr Group
##############################

# list files with agr data
directories <- list.dirs(path = "/12_CAS/6_ST/agrvate", full.names = FALSE, recursive = FALSE)

# keep only the 1479 directories data for the selected S. aureus genomes
genomes_to_grep <- genomes_df %>% filter(Species == "Staphylococcus aureus") %>% pull(Genome) 
directories <- grep(paste(genomes_to_grep, collapse="|"), directories, value = TRUE)
length(directories) # 1479

# read the agr data
agr_files <- paste0("/12_CAS/6_ST/agrvate/",directories,"/",gsub("-results","",directories),"-summary.tab")
agr <- do.call(rbind, lapply(agr_files, function(file) {
  cat(file, "\n")
  df <- read.delim(file, header = TRUE, sep = "\t", check.names = FALSE)
  df$source_file <- basename(file) 
  return(df)
}))

# tidy the data
agr <- agr %>%
  mutate(Genome = gsub(".*filt_(CAS_.*)-summary.tab", "\\1", source_file)) %>%
  mutate(Genome = gsub(".*(GCA_...........).*","\\1",Genome)) %>%
  rename(agr.Group = agr_group) %>%
  select(Genome, agr.Group)

setdiff(genomes_df %>% filter(Species == "Staphylococcus aureus") %>% pull(Genome), agr$Genome) # ok, all genomes are present

# add the agr group to the genomes_df
dim(agr) # 1479x2
genomes_df <- left_join(genomes_df, agr, by = "Genome")

# replace NA with ""
genomes_df$agr.Group[is.na(genomes_df$agr.Group)] <- ""

##############################
# Create column with spa Type
##############################

# list files with spa data
directories <- list.files(path = "/12_CAS/6_ST/spaTyper", pattern = "\\.txt$", full.names = TRUE)

# keep only the 1479 directories data for the selected S. aureus genomes
directories <- grep(paste(genomes_to_grep, collapse="|"), directories, value = TRUE)
length(directories) # 1479

# read the spa data
spa <- data.frame()
for (file in directories) {
  cat(file, "\n")
  genome_from_file <- gsub(".*filt_(CAS_.*).txt", "\\1", basename(file))
  genome_from_file <- gsub(".*(GCA_...........).*","\\1",genome_from_file)

  tmp <- read.delim(file, header = TRUE, sep = "\t", check.names = FALSE) %>% 
    mutate(Genome = genome_from_file) %>%
    rename(spa.Type = Type) %>%
    group_by(Genome) %>%
    summarise(spa.Type = paste(unique(spa.Type), collapse = ";")) %>% # collapse multiple spa types
    select(Genome, spa.Type)

  if(nrow(tmp) == 0) {
    tmp <- data.frame("Genome" = genome_from_file, "spa.Type" = "Spa type not available")
  }
  spa <- rbind(spa, tmp)
}

# add the spa type to the genomes_df
dim(spa) # 1479x2
genomes_df <- left_join(genomes_df, spa, by = "Genome")

# replace NA with ""
genomes_df$spa.Type[is.na(genomes_df$spa.Type)] <- ""

##############################
# Create column with SCCmec
##############################

# read the SCCmec data
sccmec <- read.table("/12_CAS/6_ST/staphopia_public.csv", header = FALSE) %>%
                    rbind(read.table("/12_CAS/6_ST/staphopia_CAS.csv", header = FALSE)) %>%
                    mutate(Genome = gsub("(GCA_...........).*","\\1",V1)) %>%
                    mutate(Genome = gsub(".*filt_(CAS_.*).fasta", "\\1", Genome)) %>%
                    filter(Genome %in% (genomes_df %>% filter(Species == "Staphylococcus aureus") %>% pull(Genome))) %>%
                    select(-V1)

# correct column names
colnames(sccmec)[1:20] <- c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "meca", "Ia", "IIa", "IIb", "IIIa", "IVa", "IVb", "IVc", "IVd", "IVg", "IVh")

# replace "False" with 0 and True with 1
sccmec[sccmec == "False"] <- 0
sccmec[sccmec == "True"] <- 1

# summarize the types and subtypes
summary_sccmec <- data.frame()
for (i in 1:nrow(sccmec)) {
  sccmec_row <- sccmec[i, ]
  present_types <- names(sccmec_row)[sccmec_row == 1]
  # check if mecA is present
  if ("meca" %in% present_types) {
    mecA_result <- "mecA gene found"
    present_types <- present_types[present_types != "meca"]
  } else {
    mecA_result <- "mecA gene not found"
  }
  # aggregate the types
  if (length(present_types) > 0) {
    present_types <- paste(present_types, collapse=", ")
    present_types <- paste(present_types, mecA_result, sep = "; ")
  } else {
    present_types <- "SCCmec not present"
    present_types <- paste(present_types, mecA_result, sep = "; ")
  }
  summary_sccmec <- rbind(summary_sccmec, data.frame(Genome = sccmec_row$Genome, SCCmec = present_types))
}

# add the SCCmec to the genomes_df
dim(summary_sccmec) # 1479x2
genomes_df <- left_join(genomes_df, summary_sccmec, by = "Genome")

# replace NA with ""
genomes_df$SCCmec[is.na(genomes_df$SCCmec)] <- ""

##############################
# Create column with Capsular Serotype (CAS only because it needs the raw reads)
##############################

# list files with pneumocat data
library(xml2)
pneumocat_files <- list.files(path = "/6_ST/pneumocat", pattern = "CAS_", full.names = TRUE)

# read the serotype data and produce the final result
pneumocat <- data.frame()
for (file in pneumocat_files) {
  cat(file, "\n")
  sample_name <- gsub(".*pneumocat/","",file)
  # result from coverage approach
  tmp <- read_xml(paste0(file,"/",sample_name,"_R1.results.xml"))
  serotype_value <- xml_find_first(tmp, ".//result[@type='Serotype']")
  serotype_cov <- xml_attr(serotype_value, "value")
  # result from variant approach
  variant_file <- paste0(file,"/SNP_based_serotyping/",sample_name,"_R1.results.xml")
  serotype_var <- tryCatch({
    tmp <- read_xml(variant_file)
    serotype_value <- xml_find_first(tmp, ".//result[@type='Serotype_Distinction']")
    xml_attr(serotype_value, "value")
  }, error = function(e) {
    NA  # Return NA if there's an error
  })
  result <- data.frame(Genome = sample_name, serotype_coverage_approach = serotype_cov, serotype_variant_approach = serotype_var)
  pneumocat <- rbind(pneumocat, result)
}

# the final result: the serotype_variant_approach is the final result, unless it is NA and the serotype_coverage_approach is not NA
pneumocat$final_result_is_from = ifelse(is.na(pneumocat$serotype_variant_approach), "coverage-based approach", "variant-based approach")
pneumocat$serotype_coverage_approach <- paste0( pneumocat$serotype_coverage_approach, " (coverage-based approach)")
pneumocat$serotype_variant_approach <- paste0( pneumocat$serotype_variant_approach, " (variant-based approach)")
pneumocat$Capsular.Serotype <- ifelse(pneumocat$final_result_is_from == "coverage-based approach", pneumocat$serotype_coverage_approach, pneumocat$serotype_variant_approach)
pneumocat$Capsular.Serotype <- gsub("Failed .*", "Non-typable (absence of a fully functional capsular locus)", pneumocat$Capsular.Serotype)
pneumocat <- pneumocat %>% select(Genome, Capsular.Serotype)

# add the Capsular Serotype to the genomes_df
dim(pneumocat) # 46x2
genomes_df <- left_join(genomes_df, pneumocat, by = "Genome")

# replace NA with ""
genomes_df$Capsular.Serotype[is.na(genomes_df$Capsular.Serotype)] <- ""

##############################
# Manual annotation of undefined STs using trees
##############################

# Check distribution of STs that do not match any known STs from CAS genomes
genomes_df %>% filter(ST == "Alleles do not match any known STs", Genome.Source=="CAS") %>% group_by(Species) %>% summarise(n = n())
#   Species                      n
# 1 Moraxella catarrhalis       49
# 2 Staphylococcus aureus       34
# 3 Streptococcus pneumoniae     2

# Create a table with the alleles these undefined STs present
## Staphylococcus aureus
undefined_sts_aureus <- genomes_df %>% filter(ST == "Alleles do not match any known STs", Genome.Source=="CAS", Species=="Staphylococcus aureus") %>% pull(Genome) 
alleles_file_aureus <- read.csv("/12_CAS/6_ST/saureus_231_CAS_14959_public_mlst.csv") 
to_keep_aureus <- grep(paste(undefined_sts_aureus, collapse="|"), alleles_file_aureus$FILE, value = TRUE)
alleles_file_aureus <- alleles_file_aureus %>% filter(FILE %in% to_keep_aureus) %>% mutate(Genome = gsub(".*filt_(CAS_.*).fasta", "\\1", FILE)) %>% select(-FILE)
alleles_file_aureus$SCHEME <- paste("saureus: arcC", alleles_file_aureus$arcC,                  
                                             "aroE", alleles_file_aureus$aroE, 
                                             "glpF", alleles_file_aureus$glpF, 
                                             "gmk", alleles_file_aureus$gmk, 
                                             "pta", alleles_file_aureus$pta, 
                                             "tpi", alleles_file_aureus$tpi, 
                                             "yqiL", alleles_file_aureus$yqiL) 
alleles_file_aureus <- alleles_file_aureus %>% select(SCHEME, Genome)
## Moraxella catarrhalis
undefined_sts_catarrhalis <- genomes_df %>% filter(ST == "Alleles do not match any known STs", Genome.Source=="CAS", Species=="Moraxella catarrhalis") %>% pull(Genome)
alleles_file_catarrhalis <- read.csv("/12_CAS/6_ST/3_march_mcatarrhalis_165_CAS_208_public_mlst.csv")
to_keep_catarrhalis <- grep(paste(undefined_sts_catarrhalis, collapse="|"), alleles_file_catarrhalis$FILE, value = TRUE)
alleles_file_catarrhalis <- alleles_file_catarrhalis %>% filter(FILE %in% to_keep_catarrhalis) %>% mutate(Genome = gsub(".*filt_(CAS_.*).fasta", "\\1", FILE)) %>% select(-FILE)
alleles_file_catarrhalis$SCHEME <- paste("mcatarrhalis: abcZ", alleles_file_catarrhalis$abcZ, 
                                             "adk", alleles_file_catarrhalis$adk, 
                                             "efp", alleles_file_catarrhalis$efp, 
                                             "fumC", alleles_file_catarrhalis$fumC, 
                                             "glyBeta", alleles_file_catarrhalis$glyBeta, 
                                             "mutY", alleles_file_catarrhalis$mutY, 
                                             "ppa", alleles_file_catarrhalis$ppa, 
                                             "trpE", alleles_file_catarrhalis$trpE)
alleles_file_catarrhalis <- alleles_file_catarrhalis %>% select(SCHEME, Genome)
## Streptococcus pneumoniae
undefined_sts_pneumoniae <- genomes_df %>% filter(ST == "Alleles do not match any known STs", Genome.Source=="CAS", Species=="Streptococcus pneumoniae") %>% pull(Genome)
alleles_file_pneumoniae <- read.csv("/12_CAS/6_ST/spneumoniae_46_CAS_8895_public_mlst.csv")
to_keep_pneumoniae <- grep(paste(undefined_sts_pneumoniae, collapse="|"), alleles_file_pneumoniae$FILE, value = TRUE)
alleles_file_pneumoniae <- alleles_file_pneumoniae %>% filter(FILE %in% to_keep_pneumoniae) %>% mutate(Genome = gsub(".*filt_(CAS_.*).fasta", "\\1", FILE)) %>% select(-FILE)
alleles_file_pneumoniae$SCHEME <- paste("spneumoniae: aroE", alleles_file_pneumoniae$aroE, 
                                             "gdh", alleles_file_pneumoniae$gdh, 
                                             "gki", alleles_file_pneumoniae$gki, 
                                             "recP", alleles_file_pneumoniae$recP, 
                                             "spi", alleles_file_pneumoniae$spi, 
                                             "xpt", alleles_file_pneumoniae$xpt, 
                                             "ddl", alleles_file_pneumoniae$ddl)
alleles_file_pneumoniae <- alleles_file_pneumoniae %>% select(SCHEME, Genome)
## save the three tables
alleles <- rbind(alleles_file_aureus, alleles_file_catarrhalis, alleles_file_pneumoniae)
write.csv(alleles, "alleles_undefined_CAS_STs.csv", row.names = FALSE)

# use the table to manually annotate the undefined STs using trees
alleles <- read.csv("manual_anot_ready_alleles_undefined_CAS_STs.csv")

head(alleles)
#             Genome
# 1 CAS_1537_14_CHOC
# 2 CAS_1537_17_MORA
# 3 CAS_1537_57_CHOC
# 4 CAS_1560_01_BLOO
# 5 CAS_1560_02_BLOO
# 6 CAS_1560_03_BLOO
#                                                                   SCHEME
# 1 mcatarrhalis: abcZ 2 adk 2 efp 20 fumC 3 glyBeta 2 mutY 2 ppa 2 trpE 2
# 2 mcatarrhalis: abcZ 2 adk 2 efp 20 fumC 3 glyBeta 2 mutY 2 ppa 2 trpE 2
# 3 mcatarrhalis: abcZ 2 adk 2 efp 20 fumC 3 glyBeta 2 mutY 2 ppa 2 trpE 2
# 4 mcatarrhalis: abcZ 24 adk 2 efp 2 fumC 3 glyBeta 2 mutY 2 ppa 2 trpE 2
# 5 mcatarrhalis: abcZ 24 adk 2 efp 2 fumC 3 glyBeta 2 mutY 2 ppa 2 trpE 2
# 6 mcatarrhalis: abcZ 24 adk 2 efp 2 fumC 3 glyBeta 2 mutY 2 ppa 2 trpE 2
#         ST               Clonal.Complex
# 1 New ST 1 Clonal Complex not available
# 2 New ST 1 Clonal Complex not available
# 3 New ST 1 Clonal Complex not available
# 4 New ST 2 Clonal Complex not available
# 5 New ST 2 Clonal Complex not available
# 6 New ST 2 Clonal Complex not available

# add the manual annotation to the genomes_df only for the undefined STs listed in the table
for(g in alleles$Genome) {
  genomes_df$ST[genomes_df$Genome == g] <- alleles$ST[alleles$Genome == g]
  genomes_df$Clonal.Complex[genomes_df$Genome == g] <- alleles$Clonal.Complex[alleles$Genome == g]
}

##############################
# Tidy and save the table
##############################

# order columns
genomes_df <- genomes_df %>% select(Genome, Genome.Source, Species, ST, Clonal.Complex, agr.Group, spa.Type, SCCmec, Capsular.Serotype)

# write the table
write.csv(genomes_df, "02_tableS2.csv", row.names = FALSE)
```
