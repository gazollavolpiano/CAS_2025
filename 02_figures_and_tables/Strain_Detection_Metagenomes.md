# Create Figure 2 and Table S5 and S6

Table S5: Strain Clusters Identified for CAS and Public Genomes Using StrainScan
  Species: the species being profiled (e.g. Moraxella catarrhalis)
  Cluster ID: The cluster ID assigned by StrainScan
  Cluster Size: The number of genomes in the cluster
  Cluster With CAS: If the cluster contains CAS genomes
  Genomes: The genomes in the cluster (CAS isolate or assembly accession for public genomes)

Table S6: Comparision of Strain Clusters Detected in Plate Sweeps Metagenomes vs. Isolated Strains
  Subject ID: A unique identifier for each infant in the study
  Sampling Age (months): The time point of sampling (2, 6, or 12 months)
  Species: The species being compared 
  Strains Isolated [culture medium-Patch.Plate.Number]: The strain isolated from the same infant concatenated considering different media (e.g., sau_C21[BLOO-01, BLOO-02, CHOC-03, CHOC-34] means that the S. aureus strain C21 had two isolates (01 and 02) from BLOO and two (03 and 34) from CHOC)
	Strains Detected Metagenomes [Abundance%]: The strain clusters detected in the metagenomes of the infant and their relative abundance inside the cluster (e.g. sau_C21[100]). These abundances were rescaled from 0 to 100 after removing clusters defined exclusively by public genomes.
  Observation: If the metagenome from the plate sweep is not available

Important: only the strains detected in metagenomes that were defined with a cluster that contained CAS genome should be included, as they have stronger evidence of existence.

Figure 2: Distribution of Strain Clusters from 11 Dominant Species Across Infant Age Groups

Panel A:
  - Upper annotation: hierarchical clustering dendrogram 
  - Left annotation: infant age group bar + Microbiome Profile Groups (MPG) bar
  - x-axis: strain clusters per species
  - y-axis: metagenomic species 
  - color: relative abundance whithin each cluster

Panel B (left):
  - x-axis: age group
  - y-axis: proportion of infants with strain clusters detected

Panel B (right):
  - x-axis: age group
  - y-axis: strain cluster richness top=total, bottom=per infant

The available data to create the table can be found in the following sources:
  - `/01_CAS/strain_scan/database/accolens/accolens_strainscan_db/Tree_database/hclsMap_95_recls.txt`: the strain clusters for Corynebacterium accolens
  - `/01_CAS/strain_scan/database/agalactiae/agalactiae_strainscan_db/Tree_database/hclsMap_95_recls.txt`: the strain clusters for Streptococcus agalactiae
  - `/01_CAS/strain_scan/database/aureus/aureus_strainscan_db/Tree_database/hclsMap_95_recls.txt`: the strain clusters for Staphylococcus aureus
  - `/01_CAS/strain_scan/database/catarrhalis/catarrhalis_strainscan_db/Tree_database/hclsMap_95_recls.txt`: the strain clusters for Moraxella catarrhalis
  - `/01_CAS/strain_scan/database/mitis_relaxed/mitis_relaxed_strainscan_db/Tree_database/hclsMap_95_recls.txt`: the strain clusters for Streptococcus mitis (relaxed)
  - `/01_CAS/strain_scan/database/nonliquefaciens/nonliquefaciens_strainscan_db/Tree_database/hclsMap_95_recls.txt`: the strain clusters for Moraxella nonliquefaciens
  - `/01_CAS/strain_scan/database/pigrum/pigrum_strainscan_db/Tree_database/hclsMap_95_recls.txt`: the strain clusters for Dolosigranulum pigrum
  - `/01_CAS/strain_scan/database/pneumoniae/pneumoniae_strainscan_db/Tree_database/hclsMap_95_recls.txt`: the strain clusters for Streptococcus pneumoniae
  - `/01_CAS/strain_scan/database/propinquum/propinquum_strainscan_db/Tree_database/hclsMap_95_recls.txt`: the strain clusters for Corynebacterium propinquum
  - `/01_CAS/strain_scan/database/pseudodiphtheriticum/pseudodiphtheriticum_strainscan_db/Tree_database/hclsMap_95_recls.txt`: the strain clusters for Corynebacterium pseudodiphtheriticum
  - `/01_CAS/strain_scan/database/rothia/rothia_strainscan_db/Tree_database/hclsMap_95_recls.txt`: the strain clusters for Rothia sp902373285
  - `/CAS_2025/01_tableS1/01_tableS1.csv`: it is a summary data for isolates, also including infant age and log number
  - `/CAS_2025/04_figure1_tableS4_table2/04_tableS4.csv`: has data for the metagenomes available, will be use to get the Subject.IDs not present in the isolates data
  - `/01_CAS/strain_scan/profile_plate_sweep/strain_scan_results_plate_sweep_11_databases.csv`: the strain clusters detected in the metagenomes of the infants

```R
# library
library(tidyverse)
library(ComplexHeatmap)

##############################
# Create Table S5
##############################

# list of species
dabase_names <- c("accolens", "agalactiae", "aureus", "catarrhalis", "mitis_relaxed", "nonliquefaciens", "pigrum", "pneumoniae", "propinquum", "pseudodiphtheriticum", "rothia")

table_s5 <- data.frame()
for(i in dabase_names){
 # read the database file and tidy it
 staindb <- read.table(paste0("/01_CAS/strain_scan/database/", i, "/", i, "_strainscan_db/Cluster_Result/hclsMap_95_recls.txt"), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
 colnames(staindb) <- c("Cluster.ID", "Cluster.Size", "Genomes") 
 staindb$Cluster.ID <- paste0("C", as.character(staindb$Cluster.ID))
 # column indicating if CAS genomes are present
 staindb$Cluster.With.CAS <- grepl("CAS", staindb$Genomes)
 # add species column (we will tidy this later)
  staindb$Species <- i
  # bind to the table
  table_s5 <- rbind(table_s5, staindb)
}

# tidy the species column
table_s5$Species <- case_when(
  table_s5$Species == "accolens" ~ "Corynebacterium accolens",
  table_s5$Species == "agalactiae" ~ "Streptococcus agalactiae",
  table_s5$Species == "aureus" ~ "Staphylococcus aureus",
  table_s5$Species == "catarrhalis" ~ "Moraxella catarrhalis",
  table_s5$Species == "mitis_relaxed" ~ "Streptococcus mitis",
  table_s5$Species == "nonliquefaciens" ~ "Moraxella nonliquefaciens",
  table_s5$Species == "pigrum" ~ "Dolosigranulum pigrum",
  table_s5$Species == "pneumoniae" ~ "Streptococcus pneumoniae",
  table_s5$Species == "propinquum" ~ "Corynebacterium propinquum",
  table_s5$Species == "pseudodiphtheriticum" ~ "Corynebacterium pseudodiphtheriticum",
  table_s5$Species == "rothia" ~ "Rothia sp902373285",
  TRUE ~ table_s5$Species
)

# reorder the columns
table_s5 <- table_s5 %>% select(Species, Cluster.ID, Cluster.Size, Cluster.With.CAS, Genomes)

# save the table
write.csv(table_s5, "06_tableS5.csv", row.names = FALSE)

##############################
# Create Table S6
##############################

# read summary data for isolates and filter for genomes of interest
isolates_data <- read.csv("/CAS_2025/01_tableS1/01_tableS1.csv") %>%
                    mutate(Species = gsub(".*;s__", "", Taxonomy.GTDB.r214)) %>% # remove the prefix
                    mutate(Species = gsub("Streptococcus mitis.*", "Streptococcus mitis", Species)) %>% # collapse all S. mitis species
                    filter(Species %in% table_s5$Species) %>% # filter for species of interest
                    select(Species, Isolate, Culture.Medium, Subject.ID, Sample.Log.Number, Sampling.Age) 
dim(isolates_data) # 750x6 --> correct number of CAS isolates for the 11 species

# add the Cluster ID to the isolates data
strains_isolates <- table_s5 %>%
                    separate_rows(Genomes, sep = ",") %>%
                    mutate(Genomes = gsub("filt_", "", Genomes)) %>%
                    filter(grepl("CAS", Genomes)) %>% # filter for CAS genomes
                    select(Species, Genomes, Cluster.ID) 
setdiff(isolates_data$Isolate, strains_isolates$Genomes) # no missing genomes
isolates_data <- left_join(isolates_data, strains_isolates, by = c("Isolate" = "Genomes", "Species" = "Species"))

# add a prefix to the Cluster ID to indicate the species
isolates_data <- isolates_data %>%
                    mutate(Cluster.ID = as.character(Cluster.ID)) %>%
                    mutate(Cluster.ID = case_when(
                      Species == "Corynebacterium accolens" ~ paste0("cac_", Cluster.ID),
                      Species == "Streptococcus agalactiae" ~ paste0("sag_", Cluster.ID),
                      Species == "Staphylococcus aureus" ~ paste0("sau_", Cluster.ID),
                      Species == "Moraxella catarrhalis" ~ paste0("mca_", Cluster.ID),
                      Species == "Streptococcus mitis" ~ paste0("smi_", Cluster.ID),
                      Species == "Moraxella nonliquefaciens" ~ paste0("mno_", Cluster.ID),
                      Species == "Dolosigranulum pigrum" ~ paste0("dpi_", Cluster.ID),
                      Species == "Streptococcus pneumoniae" ~ paste0("spo_", Cluster.ID),
                      Species == "Corynebacterium propinquum" ~ paste0("cpr_", Cluster.ID),
                      Species == "Corynebacterium pseudodiphtheriticum" ~ paste0("cps_", Cluster.ID),
                      Species == "Rothia sp902373285" ~ paste0("rsp_", Cluster.ID),
                      TRUE ~ Cluster.ID
                    ))

# summarize the table to get the desired format for Strains.Isolated
isolates_summary <- isolates_data %>%
                    mutate(iso_num = str_split(Isolate, "_", simplify = TRUE)[,3], # extract the number from the Isolate column (the third element after splitting by "_")
                    element = paste0(Culture.Medium, "-", iso_num)) %>% # create the individual element (e.g. "BLOO-01")
                    group_by(Species, Subject.ID, Sampling.Age, Cluster.ID) %>% # group by the desired columns
                    reframe(Strains.Isolated = paste0(Cluster.ID, "[", paste(element, collapse = ","), "]"), .groups = 'drop') %>% # concatenate the elements, prefixing with the Cluster.ID and enclosing in square brackets
                    select(Species, Subject.ID, Sampling.Age, Strains.Isolated) %>% # select the columns in the desired order
                    distinct() # remove duplicated rows

# add the Strains.Detected.Metagenomes column
## read data from plate sweeps metagenomes processed with strainscan and tidy it
strains_sweeps <- read.csv("/01_CAS/strain_scan/profile_plate_sweep/strain_scan_results_plate_sweep_11_databases.csv", header = TRUE, stringsAsFactors = FALSE) %>%
                    filter(CAS_genome_included_cluster == "TRUE") %>% # filter for clusters with CAS genomes (stronger evidence of existence)
                    filter(!database == "mitis_constrained") %>% # remove poorly performing database
                    mutate(Cluster.ID = as.character(cluster_ID)) %>%
                    mutate(Cluster.ID = case_when(
                      database == "accolens" ~ paste0("cac_", Cluster.ID),
                      database == "agalactiae" ~ paste0("sag_", Cluster.ID),
                      database == "aureus" ~ paste0("sau_", Cluster.ID),
                      database == "catarrhalis" ~ paste0("mca_", Cluster.ID),
                      database == "mitis_relaxed" ~ paste0("smi_", Cluster.ID),
                      database == "nonliquefaciens" ~ paste0("mno_", Cluster.ID),
                      database == "pigrum" ~ paste0("dpi_", Cluster.ID),
                      database == "pneumoniae" ~ paste0("spo_", Cluster.ID),
                      database == "propinquum" ~ paste0("cpr_", Cluster.ID),
                      database == "pseudodiphtheriticum" ~ paste0("cps_", Cluster.ID),
                      database == "rothia" ~ paste0("rsp_", Cluster.ID),
                      TRUE ~ Cluster.ID
                    )) %>%
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
                    select(metagenomic_sample, Species, Cluster.ID, relative_abundance) %>%
                    filter(metagenomic_sample != "CAS_2152") # remove the plate that cannot be linked to age

# rescale abundances from 0 to 100 because we removed the clusters defined exclusively by public genomes
strains_sweeps <- strains_sweeps %>%
                  group_by(metagenomic_sample, Species) %>%
                  mutate(relative_abundance = relative_abundance/sum(relative_abundance)*100) %>%
                  mutate(relative_abundance = round(relative_abundance, 2)) %>%
                  ungroup()
strains_sweeps %>% group_by(metagenomic_sample, Species) %>% summarise(total_abundance = sum(relative_abundance)) %>% pull(total_abundance) %>% table() # floating-point imprecision from can canuse small deviations, that's ok!     

## add columns for Infant Log Number and Sampling Age (months) to strains_sweeps
log_age_subjectid <- read.csv("/2024/12_CAS/4_gtdbtk/log_age_subjectid.csv") %>%
                      rename(Subject.ID = subject.id, Sampling.Age = age_months, Sample.Log.Number = log.number) %>%
                      mutate(Sampling.Age = case_when(Sampling.Age %in% c(1, 2, 3) ~ "2", Sampling.Age %in% c(6, 7) ~ "6", Sampling.Age %in% c(11, 12, 13) ~ "12", TRUE ~ as.character(Sampling.Age))) %>%
                      mutate(Sampling.Age = as.numeric(Sampling.Age)) 
strains_sweeps$Sample.Log.Number <- as.numeric(gsub("CAS_", "", strains_sweeps$metagenomic_sample))
strains_sweeps <- left_join(strains_sweeps, log_age_subjectid, by = "Sample.Log.Number")
    
## create a summary version of strains_sweeps
strains_sweeps_summary <- strains_sweeps %>%
                          group_by(Subject.ID, Sampling.Age, Species, Cluster.ID) %>%
                          reframe(Strains.Detected.Metagenomes = paste0(Cluster.ID, "[", paste(relative_abundance, collapse = ","), "]"), .groups = 'drop') %>%
                          group_by(Subject.ID, Sampling.Age, Species) %>%
                          summarise(Strains.Detected.Metagenomes = paste(Strains.Detected.Metagenomes, collapse = ","), .groups = 'drop') 

# combine the two tables and create Table S6
table_s6 <- full_join(strains_sweeps_summary, isolates_summary, by = c("Species", "Subject.ID", "Sampling.Age"))
table_s6$Strains.Detected.Metagenomes <- ifelse(is.na(table_s6$Strains.Detected.Metagenomes), "", table_s6$Strains.Detected.Metagenomes)
table_s6$Strains.Isolated <- ifelse(is.na(table_s6$Strains.Isolated), "", table_s6$Strains.Isolated)

# not every plate from which an isolate was obtained has a metagenomic sample, we need to add this information to the table
metagenomes_available <- read.csv("/CAS_2025/04_figure1_tableS4_table2/04_tableS4.csv") %>%   
                          select(Subject.ID, Sampling.Age, MPG.Cluster) %>%
                          mutate(key = paste(Subject.ID, Sampling.Age, sep = "_"))

table_s6 <- table_s6 %>% mutate(key = paste(Subject.ID, Sampling.Age, sep = "_"),
         Observation = if_else(key %in% metagenomes_available$key,
                               "",          # leave blank if available
                               "Metagenome from plate sweep not available")) %>% select(-key)  # flag as not available if not in the metagenomes_available table

# save the table
write.csv(table_s6, "06_tableS6.csv", row.names = FALSE)

##############################
# Create Figure 2, Panel A
##############################

# create matrix with strains_sweeps data
# rows = a new_key with Subject.ID+Sampling.Age (e.g.0003THCO_2), columns = Cluster.ID, values = Relative_Abundance scaled to 0-100 per species
strains_sweeps <- strains_sweeps %>% mutate(new_key = paste(Subject.ID, Sampling.Age, sep = "_")) # create new_key

strains_sweeps_matrix <- strains_sweeps %>%
   select(new_key, Cluster.ID, relative_abundance) %>%
   spread(key = Cluster.ID, value = relative_abundance, fill = 0) %>%
   column_to_rownames("new_key") 
dim(strains_sweeps_matrix) # 130x136 --> 12 Subject.ID+Sampling.Age missing (no strains detected in the metagenomes)

# add missing rows using metagenomes_available dataframe
missing <- metagenomes_available %>% filter(!key %in% rownames(strains_sweeps_matrix))
to_bind <- data.frame(matrix(0, nrow = nrow(missing), ncol = ncol(strains_sweeps_matrix)))
rownames(to_bind) <- missing$key
colnames(to_bind) <- colnames(strains_sweeps_matrix)
strains_sweeps_matrix <- rbind(strains_sweeps_matrix, to_bind) %>% as.matrix()
dim(strains_sweeps_matrix) # 142x136 --> all 142 samples are now in the matrix

# sort by age (can't use ComplexHeatmap native clustering because it cannot cluster inside the age groups)
## 2 months 
mo2 <- metagenomes_available %>% filter(Sampling.Age == "2") %>% pull(key) %>% unique()
hc <- hclust(dist(strains_sweeps_matrix[mo2,]))
mo2 <- mo2[hc$order]
## 6 months
mo6 <- metagenomes_available %>% filter(Sampling.Age == "6") %>% pull(key) %>% unique()
hc <- hclust(dist(strains_sweeps_matrix[mo6,]))
mo6 <- mo6[hc$order]
## 12 months
mo12 <- metagenomes_available %>% filter(Sampling.Age == "12") %>% pull(key) %>% unique()
hc <- hclust(dist(strains_sweeps_matrix[mo12,]))
mo12 <- mo12[hc$order]
### order the matrix according to the clustering
strains_sweeps_matrix <- strains_sweeps_matrix[c(mo2, mo6, mo12),]

# create left annotation with age and MPGs
metagenomes_available <- metagenomes_available %>% column_to_rownames("key")
metagenomes_available <- metagenomes_available[row.names(strains_sweeps_matrix),] # order the metagenomes_available according to the strains_sweeps_matrix

left_annotation <- rowAnnotation("MPGs" = as.factor(metagenomes_available$MPG.Cluster), "Infant age" = metagenomes_available$Sampling.Age,
                              col = list("MPGs" = c("1"="#009E73", "2"="#E6E6E6", "3"="#E69F00", "4"="#332288"),
                                   "Infant age" = c("2" = "#CC79A7", "6" = "#56B4E9", "12" = "#0072B2")))

# create top_annotation to show the species
top_annotation_df <- strains_sweeps %>% select(Species, Cluster.ID) %>% distinct() %>% as.data.frame()
rownames(top_annotation_df) <- top_annotation_df$Cluster.ID
top_annotation_df$Species <- gsub("^(\\w)\\w+", "\\1.", top_annotation_df$Species) # shorten the species names
top_annotation_df$Species <- gsub("R. sp902373285", "Rothia sp902373285", top_annotation_df$Species) # Rothia sp902373285 is a special case
top_annotation_df <- top_annotation_df[colnames(strains_sweeps_matrix),]

# draw the heatmap 
htmp <- Heatmap(strains_sweeps_matrix,
        column_split = factor(top_annotation_df$Species, levels = c("S. aureus", "M. catarrhalis", "C. pseudodiphtheriticum", "D. pigrum", "S. mitis","S. pneumoniae", "C. accolens","M. nonliquefaciens", "C. propinquum", "S. agalactiae", "Rothia sp902373285")),
        row_split = as.factor(metagenomes_available$Sampling.Age),
        column_title_rot = 45,
        column_title_gp = grid::gpar(fontsize = 10, fontface = "italic"),
        row_title_gp = grid::gpar(fontsize = 0),
        column_names_gp = grid::gpar(fontsize = 5),
        show_row_names = FALSE,
        show_column_names = TRUE,
        show_parent_dend_line = FALSE,
        cluster_columns = TRUE,
        cluster_column_slices = FALSE,
        cluster_rows = FALSE,
        left_annotation = left_annotation,
        border = TRUE,
        col = circlize::colorRamp2(c(0, 25, 50, 75, 100), c("white", "#9BD9F2", "#0072B2", "#E69F00", "#D73027")),
                  heatmap_legend_param = list(
              title = "Relative abundance/species (%)",  # add a title to legend
                legend_direction = "horizontal"  # set legend direction to horizontal
             ))

svg("06_Figure2_A.svg", width = 10, height = 6)
draw(htmp, heatmap_legend_side="bottom", annotation_legend_side="bottom",
padding = unit(c(2, 2, 2, 20), "mm")) #bottom, left, top, right paddings
dev.off()

# How many strain clusters per species? It will be used to annotate the heatmap later
strains_sweeps %>% select(Species, Cluster.ID) %>% distinct() %>% group_by(Species) %>% summarise(n = n()) %>% arrange(desc(n))
#    Species                                  n
#  1 Moraxella catarrhalis                   25
#  2 Dolosigranulum pigrum                   22
#  3 Corynebacterium pseudodiphtheriticum    20
#  4 Streptococcus mitis                     16
#  5 Staphylococcus aureus                   15
#  6 Streptococcus pneumoniae                12
#  7 Corynebacterium accolens                 9
#  8 Moraxella nonliquefaciens                7
#  9 Corynebacterium propinquum               5
# 10 Streptococcus agalactiae                 3
# 11 Rothia sp902373285                       2

##############################
# Create Figure 2, Panel B (right)
##############################

# calculate the proportion of infants with strains detected
strain_proportion <- strains_sweeps %>%
  select(Sampling.Age, Subject.ID, Species) %>%
  distinct() %>% # some infats have more than one strain detected, do this to avoid overcounting
  group_by(Sampling.Age, Species) %>%
  summarise(number_kids = n()) 

# we need to consider the different number of infants sampled to get the right proportion 
    # 2 months: 50 plates, 50 infants
    # 6 months: 50 plates, 50 infants
    # 12 months: 42 plates, 42 infants
strain_proportion$number_kids_prop <- ifelse(strain_proportion$Sampling.Age == 2, (strain_proportion$number_kids/50)*100, strain_proportion$number_kids)
strain_proportion$number_kids_prop <- ifelse(strain_proportion$Sampling.Age == 6, (strain_proportion$number_kids/50)*100, strain_proportion$number_kids_prop)
strain_proportion$number_kids_prop <- ifelse(strain_proportion$Sampling.Age == 12, (strain_proportion$number_kids/42)*100, strain_proportion$number_kids_prop)

# shorten the species names 
strain_proportion$Species <- gsub("^(\\w)\\w+", "\\1.", strain_proportion$Species) 
strain_proportion$Species <- gsub("R. sp902373285", "Rothia sp902373285", strain_proportion$Species) # Rothia sp902373285 is a special case

# factorize the species column
strain_proportion$Species <- factor(strain_proportion$Species, levels = c("S. aureus", "M. catarrhalis", "C. pseudodiphtheriticum", "D. pigrum", "S. mitis","S. pneumoniae", "C. accolens","M. nonliquefaciens", "C. propinquum", "S. agalactiae", "Rothia sp902373285"))

# list colors for each species
colors <- c(
  "#009E73",  # green
  "#332288",  # dark purple
  "#E69F00",  # orange
  "#56B4E9",  # light sky blue
  "#F0E442",  # yellow
  "#0072B2",  # dark blue
  "#D55E00",  # vermillion
  "#CC79A7",  # reddish purple
  "#999999",  # medium gray
  "#882255",  # dark magenta
  "#AA4499",  # pinkish purple
  "#44AA99"   # teal
)

svg("06_Figure2_B_right.svg",  width = 6, height = 4)
ggplot(strain_proportion, aes(x = Sampling.Age, y = number_kids_prop, group = Species, color = Species, shape = Species)) +
  geom_line(lwd = 1, alpha = 0.8) +
  geom_point(size = 5, alpha = 0.8) +
  labs(x = "Infant age (months)", y = "Proportion of infants with \n detected strains (%)") +
  theme_classic(14) +
  scale_x_continuous(breaks = c(2, 6, 12)) +
  theme(legend.text = element_text(face = "italic")) + # legend text italic
  scale_color_manual(values = colors, name = "") +
  scale_shape_manual(values = c(15, 16, 17, 18, 15, 16, 17, 18, 15, 16, 17, 18), name = "")
dev.off()

strain_proportion %>% as.data.frame()
#    Sampling.Age                 Species number_kids number_kids_prop
# 1             2             C. accolens          13        26.000000
# 2             2           C. propinquum           4         8.000000
# 3             2 C. pseudodiphtheriticum          16        32.000000
# 4             2               D. pigrum           8        16.000000
# 5             2          M. catarrhalis           7        14.000000
# 6             2      M. nonliquefaciens           2         4.000000
# 7             2      Rothia sp902373285          19        38.000000
# 8             2               S. aureus          32        64.000000
# 9             2           S. agalactiae           4         8.000000
# 10            2                S. mitis          15        30.000000
# 11            2           S. pneumoniae           1         2.000000
# 12            6             C. accolens           5        10.000000
# 13            6 C. pseudodiphtheriticum          12        24.000000
# 14            6               D. pigrum          11        22.000000
# 15            6          M. catarrhalis           8        16.000000
# 16            6      M. nonliquefaciens           6        12.000000
# 17            6      Rothia sp902373285           3         6.000000
# 18            6               S. aureus          30        60.000000
# 19            6           S. agalactiae           1         2.000000
# 20            6                S. mitis           7        14.000000
# 21            6           S. pneumoniae           6        12.000000
# 22           12             C. accolens           1         2.380952
# 23           12           C. propinquum          12        28.571429
# 24           12 C. pseudodiphtheriticum          14        33.333333
# 25           12               D. pigrum          13        30.952381
# 26           12          M. catarrhalis          19        45.238095
# 27           12      M. nonliquefaciens           7        16.666667
# 28           12               S. aureus          15        35.714286
# 29           12                S. mitis           1         2.380952
# 30           12           S. pneumoniae           9        21.428571

##############################
# Create Figure 2, Panel B (left)
##############################

# calculate the strain richness per kid / age group
strain_richness <- strains_sweeps %>%
   group_by(Sampling.Age, Species, Subject.ID) %>%
   summarise(strain_n = n())

# add the total number of strains detected for each species per age group
total_strain_n <- strains_sweeps %>%
  select(Sampling.Age, Cluster.ID, Species) %>%
  distinct() %>%
  group_by(Sampling.Age, Species) %>%
  summarise(number_strains_total = n()) 

strain_richness <- left_join(strain_richness, total_strain_n, by = c("Sampling.Age" = "Sampling.Age", "Species" = "Species"))

# shorten the species names
strain_richness$Species <- gsub("^(\\w)\\w+", "\\1.", strain_richness$Species)
strain_richness$Species <- gsub("R. sp902373285", "Rothia sp902373285", strain_richness$Species) # Rothia sp902373285 is a special case

# factorize the species column
strain_richness$Species <- factor(strain_richness$Species, levels = c("S. aureus", "M. catarrhalis", "C. pseudodiphtheriticum", "D. pigrum", "S. mitis","S. pneumoniae", "C. accolens","M. nonliquefaciens", "C. propinquum", "S. agalactiae", "Rothia sp902373285"))

# plot for `number_strains_total`
plot1 <- ggplot(strain_richness, aes(x = Sampling.Age, y = number_strains_total, color = Species)) +
  geom_line(aes(group = Species), lwd = 1.5) +
  geom_point(size = 3) +
  labs(x = "", y = "Total richness") +
  theme_classic(14) +
  scale_x_continuous(breaks = c(2, 6, 12)) +
  theme(legend.position = "none") +
  scale_color_manual(values = scales::alpha(colors, 0.5), name = "Species") +
  facet_wrap(~ Species, nrow = 1) +
  theme(strip.text = element_blank()) # Removes the panel titles

# plot for `strain_n`
plot2 <- ggplot(strain_richness, aes(x = Sampling.Age, y = strain_n, color = Species)) +
  geom_boxplot(outlier.size = 1.5, alpha = 0.1, aes(group = Sampling.Age), fill = NA) +
  geom_jitter(alpha = 0.8, width = 0.15, size = 1) +
  labs(x = "Infant age (months)", y = "Richness per infant") +
  theme_classic(14) +
  scale_x_continuous(breaks = c(2, 6, 12)) +
  theme(legend.position = "none") +
  scale_color_manual(values = scales::alpha(colors, 0.5), name = "Species") +
  facet_wrap(~ Species, nrow = 1) +
  theme(strip.text = element_blank()) # Removes the panel titles

# save the plots
svg("06_Figure2_B_left_top.svg",  width = 9, height = 2)
plot1
dev.off()

svg("06_Figure2_B_left_bottom.svg", width = 9, height = 3)
plot2
dev.off()

# Quick check on number of strains total 
strain_richness %>% 
select(Sampling.Age, Species, number_strains_total) %>%
distinct() %>%
as.data.frame()

#    Sampling.Age                 Species number_strains_total
# 1             2             C. accolens                    9
# 2             2           C. propinquum                    3
# 3             2 C. pseudodiphtheriticum                   14
# 4             2               D. pigrum                    7
# 5             2          M. catarrhalis                    7
# 6             2      M. nonliquefaciens                    2
# 7             2      Rothia sp902373285                    2
# 8             2               S. aureus                   13
# 9             2           S. agalactiae                    3
# 10            2                S. mitis                    9
# 11            2           S. pneumoniae                    1
# 12            6             C. accolens                    5
# 13            6 C. pseudodiphtheriticum                    9
# 14            6               D. pigrum                   10
# 15            6          M. catarrhalis                    9
# 16            6      M. nonliquefaciens                    6
# 17            6      Rothia sp902373285                    2
# 18            6               S. aureus                   13
# 19            6           S. agalactiae                    1
# 20            6                S. mitis                    8
# 21            6           S. pneumoniae                    5
# 22           12             C. accolens                    1
# 23           12           C. propinquum                    4
# 24           12 C. pseudodiphtheriticum                   10
# 25           12               D. pigrum                   10
# 26           12          M. catarrhalis                   19
# 27           12      M. nonliquefaciens                    6
# 28           12               S. aureus                    8
# 29           12                S. mitis                    1
# 30           12           S. pneumoniae                    8
```
