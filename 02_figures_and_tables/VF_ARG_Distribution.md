# Create Figures with VF and ARG distribution (Adjusted for Multiple Cultivation Media + CAS Isolate & Metagenomic Detection)

The available data to create the table can be found in the following sources:
  - `/CAS_2025/08_prevalence_strains_colonized_infants_tables7/08_tableS7.csv`: summary data of strain clusters detected 
  - `/CAS_2025/03_tableS3/03_tableS3.csv`: VF and ARG detections
  - `/CAS_2025/06_figure2_tableS5_tableS6/06_tableS5.csv`: classification of strains into strain clusters
  - `/CAS_2025/01_tableS1/01_tableS1.csv`: has the metadata for the isolates including the age of the infant and genome statistics
  - `/2025/01_CAS/strain_scan/profile_plate_sweep/strain_scan_results_plate_sweep_11_databases.csv`: the strain clusters detected in the metagenomes of the infants

```R
# library
library(tidyverse)

# increase the stdout limit
options(width = 200)

##############################
# Prepare the data
##############################

# select what strains have CAS Isolate & Metagenomic Detection
selected_clusters <- read.csv("/CAS_2025/08_prevalence_strains_colonized_infants_tables7/08_tableS7.csv") %>%
  filter(Existence.Evidence == "CAS Isolate & Metagenomic Detection") %>%
  select(Species, Cluster.ID)

# add a prefix to the Cluster ID to indicate the species
selected_clusters <- selected_clusters %>%
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

# list the genomes of the selected strains
genomes_selected_clusters <- read.csv("/CAS_2025/06_figure2_tableS5_tableS6/06_tableS5.csv") %>%
                                        separate_rows(Genomes, sep = ",") %>%
                                        mutate(Genomes = gsub("filt_", "", Genomes)) %>%
                                        filter(grepl("CAS", Genomes)) %>% # filter for CAS genomes
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
                                        )) %>%
                                        filter(Cluster.ID %in% selected_clusters$Cluster.ID) %>%
                                        select(Species, Cluster.ID, Genomes)
setdiff(selected_clusters$Cluster.ID, genomes_selected_clusters$Cluster.ID) %>% length() # 0

# add a column with "Group_Code" to the genomes_selected_clusters
# this column will be used to select a representative genome for each strain cluster that was detected in the same infant but in different media (thus adjusteding for multiple cultivation media)
## first, add more information to the genomes_selected_clusters
genomes_selected_clusters <- left_join(genomes_selected_clusters, read.csv("/CAS_2025/01_tableS1/01_tableS1.csv") %>% select(Isolate, Subject.ID, Sampling.Age, Quality.Score, Contigs.Number, N50), by = c("Genomes" = "Isolate"))

## add the Group_to_filter
genomes_selected_clusters <- genomes_selected_clusters %>%
  group_by(Species, Cluster.ID, Subject.ID, Sampling.Age) %>%
  mutate(Group_Code = paste0(Cluster.ID, "_", Subject.ID, "_", Sampling.Age)) %>%
  ungroup()

genomes_selected_clusters$Group_Code %>% unique() %>% length() # 248 genomes should be selected

## add a column Selected_Genome to select a representative for each Group_Code according to the highest Quality.Score, the lowest Contigs.Number, and the highest N50
representatives <- genomes_selected_clusters %>%
  group_by(Species, Cluster.ID, Group_Code) %>%
  arrange(desc(Quality.Score), Contigs.Number, desc(N50)) %>%
  slice(1) %>%          # pick the top row per group
  ungroup()

genomes_selected_clusters <- genomes_selected_clusters %>%
  left_join(
    representatives %>% 
      select(Species, Cluster.ID, Group_Code, Genomes) %>%
      rename(Selected_Genome = Genomes),
    by = c("Species", "Cluster.ID", "Group_Code")
  )

genomes_selected_clusters$Selected_Genome %>% unique() %>% length() # 248 genomes selected

# save the data
write.csv(genomes_selected_clusters, "genomes_selected_for_each_cluster.csv", row.names = FALSE)

# keep only selected genomes on the dataframe
isolates_selected <- genomes_selected_clusters %>% pull(Selected_Genome) %>% unique()
genomes_selected_clusters <- genomes_selected_clusters %>% filter(Genomes %in% isolates_selected)
dim(genomes_selected_clusters) # 248x10

# read the VF and ARG detections and filter for the selected genomes
vf_args <- read.csv("/CAS_2025/03_tableS3/03_tableS3.csv") %>% 
            filter(Genome %in% isolates_selected)

# add the Cluster.ID and Sampling.Age to the vf_args
vf_args <- left_join(vf_args, genomes_selected_clusters %>% select(Genomes, Cluster.ID, Sampling.Age), by = c("Genome" = "Genomes"))

# read data from plate sweeps metagenomes processed with strainscan and tidy it
strains_sweeps <- read.csv("/2025/01_CAS/strain_scan/profile_plate_sweep/strain_scan_results_plate_sweep_11_databases.csv", header = TRUE, stringsAsFactors = FALSE) %>%
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

############################################################
#                   Heatmap with VFs and ARGs
############################################################

# select taxa of interest
# taxa <- "Corynebacterium accolens"
# taxa <- "Streptococcus agalactiae"
# taxa <- "Staphylococcus aureus"
# taxa <- "Moraxella catarrhalis"
# taxa <- "Streptococcus mitis"
# taxa <- "Moraxella nonliquefaciens"
# taxa <- "Dolosigranulum pigrum"
# taxa <- "Streptococcus pneumoniae" 
# taxa <- "Corynebacterium propinquum"
# taxa <- "Corynebacterium pseudodiphtheriticum"
# taxa <- "Rothia sp902373285"

# genome information about the selected strains for the taxa of interest
genomes_selected_clusters_filt_taxa <- genomes_selected_clusters %>% filter(Species == taxa) 

strain_number <- genomes_selected_clusters_filt_taxa %>% pull(Cluster.ID) %>% unique() %>% length()
genome_number <- genomes_selected_clusters_filt_taxa %>% nrow()

# select the VFs and ARGs for the taxa of interest
vf_args_interest <- vf_args %>% filter(Genome %in% unique(genomes_selected_clusters_filt_taxa$Genomes))

# --------------------------------------------------------------
# get basic statistic before we create the heatmap
strain_number # X strains isolated and detected in sweeps
genome_number # x genomes for the X strains
vf_args_interest %>% filter(Database == "ResFinder") %>% pull(Genome) %>% unique() %>% length() # x genomes with ARGs detected
vf_args_interest %>% filter(Database == "VFDB") %>% pull(Genome) %>% unique() %>% length() # x genomes with VFs detected
vf_args_interest %>% filter(Database != "") # overview 
# --------------------------------------------------------------

# create a dataframe with Cluster.ID and the number of genomes
genome_n_clust <- genomes_selected_clusters_filt_taxa %>%
  group_by(Cluster.ID) %>%
  summarise(number_of_genomes = n())

# calculate the prevalence of VFs and ARGs in each strain cluster
vf_preval <- vf_args_interest %>%
  filter(Database == "VFDB") %>%
  select(Genome, Cluster.ID, Observation) %>%
  distinct() %>%
  group_by(Cluster.ID, Observation) %>%
  summarise(n=n()) %>%
  left_join(genome_n_clust, by = "Cluster.ID") %>%
  mutate(perc = round((n/number_of_genomes)*100,0)) 

arg_preval <- vf_args_interest %>%
  filter(Database == "ResFinder") %>%
  select(Genome, Cluster.ID, Gene) %>%
  distinct() %>%
  group_by(Cluster.ID, Gene) %>%
  summarise(n=n()) %>%
  left_join(genome_n_clust, by = "Cluster.ID") %>%
  mutate(perc = round((n/number_of_genomes)*100,0))

# transform into a matrix
vf_matrix <- genome_n_clust %>%
                    left_join(vf_preval %>% 
                            select(-n) %>%
                            pivot_wider(names_from = Observation, values_from = perc, values_fill = list(perc = 0)), by = c("Cluster.ID", "number_of_genomes")) %>% 
                            ungroup() %>% 
                            mutate(ID = paste0(Cluster.ID, " (n=", number_of_genomes, ")")) %>%
                            select(-number_of_genomes, -Cluster.ID) %>%
                            column_to_rownames("ID") %>%
                            mutate_all(~replace_na(.x, 0)) %>%
                            as.matrix()

arg_matrix <- genome_n_clust %>%
                    left_join(arg_preval %>% 
                            select(-n) %>%
                            pivot_wider(names_from = Gene, values_from = perc, values_fill = list(perc = 0)), by = c("Cluster.ID", "number_of_genomes")) %>% 
                            ungroup() %>% 
                            mutate(ID = paste0(Cluster.ID, " (n=", number_of_genomes, ")")) %>%
                            select(-number_of_genomes, -Cluster.ID) %>%
                            column_to_rownames("ID") %>%
                            mutate_all(~replace_na(.x, 0)) %>%
                            as.matrix()


  # look at core VFs and ARGs (present in 100% of the genomes of a strain cluster)
  vf_matrix %>% apply(2, function(x) sum(x == 100)) %>% as.data.frame() %>% rownames_to_column("VF") %>% arrange(desc(.)) 


# gather data to also create a heatmap with age distribution of the strains in the metagenomic data
## create a matrix with the prevalence of strains in the metagenomic data at each age
age_n_clust <- strains_sweeps %>%
    filter(Species == taxa) %>%
    filter(Cluster.ID %in% genomes_selected_clusters_filt_taxa$Cluster.ID) %>%
    group_by(Cluster.ID, Sampling.Age) %>%
    summarise(n=n()) %>%
    group_by(Cluster.ID) %>%
    mutate(number_of_detections = sum(n)) %>%
    mutate(perc = round((n/number_of_detections)*100,0)) %>%
    ungroup() %>%
    select(Cluster.ID, Sampling.Age, perc) 
## add an ID column to the age_n_clust similar to the vf_matrix and arg_matrix
age_n_clust <- left_join(age_n_clust, genome_n_clust %>% select(Cluster.ID, number_of_genomes), by = "Cluster.ID") %>%
                mutate(ID = paste0(Cluster.ID, " (n=", number_of_genomes, ")")) %>%
                select(-number_of_genomes, -Cluster.ID) %>%
                pivot_wider(names_from = Sampling.Age, values_from = perc, values_fill = 0) %>%
                column_to_rownames("ID") %>%
                as.matrix()
## ordering by age 2, then age 6, then age 12
strains_ord_by_age <- age_n_clust %>%
    as.data.frame() %>%
    mutate(most_prev_age = apply(age_n_clust, 1, which.max)) %>%
    arrange(most_prev_age) 

strains_final_order <- c()
for(i in unique(strains_ord_by_age$most_prev_age)){
  if(i==1)(
    month <- "2"
  ) else if(i==2){
    month <- "6"
  } else {
    month <- "12"
  }

  to_order <- strains_ord_by_age[strains_ord_by_age$most_prev_age == i, ] 
  to_order <- to_order %>% arrange(to_order[[month]]) # trying ascending order
  strains_final_order <- c(strains_final_order, row.names(to_order))
}

## order the matrices
vf_matrix <- vf_matrix[strains_final_order, ]
arg_matrix <- arg_matrix[strains_final_order, ]
age_matrix <- age_n_clust[strains_final_order, ]
# age_matrix <- age_matrix[, c("2", "6")]  # sometimes you might need to select only some ages
age_matrix <- age_matrix[, c("2", "6", "12")]

## simplify the row/col names 
### age_matrix
rownames(age_matrix) <- gsub(".*_(C.*) \\(.*", "\\1", rownames(age_matrix))
### vf_matrix
rownames(vf_matrix) <- gsub(".*_(C.*)", "\\1", rownames(vf_matrix))
colnames(vf_matrix) <- gsub("VF.*- (.*)", "\\1", colnames(vf_matrix))
colnames(vf_matrix) <- gsub("<alpha>", "Alpha", colnames(vf_matrix))
colnames(vf_matrix) <- gsub("<beta>", "Beta", colnames(vf_matrix))
colnames(vf_matrix) <- gsub("<gamma>", "Gamma", colnames(vf_matrix))
colnames(vf_matrix) <- gsub("<delta>", "Delta", colnames(vf_matrix))
colnames(vf_matrix) <- gsub("Effector delivery system", "EDS", colnames(vf_matrix))
colnames(vf_matrix) <- gsub("proteins", "PTNs", colnames(vf_matrix))
colnames(vf_matrix) <- gsub("protein", "PTN", colnames(vf_matrix))

# plot the heatmaps
library(ComplexHeatmap)

htmp_vf <- Heatmap(vf_matrix,
        show_column_names = TRUE,
        border = TRUE,
        row_names_side = "left",
        cluster_columns = TRUE,
        cluster_rows = FALSE,
        column_names_rot = 45,  # column names at a 45-degree angle
        row_names_gp = grid::gpar(fontsize = 10),
        column_names_gp = grid::gpar(fontsize = 8),
        col = circlize::colorRamp2(c(0, 25, 50, 75, 100), c("white", "#9BD9F2", "#0072B2", "#E69F00", "#D73027")),
        heatmap_legend_param = list(title = "Prev strains (%)", legend_direction = "horizontal"), # legend direction to horizontal
        show_heatmap_legend = TRUE) 

htmp_arg <- Heatmap(arg_matrix,
        show_column_names = TRUE,
        show_row_names = FALSE, # do not show row names
        border = TRUE,
        row_names_side = "left",
        cluster_columns = TRUE,
        cluster_rows = FALSE,
        column_names_rot = 45,  # column names at a 45-degree angle
        row_names_gp = grid::gpar(fontsize = 10),
        column_names_gp = grid::gpar(fontsize = 8),
        col = circlize::colorRamp2(c(0, 25, 50, 75, 100), c("white", "#9BD9F2", "#0072B2", "#E69F00", "#D73027")),
        heatmap_legend_param = list(title = "Prev strains (%)", legend_direction = "horizontal"), # legend direction to horizontal
        show_heatmap_legend = FALSE)  

htmp_age <- Heatmap(age_matrix,
        show_column_names = TRUE,
        border = TRUE,
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        column_names_rot = 0,
        row_names_gp = grid::gpar(fontsize = 10),
         width = unit(25, "mm"), 
        column_names_gp = grid::gpar(fontsize = 10, just = "center"),
        col = circlize::colorRamp2(c(0, 25, 50, 75, 100), c("white", "#F0E442", "#E69F00", "#CC99FF","#332288")),
        heatmap_legend_param = list(title = "Prop metag. (%)", legend_direction = "horizontal"))  # legend direction to horizontal

svg(paste0("VF_ARG_age_heatmap_",gsub(" ", "_", taxa),".svg"), width = 10, height = 6)
draw(
  htmp_vf + htmp_arg +  htmp_age,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  column_title = bquote(atop(italic(.(taxa)),
                             .(genome_number) ~ " genomes from " ~ .(strain_number) ~ " strains detected by metagenomics and with isolates available")),
  column_title_gp = gpar(fontsize = 10),  
  padding = unit(c(2, 30, 2, 20), "mm") # bottom, left, top, right
)
dev.off()

```
