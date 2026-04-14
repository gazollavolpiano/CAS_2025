# Create Figure with Microbiome Profile Groups (MPG) showing the relative abundance of dominant species vs age

Panel A:
  - Upper annotation: hierarchical clustering dendrogram + infant age group bar + Microbiome Profile Groups (MPG) bar
  - x-axis: metagenomic sample
  - y-axis: dominant species 
  - color: log10(1+relative abundance)

Panel B:
  - x-axis: infant age (months)
  - y-axis: CLR tranf. abundance
  - color: Staphylococcus, Moraxella
  - add line

The available data to create the figure can be found in the following sources:
  - `/12_CAS/11_plate_sweep_taxonomy/phyloseq.rds`: phyloseq object with the metagenomic data processed with Kraken2 and Bracken
  - `/12_CAS/CAS_425days_OTUmat_clean_culturedLog.csv`: has several columns with clinical/phenotype for each sample
  - `/12_CAS/4_gtdbtk/log_age_subjectid.csv`: has columns log.number | age_months | subject.id

Table S4 (Microbiome Profile Group (MPG) Clustering Information and Taxa Abundances), the columns will be:
  Subject ID: A unique identifier for each infant in the study
  Sample Log Number: The log number of the sample
  Sampling Age (months): The time point of sampling (2, 6, or 12 months)
  MPG Cluster: The Microbiome Profile Group (MPG) cluster to which the sample is classified
  Other Columns: Relative abundances (0 to 100%) of taxa that were used to cluster the samples into MPGs

Table 2 (Logistic regression results for the association between age and Microbiome Profile Groups), the columns will be:
  MPG: The Microbiome Profile Group cluster
  Age: The age group (6 or 12 months)
  Estimate: The odds ratio
  Confidence Interval: The confidence interval of the odds ratio
  P-value: The p-value of the odds ratio

```R
# library
library(phyloseq)
library(tidyverse)

##############################
# Organize the metagenomic and clinical data
##############################

# load phyloseq object
physeq <- readRDS("/12_CAS/11_plate_sweep_taxonomy/phyloseq.rds")
physeq
# otu_table()   OTU Table:         [ 796 taxa and 143 samples ]
# tax_table()   Taxonomy Table:    [ 796 taxa by 7 taxonomic ranks ]

# list the Sample.Log.Number in the physeq
physeq.Sample.Log.Number <- sample_names(physeq) %>% as.data.frame() %>% mutate(log.number = as.numeric(gsub("CAS_", "", .))) %>% pull(log.number)

# read clinical data
log_age_subjectid <- read.csv("/12_CAS/4_gtdbtk/log_age_subjectid.csv") %>%
  rename(Subject.ID = subject.id, Sampling.Age = age_months, Sample.Log.Number = log.number) 
dim(log_age_subjectid) # 149x3
setdiff(physeq.Sample.Log.Number, log_age_subjectid$Sample.Log.Number) # 2152 is missing from log_age_subjectid
read.csv("/12_CAS/CAS_425days_OTUmat_clean_culturedLog.csv") %>% filter(Log.No == 2152) # 2152 is missing from pheno too!

# remove plate 2152 from physeq
physeq <- prune_samples(sample_names(physeq) %>% setdiff("CAS_2152"), physeq)
physeq  
# otu_table()   OTU Table:         [ 796 taxa and 142 samples ]
# tax_table()   Taxonomy Table:    [ 796 taxa by 7 taxonomic ranks ]

# tidy log_age_subjectid (remove samples not in physeq and recategorize age)
log_age_subjectid <- log_age_subjectid %>% filter(Sample.Log.Number %in% physeq.Sample.Log.Number)
dim(log_age_subjectid) # 142x3

# recategorize Sampling.Age column
log_age_subjectid <- log_age_subjectid %>%
  mutate(Sampling.Age = case_when(
    Sampling.Age %in% c(1, 2, 3) ~ "2",
    Sampling.Age %in% c(6, 7) ~ "6",
    Sampling.Age %in% c(11, 12, 13) ~ "12",
    TRUE ~ as.character(Sampling.Age)
  )) 

# quick check on number of infants
log_age_subjectid %>% 
pull(Subject.ID) %>%
unique() %>%
length() # 58

# answer: how many infants sampled at 2 months, 6 months, and 12 months?
log_age_subjectid %>%
select(Subject.ID, Sampling.Age) %>%
distinct() %>%
group_by(Sampling.Age) %>%
summarise(n = n()) %>%
as.data.frame()
#   Sampling.Age  n
# 1           12 42
# 2            2 50
# 3            6 50

##############################
# Cluster into Microbiome Profile Groups (MPG)
##############################

# STEP 1 - DEFINE COMMON SPECIES:
# Mean relative abundance > 0.1%
# Present in >20% of samples
# Dominating (>50%) at least one sample

# calculate the mean relative abundance of each species
mean_relative_ab <- microbiome::transform(physeq, transform = "compositional") %>%
  psmelt() %>%
  group_by(Species) %>%
  summarise(mean_abundance = mean(Abundance)*100) %>%
  arrange(desc(mean_abundance))

# get candidates for common species
cand_mean_relative_ab <- mean_relative_ab %>%
  filter(mean_abundance > 0.1) %>%
  pull(Species)

# calculate the prevalence among candidates
prevalence <- physeq %>%
  psmelt() %>%
  filter(Abundance != 0) %>%
  filter(Species %in% cand_mean_relative_ab) %>%
  group_by(Species) %>%
  summarise(prevalence = n()) %>%
  arrange(desc(prevalence))

# get candidates for common species
sample_total <- 142 # number of samples

cand_prevalence <- prevalence %>%
  filter(prevalence > (0.2 * sample_total)) %>%
  pull(Species)          

# are cand_prevalence dominant in at least one sample (>50%)?
dominating <-  microbiome::transform(physeq, transform = "compositional") %>%
  psmelt() %>%
  filter(Species %in% cand_prevalence) %>%
  filter(Abundance > 0.5) %>%
  pull(Species) %>%
  unique()             

# STEP 2 - DEFINE "MAJOR" GENERA":
# Remove common species and see what is left
# Concatenate the species were mean relative abundance > 0.1% AND prevalence > 10%
# The other species will be considered "rare"

# calculate the mean relative abundance of each GENERA after REMOVING COMMON SPECIES
mean_relative_ab_genus_select <- microbiome::transform(physeq, transform = "compositional")%>%
  psmelt() %>%
  filter(!Species %in% dominating) %>%
  group_by(Genus) %>%
  mutate(Abundance = Abundance*100) %>%
  summarise(mean_abundance = mean(Abundance)) %>%
  filter(mean_abundance > 0.1) %>%
  pull(Genus)

prev_genus_select <- microbiome::transform(physeq, transform = "compositional")%>%
  psmelt() %>%
  filter(!Species %in% dominating) %>%
  filter(Abundance != 0) %>%
  select(Genus, Sample) %>%
  distinct() %>%
  group_by(Genus) %>%
  summarise(prevalence = n()) %>%
  filter(prevalence > 14.2) %>% # ---> 14.2 = 0.1 * 142 (10% of samples)
  pull(Genus)

genus_to_aggregate <- intersect(mean_relative_ab_genus_select, prev_genus_select)

# STEP 3 - FINALLY AGREGATE THE SPECIES 

# calculate relative abundance
physeq_melt <- microbiome::transform(physeq, transform = "compositional") %>% psmelt() 

# ----  Extract "common species" ----
common_species <- physeq_melt %>%
  filter(Species %in% dominating) %>%
  select(Sample, Species, Abundance) 

# ----  Extract "major genera" ----
agreggated_genus <- physeq_melt %>%
  filter(!Species %in% dominating) %>%
  filter(Genus %in% genus_to_aggregate) %>%
  group_by(Sample, Genus) %>%
  summarise(Abundance = sum(Abundance)) 

# ----  Extract "rare" ----

# agregate the rest 
rare <- physeq_melt %>%
  filter(!Species %in% dominating) %>%
  filter(!Genus %in% genus_to_aggregate) %>%
  group_by(Sample) %>%
  summarise(Abundance = sum(Abundance)) %>%
  mutate(Species = "Rare") %>% 
  select(Sample, Species, Abundance)

# unite the dataframes
colnames(rare) <- c("Sample", "Taxa", "Abundance")
colnames(common_species) <- c("Sample", "Taxa", "Abundance")
colnames(agreggated_genus) <- c("Sample", "Taxa", "Abundance")

physeq_melt_agreggated <- rbind(common_species, agreggated_genus, rare)

physeq_melt_agreggated %>%
  group_by(Sample) %>%
  summarise(abundance_summed = sum(Abundance)) # looks good!

# tranform the relative abundance from 0,1 to 0,100
physeq_melt_agreggated$Abundance <- physeq_melt_agreggated$Abundance * 100

# STEP 4 - CALCULATE THE MPG PROFILES

# reshape data from long to wide format
wide_data <- physeq_melt_agreggated %>%
  spread(Taxa, Abundance, fill = 0) # Replace Taxa and Abundance with your column names

# remove the Sample column for dissimilarity calculation
samples_wide <- wide_data$Sample
wide_data <- wide_data[, -which(names(wide_data) == "Sample")]

# calculate Bray-Curtis dissimilarity
library(vegan)
bc_dist <- vegdist(wide_data, method = "bray")

# perform hierarchical clustering
hc <- hclust(bc_dist, method = "complete")

# determine the optimal number of clusters
silhouette_widths <- sapply(2:10, function(k) {
  cluster_assignments <- cutree(hc, k = k)
  cluster::silhouette(as.vector(cluster_assignments), as.dist(bc_dist))
})

optimal_clusters <- which.max(apply(silhouette_widths, 2, median))

# cut the dendrogram at the determined number of clusters
mpg_clusters <- data.frame("MPG.Cluster"=cutree(hc, optimal_clusters), "Sample"=samples_wide)

# STEP 5 - SAVE THE MPG PROFILES

# agregate all data 
wide_data$Sample <- samples_wide
data_mpg <- left_join(mpg_clusters, wide_data)
data_mpg$Sample.Log.Number <- as.numeric(gsub("CAS_", "", data_mpg$Sample))
data_mpg <- left_join(data_mpg, log_age_subjectid, by = "Sample.Log.Number")

# reorder columns
data_mpg <- data_mpg %>% select(Subject.ID, Sample.Log.Number, Sampling.Age, MPG.Cluster, everything()) 

# save data
write.csv(data_mpg %>% select(-Sample), "04_tableS4.csv", row.names = FALSE)

##############################
# Add M. nonliquefaciens, S. agalactiae, S. pneumoniae, and S. mitis (including any GTDB alphabetic suffix) to the plot for visualization purposes
##############################

# list S. mitis to aggregate
Streptococcus_species <- physeq_melt %>% filter(Genus == "Streptococcus") %>% pull(Species) %>% unique()
Streptococcus_species <- grep("Streptococcus mitis", Streptococcus_species, value = TRUE)
Streptococcus_species # 28 species including GTDB alphabetic suffix

# replace the species name with "Streptococcus mitis"
physeq_melt$Species <-  ifelse(physeq_melt$Species %in% Streptococcus_species, "Streptococcus mitis", physeq_melt$Species)

# calculate the relative abundance of these species
relative_abundances <- physeq_melt %>%
  filter(Species %in% c("Moraxella nonliquefaciens", "Streptococcus agalactiae", "Streptococcus pneumoniae", "Streptococcus mitis")) %>%
  group_by(Sample, Species) %>%
  summarise(Abundance = sum(Abundance)*100) %>% # transform to percentage again
  spread(Species, Abundance, fill = 0) %>%
  as.data.frame()

# add to wide_data
wide_data <- left_join(wide_data, relative_abundances, by = "Sample")

# for Moraxella nonliquefaciens, substract the abundance from the "Moraxella" column
wide_data$Moraxella <- wide_data$Moraxella - wide_data$`Moraxella nonliquefaciens`

# for Streptococcus species, substract the abundance from the "Rare" column
wide_data$Rare <- wide_data$Rare - (wide_data$`Streptococcus pneumoniae` + wide_data$`Streptococcus mitis` + wide_data$`Streptococcus agalactiae`)

# check if the sum of the relative abundances is 100%
wide_data %>% select(-Sample) %>% rowSums() # perfect!

# tidy taxa names (if genus has 1 species use the species name, otherwise add spp.)
# spp. = “species” (plural)
physeq_melt %>% filter(Genus == "Escherichia") %>% pull(Species) %>% unique() # several Escherichia species
colnames(wide_data) <- gsub("^Escherichia$", "*Escherichia* spp.", colnames(wide_data))

physeq_melt %>% filter(Genus == "Moraxella") %>% pull(Species) %>% unique() # several Moraxella species
colnames(wide_data) <- gsub("^Moraxella$", "*Moraxella* spp.", colnames(wide_data))

physeq_melt %>% filter(Genus == "Pseudomonas") %>% pull(Species) %>% unique() # 2 Pseudomonas species
colnames(wide_data) <- gsub("^Pseudomonas$", "*Pseudomonas* spp.", colnames(wide_data))

physeq_melt %>% filter(Genus == "Pseudomonas_E") %>% pull(Species) %>% unique() # several Pseudomonas_E species
colnames(wide_data) <- gsub("^Pseudomonas_E$", "*Pseudomonas*_E spp.", colnames(wide_data))

physeq_melt %>% filter(Genus == "Staphylococcus") %>% pull(Species) %>% unique() # several Staphylococcus species
colnames(wide_data) <- gsub("^Staphylococcus$", "*Staphylococcus* spp.", colnames(wide_data))

# also add italics to the species names
colnames(wide_data) <- gsub("Rothia sp902373285", "*Rothia* sp902373285", colnames(wide_data))
colnames(wide_data) <- gsub("Staphylococcus aureus", "***S. aureus***", colnames(wide_data))
colnames(wide_data) <- gsub("Staphylococcus capitis", "*S. capitis*", colnames(wide_data))
colnames(wide_data) <- gsub("Staphylococcus epidermidis", "*S. epidermidis*", colnames(wide_data))
colnames(wide_data) <- gsub("Staphylococcus haemolyticus", "*S. haemolyticus*", colnames(wide_data))
colnames(wide_data) <- gsub("Moraxella nonliquefaciens", "*M. nonliquefaciens*", colnames(wide_data))
colnames(wide_data) <- gsub("Streptococcus agalactiae", "*S. agalactiae*", colnames(wide_data))
colnames(wide_data) <- gsub("Streptococcus mitis", "*S. mitis*", colnames(wide_data))
colnames(wide_data) <- gsub("Streptococcus pneumoniae", "*S. pneumoniae*", colnames(wide_data))
colnames(wide_data) <- gsub("Corynebacterium accolens", "*C. accolens*", colnames(wide_data))
colnames(wide_data) <- gsub("Corynebacterium propinquum", "*C. propinquum*", colnames(wide_data))
colnames(wide_data) <- gsub("Corynebacterium pseudodiphtheriticum", "***C. pseudodiphtheriticum***", colnames(wide_data))
colnames(wide_data) <- gsub("Dolosigranulum pigrum", "*D. pigrum*", colnames(wide_data))
colnames(wide_data) <- gsub("Moraxella catarrhalis", "***M. catarrhalis***", colnames(wide_data))
colnames(wide_data) <- gsub("Rare", "Rare taxa", colnames(wide_data))

##############################
#  Figure 1, Panel A
##############################
library(ComplexHeatmap)
library(reshape2)

# create a matrix for the heatmap
ab_data_matrix <- wide_data %>%
  select(-Sample) %>%
  t() %>%
  as.matrix()
colnames(ab_data_matrix) <- wide_data$Sample

# order the rows
row.names(ab_data_matrix)
order_rows <- c("*Staphylococcus* spp.", "*S. haemolyticus*", "*S. capitis*", "*S. epidermidis*", "***S. aureus***",
                "Rare taxa","*Rothia* sp902373285", "*Escherichia* spp.", "*Pseudomonas*_E spp.", "*Pseudomonas* spp.", "*S. pneumoniae*", "*S. agalactiae*", "*S. mitis*",
                "*C. accolens*", "***C. pseudodiphtheriticum***", "*C. propinquum*", 
                "*Moraxella* spp.", "*M. nonliquefaciens*", "*D. pigrum*", "***M. catarrhalis***")

ab_data_matrix <- ab_data_matrix[order_rows, ]

# remove Moraxella spp. and Pseudomonas spp. from the matrix (they are barely present and not informative)
to_remove <- c("*Moraxella* spp.", "*Pseudomonas* spp.")
ab_data_matrix <- ab_data_matrix[!rownames(ab_data_matrix) %in% to_remove, ]

# create annotation for age and MPGs
age_data <- data_mpg %>% select(Sample, Sampling.Age) %>% mutate(Sampling.Age = factor(Sampling.Age, levels = c("2", "6", "12")))
mpg_data <- data_mpg %>% select(Sample, MPG.Cluster) %>% mutate(MPG.Cluster = factor(MPG.Cluster, levels = c("1", "2", "3", "4")))

top_annotation <- HeatmapAnnotation("MPGs" = mpg_data$MPG.Cluster, "Infant age group" = age_data$Sampling.Age,
                         col = list("MPGs" = c("1"="#009E73", "2"="#E6E6E6", "3"="#E69F00", "4"="#332288"),
                                   "Infant age group" = c("2" = "#CC79A7", "6" = "#56B4E9", "12" = "#0072B2")))

# use the previous hierarchical clustering object
hc

# draw the heatmap and save it
htmp <- Heatmap(ab_data_matrix,
        column_title = "Individual NP metagenome sample (n=142)",
        column_title_side = "bottom",  #  column title at the bottom
        column_title_gp = grid::gpar(fontsize = 10),
        show_column_names = FALSE,
        border = TRUE,
        cluster_columns = rev(as.dendrogram(hc)),
        row_labels = gt_render(rownames(ab_data_matrix)),
        cluster_rows = FALSE,
        top_annotation = top_annotation,
        row_names_gp = grid::gpar(fontsize = 10),
        col = circlize::colorRamp2(c(0, 25, 50, 75, 100), c("white", "#9BD9F2", "#0072B2", "#E69F00", "#D73027")),
                  heatmap_legend_param = list(
              title = "Relative abundance (%)",  # add a title to legend
                legend_direction = "horizontal"  # set legend direction to horizontal
             ))

svg("04_Figure1_A.svg", width = 6, height = 5)
draw(htmp, 
     heatmap_legend_side = "bottom", 
     annotation_legend_side = "bottom",
     padding = unit(c(2, 2, 2, 20), "mm")) 
dev.off()

##############################
#  Figure 1, Panel B
##############################

# aggregate the data at the genus level
phyloseq_bac_genus <- tax_glom(physeq, taxrank = "Genus")
taxa_names(phyloseq_bac_genus) <- tax_table(phyloseq_bac_genus)[,"Genus"]
phyloseq_bac_genus # 222 taxa and 142 samples

# calculate the CLR abundance of Staphylococcus and Moraxella
selected_genus_abund <- phyloseq_bac_genus %>% 
  microbiome::transform(transform = "clr") %>% # get the CLR transformation
  psmelt() %>%
  filter(Genus %in% c("Staphylococcus", "Moraxella")) %>%
  select(Sample, Genus, Abundance) 

# add the age to the dataframe
selected_genus_abund$Sample.Log.Number <- as.numeric(gsub("CAS_", "", selected_genus_abund$Sample))
selected_genus_abund <- left_join(selected_genus_abund, log_age_subjectid, by = "Sample.Log.Number")

# convert Sampling.Age to numeric
selected_genus_abund$Sampling.Age <- as.numeric(selected_genus_abund$Sampling.Age)

# create the plot
p <- ggplot(selected_genus_abund, aes(x = Sampling.Age, y = Abundance, color = Genus, fill = Genus)) +
  geom_boxplot(aes(group = interaction(Sampling.Age, Genus)), width = 2, size = 0.4, alpha = 0.2) +  # Doubled line size (from 0.2 to 0.4)
  geom_jitter(position = position_dodge(width = 1), shape = 16, size = 4, alpha = 0.1) +
  geom_smooth(method = "lm", se = TRUE, size = 2) +  # Doubled line thickness (from 1 to 2)
  theme_classic(base_size = 15) +  # Makes text 2x bigger (default is ~14-15)
  scale_x_continuous(breaks = c(2, 6, 12)) +
  labs(x = "Infant age group (months)", y = "CLR-transformed abundances") +
  scale_color_manual(values = c("Staphylococcus" = "#009E73", "Moraxella" = "#332288", "Corynebacterium" = "#e69f00"), name = "") +
  scale_fill_manual(values = c("Staphylococcus" = "#009E73", "Moraxella" = "#332288", "Corynebacterium" = "#e69f00"), name = "") +
  theme(
    legend.position = "none",  # Remove legend
    text = element_text(size = 15),  # Doubles text size everywhere
    axis.text = element_text(size = 15),  # Enlarges axis labels
    axis.title = element_text(size = 15, face = "bold"),  # Enlarges axis titles, makes them bold
    strip.text = element_text(size = 15, face = "bold")  # Makes facet panel labels bigger and bold
  ) + 
  facet_wrap(~Genus, ncol = 1)

# Export as needed
svg("04_Figure1_B.svg", width = 4, height = 3)
p
dev.off()

##############################
#  Table 1 (old)
##############################

# create dummy variables for each MPG cluster
data_mpg$staph <- factor(ifelse(data_mpg$MPG.Cluster == 1, 1, 0))
data_mpg$mixed <- factor(ifelse(data_mpg$MPG.Cluster == 2, 1, 0))
data_mpg$cory <- factor(ifelse(data_mpg$MPG.Cluster == 3, 1, 0))
data_mpg$morax <- factor(ifelse(data_mpg$MPG.Cluster == 4, 1, 0))

# ensure Sampling.Age is a factor with "2" as the reference level
data_mpg$Sampling.Age <- factor(data_mpg$Sampling.Age, levels = c("2", "6", "12")) 

# run logistic regression for each cluster
staph <- glm(staph ~ Sampling.Age, family = binomial(), data = data_mpg) %>% broom::tidy(exponentiate = TRUE, conf.int = TRUE)
staph$cluster <- "Staphylococcus-dominated"
staph <- staph[2:3,]

mixed <- glm(mixed ~ Sampling.Age, family = binomial(), data = data_mpg) %>% broom::tidy(exponentiate = TRUE, conf.int = TRUE)
mixed$cluster <- "Mixed"
mixed <- mixed[2:3,]

cory <- glm(cory ~ Sampling.Age, family = binomial(), data = data_mpg) %>% broom::tidy(exponentiate = TRUE, conf.int = TRUE)
cory$cluster <- "Corynebacterium-dominated"
cory <- cory[2:3,]

morax <- glm(morax ~ Sampling.Age, family = binomial(), data = data_mpg) %>% broom::tidy(exponentiate = TRUE, conf.int = TRUE)
morax$cluster <- "Moraxella-dominated"
morax <- morax[2:3,]

# combine results and tidy
or_result <- rbind(staph, mixed, cory, morax) 
or_result$estimate <- round(or_result$estimate, 1)
or_result$term <- gsub("Sampling.Age", "", or_result$term)

or_result %>% select(cluster, term, estimate, conf.low, conf.high, p.value) %>% mutate(conf.low = round(conf.low, 2), conf.high = round(conf.high, 1)) 
#   cluster                   term  estimate conf.low conf.high    p.value
# 1 Staphylococcus-dominated  6          0.3     0.13       0.7 0.00497   
# 2 Staphylococcus-dominated  12         0.1     0.04       0.3 0.00000524

# 7 Moraxella-dominated       6          3.9     0.89      27.2 0.100     
# 8 Moraxella-dominated       12        18       4.69     119.  0.000236  

# 5 Corynebacterium-dominated 6          1.6     0.42       6.5 0.508     
# 6 Corynebacterium-dominated 12         2.7     0.78      10.8 0.127     

# 3 Mixed                     6          2.4     0.89       6.9 0.0909    
# 4 Mixed                     12         1       0.3        3.4 0.969   

##############################
#  Table 1: new (adjusted for subject ID as random effect)
##############################

# ensure subject ID is a factor
data_mpg$Subject.ID <- factor(data_mpg$Subject.ID)

# Mixed-effects logistic regression (glmer)
library(lme4)
staph_glmer <- glmer(staph ~ Sampling.Age + (1 | Subject.ID), data = data_mpg, family = binomial())
mixed_glmer <- glmer(mixed ~ Sampling.Age + (1 | Subject.ID), data = data_mpg, family = binomial())
cory_glmer <- glmer(cory ~ Sampling.Age + (1 | Subject.ID), data = data_mpg, family = binomial())
morax_glmer <- glmer(morax ~ Sampling.Age + (1 | Subject.ID), data = data_mpg, family = binomial())

# helper to tidy glmer output (broom.mixed)
tidy_glmer <- function(model, cluster_name) {
  broom.mixed::tidy(model, effects = "fixed", exponentiate = TRUE, conf.int = TRUE) %>%
    filter(grepl("Sampling.Age", term)) %>%
    mutate(
      cluster = cluster_name,
      term = gsub("Sampling.Age", "", term)
    ) %>%
    select(cluster, term, estimate, conf.low, conf.high, p.value)
}

or_result_glmer <- bind_rows(
  tidy_glmer(staph_glmer,  "Staphylococcus-dominated"),
  tidy_glmer(mixed_glmer,  "Mixed"),
  tidy_glmer(cory_glmer,   "Corynebacterium-dominated"),
  tidy_glmer(morax_glmer,  "Moraxella-dominated")
) %>%
  mutate(
    estimate  = round(estimate,  1),
    conf.low  = round(conf.low,  2),
    conf.high = round(conf.high, 1)
  )

print(or_result_glmer)
#   cluster                   term      estimate  conf.low conf.high    p.value
# 1 Staphylococcus-dominated  6              0.3      0.11  7   e- 1 0.00573   
# 2 Staphylococcus-dominated  12             0.1      0.03  3   e- 1 0.0000270 
# 3 Mixed                     6              2.5      0.87  7.1 e+ 0 0.0899    
# 4 Mixed                     12             1        0.3   3.4 e+ 0 0.985     
# 5 Corynebacterium-dominated 6              1.5      0.32  6.9 e+ 0 0.620     
# 6 Corynebacterium-dominated 12             3        0.66  1.35e+ 1 0.156     
# 7 Moraxella-dominated       6          57131.      41.5   7.86e+ 7 0.00297   
# 8 Moraxella-dominated       12    1383721934.  134074.    1.43e+13 0.00000805

##############################
#  Table 1: additional GEE (better than previous glmer)
##############################

# sort first (required for geeglm)
data_mpg <- data_mpg %>% arrange(Subject.ID, Sampling.Age)

tidy_gee <- function(model, cluster_name) {
  broom::tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(grepl("Sampling.Age", term)) %>%
    mutate(
      cluster = cluster_name,
      term    = gsub("Sampling.Age", "", term)
    ) %>%
    select(cluster, term, estimate, conf.low, conf.high, p.value)
}

# Convert MPG dummy variables from factor to numeric
data_mpg$staph <- as.numeric(as.character(data_mpg$staph))
data_mpg$mixed <- as.numeric(as.character(data_mpg$mixed))
data_mpg$cory  <- as.numeric(as.character(data_mpg$cory))
data_mpg$morax <- as.numeric(as.character(data_mpg$morax))

library(geepack)

staph_gee <- geeglm(staph ~ Sampling.Age, id = Subject.ID, data = data_mpg, family = binomial(), corstr = "exchangeable")
mixed_gee <- geeglm(mixed ~ Sampling.Age, id = Subject.ID, data = data_mpg, family = binomial(), corstr = "exchangeable")
cory_gee  <- geeglm(cory  ~ Sampling.Age, id = Subject.ID, data = data_mpg, family = binomial(), corstr = "exchangeable")
morax_gee <- geeglm(morax ~ Sampling.Age, id = Subject.ID, data = data_mpg, family = binomial(), corstr = "exchangeable")

or_result_gee <- bind_rows(
  tidy_gee(staph_gee, "Staphylococcus-dominated"),
  tidy_gee(mixed_gee, "Mixed"),
  tidy_gee(cory_gee,  "Corynebacterium-dominated"),
  tidy_gee(morax_gee, "Moraxella-dominated")
) %>%
  mutate(
    estimate  = round(estimate,  1),
    conf.low  = round(conf.low,  2),
    conf.high = round(conf.high, 1)
  )

print(or_result_gee)

#   cluster                   term  estimate conf.low conf.high    p.value
# 1 Staphylococcus-dominated  6          0.3     0.14       0.6 0.00214   
# 2 Staphylococcus-dominated  12         0.1     0.04       0.3 0.00000494
# 3 Mixed                     6          2.4     0.87       6.4 0.0914    
# 4 Mixed                     12         1       0.3        3.5 0.979     
# 5 Corynebacterium-dominated 6          1.4     0.45       4.4 0.552     
# 6 Corynebacterium-dominated 12         2.4     0.72       7.8 0.157     
# 7 Moraxella-dominated       6          4.5     0.97      21   0.0543    
# 8 Moraxella-dominated       12        21.1     4.34     102.  0.000155  
```
