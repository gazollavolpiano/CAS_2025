# Diversity captured by metagenomic plate-sweeps versus isolates

Create a table with:
* Subject.ID	
* Sample.Log.Number	
* Sampling.Age
* Diversity (number of species) from plate sweep metagenomics
* Species n1 .. nx abundance from plate sweep metagenomics

I will then cross the results with the `TableS1.xlsx`.

```R
# Load necessary libraries
library(phyloseq)
library(tidyverse)

# load phyloseq object
physeq <- readRDS("/12_CAS/11_plate_sweep_taxonomy/phyloseq.rds")
physeq
# otu_table()   OTU Table:         [ 796 taxa and 143 samples ]
# tax_table()   Taxonomy Table:    [ 796 taxa by 7 taxonomic ranks ]

# list the Sample.Log.Number in the physeq
physeq.Sample.Log.Number <- sample_names(physeq) %>% as.data.frame() %>% mutate(log.number = as.numeric(gsub("CAS_", "", .))) %>% pull(log.number)

# read clinical data
log_age_subjectid <- read.csv("/2024/12_CAS/4_gtdbtk/log_age_subjectid.csv") %>%
  rename(Subject.ID = subject.id, Sampling.Age = age_months, Sample.Log.Number = log.number) 
dim(log_age_subjectid) # 149x3
setdiff(physeq.Sample.Log.Number, log_age_subjectid$Sample.Log.Number) # 2152 is missing from log_age_subjectid
read.csv("/labs/workspace/users/camilagv/2024/12_CAS/CAS_425days_OTUmat_clean_culturedLog.csv") %>% filter(Log.No == 2152) # 2152 is missing from pheno too!

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

# calculate the richness (number of species) per sample
richness_per_sample <- microbiome::alpha(physeq, index = "Observed") 
richness_per_sample$Sample <- rownames(richness_per_sample)

# relative abundance of each species
matrix_rela_ab <- microbiome::transform(physeq, transform = "compositional") %>%
  psmelt() %>%
  select(Sample, OTU, Abundance) %>%
  pivot_wider(names_from = OTU, values_from = Abundance, values_fill = 0)
  
# unify all together
diversity_table <- left_join(richness_per_sample, matrix_rela_ab)
dim(diversity_table) # 142 798

# add to log_age_subjectid
diversity_table$Sample.Log.Number <- as.numeric(gsub("CAS_", "", diversity_table$Sample))
diversity_table <- left_join(log_age_subjectid, diversity_table, by = "Sample.Log.Number") 
diversity_table[1:5, 1:10]  # quick check

# save the table and export for further analysis
write.csv(diversity_table, "plate_sweep_diversity_abundance.csv", row.names = FALSE)
```

I renamed the file `plate_sweep_diversity_abundance.`to `TableSC.csv` and downloaded it locally.

I will use both data on local to compare the diversity from plate sweeps versus isolates from `TableS1.xlsx`.

```R
# ============================================================
# Inputs:
#   - TableSC.csv  (wide abundance matrix)
#   - TableS1.xlsx (isolate metadata + GTDB taxonomy)
# Output:
#   - Figure_SA.svg
#   - Table_SC.csv (agreement metrics per sample)
# ============================================================

options(width = 200, scipen = 999) # improve console readability

# Load necessary libraries
library(tidyverse)

# Go to the project directory
setwd("/CAS_V3/New_materials_table_and_Figures")

# ============================================================
# Richness scatter plot: plate-sweep vs isolates
# ============================================================

# Get dataframes with richness 
sweep_alpha <- read_csv("TableSC.csv", show_col_types = FALSE) %>%
                select(Sample, observed) %>%
                rename(rich_sweep = observed)

head(sweep_alpha, 2)
#   Sample   rich_sweep
# 1 CAS_0006          5
# 2 CAS_0007         15

iso_alpha <- readxl::read_xlsx("/CAS_V3/Supplementary_Tables/TableS1.xlsx", skip = 1) %>% 
            filter(`Retained for Analysis` == "Yes") %>%
            mutate(Species = gsub(".*s__(.*)", "\\1", `Taxonomy.GTDB.r214`),
                   Sample = gsub("CAS_(\\d+).*", "CAS_\\1", Isolate)) %>%
            filter(!is.na(Species)) %>%
            filter(Species != "") %>%
            select(Sample, Species) %>%
            distinct() %>%
            group_by(Sample) %>%
            summarise(rich_iso = n_distinct(Species), .groups = "drop")

head(iso_alpha, 2)

alpha <- left_join(iso_alpha, sweep_alpha, by = "Sample") %>% drop_na()
dim(alpha) # 130 x 3

# Check for some infant visits
head(alpha, 2)
#   Sample   rich_iso rich_sweep
# 1 CAS_0006        1          5 --> correct, only S. aureus isolated
# 2 CAS_0007        4         15 --> correct, 4 species isolated

tail(alpha, 2)
#   Sample   rich_iso rich_sweep
# 1 CAS_2105        1          8 --> correct, only S. haemolyticus isolated
# 2 CAS_2412        1         34 --> correct, only S. pneumoniae isolated

# Checked with all Samples and they all match
read_csv("plate_sweep_diversity_abundance.csv", show_col_types = FALSE) %>%
                filter(Sample == "CAS_2412") %>%
                select(-Sample, -observed) %>%
                pivot_longer(cols = everything(), names_to = "OTU", values_to = "Abundance") %>%
                filter(Abundance > 0) %>%
                as.data.frame()

# Spearman correlation with p-value
cor.test(alpha$rich_iso, alpha$rich_sweep, method = "spearman", exact = FALSE)
#         Spearman's rank correlation rho

# data:  alpha$rich_iso and alpha$rich_sweep
# S = 166400, p-value = 0.00000000001924
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.5455346 

svg("Figure_SA.svg", width = 4, height = 4)
ggplot(alpha, aes(x = rich_iso, y = rich_sweep)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  labs(x = "Species richness (pooled across media)", y = "Species richness (plate-sweep metagenome)",
    subtitle = paste0("Spearman \u03C1 = 0.55 (points are infant visits), p < 0.0001")) +
  theme_classic(10)
dev.off()

# What is the median IQR per infant visit for both methods?
summary(alpha$rich_iso)
  #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  #   1.0     1.0     2.0     2.6     3.0     7.0 

summary(alpha$rich_sweep)
  #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  #  2.00    9.00   16.00   21.87   27.00  110.00 
```
