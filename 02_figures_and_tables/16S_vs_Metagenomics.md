# Comparing 16S vs metagenomics data

We are going to use the file `CAS_425days_OTUmat_6genera.xlsx` to compare the 16S and metagenomics data. This file contains the OTU table from Teo 2018, with the same samples and genera.

We are also using the file `phyloseq.rds` which contains taxonomy with metagenomics data.

```R
# Load libraries and data
library(phyloseq)
library(readxl)
library(tidyverse)

physeq <- readRDS("phyloseq.rds")
barcode <- read_excel("CAS_425days_OTUmat_6genera.xlsx")

# Inspect the metagenomics taxonomy 
# Confirm which genus name is used for Alloiococcus/Dolosigranulum
tax_df <- as.data.frame(tax_table(physeq))
unique(tax_df$Genus[grepl("Dolosigranulum|Alloiococcus", tax_df$Genus, ignore.case = TRUE)]) # "Dolosigranulum" in GTDB

# Check 16S percents sum to ~1 across the 6 genera 
barcode %>%
  mutate(rowsum = rowSums(across(ends_with(".percent")))) %>%
  pull(rowsum) %>% summary()

# 16S values are full-community proportions!!
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.006253 0.793210 0.941716 0.844541 0.990947 0.999564 

# Collapse MG to genus, keep ALL genera (don't subset yet) 
physeq_genus <- tax_glom(physeq, taxrank = "Genus", NArm = FALSE)
taxa_names(physeq_genus) <- as.character(tax_table(physeq_genus)[, "Genus"])

mg_mat <- as(otu_table(physeq_genus), "matrix")

head(mg_mat[1:3,1:3])
#               CAS_0006 CAS_0007 CAS_0008
# Abiotrophia          0        0        0
# Acetatifactor        0     4603        0
# Acinetobacter        0        0        0

# Full-community relative abundance, then subset to 6 genera
mg_rel_full <- sweep(mg_mat, 2, colSums(mg_mat), "/")
mg_rel_full[is.nan(mg_rel_full)] <- 0

target_genera <- c("Staphylococcus", "Streptococcus", "Moraxella", "Corynebacterium", "Dolosigranulum", "Haemophilus")
present <- intersect(target_genera, rownames(mg_rel_full))
mg_rel <- mg_rel_full[present, , drop = FALSE]

head(mg_rel[,1:3])
#                     CAS_0006    CAS_0007    CAS_0008
# Staphylococcus  9.997966e-01 0.504523836 0.740503966
# Streptococcus   0.000000e+00 0.000000000 0.001628612
# Moraxella       0.000000e+00 0.000000000 0.000000000
# Corynebacterium 7.167534e-05 0.001884804 0.000000000
# Dolosigranulum  0.000000e+00 0.000000000 0.000000000
# Haemophilus     0.000000e+00 0.000000000 0.000000000

# Long-format both tables and join on Log.No
mg_long <- as.data.frame(t(mg_rel)) %>%
  rownames_to_column("sample_id") %>%
  mutate(Log.No = as.integer(sub("CAS_0*", "", sample_id))) %>%
  pivot_longer(-c(sample_id, Log.No),
               names_to = "Genus", values_to = "MG_relabund")

s16_long <- barcode %>%
  pivot_longer(-Log.No, names_to = "Genus", values_to = "S16_relabund") %>%
  mutate(Genus = sub("\\.percent$", "", Genus),
         Genus = recode(Genus, "Alloiococcus" = "Dolosigranulum"))

paired <- inner_join(mg_long, s16_long, by = c("Log.No", "Genus"))

# Quick checks
length(unique(paired$Log.No)) # 125 samples matched
table(paired$Genus) # all 6 represented? yes!

# Per-genus concordance 
paired %>%
  group_by(Genus) %>%
  summarise(
    n = n(),
    spearman = cor(MG_relabund, S16_relabund, method = "spearman"),
    pearson_clr = {
      # CLR-ish on a pseudocount to handle zeros
      a <- log(MG_relabund + 1e-4); b <- log(S16_relabund + 1e-4)
      cor(a, b, method = "pearson")
    }
  )
#     Genus               n spearman pearson_clr
# 1 Corynebacterium   125   0.440       0.422 
# 2 Dolosigranulum    125   0.474       0.516 
# 3 Haemophilus       125   0.0617      0.0540
# 4 Moraxella         125   0.507       0.649 
# 5 Staphylococcus    125   0.635       0.537 
# 6 Streptococcus     125   0.347       0.379 

# Visualize
ggplot(paired, aes(S16_relabund, MG_relabund)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  facet_wrap(~ Genus, scales = "free") +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "16S relative abundance (Teo 2018)",
       y = "Metagenomics relative abundance (this study)",
       title = "Per-sample concordance, 6 dominant genera") +
  theme_bw()
# scatterplots are heavily zero-dominated! we will try other concordance approaches 

############################################################
# 1. Top-genus agreement (sample-level "winner")
############################################################

top_genus <- paired %>%
  pivot_longer(c(MG_relabund, S16_relabund),
               names_to = "platform", values_to = "relabund") %>%
  group_by(Log.No, platform) %>%
  slice_max(relabund, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(Log.No, platform, top = Genus) %>%
  pivot_wider(names_from = platform, values_from = top)

head(top_genus)
#   Log.No MG_relabund     S16_relabund   
# 1      6 Staphylococcus  Staphylococcus 
# 2      7 Staphylococcus  Staphylococcus 
# 3      8 Staphylococcus  Streptococcus  
# 4     13 Streptococcus   Streptococcus  
# 5     24 Corynebacterium Corynebacterium
# 6     43 Corynebacterium Streptococcus  

# Confusion matrix
tab <- table(MG = top_genus$MG_relabund, S16 = top_genus$S16_relabund)

tab
#                  S16
# MG                Corynebacterium Dolosigranulum Haemophilus Moraxella
#   Corynebacterium              10              6           0         1
#   Dolosigranulum                0              2           0         0
#   Haemophilus                   0              0           0         0
#   Moraxella                     3              2           0        18
#   Staphylococcus               14             10           2         9
#   Streptococcus                 0              0           0         0
#                  S16
# MG                Staphylococcus Streptococcus
#   Corynebacterium              1             3
#   Dolosigranulum               0             0
#   Haemophilus                  1             0
#   Moraxella                    0             0
#   Staphylococcus              26            14
#   Streptococcus                0             3

# Overall agreement
agree <- sum(diag(tab)) / sum(tab)
agree #[1] 0.472

# Cohen's kappa
irr::kappa2(top_genus[, c("MG_relabund", "S16_relabund")])
#  Cohen's Kappa for 2 Raters (Weights: unweighted)

#  Subjects = 125 
#    Raters = 2 
#     Kappa = 0.324 

#         z = 8.07 
#   p-value = 6.66e-16 

# any samples with no signal at all from the 6 genera?
paired %>%
  group_by(Log.No) %>%
  summarise(mg_total = sum(MG_relabund),
            s16_total = sum(S16_relabund)) %>%
  filter(mg_total == 0 | s16_total == 0)
# no!!!

############################################################
# 2. Rank concordance per sample (Kendall's tau or Spearman across the 6 genera, within each sample)
############################################################

rank_conc <- paired %>%
  group_by(Log.No) %>%
  summarise(
    kendall = cor(MG_relabund, S16_relabund, method = "kendall"),
  #  spearman = cor(MG_relabund, S16_relabund, method = "spearman"),
    .groups = "drop"
  )

summary(rank_conc$kendall)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -0.6025  0.2000  0.4472  0.3884  0.6000  0.9661 

# hist(rank_conc$kendall, breaks = 20, xlab = "Kendall's tau (rank agreement within sample, 6 genera)")

############################################################
# 3. Presence/absence concordance at a threshold
############################################################

threshold <- 0.01  # 1% relative abundance threshold for "presence"

paired %>%
  mutate(MG_pres = MG_relabund >= threshold,
         S16_pres = S16_relabund >= threshold) %>%
  group_by(Genus) %>%
  summarise(
    n = n(),
    both_present = sum(MG_pres & S16_pres),
    only_16S     = sum(!MG_pres & S16_pres),
    only_MG      = sum(MG_pres & !S16_pres),
    both_absent  = sum(!MG_pres & !S16_pres),
    jaccard = both_present / (both_present + only_16S + only_MG)
  )

#   Genus               n both_present only_16S only_MG both_absent jaccard
# 1 Corynebacterium   125           40       55       3          27   0.408
# 2 Dolosigranulum    125           22       47       1          55   0.314
# 3 Haemophilus       125            0       22       0         103   0    
# 4 Moraxella         125           29       21       3          72   0.547
# 5 Staphylococcus    125           51        3      43          28   0.526
# 6 Streptococcus     125           13       60       3          49   0.171
```