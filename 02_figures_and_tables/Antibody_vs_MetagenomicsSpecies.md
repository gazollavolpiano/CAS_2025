# Association between detectable antibodies and MPGs 

Teo et al. 2015 used IgG antibody data (from https://doi.org/10.1136/thoraxjnl-2011-200650) to validate that the 16S rRNA sequencing was detecting clinically meaningful bacterial species (since it alone cannot resolve species level). They measured IgG antibodies against species-specific surface proteins (H. influenzae and S. pneumoniae) at 12 months of age and correlated those antibody levels with sequencing data. Results are shown in Figure S2 in the paper.

We do not detect H. influenzae reliably in our data, but we do detect S. pneumoniae. We can attempt to replicate the analysis for S. pneumoniae using antibody data and our sequencing data.

The variables used for the association between detectable IgG1 antibodies to S. pneumoniae pneumococcal surface protein A1, A2, or C (measured at 12 months of age) were:
  522.	pspa1g1_1	PspA1_IgG1 ng/ml  (Cut off 2000: <2000 assigned 1000ng/ml)
  IgG1 to PspA1 at 1 year. PspA1 is pneumococcal surface protein A family 1 of S. pneumoniae Continuous ng/ml

  527.	pspa2g1_1	PspA2_IgG1 ng/ml  (Cut off 2000: <2000 assigned 1000ng/ml)
  IgG1 to PspA2 at 1 year. PspA2 is pneumococcal surface protein A family 2 of S. pneumoniae  Continuous ng/ml

  532.	pspcg1_1	PspC_IgG1 ng/ml  (Cut off 2000: <2000 assigned 1000ng/ml)
  IgG1 to PspC at 1 year. PspC is pneumococcal surface protein C of S. pneumoniae  Continuous ng/ml

What we want to answer in a exposure-response model is: 
- Cumulative carriage during 0–12 months (any positive visit at 2, 6, or 12 months) - antibody at 1 year
- Carriage at 12 months only (the visit closest in time to the antibody measurement) - antibody at 1 year

```R
# libraries
library(tidyverse)
library(phyloseq)

# load antibody data
antibody_data <- read.csv("cas_immunology_db_131016.csv") %>%
                select(SUBJ_ID, pspa1g1_1, pspa2g1_1, pspcg1_1)

head(antibody_data, 3)
#    SUBJ_ID pspa1g1_1 pspa2g1_1 pspcg1_1
# 1 ........      1000    1000.0     1000
# 2 ........      1000  817378.5     1000
# 3 ........      1000    1000.0     1000

# prepare plate sweeps data 
# - ever_pos_Streptococcus: Streptococcus genus detected in any of 2/6/12 month 
# - pos_12mo_Streptococcus: Streptococcus genus detected in the 12-month  only
# - ever_pos_pneumoniae: Streptococcus pneumoniae species detected in any of 2/6/12 month 
# - pos_12mo_pneumoniae: Streptococcus pneumoniae species detected in the 12-month  only

log_age_subjectid <- read.csv("/12_CAS/4_gtdbtk/log_age_subjectid.csv") %>%
  rename(Subject.ID = subject.id, Sampling.Age = age_months, Sample.Log.Number = log.number) %>%
  mutate(Sampling.Age = case_when(
    Sampling.Age %in% c(1, 2, 3) ~ "2",
    Sampling.Age %in% c(6, 7) ~ "6",
    Sampling.Age %in% c(11, 12, 13) ~ "12",
    TRUE ~ as.character(Sampling.Age)
  )) 

physeq <- readRDS("/12_CAS/11_plate_sweep_taxonomy/phyloseq.rds")
physeq_melt_relab <- psmelt(physeq)
physeq_melt_relab$Sample.Log.Number <- as.numeric(gsub("CAS_", "", physeq_melt_relab$Sample))

physeq_melt_relab <- physeq_melt_relab %>% left_join(log_age_subjectid)

mg_data <- physeq_melt_relab %>%
  group_by(Subject.ID) %>% 
  summarise(
    ever_pos_Streptococcus = any(Genus == "Streptococcus" & Abundance != 0),
    pos_12mo_Streptococcus = any(Genus == "Streptococcus" & Abundance != 0 & Sampling.Age == 12),
    ever_pos_pneumoniae    = any(Species == "Streptococcus pneumoniae" & Abundance != 0),
    pos_12mo_pneumoniae    = any(Species == "Streptococcus pneumoniae" & Abundance != 0 & Sampling.Age == 12),
    has_12mo               = any(Sampling.Age == 12),
    .groups = "drop"
  ) %>%
  select(Subject.ID, ever_pos_Streptococcus, pos_12mo_Streptococcus, ever_pos_pneumoniae, pos_12mo_pneumoniae, has_12mo)

head(mg_data,2) %>% as.data.frame()
#   Subject.ID ever_pos_Streptococcus pos_12mo_Streptococcus ever_pos_pneumoniae
# 1   ........                   TRUE                   TRUE               FALSE
# 2   ........                   TRUE                   TRUE               FALSE
#   pos_12mo_pneumoniae has_12mo
# 1               FALSE     TRUE
# 2               FALSE     TRUE

# unite data
antibody_data <- antibody_data %>% rename(Subject.ID = SUBJ_ID)
combined <- antibody_data %>% left_join(mg_data, by = "Subject.ID") %>% drop_na(ever_pos_Streptococcus)
dim(combined) # 58x9

combined <- combined %>%
  filter(!is.na(pspa1g1_1)) %>%   # drops the 1 infant with NA antibody data
  mutate(seropos_pneumo = (pspa1g1_1 > 2000) | (pspa2g1_1 > 2000) | (pspcg1_1 > 2000))

# How many seropositive?
combined %>% filter(seropos_pneumo) %>% nrow() # 11 seropositive

# How many ever positive for S. pneumoniae?
combined %>% filter(ever_pos_pneumoniae) %>% nrow() # 15 ever positive for S. pneumoniae

# save combined data for future use (csv)
write.csv(combined, "combined_antibody_metagenomics.csv", row.names = FALSE)

# Helper
run_fisher <- function(exp_var, data) {
  tab <- table(exposure = data[[exp_var]], seropos = data$seropos_pneumo)
  cat("\n---", exp_var, "---\n")
  print(tab)
  ft <- fisher.test(tab)
  tibble(
    exposure = exp_var,
    n = nrow(data),
    n_carriers = sum(data[[exp_var]]),
    seropos_in_carriers    = sum(data[[exp_var]] & data$seropos_pneumo),
    seropos_in_noncarriers = sum(!data[[exp_var]] & data$seropos_pneumo),
    OR = unname(ft$estimate),
    OR_lower = ft$conf.int[1],
    OR_upper = ft$conf.int[2],
    p = ft$p.value
  )
}

# Cumulative (any visit in year 1)
res_cumulative <- run_fisher("ever_pos_pneumoniae", combined)
# --- ever_pos_pneumoniae ---
#         seropos
# exposure FALSE TRUE
#    FALSE    38    4
#    TRUE      8    7

# Contemporaneous (12-month visit only) — restrict to infants with 12mo sample
combined_12mo <- combined %>% filter(has_12mo)   # use your existing flag
res_12mo <- run_fisher("pos_12mo_pneumoniae", combined_12mo)
# --- pos_12mo_pneumoniae ---
#         seropos
# exposure FALSE TRUE
#    FALSE    28    4
#    TRUE      5    4

# Also test the broader Streptococcus genus as exploratory
res_genus_cumulative <- run_fisher("ever_pos_Streptococcus", combined)
#         seropos
# exposure FALSE TRUE
#    FALSE    15    2
#    TRUE     31    9

res_genus_12mo <- run_fisher("pos_12mo_Streptococcus", combined_12mo)
# --- pos_12mo_Streptococcus ---
#         seropos
# exposure FALSE TRUE
#    FALSE    20    2
#    TRUE     13    6

bind_rows(res_cumulative, res_12mo, res_genus_cumulative, res_genus_12mo) %>% mutate(p_BH = p.adjust(p, method = "BH")) %>% as.data.frame()
#                 exposure  n n_carriers seropos_in_carriers
# 1    ever_pos_pneumoniae 57         15                   7
# 2    pos_12mo_pneumoniae 41          9                   4
# 3 ever_pos_Streptococcus 57         40                   9
# 4 pos_12mo_Streptococcus 41         19                   6
#   seropos_in_noncarriers       OR  OR_lower OR_upper          p       p_BH
# 1                      4 7.902890 1.5831208 46.57393 0.00432813 0.01731252
# 2                      4 5.295199 0.7352930 40.46444 0.05440963 0.10881927
# 3                      2 2.150893 0.3743470 22.90310 0.47608670 0.47608670
# 4                      2 4.444323 0.6628201 51.57873 0.11524771 0.15366361
```
