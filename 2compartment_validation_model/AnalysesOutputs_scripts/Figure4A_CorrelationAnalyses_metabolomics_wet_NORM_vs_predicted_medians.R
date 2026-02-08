###############################################################################
#
# Compare in silico predicted metabolite fold changes (from both
# the 2-compartment and 6-compartment GIT models) against
# experimentally measured metabolite fold changes (cecal metabolomics,
# PQN-normalized data).
#
# OVERVIEW:
#   1. Load and normalize metabolomics data (PQN-normalized, filtered)
#   2. Load 48h predictions (2-compartment ileum-cecum model) and compute
#      fold changes relative to control (G1) medians
#   3. Load 96h predictions (6-compartment GIT model) fold changes
#   4. For each model, compute per-metabolite Spearman correlations between
#      predicted and experimental fold changes across treatment groups
#   5. Generate Figure 4A: faceted scatter plots for 7 key metabolites
#
# TREATMENT GROUPS AND SUPPLEMENTATION DOSES:
#   The experimental design tested three dietary supplements at two doses.
#   Each simulation used a specific compound_amount that maps to a treatment:
#     G1 = Control (corn-based grower diet, no supplement)
#     G2 = Cellulose low  (1 g/100g)      -> compound_amount = "Cellulose_1"
#     G3 = Cellulose high (4 g/100g)      -> compound_amount = "Cellulose_4"
#     G4 = L-Threonine low  (0.4 g/100g)  -> compound_amount = "L-Threonine_0.4"
#     G5 = L-Threonine high (2.5 g/100g)  -> compound_amount = "L-Threonine_2.5"
#     G6 = Starch low  (1.25 g/100g)      -> compound_amount = "Starch_n19_1.25"
#     G7 = Starch high (2 g/100g)         -> compound_amount = "Starch_n19_2"

###############################################################################

library(tidyverse)
library(readxl)
library(ggpubr)
library(tidyr)
library(dplyr)
library(rstatix)
library(stringr)
library(broom)
library(purrr)

# ---- Set base directory ----
base_dir <- "path/to/metabolomics_data/"
setwd(base_dir)

# ---- Paths to simulation prediction tables ----
path_2comp_predictions <- "path/to/2comp_simulation_outputs/"


###############################################################################
# DEFINE TREATMENT-TO-DOSE MAPPING (single chosen combination)
###############################################################################
# This mapping defines which simulated supplementation dose corresponds to
# each experimental treatment group. The doses were chosen to best match
# the actual experimental supplementation levels.

treatment_dose_map <- tibble(
  Treatment      = paste0("G", 2:7),
  compound_amount = c("Cellulose_1",       # G2: cellulose low dose
                      "Cellulose_4",       # G3: cellulose high dose
                      "L-Threonine_0.4",   # G4: L-threonine low dose
                      "L-Threonine_2.5",   # G5: L-threonine high dose
                      "Starch_n19_1.25",   # G6: starch low dose
                      "Starch_n19_2")      # G7: starch high dose 
)


###############################################################################
#  LOAD AND PROCESS METABOLOMICS DATA
###############################################################################

# Raw metabolomics (micromol/gram wet weight)
metabolomics <- read_excel("Metabolomics_data_clean_micromol_gram.xlsx")
treatments   <- read_excel("TreatmentGroup_labeling.xlsx")
colnames(treatments) <- c("Group", "Treatment")

# ==========================================
# MetaboAnalyst preprocessing:
#   - Filtered: 5% IQR low-variance + 5% mean intensity low-abundance
#   - PQN normalization (reference: pooled sample)
# ==========================================

### CECUM (PQN-normalized)
metabolomics_pqn_analyst <- read.csv(
  "MetaboAnalyst/cecum_WETdata_0.05IQR_0.05medianIntensity_filtered_PQNormalized.csv",
  check.names = FALSE)
colnames(metabolomics_pqn_analyst)[1:2] <- c("Sample_Name", "Treatment")
metabolite_cols <- colnames(metabolomics_pqn_analyst)[3:14]

# Extract treatment group from sample name (e.g., "G1-1" -> "G1")
metabolomics_pqn_analyst$Treatment <- sapply(
  metabolomics_pqn_analyst$Sample_Name,
  function(x) strsplit(x, "-")[[1]][1]
)

# Compute fold changes relative to control (G1) median for each metabolite
metabolomics_pqn_analyst_fc <- metabolomics_pqn_analyst %>%
  pivot_longer(cols = all_of(metabolite_cols),
               names_to = "Metabolite", values_to = "Value") %>%
  group_by(Metabolite) %>%
  mutate(
    Control_med  = median(Value[Treatment == "G1"], na.rm = TRUE),
    Control_mean = mean(Value[Treatment == "G1"], na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    FoldChange      = Value / Control_med,
    FoldChange_mean = Value / Control_mean
  )


###############################################################################
# PROCESS 48H PREDICTIONS (2-COMPARTMENT MODEL)
###############################################################################

insilico_48h_G1 <- read.table(
  paste0(path_2comp_predictions, "Tables/conc_table_hour48_metabolomics_G1_updmods_lessIlealBac_sim20250801.txt"),
  sep = '\t', header = TRUE)
insilico_48h_treatments <- read.table(
  paste0(path_2comp_predictions, "Tables/conc_table_hour48_metabolomics_treatment_updmods_lessIlealBac_sim20250802.txt"),
  sep = '\t', header = TRUE)

# Average replicates for control (G1)
insilico_48h_G1_avg <- insilico_48h_G1 %>%
  group_by(Metabolite, Compartment, Sample, compound, compound_gper100g, Metabolite_Name) %>%
  summarise(Conc = mean(Conc, na.rm = TRUE), .groups = 'drop')

# Average replicates for treatments
insilico_48h_treatments_avg <- insilico_48h_treatments %>%
  group_by(Metabolite, Compartment, Sample, compound, compound_gper100g, Metabolite_Name) %>%
  summarise(Conc = mean(Conc, na.rm = TRUE), .groups = 'drop')

# Control medians (per metabolite and compartment)
control_medians_48h <- insilico_48h_G1_avg %>%
  group_by(Metabolite_Name, Compartment) %>%
  summarise(Control_median = median(Conc, na.rm = TRUE), .groups = 'drop')

# Compute fold changes for treatment predictions
insilico_48h_processed <- insilico_48h_treatments_avg %>%
  left_join(control_medians_48h, by = c("Metabolite_Name", "Compartment")) %>%
  mutate(
    FoldChange      = Conc / Control_median,
    Treatment       = substr(Sample, 1, 2),
    Source           = "48h predictions",
    compound_amount = paste0(compound, "_", compound_gper100g)
  ) %>%
  filter(Treatment != "G1") %>%
  dplyr::select(Treatment, Compartment, Metabolite_Name, FoldChange,
                compound_amount, compound, compound_gper100g, Source, Sample)

# Append compartment suffix to sample names for unique identification
insilico_48h_processed$Sample[insilico_48h_processed$Compartment == "Cecum"] <-
  gsub("_", "_T", insilico_48h_processed$Sample[insilico_48h_processed$Compartment == "Cecum"])



###############################################################################
#   CORRELATIONS: 48h 2-COMPARTMENT MODEL vs METABOLOMICS (CECUM)
###############################################################################

compartment_focus <- "Cecum"

# Helper: safe Spearman correlation (handles small sample sizes)
safe_spear <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(list(rho = NA_real_, p = NA_real_))
  list(
    rho = suppressWarnings(cor(x[ok], y[ok], method = "spearman")),
    p   = tryCatch(suppressWarnings(cor.test(x[ok], y[ok], method = "spearman")$p.value),
                   error = function(e) NA_real_)
  )
}

# --- Prepare experimental fold changes ---
exp_fc_cecum <- metabolomics_pqn_analyst_fc %>%
  filter(Treatment != "G1") %>%
  mutate(
    Sample_Name = gsub("-", "_", Sample_Name),
    BirdID      = str_extract(Sample_Name, "\\d+$")
  ) %>%
  dplyr::select(Treatment, Metabolite, FoldChange, FoldChange_mean, Sample_Name, BirdID)

# --- Prepare predicted fold changes (filtered to chosen doses) ---
pred_fc_cecum <- insilico_48h_processed %>%
  filter(Compartment == compartment_focus, Treatment != "G1") %>%
  inner_join(treatment_dose_map, by = c("Treatment", "compound_amount")) %>%
  rename(Metabolite = Metabolite_Name, Pred_FC = FoldChange) %>%
  mutate(
    Sample_Name = Sample,
    BirdID      = str_extract(Sample_Name, "\\d+$")
  ) %>%
  dplyr::select(Treatment, Metabolite, Pred_FC, Sample_Name, BirdID, compound_amount)

# --- Compute per-metabolite correlations (group medians across treatments) ---
exp_group <- exp_fc_cecum %>%
  group_by(Treatment, Metabolite) %>%
  summarise(Exp_FC = median(FoldChange, na.rm = TRUE), .groups = "drop")

pred_group <- pred_fc_cecum %>%
  group_by(Treatment, Metabolite) %>%
  summarise(Pred_FC = median(Pred_FC, na.rm = TRUE), .groups = "drop")

group_join <- inner_join(exp_group, pred_group, by = c("Treatment", "Metabolite"))

per_metab_corr <- group_join %>%
  group_by(Metabolite) %>%
  summarise(
    n_treatments = n(),
    .corr        = list(safe_spear(Exp_FC, Pred_FC)),
    rho_group    = .corr[[1]]$rho,
    p_group      = .corr[[1]]$p,
    agree_group  = mean((Exp_FC > 1 & Pred_FC > 1) |
                         (Exp_FC < 1 & Pred_FC < 1) |
                         (Exp_FC == 1 & Pred_FC == 1)),
    .groups = "drop"
  ) %>%
  mutate(q_group = p.adjust(p_group, method = "BH"))

# Save correlation scores
write.table(per_metab_corr %>% dplyr::select(-.corr),
            "CorrScores_cecum_chosenDoses_vs_metabolomicsPQNnorm_wet_spearman.tsv",
            sep = '\t', quote = FALSE, row.names = FALSE)


###############################################################################
#  LOAD 96H PREDICTIONS (6-COMPARTMENT MODEL)
###############################################################################
path_6comp_predictions = "Tables/"
# Metabolite-to-reaction ID mapping
metabolomics_metabolites <- c(
  "Acetic acid"           = "EX_cpd00029_e0",
  "Propionic acid"        = "EX_cpd00141_e0",
  "Butyric acid"          = "EX_cpd00211_e0",
  "Lactic acid"           = "EX_cpd00159_e0",
  "Succinic acid"         = "EX_cpd00036_e0",
  "alpha-Ketoglutaric acid" = "EX_cpd00024_e0",
  "Fumaric acid"          = "EX_cpd00106_e0"
)

# Treatment label mapping for 96h simulations
treatment_mapping <- c(
  "Cellulose_2"    = "Cellulose Low",
  "Cellulose_8"    = "Cellulose High",
  "L-Threonine_1"  = "L-Threonine Low",
  "L-Threonine_5"  = "L-Threonine High",
  "Starch_n19_5"   = "Starch Low",
  "Starch_n19_10"  = "Starch High"
)

treatment_group_mapping <- c(
  "Cellulose_2"    = "G2",
  "Cellulose_8"    = "G3",
  "L-Threonine_1"  = "G4",
  "L-Threonine_5"  = "G5",
  "Starch_n19_5"   = "G6",
  "Starch_n19_10"  = "G7"
)


insilico_conc_96h <- read.table(paste0(path_6comp_predictions,
                                       "conc_table_hour96_allMets_allCompounds_diet_grower_96h_econc1_propGizzard_corn180_wheat130_CobaltX12_rmAnaero_less02_mucUrea_wCoA_moreNightFlow_corn20240721_wheat20240721.txt"),
                                sep='\t',header=T) %>% mutate(compound_amount = paste0(compound,"_",compound_gper100g)) %>%
  filter(medium == "corn") %>%
  dplyr::select(-essential_conc,-proportion_medium,-Hour)

insilico_conc_96h_control <- read.table(paste0(path_6comp_predictions,
                                               "conc_table_hour96_allMets_control_corn_diet_grower_NEW_AFTER48h_96h_econc1_propGizzard180_CobaltX12_rmAnaero_less02_mucUreaCecaColon_wCoA_moreNightFlow_fix_sim20240719.txt"),
                                        sep='\t',header=T)

selected_compounds <- c("L-Threonine_1", "L-Threonine_5",
                        "Cellulose_2", "Cellulose_8","Starch_n19_5", "Starch_n19_10")
compartments <- c("Ileum","Cecum","Colon")

# average by replicate
insilico_conc_96h_control_sub <- insilico_conc_96h_control %>%
  filter(Metabolite %in% metabolomics_metabolites & 
           Compartment %in% compartments) %>%
  group_by(Metabolite, Compartment, SampleName)%>%
  summarize(Conc_mean_control = mean(Conc))

# calculate fold change for each sample (mean of 5 reps treatment/control)
insilico_conc_96h_fc <- insilico_conc_96h %>%
  filter(Metabolite %in% metabolomics_metabolites & 
           Compartment %in% compartments & 
           compound_amount %in% selected_compounds) %>%
  left_join(insilico_conc_96h_control_sub,
            by = c("Metabolite", "Compartment","SampleName")) %>%
  group_by(Metabolite, Compartment, SampleName, compound_amount) %>%
  summarize(FoldChange = mean(Conc) / Conc_mean_control) %>%
  distinct()


insilico_conc_96h_fc$Metabolite_Name <- 
  unlist(lapply(insilico_conc_96h_fc$Metabolite, 
                function(x) names(metabolomics_metabolites)[which(metabolomics_metabolites == x)]))

# 96h predictions 
insilico_96h_fc <- insilico_conc_96h_fc %>%
  group_by(compound_amount) %>%
  mutate(Treatment = treatment_group_mapping[compound_amount]) %>%
  dplyr::select(-Metabolite) %>%
  rename(Sample_Name = SampleName, Metabolite = Metabolite_Name) %>%
  dplyr::select(Treatment, Compartment, Metabolite, FoldChange, compound_amount, Sample_Name) %>%
  mutate(Source = "96h predictions")



###############################################################################
# 96H CORRELATIONS (3-COMPARTMENT: ILEUM, CECUM, COLON)
###############################################################################

### Load PQN-normalized metabolomics for all 3 compartments
metabolomics_pqn_analyst_cecum <- read.csv(
  "MetaboAnalyst/cecum_WETdata_0.05IQR_0.05medianIntensity_filtered_PQNormalized.csv",
  check.names = FALSE)
metabolomics_pqn_analyst_colon <- read.csv(
  "MetaboAnalyst/colon_WETdata_0.05IQR_0.05medianIntensity_filtered_PQNormalized.csv",
  check.names = FALSE)
metabolomics_pqn_analyst_ileum <- read.csv(
  "MetaboAnalyst/ileum_WETdata_0.05IQR_0.05medianIntensity_filtered_PQNormalized.csv",
  check.names = FALSE)

metabolomics_pqn_analyst_3comp <- bind_rows(
  metabolomics_pqn_analyst_ileum,
  metabolomics_pqn_analyst_cecum,
  metabolomics_pqn_analyst_colon
)

colnames(metabolomics_pqn_analyst_3comp)[1:2] <- c("Sample_Name", "Treatment")
metabolomics_pqn_analyst_3comp <- metabolomics_pqn_analyst_3comp %>%
  mutate(Compartment = case_when(
    grepl("I", Sample_Name) ~ "Ileum",
    grepl("T", Sample_Name) ~ "Cecum",
    grepl("C", Sample_Name) ~ "Colon"
  ))
metabolite_cols_3comp <- colnames(metabolomics_pqn_analyst_3comp)[3:15]

metabolomics_pqn_analyst_3comp$Treatment <- sapply(
  metabolomics_pqn_analyst_3comp$Sample_Name,
  function(x) strsplit(x, "-")[[1]][1]
)

# Compute fold changes per compartment
metabolomics_pqn_analyst_3comp_fc <- metabolomics_pqn_analyst_3comp %>%
  group_by(Compartment) %>%
  pivot_longer(cols = all_of(metabolite_cols_3comp),
               names_to = "Metabolite", values_to = "Value") %>%
  group_by(Compartment, Metabolite) %>%
  mutate(
    Control_med  = median(Value[Treatment == "G1"], na.rm = TRUE),
    Control_mean = mean(Value[Treatment == "G1"], na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    FoldChange      = Value / Control_med,
    FoldChange_mean = Value / Control_mean
  )


### Correlate 96h simulated vs experimental medians
sim_medians <- insilico_96h_fc %>%
  group_by(Compartment, Metabolite, Treatment) %>%
  summarize(FoldChange_sim = median(FoldChange, na.rm = TRUE), .groups = 'drop')

exp_medians <- metabolomics_pqn_analyst_3comp_fc %>%
  group_by(Compartment, Metabolite, Treatment) %>%
  summarize(FoldChange_exp = median(FoldChange, na.rm = TRUE), .groups = 'drop')

summary_medians <- left_join(sim_medians, exp_medians,
                             by = c("Compartment", "Metabolite", "Treatment"))

# Per-metabolite/compartment Spearman and Pearson correlations
corr_stats <- summary_medians %>%
  group_by(Compartment, Metabolite) %>%
  summarize(
    rho = if (sum(!is.na(FoldChange_sim) & !is.na(FoldChange_exp)) >= 2)
      cor(FoldChange_sim, FoldChange_exp, method = "spearman") else NA_real_,
    p = if (sum(!is.na(FoldChange_sim) & !is.na(FoldChange_exp)) >= 2)
      cor.test(FoldChange_sim, FoldChange_exp, method = "spearman")$p.value else NA_real_,
    pearson = if (sum(!is.na(FoldChange_sim) & !is.na(FoldChange_exp)) >= 2)
      cor(log2(FoldChange_sim), log2(FoldChange_exp), method = "pearson") else NA_real_,
    pearson_p = if (sum(!is.na(FoldChange_sim) & !is.na(FoldChange_exp)) >= 2)
      cor.test(log2(FoldChange_sim), log2(FoldChange_exp), method = "pearson")$p.value else NA_real_,
    .groups = "drop"
  ) %>%
  filter(!is.na(rho))

# Prepare plot data with annotation labels
plot_df <- left_join(summary_medians, corr_stats, by = c("Compartment", "Metabolite"))
plot_df$Compartment <- factor(plot_df$Compartment, levels = c("Ileum", "Cecum", "Colon"))
plot_df <- plot_df %>%
  group_by(Compartment) %>%
  mutate(
    q_spearman    = p.adjust(p, method = "BH"),
    q_pearson     = p.adjust(pearson_p, method = "BH"),
    label_spear   = sprintf("rho = %.2f\np = %.3g\nq = %.3g", rho, p, q_spearman),
    label_pearson = sprintf("R = %.2f\np = %.3g\nq = %.3g", pearson, pearson_p, q_pearson)
  )

# --- 96h correlation plot: 3 compartments, Spearman ---
corr_plot96h_spearman_3comp <- ggplot(
  plot_df,
  aes(x = log2(FoldChange_sim), y = log2(FoldChange_exp), color = Treatment)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "gray80", linetype = 2) +
  facet_wrap(Compartment ~ Metabolite, scales = "free", nrow = 6) +
  theme_bw(base_size = 12) +
  labs(x = "Simulated log2(Fold Change) (96h)",
       y = "Experimental log2(Fold Change)",
       color = "Treatment") +
  theme(strip.text = element_text(size = 10, face = "bold")) +
  geom_text(
    data = plot_df %>% group_by(Compartment, Metabolite) %>% slice_head(n = 1),
    aes(label = label_spear),
    x = -Inf, y = Inf, hjust = -0.05, vjust = 1.1, size = 4, inherit.aes = FALSE
  )

# --- 96h correlation plot: Cecum only, Spearman ---
corr_plot96h_spearman <- ggplot(
  plot_df[plot_df$Compartment == "Cecum", ],
  aes(x = log2(FoldChange_sim), y = log2(FoldChange_exp), color = Treatment)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "gray80", linetype = 2) +
  facet_wrap(. ~ Metabolite, scales = "free", nrow = 6) +
  theme_bw(base_size = 12) +
  labs(x = "Simulated log2(Fold Change) (96h)",
       y = "Experimental log2(Fold Change)",
       color = "Treatment") +
  theme(strip.text = element_text(size = 10, face = "bold")) +
  geom_text(
    data = plot_df %>% group_by(Metabolite) %>% slice_head(n = 1),
    aes(label = label_spear),
    x = -Inf, y = Inf, hjust = -0.05, vjust = 1.1, size = 4, inherit.aes = FALSE
  )

# --- 96h correlation plot: Cecum only, Pearson ---
corr_plot96h_pearson <- ggplot(
  plot_df[plot_df$Compartment == "Cecum", ],
  aes(x = log2(FoldChange_sim), y = log2(FoldChange_exp), color = Treatment)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "gray80", linetype = 2) +
  facet_wrap(. ~ Metabolite, scales = "free", nrow = 6) +
  theme_bw(base_size = 12) +
  labs(x = "Simulated log2(Fold Change) (96h)",
       y = "Experimental log2(Fold Change)",
       color = "Treatment") +
  theme(strip.text = element_text(size = 10, face = "bold")) +
  geom_text(
    data = plot_df %>% group_by(Metabolite) %>% slice_head(n = 1),
    aes(label = label_pearson),
    x = -Inf, y = Inf, hjust = -0.05, vjust = 1.1, size = 4, inherit.aes = FALSE
  )

ggsave(filename = "CorrelationPlots_13mets_3compartments_96hSim_vs_metabolomics_wet_pqnNorm_medians_Spearman.pdf",
       corr_plot96h_spearman_3comp, width = 14, height = 20)

ggsave(filename = "CorrelationPlots_13mets_Cecum_96hSim_vs_metabolomics_wet_pqnNorm_medians_Spearman.pdf",
       corr_plot96h_spearman, width = 12, height = 10)


###############################################################################
# FIGURE 4A: 7 METABOLITES, CECUM, BOTH MODELS (48h + 96h)
###############################################################################

#### LOAD 48h PREDICTIONS FC DATA
insilico_48h_processed <- read.table(
  "InsilicoPredictions_48h_2compModel_FCmedians_03082025.tsv",
  sep = '\t', header = TRUE)

metabolites_plot <- c(
  "Acetic acid", "Propionic acid", "Butyric acid", "Lactic acid",
  "Succinic acid", "alpha-Ketoglutaric acid", "Fumaric acid"
)

# --- 48h predictions: filter to chosen doses and compute treatment medians ---
insilico_48h_for_plot <- insilico_48h_processed %>%
  filter(Compartment == "Cecum", Treatment != "G1",
         Metabolite_Name %in% metabolites_plot) %>%
  inner_join(treatment_dose_map, by = c("Treatment", "compound_amount")) %>%
  group_by(Treatment, Metabolite_Name) %>%
  summarise(Pred_FC = median(FoldChange, na.rm = TRUE), .groups = "drop") %>%
  rename(Metabolite = Metabolite_Name) %>%
  mutate(source = "2-compartment model")


# --- 96h predictions: cecum, 7 metabolites ---
insilico_96h_for_plot <- insilico_96h_fc %>%
  filter(Compartment == "Cecum",
         Metabolite %in% metabolites_plot) %>%
  group_by(Metabolite, Treatment) %>%
  summarize(Pred_FC = median(FoldChange, na.rm = TRUE), .groups = 'drop') %>%
  mutate(source = "6-compartment model")


# --- Experimental metabolomics (cecum, 7 metabolites) ---
metabolomics_pqn_analyst_fc_7mets <- read.csv(
  "metabolomics_cecum_0.05IQR_0.05medianIntensity_filtered_PQNorm_FoldChanges.csv") %>%
  filter(Metabolite %in% metabolites_plot, Treatment != "G1") %>%
  group_by(Treatment, Metabolite) %>%
  summarise(Exp_FC = median(FoldChange, na.rm = TRUE), .groups = "drop")


# --- Combine predictions from both models with experimental data ---
treatment_colors <- c(
  "G1" = "#2C3E50", "G2" = "#E74C3C", "G3" = "darkred",
  "G4" = "#3498DB", "G5" = "#2f2fad",
  "G6" = "#27AE60", "G7" = "#216907"
)

plot_data <- inner_join(
  rbind(insilico_48h_for_plot, insilico_96h_for_plot),
  metabolomics_pqn_analyst_fc_7mets,
  by = c("Treatment", "Metabolite")
)

# Compute facet-level correlation statistics
facet_labels <- plot_data %>%
  group_by(source, Metabolite) %>%
  summarise(
    .corr = list({
      ok <- is.finite(Exp_FC) & is.finite(Pred_FC)
      rho <- suppressWarnings(cor(Exp_FC[ok], Pred_FC[ok], method = "spearman"))
      p   <- tryCatch(cor.test(Exp_FC[ok], Pred_FC[ok], method = "spearman")$p.value,
                      error = function(e) NA_real_)
      pearson <- suppressWarnings(cor(log2(Exp_FC[ok]), log2(Pred_FC[ok]), method = "pearson"))
      pearson_p <- tryCatch(cor.test(log2(Exp_FC[ok]), log2(Pred_FC[ok]), method = "pearson")$p.value,
                            error = function(e) NA_real_)
      tibble(rho = rho, p = p, pearson = pearson, pearson_p = pearson_p)
    }),
    .groups = "drop"
  ) %>%
  tidyr::unnest(.corr) %>%
  mutate(
    q = p.adjust(p, method = "BH"),
    facet_label = sprintf("rho=%.2f\np=%.3g\nq=%.3g", rho, p, q)
  )


# --- Generate Figure 4A ---
plot <- ggplot(plot_data,
               aes(log2(Pred_FC), log2(Exp_FC), color = Treatment, label = Treatment)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "gray80", linetype = 2) +
  ggrepel::geom_text_repel(min.segment.length = 0, size = 3) +
  facet_wrap(source ~ Metabolite, scales = "free", nrow = 2) +
  labs(
    x = "Predicted log2(Fold Change)",
    y = "Experimental log2(Fold Change)"
  ) +
  theme_bw() +
  scale_color_manual(values = treatment_colors) +
  theme(
    strip.text       = element_text(face = "bold"),
    strip.background = element_rect(fill = "transparent"),
    legend.position  = "none",
    plot.title       = element_text(size = 11)
  ) +
  geom_text(
    data = facet_labels, inherit.aes = FALSE,
    aes(x = -Inf, y = Inf, label = facet_label),
    hjust = -0.05, vjust = 1.1, size = 3.3, color = "black"
  )

ggsave(filename = "Figure4A_CorrelationPlots_7mets_Cecum_SimData48h_96h_vs_metabolomics_wet_pqnNorm_medians_Spearman.pdf",
       plot, width = 15, height = 6.5)
