# Figure 3: Bubble plot of SCFA flux fold changes (lower GIT) and bar chart
# showing the number of compounds with significant increase/decrease in SCFA production.

library(cowplot)
library(ggplot2)
library(dplyr)
library(reshape2)
library(magrittr)
library(ggpp)
library(tidyr)
library(ggrepel)

## load flux tables with wilcoxon pvalues:
setwd("WithPrebiotics/2024/")

path = "Tables/"

pattern_arg_compounds = "diet_grower_96h_econc1_propGizzard_corn180_wheat130_CobaltX12_rmAnaero_less02_mucUrea_wCoA_moreNightFlow"
pattern_arg_control = "diet_grower_96h_econc1_propGizzard_corn180_wheat130_CobaltX12_rmAnaero_less02_mucUreaCecaColon_wCoA_moreNightFlow_fix"
sim_dates_compounds = "corn20240721_wheat20240721"
sim_dates_control = "corn20240719_wheat20240719"

Compartments_list <- list(c("Gizzard","Duodenum","Jejunum","Ileum","Cecum", "Colon"))

# load compound fluxes
flux_table_allReps_96h <- read.delim(paste0(path,"FluxTable_ToPlot_PairedWilcoxPvalsBH_FoldChange_allCompounds_",
                                            sim_dates_compounds,"_vs_control_", sim_dates_control, ".txt"),
                                     sep ='\t', header=T)

# subset to SCFA only:
flux_table_allReps_96h_scfa <- flux_table_allReps_96h[which(flux_table_allReps_96h$Metabolite %in% c("Acetate","Propionate","Butyrate")),]

# remove faulty L_Threonine_3 and 4:
flux_table_allReps_96h_scfa <- flux_table_allReps_96h_scfa[which(!(flux_table_allReps_96h_scfa$compound_amount %in% c("L-Threonine_3", "L-Threonine_4"))),]

scfa_colors <- setNames(c("#ab3613","#576991","#702963"), c("Acetate","Propionate","Butyrate"))


#########################################################
######    MAKE BUBBLE PLOT FOR SCFA FOLD CHANGES   ######
#########################################################
# First reshape the data to wide format for plotting
scfa_wide <- flux_table_allReps_96h_scfa %>%
  select(Compartment, compound_amount, Metabolite, median_fold_change, pval_adj) %>%
  distinct() %>%
  tidyr::pivot_wider(
    names_from = Metabolite,
    values_from = c(median_fold_change, pval_adj)
  ) %>%
  rename(
    fc_Acetate = median_fold_change_Acetate,
    fc_Butyrate = median_fold_change_Butyrate,
    fc_Propionate = median_fold_change_Propionate,
    sig_Acetate = pval_adj_Acetate,
    sig_Butyrate = pval_adj_Butyrate,
    sig_Propionate = pval_adj_Propionate
  )

# Add a "score" column to identify "best" compounds
scfa_wide <- scfa_wide %>%
  mutate(
    # Count significant increases
    score = (sig_Acetate <= 0.05 & fc_Acetate > 1 & fc_Propionate > 0.9) +
      (sig_Butyrate <= 0.05 & fc_Butyrate > 1 & fc_Propionate > 0.9) +
      (sig_Propionate <= 0.05 & fc_Propionate > 1),
    # Calculate average fold change for significant increases
    avg_fc = ((fc_Acetate * (sig_Acetate <= 0.05 & fc_Acetate > 1  & fc_Propionate > 0.9)) +
                (fc_Butyrate * (sig_Butyrate <= 0.05 & fc_Butyrate > 1  & fc_Propionate > 0.9)) +
                (fc_Propionate * (sig_Propionate <= 0.05 & fc_Propionate > 1))) /
      (score + (score == 0))  # Avoid division by zero
  ) %>%
  # Add top rank based on both criteria
  group_by(Compartment) %>%
  mutate(
    is_top = dense_rank(desc(score * avg_fc)) <= 10
  ) %>%
  ungroup()

# remove 'is_top' if score == 0 (too few score > 0 in jejunum -> all values got 'is_top')
scfa_wide$is_top[scfa_wide$score == 0] <- FALSE


# add significance column
scfa_wide <- scfa_wide %>%
  mutate(
    significance = factor(case_when(
      sig_Acetate <= 0.05 & sig_Butyrate <= 0.05 & sig_Propionate <= 0.05 ~ "All",
      sig_Acetate <= 0.05 & sig_Butyrate <= 0.05 ~ "Ac+But",
      sig_Acetate <= 0.05 & sig_Propionate <= 0.05 ~ "Ac+Prop",
      sig_Butyrate <= 0.05 & sig_Propionate <= 0.05 ~ "But+Prop",
      sig_Acetate <= 0.05 ~ "Acetate",
      sig_Butyrate <= 0.05 ~ "Butyrate",
      sig_Propionate <= 0.05 ~ "Propionate",
      TRUE ~ "None"
    ), levels = c("All", "Ac+But", "Ac+Prop", "But+Prop",
                  "Acetate", "Butyrate", "Propionate", "None"))
  )

write.table(scfa_wide, "~/Desktop/Study/UofT/ParkinsonLab/Modeling_chickens/ModelsAnalyses_outputs/SCFA_FluxFCs_significanceLevels_table_for_bubblePlot.txt",
            sep='\t', quote=F)

significance_label_colors <- c(
  "All" = "darkred",
  "Ac+But" = "#b5a900",
  "Ac+Prop" = "darkorange",
  "But+Prop" = "green",
  "Acetate" = scfa_colors["Acetate"],
  "Butyrate" = scfa_colors["Butyrate"],
  "Propionate" = scfa_colors["Propionate"],
  "None" = "darkgray"
)

scfa_wide_lowerGIT <- scfa_wide[which(scfa_wide$Compartment %in% c("Ileum","Cecum","Colon")),]

# Create the bubble plot
bubble_plot <- ggplot(scfa_wide_lowerGIT,
                      aes(x = fc_Butyrate,
                          y = fc_Propionate)) +
  geom_point(aes(size = fc_Acetate,
                 fill = significance),
             shape = 21,
             color = "gray80",
             stroke=0.5,
             alpha = 0.7) +
  # Add labels for top compounds
  geom_text_repel(
    data = . %>% filter(is_top),
    aes(label = compound_amount),
    size = 5,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.alpha = 0.7,
    segment.curvature = -0.1,
    point.padding = 0.5,
    box.padding = 0.5,
    color = sapply(scfa_wide_lowerGIT$significance[scfa_wide_lowerGIT$is_top], function(x) significance_label_colors[x])
  ) +
  # Reference lines
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  # Facet by compartment
  facet_wrap(~factor(Compartment, levels = Compartments_list[[1]]), scales = "free_x") +
  scale_size_continuous(
    name = "Acetate\nFold Change",
    range = c(2, 20),
    trans = "exp"
  ) +
  scale_fill_manual(
    name = "Significant Changes",
    values = c("red","yellow","orange","green","#ab3613","#702963","#576991","gray90"),
    drop = FALSE
  ) +
  labs(
    x = "Butyrate Fold Change",
    y = "Propionate Fold Change",
    title = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.title = element_text(size=24),
    axis.text = element_text(size=20),
    strip.text = element_text(size = 26),
    legend.text = element_text(size = 16),
    legend.title = element_text(size=20),
    panel.grid.minor = element_blank()
  )

pdf("~/Desktop/Study/UofT/ParkinsonLab/Modeling_chickens/ModelsAnalyses_outputs/Plots/FluxFC_SCFA_bubblePlot_ButX_PropY_AcSize_lowerGIT_v191125.pdf",
    width = 21, height = 7)
bubble_plot
dev.off()



############################################################################
###########            MAKE QUANTITATIVE PANELS         ####################
###  WITH N OF COMPOUNDS THAT INCREASED/DECREASED/NO CHANGE PRODUCTION  ####
############################################################################

scfa_wide <- read.table("~/Desktop/Study/UofT/ParkinsonLab/Modeling_chickens/ModelsAnalyses_outputs/SCFA_FluxFCs_significanceLevels_table_for_bubblePlot.txt",
                        sep='\t', header= T)

# create complete categorization of all compounds
compound_categories <- scfa_wide %>%
  mutate(
    Acetate_cat = case_when(
      fc_Acetate > 1.02 & sig_Acetate <= 0.05 ~ "Significant increase >2%",
      fc_Acetate < 0.98 & sig_Acetate <= 0.05 ~ "Significant decrease >2%",
      TRUE ~ "Non-significant change"
    ),
    Butyrate_cat = case_when(
      fc_Butyrate > 1.02 & sig_Butyrate <= 0.05 ~ "Significant increase >2%",
      fc_Butyrate < 0.98 & sig_Butyrate <= 0.05 ~ "Significant decrease >2%",
      TRUE ~ "Non-significant change"
    ),
    Propionate_cat = case_when(
      fc_Propionate > 1.02 & sig_Propionate <= 0.05 ~ "Significant increase >2%",
      fc_Propionate < 0.98 & sig_Propionate <= 0.05 ~ "Significant decrease >2%",
      TRUE ~ "Non-significant change"
    )
  )


compound_categories_long <- compound_categories %>%
  pivot_longer(
    cols = ends_with("_cat"),
    names_to = "SCFA",
    values_to = "Category"
  ) %>%
  mutate(SCFA = gsub("_cat", "", SCFA))


barchart_3categories <- ggplot(compound_categories_long[which(compound_categories_long$Compartment %in% Compartments_list[[1]][4:6]),],
                               aes(x = factor(SCFA, levels = c("Acetate","Propionate","Butyrate")), fill = Category)) +
  geom_bar(position = "stack") +
  facet_wrap(~factor(Compartment, levels = Compartments_list[[1]])) +
  scale_fill_manual(values = rev(wesanderson::wes_palettes$Darjeeling1[c(2,1,6)])) +
  scale_y_continuous(
    breaks = seq(0, 64, by = 10),
    limits = c(0, 64)
  ) +
  theme_classic() +
  labs(
    title = "Distribution of SCFA production changes",
    y = "Number of compounds",
    x = "SCFA"
  ) +
  theme(axis.text.x = element_text(size = 13, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 15),
        legend.text = element_text(size=12),
        axis.title.x = element_blank(),
        strip.text = element_text(size=14),
        axis.text.y = element_text(size=12))

pdf("Plots/BarChart_NofCompounds_quantified_FluxFC_SCFA_3categories_lowerGIT.pdf",
    width = 7, height = 5)
barchart_3categories
dev.off()
