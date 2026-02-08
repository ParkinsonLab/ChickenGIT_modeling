# Figure 2 and Figure S2: metabolite concentration PCoA, line/dot plots,
# concentration heatmaps, SCFA flux bar plots by species, and EC heatmaps
# annotated by KEGG/MetaCyc pathways.

library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(BacArena)
detach("package:plyr", unload = TRUE)


# =========================================================================
#  Shared configuration
# =========================================================================
compartment_levels <- c("Gizzard", "Duodenum", "Jejunum", "Ileum", "Cecum", "Colon")

path_outfiles <- "2024/"

diet_corn  <- read.table("Diets_chicken/corn_diet_starter_aggregated_2023.txt",
                          sep = '\t', header = TRUE)
diet_wheat <- read.table("Diets_chicken/wheat_diet_starter_aggregated_2023.txt",
                          sep = '\t', header = TRUE)

seed_metabolites <- readr::read_delim("all_seed_metabolites_edited.tsv",
                                       delim = '\t', col_names = TRUE)

# Reusable compartment color palette (Darjeeling palettes)
darj1 <- wesanderson::wes_palette("Darjeeling1")
darj2 <- wesanderson::wes_palette("Darjeeling2")

compartment_colors <- c(
  Gizzard  = darj2[1],
  Duodenum = "#7a603d",
  Jejunum  = darj1[4],
  Ileum    = darj1[1],
  Cecum    = darj1[5],
  Colon    = darj2[2]
)


# =========================================================================
#  Load 96h concentration data
# =========================================================================
conc_96h_raw <- read.table(
  paste0(path_outfiles,
         "Tables/conc_table_96h_allMets_control_diet_grower_96h_econc1_propGizzard_",
         "corn180_wheat130_CobaltX12_rmAnaero_less02_mucUreaCecaColon_wCoA_",
         "moreNightFlow_fix_corn20240719_wheat20240719.txt"),
  header = TRUE
)

# Map metabolite IDs to human-readable names
conc_96h_raw$MetName <- sapply(conc_96h_raw$Metabolite, function(met) {
  met_id <- gsub("EX_|_e0", "", met)
  match_idx <- match(met_id, seed_metabolites$id)
  if (!is.na(match_idx)) as.character(seed_metabolites$name[match_idx]) else NA
})

# Average across replicates and convert to mM
conc_96h_mean_mM <- conc_96h_raw %>%
  group_by(across(-c(Conc, Replicate))) %>%
  summarise(Conc = mean(Conc), .groups = "drop") %>%
  mutate(
    Conc        = Conc / 625,
    Hour        = factor(as.numeric(as.character(Hour)), levels = 2:95),
    Compartment = factor(Compartment, levels = compartment_levels)
  )


# =========================================================================
#  FIGURE 2B: Combined PCoA of metabolite concentrations at 3 timepoints
# =========================================================================
hours_pcoa <- c(25, 35, 69)

conc_pcoa_subset <- conc_96h_mean_mM %>%
  filter(Hour %in% hours_pcoa |
           (Hour == 27 & Compartment == "Gizzard")) %>%
  mutate(sample_key = paste(Compartment, SampleName, Hour, sep = "_"))

# Build wide matrix for Bray-Curtis dissimilarity
conc_wide <- conc_pcoa_subset %>%
  select(sample_key, MetName, Conc) %>%
  tidyr::pivot_wider(names_from = MetName, values_from = Conc, values_fn = list(Conc = mean))

conc_matrix <- as.data.frame(conc_wide)
rownames(conc_matrix) <- conc_matrix$sample_key
conc_matrix$sample_key <- NULL
conc_matrix <- as.matrix(conc_matrix)
conc_matrix[is.na(conc_matrix)] <- 0
conc_matrix <- conc_matrix[, colSums(conc_matrix != 0) > 0]

# PCoA via Bray-Curtis
bray_dist     <- vegan::vegdist(conc_matrix, method = "bray")
pcoa_result   <- cmdscale(bray_dist, eig = TRUE, k = 2)
pct_var       <- pcoa_result$eig / sum(pcoa_result$eig) * 100

pcoa_df <- data.frame(
  PC1 = pcoa_result$points[, 1],
  PC2 = pcoa_result$points[, 2],
  sample_key = rownames(pcoa_result$points),
  stringsAsFactors = FALSE
) %>%
  tidyr::separate(sample_key, into = c("Compartment", "Sample", "Hour"),
                  sep = "_", remove = FALSE) %>%
  mutate(
    Diet            = gsub("[0-9]|-", "", Sample),
    Compartment_Hour = paste0(Compartment, "_", Hour),
    Diet_Hour        = paste0(Diet, "_", Hour),
    Hour             = factor(Hour, levels = c("27", "25", "35", "69")),
    Compartment_display = ifelse(Compartment_Hour == "Gizzard_27",
                                 "Diet", Compartment)
  )

# Define factor levels for Compartment_Hour
comp_hour_levels <- c(
  "Gizzard_27",
  paste0(rep(compartment_levels, each = 3), "_", rep(c(25, 35, 69), 6))
)
pcoa_df$Compartment_Hour <- factor(pcoa_df$Compartment_Hour, levels = comp_hour_levels)

# Keep only Gizzard at hour 27, all compartments for 25/35/69
pcoa_df_filtered <- pcoa_df %>%
  filter(Hour != "27" | Compartment == "Gizzard")

# Colors: each compartment repeated for the three hours + diet at hour 27
color_comp_hour <- setNames(
  c("magenta", rep(unname(compartment_colors), 3)),
  c("Gizzard_27",
    paste0(names(compartment_colors), "_25"),
    paste0(names(compartment_colors), "_35"),
    paste0(names(compartment_colors), "_69"))
)

shape_diet_hour <- c(
  Corn_27 = 23, Wheat_27 = 23,   # diamond for diet
  Corn_25 = 21, Corn_35  = 22, Corn_69  = 24,
  Wheat_25 = 21, Wheat_35 = 22, Wheat_69 = 24
)

linetype_hour <- c("27" = "dashed", "25" = "dotted", "35" = "longdash", "69" = "solid")

pcoa_plot <- ggplot(pcoa_df_filtered, aes(x = PC1, y = PC2, fill = Compartment_Hour)) +
  geom_point(aes(shape = Diet_Hour), size = 3, alpha = 0.7) +
  stat_ellipse(aes(group = Compartment_Hour, color = Compartment_Hour, linetype = Hour),
               type = "t", size = 1.1, level = 0.95) +
  scale_fill_manual(values = color_comp_hour) +
  scale_color_manual(values = color_comp_hour, guide = "none") +
  scale_shape_manual(
    values = shape_diet_hour,
    labels = c(Corn_27 = "Corn diet", Wheat_27 = "Wheat diet",
               Corn_25 = "Corn 25h", Corn_35 = "Corn 35h", Corn_69 = "Corn 69h",
               Wheat_25 = "Wheat 25h", Wheat_35 = "Wheat 35h", Wheat_69 = "Wheat 69h")
  ) +
  scale_linetype_manual(
    values = linetype_hour,
    labels = c("27" = "Diet",
               "25" = "beginning of feeding period after fasting",
               "35" = "middle of the feeding period",
               "69" = "middle of the fasting period (night)")
  ) +
  labs(
    title = "Metabolic profiles across timepoints (25, 35, 69 hours)",
    x = paste0("PCoA 1 (", format(pct_var[1], digits = 2), "%)"),
    y = paste0("PCoA 2 (", format(pct_var[2], digits = 2), "%)"),
    fill = "Compartment", shape = "Diet & Time"
  ) +
  theme_bw() +
  theme(
    legend.title    = element_text(size = 14, face = "bold"),
    title           = element_text(size = 18),
    panel.grid.major = element_blank(),
    axis.text       = element_text(size = 14),
    legend.text     = element_text(size = 12),
    axis.title      = element_text(size = 17)
  )

pdf("PCoA_ConcMets_Combined_25_35_69h_compartments_BrayCurtis_dietSameShape_thick.pdf",
    width = 20, height = 11)
print(pcoa_plot)
dev.off()


# =========================================================================
#  FIGURE 2C: Line/dot plots of metabolite concentrations across compartments
# =========================================================================
hours_line <- c(25, 35, 44)

conc_line_mean <- conc_96h_mean_mM %>%
  filter(Hour %in% hours_line) %>%
  group_by(Hour, Compartment, medium, MetName) %>%
  summarise(Conc_mean = mean(Conc), Conc_sd = sd(Conc), .groups = "drop")

mets_shortlist <- c("cellulose", "D-Glucose", "Taurine", "L-Glutamine",
                     "Pantothenic acid", "Glycochenodeoxycholate", "Hexanoate",
                     "Formate", "Acetate", "GABA")

hour_colors <- c("25" = "#ad6da3", "35" = "#48a1b5", "44" = "#330901")

# --- Corn-diet only ---
plot_line_corn <- ggplot(
  conc_line_mean %>% filter(medium == "corn", MetName %in% mets_shortlist),
  aes(x = Compartment, y = Conc_mean, color = Hour,
      group = interaction(medium, Hour))
) +
  geom_point(size = 1) +
  geom_line() +
  geom_errorbar(aes(ymin = Conc_mean - Conc_sd, ymax = Conc_mean + Conc_sd), width = 0.1) +
  facet_wrap(~ factor(MetName, levels = mets_shortlist), scales = "free", ncol = 5) +
  scale_color_manual(values = hour_colors) +
  labs(y = "Concentration, mM") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.spacing    = unit(0.5, "lines"),
    axis.text.y      = element_text(size = 13),
    axis.text.x      = element_text(size = 12, angle = 45, hjust = 1, vjust = 1.1),
    title            = element_text(size = 14),
    axis.title.x     = element_blank(),
    strip.text       = element_text(size = 13)
  )

pdf("ConcentrationChanges_lineChart_6comp_h25_35_44_10mets_cornDiet.pdf",
    height = 6, width = 16)
plot_line_corn
dev.off()

# --- Both diets ---
diet_colors <- c(wheat = "#3B270C", corn = "#03AC13")

plot_line_both <- ggplot(
  conc_line_mean %>% filter(MetName %in% mets_shortlist),
  aes(x = Compartment, y = Conc_mean, color = medium,
      group = interaction(medium, Hour), linetype = Hour)
) +
  geom_point(size = 1) +
  geom_line() +
  geom_errorbar(aes(ymin = Conc_mean - Conc_sd, ymax = Conc_mean + Conc_sd), width = 0.1) +
  facet_wrap(~ factor(MetName, levels = mets_shortlist), scales = "free", ncol = 4) +
  scale_color_manual(values = diet_colors) +
  scale_linetype_manual(values = c("solid", "longdash", "dotted")) +
  labs(y = "Concentration, mM") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.spacing    = unit(0.5, "lines"),
    axis.text.y      = element_text(size = 13),
    axis.text.x      = element_text(size = 12, angle = 45, hjust = 1, vjust = 1.1),
    title            = element_text(size = 14),
    axis.title.x     = element_blank(),
    strip.text       = element_text(size = 13)
  )

pdf("ModelsAnalyses_outputs/ConcentrationChanges_lineChart_6comp_h25_35_44_12mets_bothDiets.pdf",
    height = 10, width = 13)
plot_line_both
dev.off()


# =========================================================================
#  FIGURE 2D: Concentration heatmap across compartments and time
# =========================================================================
generate_conc_heatmap <- function(conc_table, sample_name, timesteps,
                                   mets_to_plot = c("topVarMets", "topChangedMets",
                                                     "diet_mets", "custom"),
                                   custom_mets = NULL, scale_factor = 1,
                                   legend_direction = "vertical",
                                   conc_units = "mM") {
  mets_to_plot <- match.arg(mets_to_plot)
  conc_sub <- conc_table %>% filter(SampleName == sample_name)

  diet_name <- tolower(stringr::str_remove(sample_name, "[0-9]{1,2}-[0-9]"))
  diet <- if (diet_name == "wheat") diet_wheat else diet_corn

  # Select metabolites based on method
  selected_mets <- switch(mets_to_plot,
    topVarMets = {
      conc_sub %>%
        group_by(Compartment, MetName) %>%
        summarise(var = var(Conc), .groups = "drop") %>%
        top_n(30, var) %>%
        pull(MetName) %>% unique() %>% as.character()
    },
    topChangedMets = {
      conc_sub %>%
        group_by(Compartment, MetName) %>%
        summarise(var = var(Conc), ratio = max(Conc) / min(Conc + 1e-20),
                  .groups = "drop") %>%
        top_n(35, ratio) %>% top_n(15, var) %>%
        pull(MetName) %>% unique() %>% as.character()
    },
    diet_mets = as.character(diet$met_name),
    custom    = custom_mets
  )

  conc_mat <- conc_sub %>%
    filter(MetName %in% selected_mets) %>%
    reshape2::acast(MetName ~ Compartment + Hour, value.var = "Conc", fun.aggregate = mean)

  # Column annotation
  col_info <- data.frame(
    Compartment = sapply(colnames(conc_mat), function(x) strsplit(x, "_")[[1]][1]),
    Hour        = rep(2:timesteps, 6)
  )
  col_info$Compartment <- factor(col_info$Compartment, levels = compartment_levels)

  # Compartment annotation colors (Darjeeling palettes)
  annot_comps <- HeatmapAnnotation(
    compartment = col_info$Compartment,
    col = list(compartment = c(
      Gizzard = darj2[1], Duodenum = darj2[3], Jejunum = darj1[4],
      Ileum = darj1[1], Cecum = darj1[5], Colon = darj2[2]
    )),
    show_annotation_name = FALSE, show_legend = FALSE
  )

  # Concentration color scale
  green_ramp <- c("white", "#f1f2ce", "#eaf2c4", "#e0f2c4", "#d6db93",
                   "#bde6ae", "#a4edae", "#8ad4a0", "#59a885",
                   "#1a6e48", "#0b5936", "#0a3d0e", "#032105")

  if (conc_units == "mmol/cell") {
    conc_breaks <- c(0, 1e-07, 1e-06, 1e-05, 1e-04, 1e-03, 1e-02, 0.1, 1, 10, 100)
    conc_col <- colorRamp2(conc_breaks, green_ramp[c(1:4, 6:8, 9:11)])
  } else {
    conc_breaks <- scale_factor * c(0, 1e-09, 1e-08, 1e-07, 5e-06, 1e-05,
                                     5e-05, 1e-04, 5e-04, 1e-03, 5e-03, 0.05, 0.1)
    conc_col <- colorRamp2(conc_breaks, green_ramp)
  }

  # Legend
  legend_col_fun <- colorRamp2(0:12, green_ramp)
  legend_labels  <- scale_factor * c(0, 1e-09, 1e-08, 1e-07, 5e-06, 1e-05,
                                      5e-05, 1e-04, 5e-04, 1e-03, 5e-03, 0.05, 0.5)
  conc_legend <- Legend(col_fun = legend_col_fun, labels = legend_labels,
                         title = "concentration, mM", direction = legend_direction)

  hm <- Heatmap(
    conc_mat,
    show_row_names = TRUE, show_column_names = FALSE,
    column_order = rownames(col_info$Hour),
    row_order = selected_mets,
    border = TRUE, show_column_dend = FALSE,
    cluster_rows = FALSE, show_row_dend = FALSE,
    show_heatmap_legend = FALSE,
    column_split = col_info$Compartment,
    col = conc_col,
    top_annotation = annot_comps,
    name = sample_name
  )

  list(heatmap = hm, mets = selected_mets, legend = conc_legend)
}

custom_mets_heatmap <- c(
  "D-Glucose", "D-Fructose", "Sucrose",
  "starch (n=19, 3xalpha1-6, 15xalpha1-4)", "malto-11-ose",
  "malto-6-ose", "Amylotriose",
  "L-Methionine", "L-Threonine", "L-Alanine",
  "Ethanol",
  "L-Lactate", "Propionate", "Butyrate", "O2"
)

hm_result <- generate_conc_heatmap(
  conc_table = conc_96h_mean_mM, sample_name = "Corn5-3",
  timesteps = 95, scale_factor = 2.5,
  mets_to_plot = "custom", custom_mets = custom_mets_heatmap,
  conc_units = "mM"
)

sim_plots_dir <- paste0(
  "simulation_outputs/Coupling_multiple_compartments_D10/",
  "Simulations_WithAbsorption/NoPathogens_10samples/2024/Plots/"
)

pdf(paste0(sim_plots_dir,
           "HeatMap_ConcChanges96h_Corn53_control_customMets_sim190724.pdf"),
    height = 5, width = 16)
hm_result$heatmap
dev.off()

hm_legend <- generate_conc_heatmap(
  conc_table = conc_96h_mean_mM, sample_name = "Corn5-3",
  timesteps = 95, scale_factor = 2.5,
  mets_to_plot = "custom", custom_mets = custom_mets_heatmap,
  legend_direction = "horizontal", conc_units = "mM"
)$legend

pdf(paste0(sim_plots_dir,
           "LegendHorizontal_for_HeatMap_ConcChanges96h_Corn53_control_customMets_sim190724.pdf"),
    height = 4, width = 7)
draw(hm_legend)
dev.off()


# =========================================================================
#  FIGURE 2A: SCFA flux bar plots by species (Ileum, Cecum, Colon)
# =========================================================================
flux_scfa_by_species <- read.table(
  paste0(path_outfiles,
         "Tables/flux_table_sum_SCFA_bySpecie_meanByRep_control_diet_grower_96h_econc1_",
         "propGizzard_corn180_wheat130_CobaltX12_rmAnaero_less02_mucUreaCecaColon_wCoA_",
         "moreNightFlow_fix_corn20240719_wheat20240719.txt"),
  sep = '\t', header = TRUE
)

flux_scfa_by_species$MetName <- sapply(flux_scfa_by_species$Reaction, function(x) {
  as.character(seed_metabolites$name[seed_metabolites$id == gsub("EX_|_e0", "", x)])
})
flux_scfa_by_species$Compartment <- factor(flux_scfa_by_species$Compartment,
                                            levels = compartment_levels)
flux_scfa_by_species <- flux_scfa_by_species %>% filter(MetName != "L-Lactate")

species_colors <- c(
  Blautia_sp                          = "#542294",
  Mediterraneibacter_sp002161355      = "#412566",
  Lachnoclostridium_phocaeensis       = "#2b1e3d",
  Clostridium_saudiense               = colpal3[12],
  Faecalibacterium_sp002160895        = "#d95dc4",
  Agathobaculum_sp                    = "#9c3b8b",
  Flavonifractor_plautii              = "#8f1d7b",
  Oscillospiraceae_genusUBA9475       = "#73015f",
  Romboutsia_sp                       = colpal3[9],
  Borkfalkia_sp                       = "#FFD4EB",
  Lactobacillus_crispatus             = "#f72548",
  Limosilactobacillus_reuteri         = "#bd112e",
  Lactobacillus_johnsonii             = "#99021b",
  Streptococcus_equinus               = colpal3[19],
  Staphylococcus_saprophyticus        = "#f58971",
  Enterococcus_cecorum                = "#e05536",
  Enterococcus_faecalis               = "#ab4d38",
  Erysipelatoclostridium_sp002160495  = "#BA5D00",
  Massiliomicrobiota_sp002160815      = "#8a4601",
  Megasphaera_stantonii               = colpal3[15],
  Escherichia_flexneri                = "#0569ab",
  Klebsiella_pneumoniae               = "#044d7d",
  Sphingomonas_sp                     = colpal3[4],
  Acinetobacter_lwoffii               = "#e6ba7a",
  Pseudomonas_aeruginosa              = "#b08b56",
  Bifidobacterium_gallinarum          = "#bfb5ae",
  Corynebacterium_amycolatum          = "#a19892",
  Corynebacterium_kroppenstedtii      = "#87807b",
  Phocaeicola_vulgatus                = "#595653"
)

lower_git <- c("Ileum", "Cecum", "Colon")
scfa_names <- c("Acetate", "Propionate", "Butyrate")

flux_bar_plots <- lapply(scfa_names, function(scfa) {
  ggplot(
    flux_scfa_by_species %>% filter(Compartment %in% lower_git, MetName == scfa),
    aes(x = factor(Compartment, levels = rev(lower_git)), y = Flux.mean, fill = Specie)
  ) +
    geom_bar(position = "stack", stat = "identity", alpha = 0.9) +
    scale_fill_manual(values = species_colors) +
    facet_grid(SampleName + medium + MetName ~ ., scales = "free") +
    coord_flip() +
    theme_bw() +
    theme(
      legend.title    = element_blank(),
      legend.position = "bottom",
      legend.text     = element_text(size = 10),
      title           = element_text(size = 11),
      panel.grid      = element_blank(),
      axis.title      = element_blank(),
      axis.text.x     = element_text(size = 13),
      strip.text.x    = element_text(size = 16),
      strip.text.y    = element_text(size = 10),
      panel.spacing   = unit(0.1, "lines"),
      axis.ticks      = element_blank()
    )
})
names(flux_bar_plots) <- scfa_names

flux_combined_bars <- cowplot::plot_grid(
  flux_bar_plots$Acetate,
  flux_bar_plots$Propionate,
  flux_bar_plots$Butyrate + scale_y_continuous(breaks = c(0, 2e+05, 4e+05, 6e+05, 8e+05)),
  ncol = 1
)

pdf(paste0(path_outfiles,
           "Plots/SCFAFluxBySpecie_BarPlots_coordFlip_separateScales_corn20240719_wheat20240719.pdf"),
    height = 40, width = 10)
flux_combined_bars
dev.off()


# =========================================================================
#  FIGURE S2: EC number heatmap annotated by KEGG/MetaCyc pathways
# =========================================================================
models_compartments <- readxl::read_xlsx(
  "ModelGeneration_Files/SpecieModels_present_in_compartments_D10_noAGP_2022upd.xlsx",
  sheet = "Gapseq_ModelList_upd"
)

ec_by_model <- read.delim(
  "ModelsAnalyses_outputs/Table_ECnumbers_inEachModel_gapseqAdapted_cobrapyAdapted_July2023.txt",
  sep = '\t', header = TRUE
)
ec_by_model$tax_family <- sapply(ec_by_model$SpecieName, function(x) {
  unique(models_compartments$TaxonomicFamily_2024[models_compartments$Specie == x])
})

family_colors <- c(
  Lachnospiraceae       = "#2A114A",
  Clostridiaceae        = colpal3[12],
  Ruminococcaceae       = colpal3[2],
  Oscillospiraceae      = "#bd7eb2",
  Peptostreptococcaceae = colpal3[9],
  Borkfalkiaceae        = "#FFD4EB",
  Lactobacillaceae      = colpal3[5],
  Streptococcaceae      = colpal3[19],
  Staphylococcaceae     = colpal3[11],
  Enterococcaceae       = "#cc9616",
  Erysipelotrichaceae   = "#BA5D00",
  Veillonellaceae       = colpal3[15],
  Enterobacteriaceae    = colpal3[10],
  Sphingomonadaceae     = colpal3[4],
  Pseudomonadaceae      = colpal3[6],
  Moraxellaceae         = "#dec39b",
  Bacteroidaceae        = colpal3[7],
  Bifidobacteriaceae    = "#c2b8b2",
  Corynebacteriaceae    = "#267838"
)

# Expand ECs into separate rows, then build binary presence/absence matrix
ec_expanded <- ec_by_model %>%
  tidyr::separate_rows(EC_numbers, sep = ", ") %>%
  mutate(EC_numbers = trimws(EC_numbers))

ec_binary <- ec_expanded %>%
  tidyr::pivot_wider(names_from = EC_numbers, values_from = EC_numbers,
                     values_fill = list(EC_numbers = 0),
                     values_fn = list(EC_numbers = length)) %>%
  mutate(across(-c(SpecieName, tax_family), ~ ifelse(. > 0, 1, 0))) %>%
  as.data.frame()

rownames(ec_binary) <- ec_binary$SpecieName

# Load pathway annotation tables
kegg_pwy    <- read.delim("kegg_pwy_info.tbl")
metacyc_pwy <- read.delim("meta_pwy.tbl.txt")

ec_columns <- colnames(ec_binary)[-c(1:3)]

# Function to extract pathway info for EC numbers
extract_ec_info <- function(pwy_data, ec_numbers) {
  pwy_data %>%
    tidyr::separate_rows(reaEc, sep = ",") %>%
    mutate(reaEc = trimws(reaEc)) %>%
    filter(reaEc %in% ec_numbers) %>%
    group_by(reaEc) %>%
    summarise(name = list(name), hierarchy = list(hierarchy), .groups = "drop")
}

ec_info_kegg <- extract_ec_info(kegg_pwy, ec_columns)
ec_info_meta <- extract_ec_info(metacyc_pwy, ec_columns)

# Unnest and deduplicate
ec_kegg_flat <- ec_info_kegg %>%
  tidyr::unnest(cols = c(name, hierarchy)) %>%
  distinct(reaEc, hierarchy, .keep_all = TRUE)

ec_meta_flat <- ec_info_meta %>%
  tidyr::unnest(cols = c(name, hierarchy)) %>%
  distinct(reaEc, hierarchy, .keep_all = TRUE)

# Keep only MetaCyc ECs not already in KEGG
meta_only_ecs <- setdiff(unique(ec_info_meta$reaEc), unique(ec_info_kegg$reaEc))
ec_meta_flat  <- ec_meta_flat %>% filter(reaEc %in% meta_only_ecs)

# Map MetaCyc hierarchy strings to standardized pathway names
metacyc_pathway_map <- c(
  "Lipid-biosynthesis|Lipid-Biosynthesis|Lipid-Degradation" = "Lipid metabolism",
  "Amino-Acid-Degradation"    = "Amino acid metabolism",
  "Amino-Acid-Biosynthesis"   = "Amino acid metabolism",
  "Glycan"                    = "Glycan biosynthesis and metabolism",
  "Energy-Metabolism"         = "Energy metabolism",
  "Carbohydrates-Degradation" = "Carbohydrate metabolism",
  "Carbohydrates-Biosynthesis" = "Carbohydrate metabolism",
  "CARBOXYLATES-DEG"          = "Carbohydrate metabolism",
  "Cofactor"                  = "Metabolism of cofactors and vitamins",
  "AROMATIC"                  = "Xenobiotics biodegradation and metabolism",
  "CO2-Fixation"              = "Energy metabolism",
  "Polymer-Degradation"       = "Carbohydrate metabolism"
)

ec_meta_flat$Pathway <- NA
for (pattern in names(metacyc_pathway_map)) {
  matches <- grep(pattern, ec_meta_flat$hierarchy)
  ec_meta_flat$Pathway[matches] <- metacyc_pathway_map[pattern]
}

ec_meta_flat <- ec_meta_flat %>%
  filter(!is.na(Pathway)) %>%
  distinct(reaEc, Pathway, .keep_all = TRUE)

# Standardize KEGG pathway names
ec_kegg_flat$Pathway <- gsub(
  "kegg;Metabolism;|kegg;Genetic Information Processing;|kegg;Environmental Information Processing;",
  "", ec_kegg_flat$hierarchy
)
ec_kegg_flat <- ec_kegg_flat %>% filter(Pathway != "Signal transduction")

# Combine KEGG + MetaCyc
ec_pwy_combined <- bind_rows(ec_kegg_flat, ec_meta_flat)

# Handle ECs with multiple pathway annotations
multi_annot <- ec_pwy_combined %>%
  group_by(reaEc) %>% filter(n() > 1) %>% ungroup()

multi_annot_filtered <- multi_annot %>%
  filter(!(reaEc %in% reaEc[hierarchy == "kegg;Metabolism;Global and overview maps"] &
             hierarchy == "kegg;Metabolism;Global and overview maps"))

ec_pwy_clean <- ec_pwy_combined %>%
  filter(!reaEc %in% multi_annot$reaEc) %>%
  bind_rows(multi_annot_filtered)

# Assign Pathway_category1 and Pathway_category2 per EC
assign_pathway_categories <- function(pathways) {
  uq <- unique(pathways)
  if (length(uq) == 1) {
    return(data.frame(Pathway_category1 = uq, Pathway_category2 = NA))
  }
  if (length(uq) >= 3) {
    if ("Metabolism of other amino acids" %in% uq & "Amino acid metabolism" %in% uq)
      uq <- uq[uq != "Metabolism of other amino acids"]
    if ("Carbohydrate metabolism" %in% uq & "Energy metabolism" %in% uq)
      uq <- uq[uq != "Energy metabolism"]
  }
  if (length(uq) > 3) {
    return(data.frame(Pathway_category1 = "Multiple superpathways", Pathway_category2 = NA))
  }
  sorted <- sort(uq)
  data.frame(Pathway_category1 = sorted[1],
             Pathway_category2 = if (length(sorted) >= 2) sorted[2] else NA)
}

ec_pwy_categorized <- ec_pwy_clean %>%
  group_by(reaEc) %>%
  do(assign_pathway_categories(.$Pathway)) %>%
  right_join(ec_pwy_clean, by = "reaEc") %>%
  ungroup() %>%
  distinct(reaEc, Pathway_category1, Pathway_category2)

pathway_order <- c(
  "Multiple superpathways", "Carbohydrate metabolism", "Energy metabolism",
  "Amino acid metabolism", "Metabolism of other amino acids", "Lipid metabolism",
  "Metabolism of terpenoids and polyketides", "Metabolism of cofactors and vitamins",
  "Glycan biosynthesis and metabolism", "Xenobiotics biodegradation and metabolism",
  "Nucleotide metabolism", "Biosynthesis of other secondary metabolites",
  "Global and overview maps", "Translation", "Unknown"
)

# Ensure all ECs in the binary matrix are represented
all_ecs <- colnames(ec_binary)[-c(1:2)]
ec_annot <- tibble(reaEc = all_ecs) %>%
  left_join(ec_pwy_categorized, by = "reaEc") %>%
  mutate(
    Pathway_category1 = factor(ifelse(is.na(Pathway_category1), "Unknown", Pathway_category1),
                                levels = pathway_order),
    Pathway_category2 = ifelse(is.na(Pathway_category2), "Unknown", Pathway_category2)
  ) %>%
  arrange(Pathway_category1, match(reaEc, all_ecs))

# Reorder binary matrix columns by pathway
ec_binary <- ec_binary %>%
  select(tax_family, SpecieName, all_of(ec_annot$reaEc))

# Color palette for pathway categories
all_categories <- unique(c(
  as.character(ec_annot$Pathway_category1[!is.na(ec_annot$Pathway_category1)]),
  ec_annot$Pathway_category2
))
all_categories <- all_categories[!is.na(all_categories)]

if (length(all_categories) > 11) {
  pathway_colors <- setNames(
    colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(length(all_categories)),
    all_categories
  )
}
pathway_colors["Multiple superpathways"] <- "black"
pathway_colors["Unknown"] <- "white"

col_annotation <- HeatmapAnnotation(
  Pathway_category1 = ec_annot$Pathway_category1,
  Pathway_category2 = ec_annot$Pathway_category2,
  col = list(Pathway_category1 = pathway_colors, Pathway_category2 = pathway_colors),
  annotation_legend_param = list(
    Pathway_category1 = list(title = "Pathway Category 1", at = pathway_order, labels = pathway_order),
    Pathway_category2 = list(title = "Pathway Category 2", at = pathway_order, labels = pathway_order)
  ),
  show_annotation_name = FALSE
)

row_annotation <- rowAnnotation(
  tax_family = ec_binary$tax_family,
  col = list(tax_family = family_colors),
  show_annotation_name = FALSE
)

ec_heatmap_data <- ec_binary %>% select(-tax_family, -SpecieName) %>% as.matrix()

ec_hm_colors <- colorRamp2(c(0, 1), c("#0f075e", "#ed1561"))

ec_heatmap <- Heatmap(
  ec_heatmap_data, name = "ECs",
  cluster_columns = FALSE, show_column_dend = FALSE,
  top_annotation = col_annotation, right_annotation = row_annotation,
  show_row_names = TRUE, show_column_names = FALSE,
  show_heatmap_legend = TRUE, col = ec_hm_colors
)

ec_heatmap_clustered <- Heatmap(
  ec_heatmap_data, name = "ECs",
  cluster_columns = TRUE, cluster_rows = TRUE,
  show_column_dend = FALSE,
  top_annotation = col_annotation, right_annotation = row_annotation,
  show_row_names = TRUE, show_column_names = FALSE,
  show_heatmap_legend = TRUE, col = ec_hm_colors
)

pdf("ModelsAnalyses_outputs/HeatMap_allECs_KEGGMetaCyc_annotated_30MAGs_v280824.pdf",
    width = 30, height = 8)
draw(ec_heatmap)
dev.off()

# Export heatmap order as table
ht_drawn <- draw(ec_heatmap)
ro <- unlist(row_order(ht_drawn))
co <- unlist(column_order(ht_drawn))
mat_ordered <- ec_heatmap_data[ro, co]
mat_with_annot <- rbind(
  as.character(ec_annot$Pathway_category1),
  gsub("Unknown", "", as.character(ec_annot$Pathway_category2)),
  mat_ordered
)
write.csv(mat_with_annot,
          "ModelsAnalyses_outputs/Table_HeatMap_allECs_KEGGMetaCyc_annotated_30MAGs_v280824.csv")

# Clustered version
ht_clus_drawn <- draw(ec_heatmap_clustered)
ro_clus <- unlist(row_order(ht_clus_drawn))
co_clus <- unlist(column_order(ht_clus_drawn))
mat_clus <- ec_heatmap_data[ro_clus, co_clus]
mat_clus_annot <- rbind(
  as.character(ec_annot$Pathway_category1),
  gsub("Unknown", "", as.character(ec_annot$Pathway_category2)),
  mat_clus
)
write.csv(mat_clus_annot,
          "ModelsAnalyses_outputs/Table_HeatMap_Clustered_allECs_KEGGMetaCyc_annotated_30MAGs_v051225.csv")


