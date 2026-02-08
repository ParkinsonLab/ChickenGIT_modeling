# PCoA analysis comparing microbial family-level compositions across three datasets:
# initial 16S (D10 corn), 96h simulation relative abundances, and 16S from the validation trial
# for ileum and cecum compartments.

library(phyloseq)   # for microbiome data structures and analyses
library(vegan)      # for PERMANOVA (adonis/adonis2) and ordination
library(ggplot2)    # for plotting
library(dplyr)      # for data manipulation
library(ALDEx2)
library(reshape2)
library(tidyr)
library(stringr)
library(ANCOMBC)
library(ggpubr)
library(readxl)
library(BacArena)
library(tibble)

# Define custom colors for treatments
treatment_colors <- c(
  "Control" = "#2C3E50",
  "Cellulose Low (2g)" = "#E74C3C",
  "Cellulose High (8g/100g)" = "darkred",
  "L-Threonine Low (1g/100g)" = "#3498DB",
  "L-Threonine High (5g/100g)" = "#2f2fad",
  "Starch Low (5g/100g)" = "#27AE60",
  "Starch High (10g/100g)" = "#216907"
)

# Function to prepare dataset 1: Post-trial 2025 data (validation)
prepare_validation_data <- function(validation_rel_abundance) {
  # This uses the relative abundance data from your validation trials
  validation_data <- validation_rel_abundance %>%
    filter(site %in% c("Ileum", "Cecum")) %>%
    group_by(sample_id, site, group, family) %>%
    summarise(rel_abundance = sum(rel_abundance, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      data_source = "Validation_Trial",
      treatment = group,
      compartment = site
    ) %>%
    dplyr::select(sample_id, compartment, treatment, family, rel_abundance, data_source)

  return(validation_data)
}

# Function to prepare dataset 2: Initial D10 data (modeling source)
prepare_initial_data <- function(initial_rel_abundance) {
  # Filter for Ileum and Cecum from the initial 6-compartment data
  initial_data <- initial_rel_abundance %>%
    filter(description %in% c("Ileum", "Cecum")) %>%
    group_by(sample_id, description, family) %>%
    summarise(rel_abundance = sum(rel_abundance, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      data_source = "Initial_Data",
      treatment = "Control",  # All initial data treated as baseline/control
      compartment = description
    ) %>%
    dplyr::select(sample_id, compartment, treatment, family, rel_abundance, data_source)

  return(initial_data)
}

# Function to prepare dataset 3: Simulation predictions (from loaded data)
prepare_simulation_data_from_loaded <- function(simulation_rel_abundance) {
  simulation_data <- simulation_rel_abundance %>%
    mutate(
      data_source = "Simulation",
      compartment = site  # Ensure consistent naming
    ) %>%
    dplyr::select(sample_id, compartment, treatment, family, rel_abundance, data_source)

  return(simulation_data)
}

# Function to combine all datasets
combine_datasets <- function(validation_data, initial_data, simulation_data) {
  # Ensure all datasets have the same structure and create unique sample IDs
  # that include compartment information to avoid duplicates
  combined_data <- rbind(validation_data, initial_data, simulation_data) %>%
    # Create unique sample IDs by combining sample_id and compartment
    mutate(
      unique_sample_id = paste(sample_id, compartment, sep = "_"),
      # Handle missing families and standardize names
      family = ifelse(is.na(family) | family == "", "Unknown", family),
      family = str_remove(family, "^f__"),  # Remove QIIME2 prefixes if present
      # Standardize treatment names
      treatment = case_when(
        str_detect(treatment, "Cellulose.*2g|Low.*Cellulose") ~ "Cellulose Low (2g)",
        str_detect(treatment, "Cellulose.*8g|High.*Cellulose") ~ "Cellulose High (8g/100g)",
        str_detect(treatment, "Threonine.*1g|Low.*Threonine") ~ "L-Threonine Low (1g/100g)",
        str_detect(treatment, "Threonine.*5g|High.*Threonine") ~ "L-Threonine High (5g/100g)",
        str_detect(treatment, "Starch.*5g|Low.*Starch") ~ "Starch Low (5g/100g)",
        str_detect(treatment, "Starch.*10g|High.*Starch") ~ "Starch High (10g/100g)",
        TRUE ~ "Control"
      )
    )

  return(combined_data)
}

# Function to create wide format matrix for ordination
create_abundance_matrix <- function(combined_data) {
  # First, ensure no duplicates by summing any duplicate unique_sample_id + family combinations
  abundance_matrix <- combined_data %>%
    dplyr::select(unique_sample_id, family, rel_abundance) %>%
    group_by(unique_sample_id, family) %>%
    summarise(rel_abundance = sum(rel_abundance, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = family, values_from = rel_abundance, values_fill = 0) %>%
    column_to_rownames("unique_sample_id")

  # Create metadata using unique sample IDs
  metadata <- combined_data %>%
    dplyr::select(unique_sample_id, compartment, treatment, data_source) %>%
    distinct() %>%
    column_to_rownames("unique_sample_id")

  return(list(abundance = abundance_matrix, metadata = metadata))
}

# Function to perform PCoA analysis
perform_pcoa <- function(abundance_matrix, metadata, compartment_filter) {
  # Filter for specific compartment
  comp_metadata <- metadata[metadata$compartment == compartment_filter, ]
  comp_abundance <- abundance_matrix[rownames(comp_metadata), ]

  # Remove families with zero abundance across all samples
  comp_abundance <- comp_abundance[, colSums(comp_abundance) > 0]

  # Calculate Bray-Curtis distance
  bray_dist <- vegdist(comp_abundance, method = "bray")

  # Perform PCoA
  pcoa_result <- cmdscale(bray_dist, eig = TRUE, k = 2)

  # Calculate percentage of variance explained
  eigenvals <- pcoa_result$eig[pcoa_result$eig > 0]
  var_explained <- eigenvals / sum(eigenvals) * 100

  # Create PCoA dataframe
  pcoa_df <- data.frame(
    unique_sample_id = rownames(pcoa_result$points),
    PC1 = pcoa_result$points[, 1],
    PC2 = pcoa_result$points[, 2]
  ) %>%
    merge(comp_metadata, by.x = "unique_sample_id", by.y = "row.names")

  return(list(
    pcoa_df = pcoa_df,
    var_explained = var_explained[1:2],
    distance_matrix = bray_dist
  ))
}

# Function to create PCoA plot with confidence ellipses
create_pcoa_plot <- function(pcoa_result, compartment_name) {
  pcoa_df <- pcoa_result$pcoa_df
  var_explained <- pcoa_result$var_explained

  # Create base plot
  p <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = treatment, shape = data_source)) +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(aes(group = interaction(treatment, data_source)),
                 level = 0.68, type = "norm", alpha = 0.4) +
    scale_color_manual(values = treatment_colors) +
    scale_shape_manual(values = c("Validation_Trial" = 16,
                                  "Initial_Data" = 17,
                                  "Simulation" = 15)) +
    labs(
      title = paste("PCoA Analysis -", compartment_name),
      subtitle = "Comparison of microbial compositions (family level) across data sources",
      x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
      y = paste0("PC2 (", round(var_explained[2], 1), "%)"),
      color = "Treatment",
      shape = "Data Source"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 15, face = "bold"),
      plot.subtitle = element_text(size = 14),
      axis.text = element_text(size=13),
      axis.title = element_text(size=15),
      legend.position = "bottom",
      legend.box = "horizontal",
      panel.grid.minor = element_blank()
    ) +
    guides(
      color = guide_legend(override.aes = list(size = 4)),
      shape = guide_legend(override.aes = list(size = 4))
    )

  return(p)
}

# Function to create faceted plot by treatment
create_faceted_pcoa_plot <- function(pcoa_result, compartment_name) {
  pcoa_df <- pcoa_result$pcoa_df
  var_explained <- pcoa_result$var_explained

  p <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = data_source, shape = data_source)) +
    geom_point(size = 3, alpha = 0.7) +
    stat_ellipse(aes(group = data_source), level = 0.68, type = "norm", alpha = 0.3) +
    scale_color_manual(values = c("Validation_Trial" = "#E74C3C",
                                  "Initial_Data" = "#3498DB",
                                  "Simulation" = "#27AE60")) +
    scale_shape_manual(values = c("Validation_Trial" = 16,
                                  "Initial_Data" = 17,
                                  "Simulation" = 15)) +
    facet_wrap(~ treatment, scales = "free") +
    labs(
      title = paste("PCoA Analysis by Treatment -", compartment_name),
      x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
      y = paste0("PC2 (", round(var_explained[2], 1), "%)"),
      color = "Data Source",
      shape = "Data Source"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "lightgray")
    )

  return(p)
}

# Function to perform PERMANOVA test
perform_permanova <- function(pcoa_result) {
  # Extract distance matrix and metadata
  dist_matrix <- pcoa_result$distance_matrix
  metadata <- pcoa_result$pcoa_df

  # PERMANOVA test
  permanova_result <- adonis2(dist_matrix ~ treatment * data_source,
                              data = metadata, permutations = 999)

  return(permanova_result)
}


## =====================================================================
##  LOAD 16S data post-trial (ileum and cecum)
## =====================================================================

load_validation_data <- function() {
  setwd("20241217_Parkinson_Results_16S_chicken_133samples")

  # Load ASV table
  asv_table <- read.table("qiime_2023.2_reprocessing/asv_table_w_tax_112samples_minDepth500.txt", sep='\t', header=T)
  colnames(asv_table) <- gsub("X","",colnames(asv_table))

  # Load taxonomy
  taxonomy <- read.csv("qiime_2023.2_reprocessing/taxonomy.tsv", sep='\t')
  taxonomy_split <- separate(data = taxonomy,
                             col = taxonomy,
                             into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                             sep = "\\; ")
  taxonomy_split <- taxonomy_split[,-c(1,9)]
  rownames(taxonomy_split) <- taxonomy[,1]

  # Load sample and treatment info
  sample_info <- read_excel("20241217_Parkinson_PrepLog.xlsx", sheet = "Parkinson subm", col_names = F)
  treatment_info <- read_excel("../../SampleData_chicken_trials_2024/TreatmentGroup_labeling.xlsx")

  # Process sample metadata
  treatment_metadata <- treatment_info %>%
    mutate(group = gsub("G\\d+\\- ","", `Treatment Groups`)) %>%
    mutate(group_short = str_extract(`Treatment Groups`, "G\\d+"))

  sample_metadata <- sample_info %>%
    dplyr::select(1:3) %>%
    rename(sample_id = 1, sample_name = 2, description = 3) %>%
    mutate(
      group_short = str_extract(sample_name, "G\\d+"),
      site = case_when(
        str_detect(sample_name, "_I") ~ "Ileum",
        str_detect(sample_name, "_T") ~ "Cecum",
        TRUE ~ NA_character_
      )
    ) %>%
    left_join(treatment_metadata, by = "group_short") %>%
    mutate(sample_id = as.character(sample_id))

  # Create abundance data
  asv_table_num <- data.matrix(asv_table[,2:113])
  rownames(asv_table_num) <- asv_table[,1]

  abundance_data <- asv_table %>%
    dplyr::select(-taxonomy) %>%
    pivot_longer(cols = -c(1), names_to = "sample_id", values_to = "count") %>%
    rename(OTU_id = 1)

  # Calculate relative abundance with family-level aggregation
  relative_abundance <- abundance_data %>%
    group_by(sample_id) %>%
    mutate(rel_abundance = count / sum(count)) %>%
    ungroup() %>%
    left_join(data.frame(OTU_id = rownames(taxonomy_split), taxonomy_split), by = "OTU_id") %>%
    left_join(sample_metadata, by = "sample_id") %>%
    filter(!is.na(site)) %>%  # Only keep Ileum and Cecum samples
    group_by(sample_id, site, group, Family) %>%
    summarise(rel_abundance = sum(rel_abundance, na.rm = TRUE), .groups = "drop") %>%
    rename(family = Family) %>%
    mutate(family = ifelse(is.na(family) | family == "", "Unknown", family))

  return(relative_abundance)
}

validation_data_rel_abundance <- load_validation_data()


## =====================================================================
##  LOAD 16S initial data (6 compartments) used for modeling
## =====================================================================

load_initial_data <- function() {
  setwd("Chicken_16S_Angela_D10_noAGP_qiime2_rerun_2025/")

  # Read ASV table
  asv_table <- read.csv("exports/feature-table_w_tax.txt", sep='\t', skip=1, header=T)
  colnames(asv_table)[1] <- "OTU_id"
  colnames(asv_table) <- gsub("\\.","-",colnames(asv_table))

  # Read sample info
  sample_info <- read.table("../ChickenGIT_modeling_project/ModelGeneration_Files/metadata.txt", header = T)

  # Process taxonomy
  taxonomy <- asv_table %>%
    dplyr::select(OTU_id, taxonomy) %>%
    separate(taxonomy,
             into = c("domain", "phylum", "class", "order", "family", "genus", "species"),
             sep = "; ",
             fill = "right") %>%
    mutate(across(everything(), ~str_remove(., "^[a-z]__")))

  # Convert to long format
  abundance_data <- asv_table %>%
    dplyr::select(-taxonomy) %>%
    pivot_longer(cols = -OTU_id, names_to = "sample_id", values_to = "count")

  # Process sample metadata
  sample_metadata <- sample_info %>%
    dplyr::select(1,3,6,9) %>%
    rename(sample_id = 2, sample_name = 1, diet = 3, description = 4)

  # Calculate relative abundance with family-level aggregation
  relative_abundance <- abundance_data %>%
    group_by(sample_id) %>%
    mutate(rel_abundance = count / sum(count)) %>%
    ungroup() %>%
    left_join(taxonomy, by = "OTU_id") %>%
    left_join(sample_metadata, by = "sample_id") %>%
    filter(description %in% c("Ileum", "Ceca")) %>%  # Only keep Ileum and Cecum
    group_by(sample_id, description, family) %>%
    summarise(rel_abundance = sum(rel_abundance, na.rm = TRUE), .groups = "drop") %>%
    mutate(family = ifelse(is.na(family) | family == "", "Unknown", family))

  return(relative_abundance)
}

initial_data_rel_abundance <- load_initial_data()

# keep only corn diet:
initial_data_rel_abundance_corn <- initial_data_rel_abundance[
  grep("D10-5-1|D10-5-2|D10-5-3|D10-10-1|D10-10-2", initial_data_rel_abundance$sample_id),
]
initial_data_rel_abundance_corn$description <- gsub("Ceca","Cecum",
                                                    initial_data_rel_abundance_corn$description)


## =====================================================================
##  LOAD relative abundances resulting from 96h simulations with/without compounds
## =====================================================================
load_simulation_data <- function() {
  setwd("WithPrebiotics/2024/")
  path <- "Tables/"
  pattern <- "corn_diet_grower_NEW_AFTER48h_96h_econc1_propGizzard180_CobaltX12_rmAnaero_less02_mucUrea_wCoA_moreNightFlow_fix_allCompounds_20240721"

  # Read the simulation relative abundance table
  relab_file <- paste0(path, "RelAbund_tables_96h_", pattern, "_meanByRep.txt")

  if (!file.exists(relab_file)) {
    stop(paste("Simulation file not found:", relab_file))
  }

  relab_table <- read.table(relab_file, sep='\t', header=T)

  # Filter for only Ileum and Cecum, and final timepoint (96h)
  relative_abundance <- relab_table %>%
    filter(Compartment %in% c("Ileum", "Cecum"),
           Hour == 96) %>%  # Use final timepoint
    # Clean species names and map to families
    mutate(
      Species_clean = gsub("-gapseqAdRm|_cobrapy_adapted", "", Specie),
      family = case_when(
        str_detect(Species_clean, "Blautia|Mediterraneibacter|Lachnoclostridium") ~ "Lachnospiraceae",
        str_detect(Species_clean, "Clostridium_saudiense") ~ "Clostridiaceae",
        str_detect(Species_clean, "Faecalibacterium|Agathobaculum|Flavonifractor|Oscillospiraceae") ~ "Ruminococcaceae",
        str_detect(Species_clean, "Romboutsia") ~ "Peptostreptococcaceae",
        str_detect(Species_clean, "Borkfalkia") ~ "Lachnospiraceae",
        str_detect(Species_clean, "Lactobacillus|Limosilactobacillus") ~ "Lactobacillaceae",
        str_detect(Species_clean, "Streptococcus") ~ "Streptococcaceae",
        str_detect(Species_clean, "Staphylococcus|Enterococcus") ~ "Enterococcaceae",
        str_detect(Species_clean, "Erysipelatoclostridium") ~ "Erysipelatoclostridiaceae",
        str_detect(Species_clean, "Massiliomicrobiota") ~ "Erysipelotrichaceae",
        str_detect(Species_clean, "Megasphaera") ~ "Veillonellaceae",
        str_detect(Species_clean, "Escherichia|Klebsiella") ~ "Enterobacteriaceae",
        str_detect(Species_clean, "Sphingomonas") ~ "Sphingomonadaceae",
        str_detect(Species_clean, "Acinetobacter") ~ "Moraxellaceae",
        str_detect(Species_clean, "Pseudomonas") ~ "Pseudomonadaceae",
        str_detect(Species_clean, "Bifidobacterium") ~ "Bifidobacteriaceae",
        str_detect(Species_clean, "Corynebacterium") ~ "Corynebacteriaceae",
        str_detect(Species_clean, "Phocaeicola") ~ "Bacteroidaceae",
        TRUE ~ "Other"
      )
    ) %>%
    # Aggregate by family and create consistent structure
    group_by(Sample, compound, compound_gper100g, Compartment, medium, family) %>%
    summarise(rel_abundance = sum(RelAbundance, na.rm = TRUE), .groups = "drop") %>%
    # Create consistent sample IDs and treatment names
    mutate(
      sample_id = paste0("sim_", Sample, "_", compound, "_", compound_gper100g),
      treatment = case_when(
        compound == "None" ~ "Control",
        compound == "Cellulose" & compound_gper100g == "2" ~ "Cellulose Low (2g)",
        compound == "Cellulose" & compound_gper100g == "8" ~ "Cellulose High (8g/100g)",
        compound == "L-Threonine" & compound_gper100g == "1" ~ "L-Threonine Low (1g/100g)",
        compound == "L-Threonine" & compound_gper100g == "5" ~ "L-Threonine High (5g/100g)",
        compound == "Starch_n19" & compound_gper100g == "5" ~ "Starch Low (5g/100g)",
        compound == "Starch_n19" & compound_gper100g == "10" ~ "Starch High (10g/100g)",
        TRUE ~ paste(compound, compound_gper100g)
      ),
      site = Compartment  # Rename to match validation data structure
    ) %>%
    # Final family-level aggregation and cleanup
    group_by(sample_id, site, treatment, family) %>%
    summarise(rel_abundance = sum(rel_abundance, na.rm = TRUE), .groups = "drop") %>%
    mutate(family = ifelse(is.na(family) | family == "", "Unknown", family))

  return(relative_abundance)
}

simulation_rel_abundance <- load_simulation_data()

# subset to 3 treatments only:
simulation_rel_abundance <-
  simulation_rel_abundance[grep(
    "Threonine_1|Threonine_5|Cellulose_2|Cellulose_8|Starch_n19_5|Starch_n19_10", simulation_rel_abundance$sample_id),]


## =====================================================
## Combine datasets and run PCoA - make plots
## =====================================================

# Main analysis function for three loaded datasets
run_pcoa_analysis_three_datasets <- function(validation_rel_abundance, initial_rel_abundance, simulation_rel_abundance) {

  # Prepare datasets
  cat("Preparing validation data...\n")
  validation_data <- prepare_validation_data(validation_rel_abundance)

  cat("Preparing initial data...\n")
  initial_data <- prepare_initial_data(initial_rel_abundance)

  cat("Preparing simulation data...\n")
  simulation_data <- prepare_simulation_data_from_loaded(simulation_rel_abundance)

  cat("Combining datasets...\n")
  combined_data <- combine_datasets(validation_data, initial_data, simulation_data)

  cat("Creating abundance matrix...\n")
  matrix_data <- create_abundance_matrix(combined_data)

  # Analyze only Ileum and Cecum
  compartments <- c("Ileum", "Cecum")
  results <- list()

  for (comp in compartments) {
    cat(paste("Analyzing", comp, "...\n"))

    # Perform PCoA
    pcoa_result <- perform_pcoa(matrix_data$abundance, matrix_data$metadata, comp)

    # Create plots
    plot_combined <- create_pcoa_plot(pcoa_result, comp)
    plot_faceted <- create_faceted_pcoa_plot(pcoa_result, comp)

    # Perform PERMANOVA
    permanova_result <- perform_permanova(pcoa_result)

    # Store results
    results[[comp]] <- list(
      pcoa_result = pcoa_result,
      plot_combined = plot_combined,
      plot_faceted = plot_faceted,
      permanova = permanova_result
    )
  }

  return(results)
}

results <- run_pcoa_analysis_three_datasets(validation_rel_abundance = validation_data_rel_abundance,
                                            simulation_rel_abundance = simulation_rel_abundance,
                                            initial_rel_abundance = initial_data_rel_abundance_corn)

# Function to save plots
save_pcoa_plots <- function(results,
                            output_dir = "PCoA_plots_3datasets/") {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  for (comp in names(results)) {
    # Save combined plot
    ggsave(
      filename = paste0(output_dir, "PCoA_", comp, "_combined.pdf"),
      plot = results[[comp]]$plot_combined,
      width = 10, height = 7,
    )

    # Save faceted plot
    ggsave(
      filename = paste0(output_dir, "PCoA_", comp, "_faceted_by_treatment.pdf"),
      plot = results[[comp]]$plot_faceted,
      width = 15, height = 10, units = "in"
    )

    # Save PERMANOVA results
    write.table(results[[comp]]$permanova,
                file = paste0(output_dir, "PERMANOVA_", comp, ".txt"),
                sep = "\t", quote = FALSE)
  }

}

save_pcoa_plots(results)


## =============================================================================
## compare relative abundances of certain taxa in 16S datasets
## =============================================================================
# load original 16S (6-comp) dataset
initial_data_rel_abundance <- load_initial_data()
# keep only corn diet:
initial_data_rel_abundance_corn <- initial_data_rel_abundance[
  grep("D10-5-1|D10-5-2|D10-5-3|D10-10-1|D10-10-2", initial_data_rel_abundance$sample_id),
]

# load validation 16S
validation_data_rel_abundance <- load_validation_data()

mean_relabund_families_initial <- initial_data_rel_abundance_corn %>%
  dplyr::filter(description == "Ceca") %>%
  group_by(family) %>%
  mutate(mean = mean(rel_abundance), sd = sd(rel_abundance)) %>%
  dplyr::select(-c(sample_id,rel_abundance)) %>%
  distinct()

mean_relabund_families_trial <- validation_data_rel_abundance %>%
  dplyr::filter(site == "Cecum") %>%
  group_by(family) %>%
  mutate(mean = mean(rel_abundance), sd = sd(rel_abundance)) %>%
  dplyr::select(-c(sample_id,rel_abundance,group)) %>%
  distinct()

mean_relabund_families_trial$family <- gsub("f__","",mean_relabund_families_trial$family)

families_of_interest = c("Enterobacteriaceae", "Lactobacillaceae", "Ruminococcaceae",
                         "Streptococcaceae")
mean_relabund_families_initial %>% filter(family %in% families_of_interest)
# description family                mean      sd
# 1 Ceca        Enterobacteriaceae 0.213   0.0784
# 2 Ceca        Lactobacillaceae   0.0530  0.0389
# 3 Ceca        Ruminococcaceae    0.0423  0.0302
# 4 Ceca        Streptococcaceae   0.00587 0.00422
mean_relabund_families_trial %>% filter(family %in% families_of_interest)
# site  family                mean      sd
# 1 Cecum Enterobacteriaceae 0.00850 0.0175
# 2 Cecum Lactobacillaceae   0.00514 0.0188
# 3 Cecum Ruminococcaceae    0.294   0.0835
# 4 Cecum Streptococcaceae   0.00217 0.00204
