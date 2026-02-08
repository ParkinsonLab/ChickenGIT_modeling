# Paired Wilcoxon signed-rank tests on fold changes (compound vs control) for
# total metabolite fluxes and 96h concentrations across 6 GIT compartments.
# Includes BH correction and a by-diet subsampled analysis (10 reps per diet).

library(dplyr)
library(magrittr)
library(tidyr)
library(purrr)

myargs <- commandArgs(trailingOnly = TRUE)
print(myargs)

# ===== Loading input parameters =============================

pattern_control=as.character(myargs[1])

pattern_compounds=as.character(myargs[2])

# in the format corn20240707_wheat20240708
sim_dates_control=as.character(myargs[3])

# in the format corn20240709_wheat20240709
sim_dates_compounds=as.character(myargs[4])

# ===== Paths =============================

base_dir <- "Chicken_10D_simulation/"
path <- file.path(base_dir, "arena_sim_objects/CoupledSixComps_D10/WithPrebiotics/2024/gapseq_models/OUT_files/Tables/")
path_control <- file.path(base_dir, "arena_sim_objects/CoupledSixComps_D10/NoPathogens/2024/gapseq_models/OUT_files/Tables/")

seed_metabolites <- data.table::fread(file.path(base_dir, "all_seed_metabolites_edited.tsv"), sep = '\t', data.table=F)
seed_metabolites <- as.data.frame(seed_metabolites)

# ----------------------------------------------------------------------------------------
##########################################################################################
################        FLUXES (PRODUCTION) WILCOX TEST             ######################
##########################################################################################
# ----------------------------------------------------------------------------------------

fluxes_production_control <- read.table(paste0(path_control,"flux_table_96h_sum_allMets_ByCommunity_control_",
                                               pattern_control,
                                               "_", sim_dates_control, ".txt"),
                                        sep='\t',header=T)

fluxes_production_control$compound_amount <- "control"


flux_total_sum_compounds <- read.table(paste0(path,"flux_table_96h_sum_allMets_ByCommunity_",
                                              pattern_compounds,
                                              "_", sim_dates_compounds, ".txt"),
                                       sep='\t',header=T)


flux_total_sum_compounds$compound_amount <- paste0(flux_total_sum_compounds$compound, "_", flux_total_sum_compounds$compound_gper100g)
flux_total_sum_compounds$compound <- NULL
flux_total_sum_compounds$compound_gper100g <- NULL


# leave only non-negative fluxes (metabolite production):
fluxes_total_sum_compounds_production <- flux_total_sum_compounds %>%
  filter(Flux >= 0)

# Remove biomass
fluxes_total_sum_compounds_production <- fluxes_total_sum_compounds_production %>%
  filter(Reaction != "EX_cpd11416_c0")

# Prepare control data by averaging replicates
control_means <- fluxes_production_control %>%
  group_by(SampleName, Compartment, Reaction) %>%
  summarise(mean_Flux_control = mean(Flux, na.rm = TRUE), .groups = 'drop')

# Join control means with compounds data
combined_compounds_control <- fluxes_total_sum_compounds_production %>%
  inner_join(control_means, by = c("SampleName", "Compartment", "Reaction"))

# Calculate fold changes for each replicate
combined_compounds_control <- combined_compounds_control %>%
  mutate(FoldChange = Flux / mean_Flux_control)

# Flatten the data for easier manipulation
combined_data_long <- combined_compounds_control %>%
  select(SampleName, Compartment, medium, Reaction, Replicate, compound_amount, FoldChange)

# remove FoldChanges == NaN (in case when reaction flux is non-zero is in compounds data and zero in control)
# or Inf (if both fluxes are zero)
combined_data_long <- combined_data_long[which(!is.na(combined_data_long$FoldChange) & !(combined_data_long$FoldChange == Inf)),]



# Function to perform paired Wilcoxon test for each combination of compound_amount, Compartment, and Reaction
perform_paired_wilcox_test <- function(data, compound_amount) {
  compound_filtered <- data %>%
    filter(compound_amount == !!compound_amount)

  combined_data_filtered <- compound_filtered %>%
    group_by(Compartment, Reaction) %>%
    filter(n() >= 3 & sum(!is.na(FoldChange)) >= 3) %>%
    summarise(
      median_fold_change = median(FoldChange, na.rm = TRUE),
      mean_fold_change = mean(FoldChange, na.rm = TRUE),
      pval = ifelse(sum(!is.na(FoldChange)) >= 3,
                    wilcox.test(FoldChange, mu = 1, alternative = "two.sided")$p.value,
                    NA),
      .groups = 'drop'
    ) %>%
    filter(!is.na(pval)) %>%
    mutate(
      compound_amount = !!compound_amount
    )

  return(combined_data_filtered)
}

# Apply the function to each unique compound_amount
results_list <- lapply(unique(combined_data_long$compound_amount), perform_paired_wilcox_test,
                       data = combined_data_long)


# Combine results into a single data frame
flux_table_wilcox_control_compound <- bind_rows(results_list)


# Apply the Benjamini-Hochberg correction for multiple comparisons
flux_table_wilcox_control_compound_corrected <- flux_table_wilcox_control_compound %>%
  mutate(pval_adj = p.adjust(pval, method = "BH"))

# Add Metabolite names
flux_table_wilcox_control_compound_corrected <- flux_table_wilcox_control_compound_corrected %>%
  mutate(
    Metabolite = sapply(Reaction, function(x) {
      met_id <- gsub("EX_|_e0", "", x)
      if (met_id %in% seed_metabolites$id) {
        return(seed_metabolites$name[seed_metabolites$id == met_id])
      } else {
        return(NA)
      }
    })
  )

# Save the table
write.table(flux_table_wilcox_control_compound_corrected, paste0(path, "FluxTable_PairedWilcoxPvalsBH_FoldChange_allCompounds_", sim_dates_compounds,
                                                                 "_vs_control_", sim_dates_control, ".txt"), sep = '\t', quote = F)

# combine with combined_data_long for plotting the boxplots
flux_table_wilcox_control_compound_corrected_to_plot <- merge(combined_data_long, flux_table_wilcox_control_compound_corrected,
                                                              by = c("Compartment","Reaction","compound_amount"))

write.table(flux_table_wilcox_control_compound_corrected_to_plot, paste0(path, "FluxTable_ToPlot_PairedWilcoxPvalsBH_FoldChange_allCompounds_", sim_dates_compounds,
                                                                 "_vs_control_", sim_dates_control, ".txt"), sep = '\t', quote = F)




##################################################################################################
###############     REPEAT - BUT SPLIT BY DIET AND ONLY 10 REPS FROM EACH DIET (2 per sample):
##################################################################################################
# Function to sample exactly 10 replicates for each combination of Compartment and Reaction within each medium
# (groups by Reaction; the concentration version below groups by Metabolite instead)
sample_replicates <- function(data, n = 10) {
  data %>%
    group_by(compound_amount, medium, Compartment, Reaction) %>%
    group_modify(~ {
      # Check how many samples are available for the group
      available_samples <- unique(.x$SampleName)

      # Calculate the maximum number of replicates per sample to collect
      reps_per_sample <- ceiling(n / length(available_samples))

      # Sample the rows ensuring a total of `n` replicates are collected
      sampled_rows <- .x %>%
        group_by(SampleName) %>%
        sample_n(min(n(), reps_per_sample)) %>%
        ungroup()

      # If we sampled more than `n` rows (because of rounding up), downsample
      if (nrow(sampled_rows) > n) {
        sampled_rows <- sampled_rows %>% sample_n(n)
      }

      return(sampled_rows)
    }) %>%
    ungroup()
}


# Sample replicates in advance for control and compound data
combined_data_long_sampled_10reps <- sample_replicates(combined_data_long)


# Apply the function to each compound amount
results_list_byDiet_10reps <- lapply(unique(combined_data_long_sampled_10reps$compound_amount),
                                     perform_paired_wilcox_test_byDiet,
                                     data=combined_data_long_sampled_10reps)

# Combine results into a single data frame
flux_table_wilcox_control_compound_byDiet_10reps <- bind_rows(results_list_byDiet_10reps)

# Apply the Benjamini-Hochberg correction for multiple comparisons
flux_table_wilcox_control_compound_byDiet_10reps_corrected <- flux_table_wilcox_control_compound_byDiet_10reps %>%
  mutate(pval_adj = p.adjust(pval, method = "BH"))

# Add Metabolite names
flux_table_wilcox_control_compound_byDiet_10reps_corrected <- flux_table_wilcox_control_compound_byDiet_10reps_corrected %>%
  mutate(
    Metabolite = sapply(Reaction, function(x) {
      met_id <- gsub("EX_|_e0", "", x)
      if (met_id %in% seed_metabolites$id) {
        return(seed_metabolites$name[seed_metabolites$id == met_id])
      } else {
        return(NA)
      }
    })
  )

# Save the table
write.table(flux_table_wilcox_control_compound_byDiet_10reps_corrected,
            paste0(path, "FluxTable_PairedWilcoxPvalsBH_FoldChange_byDiet_10reps_allCompounds_", sim_dates_compounds,
                   "_vs_control_", sim_dates_control, ".txt"), sep = '\t', quote = F)

# combine with combined_data_long for plotting
flux_table_wilcox_control_compound_byDiet_10reps_corrected_to_plot <- merge(combined_data_long_sampled_10reps,
                                                                            flux_table_wilcox_control_compound_byDiet_10reps_corrected,
                                                              by = c("Compartment","Reaction","compound_amount","medium"))

write.table(flux_table_wilcox_control_compound_byDiet_10reps_corrected_to_plot,
            paste0(path, "FluxTable_ToPlot_PairedWilcoxPvalsBH_FoldChange_byDiet_10reps_allCompounds_", sim_dates_compounds,
                   "_vs_control_", sim_dates_control, ".txt"), sep = '\t', quote = F)



# ----------------------------------------------------------------------------------------
##########################################################################################
################        CONCENTRATIONS AT 96h WILCOX TEST           ######################
##########################################################################################
# ----------------------------------------------------------------------------------------

conc_table_control_h96 <- read.table(paste0(path_control,"conc_table_hour96_allMets_control_",
                                            pattern_control,
                                            "_", sim_dates_control, ".txt"),
                                     sep='\t',header=T)


conc_table_control_h96$compound_amount <- "control"
conc_table_control_h96 <- conc_table_control_h96 %>%
  filter(Metabolite != "EX_cpd11416_c0")


conc_table_compounds_h96 <- read.table(paste0(path,"conc_table_hour96_allMets_",
                                              pattern_compounds,
                                              "_", sim_dates_compounds, ".txt"),
                                       sep='\t',header=T)


# Create a new compound_amount column
conc_compounds <- conc_table_compounds_h96 %>%
  mutate(compound_amount = paste0(compound, "_", compound_gper100g)) %>%
  select(-compound, -compound_gper100g)

# Remove biomass
conc_compounds <- conc_compounds %>%
  filter(Metabolite != "EX_cpd11416_c0")

# Prepare control data by averaging replicates
conc_control_means <- conc_table_control_h96 %>%
  group_by(SampleName, Compartment, Metabolite) %>%
  summarise(mean_Conc_control = mean(Conc, na.rm = TRUE), .groups = 'drop')

# Join control means with compounds data
combined_compounds_control <- conc_compounds %>%
  inner_join(conc_control_means, by = c("SampleName", "Compartment", "Metabolite"))

# Calculate fold changes for each replicate
combined_compounds_control <- combined_compounds_control %>%
  mutate(FoldChange = Conc / mean_Conc_control)

# Flatten the data for easier manipulation
combined_data_long <- combined_compounds_control %>%
  select(SampleName, Compartment, medium, Metabolite, Replicate, compound_amount, FoldChange)

# remove FoldChanges == NaN (in case when reaction flux is non-zero is in compounds data and zero in control)
# or Inf (if both fluxes are zero)
combined_data_long <- combined_data_long[which(!is.na(combined_data_long$FoldChange) & !(combined_data_long$FoldChange == Inf)),]


# Function to perform paired Wilcoxon test for each combination of compound_amount, Compartment, and Metabolite
perform_paired_wilcox_test_conc <- function(data, compound_amount) {
  compound_filtered <- data %>%
    filter(compound_amount == !!compound_amount)

  combined_data_filtered <- compound_filtered %>%
    group_by(Compartment, Metabolite) %>%
    filter(n() >= 3 & sum(!is.na(FoldChange)) >= 3) %>%
    summarise(
      median_fold_change = median(FoldChange, na.rm = TRUE),
      mean_fold_change = mean(FoldChange, na.rm = TRUE),
      pval = ifelse(sum(!is.na(FoldChange)) >= 3,
                    wilcox.test(FoldChange, mu = 1, alternative = "two.sided")$p.value,
                    NA),
      .groups = 'drop'
    ) %>%
    filter(!is.na(pval)) %>%
    mutate(
      compound_amount = !!compound_amount
    )

  return(combined_data_filtered)
}

# Apply the function to each unique compound_amount
results_list <- lapply(unique(combined_data_long$compound_amount), perform_paired_wilcox_test_conc,
                       data = combined_data_long)


# Combine results into a single data frame
conc_table_wilcox_control_compound <- bind_rows(results_list)

# Apply the Benjamini-Hochberg correction for multiple comparisons
conc_table_wilcox_control_compound_corrected <- conc_table_wilcox_control_compound %>%
  group_by(medium) %>%
  mutate(pval_adj = p.adjust(pval, method = "BH"))

# Add Metabolite names
conc_table_wilcox_control_compound_corrected <- conc_table_wilcox_control_compound_corrected %>%
  mutate(
    MetName = sapply(Metabolite, function(x) {
      met_id <- gsub("EX_|_e0", "", x)
      if (met_id %in% seed_metabolites$id) {
        return(seed_metabolites$name[seed_metabolites$id == met_id])
      } else {
        return(NA)
      }
    })
  )


# Save the table
write.table(conc_table_wilcox_control_compound_corrected, paste0(path, "ConcTable_PairedWilcoxPvalsBH_FoldChange_allCompounds_",
                                                                 sim_dates_compounds,
                                                                 "_vs_control_", sim_dates_control, ".txt"), sep = '\t', quote = F)



##################################################################################################
###############     REPEAT - BUT SPLIT BY DIET AND ONLY 10 REPS FROM EACH DIET (2 per sample):
##################################################################################################
# Function to sample exactly 10 replicates for each combination of Compartment and Metabolite within each medium
# (groups by Metabolite; the flux version above groups by Reaction instead)
sample_replicates <- function(data, n = 10) {
  data %>%
    group_by(compound_amount, medium, Compartment, Metabolite) %>%
    group_modify(~ {
      # Check how many samples are available for the group
      available_samples <- unique(.x$SampleName)

      # Calculate the maximum number of replicates per sample to collect
      reps_per_sample <- ceiling(n / length(available_samples))

      # Sample the rows ensuring a total of `n` replicates are collected
      sampled_rows <- .x %>%
        group_by(SampleName) %>%
        sample_n(min(n(), reps_per_sample)) %>%
        ungroup()

      # If we sampled more than `n` rows (because of rounding up), downsample
      if (nrow(sampled_rows) > n) {
        sampled_rows <- sampled_rows %>% sample_n(n)
      }

      return(sampled_rows)
    }) %>%
    ungroup()
}

### DO FOR CORN ONLY:
combined_data_long <- combined_data_long[which(combined_data_long$medium == "corn"),]

for (sampling_time in c(1:5)) {

  # Sample replicates in advance for control and compound data
  combined_data_long_sampled_10reps <- sample_replicates(combined_data_long)


  # Apply the function to each compound amount
  results_list_byDiet_10reps <- lapply(unique(combined_data_long_sampled_10reps$compound_amount),
                                       perform_paired_wilcox_test_conc_byDiet,
                                       data=combined_data_long_sampled_10reps)

  # Combine results into a single data frame
  conc_table_wilcox_control_compound_byDiet_10reps <- bind_rows(results_list_byDiet_10reps)


  # Apply the Benjamini-Hochberg correction for multiple comparisons
  conc_table_wilcox_control_compound_byDiet_10reps_corrected <- conc_table_wilcox_control_compound_byDiet_10reps %>%
    group_by(medium) %>%
    mutate(pval_adj = p.adjust(pval, method = "BH"))


  # Add Metabolite names
  conc_table_wilcox_control_compound_byDiet_10reps_corrected <- conc_table_wilcox_control_compound_byDiet_10reps_corrected %>%
    mutate(
      MetName = sapply(Metabolite, function(x) {
        met_id <- gsub("EX_|_e0", "", x)
        if (met_id %in% seed_metabolites$id) {
          return(seed_metabolites$name[seed_metabolites$id == met_id])
        } else {
          return(NA)
        }
      })
    )

  # Save the table
  write.table(conc_table_wilcox_control_compound_byDiet_10reps_corrected,
              paste0(path, "ConcTable_PairedWilcoxPvalsBH_FoldChange_byDiet_10reps_allCompounds_", sim_dates_compounds,
                     "_vs_control_", sim_dates_control, "_CORN_sampled", sampling_time ,".txt"), sep = '\t', quote = F)
}
