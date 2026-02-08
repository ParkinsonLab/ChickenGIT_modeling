# Extract abundance, flux, and concentration data from 2-compartment (Ileum/Cecum)
# simulation .rds files and write summary tables for downstream analysis.

library(data.table)
library(reshape2)
library(parallel)

myargs <- commandArgs(trailingOnly = TRUE)
if (length(myargs) < 2) stop("Usage: Rscript script.R <directory_pattern> <output_prefix> [num_cores]")

directory_pattern <- myargs[1]
output_prefix <- myargs[2]
num_cores <- ifelse(length(myargs) >= 3, as.numeric(myargs[3]), 6)

base_dir <- "Chicken_IleumCecum_2024Trial_simulation/simulation_outputs/gapseq_models/"
input_path <- file.path(base_dir, directory_pattern)
output_path <- file.path(base_dir, "OUT_files", "Tables")
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

compartments <- c("Ileum", "Cecum")

# Metabolomics metabolites
metabolomics_metabolites <- c(
  "Acetic acid"             = "EX_cpd00029_e0",
  "Propionic acid"          = "EX_cpd00141_e0",
  "Butyric acid"            = "EX_cpd00211_e0",
  "Lactic acid"             = "EX_cpd00159_e0",
  "Succinic acid"           = "EX_cpd00036_e0",
  "alpha-Ketoglutaric acid" = "EX_cpd00024_e0",
  "Fumaric acid"            = "EX_cpd00106_e0"
)

metabolomics_reactions <- as.character(metabolomics_metabolites)

# Extract additional parameters from filename (e.g., suffix after "L2reg_ileO2_")
extract_additional_params <- function(filepath) {
  filename <- basename(filepath)
  if (grepl("L2reg_ileO2_", filename)) {
    additional_part <- sub(".*L2reg_ileO2_(.*)\\.rds$", "\\1", filename)
    return(additional_part)
  } else if (grepl("L2reg_ileO2\\.rds$", filename)) {
    return("base")
  } else {
    return("unknown")
  }
}

process_single_file <- function(filepath) {
  simlist <- readRDS(filepath)
  n_reps <- length(simlist)

  # Extract parameters
  params <- simlist[[1]][[5]]
  essential_conc <- params[[1]]
  proportion_medium <- params[[2]]
  timesteps <- params[[3]]
  sample_name <- params[[4]]
  compound <- ifelse(is.null(params[[5]]), "None", params[[5]])
  compound_gper100g <- ifelse(is.null(params[[6]]), 0, params[[6]])

  additional_params <- extract_additional_params(filepath)

  all_abundance <- list()
  all_fluxes <- list()
  all_concentrations <- list()
  all_commsize <- list()

  for (rep in 1:n_reps) {
    for (comp_idx in 1:2) {
      compartment <- compartments[comp_idx]

      # Population data
      pop_data_hours <- list()
      for (hour in 1:timesteps) {
        if (!is.null(simlist[[rep]][[comp_idx]][[hour]]$Population)) {
          pop_matrix <- simlist[[rep]][[comp_idx]][[hour]]$Population
          pop_data_hours[[hour]] <- data.table(
            Specie = rownames(pop_matrix),
            AbsAbundance = pop_matrix[, 2],
            Hour = hour
          )
        }
      }

      if (length(pop_data_hours) > 0) {
        pop_combined <- rbindlist(pop_data_hours, fill = TRUE)
        pop_combined[is.na(AbsAbundance), AbsAbundance := 0]

        pop_combined[, RelAbundance := AbsAbundance / sum(AbsAbundance), by = Hour]
        pop_combined[is.na(RelAbundance), RelAbundance := 0]

        pop_combined[, `:=`(
          Replicate = rep, Compartment = compartment, Sample = sample_name,
          essential_conc = essential_conc, proportionMedium = proportion_medium,
          compound = compound, compound_gper100g = compound_gper100g,
          additional_params = additional_params
        )]

        all_abundance[[length(all_abundance) + 1]] <- pop_combined

        commsize_dt <- pop_combined[, .(CommSize = sum(AbsAbundance)),
                                    by = .(Hour, Replicate, Compartment, Sample, essential_conc,
                                           proportionMedium, compound, compound_gper100g, additional_params)]
        all_commsize[[length(all_commsize) + 1]] <- commsize_dt
      }

      # Fluxes
      for (hour in 2:timesteps) {
        if (!is.null(simlist[[rep]][[comp_idx]][[hour]]$Fluxlist)) {
          fluxlist <- simlist[[rep]][[comp_idx]][[hour]]$Fluxlist[[2]]

          flux_data <- rbindlist(lapply(names(fluxlist), function(bac) {
            if (length(names(fluxlist[[bac]])) != 0) {
              bac_fluxes <- fluxlist[[bac]]
              ex_reactions <- names(bac_fluxes)[startsWith(names(bac_fluxes), "EX_")]

              if (length(ex_reactions) > 0) {
                flux_values <- bac_fluxes[ex_reactions]
                non_zero <- abs(flux_values) > 1e-6
                if (sum(non_zero) > 0) {
                  data.table(
                    Reaction = ex_reactions[non_zero],
                    Flux = flux_values[non_zero],
                    Specie = bac,
                    Hour = hour, Replicate = rep, Compartment = compartment
                  )
                }
              }
            }
          }), fill = TRUE)

          if (!is.null(flux_data) && nrow(flux_data) > 0) {
            flux_data[, `:=`(
              Sample = sample_name, essential_conc = essential_conc,
              proportion_medium = proportion_medium, compound = compound,
              compound_gper100g = compound_gper100g, additional_params = additional_params
            )]
            all_fluxes[[length(all_fluxes) + 1]] <- flux_data
          }
        }

        # Concentrations
        if (!is.null(simlist[[rep]][[comp_idx]][[hour]]$Substances)) {
          substances <- simlist[[rep]][[comp_idx]][[hour]]$Substances
          conc_dt <- data.table(
            Metabolite = rownames(substances),
            Conc = substances[, 2],
            Hour = hour, Replicate = rep, Compartment = compartment,
            Sample = sample_name, essential_conc = essential_conc,
            proportion_medium = proportion_medium, compound = compound,
            compound_gper100g = compound_gper100g, additional_params = additional_params
          )
          all_concentrations[[length(all_concentrations) + 1]] <- conc_dt
        }
      }
    }
  }

  return(list(
    abundance = rbindlist(all_abundance, fill = TRUE),
    fluxes = rbindlist(all_fluxes, fill = TRUE),
    concentrations = rbindlist(all_concentrations, fill = TRUE),
    commsize = rbindlist(all_commsize, fill = TRUE)
  ))
}

# Main execution
files <- list.files(input_path, pattern = "\\.rds$", full.names = TRUE)
if (length(files) == 0) stop("No .rds files found in: ", input_path)

cat("Processing", length(files), "files with", num_cores, "cores\n")

cl <- makeCluster(num_cores, type = "PSOCK")
clusterExport(cl, c("process_single_file", "compartments", "metabolomics_reactions",
                     "metabolomics_metabolites", "extract_additional_params"))
clusterEvalQ(cl, { library(data.table); library(reshape2) })

file_results <- parLapply(cl, files, process_single_file)
stopCluster(cl)

# Combine all results
cat("Combining results and generating tables...\n")
all_abundance <- rbindlist(lapply(file_results, function(x) x$abundance), fill = TRUE)
all_fluxes <- rbindlist(lapply(file_results, function(x) x$fluxes), fill = TRUE)
all_concentrations <- rbindlist(lapply(file_results, function(x) x$concentrations), fill = TRUE)
all_commsize <- rbindlist(lapply(file_results, function(x) x$commsize), fill = TRUE)

date_str <- gsub("-", "", directory_pattern)

cat("Generating analysis tables...\n")

# 1. Relative abundance (mean across replicates)
if (nrow(all_abundance) > 0) {
  relabund_mean <- all_abundance[, .(RelAbundance = mean(RelAbundance, na.rm = TRUE)),
                                 by = .(Specie, Hour, Compartment, Sample, essential_conc,
                                        proportionMedium, compound, compound_gper100g, additional_params)]

  write.table(relabund_mean,
              file.path(output_path, paste0("RelAbund_tables_", output_prefix, "_", date_str, "_meanByRep.txt")),
              sep = '\t', quote = FALSE, row.names = FALSE)
}

# 2. Community size
if (nrow(all_commsize) > 0) {
  commsize_mean <- all_commsize[, .(CommSize = mean(CommSize, na.rm = TRUE)),
                                by = .(Hour, Compartment, Sample, essential_conc,
                                       proportionMedium, compound, compound_gper100g, additional_params)]

  write.table(commsize_mean,
              file.path(output_path, paste0("CommSize_tables_", output_prefix, "_", date_str, "_meanByRep.txt")),
              sep = '\t', quote = FALSE, row.names = FALSE)
}

# 3. Absolute abundance (all replicates)
if (nrow(all_abundance) > 0) {
  write.table(all_abundance,
              file.path(output_path, paste0("absab_table_", output_prefix, "_", date_str, "_allReplicates.txt")),
              sep = '\t', quote = FALSE, row.names = FALSE)
}

# 4. Metabolomics flux tables
if (nrow(all_fluxes) > 0) {
  metabolomics_fluxes <- all_fluxes[Reaction %in% metabolomics_reactions]

  if (nrow(metabolomics_fluxes) > 0) {
    metabolomics_fluxes[, Metabolite_Name := names(metabolomics_metabolites)[match(Reaction, metabolomics_metabolites)]]

    # By community (summed across species)
    metabolomics_community <- metabolomics_fluxes[, .(Flux = sum(Flux, na.rm = TRUE)),
                                                  by = .(Reaction, Metabolite_Name, Replicate, Hour, Compartment, Sample,
                                                         essential_conc, proportion_medium, compound, compound_gper100g, additional_params)]

    write.table(metabolomics_community,
                file.path(output_path, paste0("flux_table_sum_metabolomics_ByCommunity_", output_prefix, "_sim", date_str, ".txt")),
                sep = '\t', quote = FALSE, row.names = FALSE)

    # By species (mean across replicates)
    metabolomics_species <- metabolomics_fluxes[, .(Flux.mean = mean(Flux, na.rm = TRUE), Flux.sd = sd(Flux, na.rm = TRUE)),
                                                by = .(Reaction, Metabolite_Name, Specie, Compartment, Sample, essential_conc,
                                                       proportion_medium, compound, compound_gper100g, additional_params)]

    write.table(metabolomics_species,
                file.path(output_path, paste0("flux_table_sum_metabolomics_bySpecie_meanByRep_", output_prefix, "_sim", date_str, ".txt")),
                sep = '\t', quote = FALSE, row.names = FALSE)
  }

  # All metabolites flux by community
  all_flux_community <- all_fluxes[, .(Flux = sum(Flux, na.rm = TRUE)),
                                   by = .(Reaction, Replicate, Hour, Compartment, Sample,
                                          essential_conc, proportion_medium, compound, compound_gper100g, additional_params)]

  write.table(all_flux_community,
              file.path(output_path, paste0("flux_table_sum_allMets_ByCommunity_", output_prefix, "_sim", date_str, ".txt")),
              sep = '\t', quote = FALSE, row.names = FALSE)
}

# 5. Metabolomics concentration tables
if (nrow(all_concentrations) > 0) {
  metabolomics_conc <- all_concentrations[Metabolite %in% metabolomics_reactions]

  if (nrow(metabolomics_conc) > 0) {
    metabolomics_conc[, Metabolite_Name := names(metabolomics_metabolites)[match(Metabolite, metabolomics_metabolites)]]

    write.table(metabolomics_conc,
                file.path(output_path, paste0("conc_table_metabolomics_", output_prefix, "_sim", date_str, ".txt")),
                sep = '\t', quote = FALSE, row.names = FALSE)

    # Final timepoint concentrations
    max_hour <- max(metabolomics_conc$Hour, na.rm = TRUE)
    metabolomics_final <- metabolomics_conc[Hour == max_hour]

    write.table(metabolomics_final,
                file.path(output_path, paste0("conc_table_hour", max_hour, "_metabolomics_", output_prefix, "_sim", date_str, ".txt")),
                sep = '\t', quote = FALSE, row.names = FALSE)
  }
}

cat("Completed! Processed", length(files), "files\n")
cat("Tables saved to:", output_path, "\n")

cat("\n=== PROCESSING SUMMARY ===\n")
if (nrow(all_abundance) > 0) cat("Abundance records:", nrow(all_abundance), "\n")
if (nrow(all_commsize) > 0) cat("Community size records:", nrow(all_commsize), "\n")
if (nrow(all_fluxes) > 0) cat("Flux records:", nrow(all_fluxes), "\n")
if (nrow(all_concentrations) > 0) cat("Concentration records:", nrow(all_concentrations), "\n")

if (nrow(all_abundance) > 0) {
  unique_params <- unique(all_abundance$additional_params)
  cat("Additional parameter sets found:", paste(unique_params, collapse = ", "), "\n")
}
