###############################################################################
# Extract_OUTfiles_from_simulations_Controls.R
#
# Extract simulation output data (relative abundances, community sizes,
# fluxes, metabolite concentrations, shadow prices) from raw .rds
# simulation output files for CONTROL simulations (no added compounds).
#
# USAGE (command line):
#   Rscript Extract_OUTfiles_from_simulations_Controls.R <pattern_arg>
#
#   pattern_arg: regex pattern to match simulation .rds files
#                (e.g., "10samples_grower" to match all grower diet samples)
#
# OUTPUT: Two .rds files per simulation batch, saved to OUT_files/<sim_date>/:
#   - OUTfile_RelAbund_CommSize_*  (relative/absolute abundances + community size)
#   - OUTfile_Fluxes_Substances_Shadow_*  (fluxes, concentrations, shadow prices)
#
# This script is designed to run on a computing cluster
###############################################################################

library(data.table)
library(dplyr)
library(reshape2)
library(stringr)
library(parallel)

myargs <- commandArgs(trailingOnly = TRUE)
print(myargs)

# ---- Input arguments ----
pattern_arg <- as.character(myargs[1])

# ---- Paths ----
base_dir        <- "path/to/simulation_outputs/"
path            <- paste0(base_dir, "CoupledSixComps_D10/NoPathogens/2024/gapseq_models/")
ExchangeConstraints <- read.table(paste0(base_dir, "ModelGeneration_Files/Compounds_with_uptake_lbnds_gapseq.txt"),
                                  sep = '\t', header = TRUE)
AddedReactions  <- read.table(paste0(base_dir, "Reactions_toAdd_toModels_dbCAN3.txt"),
                              sep = '\t', header = TRUE)
cluster_log_dir <- paste0(base_dir, "cluster_logs/")

Compartments_list <- list("Gizzard", "Duodenum", "Jejunum", "Ileum", "Cecum", "Colon")


###############################################################################
# HELPER FUNCTIONS
###############################################################################

#' cbind data frames/matrices with different numbers of rows (pad with NA)
custom_cbind.fill <- function(...) {
  nm <- lapply(list(...), as.matrix)
  n  <- max(sapply(nm, nrow))
  do.call(cbind, lapply(nm, function(x) rbind(x, matrix(NA, n - nrow(x), ncol(x)))))
}

#' Add simulation parameter columns to a data table
#' @param dt         data.table to annotate
#' @param params     named list of simulation parameters
#' @param is_control logical: if TRUE, omit compound columns
add_sim_params <- function(dt, params, is_control = TRUE) {
  setDT(dt)
  dt[, c('essential_conc', 'proportion_medium', 'medium', 'SampleName') :=
       .(params$essential_conc, params$proportion_medium, params$medium, params$sample_name)]
  return(dt)
}


###############################################################################
# FUNCTION: Extract relative abundances and community sizes
###############################################################################

getRelAbundance_Commsize <- function(files) {
  final_relabund      <- data.table()
  final_relabund_mean <- data.table()
  commsize_total      <- data.table()
  final_absabund      <- data.table()

  for (file in files) {
    print(file)

    simlist <- readRDS(paste0(path, file))
    replicates_number <- length(simlist)

    # Extract simulation parameters from metadata slot
    params <- list(
      essential_conc    = simlist[[1]][[11]][[1]],
      proportion_medium = simlist[[1]][[11]][[2]],
      Sugars_X          = simlist[[1]][[11]][[3]] %||% "None",
      AminoAcids_X      = simlist[[1]][[11]][[4]] %||% "None",
      timesteps         = simlist[[1]][[11]][[5]],
      sample_name       = simlist[[1]][[11]][[6]],
      medium            = simlist[[1]][[11]][[7]]
    )
    timesteps <- params$timesteps

    # Helper to annotate data with simulation parameters
    annotate <- function(dt) {
      dt$essential_conc    <- params$essential_conc
      dt$proportionMedium  <- params$proportion_medium
      dt$medium            <- params$medium
      dt$Sample            <- params$sample_name
      dt$Sugars_X          <- params$Sugars_X
      dt$AminoAcids_X      <- params$AminoAcids_X
      return(dt)
    }

    tax_abund_total <- data.table()
    tax_relabund    <- data.table()
    commsize        <- data.table()

    for (rep in 1:replicates_number) {
      for (compartment_number in 1:6) {
        compartment <- Compartments_list[[compartment_number]]
        # Build population matrix across all time steps
        poplist <- as.data.frame(simlist[[rep]][[compartment_number]][[1]]$Population[, 2])

        for (hour in 2:timesteps) {
          poplist2 <- as.data.frame(simlist[[rep]][[compartment_number]][[hour]]$Population[, 2])
          # Handle case where new species appear in later time steps
          if (nrow(poplist2) > nrow(poplist)) {
            poplist_new <- custom_cbind.fill(poplist, poplist2)
            rownames(poplist_new)[(nrow(poplist) + 1):nrow(poplist2)] <-
              rownames(poplist2)[(nrow(poplist) + 1):nrow(poplist2)]
            poplist <- poplist_new
          } else {
            poplist <- cbind(poplist, poplist2)
          }
        }

        colnames(poplist) <- 1:timesteps
        poplist[is.na(poplist)] <- 0

        # Absolute abundances with replicate and compartment labels
        tax_abund_comp_rep <- cbind(poplist,
                                    rep(rep, nrow(poplist)),
                                    rep(compartment, nrow(poplist)))
        tax_abund_comp_rep <- as.data.frame(tax_abund_comp_rep)

        # Community size per time step
        commsize_comp_rep <- apply(tax_abund_comp_rep[, 1:timesteps], 2,
                                   function(x) sum(as.numeric(x)))
        commsize_comp_rep$Compartment <- compartment

        tax_abund_comp_rep$Specie <- rownames(poplist)
        colnames(tax_abund_comp_rep) <- colnames(tax_abund_total)

        tax_abund_total <- rbind(tax_abund_total, tax_abund_comp_rep)
        tax_abund_total[is.na(tax_abund_total)] <- 0

        # Relative abundances
        tax_abund_comp_rep[, 1:timesteps] <- apply(
          tax_abund_comp_rep[, 1:timesteps], 2,
          function(x) { as.numeric(x) / sum(as.numeric(x)) })
        tax_relabund <- rbind(tax_relabund, tax_abund_comp_rep)
        commsize     <- rbind(commsize, commsize_comp_rep)
        colnames(commsize) <- c(1:timesteps, "Compartment")
      }
    }

    colnames(tax_relabund) <- c(1:timesteps, "Replicate", "Compartment", "Specie")

    # Mean relative abundance across replicates
    tax_relabund_mean <- aggregate(. ~ Specie + Compartment, data = tax_relabund, FUN = mean)
    tax_relabund_mean$Replicate <- NULL

    tax_relabund      <- annotate(tax_relabund)
    tax_relabund_mean <- annotate(tax_relabund_mean)

    # Community size statistics (mean and sd across replicates)
    setDT(commsize)
    commsize_stat <- commsize[, as.list(unlist(lapply(.SD, function(x)
      list(mean = mean(x), sd = sd(x))))),
      by = "Compartment", .SDcols = 1:timesteps]
    commsize_stat_melt <- melt(commsize_stat)
    commsize_stat_melt[, stat := gsub("^[0-9]*\\.", "", variable)]
    commsize_stat_melt[, variable := gsub("\\.mean|\\.sd", "", variable)]
    commsize_stat_wide <- tidyr::spread(commsize_stat_melt, key = stat, value = value)
    colnames(commsize_stat_wide) <- c("Compartment", "Hour", "Commsize_mean", "Commsize_sd")
    commsize_stat_wide <- annotate(commsize_stat_wide)

    tax_abund_total <- annotate(tax_abund_total)

    final_absabund      <- rbind(final_absabund, tax_abund_total)
    final_relabund      <- rbind(final_relabund, tax_relabund)
    final_relabund_mean <- rbind(final_relabund_mean, tax_relabund_mean)
    commsize_total      <- rbind(commsize_total, commsize_stat_wide)

    setDT(final_absabund); setDT(final_relabund)
    setDT(final_relabund_mean); setDT(commsize_total)
  }

  return(list(final_relabund, final_relabund_mean, commsize_total, final_absabund))
}


###############################################################################
# FUNCTION: Extract fluxes, concentrations, and shadow prices (per compartment)
###############################################################################

getFluxes_Conc_Shadows <- function(files, compartment_number) {
  bac_fluxes_total_df          <- data.table()
  bac_fluxes_total_df_sum      <- data.table()
  bac_fluxes_total_sum_mean    <- data.table()
  substance_table_df_total     <- data.table()
  substance_table_df_total_mean <- data.table()
  bac_shadowprices_total_df    <- data.table()

  compartment <- Compartments_list[[compartment_number]]
  print(compartment)

  for (file in files) {
    start_time <- Sys.time()
    print(file)

    simlist <- readRDS(paste0(path, file))
    replicates_number <- length(simlist)

    # Extract parameters
    essential_conc    <- simlist[[1]][[11]][[1]]
    proportion_medium <- simlist[[1]][[11]][[2]]
    Sugars_X          <- simlist[[1]][[11]][[3]] %||% "None"
    AminoAcids_X      <- simlist[[1]][[11]][[4]] %||% "None"
    timesteps         <- simlist[[1]][[11]][[5]]
    sample_name       <- simlist[[1]][[11]][[6]]
    medium            <- simlist[[1]][[11]][[7]]

    process_data <- function(dt) {
      setDT(dt)
      dt[, c('essential_conc', 'proportion_medium', 'medium', 'SampleName') :=
           .(essential_conc, proportion_medium, medium, sample_name)]
      return(dt)
    }

    bac_fluxes_total_file       <- data.table()
    substance_table_total_file  <- data.table()
    bac_shadowprices_total_file <- data.table()

    for (rep in 1:replicates_number) {
      for (hour in 2:timesteps) {
        fluxlist     <- simlist[[rep]][[compartment_number]][[hour]]$Fluxlist[[2]]
        subconclist  <- simlist[[rep]][[compartment_number]][[hour]]$Substances[, 2]
        shadowlist_rep <- simlist[[rep]][[compartment_number]][[hour]]$Shadowlist

        # Extract exchange and added reaction fluxes for each bacterium
        bac_fluxes_total <- rbindlist(lapply(names(fluxlist), function(bac) {
          if (length(names(fluxlist[[bac]])) == 0) return(NULL)
          bac_reactions <- names(fluxlist[[bac]])
          bac_name <- gsub("(-gapseqAdRm|_cobrapy_adapted)", "", bac)
          is_exchange <- startsWith(bac_reactions, "EX_")
          is_added    <- bac_reactions %in% paste0(
            AddedReactions$reaction_id[grep(bac_name, AddedReactions$specie_models)], "_c0")
          keep <- is_exchange | is_added

          data.table(
            Reaction    = bac_reactions[keep],
            Flux        = fluxlist[[bac]][keep],
            Specie      = bac_name,
            Replicate   = rep,
            Hour        = hour,
            Compartment = compartment
          )
        }))

        # Remove zero-flux rows to save memory
        if (nrow(bac_fluxes_total) > 0) {
          bac_fluxes_total <- bac_fluxes_total[Flux != 0, ]
        }
        bac_fluxes_total_file <- rbind(bac_fluxes_total_file, bac_fluxes_total)

        # Shadow prices: top 5 most limiting metabolites per species
        bac_shadowprices <- shadowlist_rep %>% group_by(Specie) %>% top_n(-5, value)
        bac_shadowprices <- data.table(
          Metabolite  = bac_shadowprices$Metabolite,
          Shadow      = bac_shadowprices$value,
          Specie      = bac_shadowprices$Specie,
          Replicate   = rep,
          Hour        = hour,
          Compartment = compartment
        )
        bac_shadowprices_total_file <- rbind(bac_shadowprices_total_file, bac_shadowprices)

        # Substance (metabolite) concentrations
        substance_table <- data.table(
          Metabolite  = names(subconclist),
          Conc        = subconclist,
          Replicate   = rep,
          Hour        = hour,
          Compartment = compartment
        )
        substance_table_total_file <- rbindlist(list(substance_table_total_file, substance_table))
      }
    }

    end_time <- Sys.time()
    print(end_time - start_time)

    # Summarize: sum fluxes over all hours per reaction/specie/replicate
    bac_fluxes_total_file[, Flux := as.numeric(as.character(Flux))]

    bac_fluxes_total_file_sum <- aggregate(
      Flux ~ Reaction + Specie + Replicate + Compartment,
      data = bac_fluxes_total_file, FUN = sum)
    setDT(bac_fluxes_total_file_sum)

    bac_fluxes_total_file_sum <- process_data(bac_fluxes_total_file_sum)
    bac_fluxes_total_file     <- process_data(bac_fluxes_total_file)

    # Mean and sd of total fluxes across replicates
    bac_fluxes_total_file_sum_mean <- bac_fluxes_total_file_sum[,
      list(Flux.mean = mean(Flux), Flux.sd = sd(Flux)),
      by = c("Reaction", "Specie", "essential_conc", "Compartment",
             "proportion_medium", "medium", "SampleName")]
    setDT(bac_fluxes_total_file_sum_mean)

    substance_table_total_file <- process_data(substance_table_total_file)
    substance_table_total_file_mean <- substance_table_total_file[,
      list(Conc.mean = mean(Conc), Conc.sd = sd(Conc)),
      by = c("Metabolite", "Hour", "essential_conc", "Compartment",
             "proportion_medium", "medium", "SampleName")]

    bac_shadowprices_total_file <- process_data(bac_shadowprices_total_file)

    # Accumulate across files
    bac_fluxes_total_df       <- rbind(bac_fluxes_total_df, bac_fluxes_total_file)
    bac_fluxes_total_df_sum   <- rbind(bac_fluxes_total_df_sum, bac_fluxes_total_file_sum)
    bac_fluxes_total_sum_mean <- rbind(bac_fluxes_total_sum_mean, bac_fluxes_total_file_sum_mean)
    substance_table_df_total  <- rbind(substance_table_df_total, substance_table_total_file)
    substance_table_df_total_mean <- rbind(substance_table_df_total_mean, substance_table_total_file_mean)
    bac_shadowprices_total_df <- rbind(bac_shadowprices_total_df, bac_shadowprices_total_file)
  }

  return(list(
    bac_fluxes_total_df, bac_fluxes_total_df_sum, bac_fluxes_total_sum_mean,
    substance_table_df_total, substance_table_df_total_mean,
    bac_shadowprices_total_df
  ))
}


###############################################################################
# MAIN: Run extraction
###############################################################################

files <- list.files(path, pattern = pattern_arg)
date_of_simulation <- strsplit(
  strsplit(files[1], "_")[[1]][length(strsplit(files[1], "_")[[1]]) - 1],
  "\\.rds")[[1]][1]
sim_date <- gsub("-", "", date_of_simulation)
dir.create(paste0(path, "OUT_files/", sim_date), showWarnings = FALSE)

# ---- Step 1: Extract relative abundances and community sizes ----
final_relabund_commsize <- getRelAbundance_Commsize(files = files)
saveRDS(final_relabund_commsize,
        paste0(path, "OUT_files/", sim_date, "/OUTfile_RelAbund_CommSize_",
               pattern_arg, "_allReplicates_andMean.rds"))

# ---- Step 2: Extract fluxes, concentrations, shadows (parallelized by compartment) ----
cores <- 6
cl <- makeCluster(cores, type = "PSOCK",
                  outfile = paste0(cluster_log_dir, "Retrieve_tables_", pattern_arg, ".log"))
clusterExport(cl, c("getFluxes_Conc_Shadows", "path", "files",
                     "ExchangeConstraints", "AddedReactions", "Compartments_list"))
clusterEvalQ(cl, library("data.table"))
clusterEvalQ(cl, library("dplyr"))

print(system.time(outlist <- parLapply(cl, 1:6, function(i) {
  getFluxes_Conc_Shadows(files = files, compartment_number = i)
})))
stopCluster(cl)

saveRDS(outlist,
        paste0(path, "OUT_files/", sim_date, "/OUTfile_Fluxes_Substances_Shadow_",
               pattern_arg, "_allReplicates_andMean.rds"))
