###############################################################################
# CompileTables_fromOUTfiles_withCompounds.R
#
# Compile extracted simulation outputs (from Extract_OUTfiles scripts)
# into analysis-ready tables for COMPOUND/PREBIOTIC SUPPLEMENTATION
# simulations. Concatenates the 0-48h and 48-96h simulation halves,
# extracts SCFA fluxes and concentrations, and writes summary tables
# for downstream analyses (PCoA, fold change plots, etc.).
#
# USAGE (command line):
#   Rscript CompileTables_fromOUTfiles_withCompounds.R \
#       <pattern_48h> <pattern_48_96h> <sim_date_48h> <sim_date_48_96h>
#
#   pattern_48h:       file pattern for 0-48h simulation outputs
#   pattern_48_96h:    file pattern for 48-96h simulation outputs
#   sim_date_48h:      date subdirectory for 0-48h outputs
#   sim_date_48_96h:   date subdirectory for 48-96h outputs
###############################################################################

library(data.table)

myargs <- commandArgs(trailingOnly = TRUE)
print(myargs)

# ---- Input parameters ----
pattern_arg_48h    <- as.character(myargs[1])
pattern_arg_48_96h <- as.character(myargs[2])
sim_date_48h       <- as.character(myargs[3])
sim_date_48_96h    <- as.character(myargs[4])

# ---- Paths (adjust to your environment) ----
base_out_dir   <- "path/to/simulation_outputs/CoupledSixComps_D10/WithPrebiotics/2024/gapseq_models/OUT_files/"
path_to_tables <- paste0(base_out_dir, "Tables/")


###############################################################################
# HELPER: Melt relative abundance table from wide to long format
###############################################################################

meltTheRelabTable <- function(relab_table, control = TRUE) {
  if (!control) {
    relab_table <- reshape2::melt(relab_table,
      id.vars = c("Specie", "Sample", "medium", "essential_conc", "proportionMedium",
                   "Compartment", "Sugars_X", "AminoAcids_X", "compound", "compound_gper100g"))
    colnames(relab_table)[c(11, 12)] <- c("Hour", "RelAbundance")
  } else {
    relab_table <- reshape2::melt(relab_table,
      id.vars = c("Specie", "Sample", "medium", "essential_conc", "proportionMedium",
                   "Compartment", "Sugars_X", "AminoAcids_X"))
    colnames(relab_table)[c(9, 10)] <- c("Hour", "RelAbundance")
  }
  relab_table$Hour <- as.numeric(gsub('X', '', as.character(relab_table$Hour)))
  return(relab_table)
}


###############################################################################
# PART 1: CONCATENATE RELATIVE & ABSOLUTE ABUNDANCES (0-96h)
###############################################################################

setwd(base_out_dir)

relab_tables_48h    <- list.files(path = sim_date_48h,
                                  pattern = paste0("OUTfile_RelAbund_CommSize.*", pattern_arg_48h))
relab_tables_48_96h <- list.files(path = sim_date_48_96h,
                                  pattern = paste0("OUTfile_RelAbund_CommSize.*", pattern_arg_48_96h))

# ---- Load and concatenate 0-48h outputs ----
relab_tables_grower_48h    <- NULL
commsize_tables_grower_48h <- NULL

for (file in relab_tables_48h) {
  out_tables     <- readRDS(file.path(sim_date_48h, file))
  relab_table    <- meltTheRelabTable(out_tables[[2]], control = FALSE)
  commsize_table <- out_tables[[3]]

  relab_tables_grower_48h    <- rbind(relab_tables_grower_48h, relab_table)
  commsize_tables_grower_48h <- rbind(commsize_tables_grower_48h, commsize_table)
  commsize_tables_grower_48h$Hour <- as.numeric(commsize_tables_grower_48h$Hour)
}

# ---- Load 48-96h outputs, shift hours, extract hour 96 ----
absab_table_grower_h96       <- NULL
relab_table_grower_h96       <- NULL
relab_tables_grower_48_96h   <- NULL
commsize_tables_grower_48_96h <- NULL

for (file in relab_tables_48_96h) {
  out_tables     <- readRDS(file.path(sim_date_48_96h, file))
  relab_table    <- out_tables[[2]]
  commsize_table <- out_tables[[3]]
  absab_table    <- out_tables[[4]]

  # Process absolute abundances
  colnames(absab_table)[1:51] <- c(1:48, "Replicate", "Compartment", "Specie")
  absab_table_melt <- reshape2::melt(absab_table,
    id.vars = c("Specie", "Replicate", "Sample", "medium",
                "essential_conc", "proportionMedium", "Compartment",
                "Sugars_X", "AminoAcids_X", "compound", "compound_gper100g"))
  colnames(absab_table_melt)[c(12, 13)] <- c("Hour", "AbsAbundance")
  absab_table_melt$Replicate <- NULL
  absab_table_melt_mean <- aggregate(AbsAbundance ~ ., absab_table_melt, mean)

  # Process relative abundances
  relab_table <- meltTheRelabTable(relab_table, control = FALSE)

  # Shift hours: 1-48 -> 49-96
  relab_table$Hour <- relab_table$Hour + 48
  absab_table_melt_mean$Hour <- as.numeric(as.character(absab_table_melt_mean$Hour)) + 48

  # Extract hour 96 snapshot
  relab_table_h96      <- relab_table[relab_table$Hour == 96, ]
  absab_table_mean_h96 <- absab_table_melt_mean[absab_table_melt_mean$Hour == 96, ]

  commsize_table$Hour <- as.numeric(commsize_table$Hour) + 48

  relab_table_grower_h96        <- rbind(relab_table_grower_h96, relab_table_h96)
  absab_table_grower_h96        <- rbind(absab_table_grower_h96, absab_table_mean_h96)
  relab_tables_grower_48_96h    <- rbind(relab_tables_grower_48_96h, relab_table)
  commsize_tables_grower_48_96h <- rbind(commsize_tables_grower_48_96h, commsize_table)
}

# ---- Combine full 96h time series ----
relab_tables_grower_96h    <- rbind(relab_tables_grower_48h, relab_tables_grower_48_96h)
commsize_tables_grower_96h <- rbind(commsize_tables_grower_48h, commsize_tables_grower_48_96h)

# Write abundance tables
write.table(relab_table_grower_h96,
  paste0(path_to_tables, "relab_table_", pattern_arg_48_96h, "_h96_allCompounds_concatenated_", sim_date_48_96h, "_meanByRep.txt"),
  sep = '\t', quote = FALSE)
write.table(absab_table_grower_h96,
  paste0(path_to_tables, "absab_table_", pattern_arg_48_96h, "_h96_allCompounds_concatenated_", sim_date_48_96h, "_meanByRep.txt"),
  sep = '\t', quote = FALSE)
write.table(relab_tables_grower_96h,
  paste0(path_to_tables, "RelAbund_tables_96h_", pattern_arg_48_96h, "_allCompounds_", sim_date_48_96h, "_meanByRep.txt"),
  sep = '\t', quote = FALSE)
write.table(commsize_tables_grower_96h,
  paste0(path_to_tables, "CommSize_tables_96h_", pattern_arg_48_96h, "_allCompounds_", sim_date_48_96h, "_meanByRep.txt"),
  sep = '\t', quote = FALSE)


###############################################################################
# PART 2: COMPILE FLUX AND CONCENTRATION TABLES (0-96h)
###############################################################################

setwd(base_out_dir)

out_tables_grower_48h <- list.files(
  path = file.path(getwd(), sim_date_48h),
  pattern = paste0("OUTfile_Fluxes_Substances_Shadow.*", pattern_arg_48h))
out_tables_grower_48_96h <- list.files(
  path = file.path(getwd(), sim_date_48_96h),
  pattern = paste0("OUTfile_Fluxes_Substances_Shadow.*", pattern_arg_48_96h))

# SCFA reaction IDs (butyrate, acetate, propionate, lactate)
scfa_reactions <- c("EX_cpd00211_e0", "EX_cpd00029_e0", "EX_cpd00141_e0", "EX_cpd00159_e0")


#' Load and process flux/concentration OUTfiles for one time period
#' (compound supplementation version - includes compound and compound_gper100g columns)
get_flux_conc_tables <- function(sim_date, out_tables_files) {
  flux_scfa_byCommunity   <- data.table()
  flux_scfa_bySpecie_mean <- data.table()
  flux_allMets_byCommunity_sum      <- data.table()
  flux_allMets_byCommunity_eachHour <- data.table()
  conc_scfa    <- data.table()
  conc_allmets  <- data.table()

  for (out_table_file in out_tables_files) {
    start_time <- Sys.time()
    print(out_table_file)
    flux_tables <- readRDS(file.path(sim_date, out_table_file))
    # Structure: list of 6 compartments, each with 6 slots:
    #   1) all fluxes  2) sum of fluxes  3) mean by replicate
    #   4) all concentrations  5) mean concentrations  6) shadow prices

    flux_scfa_byCommunity_file   <- data.table()
    flux_scfa_bySpecie_mean_file <- data.table()
    conc_scfa_file    <- data.table()
    flux_allMets_eachHour_file <- data.table()
    flux_allMets_sum_file      <- data.table()

    for (comp_n in 1:6) {
      flux_table          <- as.data.table(flux_tables[[comp_n]][[1]])
      flux_table_sum      <- as.data.table(flux_tables[[comp_n]][[2]])
      flux_table_sum_mean <- as.data.table(flux_tables[[comp_n]][[3]])
      conc_table_all      <- as.data.table(flux_tables[[comp_n]][[4]])

      # SCFA-specific extractions
      flux_scfa_sum    <- flux_table_sum[Reaction %in% scfa_reactions]
      flux_scfa_mean   <- flux_table_sum_mean[Reaction %in% scfa_reactions]
      conc_scfa_allRep <- conc_table_all[Metabolite %in% scfa_reactions]

      # Community-level flux sums (including compound grouping columns)
      flux_community_eachHour <- flux_table[, .(Flux = sum(Flux)),
        by = .(Reaction, Replicate, Hour, compound, compound_gper100g,
               proportion_medium, medium, SampleName, Compartment)]
      flux_community_sum <- flux_table_sum[, .(Flux = sum(Flux)),
        by = .(Reaction, Replicate, compound, compound_gper100g,
               proportion_medium, medium, SampleName, Compartment)]
      flux_scfa_community <- flux_scfa_sum[, .(Flux = sum(Flux)),
        by = .(Reaction, Replicate, compound, compound_gper100g,
               proportion_medium, medium, SampleName, Compartment)]

      flux_scfa_byCommunity_file   <- rbindlist(list(flux_scfa_byCommunity_file, flux_scfa_community))
      flux_scfa_bySpecie_mean_file <- rbindlist(list(flux_scfa_bySpecie_mean_file, flux_scfa_mean))
      conc_scfa_file               <- rbindlist(list(conc_scfa_file, conc_scfa_allRep))
      flux_allMets_eachHour_file   <- rbindlist(list(flux_allMets_eachHour_file, flux_community_eachHour))
      flux_allMets_sum_file        <- rbindlist(list(flux_allMets_sum_file, flux_community_sum))
    }
    print(Sys.time() - start_time)

    flux_scfa_byCommunity   <- rbindlist(list(flux_scfa_byCommunity, flux_scfa_byCommunity_file))
    flux_scfa_bySpecie_mean <- rbindlist(list(flux_scfa_bySpecie_mean, flux_scfa_bySpecie_mean_file))
    conc_scfa               <- rbindlist(list(conc_scfa, conc_scfa_file))
    flux_allMets_byCommunity_eachHour <- rbindlist(list(flux_allMets_byCommunity_eachHour, flux_allMets_eachHour_file))
    flux_allMets_byCommunity_sum      <- rbindlist(list(flux_allMets_byCommunity_sum, flux_allMets_sum_file))
  }

  return(list(flux_scfa_byCommunity, flux_scfa_bySpecie_mean,
              conc_scfa, conc_allmets, flux_allMets_byCommunity_sum,
              flux_allMets_byCommunity_eachHour))
}


# ---- Process both time periods ----
tables_48h    <- get_flux_conc_tables(sim_date_48h, out_tables_grower_48h)
tables_48_96h <- get_flux_conc_tables(sim_date_48_96h, out_tables_grower_48_96h)

# ---- Combine 0-48h + 48-96h, summing fluxes across the two periods ----

# SCFA fluxes by community (sum over 96h)
flux_scfa_96h <- rbindlist(list(tables_48h[[1]], tables_48_96h[[1]]))[,
  .(Flux = sum(Flux)),
  by = .(Reaction, Replicate, compound, compound_gper100g,
         proportion_medium, medium, SampleName, Compartment)]

# SCFA fluxes by species, mean across replicates (sum over 96h)
flux_scfa_bySpecie_96h <- rbindlist(list(tables_48h[[2]], tables_48_96h[[2]]))[,
  .(Flux.mean = sum(Flux.mean)),
  by = .(Reaction, Specie, compound, compound_gper100g,
         proportion_medium, medium, SampleName, Compartment)]

# All metabolite fluxes by community (sum over 96h)
flux_allMets_96h <- rbindlist(list(tables_48h[[5]], tables_48_96h[[5]]))[,
  .(Flux = sum(Flux)),
  by = .(Reaction, Replicate, compound, compound_gper100g,
         proportion_medium, medium, SampleName, Compartment)]

# SCFA concentrations: shift 48-96h hours by +47
conc_scfa_48_96h      <- tables_48_96h[[3]]
conc_scfa_48_96h$Hour <- conc_scfa_48_96h$Hour + 47
conc_scfa_96h         <- rbindlist(list(tables_48h[[3]], conc_scfa_48_96h))

# Hourly flux sums: shift 48-96h hours by +48
flux_eachHour_48_96h      <- tables_48_96h[[6]]
flux_eachHour_48_96h$Hour <- flux_eachHour_48_96h$Hour + 48
flux_eachHour_96h         <- rbindlist(list(tables_48h[[6]], flux_eachHour_48_96h))


# ---- Write output tables ----
write.table(conc_scfa_96h[Hour == 95, ],
  paste0(path_to_tables, "conc_table_hour96_SCFA_allCompounds_", pattern_arg_48_96h, "_sim", sim_date_48_96h, ".txt"),
  sep = '\t', quote = FALSE)

write.table(flux_scfa_96h,
  paste0(path_to_tables, "flux_table_sum_SCFA_ByCommunity_allCompounds_", pattern_arg_48_96h, "_sim", sim_date_48_96h, ".txt"),
  sep = '\t', quote = FALSE)

write.table(flux_scfa_bySpecie_96h,
  paste0(path_to_tables, "flux_table_sum_SCFA_bySpecie_meanByRep_allCompounds_", pattern_arg_48_96h, "_sim", sim_date_48_96h, ".txt"),
  sep = '\t', quote = FALSE)

write.table(conc_scfa_96h,
  paste0(path_to_tables, "conc_table_96h_SCFA_allCompounds_", pattern_arg_48_96h, "_sim", sim_date_48_96h, ".txt"),
  sep = '\t', quote = FALSE)

write.table(flux_allMets_96h,
  paste0(path_to_tables, "flux_table_96h_sum_allMets_ByCommunity_allCompounds_", pattern_arg_48_96h, "_sim", sim_date_48_96h, ".txt"),
  sep = '\t', quote = FALSE)

write.table(flux_eachHour_96h,
  paste0(path_to_tables, "flux_table_eachHour_96h_sum_allMets_ByCommunity_allCompounds_", pattern_arg_48_96h, "_sim", sim_date_48_96h, ".txt"),
  sep = '\t', quote = FALSE)
