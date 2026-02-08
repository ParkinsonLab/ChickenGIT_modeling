# =============================================================================
# Simulation of 6-compartment chicken GIT model (48-96h continuation, with compounds)
# =============================================================================
# Continues a BacArena simulation from hour 48 to 96 for conditions
# with dietary compound/prebiotic supplementation. Loads the saved
# 48h output with compound, recreates arenas, and runs additional
# timesteps with the same simulation logic.
#
# USAGE:   Rscript ...continueAfter48h_addCompounds.R <sample> <timesteps>
#            <pattern> <sim_date> <compound> <compound_gper100g>
#
# INPUTS:  Previously saved 48h compound output, CompoundsTable, same base inputs
# =============================================================================

myargs <- commandArgs(trailingOnly = TRUE)

library(plyr)
library(dplyr)
library(magrittr)
library(BacArena)
library(glpkAPI)
library(igraph)
library(reshape2)
library(libSBML)
library(sybilSBML)
library(stringr)
library(parallel)

# ===== Base directory =====
base_dir <- "Chicken_10D_simulation/"
cluster_log_dir <- file.path(base_dir, "coupling_cluster_logs")

# ===== Loading input parameters =============================

# sample name, i.e. Corn10-1
sample <- as.character(myargs[1])

# number of timesteps
timesteps <- as.numeric(myargs[2])

# pattern of the filename of the 48h output file
pattern <- as.character(myargs[3])

# simulation date in the YYYY-MM-DD format
sim_date <- as.character(myargs[4])

# name of the compound/prebiotic to be tested (i.e., Folate)
compound <- as.character(myargs[5])

# the amount (in g per 100g) to be added to the feeding
compound_gper100g <- as.numeric(myargs[6])

#### ===== Load input files ================================
diet_corn <- read.table(file.path(base_dir, "corn_diet_grower_aggregated_2023.txt"),
                        sep = '\t', header = TRUE)
diet_wheat <- read.table(file.path(base_dir, "wheat_diet_grower_aggregated_2023.txt"),
                         sep = '\t', header = TRUE)
OTU_table <- read.table(file.path(base_dir, "OTUtable.tsv"))
metadata <- read.table(file.path(base_dir, "metadata.txt"), header = TRUE)
Parameters_table <- readxl::read_xlsx(file.path(base_dir, "ModelGeneration_Files/Table_parameters_coupling_24h_2024.xlsx"),
                                      sheet = "Table_parameters_16h_2023")
Parameters_table_night <- readxl::read_xlsx(file.path(base_dir, "ModelGeneration_Files/Table_parameters_coupling_24h_2024.xlsx"),
                                            sheet = "Table_parameters_nightmode")
Metabolite_conc_h24 <- read.delim(file.path(base_dir, "ModelGeneration_Files/MetaboliteConcentrations_h24_10Samples_growerdiet_econc1_propGizzard100_CobaltX12_rmAnaero_less02_mucUreaCecaColon_ordered_fromsim20240701.txt"),
                                  sep = '\t', header = TRUE)

CompoundsToTest_table <- read.table(file.path(base_dir, "ModelGeneration_Files/CompoundsTable_to_add_withPrebiotics_v2024.txt"),
                                    sep = '\t', header = TRUE)

ExchangeConstraints <- read.table(file.path(base_dir, "ModelGeneration_Files/Compounds_with_uptake_lbnds_gapseq.txt"),
                                  sep = '\t', header = TRUE)

# load all sbml models
model_list_all_compartments <-
  readRDS(file.path(base_dir, "ModelGeneration_Files/ChickenAllSpecies_noAGP_D10_gapseq_ModelList_MAGs_gapseqCobraAdapt_Aug2023.rds"))

# load sbml models to construct a community
model_files <- list.files(file.path(base_dir, "ModelGeneration_Files"),
                          pattern = "Chicken_.*_ModelList_MAGs_gapseqCobraAdapt_Aug2023.rds",
                          full.names = TRUE)
# rearrange in the GIT order
model_files <- model_files[c(4, 3, 6, 5, 1, 2)]
model_lists <- lapply(model_files, readRDS)
names(model_lists) <- c("Gizzard", "Duodenum", "Jejunum", "Ileum", "Ceca", "Colon")

# Load names of organisms in the community model and corresponding tax category
organism_files <- list.files(file.path(base_dir, "ModelGeneration_Files"),
                             pattern = "Chicken_.*_organisms_tax-category_table_gapseq_MAGs_gapseqCobraAdapt_Aug2023.rds",
                             full.names = TRUE)
# rearrange in the GIT order
organism_files <- organism_files[c(4, 3, 6, 5, 1, 2)]
organism_lists <- lapply(organism_files, readRDS)

# load all the functions
functions_scripts <- list.files(file.path(base_dir, "ModelGeneration_scripts"),
                                pattern = "Function.*.R", full.names = TRUE)

# Select diet based on sample name
diet_name <- tolower(str_remove(sample, '[0-9]{1,2}-[0-9]'))
if (diet_name == "wheat") {
  diet <- diet_wheat
} else {
  diet <- diet_corn
}

# ===== Compound modification of diet =============================
# get only the row of the table corresponding to the compound chosen
CompoundsToTest_table_sub <- CompoundsToTest_table[which(
  CompoundsToTest_table$variable == compound &
    CompoundsToTest_table$value == compound_gper100g), ]

# now add this row to the diet table, so that the compound is introduced into feeding:
diet <- rbind(diet, CompoundsToTest_table_sub)
# and sum the values for the compound if it was already present in the diet:
diet <- aggregate(cbind(value, met_mM) ~ unit + met_name + met_id + category, diet, sum)

# ===== Load corresponding output file after 48h of simulation =====
year <- strsplit(sim_date, '-')[[1]][1]
out_path <- file.path(base_dir, "arena_sim_objects/CoupledSixComps_D10/WithPrebiotics", year, "gapseq_models")
outfile <- readRDS(file.path(out_path,
                             paste0(sample, "_", diet_name, "_",
                                    pattern, "_", compound, "_",
                                    as.character(compound_gper100g), "gper100g_",
                                    sim_date, "_gapseq.rds")))

econc_proportion <- outfile[[1]][[11]]$econc_proportion
prop_gizzard <- outfile[[1]][[11]]$prop_gizzard
Sugars_X <- outfile[[1]][[11]]$Sugars_X
AminoAcid_X <- outfile[[1]][[11]]$AminoAcids_X

##########################################################
############ UPLOAD DATA ON CLUSTER ######################
############# RUN MODELS WITH REPLICATES #################
##########################################################

# Initialize cluster
cores <- 5
replicates <- length(outfile)

now <- Sys.time()
today <- strsplit(as.character(as.POSIXlt(now)[1]), "\\ ")[[1]][1]

cl <- makeCluster(cores, type = "PSOCK",
                  outfile = file.path(cluster_log_dir,
                                      paste0("cluster_coupling_6comps_10d_continueAfter48h_GAPSEQ_",
                                             sample, "_growerdiet_",
                                             diet_name, "_propEconc", econc_proportion, "_",
                                             compound, "_", as.character(compound_gper100g), "gper100g",
                                             "_propGizzard_mucusUrea_CoAupd", prop_gizzard, "_",
                                             today,
                                             "_uptakeConstraints_sugarsX", Sugars_X,
                                             "_AminoAcidsX", AminoAcid_X, ".log")))

clusterExport(cl, c("functions_scripts"))
# load all the scripts with functions
clusterCall(cl, function() { lapply(functions_scripts, source) })

# export all the parameters and input files
clusterExport(cl, c("prop_gizzard", "Parameters_table", "Parameters_table_night",
                     "econc_proportion", "ExchangeConstraints"))
clusterExport(cl, list("organism_lists", "model_lists", "model_list_all_compartments",
                        "diet", "diet_name", "diet_corn", "diet_wheat", "sample",
                        "Metabolite_conc_h24"))
clusterExport(cl, c("OTU_table", "metadata"))
clusterExport(cl, c("timesteps", "Sugars_X", "AminoAcid_X", "outfile",
                     "compound", "compound_gper100g"))

print(system.time(simlist <- parLapply(cl, 1:replicates, function(i) {

  # ===== Helper functions =====

  # Save simulation outputs for one compartment at one timestep
  save_sim_outputs <- function(sim_object, out_list, timestep) {
    curves <- BacArena::plotCurves(sim_object, retdata = TRUE, graph = FALSE)
    out_list[[timestep]]$Population <- curves$Population
    out_list[[timestep]]$Substances <- curves$Substances
    out_list[[timestep]]$Orgdat <- sim_object@simlist
    out_list[[timestep]]$Fluxlist <- sim_object@mfluxlist
    # save only negative shadow prices
    shadowlist <- reshape2::melt(unlist(lapply(sim_object@shadowlist[[2]], function(x) x[x < -0.1])))
    shadowlist$Metabolite <- unlist(lapply(rownames(shadowlist), function(x) { strsplit(x, "EX_")[[1]][2] }))
    shadowlist$Specie <- unlist(lapply(rownames(shadowlist), function(x) { gsub("\\.$", "", strsplit(x, "EX_")[[1]][1]) }))
    out_list[[timestep]]$Shadowlist <- shadowlist
    return(out_list)
  }

  # Update type-specie association table with any new species in an arena
  update_type_specie_table <- function(arena, type_specie_table) {
    new_types <- setdiff(unique(arena@orgdat$type), type_specie_table$Type)
    if (length(new_types) > 0) {
      for (bactype in new_types) {
        type_specie_new <- data.frame(
          Name = arena@specs[[bactype]]@model@mod_desc,
          Type = bactype
        )
        type_specie_table <- rbind(type_specie_table, type_specie_new)
      }
    }
    return(type_specie_table)
  }

  # Remove specified organisms (by name) from an arena
  remove_organisms_from_arena <- function(arena, organism_names, type_specie_table, compartment_label) {
    for (org_name in organism_names) {
      if (org_name %in% as.character(type_specie_table$Name)) {
        org_type <- as.numeric(type_specie_table$Type[as.character(type_specie_table$Name) == org_name])
        if (length(which(arena@orgdat$type == org_type)) > 0) {
          arena@orgdat <- arena@orgdat[-which(arena@orgdat$type == org_type), ]
        }
      }
    }
    return(arena)
  }

  # Split removal results from ileum into cecum and colon fractions
  split_removal_to_cecum_colon <- function(removefunction_result, fraction_to_cecum, fraction_to_colon) {
    result_to_cecum <- removefunction_result
    result_to_colon <- removefunction_result
    # sample fraction_to_cecum of orgdat removed from ileum
    result_to_cecum[[2]] <- result_to_cecum[[2]][sample(nrow(result_to_cecum[[2]]),
                                                         ceiling(fraction_to_cecum * nrow(result_to_cecum[[2]]))), ]
    # remaining goes to colon
    result_to_colon[[2]] <- dplyr::setdiff(removefunction_result[[2]], result_to_cecum[[2]])
    # do same for sublb (metabolites) table
    result_to_cecum[[3]] <- result_to_cecum[[3]][sample(nrow(result_to_cecum[[3]]),
                                                         round(fraction_to_colon * nrow(result_to_cecum[[3]]))), ]
    result_to_colon[[3]] <- dplyr::setdiff(removefunction_result[[3]], result_to_cecum[[3]])
    # update slot [[4]] - species_rm table
    type_number_to_cecum <- table(result_to_cecum[[2]]$type)
    type_number_to_colon <- table(result_to_colon[[2]]$type)
    result_to_cecum[[4]] <- subset(result_to_cecum[[4]], result_to_cecum[[4]]$Type %in% names(type_number_to_cecum))
    result_to_colon[[4]] <- subset(result_to_colon[[4]], result_to_colon[[4]]$Type %in% names(type_number_to_colon))
    result_to_cecum[[4]]$Number <- type_number_to_cecum
    result_to_colon[[4]]$Number <- type_number_to_colon
    return(list(to_cecum = result_to_cecum, to_colon = result_to_colon))
  }

  # ===== Replicate-specific setup =====

  # Recover seed from previous 48h simulation
  seedN <- outfile[[i]][[7]]

  Reflux <- function(sim_object_n1, perc_to_remove_bac_n1, perc_to_remove_mets_n1 = 0,
                     type_specie_table_n1,
                     sim_object_n0, arena2_n0, model_list_n0, model_list_n1) {
    # if perc_to_remove_mets_n1 not specified set it equal to perc_to_remove_bac
    if (perc_to_remove_mets_n1 == 0) {
      perc_to_remove_mets_n1 <- perc_to_remove_bac_n1
    }

    # REVERSE PERISTALSIS: COMP N+1 - remove; COMP N  - add
    # update arena of compartment N+1
    removefunction_result_n1 <- RemoveBacMetabolites(sim_object = sim_object_n1,
                                                     perc_to_remove_bac = perc_to_remove_bac_n1,
                                                     perc_to_remove_mets = perc_to_remove_mets_n1,
                                                     type_specie_table = type_specie_table_n1)
    arena2_n1_upd <- removefunction_result_n1[[1]]

    # then add retrieved bacs and mets to the compartment N:
    arena2_n0_upd <- AddBacMets_fromDownstreamCompartment(removefunction_result_n1,
                                                          sim_object_n0, arena2_n0, model_list_n0, model_list_n1)
    return(list("arena2_n0_upd" = arena2_n0_upd, "arena2_n1_upd" = arena2_n1_upd))
  }

  set_pH <- function(arena, pH) {
    arena <- BacArena::addSubs(arena, smax = 1000 / (10^pH), mediac = "EX_cpd00067_e0", unit = "mM", add = FALSE)
    return(arena)
  }

  set_o2 <- function(arena, max_conc) {
    arena <- BacArena::addSubs(arena, smax = max_conc, mediac = "EX_cpd00007_e0", unit = "mM", add = FALSE)
    return(arena)
  }

  if (Sugars_X == "None") {
    Sugars_X <- NULL
  }
  if (AminoAcid_X == "None") {
    AminoAcid_X <- NULL
  }

  # Update exchange reactions constraints based on the input constraints table for each modelList of each compartment
  model_lists <- lapply(model_lists, ChangeUptakeConstraints,
                        constraints_table = ExchangeConstraints,
                        Sugars_X = Sugars_X, AA_X = AminoAcid_X)

  stir_TF <- Parameters_table$stir[1]
  rich_TF <- FALSE
  diet_added_perc <- as.numeric(Parameters_table$diet_added_perc[Parameters_table$Compartment == "Gizzard"])
  diet_concentration <- as.numeric(Parameters_table$diet_proportion[Parameters_table$Compartment == "Gizzard"])

  # strict anaerobes that are not expected to survive in high oxygen environment like gizzard
  strict_anaerobes <- c("Clostridium_saudiense_cobrapy_adapted",
                        "Lachnoclostridium_phocaeensis_cobrapy_adapted",
                        "Blautia_sp_cobrapy_adapted",
                        "Mediterraneibacter_sp002161355_cobrapy_adapted",
                        "Faecalibacterium_sp002160895_cobrapy_adapted",
                        "Oscillospiraceae_genusUBA9475_cobrapy_adapted",
                        "Flavonifractor_plautii-gapseqAdRm")

  # obligate aerobes that are not expected to survive in nearly anaerobic environment (cecum and colon)
  obligate_aerobes <- c("Acinetobacter_lwoffii-gapseqAdRm")

  # ===== Recreate arenas from 48h output =====

  # Recreate arenas from last timesteps of the simulation for compartments
  RecreateModelArena <- function(compartment_number, sim_outfile, hour = 48,
                                 seedN, stir_var, tstep = 1, arena_size = 100) {
    # create the list with compartment names
    Compartments_list <- list("Gizzard", "Duodenum", "Jejunum", "Ileum", "Cecum", "Colon")
    # get the abundances after 48h
    AbundanceTable <- as.data.frame(sim_outfile[[compartment_number]][[hour]]$Population)
    AbundanceTable$Specie <- gsub("_cobrapy_adapted|-gapseqAdRm", "", rownames(AbundanceTable))

    # for gizzard:
    essential_metabolites <- c()

    # initialize arena
    arena <- BacArena::Arena(n = arena_size, m = arena_size, stir = stir_var,
                             Lx = 0.025 * (arena_size / 100), Ly = 0.025 * (arena_size / 100),
                             tstep = tstep, seed = seedN)
    # add species in the abundances they were present after 48h:
    # go through all species present in GIT:
    for (m in 1:length(model_list_all_compartments)) {
      # if abundance > 0:
      if (model_list_all_compartments[[m]]@mod_name %in% AbundanceTable$Specie) {
        if (AbundanceTable$V2[AbundanceTable$Specie == model_list_all_compartments[[m]]@mod_name] > 0) {

          bac <- BacArena::Bac(model = model_list_all_compartments[[m]], growtype = "exponential")
          arena <- BacArena::addOrg(arena, bac,
                                    amount = AbundanceTable$V2[AbundanceTable$Specie == model_list_all_compartments[[m]]@mod_name])

          ### only for gizzard get list of essential metabolites to be added along with the diet:
          if (compartment_number == 1) {
            essential_metabolites <- unique(c(essential_metabolites, BacArena::addEssentialMed(arena, bac, only_return = TRUE)))
          }
        }
      }
    }

    ###  now add metabolites in the concentrations at 48h:
    ConcTable <- as.data.frame(sim_outfile[[compartment_number]][[hour]]$Substances)
    ConcTable$Metabolite <- rownames(ConcTable)
    ConcTable <- ConcTable[which(ConcTable$V2 > 0), ]

    # remove biomass:
    ConcTable <- ConcTable[-which(ConcTable$Metabolite == "EX_cpd11416_c0"), ]

    # convert mmol/cell to mM to be added to arena:
    ConcTable$Conc_mM <- ConcTable$V2 / 625
    # add to arena:
    arena <- BacArena::addSubs(arena, smax = ConcTable$Conc_mM, unit = "mM",
                               mediac = ConcTable$Metabolite, add = FALSE, addAnyway = TRUE)

    # add nucleobases and essential CoA:
    # gapseq: add cytosine, uracil, cytidine, adenosine, guanosine, CoA
    arena <- BacArena::addSubs(arena, smax = 1e-04, unit = "mM",
                               mediac = c("EX_cpd00307_e0", "EX_cpd00092_e0", "EX_cpd00367_e0",
                                          "EX_cpd00182_e0", "EX_cpd00311_e0", "EX_cpd00010_e0"),
                               add = TRUE, addAnyway = TRUE)

    # add mucins and urea for cecum and colon:
    if (compartment_number %in% c(5, 6)) {
      mucins_urea <- c("EX_cpd00122_e0", "EX_cpd00232_e0", "EX_cpd00832_e0", "EX_cpd02992_e0",
                       "EX_cpd11842_e0", "EX_cpd21520_e0", "EX_cpd00073_e0")

      arena <- BacArena::addSubs(arena, smax = 1e-05, unit = "mM",
                                 mediac = mucins_urea, add = TRUE, addAnyway = TRUE)
    }

    return(list(arena, essential_metabolites))
  }

  compartment_arenas <- lapply(1:6, function(j) {
    RecreateModelArena(compartment_number = j, sim_outfile = outfile[[i]], hour = 48,
                       seedN = seedN, stir_var = stir_TF)[[1]]
  })

  essential_mets_gizzard <- RecreateModelArena(compartment_number = 1, sim_outfile = outfile[[i]],
                                               seedN = seedN, stir_var = stir_TF, tstep = 1)[[2]]
  # save only essential metabolites that aren't in the diet:
  essential_mets_gizzard <- essential_mets_gizzard[which(!(essential_mets_gizzard %in% diet$met_id))]

  # Initialize output lists
  out_lists <- lapply(1:6, function(x) vector("list", length = timesteps + 1))

  for (timestep in 1:timesteps) {

    names(compartment_arenas) <- c("Gizzard", "Duodenum", "Jejunum", "Ileum", "Cecum", "Colon")

    ########################################################################
    ############                  DAY TIME                     #############
    ########################################################################

    if (timestep %in% c(1:16, 25:40, 49:64)) {

      # set pH and oxygen
      for (j in 1:6) {
        compartment_arenas[[j]] <- set_pH(compartment_arenas[[j]], Parameters_table$pH[Parameters_table$Compartment == names(compartment_arenas)[j]])
        compartment_arenas[[j]] <- set_o2(compartment_arenas[[j]], max_conc = as.numeric(Parameters_table$oxygen_mM[Parameters_table$Compartment == names(compartment_arenas)[j]]))
      }

      ####################################
      ############## GIZZARD #############
      ####################################
      ### simulate gizzard for 1 h
      sim_gizzard <- BacArena::simEnv(compartment_arenas[[1]], time = 1)
      out_lists[[1]] <- save_sim_outputs(sim_gizzard, out_lists[[1]], timestep)

      arena2_gizzard <- BacArena::getArena(sim_gizzard, 1)

      # get type_specie_table - species and corresponding "type" in arena@orgdat
      if (timestep == 1) {
        type_specie_table_gizzard <- CreateTypeSpecieAssoc_table(sim_gizzard)
      }

      reflux_duo_frequency <- as.numeric(as.character(Parameters_table$reflux_freq_hrs[Parameters_table$Compartment == "Duodenum"]))

      # once in X hrs add fraction of diet (incl. essential mets) + initial composition to the grid:
      if (timestep %% as.numeric(Parameters_table$diet_added_hrs[1]) == 0) {
        arena2_gizzard <- AddBacDiet_gizzard(Chicken_organisms_gizzard = organism_lists[[1]],
                                             model_list_gizzard = model_lists[[1]], sampleName = sample,
                                             diet_mets = diet,
                                             arena_gizzard = arena2_gizzard, prop_gizzard = prop_gizzard,
                                             diet_added_perc = diet_added_perc,
                                             initial_commsize_gizzard = Parameters_table$initial_community_size[1],
                                             econc = econc_proportion * Parameters_table$essential_conc_mM[1],
                                             all_conc = Parameters_table$all_conc_mM[1],
                                             essential_mets = essential_mets_gizzard)
      }

      # select % of the grid randomly that will be removed and transferred to duodenum
      # although once every N hrs reflux happens --> reverse peristalsis duodenum->gizzard:
      if (!(timestep %% reflux_duo_frequency == 0)) {
        perc_to_remove_gizzard_bac <- Parameters_table$remove_bac_perc[Parameters_table$Compartment == "Gizzard"]
        perc_to_remove_gizzard_mets <- Parameters_table$remove_sublb_perc[Parameters_table$Compartment == "Gizzard"]

        removefunction_result_gizzard <- RemoveBacMetabolites(sim_gizzard, arena2_object = arena2_gizzard,
                                                              perc_to_remove_bac = perc_to_remove_gizzard_bac,
                                                              perc_to_remove_mets = perc_to_remove_gizzard_mets,
                                                              type_specie_table = type_specie_table_gizzard)
        arena2_gizzard <- removefunction_result_gizzard[[1]]
      }

      # update gizzard arena
      compartment_arenas[[1]] <- arena2_gizzard

      ####################################
      ############# DUODENUM #############
      ####################################
      if (timestep > 1) {
        # imitate absorption of sugars and conversion of starches in duodenum:
        compartment_arenas[[2]] <- AbsorptionConversion_duodenum(arena2 = compartment_arenas[[2]])
      }

      ### simulate duodenum for 1 h
      sim_duodenum <- BacArena::simEnv(compartment_arenas[[2]], time = 1)
      out_lists[[2]] <- save_sim_outputs(sim_duodenum, out_lists[[2]], timestep)

      # get type_specie_table - species and corresponding "type" in arena@orgdat
      if (timestep == 1) {
        type_specie_table_duodenum <- CreateTypeSpecieAssoc_table(sim_duodenum)
      }

      # NO REFLUX
      reflux_duo <- FALSE
      if (!(timestep %% reflux_duo_frequency == 0)) {
        perc_to_remove_duodenum_bac <- Parameters_table$remove_bac_perc[Parameters_table$Compartment == "Duodenum"]
        perc_to_remove_duodenum_mets <- Parameters_table$remove_sublb_perc[Parameters_table$Compartment == "Duodenum"]

        arena2_duodenum <- getArena(sim_duodenum, 1)

        removefunction_result_duodenum <- RemoveBacMetabolites(sim_duodenum, arena2_object = arena2_duodenum,
                                                               perc_to_remove_bac = perc_to_remove_duodenum_bac,
                                                               perc_to_remove_mets = perc_to_remove_duodenum_mets,
                                                               type_specie_table = type_specie_table_duodenum)
        arena2_duodenum <- removefunction_result_duodenum[[1]]

        # add mets and bacs that were removed from gizzard
        arena2_duodenum <- AddBacMets_fromUpstreamCompartment(removefunction_result_n0 = removefunction_result_gizzard,
                                                              sim_object_n1 = sim_duodenum, arena2_n1 = arena2_duodenum,
                                                              model_list_n1 = model_lists$Duodenum,
                                                              model_list_n0 = model_lists$Gizzard)

        # REFLUX
      } else if (timestep %% reflux_duo_frequency == 0) {
        reflux_duo <- TRUE
        perc_to_reflux_duodenum <- as.numeric(Parameters_table$reflux_perc[Parameters_table$Compartment == "Duodenum"])
        compartment_arenas[[1]] <- Reflux(sim_object_n1 = sim_duodenum, perc_to_remove_bac_n1 = perc_to_reflux_duodenum,
                                          type_specie_table_n1 = type_specie_table_duodenum,
                                          sim_object_n0 = sim_gizzard, arena2_n0 = arena2_gizzard,
                                          model_list_n0 = model_lists$Gizzard, model_list_n1 = model_lists$Duodenum)[[1]]
        arena2_duodenum <- Reflux(sim_object_n1 = sim_duodenum, perc_to_remove_bac_n1 = perc_to_reflux_duodenum,
                                  type_specie_table_n1 = type_specie_table_duodenum,
                                  sim_object_n0 = sim_gizzard, arena2_n0 = arena2_gizzard,
                                  model_list_n0 = model_lists$Gizzard, model_list_n1 = model_lists$Duodenum)[[2]]

        # if any new species happen to be in gizzard - add them to the type_specie_table_gizzard:
        type_specie_table_gizzard <- update_type_specie_table(compartment_arenas[[1]], type_specie_table_gizzard)

        # remove strict anaerobes (if any happened to be transferred during reflux):
        compartment_arenas[[1]] <- remove_organisms_from_arena(compartment_arenas[[1]], strict_anaerobes,
                                                               type_specie_table_gizzard, "gizzard")
      }

      # update duodenal arena
      compartment_arenas[[2]] <- arena2_duodenum

      ####################################
      ############# JEJUNUM ##############
      ####################################
      if (timestep > 2) {
        ## imitate absorption of sugars, amino acids and fats in jejunum:
        compartment_arenas[[3]] <- AbsorptionConversion_jejunum(arena2 = compartment_arenas[[3]], diet)
      }

      ### simulate jejunum for 1 h
      sim_jejunum <- BacArena::simEnv(compartment_arenas[[3]], time = 1)
      out_lists[[3]] <- save_sim_outputs(sim_jejunum, out_lists[[3]], timestep)

      # get type_specie_table
      if (timestep == 1) {
        type_specie_table_jejunum <- CreateTypeSpecieAssoc_table(sim_jejunum)
      }

      arena2_jejunum <- getArena(sim_jejunum, 1)

      perc_to_remove_jejunum_bac <- Parameters_table$remove_bac_perc[Parameters_table$Compartment == "Jejunum"]
      perc_to_remove_jejunum_mets <- Parameters_table$remove_sublb_perc[Parameters_table$Compartment == "Jejunum"]
      removefunction_result_jejunum <- RemoveBacMetabolites(sim_jejunum, arena2_object = arena2_jejunum,
                                                            perc_to_remove_bac = perc_to_remove_jejunum_bac,
                                                            perc_to_remove_mets = perc_to_remove_jejunum_mets,
                                                            type_specie_table = type_specie_table_jejunum)
      arena2_jejunum <- removefunction_result_jejunum[[1]]

      # if no reflux - add bacs and mets from duodenum:
      if (reflux_duo == FALSE) {
        arena2_jejunum <- AddBacMets_fromUpstreamCompartment(removefunction_result_n0 = removefunction_result_duodenum,
                                                             sim_object_n1 = sim_jejunum, arena2_n1 = arena2_jejunum,
                                                             model_list_n1 = model_lists$Jejunum,
                                                             model_list_n0 = model_lists$Duodenum)
      }
      # update jejunum arena
      compartment_arenas[[3]] <- arena2_jejunum

      ###################################
      ############## ILEUM ##############
      ###################################
      if (timestep > 2) {
        ## imitate absorption of sugars, amino acids and fats in ileum:
        compartment_arenas[[4]] <- AbsorptionConversion_ileum(arena2 = compartment_arenas[[4]], diet)
      }

      ### simulate ileum for 1 h
      sim_ileum <- BacArena::simEnv(compartment_arenas[[4]], time = 1)
      out_lists[[4]] <- save_sim_outputs(sim_ileum, out_lists[[4]], timestep)

      # get type_specie_table
      if (timestep == 1) {
        type_specie_table_ileum <- CreateTypeSpecieAssoc_table(sim_ileum)
      }

      arena2_ileum <- getArena(sim_ileum, 1)

      ### select % of the grid randomly that will be removed in the next iteration
      ### split between cecum and colon
      perc_to_remove_ileum_bac <- Parameters_table$remove_bac_perc[Parameters_table$Compartment == "Ileum"]
      perc_to_remove_ileum_mets <- Parameters_table$remove_sublb_perc[Parameters_table$Compartment == "Ileum"]

      perc_to_remove_ileum_to_cecum <- as.numeric(Parameters_table$remove_from_ileum[Parameters_table$Compartment == "Cecum"])
      fraction_to_cecum <- perc_to_remove_ileum_to_cecum / perc_to_remove_ileum_bac
      perc_to_remove_ileum_to_colon <- as.numeric(Parameters_table$remove_from_ileum[Parameters_table$Compartment == "Colon"])
      fraction_to_colon <- perc_to_remove_ileum_to_colon / perc_to_remove_ileum_bac

      removefunction_result_ileum <- RemoveBacMetabolites(sim_ileum, arena2_object = arena2_ileum,
                                                          perc_to_remove_bac = perc_to_remove_ileum_bac,
                                                          perc_to_remove_mets = perc_to_remove_ileum_mets,
                                                          type_specie_table = type_specie_table_ileum)
      arena2_ileum <- removefunction_result_ileum[[1]]

      # split removal results between cecum and colon
      split_result <- split_removal_to_cecum_colon(removefunction_result_ileum, fraction_to_cecum, fraction_to_colon)
      removefunction_result_ileum_to_cecum <- split_result$to_cecum
      removefunction_result_ileum_to_colon <- split_result$to_colon

      # add bacs and mets from jejunum:
      arena2_ileum <- AddBacMets_fromUpstreamCompartment(removefunction_result_n0 = removefunction_result_jejunum,
                                                         sim_object_n1 = sim_ileum, arena2_n1 = arena2_ileum,
                                                         model_list_n1 = model_lists$Ileum,
                                                         model_list_n0 = model_lists$Jejunum)

      # update ileal arena
      compartment_arenas[[4]] <- arena2_ileum

      ####################################
      ############  CECUM  ###############
      ####################################
      ### simulate cecum for 1 h
      sim_cecum <- BacArena::simEnv(compartment_arenas[[5]], time = 1)
      out_lists[[5]] <- save_sim_outputs(sim_cecum, out_lists[[5]], timestep)

      if (timestep == 1) {
        type_specie_table_cecum <- CreateTypeSpecieAssoc_table(sim_cecum)
      }

      # select pre-set % of the grid randomly that will be removed once in 4 hours (released to colon)
      if (timestep %% 4 == 0) {
        perc_to_remove_cecum_bac <- Parameters_table$remove_bac_perc[Parameters_table$Compartment == "Cecum"]
        perc_to_remove_cecum_mets <- Parameters_table$remove_sublb_perc[Parameters_table$Compartment == "Cecum"]
        removefunction_result_cecum <- RemoveBacMetabolites(sim_cecum,
                                                            perc_to_remove_bac = perc_to_remove_cecum_bac,
                                                            perc_to_remove_mets = perc_to_remove_cecum_mets,
                                                            type_specie_table = type_specie_table_cecum)
        arena2_cecum <- removefunction_result_cecum[[1]]
      } else {
        arena2_cecum <- BacArena::getArena(sim_cecum, 1)
      }

      # add metabolites and bacteria removed from the ileum
      arena2_cecum <- AddBacMets_fromUpstreamCompartment(removefunction_result_n0 = removefunction_result_ileum_to_cecum,
                                                         sim_object_n1 = sim_cecum, arena2_n1 = arena2_cecum,
                                                         model_list_n1 = model_lists$Ceca,
                                                         model_list_n0 = model_lists$Ileum)

      if (timestep > 1) {
        # if any new species happen to be transferred in cecum - add them to the type_specie_table:
        type_specie_table_cecum <- update_type_specie_table(arena2_cecum, type_specie_table_cecum)
      }

      # remove obligate aerobes (if any happened to be transferred):
      arena2_cecum <- remove_organisms_from_arena(arena2_cecum, obligate_aerobes, type_specie_table_cecum, "cecum")

      # update cecum arena
      compartment_arenas[[5]] <- arena2_cecum

      ####################################
      ############  COLON  ###############
      ####################################
      ### simulate colon for 1 h
      sim_colon <- BacArena::simEnv(compartment_arenas[[6]], time = 1)
      out_lists[[6]] <- save_sim_outputs(sim_colon, out_lists[[6]], timestep)

      # get type_specie_table
      if (timestep == 1) {
        type_specie_table_colon <- CreateTypeSpecieAssoc_table(sim_colon)
      }

      perc_to_remove_colon_bac <- Parameters_table$remove_bac_perc[Parameters_table$Compartment == "Colon"]
      perc_to_remove_colon_mets <- Parameters_table$remove_sublb_perc[Parameters_table$Compartment == "Colon"]
      removefunction_result_colon <- RemoveBacMetabolites(sim_colon,
                                                          perc_to_remove_bac = perc_to_remove_colon_bac,
                                                          perc_to_remove_mets = perc_to_remove_colon_mets,
                                                          type_specie_table = type_specie_table_colon)
      arena2_colon <- removefunction_result_colon[[1]]

      # once in 4 hrs cecal contents are released to the colon:
      if (timestep %% 4 == 0) {
        arena2_colon <- AddBacMets_fromUpstreamCompartment(removefunction_result_n0 = removefunction_result_cecum,
                                                           sim_object_n1 = sim_colon, arena2_n1 = arena2_colon,
                                                           model_list_n1 = model_lists$Colon,
                                                           model_list_n0 = model_lists$Ceca)
      }
      # add metabolites and bacteria removed from the ileum
      arena2_colon <- AddBacMets_fromUpstreamCompartment(removefunction_result_n0 = removefunction_result_ileum_to_colon,
                                                         sim_object_n1 = sim_colon, arena2_n1 = arena2_colon,
                                                         model_list_n1 = model_lists$Colon,
                                                         model_list_n0 = model_lists$Ileum)

      if (timestep > 1) {
        # if any new species happen to be transferred in colon - add them to the type_specie_table:
        type_specie_table_colon <- update_type_specie_table(arena2_colon, type_specie_table_colon)
      }

      # remove obligate aerobes (if any happened to be transferred):
      arena2_colon <- remove_organisms_from_arena(arena2_colon, obligate_aerobes, type_specie_table_colon, "colon")

      ### absorb water:
      sublb_colon <- as.data.frame(BacArena::getSublb(arena2_colon))
      # leave 20% of water:
      sublb_colon[, "EX_cpd00001_e0"] <- 0.2 * sublb_colon[1, "EX_cpd00001_e0"]
      arena2_colon@sublb <- as.matrix(sublb_colon)

      # update colon arena
      compartment_arenas[[6]] <- arena2_colon
    }


    ########################################################################
    ############                  NIGHT TIME                   #############
    ########################################################################

    if (timestep %in% c(17:24, 41:48, 65:72)) {

      # set pH and oxygen
      for (j in 1:6) {
        compartment_arenas[[j]] <- set_pH(compartment_arenas[[j]], Parameters_table$pH[Parameters_table$Compartment == names(compartment_arenas)[j]])
        compartment_arenas[[j]] <- set_o2(compartment_arenas[[j]], as.numeric(Parameters_table$oxygen_mM[Parameters_table$Compartment == names(compartment_arenas)[j]]))
      }

      ####################################
      ############## GIZZARD #############
      ####################################
      ### simulate gizzard for 1 h
      sim_gizzard <- BacArena::simEnv(compartment_arenas[[1]], time = 1)
      out_lists[[1]] <- save_sim_outputs(sim_gizzard, out_lists[[1]], timestep)

      reflux_freq_hrs_duo <- as.numeric(Parameters_table_night$reflux_freq_hrs[Parameters_table_night$Compartment == "Duodenum"])

      # reflux happens on every 3rd hour of the night mode 
      ## NO REFLUX
      if (!(timestep %% reflux_freq_hrs_duo == 1)) {
        perc_to_remove_gizzard_bac <- Parameters_table_night$remove_bac_perc[Parameters_table_night$Compartment == "Gizzard"]
        perc_to_remove_gizzard_mets <- Parameters_table_night$remove_sublb_perc[Parameters_table_night$Compartment == "Gizzard"]
        removefunction_result_gizzard <- RemoveBacMetabolites(sim_gizzard,
                                                              perc_to_remove_bac = perc_to_remove_gizzard_bac,
                                                              perc_to_remove_mets = perc_to_remove_gizzard_mets,
                                                              type_specie_table = type_specie_table_gizzard)
        arena2_gizzard <- removefunction_result_gizzard[[1]]
      } else {
        # IF REFLUX
        arena2_gizzard <- BacArena::getArena(sim_gizzard, 1)
      }

      # update gizzard arena
      compartment_arenas[[1]] <- arena2_gizzard

      ####################################
      ############# DUODENUM #############
      ####################################
      # imitate absorption of sugars and conversion of starches in duodenum:
      compartment_arenas[[2]] <- AbsorptionConversion_duodenum(arena2 = compartment_arenas[[2]])

      ### simulate duodenum for 1 h
      sim_duodenum <- BacArena::simEnv(compartment_arenas[[2]], time = 1)
      out_lists[[2]] <- save_sim_outputs(sim_duodenum, out_lists[[2]], timestep)

      # NO REFLUX
      if (!(timestep %% as.numeric(Parameters_table_night$reflux_freq_hrs[Parameters_table_night$Compartment == "Duodenum"]) == 1)) {
        perc_to_remove_duodenum_bac <- Parameters_table_night$remove_bac_perc[Parameters_table_night$Compartment == "Duodenum"]
        perc_to_remove_duodenum_mets <- Parameters_table_night$remove_sublb_perc[Parameters_table_night$Compartment == "Duodenum"]

        arena2_duodenum <- getArena(sim_duodenum, 1)

        removefunction_result_duodenum <- RemoveBacMetabolites(sim_duodenum, arena2_object = arena2_duodenum,
                                                               perc_to_remove_bac = perc_to_remove_duodenum_bac,
                                                               perc_to_remove_mets = perc_to_remove_duodenum_mets,
                                                               type_specie_table = type_specie_table_duodenum)
        arena2_duodenum <- removefunction_result_duodenum[[1]]

        # add mets and bacs that were removed from gizzard
        arena2_duodenum <- AddBacMets_fromUpstreamCompartment(removefunction_result_n0 = removefunction_result_gizzard,
                                                              sim_object_n1 = sim_duodenum, arena2_n1 = arena2_duodenum,
                                                              model_list_n1 = model_lists$Duodenum,
                                                              model_list_n0 = model_lists$Gizzard)

        # REFLUX
      } else if (timestep %% as.numeric(Parameters_table_night$reflux_freq_hrs[Parameters_table_night$Compartment == "Duodenum"]) == 1) {
        reflux_duo <- TRUE
        perc_to_reflux_duodenum <- as.numeric(Parameters_table_night$reflux_perc[Parameters_table_night$Compartment == "Duodenum"])
        compartment_arenas[[1]] <- Reflux(sim_object_n1 = sim_duodenum, perc_to_remove_bac_n1 = perc_to_reflux_duodenum,
                                          perc_to_remove_mets_n1 = perc_to_reflux_duodenum,
                                          type_specie_table_n1 = type_specie_table_duodenum,
                                          sim_object_n0 = sim_gizzard, arena2_n0 = arena2_gizzard,
                                          model_list_n0 = model_lists$Gizzard, model_list_n1 = model_lists$Duodenum)[[1]]
        arena2_duodenum <- Reflux(sim_object_n1 = sim_duodenum, perc_to_remove_bac_n1 = perc_to_reflux_duodenum,
                                  perc_to_remove_mets_n1 = perc_to_reflux_duodenum,
                                  type_specie_table_n1 = type_specie_table_duodenum,
                                  sim_object_n0 = sim_gizzard, arena2_n0 = arena2_gizzard,
                                  model_list_n0 = model_lists$Gizzard, model_list_n1 = model_lists$Duodenum)[[2]]

        # if any new species happen to be in gizzard - add them to the type_specie_table_gizzard:
        type_specie_table_gizzard <- update_type_specie_table(compartment_arenas[[1]], type_specie_table_gizzard)

        # remove strict anaerobes (if any happened to be transferred during reflux):
        compartment_arenas[[1]] <- remove_organisms_from_arena(compartment_arenas[[1]], strict_anaerobes,
                                                               type_specie_table_gizzard, "gizzard")
      }

      # update duodenal arena
      compartment_arenas[[2]] <- arena2_duodenum

      ####################################
      ############# JEJUNUM ##############
      ####################################
      # imitate absorption of sugars, amino acids and fats and conversion of starches in jejunum:
      compartment_arenas[[3]] <- AbsorptionConversion_jejunum(arena2 = compartment_arenas[[3]], diet)

      ### simulate jejunum for 1 h
      sim_jejunum <- BacArena::simEnv(compartment_arenas[[3]], time = 1)
      out_lists[[3]] <- save_sim_outputs(sim_jejunum, out_lists[[3]], timestep)

      # get type_specie_table
      if (timestep == 1) {
        type_specie_table_jejunum <- CreateTypeSpecieAssoc_table(sim_jejunum)
      }

      arena2_jejunum <- getArena(sim_jejunum, 1)

      # NO REFLUX
      if (!(timestep %% as.numeric(Parameters_table_night$reflux_freq_hrs[Parameters_table_night$Compartment == "Jejunum"]) == 1)) {

        perc_to_remove_jejunum_bac <- Parameters_table_night$remove_bac_perc[Parameters_table_night$Compartment == "Jejunum"]
        perc_to_remove_jejunum_mets <- Parameters_table_night$remove_sublb_perc[Parameters_table_night$Compartment == "Jejunum"]
        removefunction_result_jejunum <- RemoveBacMetabolites(sim_jejunum, arena2_object = arena2_jejunum,
                                                              perc_to_remove_bac = perc_to_remove_jejunum_bac,
                                                              perc_to_remove_mets = perc_to_remove_jejunum_mets,
                                                              type_specie_table = type_specie_table_jejunum)
        arena2_jejunum <- removefunction_result_jejunum[[1]]

        # if no reflux - add bacs and mets from duodenum:
        arena2_jejunum <- AddBacMets_fromUpstreamCompartment(removefunction_result_n0 = removefunction_result_duodenum,
                                                             sim_object_n1 = sim_jejunum, arena2_n1 = arena2_jejunum,
                                                             model_list_n1 = model_lists$Jejunum,
                                                             model_list_n0 = model_lists$Duodenum)

        # REFLUX
      } else if (timestep %% as.numeric(Parameters_table_night$reflux_freq_hrs[Parameters_table_night$Compartment == "Jejunum"]) == 1) {
        perc_to_reflux_jejunum <- as.numeric(Parameters_table_night$reflux_perc[Parameters_table_night$Compartment == "Jejunum"])
        compartment_arenas[[2]] <- Reflux(sim_object_n1 = sim_jejunum, perc_to_remove_bac_n1 = perc_to_reflux_jejunum,
                                          perc_to_remove_mets_n1 = perc_to_reflux_jejunum,
                                          type_specie_table_n1 = type_specie_table_jejunum,
                                          sim_object_n0 = sim_duodenum, arena2_n0 = arena2_duodenum,
                                          model_list_n0 = model_lists$Duodenum, model_list_n1 = model_lists$Jejunum)[[1]]

        arena2_jejunum <- Reflux(sim_object_n1 = sim_jejunum, perc_to_remove_bac_n1 = perc_to_reflux_jejunum,
                                 perc_to_remove_mets_n1 = perc_to_reflux_jejunum,
                                 type_specie_table_n1 = type_specie_table_jejunum,
                                 sim_object_n0 = sim_duodenum, arena2_n0 = arena2_duodenum,
                                 model_list_n0 = model_lists$Duodenum, model_list_n1 = model_lists$Jejunum)[[2]]
      }

      # update jejunum arena
      compartment_arenas[[3]] <- arena2_jejunum

      ###################################
      ############## ILEUM ##############
      ###################################
      # imitate absorption of sugars, amino acids and fats; and conversion of starches in ileum:
      compartment_arenas[[4]] <- AbsorptionConversion_ileum(arena2 = compartment_arenas[[4]], diet)

      ### simulate ileum for 1 h
      sim_ileum <- BacArena::simEnv(compartment_arenas[[4]], time = 1)
      out_lists[[4]] <- save_sim_outputs(sim_ileum, out_lists[[4]], timestep)

      # get type_specie_table
      if (timestep == 1) {
        type_specie_table_ileum <- CreateTypeSpecieAssoc_table(sim_ileum)
      }

      arena2_ileum <- getArena(sim_ileum, 1)

      # NO REFLUX
      if (!(timestep %% as.numeric(Parameters_table_night$reflux_freq_hrs[Parameters_table_night$Compartment == "Ileum"]) == 1)) {

        ### select % of the grid randomly that will be removed, split between cecum and colon
        perc_to_remove_ileum_bac <- Parameters_table_night$remove_bac_perc[Parameters_table_night$Compartment == "Ileum"]
        perc_to_remove_ileum_mets <- Parameters_table_night$remove_sublb_perc[Parameters_table_night$Compartment == "Ileum"]

        perc_to_remove_ileum_to_cecum <- as.numeric(Parameters_table_night$remove_from_ileum[Parameters_table_night$Compartment == "Cecum"])
        fraction_to_cecum <- perc_to_remove_ileum_to_cecum / perc_to_remove_ileum_bac
        perc_to_remove_ileum_to_colon <- as.numeric(Parameters_table_night$remove_from_ileum[Parameters_table_night$Compartment == "Colon"])
        fraction_to_colon <- perc_to_remove_ileum_to_colon / perc_to_remove_ileum_bac

        removefunction_result_ileum <- RemoveBacMetabolites(sim_ileum, arena2_object = arena2_ileum,
                                                            perc_to_remove_bac = perc_to_remove_ileum_bac,
                                                            perc_to_remove_mets = perc_to_remove_ileum_mets,
                                                            type_specie_table = type_specie_table_ileum)
        arena2_ileum <- removefunction_result_ileum[[1]]

        # split removal results between cecum and colon
        split_result <- split_removal_to_cecum_colon(removefunction_result_ileum, fraction_to_cecum, fraction_to_colon)
        removefunction_result_ileum_to_cecum <- split_result$to_cecum
        removefunction_result_ileum_to_colon <- split_result$to_colon

        # add bacs and mets from jejunum:
        arena2_ileum <- AddBacMets_fromUpstreamCompartment(removefunction_result_n0 = removefunction_result_jejunum,
                                                           sim_object_n1 = sim_ileum, arena2_n1 = arena2_ileum,
                                                           model_list_n1 = model_lists$Ileum,
                                                           model_list_n0 = model_lists$Jejunum)

        ### REFLUX (transfer part of ileal arena up to jejunum,
        # in the 'colon' iteration, part of colonic arena will be transferred up to ileum)
      } else if (timestep %% as.numeric(Parameters_table_night$reflux_freq_hrs[Parameters_table_night$Compartment == "Ileum"]) == 1) {
        perc_to_reflux_ileum <- as.numeric(Parameters_table_night$reflux_perc[Parameters_table_night$Compartment == "Ileum"])
        compartment_arenas[[3]] <- Reflux(sim_object_n1 = sim_ileum, perc_to_remove_bac_n1 = perc_to_reflux_ileum,
                                          perc_to_remove_mets_n1 = perc_to_reflux_ileum,
                                          type_specie_table_n1 = type_specie_table_ileum,
                                          sim_object_n0 = sim_jejunum, arena2_n0 = arena2_jejunum,
                                          model_list_n0 = model_lists$Jejunum, model_list_n1 = model_lists$Ileum)[[1]]

        arena2_ileum <- Reflux(sim_object_n1 = sim_ileum, perc_to_remove_bac_n1 = perc_to_reflux_ileum,
                               perc_to_remove_mets_n1 = perc_to_reflux_ileum,
                               type_specie_table_n1 = type_specie_table_ileum,
                               sim_object_n0 = sim_jejunum, arena2_n0 = arena2_jejunum,
                               model_list_n0 = model_lists$Jejunum, model_list_n1 = model_lists$Ileum)[[2]]
      }

      # update ileal arena
      compartment_arenas[[4]] <- arena2_ileum

      ####################################
      ############  CECUM  ###############
      ####################################
      ### simulate cecum for 1 h
      sim_cecum <- BacArena::simEnv(compartment_arenas[[5]], time = 1)
      out_lists[[5]] <- save_sim_outputs(sim_cecum, out_lists[[5]], timestep)

      # get type_specie_table
      if (timestep == 1) {
        type_specie_table_cecum <- CreateTypeSpecieAssoc_table(sim_cecum)
      }

      # select pre-set % of the grid randomly that will be removed once in 4 hours (released to colon)
      if (timestep %% 4 == 0) {
        perc_to_remove_cecum_bac <- Parameters_table_night$remove_bac_perc[Parameters_table_night$Compartment == "Cecum"]
        perc_to_remove_cecum_mets <- Parameters_table_night$remove_sublb_perc[Parameters_table_night$Compartment == "Cecum"]
        removefunction_result_cecum <- RemoveBacMetabolites(sim_cecum,
                                                            perc_to_remove_bac = perc_to_remove_cecum_bac,
                                                            perc_to_remove_mets = perc_to_remove_cecum_mets,
                                                            type_specie_table = type_specie_table_cecum)
        arena2_cecum <- removefunction_result_cecum[[1]]
      } else {
        arena2_cecum <- BacArena::getArena(sim_cecum, 1)
      }

      ## NO REFLUX - stuff transferred from ileum downstream
      if (!(timestep %% as.numeric(Parameters_table_night$reflux_freq_hrs[Parameters_table_night$Compartment == "Ileum"]) == 1)) {
        # add metabolites and bacteria removed from the ileum
        arena2_cecum <- AddBacMets_fromUpstreamCompartment(removefunction_result_n0 = removefunction_result_ileum_to_cecum,
                                                           sim_object_n1 = sim_cecum, arena2_n1 = arena2_cecum,
                                                           model_list_n1 = model_lists$Ceca,
                                                           model_list_n0 = model_lists$Ileum)
      }
      ## IF REFLUX - then nothing is transferred from ileum (fraction of colon will be transferred upstream in 'colon' section)

      # if any new species happen to be in cecum - add them to the type_specie_table:
      type_specie_table_cecum <- update_type_specie_table(arena2_cecum, type_specie_table_cecum)

      # remove obligate aerobes (if any happened to be transferred):
      arena2_cecum <- remove_organisms_from_arena(arena2_cecum, obligate_aerobes, type_specie_table_cecum, "cecum")

      # update cecum arena
      compartment_arenas[[5]] <- arena2_cecum

      ####################################
      ############  COLON  ###############
      ####################################
      ### simulate colon for 1 h
      sim_colon <- BacArena::simEnv(compartment_arenas[[6]], time = 1)
      out_lists[[6]] <- save_sim_outputs(sim_colon, out_lists[[6]], timestep)

      # get type_specie_table
      if (timestep == 1) {
        type_specie_table_colon <- CreateTypeSpecieAssoc_table(sim_colon)
      }

      ## during the night nothing gets removed from colon

      # once in 4 hrs cecal contents are released to the colon:
      if (timestep %% 4 == 0) {
        arena2_colon <- AddBacMets_fromUpstreamCompartment(removefunction_result_n0 = removefunction_result_cecum,
                                                           sim_object_n1 = sim_colon, arena2_n1 = arena2_colon,
                                                           model_list_n1 = model_lists$Colon,
                                                           model_list_n0 = model_lists$Ceca)
      }

      ## NO REFLUX
      if (!(timestep %% as.numeric(Parameters_table_night$reflux_freq_hrs[Parameters_table_night$Compartment == "Ileum"]) == 1)) {
        # add metabolites and bacteria removed from the ileum
        arena2_colon <- AddBacMets_fromUpstreamCompartment(removefunction_result_n0 = removefunction_result_ileum_to_colon,
                                                           sim_object_n1 = sim_colon, arena2_n1 = arena2_colon,
                                                           model_list_n1 = model_lists$Colon,
                                                           model_list_n0 = model_lists$Ileum)
      } else {
        ## IF REFLUX, then nothing is transferred from ileum downstream,
        # however, fraction of colon is being transferred upstream to cecum and ileum:

        ### select % of the grid randomly that will be removed
        ### split between cecum and ileum
        perc_to_remove_colon_bac <- as.numeric(Parameters_table_night$reflux_perc[Parameters_table_night$Compartment == "Colon"])
        perc_to_remove_colon_mets <- as.numeric(Parameters_table_night$reflux_perc[Parameters_table_night$Compartment == "Colon"])

        perc_to_remove_colon_to_ileum <- as.numeric(Parameters_table_night$reflux_remove_from_colon[Parameters_table_night$Compartment == "Ileum"])
        fraction_to_ileum <- perc_to_remove_colon_to_ileum / perc_to_remove_colon_bac
        perc_to_remove_colon_to_cecum <- as.numeric(Parameters_table_night$reflux_remove_from_colon[Parameters_table_night$Compartment == "Cecum"])
        fraction_to_cecum <- perc_to_remove_colon_to_cecum / perc_to_remove_colon_bac

        removefunction_result_colon <- RemoveBacMetabolites(sim_colon, arena2_object = arena2_colon,
                                                            perc_to_remove_bac = perc_to_remove_colon_bac,
                                                            perc_to_remove_mets = perc_to_remove_colon_mets,
                                                            type_specie_table = type_specie_table_colon)
        arena2_colon <- removefunction_result_colon[[1]]

        # get the part that will be transferred to cecum and to ileum:
        removefunction_result_colon_to_ileum <- removefunction_result_colon
        removefunction_result_colon_to_cecum <- removefunction_result_colon
        # sample fraction_to_ileum of orgdat removed from colon:
        removefunction_result_colon_to_ileum[[2]] <- removefunction_result_colon_to_ileum[[2]][sample(nrow(removefunction_result_colon_to_ileum[[2]]),
                                                                                                      ceiling(fraction_to_ileum * nrow(removefunction_result_colon_to_ileum[[2]]))), ]
        # remaining goes to cecum:
        removefunction_result_colon_to_cecum[[2]] <- dplyr::setdiff(removefunction_result_colon[[2]], removefunction_result_colon_to_ileum[[2]])

        # do same for sublb (metabolites) table:
        removefunction_result_colon_to_ileum[[3]] <- removefunction_result_colon_to_ileum[[3]][sample(nrow(removefunction_result_colon_to_ileum[[3]]),
                                                                                                      round(fraction_to_ileum * nrow(removefunction_result_colon_to_ileum[[3]]))), ]
        removefunction_result_colon_to_cecum[[3]] <- dplyr::setdiff(removefunction_result_colon[[3]], removefunction_result_colon_to_ileum[[3]])

        # update slot [[4]] - species_rm table
        type_number_colon_to_ileum <- table(removefunction_result_colon_to_ileum[[2]]$type)
        type_number_colon_to_cecum <- table(removefunction_result_colon_to_cecum[[2]]$type)
        removefunction_result_colon_to_ileum[[4]] <- subset(removefunction_result_colon_to_ileum[[4]],
                                                            removefunction_result_colon_to_ileum[[4]]$Type %in% names(type_number_colon_to_ileum))
        removefunction_result_colon_to_cecum[[4]] <- subset(removefunction_result_colon_to_cecum[[4]],
                                                            removefunction_result_colon_to_cecum[[4]]$Type %in% names(type_number_colon_to_cecum))
        removefunction_result_colon_to_ileum[[4]]$Number <- type_number_colon_to_ileum
        removefunction_result_colon_to_cecum[[4]]$Number <- type_number_colon_to_cecum
      }

      # if any new species happen to be in colon - add them to the type_specie_table:
      type_specie_table_colon <- update_type_specie_table(arena2_colon, type_specie_table_colon)

      # remove obligate aerobes (if any happened to be transferred):
      arena2_colon <- remove_organisms_from_arena(arena2_colon, obligate_aerobes, type_specie_table_colon, "colon")

      # update colon arena
      compartment_arenas[[6]] <- arena2_colon

      #####################################################
      ############   ILEUM AND CECUM AGAIN  ###############
      ##########  TO ADD REFLUXES FROM COLON  #############
      #####################################################
      ### REFLUX (transfer part of colonic arena up to ileum and cecum)
      if (timestep %% as.numeric(Parameters_table_night$reflux_freq_hrs[Parameters_table_night$Compartment == "Colon"]) == 1) {
        # add metabolites and bacteria removed from the colon

        compartment_arenas[[4]] <- AddBacMets_fromDownstreamCompartment(removefunction_result_colon_to_ileum,
                                                                        sim_object_n0 = sim_ileum,
                                                                        arena2_n0 = compartment_arenas[[4]],
                                                                        model_list_n0 = model_lists$Ileum,
                                                                        model_list_n1 = model_lists$Colon)

        # imitate absorption of sugars, amino acids and fats; and conversion of starches in ileum:
        compartment_arenas[[4]] <- AbsorptionConversion_ileum(arena2 = compartment_arenas[[4]], diet)

        compartment_arenas[[5]] <- AddBacMets_fromDownstreamCompartment(removefunction_result_colon_to_cecum,
                                                                        sim_object_n0 = sim_cecum,
                                                                        arena2_n0 = compartment_arenas[[5]],
                                                                        model_list_n0 = model_lists$Ceca,
                                                                        model_list_n1 = model_lists$Colon)
      }

    }
  }

  # ===== Compile output =====
  if (is.null(AminoAcid_X)) { AminoAcid_X <- "None" }
  if (is.null(Sugars_X)) { Sugars_X <- "None" }

  simlist <- c(
    out_lists[1:6],
    list(
      seedN,
      Parameters_table,
      Parameters_table_night,
      ExchangeConstraints,
      list(econc_proportion = econc_proportion, prop_gizzard = prop_gizzard,
           Sugars_X = Sugars_X, AminoAcids_X = AminoAcid_X,
           timesteps = timesteps, sample = sample, diet_name = diet_name,
           compound = compound, compound_gper100g = compound_gper100g)
    )
  )
  return(simlist)
})))

stopCluster(cl)

# ===== Save output =====
now <- Sys.time()
today <- strsplit(as.character(as.POSIXlt(now)[1]), "\\ ")[[1]][1]
stir_TF <- Parameters_table$stir[1]

output_dir <- file.path(base_dir, "arena_sim_objects/CoupledSixComps_D10/WithPrebiotics/2024/gapseq_models")

saveRDS(simlist, file.path(output_dir,
                           paste0(sample, "_", as.character(diet_name),
                                  "_diet_grower_NEW_AFTER48h_", as.character(timesteps + 48),
                                  "h_econc", econc_proportion,
                                  "_propGizzard", as.character(prop_gizzard),
                                  "_CobaltX12_rmAnaero_less02_mucUrea_wCoA_moreNightFlow_",
                                  compound, "_", as.character(compound_gper100g), "gper100g_",
                                  today,
                                  "_gapseq.rds")))
