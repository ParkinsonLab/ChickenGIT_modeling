# =============================================================================
# Simulation of 6-compartment chicken GIT model (48-96h continuation, control)
# =============================================================================
# Continues a BacArena simulation from hour 48 to 96 for the control
# (no supplements) condition. Loads the saved 48h output, recreates
# arenas from the final state, and runs additional timesteps with
# the same day/night cycling, peristalsis, and reflux logic.
#
# USAGE:   Rscript ...continueAfter48h.R <sample> <timesteps> <pattern> <sim_date>
#
# INPUTS:  Previously saved 48h output .rds file, plus same input files as 0-48h script
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

# ---- Base directory (adjust to your environment) ----
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


#### ===== Load input files ================================
diet_corn  <- read.table(file.path(base_dir, "corn_diet_grower_aggregated_2023.txt"),
                         sep = '\t', header = TRUE)
diet_wheat <- read.table(file.path(base_dir, "wheat_diet_grower_aggregated_2023.txt"),
                         sep = '\t', header = TRUE)
OTU_table  <- read.table(file.path(base_dir, "OTUtable.tsv"))
metadata   <- read.table(file.path(base_dir, "metadata.txt"), header = TRUE)

Parameters_table <- readxl::read_xlsx(file.path(base_dir, "ModelGeneration_Files/Table_parameters_coupling_24h_2024.xlsx"),
                                      sheet = "Table_parameters_16h_2023")
Parameters_table_night <- readxl::read_xlsx(file.path(base_dir, "ModelGeneration_Files/Table_parameters_coupling_24h_2024.xlsx"),
                                            sheet = "Table_parameters_nightmode")

Metabolite_conc_h24 <- read.delim(file.path(base_dir, "ModelGeneration_Files/MetaboliteConcentrations_h24_10Samples_growerdiet_econc1_propGizzard100_CobaltX12_rmAnaero_less02_mucUreaCecaColon_ordered_fromsim20240701.txt"),
                                  sep = '\t', header = TRUE)

ExchangeConstraints <- read.table(file.path(base_dir, "ModelGeneration_Files/Compounds_with_uptake_lbnds_gapseq.txt"),
                                  sep = '\t', header = TRUE)

file_path <- base_dir

# Load all SBML models for all GIT compartments (combined model list)
model_list_all_compartments <-
  readRDS(file.path(base_dir, "ModelGeneration_Files/ChickenAllSpecies_noAGP_D10_gapseq_ModelList_MAGs_gapseqCobraAdapt_Aug2023.rds"))

# Load compartment-specific SBML model lists
model_files <- list.files(file.path(file_path, "ModelGeneration_Files"),
                          pattern = "Chicken_.*_ModelList_MAGs_gapseqCobraAdapt_Aug2023.rds",
                          full.names = TRUE)
# rearrange in the GIT order
model_files <- model_files[c(4, 3, 6, 5, 1, 2)]
model_lists <- lapply(model_files, readRDS)
names(model_lists) <- c("Gizzard", "Duodenum", "Jejunum", "Ileum", "Ceca", "Colon")

# Load names of organisms in each compartment and corresponding taxonomic category
organism_files <- list.files(file.path(file_path, "ModelGeneration_Files"),
                             pattern = "Chicken_.*_organisms_tax-category_table_gapseq_MAGs_gapseqCobraAdapt_Aug2023.rds",
                             full.names = TRUE)
# rearrange in the GIT order
organism_files <- organism_files[c(4, 3, 6, 5, 1, 2)]
organism_lists <- lapply(organism_files, readRDS)

# Load all function scripts
functions_scripts <- list.files(file.path(file_path, "ModelGeneration_scripts"),
                                pattern = "Function.*.R", full.names = TRUE)


# Select diet based on sample name
diet_name <- tolower(str_remove(sample, '[0-9]{1,2}-[0-9]'))
if (diet_name == "wheat") {
  diet <- diet_wheat
} else {
  diet <- diet_corn
}

### Load corresponding output file after 48h of simulation:
year <- strsplit(sim_date, '-')[[1]][1]
out_path <- file.path(base_dir, "arena_sim_objects/CoupledSixComps_D10/NoPathogens", year, "gapseq_models")
outfile <- readRDS(file.path(out_path,
                             paste0(sample, "_", diet_name, "_",
                                    pattern, "_", sim_date, "_gapseq.rds")))

# Extract simulation parameters from the saved 48h output
econc_proportion <- outfile[[1]][[11]]$econc_proportion
prop_gizzard     <- outfile[[1]][[11]]$prop_gizzard
Sugars_X         <- outfile[[1]][[11]]$Sugars_X
AminoAcid_X      <- outfile[[1]][[11]]$AminoAcids_X


##########################################################
############ PARALLEL SIMULATION ON CLUSTER ##############
############# RUN MODELS WITH REPLICATES #################
##########################################################

# Initialize cluster
cores <- 5
replicates <- length(outfile)

now   <- Sys.time()
today <- strsplit(as.character(as.POSIXlt(now)[1]), "\\ ")[[1]][1]

cl <- makeCluster(cores, type = "PSOCK",
                  outfile = file.path(cluster_log_dir,
                                      paste0("cluster_coupling_6comps_10d_continueAfter48h_GAPSEQ_",
                                             sample, "_growerdiet_",
                                             diet_name, "_propEconc", econc_proportion,
                                             "_propGizzard_mucusUrea_CoAupd_", prop_gizzard, "_",
                                             today,
                                             "_uptakeConstraints_sugarsX", Sugars_X,
                                             "_AminoAcidsX", AminoAcid_X, ".log")))

clusterExport(cl, c("functions_scripts"))
# Load all function scripts on each worker
clusterCall(cl, function() { lapply(functions_scripts, source) })

# Export all parameters and input files to workers
clusterExport(cl, c("prop_gizzard", "Parameters_table", "Parameters_table_night",
                     "econc_proportion", "ExchangeConstraints"))
clusterExport(cl, list("organism_lists", "model_lists", "model_list_all_compartments",
                        "diet", "diet_name", "diet_corn", "diet_wheat", "sample",
                        "Metabolite_conc_h24"))
clusterExport(cl, c("OTU_table", "metadata"))
clusterExport(cl, c("timesteps", "Sugars_X", "AminoAcid_X", "outfile"))


print(system.time(simlist <- parLapply(cl, 1:replicates, function(i) {

  # ---------------------------------------------------------------------------
  # Helper: save_sim_outputs
  # Extracts and stores population, substance, orgdat, flux, and shadow price

  # data from a BacArena simulation object into the output list.
  # ---------------------------------------------------------------------------
  save_sim_outputs <- function(sim_object, out_list, comp_index, timestep) {
    curves <- BacArena::plotCurves(sim_object, retdata = TRUE, graph = FALSE)
    out_list[[comp_index]][[timestep]]$Population <- curves$Population
    out_list[[comp_index]][[timestep]]$Substances <- curves$Substances
    out_list[[comp_index]][[timestep]]$Orgdat     <- sim_object@simlist
    out_list[[comp_index]][[timestep]]$Fluxlist   <- sim_object@mfluxlist
    # Save only negative shadow prices (limiting metabolites)
    shadowlist <- reshape2::melt(unlist(lapply(sim_object@shadowlist[[2]], function(x) x[x < -0.1])))
    shadowlist$Metabolite <- unlist(lapply(rownames(shadowlist), function(x) { strsplit(x, "EX_")[[1]][2] }))
    shadowlist$Specie     <- unlist(lapply(rownames(shadowlist), function(x) { gsub("\\.$", "", strsplit(x, "EX_")[[1]][1]) }))
    out_list[[comp_index]][[timestep]]$Shadowlist <- shadowlist
    return(out_list)
  }

  # ---------------------------------------------------------------------------
  # Helper: update_type_specie_table
  # Adds any new species types that appeared in an arena to the type-specie
  # lookup table (needed after peristalsis/reflux transfers).
  # ---------------------------------------------------------------------------
  update_type_specie_table <- function(arena, type_specie_table) {
    new_types <- setdiff(unique(arena@orgdat$type), type_specie_table$Type)
    if (length(new_types) > 0) {
      for (bactype in new_types) {
        type_specie_new <- data.frame(
          Name = arena@specs[[bactype]]@model@mod_desc,
          Type = bactype,
          stringsAsFactors = FALSE
        )
        type_specie_table <- rbind(type_specie_table, type_specie_new)
      }
    }
    return(type_specie_table)
  }

  # ---------------------------------------------------------------------------
  # Helper: remove_organisms_from_arena
  # Removes organisms of given names from an arena's orgdat, using the
  # type-specie lookup table to find matching types.
  # ---------------------------------------------------------------------------
  remove_organisms_from_arena <- function(arena, organism_names, type_specie_table, compartment_label = "") {
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

  # ---------------------------------------------------------------------------
  # Helper: split_removal_to_cecum_colon
  # Splits the RemoveBacMetabolites result from ileum into two portions:
  # one destined for cecum and one for colon, based on fraction parameters.
  # ---------------------------------------------------------------------------
  split_removal_to_cecum_colon <- function(removefunction_result, fraction_to_cecum, fraction_to_colon) {
    result_to_cecum <- removefunction_result
    result_to_colon <- removefunction_result

    # Split orgdat (bacteria)
    result_to_cecum[[2]] <- result_to_cecum[[2]][sample(nrow(result_to_cecum[[2]]),
                                                         ceiling(fraction_to_cecum * nrow(result_to_cecum[[2]]))), ]
    result_to_colon[[2]] <- dplyr::setdiff(removefunction_result[[2]], result_to_cecum[[2]])

    # Split sublb (metabolites)
    result_to_cecum[[3]] <- result_to_cecum[[3]][sample(nrow(result_to_cecum[[3]]),
                                                         round(fraction_to_colon * nrow(result_to_cecum[[3]]))), ]
    result_to_colon[[3]] <- dplyr::setdiff(removefunction_result[[3]], result_to_cecum[[3]])

    # Update species_rm slot [[4]] with correct types and counts
    type_number_to_cecum <- table(result_to_cecum[[2]]$type)
    type_number_to_colon <- table(result_to_colon[[2]]$type)
    result_to_cecum[[4]] <- subset(result_to_cecum[[4]], result_to_cecum[[4]]$Type %in% names(type_number_to_cecum))
    result_to_colon[[4]] <- subset(result_to_colon[[4]], result_to_colon[[4]]$Type %in% names(type_number_to_colon))
    result_to_cecum[[4]]$Number <- type_number_to_cecum
    result_to_colon[[4]]$Number <- type_number_to_colon

    return(list(to_cecum = result_to_cecum, to_colon = result_to_colon))
  }

  # ---------------------------------------------------------------------------
  # Recover seed from the original 48h run (ensures reproducibility)
  # ---------------------------------------------------------------------------
  seedN <- outfile[[i]][[7]]

  # ---------------------------------------------------------------------------
  # Reflux function: reverse peristalsis from compartment N+1 back to N
  # ---------------------------------------------------------------------------
  Reflux <- function(sim_object_n1, perc_to_remove_bac_n1, perc_to_remove_mets_n1 = 0,
                     type_specie_table_n1,
                     sim_object_n0, arena2_n0, model_list_n0, model_list_n1) {
    # if perc_to_remove_mets_n1 not specified, set it equal to perc_to_remove_bac
    if (perc_to_remove_mets_n1 == 0) {
      perc_to_remove_mets_n1 <- perc_to_remove_bac_n1
    }

    # REVERSE PERISTALSIS: COMP N+1 - remove; COMP N - add
    # Update arena of compartment N+1
    removefunction_result_n1 <- RemoveBacMetabolites(sim_object = sim_object_n1,
                                                      perc_to_remove_bac = perc_to_remove_bac_n1,
                                                      perc_to_remove_mets = perc_to_remove_mets_n1,
                                                      type_specie_table = type_specie_table_n1)
    arena2_n1_upd <- removefunction_result_n1[[1]]

    # Then add retrieved bacs and mets to the compartment N:
    arena2_n0_upd <- AddBacMets_fromDownstreamCompartment(removefunction_result_n1,
                                                           sim_object_n0, arena2_n0,
                                                           model_list_n0, model_list_n1)
    return(list("arena2_n0_upd" = arena2_n0_upd, "arena2_n1_upd" = arena2_n1_upd))
  }

  # ---------------------------------------------------------------------------
  # set_pH: set proton concentration in arena based on pH value
  # ---------------------------------------------------------------------------
  set_pH <- function(arena, pH) {
    arena <- BacArena::addSubs(arena, smax = 1000 / (10^pH),
                               mediac = "EX_cpd00067_e0", unit = "mM", add = FALSE)
    return(arena)
  }

  # ---------------------------------------------------------------------------
  # set_o2: set dissolved oxygen concentration
  # ---------------------------------------------------------------------------
  set_o2 <- function(arena, max_conc) {
    arena <- BacArena::addSubs(arena, smax = max_conc,
                               mediac = "EX_cpd00007_e0", unit = "mM", add = FALSE)
    return(arena)
  }

  # Handle "None" string for optional multipliers
  if (Sugars_X == "None") {
    Sugars_X <- NULL
  }
  if (AminoAcid_X == "None") {
    AminoAcid_X <- NULL
  }

  # Update exchange reactions constraints based on the input constraints table
  model_lists <- lapply(model_lists, ChangeUptakeConstraints,
                        constraints_table = ExchangeConstraints,
                        Sugars_X = Sugars_X, AA_X = AminoAcid_X)

  stir_TF            <- Parameters_table$stir[1]
  rich_TF            <- FALSE
  diet_added_perc    <- as.numeric(Parameters_table$diet_added_perc[Parameters_table$Compartment == "Gizzard"])
  diet_concentration <- as.numeric(Parameters_table$diet_proportion[Parameters_table$Compartment == "Gizzard"])

  # Strict anaerobes: not expected to survive in high-oxygen environments like gizzard
  strict_anaerobes <- c("Clostridium_saudiense_cobrapy_adapted",
                        "Lachnoclostridium_phocaeensis_cobrapy_adapted",
                        "Blautia_sp_cobrapy_adapted",
                        "Mediterraneibacter_sp002161355_cobrapy_adapted",
                        "Faecalibacterium_sp002160895_cobrapy_adapted",
                        "Oscillospiraceae_genusUBA9475_cobrapy_adapted",
                        "Flavonifractor_plautii-gapseqAdRm")

  # Obligate aerobes: not expected to survive in nearly anaerobic cecum and colon
  obligate_aerobes <- c("Acinetobacter_lwoffii-gapseqAdRm")

  # ---------------------------------------------------------------------------
  # RecreateModelArena: rebuild an arena from the hour-48 state
  #
  # Reads population abundances and metabolite concentrations from the saved
  # 48h output, then reconstructs the BacArena arena for a given compartment.
  # Also adds nucleobases/CoA and mucins/urea for cecum/colon.
  # ---------------------------------------------------------------------------
  RecreateModelArena <- function(compartment_number, sim_outfile, hour = 48,
                                 seedN, stir_var, tstep = 1, arena_size = 100) {
    Compartments_list <- list("Gizzard", "Duodenum", "Jejunum", "Ileum", "Cecum", "Colon")

    # Get the abundances after 48h
    AbundanceTable <- as.data.frame(sim_outfile[[compartment_number]][[hour]]$Population)
    AbundanceTable$Specie <- gsub("_cobrapy_adapted|-gapseqAdRm", "", rownames(AbundanceTable))

    # For gizzard: collect essential metabolites
    essential_metabolites <- c()

    # Initialize arena with same seed as original run
    arena <- BacArena::Arena(n = arena_size, m = arena_size, stir = stir_var,
                             Lx = 0.025 * (arena_size / 100), Ly = 0.025 * (arena_size / 100),
                             tstep = tstep, seed = seedN)

    # Add species in the abundances they were present after 48h
    # Go through all species present in GIT:
    for (m in 1:length(model_list_all_compartments)) {
      if (model_list_all_compartments[[m]]@mod_name %in% AbundanceTable$Specie) {
        if (AbundanceTable$V2[AbundanceTable$Specie == model_list_all_compartments[[m]]@mod_name] > 0) {

          bac <- BacArena::Bac(model = model_list_all_compartments[[m]], growtype = "exponential")
          arena <- BacArena::addOrg(arena, bac,
                                    amount = AbundanceTable$V2[AbundanceTable$Specie == model_list_all_compartments[[m]]@mod_name])

          # Only for gizzard: get list of essential metabolites to be added along with the diet
          if (compartment_number == 1) {
            essential_metabolites <- unique(c(essential_metabolites,
                                             BacArena::addEssentialMed(arena, bac, only_return = TRUE)))
          }
        }
      }
    }

    # Add metabolites in the concentrations at 48h
    ConcTable <- as.data.frame(sim_outfile[[compartment_number]][[hour]]$Substances)
    ConcTable$Metabolite <- rownames(ConcTable)
    ConcTable <- ConcTable[which(ConcTable$V2 > 0), ]

    # Remove biomass
    ConcTable <- ConcTable[-which(ConcTable$Metabolite == "EX_cpd11416_c0"), ]

    # Convert mmol/cell to mM to be added to arena
    ConcTable$Conc_mM <- ConcTable$V2 / 625
    arena <- BacArena::addSubs(arena, smax = ConcTable$Conc_mM, unit = "mM",
                               mediac = ConcTable$Metabolite, add = FALSE, addAnyway = TRUE)

    # Add nucleobases and essential CoA:
    # gapseq: cytosine, uracil, cytidine, adenosine, guanosine, CoA
    arena <- BacArena::addSubs(arena, smax = 1e-04, unit = "mM",
                               mediac = c("EX_cpd00307_e0", "EX_cpd00092_e0", "EX_cpd00367_e0",
                                          "EX_cpd00182_e0", "EX_cpd00311_e0", "EX_cpd00010_e0"),
                               add = TRUE, addAnyway = TRUE)

    # Add mucins and urea for cecum and colon
    if (compartment_number %in% c(5, 6)) {
      mucins_urea <- c("EX_cpd00122_e0", "EX_cpd00232_e0", "EX_cpd00832_e0", "EX_cpd02992_e0",
                        "EX_cpd11842_e0", "EX_cpd21520_e0", "EX_cpd00073_e0")
      arena <- BacArena::addSubs(arena, smax = 1e-05, unit = "mM",
                                 mediac = mucins_urea, add = TRUE, addAnyway = TRUE)
    }

    return(list(arena, essential_metabolites))
  }

  # ---------------------------------------------------------------------------
  # Recreate arenas from hour-48 state for all 6 compartments
  # ---------------------------------------------------------------------------
  compartment_arenas <- lapply(1:6, function(j) {
    RecreateModelArena(compartment_number = j, sim_outfile = outfile[[i]], hour = 48,
                       seedN = seedN, stir_var = stir_TF)[[1]]
  })

  # Get essential metabolites for gizzard (filtered to exclude those already in diet)
  essential_mets_gizzard <- RecreateModelArena(compartment_number = 1, sim_outfile = outfile[[i]],
                                               seedN = seedN, stir_var = stir_TF, tstep = 1)[[2]]
  essential_mets_gizzard <- essential_mets_gizzard[which(!(essential_mets_gizzard %in% diet$met_id))]

  # Initialize output lists
  out_lists <- lapply(1:6, function(x) vector("list", length = timesteps + 1))

  # ===========================================================================
  # MAIN SIMULATION LOOP
  # ===========================================================================
  for (timestep in 1:timesteps) {

    names(compartment_arenas) <- c("Gizzard", "Duodenum", "Jejunum", "Ileum", "Cecum", "Colon")

    # =========================================================================
    # DAY TIME (hours 1-16, 25-40, 49-64 of each 24h cycle)
    # =========================================================================
    if (timestep %in% c(1:16, 25:40, 49:64)) {

      # Set pH and oxygen for all compartments
      for (j in 1:6) {
        compartment_arenas[[j]] <- set_pH(compartment_arenas[[j]],
                                          Parameters_table$pH[Parameters_table$Compartment == names(compartment_arenas)[j]])
        compartment_arenas[[j]] <- set_o2(compartment_arenas[[j]],
                                          max_conc = as.numeric(Parameters_table$oxygen_mM[Parameters_table$Compartment == names(compartment_arenas)[j]]))
      }


      ####################################
      ############## GIZZARD #############
      ####################################
      sim_gizzard <- BacArena::simEnv(compartment_arenas[[1]], time = 1)
      out_lists <- save_sim_outputs(sim_gizzard, out_lists, 1, timestep)

      arena2_gizzard <- BacArena::getArena(sim_gizzard, 1)

      # Build type-specie lookup table on first timestep
      if (timestep == 1) {
        type_specie_table_gizzard <- CreateTypeSpecieAssoc_table(sim_gizzard)
      }

      reflux_duo_frequency <- as.numeric(as.character(Parameters_table$reflux_freq_hrs[Parameters_table$Compartment == "Duodenum"]))

      # Once every N hrs add fraction of diet (incl. essential mets) + initial composition to the grid
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

      # Peristalsis: remove fraction from gizzard and transfer to duodenum
      # (unless reflux is happening this hour)
      if (!(timestep %% reflux_duo_frequency == 0)) {
        perc_to_remove_gizzard_bac  <- Parameters_table$remove_bac_perc[Parameters_table$Compartment == "Gizzard"]
        perc_to_remove_gizzard_mets <- Parameters_table$remove_sublb_perc[Parameters_table$Compartment == "Gizzard"]

        removefunction_result_gizzard <- RemoveBacMetabolites(sim_gizzard, arena2_object = arena2_gizzard,
                                                              perc_to_remove_bac = perc_to_remove_gizzard_bac,
                                                              perc_to_remove_mets = perc_to_remove_gizzard_mets,
                                                              type_specie_table = type_specie_table_gizzard)
        arena2_gizzard <- removefunction_result_gizzard[[1]]
      }

      compartment_arenas[[1]] <- arena2_gizzard


      ####################################
      ############# DUODENUM #############
      ####################################
      if (timestep > 1) {
        # Imitate absorption of sugars and conversion of starches in duodenum
        compartment_arenas[[2]] <- AbsorptionConversion_duodenum(arena2 = compartment_arenas[[2]])
      }

      sim_duodenum <- BacArena::simEnv(compartment_arenas[[2]], time = 1)
      out_lists <- save_sim_outputs(sim_duodenum, out_lists, 2, timestep)

      if (timestep == 1) {
        type_specie_table_duodenum <- CreateTypeSpecieAssoc_table(sim_duodenum)
      }

      # NO REFLUX: transfer from gizzard to duodenum, remove fraction from duodenum
      reflux_duo <- FALSE
      if (!(timestep %% reflux_duo_frequency == 0)) {
        perc_to_remove_duodenum_bac  <- Parameters_table$remove_bac_perc[Parameters_table$Compartment == "Duodenum"]
        perc_to_remove_duodenum_mets <- Parameters_table$remove_sublb_perc[Parameters_table$Compartment == "Duodenum"]

        arena2_duodenum <- getArena(sim_duodenum, 1)

        removefunction_result_duodenum <- RemoveBacMetabolites(sim_duodenum, arena2_object = arena2_duodenum,
                                                               perc_to_remove_bac = perc_to_remove_duodenum_bac,
                                                               perc_to_remove_mets = perc_to_remove_duodenum_mets,
                                                               type_specie_table = type_specie_table_duodenum)
        arena2_duodenum <- removefunction_result_duodenum[[1]]

        # Add mets and bacs that were removed from gizzard
        arena2_duodenum <- AddBacMets_fromUpstreamCompartment(removefunction_result_n0 = removefunction_result_gizzard,
                                                              sim_object_n1 = sim_duodenum, arena2_n1 = arena2_duodenum,
                                                              model_list_n1 = model_lists$Duodenum,
                                                              model_list_n0 = model_lists$Gizzard)

      # REFLUX: reverse peristalsis duodenum -> gizzard
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

        # Update gizzard type-specie table with any new species from reflux
        type_specie_table_gizzard <- update_type_specie_table(compartment_arenas[[1]], type_specie_table_gizzard)

        # Remove strict anaerobes from gizzard (if transferred during reflux)
        compartment_arenas[[1]] <- remove_organisms_from_arena(compartment_arenas[[1]], strict_anaerobes,
                                                               type_specie_table_gizzard, "gizzard")
      }

      compartment_arenas[[2]] <- arena2_duodenum


      ####################################
      ############# JEJUNUM ##############
      ####################################
      if (timestep > 2) {
        # Imitate absorption of sugars, amino acids and fats in jejunum
        compartment_arenas[[3]] <- AbsorptionConversion_jejunum(arena2 = compartment_arenas[[3]], diet)
      }

      sim_jejunum <- BacArena::simEnv(compartment_arenas[[3]], time = 1)
      out_lists <- save_sim_outputs(sim_jejunum, out_lists, 3, timestep)

      if (timestep == 1) {
        type_specie_table_jejunum <- CreateTypeSpecieAssoc_table(sim_jejunum)
      }

      arena2_jejunum <- getArena(sim_jejunum, 1)

      perc_to_remove_jejunum_bac  <- Parameters_table$remove_bac_perc[Parameters_table$Compartment == "Jejunum"]
      perc_to_remove_jejunum_mets <- Parameters_table$remove_sublb_perc[Parameters_table$Compartment == "Jejunum"]
      removefunction_result_jejunum <- RemoveBacMetabolites(sim_jejunum, arena2_object = arena2_jejunum,
                                                            perc_to_remove_bac = perc_to_remove_jejunum_bac,
                                                            perc_to_remove_mets = perc_to_remove_jejunum_mets,
                                                            type_specie_table = type_specie_table_jejunum)
      arena2_jejunum <- removefunction_result_jejunum[[1]]

      # If no reflux happened, add bacs and mets from duodenum
      if (reflux_duo == FALSE) {
        arena2_jejunum <- AddBacMets_fromUpstreamCompartment(removefunction_result_n0 = removefunction_result_duodenum,
                                                             sim_object_n1 = sim_jejunum, arena2_n1 = arena2_jejunum,
                                                             model_list_n1 = model_lists$Jejunum,
                                                             model_list_n0 = model_lists$Duodenum)
      }

      compartment_arenas[[3]] <- arena2_jejunum


      ###################################
      ############## ILEUM ##############
      ###################################
      if (timestep > 2) {
        # Imitate absorption of sugars, amino acids and fats in ileum
        compartment_arenas[[4]] <- AbsorptionConversion_ileum(arena2 = compartment_arenas[[4]], diet)
      }

      sim_ileum <- BacArena::simEnv(compartment_arenas[[4]], time = 1)
      out_lists <- save_sim_outputs(sim_ileum, out_lists, 4, timestep)

      if (timestep == 1) {
        type_specie_table_ileum <- CreateTypeSpecieAssoc_table(sim_ileum)
      }

      arena2_ileum <- getArena(sim_ileum, 1)

      # Remove fraction from ileum: split between cecum and colon
      perc_to_remove_ileum_bac  <- Parameters_table$remove_bac_perc[Parameters_table$Compartment == "Ileum"]
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

      # Split removed contents between cecum and colon
      split_result <- split_removal_to_cecum_colon(removefunction_result_ileum, fraction_to_cecum, fraction_to_colon)
      removefunction_result_ileum_to_cecum <- split_result$to_cecum
      removefunction_result_ileum_to_colon <- split_result$to_colon

      # Add bacs and mets from jejunum
      arena2_ileum <- AddBacMets_fromUpstreamCompartment(removefunction_result_n0 = removefunction_result_jejunum,
                                                         sim_object_n1 = sim_ileum, arena2_n1 = arena2_ileum,
                                                         model_list_n1 = model_lists$Ileum,
                                                         model_list_n0 = model_lists$Jejunum)

      compartment_arenas[[4]] <- arena2_ileum


      ####################################
      ############  CECUM  ###############
      ####################################
      sim_cecum <- BacArena::simEnv(compartment_arenas[[5]], time = 1)
      out_lists <- save_sim_outputs(sim_cecum, out_lists, 5, timestep)

      if (timestep == 1) {
        type_specie_table_cecum <- CreateTypeSpecieAssoc_table(sim_cecum)
      }

      # Every 4 hours, release fraction of cecal contents to colon
      if (timestep %% 4 == 0) {
        perc_to_remove_cecum_bac  <- Parameters_table$remove_bac_perc[Parameters_table$Compartment == "Cecum"]
        perc_to_remove_cecum_mets <- Parameters_table$remove_sublb_perc[Parameters_table$Compartment == "Cecum"]
        removefunction_result_cecum <- RemoveBacMetabolites(sim_cecum,
                                                            perc_to_remove_bac = perc_to_remove_cecum_bac,
                                                            perc_to_remove_mets = perc_to_remove_cecum_mets,
                                                            type_specie_table = type_specie_table_cecum)
        arena2_cecum <- removefunction_result_cecum[[1]]
      } else {
        arena2_cecum <- BacArena::getArena(sim_cecum, 1)
      }

      # Add metabolites and bacteria removed from the ileum
      arena2_cecum <- AddBacMets_fromUpstreamCompartment(removefunction_result_n0 = removefunction_result_ileum_to_cecum,
                                                         sim_object_n1 = sim_cecum, arena2_n1 = arena2_cecum,
                                                         model_list_n1 = model_lists$Ceca,
                                                         model_list_n0 = model_lists$Ileum)

      if (timestep > 1) {
        type_specie_table_cecum <- update_type_specie_table(arena2_cecum, type_specie_table_cecum)
      }

      # Remove obligate aerobes from cecum
      arena2_cecum <- remove_organisms_from_arena(arena2_cecum, obligate_aerobes,
                                                  type_specie_table_cecum, "cecum")

      compartment_arenas[[5]] <- arena2_cecum


      ####################################
      ############  COLON  ###############
      ####################################
      sim_colon <- BacArena::simEnv(compartment_arenas[[6]], time = 1)
      out_lists <- save_sim_outputs(sim_colon, out_lists, 6, timestep)

      if (timestep == 1) {
        type_specie_table_colon <- CreateTypeSpecieAssoc_table(sim_colon)
      }

      # Remove fraction from colon (excreted)
      perc_to_remove_colon_bac  <- Parameters_table$remove_bac_perc[Parameters_table$Compartment == "Colon"]
      perc_to_remove_colon_mets <- Parameters_table$remove_sublb_perc[Parameters_table$Compartment == "Colon"]
      removefunction_result_colon <- RemoveBacMetabolites(sim_colon,
                                                          perc_to_remove_bac = perc_to_remove_colon_bac,
                                                          perc_to_remove_mets = perc_to_remove_colon_mets,
                                                          type_specie_table = type_specie_table_colon)
      arena2_colon <- removefunction_result_colon[[1]]

      # Once every 4 hrs cecal contents are released to the colon
      if (timestep %% 4 == 0) {
        arena2_colon <- AddBacMets_fromUpstreamCompartment(removefunction_result_n0 = removefunction_result_cecum,
                                                           sim_object_n1 = sim_colon, arena2_n1 = arena2_colon,
                                                           model_list_n1 = model_lists$Colon,
                                                           model_list_n0 = model_lists$Ceca)
      }
      # Add metabolites and bacteria removed from the ileum
      arena2_colon <- AddBacMets_fromUpstreamCompartment(removefunction_result_n0 = removefunction_result_ileum_to_colon,
                                                         sim_object_n1 = sim_colon, arena2_n1 = arena2_colon,
                                                         model_list_n1 = model_lists$Colon,
                                                         model_list_n0 = model_lists$Ileum)

      if (timestep > 1) {
        type_specie_table_colon <- update_type_specie_table(arena2_colon, type_specie_table_colon)
      }

      # Remove obligate aerobes from colon
      arena2_colon <- remove_organisms_from_arena(arena2_colon, obligate_aerobes,
                                                  type_specie_table_colon, "colon")

      # Absorb water in colon: retain 20%
      sublb_colon <- as.data.frame(BacArena::getSublb(arena2_colon))
      sublb_colon[, "EX_cpd00001_e0"] <- 0.2 * sublb_colon[1, "EX_cpd00001_e0"]
      arena2_colon@sublb <- as.matrix(sublb_colon)

      compartment_arenas[[6]] <- arena2_colon
    }


    # =========================================================================
    # NIGHT TIME (hours 17-24, 41-48, 65-72 of each 24h cycle)
    # =========================================================================
    if (timestep %in% c(17:24, 41:48, 65:72)) {

      # Set pH and oxygen for all compartments
      for (j in 1:6) {
        compartment_arenas[[j]] <- set_pH(compartment_arenas[[j]],
                                          Parameters_table$pH[Parameters_table$Compartment == names(compartment_arenas)[j]])
        compartment_arenas[[j]] <- set_o2(compartment_arenas[[j]],
                                          as.numeric(Parameters_table$oxygen_mM[Parameters_table$Compartment == names(compartment_arenas)[j]]))
      }


      ####################################
      ############## GIZZARD #############
      ####################################
      sim_gizzard <- BacArena::simEnv(compartment_arenas[[1]], time = 1)
      out_lists <- save_sim_outputs(sim_gizzard, out_lists, 1, timestep)

      reflux_freq_hrs_duo <- as.numeric(Parameters_table_night$reflux_freq_hrs[Parameters_table_night$Compartment == "Duodenum"])

      # Reflux happens on every Nth hour of the night mode (reverse peristalsis duodenum->gizzard)
      # NO REFLUX: remove fraction from gizzard for transfer to duodenum
      if (!(timestep %% reflux_freq_hrs_duo == 1)) {
        perc_to_remove_gizzard_bac  <- Parameters_table_night$remove_bac_perc[Parameters_table_night$Compartment == "Gizzard"]
        perc_to_remove_gizzard_mets <- Parameters_table_night$remove_sublb_perc[Parameters_table_night$Compartment == "Gizzard"]
        removefunction_result_gizzard <- RemoveBacMetabolites(sim_gizzard,
                                                              perc_to_remove_bac = perc_to_remove_gizzard_bac,
                                                              perc_to_remove_mets = perc_to_remove_gizzard_mets,
                                                              type_specie_table = type_specie_table_gizzard)
        arena2_gizzard <- removefunction_result_gizzard[[1]]
      } else {
        # IF REFLUX: no removal from gizzard (will receive from duodenum instead)
        arena2_gizzard <- BacArena::getArena(sim_gizzard, 1)
      }

      compartment_arenas[[1]] <- arena2_gizzard


      ####################################
      ############# DUODENUM #############
      ####################################
      # Imitate absorption of sugars and conversion of starches in duodenum
      compartment_arenas[[2]] <- AbsorptionConversion_duodenum(arena2 = compartment_arenas[[2]])

      sim_duodenum <- BacArena::simEnv(compartment_arenas[[2]], time = 1)
      out_lists <- save_sim_outputs(sim_duodenum, out_lists, 2, timestep)

      # NO REFLUX: stuff transferred from duodenum downstream
      if (!(timestep %% as.numeric(Parameters_table_night$reflux_freq_hrs[Parameters_table_night$Compartment == "Duodenum"]) == 1)) {
        perc_to_remove_duodenum_bac  <- Parameters_table_night$remove_bac_perc[Parameters_table_night$Compartment == "Duodenum"]
        perc_to_remove_duodenum_mets <- Parameters_table_night$remove_sublb_perc[Parameters_table_night$Compartment == "Duodenum"]

        arena2_duodenum <- getArena(sim_duodenum, 1)

        removefunction_result_duodenum <- RemoveBacMetabolites(sim_duodenum, arena2_object = arena2_duodenum,
                                                               perc_to_remove_bac = perc_to_remove_duodenum_bac,
                                                               perc_to_remove_mets = perc_to_remove_duodenum_mets,
                                                               type_specie_table = type_specie_table_duodenum)
        arena2_duodenum <- removefunction_result_duodenum[[1]]

        # Add mets and bacs that were removed from gizzard
        arena2_duodenum <- AddBacMets_fromUpstreamCompartment(removefunction_result_n0 = removefunction_result_gizzard,
                                                              sim_object_n1 = sim_duodenum, arena2_n1 = arena2_duodenum,
                                                              model_list_n1 = model_lists$Duodenum,
                                                              model_list_n0 = model_lists$Gizzard)

      # REFLUX: reverse peristalsis duodenum -> gizzard
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

        # Update gizzard type-specie table with any new species from reflux
        type_specie_table_gizzard <- update_type_specie_table(compartment_arenas[[1]], type_specie_table_gizzard)

        # Remove strict anaerobes from gizzard (if transferred during reflux)
        compartment_arenas[[1]] <- remove_organisms_from_arena(compartment_arenas[[1]], strict_anaerobes,
                                                               type_specie_table_gizzard, "gizzard")
      }

      compartment_arenas[[2]] <- arena2_duodenum


      ####################################
      ############# JEJUNUM ##############
      ####################################
      # Imitate absorption of sugars, amino acids and fats and conversion of starches in jejunum
      compartment_arenas[[3]] <- AbsorptionConversion_jejunum(arena2 = compartment_arenas[[3]], diet)

      sim_jejunum <- BacArena::simEnv(compartment_arenas[[3]], time = 1)
      out_lists <- save_sim_outputs(sim_jejunum, out_lists, 3, timestep)

      if (timestep == 1) {
        type_specie_table_jejunum <- CreateTypeSpecieAssoc_table(sim_jejunum)
      }

      arena2_jejunum <- getArena(sim_jejunum, 1)

      # NO REFLUX: transfer from jejunum downstream
      if (!(timestep %% as.numeric(Parameters_table_night$reflux_freq_hrs[Parameters_table_night$Compartment == "Jejunum"]) == 1)) {
        perc_to_remove_jejunum_bac  <- Parameters_table_night$remove_bac_perc[Parameters_table_night$Compartment == "Jejunum"]
        perc_to_remove_jejunum_mets <- Parameters_table_night$remove_sublb_perc[Parameters_table_night$Compartment == "Jejunum"]
        removefunction_result_jejunum <- RemoveBacMetabolites(sim_jejunum, arena2_object = arena2_jejunum,
                                                              perc_to_remove_bac = perc_to_remove_jejunum_bac,
                                                              perc_to_remove_mets = perc_to_remove_jejunum_mets,
                                                              type_specie_table = type_specie_table_jejunum)
        arena2_jejunum <- removefunction_result_jejunum[[1]]

        # Add bacs and mets from duodenum
        arena2_jejunum <- AddBacMets_fromUpstreamCompartment(removefunction_result_n0 = removefunction_result_duodenum,
                                                             sim_object_n1 = sim_jejunum, arena2_n1 = arena2_jejunum,
                                                             model_list_n1 = model_lists$Jejunum,
                                                             model_list_n0 = model_lists$Duodenum)

      # REFLUX: reverse peristalsis jejunum -> duodenum
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

      compartment_arenas[[3]] <- arena2_jejunum


      ###################################
      ############## ILEUM ##############
      ###################################
      # Imitate absorption of sugars, amino acids and fats; and conversion of starches in ileum
      compartment_arenas[[4]] <- AbsorptionConversion_ileum(arena2 = compartment_arenas[[4]], diet)

      sim_ileum <- BacArena::simEnv(compartment_arenas[[4]], time = 1)
      out_lists <- save_sim_outputs(sim_ileum, out_lists, 4, timestep)

      if (timestep == 1) {
        type_specie_table_ileum <- CreateTypeSpecieAssoc_table(sim_ileum)
      }

      arena2_ileum <- getArena(sim_ileum, 1)

      # NO REFLUX: transfer from ileum to cecum/colon
      if (!(timestep %% as.numeric(Parameters_table_night$reflux_freq_hrs[Parameters_table_night$Compartment == "Ileum"]) == 1)) {
        perc_to_remove_ileum_bac  <- Parameters_table_night$remove_bac_perc[Parameters_table_night$Compartment == "Ileum"]
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

        # Split removed contents between cecum and colon
        split_result <- split_removal_to_cecum_colon(removefunction_result_ileum, fraction_to_cecum, fraction_to_colon)
        removefunction_result_ileum_to_cecum <- split_result$to_cecum
        removefunction_result_ileum_to_colon <- split_result$to_colon

        # Add bacs and mets from jejunum
        arena2_ileum <- AddBacMets_fromUpstreamCompartment(removefunction_result_n0 = removefunction_result_jejunum,
                                                           sim_object_n1 = sim_ileum, arena2_n1 = arena2_ileum,
                                                           model_list_n1 = model_lists$Ileum,
                                                           model_list_n0 = model_lists$Jejunum)

      # REFLUX: transfer part of ileal arena up to jejunum
      # (in the 'colon' section, part of colonic arena will be transferred up to ileum)
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

      compartment_arenas[[4]] <- arena2_ileum


      ####################################
      ############  CECUM  ###############
      ####################################
      sim_cecum <- BacArena::simEnv(compartment_arenas[[5]], time = 1)
      out_lists <- save_sim_outputs(sim_cecum, out_lists, 5, timestep)

      if (timestep == 1) {
        type_specie_table_cecum <- CreateTypeSpecieAssoc_table(sim_cecum)
      }

      # Every 4 hours, release fraction of cecal contents to colon
      if (timestep %% 4 == 0) {
        perc_to_remove_cecum_bac  <- Parameters_table_night$remove_bac_perc[Parameters_table_night$Compartment == "Cecum"]
        perc_to_remove_cecum_mets <- Parameters_table_night$remove_sublb_perc[Parameters_table_night$Compartment == "Cecum"]
        removefunction_result_cecum <- RemoveBacMetabolites(sim_cecum,
                                                            perc_to_remove_bac = perc_to_remove_cecum_bac,
                                                            perc_to_remove_mets = perc_to_remove_cecum_mets,
                                                            type_specie_table = type_specie_table_cecum)
        arena2_cecum <- removefunction_result_cecum[[1]]
      } else {
        arena2_cecum <- BacArena::getArena(sim_cecum, 1)
      }

      # NO REFLUX: add metabolites and bacteria removed from ileum
      if (!(timestep %% as.numeric(Parameters_table_night$reflux_freq_hrs[Parameters_table_night$Compartment == "Ileum"]) == 1)) {
        arena2_cecum <- AddBacMets_fromUpstreamCompartment(removefunction_result_n0 = removefunction_result_ileum_to_cecum,
                                                           sim_object_n1 = sim_cecum, arena2_n1 = arena2_cecum,
                                                           model_list_n1 = model_lists$Ceca,
                                                           model_list_n0 = model_lists$Ileum)
      }
      # IF REFLUX: nothing is transferred from ileum (fraction of colon will be transferred upstream in colon section)

      # Update type-specie table for cecum
      type_specie_table_cecum <- update_type_specie_table(arena2_cecum, type_specie_table_cecum)

      # Remove obligate aerobes from cecum
      arena2_cecum <- remove_organisms_from_arena(arena2_cecum, obligate_aerobes,
                                                  type_specie_table_cecum, "cecum")

      compartment_arenas[[5]] <- arena2_cecum


      ####################################
      ############  COLON  ###############
      ####################################
      sim_colon <- BacArena::simEnv(compartment_arenas[[6]], time = 1)
      out_lists <- save_sim_outputs(sim_colon, out_lists, 6, timestep)

      if (timestep == 1) {
        type_specie_table_colon <- CreateTypeSpecieAssoc_table(sim_colon)
      }

      ## During the night nothing gets removed from colon (no excretion)

      # Once every 4 hrs cecal contents are released to the colon
      if (timestep %% 4 == 0) {
        arena2_colon <- AddBacMets_fromUpstreamCompartment(removefunction_result_n0 = removefunction_result_cecum,
                                                           sim_object_n1 = sim_colon, arena2_n1 = arena2_colon,
                                                           model_list_n1 = model_lists$Colon,
                                                           model_list_n0 = model_lists$Ceca)
      }

      # NO REFLUX: add metabolites and bacteria from ileum
      if (!(timestep %% as.numeric(Parameters_table_night$reflux_freq_hrs[Parameters_table_night$Compartment == "Ileum"]) == 1)) {
        arena2_colon <- AddBacMets_fromUpstreamCompartment(removefunction_result_n0 = removefunction_result_ileum_to_colon,
                                                           sim_object_n1 = sim_colon, arena2_n1 = arena2_colon,
                                                           model_list_n1 = model_lists$Colon,
                                                           model_list_n0 = model_lists$Ileum)
      } else {
        ## IF REFLUX: nothing transferred from ileum downstream;
        # instead, fraction of colon is transferred upstream to cecum and ileum

        perc_to_remove_colon_bac  <- as.numeric(Parameters_table_night$reflux_perc[Parameters_table_night$Compartment == "Colon"])
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

        # Split removed contents between ileum and cecum
        removefunction_result_colon_to_ileum <- removefunction_result_colon
        removefunction_result_colon_to_cecum <- removefunction_result_colon
        # Sample fraction_to_ileum of orgdat removed from colon
        removefunction_result_colon_to_ileum[[2]] <- removefunction_result_colon_to_ileum[[2]][
          sample(nrow(removefunction_result_colon_to_ileum[[2]]),
                 ceiling(fraction_to_ileum * nrow(removefunction_result_colon_to_ileum[[2]]))), ]
        # Remaining goes to cecum
        removefunction_result_colon_to_cecum[[2]] <- dplyr::setdiff(removefunction_result_colon[[2]],
                                                                     removefunction_result_colon_to_ileum[[2]])

        # Split sublb (metabolites)
        removefunction_result_colon_to_ileum[[3]] <- removefunction_result_colon_to_ileum[[3]][
          sample(nrow(removefunction_result_colon_to_ileum[[3]]),
                 round(fraction_to_ileum * nrow(removefunction_result_colon_to_ileum[[3]]))), ]
        removefunction_result_colon_to_cecum[[3]] <- dplyr::setdiff(removefunction_result_colon[[3]],
                                                                     removefunction_result_colon_to_ileum[[3]])

        # Update species_rm slot [[4]]
        type_number_colon_to_ileum <- table(removefunction_result_colon_to_ileum[[2]]$type)
        type_number_colon_to_cecum <- table(removefunction_result_colon_to_cecum[[2]]$type)
        removefunction_result_colon_to_ileum[[4]] <- subset(removefunction_result_colon_to_ileum[[4]],
                                                            removefunction_result_colon_to_ileum[[4]]$Type %in% names(type_number_colon_to_ileum))
        removefunction_result_colon_to_cecum[[4]] <- subset(removefunction_result_colon_to_cecum[[4]],
                                                            removefunction_result_colon_to_cecum[[4]]$Type %in% names(type_number_colon_to_cecum))
        removefunction_result_colon_to_ileum[[4]]$Number <- type_number_colon_to_ileum
        removefunction_result_colon_to_cecum[[4]]$Number <- type_number_colon_to_cecum
      }

      # Update type-specie table for colon
      type_specie_table_colon <- update_type_specie_table(arena2_colon, type_specie_table_colon)

      # Remove obligate aerobes from colon
      arena2_colon <- remove_organisms_from_arena(arena2_colon, obligate_aerobes,
                                                  type_specie_table_colon, "colon")

      compartment_arenas[[6]] <- arena2_colon


      #####################################################
      ############   ILEUM AND CECUM AGAIN  ###############
      ##########  TO ADD REFLUXES FROM COLON  #############
      #####################################################
      # REFLUX: transfer part of colonic arena up to ileum and cecum
      if (timestep %% as.numeric(Parameters_table_night$reflux_freq_hrs[Parameters_table_night$Compartment == "Colon"]) == 1) {
        # Add metabolites and bacteria removed from colon to ileum
        compartment_arenas[[4]] <- AddBacMets_fromDownstreamCompartment(removefunction_result_colon_to_ileum,
                                                                        sim_object_n0 = sim_ileum,
                                                                        arena2_n0 = compartment_arenas[[4]],
                                                                        model_list_n0 = model_lists$Ileum,
                                                                        model_list_n1 = model_lists$Colon)

        # Imitate absorption of sugars, amino acids and fats in ileum after reflux
        compartment_arenas[[4]] <- AbsorptionConversion_ileum(arena2 = compartment_arenas[[4]], diet)

        # Add metabolites and bacteria removed from colon to cecum
        compartment_arenas[[5]] <- AddBacMets_fromDownstreamCompartment(removefunction_result_colon_to_cecum,
                                                                        sim_object_n0 = sim_cecum,
                                                                        arena2_n0 = compartment_arenas[[5]],
                                                                        model_list_n0 = model_lists$Ceca,
                                                                        model_list_n1 = model_lists$Colon)
      }

    } # end night time block

  } # end timestep loop


  # ---- Assemble output list ----
  if (is.null(AminoAcid_X)) { AminoAcid_X <- "None" }
  if (is.null(Sugars_X))    { Sugars_X <- "None" }

  simlist <- c(
    out_lists,                      # slots 1-6: compartment outputs
    list(seedN),                    # slot 7: seed
    list(Parameters_table),         # slot 8: day parameters
    list(Parameters_table_night),   # slot 9: night parameters
    list(ExchangeConstraints),      # slot 10: exchange constraints
    list(list(                      # slot 11: all simulation arguments
      econc_proportion = econc_proportion,
      prop_gizzard     = prop_gizzard,
      Sugars_X         = Sugars_X,
      AminoAcids_X     = AminoAcid_X,
      timesteps        = timesteps,
      sample           = sample,
      diet_name        = diet_name
    ))
  )
  return(simlist)
})))

stopCluster(cl)

# ---- Save output ----
now   <- Sys.time()
today <- strsplit(as.character(as.POSIXlt(now)[1]), "\\ ")[[1]][1]
stir_TF <- Parameters_table$stir[1]

output_dir <- file.path(base_dir, "arena_sim_objects/CoupledSixComps_D10/NoPathogens/2024/gapseq_models")
saveRDS(simlist, file.path(output_dir,
                           paste0(sample, "_", as.character(diet_name),
                                  "_diet_grower_NEW_AFTER48h_", as.character(timesteps + 48),
                                  "h_econc", econc_proportion,
                                  "_propGizzard", as.character(prop_gizzard),
                                  "_CobaltX12_rmAnaero_less02_mucUreaCecaColon_wCoA_moreNightFlow_",
                                  today,
                                  "_gapseq.rds")))
