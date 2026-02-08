
myargs <- commandArgs(trailingOnly = TRUE)
print(myargs)

library(plyr)
library(dplyr)
library(magrittr)
library(BacArena)
library(glpkAPI)
library(Rcpp)
library(igraph)
library(reshape2)
library(libSBML)
library(sybilSBML)
library(stringr)
library(parallel)

# ===== Loading input parameters =============================

# proportion of concentration of essential mets 
econc_proportion = as.numeric(myargs[1])

# proportion of diet being added to the
diet_proportion = as.numeric(myargs[2])

# bird identifier (e.g., G1_1 for Group G1, Bird 1)
bird_id = as.character(myargs[3])

# number of timesteps (hours)
timesteps = as.numeric(myargs[4])

# simulation date in YYYY-MM-DD format
sim_date = as.character(myargs[5])

# optional: supplement name (if testing supplements)
compound = ifelse(length(myargs) >= 6, as.character(myargs[6]), "none")

# optional: supplement concentration (g per 100g)
compound_conc = ifelse(length(myargs) >= 7, as.numeric(myargs[7]), 0)

if (is.na(compound_conc)) { compound_conc <- 0}

#### ===== Load input files ================================
file_path <- "Chicken_IleumCecum_2024Trial_simulation/"
input_files_path = paste0(file_path, "ModelGeneration_Files/")
functions_scripts_path = paste0(file_path, "ModelGeneration_scripts/")

diet_corn_ileum <- read.table(paste0(input_files_path, "corn_diet_grower_ileum.txt"),
                        sep = '\t', header = T)

# Load parameters for 2-compartment model
Parameters_table <- read.table(paste0(input_files_path,"Parameters_IleumCecum_Model_2025.txt"),sep='\t',header=T)

otu_table_metadata <- read.table(paste0(input_files_path, "RepresentativeOTUs_matchedMAGs_with_metadata_paired_birds_corrected.tsv"),sep ='\t', header=T)

CompoundsToTest_table <- read.table(paste0(input_files_path, "CompoundsTable_to_test.txt"),sep ='\t', header=T)

# load all the functions 
functions_scripts <- list.files(path = functions_scripts_path, pattern = "Function.*.R", full.names = TRUE)

# Load bird-specific models
model_list_ileum <- readRDS(file.path(input_files_path, paste0("model_lists/",bird_id, "_Ileum_ModelList_July2025_v2.rds")))
model_list_cecum <- readRDS(file.path(input_files_path, paste0("model_lists/",bird_id, "_Cecum_ModelList_July2025_v2.rds")))

otu_table_metadata$species_ID <- as.character(otu_table_metadata$species_ID)
otu_table_metadata$bird_id <- as.character(otu_table_metadata$bird_id)
otu_table_metadata$site <- as.character(otu_table_metadata$site)

# Clean collision detection (when the same specie_id in ileum and cecum were reconstructed from different MAGs)
fix_species_collisions <- function(model_list_ileum, model_list_cecum, metadata, bird_id) {
  
  # Filter metadata for current bird
  ileum_meta <- metadata[metadata$bird_id == bird_id & metadata$site == "Ileum", ]
  cecum_meta <- metadata[metadata$bird_id == bird_id & metadata$site == "Cecum", ]
  
  # Get species names
  ileum_species <- sapply(model_list_ileum, function(m) m@mod_desc)
  cecum_species <- sapply(model_list_cecum, function(m) m@mod_desc)
  common_species <- intersect(ileum_species, cecum_species)
  cat("Found", length(common_species)," common species in both compartments\n")

  # Track renames for metadata update
  renames <- list()
  
  for (species_name in common_species) {
    # Get MAG IDs from metadata
    ileum_mag <- as.character(ileum_meta$mag_id_clean[ileum_meta$species_ID == species_name][1])
    cecum_mag <-as.character(cecum_meta$mag_id_clean[cecum_meta$species_ID == species_name][1])
    
    if (!(ileum_mag == cecum_mag)) {
      cat("Collision:", species_name, "- different MAGs (", ileum_mag, "vs", cecum_mag, ")\n")
      
      # Find and rename models
      ileum_idx <- which(ileum_species == species_name)[1]
      cecum_idx <- which(cecum_species == species_name)[1]
      
      # Rename with MAG suffix
      new_ileum_name <- paste0(species_name, "_ileum_str")
      new_cecum_name <- paste0(species_name, "_cecum_str")
      
      model_list_ileum[[ileum_idx]]@mod_desc <- new_ileum_name
      model_list_ileum[[ileum_idx]]@mod_name <- new_ileum_name
      model_list_cecum[[cecum_idx]]@mod_desc <- new_cecum_name  
      model_list_cecum[[cecum_idx]]@mod_name <- new_cecum_name
      
      # Track renames for metadata update
      renames[[species_name]] <- list(
        ileum = new_ileum_name,
        cecum = new_cecum_name
      )
      
      cat("  Renamed to:", new_ileum_name, "and", new_cecum_name, "\n")
    }
  }
  
  # UPDATE METADATA with new species names
  updated_metadata <- metadata
  for (old_name in names(renames)) {
    # Update ileum entries
    ileum_mask <- (updated_metadata$species_ID == old_name) & (updated_metadata$site == "Ileum") & (updated_metadata$bird_id == bird_id)
    updated_metadata$species_ID[ileum_mask] <- renames[[old_name]]$ileum
    
    # Update cecum entries  
    cecum_mask <- (updated_metadata$species_ID == old_name) & (updated_metadata$site == "Cecum") & (updated_metadata$bird_id == bird_id)
    updated_metadata$species_ID[cecum_mask] <- renames[[old_name]]$cecum
  }
  
  return(list(
    ileum = model_list_ileum, 
    cecum = model_list_cecum,
    metadata = updated_metadata
  ))
}

fixed_models <- fix_species_collisions(model_list_ileum, model_list_cecum, 
                                       otu_table_metadata, bird_id)
model_list_ileum <- fixed_models$ileum
model_list_cecum <- fixed_models$cecum
otu_table_metadata <- fixed_models$metadata
model_lists <- list(model_list_ileum, model_list_cecum)



##########################################################
############   UPLOAD DATA ON CLUSTER  ###################
##########   RUN MODELS WITH REPLICATES  #################
##########################################################

if (!(compound == "none")) {
  # get only the row of the table corresponding to the compound chosen
  CompoundsToTest_table_sub <- CompoundsToTest_table[which(
    CompoundsToTest_table$variable == compound & 
      CompoundsToTest_table$value == compound_conc),]
  
  # now add this row to a diet table,  so that compound is introduced into feeding:
  diet_corn_ileum <- rbind(diet_corn_ileum, CompoundsToTest_table_sub[,colnames(diet_corn_ileum)])
  # and sum the values for the compound if it was present in the diet:
  diet_corn_ileum <- aggregate(met_mM ~  unit + met_name + met_id + category, diet_corn_ileum, sum)
  
}


# Initialize cluster
cores <- 5
replicates <- 5

now <- Sys.time()
today <- strsplit(as.character(as.POSIXlt(now)[1]), "\\ ")[[1]][1]

cl <- makeCluster(cores, type="PSOCK",
                  outfile=paste0(file_path, "cluster_logs/cluster_coupling_ileum_cecum_",
                                 timesteps,"h_",
                                 bird_id, "_propEconc", econc_proportion,"_",
                                 compound, "_", as.character(compound_conc), "_pant20x_dca_rm0.8met_updmods_lessIlealBac",
                                 "_diet_proportion", diet_proportion,"_", today,".log"))

clusterExport(cl, c("functions_scripts"))
# load all the scripts with functions
clusterCall(cl, function() { lapply(functions_scripts, source) })

# export all the parameters and input files
clusterExport(cl, c("bird_id", "diet_proportion", "Parameters_table", "econc_proportion"))
clusterExport(cl, c("model_lists", "diet_corn_ileum",  "otu_table_metadata")) # export model lists, diet, abundances and metadata
clusterExport(cl, c("timesteps", "compound", "compound_conc"))

print(system.time(simlist <- parLapply(cl, 1:replicates, function(i){ #do 5 replicate simulations
  
  # set seed for this replicate  (it will be used in the creation of initial arena (CreateModelArena function), 239 + seedN)
  N <- sample(c(1:100), 1)
  seedN <- as.integer(239 + N)
  print(paste0("SEED FOR THIS REPLICATE IS ", as.character(seedN)))
  
  
  ## ------------------------
  ## HELPER FUNCTIONS
  ## ------------------------
  set_pH <- function(arena, pH) {
    arena <- BacArena::addSubs(arena, smax=1000/(10^pH), mediac="EX_cpd00067_e0", unit="mM", add=F) 
    return(arena)
  }
  
  set_o2 <- function(arena, max_conc) {
    arena <- BacArena::addSubs(arena, smax = max_conc, mediac="EX_cpd00007_e0", unit="mM", add=F) 
    return(arena)
  }
  
  ### Set pH and oxygen for both compartments
  set_compartment_conditions <- function(arenas, compartment_names, params_table) {
    for (j in 1:2) {
      arenas[[j]] <- set_pH(arenas[[j]], params_table$pH[params_table$Compartment == compartment_names[j]])
      arenas[[j]] <- set_o2(arenas[[j]], as.numeric(params_table$oxygen_mM[params_table$Compartment == compartment_names[j]]))
    }
    return(arenas)
  }
  
  ### Save simulation outputs (identical for both compartments)
  save_simulation_outputs <- function(sim_obj, timestep) {
    outputs <- list()
    outputs$Population <- BacArena::plotCurves(sim_obj, retdata = T, graph = F)$Population
    outputs$Substances <- BacArena::plotCurves(sim_obj, retdata = T, graph = F)$Substances
    outputs$Orgdat <- sim_obj@simlist
    outputs$Fluxlist <- sim_obj@mfluxlist
    
    # save only negative shadow prices
    shadowlist <- reshape2::melt(unlist(lapply(sim_obj@shadowlist[[2]], function(x) x[x < -0.1])))
    shadowlist$Metabolite <- unlist(lapply(rownames(shadowlist), function(x) {strsplit(x, "EX_")[[1]][2]}))
    shadowlist$Specie <- unlist(lapply(rownames(shadowlist), function(x) { gsub("\\.$", "",  strsplit(x, "EX_")[[1]][1]) }))
    outputs$Shadowlist <- shadowlist
    
    return(outputs)
  }
  
  ### Add diet and essential metabolites to ileum
  add_diet_to_ileum <- function(arena, diet, diet_proportion, econc_proportion, feeding_rate = 1, 
                                essential_metabolites) {
    # add metabolites (diet without metabolites that are absorbed upstream)
    arena <- BacArena::addSubs(arena,
                               smax = feeding_rate * diet_proportion * diet$met_mM,
                               unit="mM",
                               mediac = as.character(diet$met_id), add = T, addAnyway = T)
    
    # add essential metabolites that aren't in the diet:
    essential_metabolites_not_in_diet <- setdiff(essential_metabolites, as.character(diet$met_id))
    arena <- BacArena::addSubs(arena, smax = feeding_rate * econc_proportion * Parameters_table$essential_conc_mM[1],
                               unit="mM",
                               mediac = as.character(essential_metabolites_not_in_diet), add = T, addAnyway = T)
    return(arena)
  }
  
  ### Process ileum removal and split transfer to cecum/colon
  process_ileum_removal <- function(sim_ileum, arena2_ileum, type_specie_table_ileum, params_table, is_night = FALSE) {
    # select % of the grid randomly that will be removed in the next iteration
    # one part will be transferred to colon, and remaining part to cecum
    # assuming that majority of bacteria from ileum transit directly to colon / do not survive in anaerobic cecum
    # while half of the contents leaving ileum (metabolites) ends up in cecum  (half in colon)
    if (is_night) {
      perc_to_remove_ileum_bac <- params_table$remove_bac_perc_night[params_table$Compartment == "Ileum"]
      perc_to_remove_ileum_mets <- params_table$remove_sublb_perc_night[params_table$Compartment == "Ileum"]
      perc_to_remove_ileum_to_cecum <- as.numeric(params_table$remove_from_ileum_night[params_table$Compartment == "Cecum"])
    } else {
      perc_to_remove_ileum_bac <- params_table$remove_bac_perc[params_table$Compartment == "Ileum"]
      perc_to_remove_ileum_mets <- params_table$remove_sublb_perc[params_table$Compartment == "Ileum"]
      perc_to_remove_ileum_to_cecum <- as.numeric(params_table$remove_from_ileum[params_table$Compartment == "Cecum"])
    }
    
    fraction_to_cecum <- perc_to_remove_ileum_to_cecum / perc_to_remove_ileum_bac
    fraction_to_cecum_mets <- 0.5
    # get specs, orgdat and metabolites removed from ileum in removefunction_result_ileum list
    removefunction_result_ileum <- RemoveBacMetabolites(sim_ileum, arena2_object=arena2_ileum,
                                                        perc_to_remove_bac=perc_to_remove_ileum_bac,
                                                        perc_to_remove_mets=perc_to_remove_ileum_mets,
                                                        type_specie_table=type_specie_table_ileum)
    
    # get the part that will be transferred to cecum (note: slot [[1]], i.e. arena_ileum doesn't matter here):
    removefunction_result_ileum_to_cecum <- removefunction_result_ileum
    
    # first, sample fraction_to_cecum of orgdat removed from ileum:
    removefunction_result_ileum_to_cecum[[2]] <- removefunction_result_ileum_to_cecum[[2]][
      sample(nrow(removefunction_result_ileum_to_cecum[[2]]),
             ceiling(fraction_to_cecum * nrow(removefunction_result_ileum_to_cecum[[2]]))), ]
    
    # secondly, do same for sublb (metabolites) table:
    removefunction_result_ileum_to_cecum[[3]] <- removefunction_result_ileum_to_cecum[[3]][
      sample(nrow(removefunction_result_ileum_to_cecum[[3]]),
             round(fraction_to_cecum_mets * nrow(removefunction_result_ileum_to_cecum[[3]]))), ]
    
    # last, change the slot [[4]] - species_rm (table with names of specie models and numbers of them that should be transferred)
    type_number_ileum_to_cecum <- table(removefunction_result_ileum_to_cecum[[2]]$type) # get type and numbers of each type from orgdat
    
    # leave in [[4]] slot only types of bacteria being transferred to cecum/colon
    # make sure you don't try to transfer types of bacteria not present
    removefunction_result_ileum_to_cecum[[4]] <- subset(removefunction_result_ileum_to_cecum[[4]], 
                                                        removefunction_result_ileum_to_cecum[[4]]$Type %in% names(type_number_ileum_to_cecum))
    # update numbers of each type to be transferred to cecum/colon:
    removefunction_result_ileum_to_cecum[[4]]$Number <- type_number_ileum_to_cecum
    
    return(list(
      updated_arena = removefunction_result_ileum[[1]],
      transfer_to_cecum = removefunction_result_ileum_to_cecum
    ))
  }

  
  ### Process cecum simulation and add transferred material from ileum
  process_cecum <- function(sim_cecum, timestep, type_specie_table_cecum, removefunction_result_ileum_to_cecum, 
                            model_lists, params_table, is_night = FALSE) {
    
    # select pre-set (in Parameters_table) % of the grid randomly that will be removed once in 4 hours (released to colon)
    if (timestep %% 4 == 0) {
      if (is_night) {
        perc_to_remove_cecum_bac <- params_table$remove_bac_perc_night[params_table$Compartment == "Cecum"]
        perc_to_remove_cecum_mets <- params_table$remove_sublb_perc_night[params_table$Compartment == "Cecum"]
      } else {
        perc_to_remove_cecum_bac <- params_table$remove_bac_perc[params_table$Compartment == "Cecum"]
        perc_to_remove_cecum_mets <- params_table$remove_sublb_perc[params_table$Compartment == "Cecum"]
      }
      
      # get specs, orgdat and metabolites removed from cecum in removefunction_result_cecum list
      removefunction_result_cecum <- RemoveBacMetabolites(sim_cecum, 
                                                          perc_to_remove_bac=perc_to_remove_cecum_bac,
                                                          perc_to_remove_mets=perc_to_remove_cecum_mets,
                                                          type_specie_table=type_specie_table_cecum)
      arena2_cecum <- removefunction_result_cecum[[1]]
    } else {
      arena2_cecum <- BacArena::getArena(sim_cecum, 1)
    }
    
    # add metabolites and bacteria removed from the ileum (i.e. X% of ileal contents)
    arena2_cecum <- AddBacMets_fromUpstreamCompartment(removefunction_result_n0 = removefunction_result_ileum_to_cecum,
                                                       sim_object_n1 = sim_cecum, arena2_n1 = arena2_cecum,
                                                       model_list_n0 = model_lists[[1]])
    
    # if any new species happen to be transferred in cecum - add them to the type_specie_table:
    new_bactype_arena2_cecum <- setdiff(unique(arena2_cecum@orgdat$type), type_specie_table_cecum$Type)
    if (length(new_bactype_arena2_cecum) > 0) {
      for (bactype in new_bactype_arena2_cecum) {
        type_specie_new <- data.frame(cbind(arena2_cecum@specs[[bactype]]@model@mod_desc, bactype))
        colnames(type_specie_new) <- c("Name", "Type")
        type_specie_table_cecum <- rbind(type_specie_table_cecum, type_specie_new)
      }
    }
    
    
    return(list(
      updated_arena = arena2_cecum,
      updated_type_table = type_specie_table_cecum
    ))
  }

  

  ## ---------------------------------------------------------------
  ## Pick the highest trade-off that satisfies the fraction-growing criterion and
  ## keeps mu_star in a plausible range: (1) growers ≥ min_frac  and  (2) max mu_star ≤ max_mu
  ## ---------------------------------------------------------------
  choose_tau <- function(comm_mod,
                         taus       = seq(0.2, 0.9, 0.05),
                         abund_vec,
                         min_frac   = 0.90,
                         max_mu     = 3,         # ~ realistic upper limit (h-1)
                         thresh     = 1e-6,
                         verbose    = TRUE) {

    # Store both criteria for every candidate τ
    grow_frac <- numeric(length(taus))
    mu_max    <- numeric(length(taus))

    for (j in seq_along(taus)) {
      mu <- ctFBA(comm_mod, tau = taus[j], abund_vec = abund_vec, verbose = TRUE)
      grow_frac[j] <- mean(mu > thresh)
      mu_max[j]    <- max(mu)
      if (verbose)
        message(sprintf("τ = %.2f  -> growers %.2f  max μ %.3f",
                        taus[j], grow_frac[j], mu_max[j]))
    }

    ok <- which(grow_frac >= min_frac & mu_max <= max_mu)

    if (length(ok)) {
      tau_star <- taus[max(ok)]            # **highest** τ that passes both
    } else {
      # fallback: use default tau=0.5
      tau_star = 0.5
      # ok2 <- which(grow_frac >= min_frac)
      # tau_star <- if (length(ok2)) taus[max(ok2)] else taus[which.max(grow_frac)]
      warning("No tau satisfied the mu_max criterion – using ", tau_star)
    }

    list(tau   = tau_star,
         table = data.frame(tau = taus,
                            growing = grow_frac,
                            mu_max  = mu_max))
  }
  
  ## ---------------------------------------------------
  ## --------------- MAIN SIMULATION SETUP -------------
  ## ---------------------------------------------------
  compartments = c("Ileum","Cecum")
  # Initialize arenas for compartments
  setup_lists <- list(
    CreateModelArena(model_list  = model_lists[[1]],  # Ileum
                     bird_id     = bird_id,
                     compartment = "Ileum",
                     diet        = diet_corn_ileum,
                     diet_proportion = diet_proportion,          # ileum gets feed
                     initial_commsize = Parameters_table$initial_community_size[1],
                     econc       = econc_proportion * Parameters_table$essential_conc_mM[1],
                     stir_var    = Parameters_table$stir[1],
                     arena_size  = 100, speed = 5, tstep = 1,
                     seed        = seedN,
                     metadata_file = otu_table_metadata),
    
    CreateModelArena(model_list  = model_lists[[2]],  # Cecum
                     bird_id     = bird_id,
                     compartment = "Cecum",
                     diet        = diet_corn_ileum,
                     diet_proportion = 0,                       # cecum not fed directly
                     initial_commsize = Parameters_table$initial_community_size[2],
                     econc       = econc_proportion * Parameters_table$essential_conc_mM[2],
                     stir_var    = Parameters_table$stir[2],
                     arena_size  = 100, speed = 5, tstep = 1,
                     seed        = seedN,
                     metadata_file = otu_table_metadata)
  )
  
  # get arenas
  compartment_arenas          <- lapply(setup_lists, `[[`, 1)  # keep only arenas
  names(compartment_arenas)   <- compartments
  # get essential mets for ileum
  essential_metabolites_ileum <- setup_lists[[1]][[2]] 
  
  # remove oxygen from essential_metabolites: 
  essential_metabolites_ileum <- essential_metabolites_ileum[which(!(essential_metabolites_ileum == "EX_cpd00007_e0"))]

  
  
  # get type_specie_table - species and corresponding "type" in arena@orgdat
  type_specie_table_ileum <- CreateTypeSpecieAssoc_table(compartment_arenas[[1]])
  type_specie_table_cecum <- CreateTypeSpecieAssoc_table(compartment_arenas[[2]])
  
  
  ################################################################
  ##  Calibrate trade-off tau for each compartment 
  ################################################################
  # build_comm_model() and ctFBA() are defined in Function_ctFBA_BacArena.R
  # At this point each arena already contains (i) the agent counts seeded from OTU table;
  # (ii) the initial medium (diet for ileum; none for cecum).
  # That matches the “initial in vivo snapshot” MICOM uses to calibrate tau
  # Choosing tay once keeps the whole simulation deterministic across the five replicates that differ only by random spatial seeding.
  calibrated_ileum  <- choose_tau(build_comm_model(compartment_arenas[[1]]),
                           taus = seq(0.2, 0.9, 0.05), 
                           abund_vec = as.numeric(table(compartment_arenas[[1]]@orgdat$type)),
                           max_mu = 3.5,
                           min_frac = 0.90)
  tau_star_ileum  <- calibrated_ileum$tau
  
  calibrated_cecum  <- choose_tau(build_comm_model(compartment_arenas[[2]]),
                           taus = seq(0.2, 0.9, 0.05),
                           abund_vec = as.numeric(table(compartment_arenas[[2]]@orgdat$type)),
                           max_mu = 5,
                           min_frac = 0.90)
  tau_star_cecum <- calibrated_cecum$tau
  
  cat("Chosen trade-off tau values  – Ileum:",  tau_star_ileum,
      " | Cecum:", tau_star_cecum, "\n")
  
  ## add metabolites from hour 24 in concentrations from hour 24 instead of 'essential' metabolites 
  # the mets were collected at time point 24 when the proportion of the medium added was X90, that's why conc_multiply = prop_gizzard/90, 
  # i.e. scaling the added concentrations (instead of essential concentrations) to the medium added
  # compartment_arenas <- lapply(1:2, function(j) {
  #   addMets_from24h(compartment_arenas[[j]],sampleName = sample,compartment = names(compartment_arenas)[j],
  #                   conc_multiply = prop_gizzard/unique(Metabolite_conc_h24$proportion_medium))
  # })
  
  # Initialize output lists
  out_lists <- lapply(1:2, function(x) vector("list", length = timesteps + 1))
  
  for (timestep in 1:timesteps) {
    
    names(compartment_arenas) <- c("Ileum","Cecum")
    
    # Determine if day or night - DAY TIME: hours 1-16, 25-40 | NIGHT TIME: hours 17-24, 41-48
    is_day_time <- timestep %in% c(1:16, 25:40)
    is_night_time <- timestep %in% c(17:24, 41:48)
    
    if (is_day_time || is_night_time) {
      print(paste0('timestep ', as.character(timestep)))
      
      # set pH and oxygen
      compartment_arenas <- set_compartment_conditions(compartment_arenas, compartments, Parameters_table)
      
      ########################################################################
      ############                  DAY/NIGHT TIME               #############
      ########################################################################
      
      ###################################
      ############## ILEUM ##############
      ###################################

      # add night rate of metabolites during night time (night_rate = 0.25)
      feeding_rate <- if (is_night_time) 0.25 else 1.0
      compartment_arenas[[1]] <- add_diet_to_ileum(compartment_arenas[[1]], diet_corn_ileum, diet_proportion, 
                                                   econc_proportion, feeding_rate, essential_metabolites_ileum)
      
      
      # Apply cooperative trade-off at each timestep to prevent competitive exclusion
      # Maintain community diversity by ensuring min_fraction of species can grow
      # Use the calibrated tau:
      compartment_arenas[[1]] <- apply_ctFBA_to_arena(
        arena  = compartment_arenas[[1]],
        tau    = tau_star_ileum,
        solver = "glpkAPI", 
        type_specie_table = type_specie_table_ileum)
      
      ### simulate ileum for 1 h
      sim_ileum <- BacArena::simEnv(compartment_arenas[[1]], time = 1) #simulate for 1h
      # save outputs of this iteration
      out_lists[[1]][[timestep]] <- save_simulation_outputs(sim_ileum, timestep)
      
      arena2_ileum <- BacArena::getArena(sim_ileum, 1)
      
      ## imitate absorption of sugars, amino acids and fats in ileum:
      arena2_ileum <- AbsorptionConversion_ileum(arena2 = arena2_ileum, diet_corn_ileum,  diet_proportion)
      
      # NO REFLUX - Process removal and transfer
      ileum_results <- process_ileum_removal(sim_ileum, arena2_ileum, type_specie_table_ileum, 
                                             Parameters_table, is_night = is_night_time)
      
      # update ileal arena
      compartment_arenas[[1]] <- ileum_results$updated_arena
      print(paste0('arena_ileum after ', timestep, ' hours: ', as.character(nrow(compartment_arenas[[1]]@orgdat))))
      
      ####################################
      ############  CECUM  ###############
      ####################################
      # get type_specie_table - species and corresponding "type" in arena@orgdat
      if (timestep == 1) {
        type_specie_table_cecum <- CreateTypeSpecieAssoc_table(compartment_arenas[[2]])
      }
      
      if (timestep > 1) {
        # Newcomers from ileum enter cecum with their caps already set (they inherit them from the model), 
        # so ensuring consistency across compartments
        # Apply cooperative trade-off at each timestep to prevent competitive exclusion
        # Maintain community diversity by ensuring min_fraction of species can grow
        # calibrate tau trade-off (MICOM treats every LP on the whole community each time it changes)
        compartment_arenas[[2]] <- apply_ctFBA_to_arena(
          arena  = compartment_arenas[[2]],
          tau    = tau_star_cecum,
          solver = "glpkAPI",
          type_specie_table = type_specie_table_cecum)
      }
      
     
      
      ### simulate cecum for 1 h
      sim_cecum <- BacArena::simEnv(compartment_arenas[[2]], time = 1) #simulate for 1h
      # save outputs of this step
      out_lists[[2]][[timestep]] <- save_simulation_outputs(sim_cecum, timestep)
    
      
      ## NO REFLUX - Process cecum
      cecum_results <- process_cecum(sim_cecum, timestep, type_specie_table_cecum, 
                                     ileum_results$transfer_to_cecum, model_lists, 
                                     Parameters_table, is_night = is_night_time)
      
      # update cecum arena
      compartment_arenas[[2]] <- cecum_results$updated_arena
      type_specie_table_cecum <- cecum_results$updated_type_table
      print(paste0('arena_cecum after ', timestep, ' hours: ', as.character(nrow(compartment_arenas[[2]]@orgdat))))
    }
  }
  
  
 
  simlist = list()
  simlist[[1]] = out_lists[[1]]
  simlist[[2]] = out_lists[[2]]
  simlist[[3]] = seedN
  simlist[[4]] = Parameters_table
  # save the list with all args:
  simlist[[5]] = list('econc_proportion' = econc_proportion,'diet_proportion' = diet_proportion,
                       'timesteps' = timesteps, 'bird_id' = bird_id, 
                       'compound' = compound, 'compound_conc' = compound_conc)
  return(simlist)
} ))) 

stopCluster(cl)


# Create filename
if (compound == "none") {
  output_filename <- paste0(bird_id, "_2comp_", timesteps, "h_",
                            "econc", econc_proportion, "_diet", as.character(diet_proportion), "_",
                            sim_date, "_L2reg_ileO2_pant20x_dca_rm0.8met_updmods_lessIlealBac.rds")
} else {
  output_filename <- paste0(bird_id, "_2comp_", timesteps, "h_", 
                            "econc", econc_proportion, "_diet", as.character(diet_proportion), "_",
                            compound, "_", compound_conc, "g100g_", sim_date, "_L2reg_ileO2_pant20x_dca_rm0.8met_updmods_lessIlealBac.rds")
}

# create output dir
output_dir <- "Chicken_IleumCecum_2024Trial_simulation/simulation_outputs/gapseq_models/"
# dir.create(output_dir, showWarnings = FALSE)
dir.create(file.path(output_dir, sim_date), showWarnings = FALSE)

saveRDS(simlist, file.path(output_dir, sim_date, output_filename))



