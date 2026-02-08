library(libSBML)
library(sybilSBML)

################################################################################
####################    ABSORPTION RATES DURING THE DAY     ####################
################################################################################

# HELPER FUNCTIONS TO REDUCE REPETITION

# Get current metabolite concentrations from arena
get_metabolite_concentrations <- function(arena2, metabolite_ids) {
  sublb <- BacArena::getSublb(arena2)
  sublb_df <- as.data.frame(sublb)
  
  concentrations <- list()
  for (met_id in metabolite_ids) {
    if (met_id %in% colnames(sublb_df)) {
      concentrations[[met_id]] <- mean(sublb_df[, met_id])
    } else {
      concentrations[[met_id]] <- 0
    }
  }
  return(concentrations)
}

# Apply absorption with the complex diet concentration formula
# diet_concentration * prop_gizzard * 0.5 * 625 * diet$met_mM[diet$met_id == X]/10^12
apply_absorption <- function(arena2, met_id, current_conc, retention_fraction, diet, diet_concentration, prop_gizzard) {
  if (length(diet$met_mM[diet$met_id == met_id]) == 0) {
    # If metabolite not in diet, just return arena unchanged
    return(arena2)
  }
  
  # Calculate target concentration: retention_fraction of diet intake per hour
  # 0.5 - because diet added every 2hrs, 625 - conversion from mM to mmol/cell, 10^12 - from mmol/cell to fmol/cell
  target_conc <- retention_fraction * diet_concentration * prop_gizzard * 0.5 * 625 * diet$met_mM[diet$met_id == met_id] / 10^12
  
  # Set concentration to minimum of current or target (absorption effect)
  arena2 <- BacArena::addSubs(arena2, 
                              smax = min(current_conc, target_conc),
                              mediac = met_id, 
                              add = FALSE, 
                              unit = "mmol/cell")
  return(arena2)
}

# Process amino acid absorption (used in all three compartments)
process_amino_acid_absorption <- function(arena2, diet, retention_fraction, diet_concentration, prop_gizzard, 
                                          category_filter = c("Proteins", "Proteins_shadow")) {
  sublb <- BacArena::getSublb(arena2)
  sublb_df <- as.data.frame(sublb)
  
  # Get amino acids from diet
  AA_diet_concentration <- subset(diet, diet$category %in% category_filter)
  
  # 0.5*diet_concentration - approximately how much diet-metabolites is ingested per hour given that
  # total amount of metabolites ingested during 16hours is 8*diet 
  # (8 * mets concentrations given as input, 8X over 16 hrs, i.e. approximately ~0.5X per hour)
  # *625 - convert mM (in diet table) to mM/cell (as in arena)
  for (AminoAcid in AA_diet_concentration$met_id) {
    if (AminoAcid %in% colnames(sublb_df)) {
      AminoAcid_diet_conc_mM <- diet_concentration * prop_gizzard * AA_diet_concentration$met_mM[AA_diet_concentration$met_id == AminoAcid]
      AminoAcid_diet_conc_mMcell <- AminoAcid_diet_conc_mM * 625 
      AminoAcid_diet_conc_mMcell_perHour <- 0.5 * AminoAcid_diet_conc_mMcell / 10^12
      
      # current concentration of AminoAcid in arena:
      AminoAcid_arena_conc <- mean(sublb_df[, AminoAcid])
      # leave conc of AminoAcid not higher than retention_fraction of conc_diet_perHour
      arena2 <- BacArena::addSubs(arena2, 
                                  smax = min(AminoAcid_arena_conc, retention_fraction * AminoAcid_diet_conc_mMcell_perHour),
                                  mediac = AminoAcid, 
                                  add = FALSE, 
                                  unit = "mmol/cell")
    }
  }
  return(arena2)
}

# Process fatty acid absorption (used in jejunum and ileum)
process_fatty_acid_absorption <- function(arena2, diet, retention_fractions, diet_concentration, prop_gizzard) {
  # REMOVE/ABSORB 4 TAGs: palmitic (hdca), linoleic (lnlc), oleic (ocdcea) and stearic (ocdca) acids
  sublb <- BacArena::getSublb(arena2)
  sublb_df <- as.data.frame(sublb)
  
  fatty_acids <- list(
    "EX_cpd00214_e0" = retention_fractions$palmitic,   # palmitic acid
    "EX_cpd15269_e0" = retention_fractions$oleic,      # oleic acid  
    "EX_cpd01122_e0" = retention_fractions$linoleic,    # linoleic acid
    "EX_cpd01080_e0" = retention_fractions$stearic     # stearic acid
  )
  
  for (fatty_acid_id in names(fatty_acids)) {
    if (fatty_acid_id %in% colnames(sublb_df)) {
      current_conc <- mean(sublb_df[, fatty_acid_id])
      retention_fraction <- fatty_acids[[fatty_acid_id]]
      
      arena2 <- apply_absorption(arena2, fatty_acid_id, current_conc, retention_fraction, 
                                 diet, diet_concentration, prop_gizzard)
    }
  }
  
  return(arena2)
}

################################################################################
####################           MAIN FUNCTIONS              ####################
################################################################################

AbsorptionConversion_duodenum <- function(arena2) {
  
  diet_concentration <- as.numeric(Parameters_table$diet_proportion[Parameters_table$Compartment == "Gizzard"])
  
  # Get current concentrations for key metabolites
  # starch_n19, starch_n27, glucose, fructose, maltose, sucrose, maltotetraose
  metabolites <- c("EX_cpd90003_e0", "EX_cpd90004_e0", "EX_cpd00027_e0", "EX_cpd00082_e0", 
                   "EX_cpd00179_e0", "EX_cpd00076_e0", "EX_cpd01262_e0")
  current_concs <- get_metabolite_concentrations(arena2, metabolites)
  
  # leave conc of starches not higher than 40% of conc_diet_perHour
  arena2 <- apply_absorption(arena2, "EX_cpd90003_e0", current_concs[["EX_cpd90003_e0"]], 0.4, diet, diet_concentration, prop_gizzard)
  arena2 <- apply_absorption(arena2, "EX_cpd90004_e0", current_concs[["EX_cpd90004_e0"]], 0.4, diet, diet_concentration, prop_gizzard)
  
  ##### ABSORB SUGARS
  # absorb/remove 80% of glucose and fructose: leave conc not higher than 20% of conc_diet_perHour
  arena2 <- apply_absorption(arena2, "EX_cpd00027_e0", current_concs[["EX_cpd00027_e0"]], 0.2, diet, diet_concentration, prop_gizzard)
  arena2 <- apply_absorption(arena2, "EX_cpd00082_e0", current_concs[["EX_cpd00082_e0"]], 0.2, diet, diet_concentration, prop_gizzard)
  
  ##### ABSORB AMINO ACIDS
  ## leave concentrations of AA not higher than 60% of ingested with the diet
  arena2 <- process_amino_acid_absorption(arena2, diet, 0.6, diet_concentration, prop_gizzard)
  
  return(arena2)
}

AbsorptionConversion_jejunum <- function(arena2, diet) {
  
  diet_concentration <- as.numeric(Parameters_table$diet_proportion[Parameters_table$Compartment == "Gizzard"])
  
  # Get current concentrations for key metabolites
  # starch_n19, starch_n27, glucose, fructose, maltose, sucrose, maltotetraose
  metabolites <- c("EX_cpd90003_e0", "EX_cpd90004_e0", "EX_cpd00027_e0", "EX_cpd00082_e0", 
                   "EX_cpd00076_e0", "EX_cpd00179_e0", "EX_cpd01262_e0")
  current_concs <- get_metabolite_concentrations(arena2, metabolites)
  
  ##### CONVERT STARCH
  # leave conc of starches not higher than 25% of conc_diet_perHour
  # prop_gizzard - proportion of diet added, 0.5 - because diet added every 2hrs,
  # 625 - conversion from mM to mmol/cell, 10^12 - from mmol/cell to fmol/cell 
  arena2 <- apply_absorption(arena2, "EX_cpd90003_e0", current_concs[["EX_cpd90003_e0"]], 0.25, diet, diet_concentration, prop_gizzard)
  arena2 <- apply_absorption(arena2, "EX_cpd90004_e0", current_concs[["EX_cpd90004_e0"]], 0.25, diet, diet_concentration, prop_gizzard)
  
  ##### ABSORB SUGARS
  # absorb/remove glucose and fructose and leave their conc not higher than 10% of conc_diet_perHour
  arena2 <- apply_absorption(arena2, "EX_cpd00027_e0", current_concs[["EX_cpd00027_e0"]], 0.1, diet, diet_concentration, prop_gizzard)
  arena2 <- apply_absorption(arena2, "EX_cpd00082_e0", current_concs[["EX_cpd00082_e0"]], 0.1, diet, diet_concentration, prop_gizzard)
  
  # digest 85% of maltose, malttr and sucrose ---> glucose (being absorbed) as well:
  arena2 <- apply_absorption(arena2, "EX_cpd00076_e0", current_concs[["EX_cpd00076_e0"]], 0.15, diet, diet_concentration, prop_gizzard)
  arena2 <- apply_absorption(arena2, "EX_cpd00179_e0", current_concs[["EX_cpd00179_e0"]], 0.15, diet, diet_concentration, prop_gizzard)
  
  ##### ABSORB AMINO ACIDS
  ## leave concentrations of AA not higher than 30% of ingested with the diet
  arena2 <- process_amino_acid_absorption(arena2, diet, 0.3, diet_concentration, prop_gizzard)
  
  ##### REMOVE/ABSORB FATTY ACIDS
  # leave concentrations of palmitic (hdca), oleic acid (ocdcea) not higher than 40% of ingested
  # stearic (ocdca) - not higher than 50% of ingested
  fatty_acid_retention <- list(palmitic = 0.4, oleic = 0.4,  linoleic = 0.4, stearic = 0.5)
  arena2 <- process_fatty_acid_absorption(arena2, diet, fatty_acid_retention, diet_concentration, prop_gizzard)
  
  return(arena2)
}

AbsorptionConversion_ileum <- function(arena2, diet) {
  
  diet_concentration <- as.numeric(Parameters_table$diet_proportion[Parameters_table$Compartment == "Gizzard"])
  
  # Get current concentrations for key metabolites
  # starch_n19, starch_n27, glucose, fructose, maltose, sucrose, maltotetraose
  metabolites <- c("EX_cpd90003_e0", "EX_cpd90004_e0", "EX_cpd00027_e0", "EX_cpd00082_e0", 
                   "EX_cpd00076_e0", "EX_cpd00179_e0")
  current_concs <- get_metabolite_concentrations(arena2, metabolites)
  
  ##### CONVERT STARCH
  # adjust starch1-2 conc to 10% of diet intake:
  arena2 <- apply_absorption(arena2, "EX_cpd90003_e0", current_concs[["EX_cpd90003_e0"]], 0.1, diet, diet_concentration, prop_gizzard)
  arena2 <- apply_absorption(arena2, "EX_cpd90004_e0", current_concs[["EX_cpd90004_e0"]], 0.1, diet, diet_concentration, prop_gizzard)
  
  ###### ABSORB SUGARS
  # absorb/remove 95% of all ingested glucose, fructose: leave only 5%
  arena2 <- apply_absorption(arena2, "EX_cpd00027_e0", current_concs[["EX_cpd00027_e0"]], 0.05, diet, diet_concentration, prop_gizzard)
  arena2 <- apply_absorption(arena2, "EX_cpd00082_e0", current_concs[["EX_cpd00082_e0"]], 0.05, diet, diet_concentration, prop_gizzard)
  
  # digest 90% of maltose and sucrose ---> glucose (being absorbed) as well: leave only 10%
  arena2 <- apply_absorption(arena2, "EX_cpd00076_e0", current_concs[["EX_cpd00076_e0"]], 0.1, diet, diet_concentration, prop_gizzard)
  arena2 <- apply_absorption(arena2, "EX_cpd00179_e0", current_concs[["EX_cpd00179_e0"]], 0.1, diet, diet_concentration, prop_gizzard)
  
  ##### ABSORB AMINO ACIDS
  ## leave concentrations of AA not higher than 20% of ingested with the diet
  arena2 <- process_amino_acid_absorption(arena2, diet, 0.2, diet_concentration, prop_gizzard, category_filter = "Proteins")
  
  ##### REMOVE/ABSORB FATTY ACIDS
  # leave concentrations of palmitic (hdca), oleic acid (ocdcea) not higher than 20% of ingested
  # stearic (ocdca) - not higher than 30% of ingested
  fatty_acid_retention <- list(palmitic = 0.2, oleic = 0.2, linoleic = 0.2, stearic = 0.3)
  arena2 <- process_fatty_acid_absorption(arena2, diet, fatty_acid_retention, diet_concentration, prop_gizzard)
  
  return(arena2)
}