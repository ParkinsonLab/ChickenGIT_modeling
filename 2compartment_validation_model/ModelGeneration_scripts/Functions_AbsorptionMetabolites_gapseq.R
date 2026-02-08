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
# diet_concentration * 625 * diet$met_mM[diet$met_id == X]/10^12
apply_absorption <- function(arena2, met_id, current_conc, retention_fraction, diet_ileum, diet_concentration) {
  if (length(diet_ileum$met_mM[diet_ileum$met_id == met_id]) == 0) {
    # If metabolite not in diet, just return arena unchanged
    return(arena2)
  }
  
  # Calculate target concentration: retention_fraction of diet intake per hour
  #  625 - conversion from mM to mmol/cell, 10^12 - from mmol/cell to fmol/cell
  target_conc <- retention_fraction * diet_concentration * 625 * diet_ileum$met_mM[diet_ileum$met_id == met_id] / 10^12
  
  # Set concentration to minimum of current or target (absorption effect)
  arena2 <- BacArena::addSubs(arena2, 
                              smax = min(current_conc, target_conc),
                              mediac = met_id, 
                              add = FALSE, 
                              unit = "mmol/cell")
  return(arena2)
}

# Process amino acid absorption (used in all three compartments)
process_amino_acid_absorption <- function(arena2, diet_ileum, retention_fraction, diet_concentration,
                                          category_filter = c("Proteins", "Proteins_shadow")) {
  sublb <- BacArena::getSublb(arena2)
  sublb_df <- as.data.frame(sublb)
  
  # Get amino acids from diet
  AA_diet_concentration <- subset(diet_ileum, diet_ileum$category %in% category_filter)
  
  # diet_concentration - approximately how much diet-metabolites is ingested per hour 
  # *625 - convert mM (in diet table) to mM/cell (as in arena)
  for (AminoAcid in AA_diet_concentration$met_id) {
    if (AminoAcid %in% colnames(sublb_df)) {
      AminoAcid_diet_conc_mM <- diet_concentration  * AA_diet_concentration$met_mM[AA_diet_concentration$met_id == AminoAcid]
      AminoAcid_diet_conc_mMcell_perHour <- AminoAcid_diet_conc_mM * 625 / 10^12

      
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
process_fatty_acid_absorption <- function(arena2, diet_ileum, retention_fractions, diet_concentration) {
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
                                 diet_ileum, diet_concentration)
    }
  }
  
  return(arena2)
}

################################################################################
####################           MAIN FUNCTION                ####################
################################################################################

# rescaling the ileal absorption factors using the ileum diet instead of full diet (original input for 6-compartment model)
# F - Fraction of the original (full) feed we want to remain after ileum
# S - Fraction of the original feed that enters ileum - in the ileum diet file
# R=F/S - factor to be applied for absorption below
# Class	                         S     F  	 R = F/S  
# Starch (total)	              0.15	0.03	0.20
# Glucose / fructose          	0.05	0.02	0.40
# Maltose, sucrose, etc.      	0.10	0.03	0.33
# Amino acids (total)         	0.25	0.10	0.4
# Palmitic                     	0.25	0.15	0.60
# Stearic	                      0.40	0.30	0.75

AbsorptionConversion_ileum <- function(arena2, diet_ileum,  diet_proportion) {
  
   # Get current concentrations for key metabolites
  # starch_n19, starch_n27, glucose, fructose,  sucrose, maltose
  metabolites <- c("EX_cpd90003_e0", "EX_cpd90004_e0", "EX_cpd00027_e0", "EX_cpd00082_e0", 
                   "EX_cpd00076_e0", "EX_cpd00179_e0")
  current_concs <- get_metabolite_concentrations(arena2, metabolites)
  
  ##### STARCH
  arena2 <- apply_absorption(arena2, "EX_cpd90003_e0", current_concs[["EX_cpd90003_e0"]], 0.2, diet_ileum, diet_proportion)
  arena2 <- apply_absorption(arena2, "EX_cpd90004_e0", current_concs[["EX_cpd90004_e0"]], 0.2, diet_ileum, diet_proportion)
  
  ###### ABSORB SUGARS
  # absorb/remove 95% of all ingested glucose, fructose: leave only 5%
  arena2 <- apply_absorption(arena2, "EX_cpd00027_e0", current_concs[["EX_cpd00027_e0"]], 0.4, diet_ileum, diet_proportion)
  arena2 <- apply_absorption(arena2, "EX_cpd00082_e0", current_concs[["EX_cpd00082_e0"]], 0.4, diet_ileum, diet_proportion)
  
  # digest 90% of maltose and sucrose ---> glucose (being absorbed) as well: leave only 10%
  arena2 <- apply_absorption(arena2, "EX_cpd00076_e0", current_concs[["EX_cpd00076_e0"]], 0.33, diet_ileum, diet_proportion)
  arena2 <- apply_absorption(arena2, "EX_cpd00179_e0", current_concs[["EX_cpd00179_e0"]], 0.33, diet_ileum, diet_proportion)
  
  ##### ABSORB AMINO ACIDS
  arena2 <- process_amino_acid_absorption(arena2, diet_ileum, 0.4, diet_proportion, category_filter = "Proteins")
  
  ##### REMOVE/ABSORB FATTY ACIDS
  # leave concentrations of palmitic (hdca), oleic acid (ocdcea) not higher than 20% of ingested
  # stearic (ocdca) - not higher than 30% of ingested
  fatty_acid_retention <- list(palmitic = 0.6, oleic = 1, linoleic = 1, stearic = 0.75)
  arena2 <- process_fatty_acid_absorption(arena2, diet_ileum, fatty_acid_retention, diet_proportion)
  
  return(arena2)
}