library(dplyr)
library(magrittr)
library(BacArena)
library(glpkAPI)
library(Rcpp)
library(reshape2)
library(libSBML)
library(sybilSBML)
library(stringr)



RenameModelReactions = function(model) {
  model@react_id <- gsub('__', '_',model@react_id) 
  model@met_id <- gsub('__', '_',model@met_id) 
  model@met_id <- gsub('\\[', '(',model@met_id) 
  model@met_id <- gsub('\\]', ')',model@met_id) 
  return(model)
}

##############################################################################################################
####################       FUNCTION TO CONSTRUCT ARENA FOR CECUM. MIXING, DIET   #############################
##############################################################################################################

CreateModelArena = function(model_list, bird_id, compartment, 
                            diet, diet_proportion, initial_commsize=500, econc, 
                            stir_var = T, arena_size=100, speed=5, tstep=1, seed=239,
                            metadata_file) {
  
  otu_mag_metadata_paired <- metadata_file
  
  # Parse bird_id to find the actual sample names
  # bird_id format: G1_1 -> need to find G1_I1 (ileum) and G1_T1 (cecum)
  group_id <- stringr::str_extract(bird_id, "^[^_]+")
  bird_number <- stringr::str_extract(bird_id, "\\d+$")
  
  if (compartment == "Ileum") {
    target_sample <- paste0(group_id, "_I", bird_number)
  } else if (compartment == "Cecum") {
    target_sample <- paste0(group_id, "_T", bird_number)
  } 
  
  # Get abundance data for this specific sample
  abundance_data <- otu_mag_metadata_paired %>%
    dplyr::filter(sample_name == target_sample, total_rel_abundance > 0) %>%
    dplyr::select(mag_id_clean, total_rel_abundance, taxonomy, species_ID) %>%
    dplyr::rename(MAG_ID = mag_id_clean, Relative_Abundance = total_rel_abundance) 
  
  if (nrow(abundance_data) == 0) {
    stop(paste("No abundance data found for", target_sample))
  }
  
  cat("Creating arena for", bird_id, "-", compartment, "(sample:", target_sample, ") with", 
      length(model_list), "MAGs\n")
  
  arena <- BacArena::Arena(n=arena_size, m=arena_size,
                           stir = stir_var, Lx=0.025*(arena_size/100), Ly=0.025*(arena_size/100),
                           tstep = tstep,
                           seed = seed)
  
  essential_meds <- c()
  for (i in 1:length(model_list)) {
    # find necessary models in the ModelList and upload them through Bac function
    # and then add to arena in the amount defined above (abundance)
    model_list[[i]] <- RenameModelReactions(model_list[[i]])
    abundance = abundance_data$Relative_Abundance[abundance_data$species_ID == model_list[[i]]@mod_name]
    
    if (ceiling(abundance*initial_commsize)>0){
      
      bac = BacArena::Bac(model=model_list[[i]], growtype="exponential")
      arena = BacArena::addOrg(arena, bac, amount=ceiling(abundance*initial_commsize))
      
      # save essential metabolites for each bac-organism
      # 'only_return=T' if essential metabolites should only be returned but not added to arena
      essential_meds <- unique(c(essential_meds,
                                 BacArena::addEssentialMed(arena, bac, limit = 100, only_return = T)))
      
    }
  }
  
  # remove acetate and propionate from essential metabolites if present!
  essential_meds <- essential_meds[which(!(essential_meds %in% c("EX_cpd00029_e0", "EX_cpd00141_e0")))]
  
  # save only essential metabolites that aren't in the diet:
  essential_meds <- essential_meds[which(!(essential_meds %in% diet$met_id))]
  
  # add essential metabolites 
  if (length(essential_meds) > 0) {
    arena <- BacArena::addSubs(arena,smax=econc,unit="mM", mediac=essential_meds,add=F, addAnyway = T)
  }
  
  # add nucleobases and gam to essential mets
  # gapseq: add cytosine, uracil, cytidine, adenosine, guanosine and CoA
  arena <- BacArena::addSubs(arena,smax=econc,unit="mM", 
                             mediac=c("EX_cpd00307_e0",  "EX_cpd00092_e0","EX_cpd00367_e0",
                                      "EX_cpd00182_e0", "EX_cpd00311_e0","EX_cpd00010_e0","EX_cpd00276"), add=T, addAnyway = T)
  
  essential_meds <- unique(c(essential_meds, c("EX_cpd00307_e0",  "EX_cpd00092_e0","EX_cpd00367_e0",
                                               "EX_cpd00182_e0", "EX_cpd00311_e0","EX_cpd00010_e0","EX_cpd00276")))
  
  # add metabolite concentrations from input diet
  arena <- BacArena::addSubs(arena,
                             smax=diet_proportion * diet$met_mM,
                             unit="mM", mediac=diet$met_id, 
                             add = F, addAnyway = T)
  
  
  return(list(arena, essential_meds))
}
