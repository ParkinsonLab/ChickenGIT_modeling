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


AddBacDiet_gizzard <- function(Chicken_organisms_gizzard, model_list_gizzard, sampleName, diet_mets, arena_gizzard,
                               prop_gizzard, diet_added_perc, initial_commsize_gizzard, econc, all_conc, essential_mets) {
  
  diet_concentration <- as.numeric(Parameters_table$diet_proportion[Parameters_table$Compartment == "Gizzard"])
  
  model_list <- model_list_gizzard
  TableRelativeAbundance <- as.data.frame(Get16SAbundance_table(OTU_table, Day = "D10", Compartment = "Gizzard"))
  ########## ADDING 20% OF INITIAL ARENA (BACTERIA) TO GIZZARD #########
  # essential_meds <- c()
  for (taxgroup in unique(Chicken_organisms_gizzard$TaxCategory)) {
    # if sample is missing - reconstruct its rel. abundance as an average of rel. abundances in samples from the same day+diet
    if (!(sampleName %in% TableRelativeAbundance$variable)) {
      abundance = mean(TableRelativeAbundance$value[TableRelativeAbundance$taxonomy %in% taxgroup & 
                                                      TableRelativeAbundance$Diet == stringr::str_remove(gsub('[0-9]','', sampleName), '-')]) /
        length(Chicken_organisms_gizzard$Specie[Chicken_organisms_gizzard$TaxCategory %in% taxgroup]) # divide by the number of organisms
      # chosen to represent this taxgroup
    }
    abundance = TableRelativeAbundance$value[(TableRelativeAbundance$taxonomy %in% taxgroup) &
                                               (TableRelativeAbundance$variable %in% sampleName)] /
      length(Chicken_organisms_gizzard$Specie[Chicken_organisms_gizzard$TaxCategory %in% taxgroup]) # divide by the number of organisms
    # chosen to represent this taxgroup
    
    # get the coordinates of unoccupied (free) grid cells:
    y <- c(); for (x in 1:100) {y <- c(y, rep(x,100))} 
    grid_coord <- data.frame(cbind(x = rep(c(1:100),100), y = y)) # created df with coordinates x,y - (1:100)x(1:100)
    grid_coord_unoccupied <- dplyr::setdiff(grid_coord, arena_gizzard@orgdat[, c("x", "y")]) # get coordinates of unoccupied (free) cells
    
    for (model_name in Chicken_organisms_gizzard$Specie[Chicken_organisms_gizzard$TaxCategory %in% taxgroup]) {
      for (i in 1:length(model_list)) {
        # find necessary models in the ModelList and upload them through Bac function
        # and then add to arena in the amount defined above (abundance)
        if (model_list[[i]]@mod_name == model_name) {
          model_list[[i]] <- RenameModelReactions(model_list[[i]])
          # skip those bacs which abundance*0.2*commsize < 0.1
          if ((ceiling(abundance*initial_commsize_gizzard)>0)) { 
            bac = BacArena::Bac(model=model_list[[i]], growtype="exponential")
            
            # check if bac name is already used by other organism to avoid error:
            # "Organism of the same type but with different model already present, added a new one"
            idx.dupl <- which(names(arena_gizzard@specs) == bac@type)  # returns an empty vector if no bac@type among arena_gizzard@specs
            if( length(idx.dupl) > 0 ){
              if( !identical( bac@model, arena_gizzard@specs[[idx.dupl]]@model ) ){
                bac@model <- arena_gizzard@specs[[idx.dupl]]@model
              } 
            }
            # adding 20% of initial bac number back if 0.2*that number doesn't exceed number of onoccupied cells:
            number_of_bac_to_add <- min(ceiling(abundance*0.2*initial_commsize_gizzard), floor(abundance*nrow(grid_coord_unoccupied) ))
            if (number_of_bac_to_add > 0) {
              arena_gizzard = BacArena::addOrg(arena_gizzard, bac, amount=number_of_bac_to_add )
              # update grid_coord_unoccupied after filling some part of free grid with bac
              grid_coord_unoccupied <- dplyr::setdiff(grid_coord_unoccupied, arena_gizzard@orgdat[, c("x", "y")])
               }
            
            
           }
        }
      }
    }
  }
  
  ########## ADDING % OF INITIAL DIET (METABOLITES) ################
  # add all mets with low concentrations
  # arena_gizzard <- BacArena::addSubs(arena_gizzard, smax = all_conc,  unit = "mM")
  # add only those essential metabolites that aren't present in the diet:
  # essential_meds_nondiet <- setdiff(essential_meds, as.character(diet_mets$met_id))
  essential_meds_nondiet <- setdiff(essential_mets, as.character(diet_mets$met_id))
  # add the essential diets on top
  arena_gizzard <- BacArena::addSubs(arena_gizzard,smax = diet_added_perc*econc,unit="mM",
                                     mediac= essential_meds_nondiet, add=T, addAnyway = T)
  # add metabolite fluxes generated bu diet design (vmh)
  # diet_added_perc - how much added every X hrs; prop_gizzard * diet_concentration - fraction of 100g diet
  # diet_mets$met_mM <- prop_gizzard * diet_added_perc * diet_mets$met_mM
  arena_gizzard <- BacArena::addSubs(arena_gizzard,smax =  prop_gizzard * diet_concentration *  diet_added_perc * diet_mets$met_mM,
                                     unit="mM",
                                     mediac = as.character(diet_mets$met_id), add = T, addAnyway = T)
  
  arena_gizzard_to_return <- arena_gizzard
  return(arena_gizzard_to_return)
}

AddBacDiet_gizzard_random <- function(Chicken_organisms_gizzard, model_list_gizzard, sampleName, diet_mets, arena_gizzard,
                               prop_gizzard, diet_added_perc, initial_commsize_gizzard, econc, all_conc, relabundance_arg,
                               abundance_gizzard_df=NULL) {
  
  diet_concentration <- as.numeric(Parameters_table$diet_proportion[Parameters_table$Compartment == "Gizzard"])
  
  model_list <- model_list_gizzard
  TableRelativeAbundance <- as.data.frame(Get16SAbundance_table(OTU_table, Day = "D10", Compartment = "Gizzard"))
  # make a table with taxonomic categories present in a sample
  TableRelativeAbundance_sample_presence <- TableRelativeAbundance[which(TableRelativeAbundance$value > 0.001 &
                                                                           TableRelativeAbundance$variable == sampleName),]
  
  # then subset the models from model_list to be used for the respective sample's community to keep only the ones present in this sample
  Chicken_organisms_gizzard <- Chicken_organisms_gizzard[which(Chicken_organisms_gizzard$TaxCategory %in% TableRelativeAbundance_sample_presence$taxonomy), ] 
  
  N_species <- nrow(Chicken_organisms_gizzard)
  
  # bring back the random relabundance that was used to build gizzard's arena at timestep=0
  if (relabundance_arg == "random") {
    abundance_df <- abundance_gizzard_df
  }
  
  ########## ADDING 20% OF INITIAL ARENA (BACTERIA) TO GIZZARD #########
  essential_meds <- c()
  
  # get the coordinates of unoccupied (free) grid cells:
  y <- c(); for (x in 1:100) {y <- c(y, rep(x,100))} 
  grid_coord <- data.frame(cbind(x = rep(c(1:100),100), y = y)) # created df with coordinates x,y - (1:100)x(1:100)
  grid_coord_unoccupied <- dplyr::setdiff(grid_coord, arena_gizzard@orgdat[, c("x", "y")]) # get coordinates of unoccupied (free) cells
  
  for (model_name in Chicken_organisms_gizzard$Specie) {
    
    if (relabundance_arg == "equal") {
      abundance <- 1/N_species 
    } 
    if (relabundance_arg == "random") {
      abundance <- as.numeric(as.character(abundance_df$relabundance[abundance_df$Specie == model_name]))
    }
    
    for (i in 1:length(model_list)) {
      # find necessary models in the ModelList and upload them through Bac function
      # and then add to arena in the amount defined above (abundance)
      if (model_list[[i]]@mod_name == model_name) {
        model_list[[i]] <- RenameModelReactions(model_list[[i]])
        # skip those bacs which abundance*0.2*commsize < 0.1
        if ((ceiling(abundance*initial_commsize_gizzard)>0)) { 
          bac = BacArena::Bac(model=model_list[[i]], growtype="exponential")
          
          # check if bac name is already used by other organism to avoid error:
          # "Organism of the same type but with different model already present, added a new one"
          idx.dupl <- which(names(arena_gizzard@specs) == bac@type)  # returns an empty vector if no bac@type among arena_gizzard@specs
          if( length(idx.dupl) > 0 ){
            if( !identical( bac@model, arena_gizzard@specs[[idx.dupl]]@model ) ){
              bac@model <- arena_gizzard@specs[[idx.dupl]]@model
            } 
          }
          # adding 20% of initial bac number back if 0.2*that number doesn't exceed number of onoccupied cells:
          number_of_bac_to_add <- min(ceiling(abundance*0.2*initial_commsize_gizzard), floor(abundance*nrow(grid_coord_unoccupied) ))
          if (number_of_bac_to_add > 0) {
            arena_gizzard = BacArena::addOrg(arena_gizzard, bac, amount=number_of_bac_to_add )
            # update grid_coord_unoccupied after filling some part of free grid with bac
            grid_coord_unoccupied <- dplyr::setdiff(grid_coord_unoccupied, arena_gizzard@orgdat[, c("x", "y")])
            
            # to get essential metabolites list for this bacteria, create empty arena with this bac only:
            arena_empty <- BacArena::Arena(100, 100, stir = T, Lx=0.025, Ly=0.025, tstep=1, seed = 239)
            arena_empty <- BacArena::addOrg(arena_empty, bac, amount=ceiling(abundance*0.2*initial_commsize_gizzard))
          }
          
          
          # save essential metabolites for each bac-organism
          # only_return=T if essential metabolites should only be returned but not added to arena
          essential_meds <- unique(c(essential_meds, BacArena::addEssentialMed(arena_empty, bac, only_return = T)))
        }
      }
    }
  }

  
  ########## ADDING % OF INITIAL DIET (METABOLITES) ################
  # add all mets with low concentrations
  # arena_gizzard <- BacArena::addSubs(arena_gizzard, smax = all_conc,  unit = "mM")
  # add only those essential metabolites that aren't present in the diet:
  essential_meds_nondiet <- setdiff(essential_meds, as.character(diet_mets$met_id))
  # add the essential diets on top
  arena_gizzard <- BacArena::addSubs(arena_gizzard,smax = diet_added_perc*econc,unit="mM",
                                     mediac= essential_meds_nondiet, add=T, addAnyway = T)
  # add metabolite fluxes generated bu diet design (vmh)
  # diet_added_perc - how much added every X hrs; prop_gizzard * diet_concentration - fraction of 100g diet
  # diet_mets$met_mM <- prop_gizzard * diet_added_perc * diet_mets$met_mM
  arena_gizzard <- BacArena::addSubs(arena_gizzard,smax =  prop_gizzard * diet_concentration *  diet_added_perc * diet_mets$met_mM,
                                     unit="mM",
                                     mediac = as.character(diet_mets$met_id), add = T, addAnyway = T)
  
  arena_gizzard_to_return <- arena_gizzard
  return(arena_gizzard_to_return)
}



AddBacMets_fromUpstreamCompartment <- function(removefunction_result_n0,
                                               sim_object_n1, arena2_n1, model_list_n1, model_list_n0) {
  
  # get substances, species and orgdat (positions with species) removed from compartment N:
  rm_orgdat_n0 <- removefunction_result_n0[[2]]
  sublb_n0_rm <- removefunction_result_n0[[3]]
  species_rm_n0 <- removefunction_result_n0[[4]]
  # get substance matrix (concentrations of each metabolite in each occupied cell on a compartment N+1' grid)
  sublb_n1 <- BacArena::getSublb(arena2_n1)
  sublb_n1 <- as.data.frame(sublb_n1)
  # get the coordinates of unoccupied (free) grid cells:
  y <- c(); for (x in 1:100) {y <- c(y, rep(x,100))} 
  grid_coord <- data.frame(cbind(x = rep(c(1:100),100), y = y)) # created df with coordinates x,y - (1:100)x(1:100)
  grid_coord_unoccupied <- dplyr::setdiff(grid_coord, arena2_n1@orgdat[, c("x", "y")]) # get coordinates of unoccupied (free) cells
  
  
  
  ### add bacteria deleted from compartment N in this iteration to the compartment N+1 grid: 
  # subsample the free space on the grid to get sub-grid of the size of rm_orgdat_n0 (or all of the free space, 
  # whatever is smaller) and fill it with compartment N community and metabolites to known positions (x,y):
  add_orgdat_n1 <- grid_coord_unoccupied[sample(nrow(grid_coord_unoccupied), min(nrow(rm_orgdat_n0), nrow(grid_coord_unoccupied))), ] 
  
  model_list <- model_list_n0
  # go through each specie being transferred from n0 to n1
  for (j in 1:nrow(species_rm_n0)) {
    for (i in 1:length(model_list)) {
      if (nrow(add_orgdat_n1) > 1) {
        # find necessary models in the ModelList and upload them through Bac function
        # and then add to arena in the amount defined above (abundance)
        if (model_list[[i]]@mod_desc == as.character(species_rm_n0$Name[j])) {
          model_list[[i]] <- RenameModelReactions(model_list[[i]])
          number = species_rm_n0$Number[j] # how many individuals of this specie to add
          positions_to_fill <- add_orgdat_n1[sample(nrow(add_orgdat_n1), min(number, nrow(add_orgdat_n1))), ]
          bac = BacArena::Bac(model=model_list[[i]], growtype="exponential")
          
          # check if bac name is already used by other organism to avoid "Organism of the same type but with different model already present, added a new one"
          idx.dupl <- which(names(arena2_n1@specs) == bac@type)  # returns an empty vector if no bac@type among arena2_n1@specs
          if( length(idx.dupl) > 0 ){
            if( !identical( bac@model, arena2_n1@specs[[idx.dupl]]@model ) ){
              bac@model <- arena2_n1@specs[[idx.dupl]]@model
            } 
          }
          
          arena2_n1 <- BacArena::addOrg(arena2_n1, bac, amount = number, x = positions_to_fill$x, y = positions_to_fill$y)
          # now remove positions from add_orgdat_n1 that have been just filled with individuals so that other species from
          # species_rm_n0 list won't be added to the same grid cell
          add_orgdat_n1 <- dplyr::setdiff(add_orgdat_n1, positions_to_fill)
        }
      }
    }
  }
  
  ####  ADD METABOLITES FROM COMPARTMENT N: UPDATED AUGUST 2023, ANOTHER APPROACH via addSubs #####
  # get a mean of all concentrations from sublb_n0_rm:
  sublb_n0_rm_meanconc <- colMeans(sublb_n0_rm[, c(-1,-2)])
  # then get mean of all concentrations from sublb_n1:
  sublb_n1_meanconc <- colMeans(sublb_n1[, c(-1,-2)])
  # the resulting concentration after adding sublb_n0_rm to sublb_n1 for metabolite i would be:
  # (sublb_n0_rm_meanconc[i]*nrow(sublb_n0_rm) + sublb_n1_meanconc[i]*nrow(sublb_n1)) / (nrow(sublb_n0_rm) + nrow(sublb_n1))
  # add this concentrations via addSubs using add=F to replace the existing concentration value with the new one! 
  # if the metabolite from sublb_n0_rm isn't present in the n1 compartment, simply consider sublb_n1_meanconc[metabolite] = 0
  for (met in names(sublb_n0_rm_meanconc)) {
    if (met %in% names(sublb_n1_meanconc)) {
      new_conc <- (sublb_n0_rm_meanconc[met]*nrow(sublb_n0_rm) + sublb_n1_meanconc[met]*nrow(sublb_n1)) / (nrow(sublb_n0_rm) + nrow(sublb_n1))
    } else {
      new_conc <- (sublb_n0_rm_meanconc[met]*nrow(sublb_n0_rm) + 0) / (nrow(sublb_n0_rm) + nrow(sublb_n1))
    }
    new_conc <- as.numeric(new_conc)
    # add=F to replace the existing concentration value with the new one
    # by default the conc value is considered in fmol/cell => divide by 10^12 for that to be mmol/cell
    arena2_n1 <- BacArena::addSubs(arena2_n1, smax = new_conc/10^12, mediac = met, unit = "mmol/cell", add = FALSE, addAnyway = T)
  }
  

  return(arena2_n1)
}

AddBacMets_fromDownstreamCompartment <- function(removefunction_result_n1,
                                                 sim_object_n0, arena2_n0,
                                                 model_list_n0, model_list_n1){
  
  arena2_n0_upd <- AddBacMets_fromUpstreamCompartment(removefunction_result_n0 = removefunction_result_n1,
                                                      sim_object_n1=sim_object_n0, arena2_n1=arena2_n0,
                                                      model_list_n1=model_list_n0, model_list_n0=model_list_n1)
  return(arena2_n0_upd)
}