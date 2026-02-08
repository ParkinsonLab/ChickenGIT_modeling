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

CreateTypeSpecieAssoc_table <- function(arena) {
  # arena <- getArena(sim_object, 1)
  specs_names <- c()
  for (i in 1:length(arena@models)) {
    specs_names <- c(specs_names, arena@models[[i]]@mod_desc) 
  }
  type_specie_table <- as.data.frame(cbind(unique(specs_names), unique(sort(arena@orgdat$type))))
  colnames(type_specie_table) <- c("Name", "Type")
  return(type_specie_table)
}


RemoveBacMetabolites <- function(sim_object, arena2_object, perc_to_remove_bac, 
                                 perc_to_remove_mets = 0, type_specie_table) {
  # if perc_to_remove_mets not specified (i.e. in older scripts) set it equal to perc_to_remove_bac
  if (perc_to_remove_mets == 0) {
    perc_to_remove_mets <- perc_to_remove_bac
  }
  # For cecum (no absorption/conversion there), for ileum - arena2 already after absorption function
  if (missing(arena2_object) ) {
    arena2 <- BacArena::getArena(sim_object,1) 
  } else { arena2 <- arena2_object }
  
  # select perc_to_remove_bac % of the grid randomly that will be removed in the next iteration and transferred to cecum
  rm_orgdat <- arena2@orgdat[sample(nrow(arena2@orgdat), 
                                    round(perc_to_remove_bac*nrow(arena2@orgdat))), c("x", "y", "type") ]
  
  
  # connect the numeric "type" in orgdat with the name of corresponding specie:
  # problem here: when some species get removed completely from the grid. solution: create and type_specie_table after 1st iteration
  type_specie <- subset(type_specie_table, type_specie_table$Type %in% unique(arena2@orgdat$type))
  # another problem: when species not present in comp initially get transferred from other compartment, so they aren't in type_specie_table
  # solution: add this bacteria to type_specie table separately:
  new_bactype_arena2 <- setdiff(unique(arena2@orgdat$type), type_specie_table$Type)
  if (length(new_bactype_arena2) > 0) {
    for (bactype in new_bactype_arena2) {
      type_specie_new <- data.frame(cbind(arena2@specs[[bactype]]@model@mod_desc, bactype))
      colnames(type_specie_new) <- c("Name", "Type")
      type_specie <- rbind(type_specie, type_specie_new)
    }
  }
  
  # and create a table with numbers of each specie individuals being removed / transferred to compartment N+1
  species_numbers_rm <- data.frame(table(rm_orgdat$type)) # Var 1 = type, Freq = how many individuals of this type
  colnames(species_numbers_rm) <- c("Type", "Number")
  species_rm <- dplyr::inner_join(species_numbers_rm, type_specie , by = "Type") # this table will be used to add to compartment N+1
  
  # get substance matrix (concentrations of each metabolite in each occupied grid cell)
  sublb <- BacArena::getSublb(arena2)
  sublb <- as.data.frame(sublb)
  # substances in grid cells to be deleted in the next iteration from compartment N and transferred to compartment N+1:
  # IMPORTANT: do this BEFORE deleting bacteria (from orgdat), as the substances (might*) get deleted as well!
  if (perc_to_remove_bac == perc_to_remove_mets) {
    sublb_rm <- plyr::match_df(sublb, rm_orgdat) # before introducing perc_to_remove_mets
  } else if (perc_to_remove_bac < perc_to_remove_mets){
    # FIRST, save sublb_rm_1 as metabolites in rm_orgdat cells (they will be removed together with organisms from these cells)
    sublb_rm_1 <- plyr::match_df(sublb, rm_orgdat)
    # save remaining metabolites after organisms removal sublb_upd:
    sublb_upd <- dplyr::setdiff(sublb, sublb_rm_1)
    # SECONDLY, save sublb_rm_2 which will consist of metabolites from the % of randomly sampled cells from
    # sublb_upd =   (%perc_to_remove_mets - %perc_to_remove_bac)*nrow(sublb) 
    # NOT sublb_upd, because number of rows in sublb_upd is already less by %perc_to_remove_bac compared with initial sublb:
    sublb_rm_2 <- sublb_upd[sample(nrow(sublb_upd), 
                                   round((perc_to_remove_mets - perc_to_remove_bac)*nrow(sublb))), ] # with new parameter perc_to_remove_mets
    # FINALLY, combine sublb_rm_1 and sublb_rm_2 into sublb_rm - mets removed from the grid after iteration:
    sublb_rm <- rbind(sublb_rm_1, sublb_rm_2)
    
  } else if (perc_to_remove_bac > perc_to_remove_mets){
    # FIRST, save sublb_rm_1 as metabolites in rm_orgdat cells (they are removed together with organisms from these cells)
    sublb_rm_1 <- plyr::match_df(sublb, rm_orgdat)
    # since % of metabolites to be removed should have been smaller than what was removed (perc_to_remove_mets < perc_to_remove_bac)
    # we have to bring back part of removed sublb_rm_1, specifically:  (%perc_to_remove_bac - %perc_to_remove_mets)*nrow(sublb_rm_1) 
    sublb_rm_1_back <- sublb_rm_1[sample(nrow(sublb_rm_1), 
                                         round((perc_to_remove_bac - perc_to_remove_mets)*nrow(sublb_rm_1))), ] 
    
    # therefore, sublb_rm - mets removed from the grid after iteration - are:
    sublb_rm <- dplyr::setdiff(sublb_rm_1, sublb_rm_1_back)
    
  }
 
  # update arena2 (remove those pre-selected random perc_to_remove_bac %)
  # keep in mind that sublb_rm_1 (might* ?) be removed together with organisms in corresponding cells (rm_orgdat)
  for (row in 1:nrow(rm_orgdat)) {
    arena2@orgdat <- arena2@orgdat[-which((arena2@orgdat$x == rm_orgdat[row, 1]) &
                                            (arena2@orgdat$y == rm_orgdat[row, 2])), ]
  }
  
  
  # update all metabolites concentrations in emptied cells: set them to zero;
  # then re-mix the environment (stirEnv) to evenly distribute mets
  sublb_remaining <- dplyr::setdiff(sublb, sublb_rm)
  sublb_rm_zeros <- plyr::match_df(sublb, sublb_rm)
  sublb_rm_zeros[, 3:ncol(sublb_rm_zeros)] <- 0
  sublb_new  <- rbind(sublb_remaining, sublb_rm_zeros)
  sublb_new <- stirEnv(arena2, sublb_new)
  arena2@sublb <- as.matrix(sublb_new)
  
  
  # for some reason after manipulations arena2@models slot becomes empty and the following sim_object@models as well, so:
  arena2@models <- sim_object@models
  
  return(list("arena2" = arena2, "rm_orgdat" = rm_orgdat, "sublb_rm" = sublb_rm, "species_rm" = species_rm))
}
