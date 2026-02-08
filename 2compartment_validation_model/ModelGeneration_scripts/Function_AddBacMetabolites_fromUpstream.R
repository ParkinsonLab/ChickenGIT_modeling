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



AddBacMets_fromUpstreamCompartment <- function(removefunction_result_n0,
                                               sim_object_n1, arena2_n1, model_list_n0) {
  
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
              # Instead of overwriting model, rename the transferred species
              # bac@type <- paste0(bac@type, "_transferred")
              # cat("Renamed transferred species to:", bac@type, "\n")
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
  
  ####  ADD METABOLITES FROM COMPARTMENT N: 
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
  
  # Force rebuild @models from @specs (as some samples fail to keep the @models slot)
  arena2_n1@models <- lapply(arena2_n1@specs, function(spec) spec@model)
  names(arena2_n1@models) <- names(arena2_n1@specs)
  
  return(arena2_n1)
}

