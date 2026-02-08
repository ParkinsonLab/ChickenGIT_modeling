library(BacArena)
library(parallel)
library(glpkAPI)
library(Rcpp)
library(igraph)
library(reshape2)
library(libSBML)
library(sybilSBML)
library(ggplot2)
library(hash)
library(scales)


Compartments_models_noAGP_D10_gapseq <- readxl::read_xlsx("/Users/irina/Desktop/Study/UofT/ParkinsonLab/Modeling_chickens/ModelGeneration_Files/SpecieModels_present_in_compartments_D10_noAGP_2022upd.xlsx",
                                                          sheet = "Gapseq_ModelList_upd")
Compartments_models_noAGP_D10_gapseq <- as.data.frame(Compartments_models_noAGP_D10_gapseq)


####################################################################################
###############        GAPSEQ-ADAPTED MODELS            ############################
####################################################################################
# save list of pre-loaded models for each compartment and the table with taxonomic category for each specie
for (compartment in c( "Gizzard",  "Duodenum","Jejunum", "Ileum", "Ceca", "Colon")) {
  model_list <- c()
  models <- subset(Compartments_models_noAGP_D10_gapseq[,  c(1:3,  which(
    colnames(Compartments_models_noAGP_D10_gapseq) %in% compartment))], 
                   Compartments_models_noAGP_D10_gapseq[, compartment] == "yes") 
  
  adapted_models <- list.files("~/Desktop/Study/UofT/ParkinsonLab/gapseq/gapseq_models/2023/final_models_29taxa_6comp/")
  
  for (specie in models$Specie) {
    print(specie)

    bac_gapseq <- readSBMLmod(paste0("~/Desktop/Study/UofT/ParkinsonLab/gapseq/gapseq_models/2023/final_models_29taxa_6comp/",
                                     adapted_models[grep(specie,adapted_models)]))
    bac_gapseq@mod_name <- specie
    model_list <- c(model_list, bac_gapseq)
  }
  
  # save specie-taxcategory relation:
  saveRDS(models[, c(1,2)],  # save only Specie, TaxCategory
          paste0("/Users/irina/Desktop/Study/UofT/ParkinsonLab/Modeling_chickens/ModelGeneration_Files/Chicken_",
                 compartment, "_D10_organisms_tax-category_table_gapseq_MAGs_gapseqCobraAdapt_Aug2023.rds"))

  # save list of pre-loaded models for this compartment:
  saveRDS(model_list,
          paste0("/Users/irina/Desktop/Study/UofT/ParkinsonLab/Modeling_chickens/ModelGeneration_Files/Chicken_",
                 compartment,
                 "_noAGP_D10_gapseq_ModelList_MAGs_gapseqCobraAdapt_Aug2023.rds"))
}


##############################################################################################################
###########     generate a model list with all species from all compartments          ########################
##############################################################################################################

adapted_models <- list.files("~/Desktop/Study/UofT/ParkinsonLab/gapseq/gapseq_models/2023/final_models_29taxa_6comp/")
model_list_all <- c()

for (specie in gsub("-gapseqAdRm\\.xml|_cobrapy_adapted\\.xml","", adapted_models)) {
  bac_gapseq <- readSBMLmod(paste0("~/Desktop/Study/UofT/ParkinsonLab/gapseq/gapseq_models/2023/final_models_29taxa_6comp/",
                                   adapted_models[grep(specie,adapted_models)]))
  bac_gapseq@mod_name <- specie
  model_list_all <- c(model_list_all, bac_gapseq)
}

# save list of all pre-loaded models:
saveRDS(model_list_all,
        paste0("/Users/irina/Desktop/Study/UofT/ParkinsonLab/Modeling_chickens/ModelGeneration_Files/ChickenAllSpecies",
               "_noAGP_D10_gapseq_ModelList_MAGs_gapseqCobraAdapt_Aug2023.rds"))



