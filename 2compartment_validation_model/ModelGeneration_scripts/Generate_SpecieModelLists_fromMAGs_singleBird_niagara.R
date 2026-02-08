library(BacArena)
library(libSBML)
library(sybilSBML)
library(dplyr)

myargs <- commandArgs(trailingOnly = TRUE)
print(myargs)

bird = as.character(myargs[1])

# ========================================================================
# LOAD OTU-MAG MAPPING AND SAMPLE METADATA
# ============================================================================

final_models_path = "/scratch/j/jparkin/utkinair/gapseq1.4.0/IleumCecum_2024trial/results/final_models/"
modeling_path = "/scratch/j/jparkin/utkinair/Modeling/Chicken_IleumCecum_2024Trial_simulation/ModelGeneration_Files/"

# ============================================================================
# GENERATE MODEL LISTS FOR EACH BIRD'S ILEUM-CECUM PAIR
# ============================================================================
otu_mag_metadata_paired <- read.table(paste0(modeling_path, "RepresentativeOTUs_matchedMAGs_with_metadata_paired_birds_corrected.tsv"),
                                      sep='\t',header=T)

for (compartment in c("Ileum","Cecum")) {
  
  sample_id = unique(otu_mag_metadata_paired$sample_name[
    otu_mag_metadata_paired$site == compartment &
      otu_mag_metadata_paired$bird_id == bird])
  
  model_list <- c()
  # Get MAGs present in this sample
  sample_mags <- otu_mag_metadata_paired %>%
    dplyr::filter(sample_name == sample_id, total_rel_abundance > 0) %>%
    pull(mag_id_clean)
  
  for (mag_id in sample_mags) {
    specie_id <- as.character(unique(otu_mag_metadata_paired$species_ID[otu_mag_metadata_paired$mag_id_clean == mag_id  & 
                                                                          otu_mag_metadata_paired$sample_name == sample_id]))
   
    # Load gapseq model for this MAG
    model_file <- file.path(paste0(final_models_path, mag_id, ".xml"))
    
    if (file.exists(model_file)) {
      tryCatch({
        bac_gapseq <- readSBMLmod(model_file)
        bac_gapseq@mod_name <- specie_id
        bac_gapseq@mod_desc <- specie_id
        
        model_list <- c(model_list, bac_gapseq)
      }, error = function(e) {
        warning(paste("Failed to load model for", mag_id, ":", e$message))
      })
    } else {
      warning(paste("Model file not found for MAG:", mag_id))
    }
  }
  
  # Save model list for this bird-compartment using bird_id
  saveRDS(model_list,
          paste0(modeling_path,"model_lists/", bird, "_", compartment, 
                 "_ModelList_July2025_v2.rds"))
  
  print(paste0("Saved ", length(model_list), " models for ", bird, "-", compartment))
}


