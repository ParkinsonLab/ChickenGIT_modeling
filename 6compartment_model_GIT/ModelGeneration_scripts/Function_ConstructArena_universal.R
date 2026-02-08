library(data.table)
library(magrittr)
library(BacArena)
library(glpkAPI)
library(libSBML)
library(sybilSBML)



############################################################################################
####### FUNCTION TO GET 16S RELATIVE ABUNDANCE in DAY-COMPARTMENT CATEGORY #################
############################################################################################

Get16SAbundance_table <-function(OTU_table, Compartment, Day, With_Antibiotics = F) {
  p <- regexpr(Compartment, colnames(OTU_table))
  OTU_table_Comp <- cbind(OTU_table[,c(1:8)], OTU_table[, p != -1])
  q <- regexpr(Day, colnames(OTU_table_Comp))
  OTU_table_Comp_Day <- cbind(OTU_table_Comp[,c(1:8)], OTU_table_Comp[, q != -1])
  
  OTU_table_Comp_Day <- stats::aggregate(.~taxonomy, data = OTU_table_Comp_Day[,c(1,9:ncol(OTU_table_Comp_Day))], FUN=sum)
  OTU_table_Comp_Day <- reshape2::melt(OTU_table_Comp_Day)
  
  
  if (With_Antibiotics == F) {
    metadata <- subset(metadata, metadata$Antibiotics %in% "None")
  }
  metadata$SampleID <- gsub('-','.', metadata$SampleID)
  
  for (var in OTU_table_Comp_Day$variable) {
    if (var %in% metadata$SampleID) {
      OTU_table_Comp_Day$variable <- gsub(var, paste0(metadata$Treatment[metadata$SampleID == var],
                                                      metadata$SampleNum[metadata$SampleID == var],
                                                      "-",
                                                      metadata$SampleRep[metadata$SampleID == var]),
                                          OTU_table_Comp_Day$variable)
    }
    else {
      OTU_table_Comp_Day <- subset(OTU_table_Comp_Day, OTU_table_Comp_Day$variable != var)
    }
  }
  
  OTU_table_Comp_Day <- OTU_table_Comp_Day %>% dplyr::group_by(variable)  %>% dplyr::mutate(value = value/sum(value))
  OTU_table_Comp_Day$taxonomy <- factor(OTU_table_Comp_Day$taxonomy, 
                                        levels = c("Lachnospiraceae", "Clostridiaceae","Ruminococcaceae",
                                                   "Eubacteriaceae","Peptostreptococcaceae","Other Clostridiales",
                                                   "Lactobacillaceae","Streptococcaceae","Other Bacilli",
                                                   "Erysipelotrichaceae","Other Firmicutes","Enterobacteriaceae",
                                                   "Alphaproteobacteria","Other Proteobacteria","Chloroplast",
                                                   "Other Bacteria"))
  OTU_table_Comp_Day$Diet <- stringr::str_remove(gsub('[0-9]','', OTU_table_Comp_Day$variable), '-')
  
  return(OTU_table_Comp_Day)
}



##############################################################################################################
#################### FUNCTION TO CONSTRUCT ARENA. MIXING, MUCUS, PATHOGENS, DIET #############################
##############################################################################################################

CreateModelArena = function(modelList, organisms,  sampleName, 
                            day, compartment, diet_name = c("wheat", "corn"), proportion, initial_commsize=500, econc, 
                            all_conc = 1e-05, bacteroides_perc = 0,
                            rich_medium = F, stir_var = T, mucus = F,
                            arena_size=100, speed=5, tstep=1, seed = 239, cecum_only=F) {
  
  TableRelativeAbundance <- as.data.frame(Get16SAbundance_table(OTU_table, Day = day, Compartment = compartment))
  
  if (cecum_only == F) {
    diet_concentration <- as.numeric(Parameters_table$diet_proportion[Parameters_table$Compartment == "Gizzard"])
  } else {
    # if only cecum is simulated, then proportion parameter already represents concentration of the diet
    diet_concentration <- 1
  }
  
  

  model_list <- modelList
  Chicken_organisms <- organisms

  if (diet_name == "wheat"){
    diet <- diet_wheat
  } else {
    diet <- diet_corn
  }
    
  arena <- BacArena::Arena(n=arena_size, m=arena_size, stir = stir_var, Lx=0.025*(arena_size/100), Ly=0.025*(arena_size/100), tstep = tstep,
                           seed = seed)
  
  essential_meds <- c()
  for (taxgroup in unique(Chicken_organisms$TaxCategory)) { 
    # if sample is missing - reconstruct its rel. abundance as an average of rel. abundances in samples from the same day+diet
    if (!(sampleName %in% TableRelativeAbundance$variable)) {
      abundance = mean(TableRelativeAbundance$value[TableRelativeAbundance$taxonomy %in% taxgroup & 
                                                       TableRelativeAbundance$Diet == stringr::str_remove(gsub('[0-9]','', sampleName), '-')]) /
        length(Chicken_organisms$Specie[Chicken_organisms$TaxCategory %in% taxgroup]) # divide by the number of organisms
      # selected to represent this taxgroup
    } else {
    abundance = TableRelativeAbundance$value[(TableRelativeAbundance$taxonomy %in% taxgroup) &
                                               (TableRelativeAbundance$variable %in% sampleName)] /
      length(Chicken_organisms$Specie[Chicken_organisms$TaxCategory %in% taxgroup]) # divide by the number of organisms
    # selected to represent this taxgroup
    }
    
    for (model_name in Chicken_organisms$Specie[Chicken_organisms$TaxCategory %in% taxgroup]) {
     
       if (bacteroides_perc > 0) {
        if (model_name == "Phocaeicola_vulgatus") {
          abundance = bacteroides_perc
        } else {
        # scale relative abundances of other species considering that bacteroides will take up bacteroides_perc %
        abundance = abundance * (1-bacteroides_perc)
        }
      }
     
      for (i in 1:length(model_list)) {
        # find necessary models in the ModelList and upload them through Bac function
        # and then add to arena in the amount defined above (abundance)
        if (model_list[[i]]@mod_name == model_name) {
          
          model_list[[i]] <- RenameModelReactions(model_list[[i]])
          
          print(model_name)
          if (ceiling(abundance*initial_commsize)>0){
            
            bac = BacArena::Bac(model=model_list[[i]], growtype="exponential")
            arena = BacArena::addOrg(arena, bac, amount=ceiling(abundance*initial_commsize))
            
            # save essential metabolites for each bac-organism
            # 'only_return=T' if essential metabolites should only be returned but not added to arena
            essential_meds <- unique(c(essential_meds, BacArena::addEssentialMed(arena, bac, only_return = T)))
            
          }
        }
      }
    }
  }
  
  
  # remove oxygen from essential metabolites if present!
  essential_meds <- essential_meds[which(!(essential_meds %in% c("EX_cpd00007_e0")))]
  

  if (rich_medium == T) {
    # add all metabolites
    arena <- BacArena::addSubs(arena, smax = all_conc, unit = "mM")
  }
  # add essential metabolites 
  if (length(essential_meds) > 0) {
    arena <- BacArena::addSubs(arena,smax=econc,unit="mM", mediac=essential_meds,add=T, addAnyway = T)
  }
  # add nucleobases and essential CoA:
  # gapseq: add cytosine, uracil, cytidine, adenosine, guanosine, CoA
  arena <- BacArena::addSubs(arena,smax=econc,unit="mM", 
                             mediac=c("EX_cpd00307_e0", "EX_cpd00092_e0","EX_cpd00367_e0","EX_cpd00182_e0", "EX_cpd00311_e0",
                                      "EX_cpd00010_e0"),
                             add=T, addAnyway = T)
  # add mucins and urea for cecum and colon:
  if (compartment %in% c("Ceca", "Colon")) {
    mucins_urea <- c("EX_cpd00122_e0", "EX_cpd00232_e0","EX_cpd00832_e0","EX_cpd02992_e0",
                     "EX_cpd11842_e0", "EX_cpd21520_e0","EX_cpd00073_e0")
    
    arena <- BacArena::addSubs(arena,smax=1e-03,unit="mM", 
                               mediac=mucins_urea,
                               add=T, addAnyway = T)
  }
  
  # add metabolite concentrations generated by diet design (vmh, feed tables (inrae) and foodb.ca)
  # add = F to avoid summing the amounts for metabolites present in both diet and essential (i.e., add dietary mets only in conc in the diet, not as essential)
  # if this is initial arena setup (iteration 1) in dowsntream compartments, add only essential metabolites 
  # to avoid replacing essential concentrations of overlapping metabolites by 0:
  if (proportion != 0) {
    arena <- BacArena::addSubs(arena,smax=proportion * diet_concentration * diet$met_mM,
                               unit="mM", mediac=diet$met_id, 
                               add = F, addAnyway = T)
  }
  
  if (mucus == T) {
    #add a mucus gradient
    arena <- BacArena::createGradient(arena,smax=1e-03,mediac=mucin,position='bottom',steep=0.5,add=F,unit="mM")
  }
  
  # save only essential metabolites that aren't in the diet:
   essential_meds <- essential_meds[which(!(essential_meds %in% diet$met_id))]
  
  
  return(list(arena, essential_meds))
}
