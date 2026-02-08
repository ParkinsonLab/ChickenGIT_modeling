
## Change exchange reactions constraints
# constraints table is supposed to have columns named  'Category' (with 'Inorganics', 'Other', 'Carbohydrates' and 'Amino acids' categories)
# 'met_cpd_id' - cpd_id of a reaction, i.e. EX_cpd00007_e0
# and 'uptake_lbnd' with the value of the maximum uptake flux for the reactions
# Note: the value is positive in the table, when changing lower bounds multiply by (-1) otherwise production of the metabolite will be enforced

ChangeUptakeConstraints <- function(modelList, constraints_table, Sugars_X = NULL, AA_X = NULL ) {
  for (i in 1:length(modelList)) {
    # water and inorganics
    for (compound in as.character(constraints_table$met_cpd_id[constraints_table$Category == "Other"])) {
      modelList[[i]]@lowbnd[which(modelList[[i]]@react_id == compound)] <- (-1)*constraints_table$uptake_lbnd[constraints_table$met_cpd_id == compound]
    }
   
    if (!is.null(Sugars_X)) {
      # saccharides
      sachharides <- as.character(constraints_table$met_cpd_id[constraints_table$Category == "Carbohydrates"])
      for (sachharide in sachharides) { 
          modelList[[i]]@lowbnd[which(modelList[[i]]@react_id == sachharide)] <- (-1) * as.numeric(Sugars_X) * constraints_table$uptake_lbnd[constraints_table$met_cpd_id == sachharide]
        }
    }
    
    if (!is.null(AA_X)) {
      # amino acids
      amino_acids <- as.character(constraints_table$met_cpd_id[constraints_table$Category == "Amino acids"])
      for (amino_acid in  amino_acids) { 
        modelList[[i]]@lowbnd[which(modelList[[i]]@react_id == amino_acid)] <- (-1) * as.numeric(AA_X) * constraints_table$uptake_lbnd[constraints_table$met_cpd_id == amino_acid]
      }
    }
  }
  return(modelList)
}

