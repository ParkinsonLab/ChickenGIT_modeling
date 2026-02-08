library(BacArena)
library(glpkAPI)
library(sybilSBML)


RenameModelReactions = function(model) {
  model@react_id <- gsub('__', '_',model@react_id) 
  model@met_id <- gsub('__', '_',model@met_id) 
  model@met_id <- gsub('\\[', '(',model@met_id) 
  model@met_id <- gsub('\\]', ')',model@met_id) 
  return(model)
}