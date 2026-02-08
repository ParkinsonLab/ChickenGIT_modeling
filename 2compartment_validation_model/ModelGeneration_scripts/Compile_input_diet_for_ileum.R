library(dplyr)
library(reshape2)

setwd("/Users/irina/Desktop/Study/UofT/ParkinsonLab/Modeling_chickens/Chicken_IleumCecum_2024trial_modeling")

# load the chicken feed corn diet
diet_corn_grower <- read.table("ModelGeneration_Files/corn_diet_grower_aggregated_2023.txt", sep='\t',header=T)
# Create a copy and remove the 'value' column
diet_corn_grower_ileum <- diet_corn_grower
diet_corn_grower_ileum$value <- NULL

# remove *12 scale for cobalt introduced initially due to shadow prices - keep it at *4:
diet_corn_grower_ileum$met_mM[diet_corn_grower_ileum$met_name == "Cobalt"] <- 
  diet_corn_grower_ileum$met_mM[diet_corn_grower_ileum$met_name == "Cobalt"]/3

##########################################################################################################################
##########     Adjust concentrations to imitate absorption/breakdown of metabolites before reaching the ileum   ##########
##########################################################################################################################

### 0) Minerals 
## Phosphorus: ileal digestibility 67% for corn diet (10.3382/ps.2013-03419); 
# ileal digestibility ~55-60% for grower diets with soybean meal (10.1016/j.psj.2024.103602)
# ileal digestibility 48-70% for different corn-based diets with supplementation (doi: 10.5187/jast.2019.61.4.192)
# total tract retention for corn-based diet is ~0.5 (10.5713/ajas.2010.90129)
## Calcium: ~55-65% ileal digestibility for grower diets with soybean meal (10.1016/j.psj.2024.103602);
# ileal digestibility coefficients of limestone, MCP, and DCP - 0.51, 0.43, and 0.32, respectively (10.3382/ps/pez314)
# ileal digestibility 50-75% for different corn-based diets with supplementation (doi: 10.5187/jast.2019.61.4.192)
## Magnesium: ~30-35% ileal digestibility for different corn-based diets with supplementation (doi: 10.5187/jast.2019.61.4.192)
# ~50% dietary Mg absorbed until ileum (https://doi.org/10.1093/jn/103.6.875)
## Potassium: ~80-90% ileal digestibility for different corn-based diets with supplementation (doi: 10.5187/jast.2019.61.4.192)
# total tract retention in corn-based fed chickens - 30-35% (10.5713/ajas.2010.90129)
## Manganese: ~25-40% ileal digestibility for different corn-based diets with supplementation (doi: 10.5187/jast.2019.61.4.192)
## Zinc: ~30-40% ileal digestibility for different corn-based diets with supplementation (doi: 10.5187/jast.2019.61.4.192)
# ileum is main site of zinc absorption
## Iron: ~20-50% ileal digestibility for different corn-based diets with supplementation (doi: 10.5187/jast.2019.61.4.192)
## Copper: Copper absorption occurs mainly in the duodenum of chickens (doi: 10.3382/ps.0740502) - keep ~ 10%?
## Sodium: little info, but total tract retention for corn-based fed chickens on d14 is 0.68 (10.5713/ajas.2010.90129)

# Adjust percentages of intake amounts for specific minerals (retention %), 
# considering they should provide an average of digestibility measurements between distal jejunum and proximal ileum
# Ca - 50% (A conservative mean of 0.5 (ingredient-level));
# Phosphorus - 45%; Potassium - 25%; Sodium - 60%; Magnesium - 60%; Iron - 70%; Manganese - 75%; Zinc - 70%, Copper - 10%
mineral_names <- c("Calcium", "Phosphorus", "Potassium", "Sodium", "Magnesium", "Fe3+", "Manganese", "Zinc", "Copper")
mineral_factors <- c(0.5, 0.45, 0.25, 0.6, 0.6, 0.7, 0.75, 0.7, 0.1)
names(mineral_factors) <- mineral_names

diet_corn_grower_ileum <- diet_corn_grower_ileum |>
  mutate(met_mM = ifelse(met_name %in% names(mineral_factors),
                         met_mM * mineral_factors[met_name], met_mM))


###  1) Carbohydrates 
# starches: 10.1016/j.aninu.2022.01.003 (~10-15% is retained by distal jejunum)
# ~80% of starch is figested by distal jejunum and ~90% by proximal ileum 
# during grower corn-based diet phase (doi:10.1017/S0007114517003257, doi:10.3382/ps/pev244)

# retain 15% of starches
starch_ids <- c("EX_cpd90003_e0", "EX_cpd90004_e0")
diet_corn_grower_ileum$met_mM[diet_corn_grower_ileum$met_id %in% starch_ids] <- 
  0.15 * diet_corn_grower_ileum$met_mM[diet_corn_grower_ileum$met_id %in% starch_ids]


# ~5% glucose is retained in distal jejunum (https://doi.org/10.1016/j.psj.2020.08.072)
# keep 5% of glucose and fructose (mono-), and 10% of sucrose, maltose, stachyose, raffinose (di-, trisaccharides)
diet_corn_grower_ileum$met_mM <- as.numeric(diet_corn_grower_ileum$met_mM)
diet_corn_grower_ileum$met_mM[diet_corn_grower_ileum$met_name %in% c("D-Glucose", "D-Fructose")] <- 
  0.05 * diet_corn_grower_ileum$met_mM[diet_corn_grower_ileum$met_name %in% c("D-Glucose", "D-Fructose")]
diet_corn_grower_ileum$met_mM[diet_corn_grower_ileum$met_name %in%  c("Galactose", "Melitose","Sucrose", "Maltose","Stachyose")] <- 
  0.1 * diet_corn_grower$met_mM[diet_corn_grower$met_name %in% c("Galactose", "Melitose","Sucrose", "Maltose","Stachyose")]


###  2) Proteins (amino acids) 
# ~70-80% is absorbed in distal jejunum (10.1016/j.aninu.2022.01.003)
# ~66% of protein digested in distal jejunum during grower corn-based diet,
# ~75% - in proximal ileum during grower corn-based diet (doi:10.1017/S0007114517003257)
# retain 25% of AAs 
diet_corn_grower_ileum$met_mM[diet_corn_grower_ileum$category == "Proteins"] <-
  0.25 * diet_corn_grower$met_mM[diet_corn_grower$category == "Proteins"]

###  3) Fats
## lower jejunum: total fat ~65% digested
# 71% palmitic, 44% stearic, 80% of oleic and 90% of linoleic acids digested (https://doi.org/10.3382/ps.2013-03344) 
## upper ileum: total fat ~80% digested
# 83% palmitic, 72% stearic, 89% of oleic and 94% of linoleic acids digested (https://doi.org/10.3382/ps.2013-03344) 

# Retain 25% of palmitic, 40% of stearic, 15% of oleic and 8% of linoleic acids, and 30% of remaining fats:
diet_corn_grower_ileum$met_mM[diet_corn_grower_ileum$met_name %in% c("Palmitic acid")] <- 
  0.25 * diet_corn_grower_ileum$met_mM[diet_corn_grower_ileum$met_name %in% c("Palmitic acid")]
diet_corn_grower_ileum$met_mM[diet_corn_grower_ileum$met_name %in% c("Stearic acid")] <- 
  0.4 * diet_corn_grower_ileum$met_mM[diet_corn_grower_ileum$met_name %in% c("Stearic acid")]
diet_corn_grower_ileum$met_mM[diet_corn_grower_ileum$met_name %in% c("Oleic acid")] <- 
  0.15 * diet_corn_grower_ileum$met_mM[diet_corn_grower_ileum$met_name %in% c("Oleic acid")]
diet_corn_grower_ileum$met_mM[diet_corn_grower_ileum$met_name %in% c("Linoleic acid")] <- 
  0.08 * diet_corn_grower_ileum$met_mM[diet_corn_grower_ileum$met_name %in% c("Linoleic acid")]


diet_corn_grower_ileum$met_mM[
  diet_corn_grower_ileum$category == "Fats" & 
    !(diet_corn_grower_ileum$met_name %in% c("Palmitic acid","Stearic acid","Oleic acid","Linoleic acid"))] <-
  0.3 * diet_corn_grower$met_mM[
    diet_corn_grower$category == "Fats" & 
      !(diet_corn_grower_ileum$met_name %in% c("Palmitic acid","Stearic acid","Oleic acid","Linoleic acid"))]

# post-simulation addition:
# from shadow prices add dodecanoic acid (ddca) and capric acid (dca) 
# at roughly similar concentration as other fatty acids (on the lower side):
diet_corn_grower_ileum <- rbind(diet_corn_grower_ileum,
                                c("Dodecanoic acid", "g_per100g","Dodecanoic acid", "EX_cpd01741_e0",
                                  "0.01", "Fats"),
                                c("Capric acid", "g_per100g","Decanoate", "EX_cpd01107_e0",
                                  "0.01", "Fats"))

### 4) Fibre cellulose 
# retain 50% of cellulose:
diet_corn_grower_ileum$met_mM[diet_corn_grower_ileum$met_name == "cellulose"] <- 
  0.5 * diet_corn_grower$met_mM[diet_corn_grower$met_name == "cellulose"]


# Based on shadow prices, increase concentration of Riboflavin by 3x and Pantothenic acid by 20x
diet_corn_grower_ileum$met_mM[diet_corn_grower_ileum$met_name %in% c("Riboflavin")] <- 
   3 * as.numeric(diet_corn_grower_ileum$met_mM[diet_corn_grower_ileum$met_name %in% c("Riboflavin")])
diet_corn_grower_ileum$met_mM[diet_corn_grower_ileum$met_name %in% c("Pantothenic acid")] <- 
   20 * as.numeric(diet_corn_grower_ileum$met_mM[diet_corn_grower_ileum$met_name %in% c("Pantothenic acid")])

# Save tables for input:
write.table(diet_corn_grower_ileum, "ModelGeneration_Files/corn_diet_grower_ileum.txt",
            sep = '\t',quote=F)


