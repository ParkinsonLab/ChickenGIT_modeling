###############################################################################
# Diet_feedtables_design_2023.R
#
# Compile metabolic representations of corn-based and wheat-based
# chicken diets (starter and grower phases) from ingredient-level
# nutritional data (FeedTables), map each nutrient to its ModelSEED
# metabolite identifier, and convert concentrations to mM
#
# INPUT FILES 
#   - Diet_FeedTables_ingredients.xlsx   (sheets: FeedTable, Sugars_content,
#                                         DietsComposition_starter,
#                                         DietsComposition_grower)
#   - all_seed_metabolites_edited.tsv    (ModelSEED metabolite database)

###############################################################################

library(ggplot2)
library(reshape2)
library(dplyr)
library(stringr)

base_dir <- "Input_diets_chicken/"
setwd(base_dir)


###############################################################################
#  LOAD RAW DATA
###############################################################################

# Feed ingredient nutritional composition (from published feed tables)
FeedTable_diet <- as.data.frame(
  readxl::read_xlsx("Diet_FeedTables_ingredients.xlsx", sheet = "FeedTable"))

# Diet formulations (% of each ingredient in the final feed)
Diet_components_starter <- as.data.frame(
  readxl::read_xlsx("Diet_FeedTables_ingredients.xlsx", sheet = "DietsComposition_starter"))
Diet_components_grower <- as.data.frame(
  readxl::read_xlsx("Diet_FeedTables_ingredients.xlsx", sheet = "DietsComposition_grower"))

# ModelSEED metabolite reference database (for mapping nutrient names to cpd IDs)
seed_metabolites <- readr::read_delim("all_seed_metabolites_edited.tsv",
                                      delim = '\t', col_names = TRUE)


###############################################################################
#  RESHAPE FEED TABLE: wide -> long format, separate units from metabolite names
###############################################################################

FeedTable_diet <- melt(FeedTable_diet, id.vars = c("DietTable_name", "FeedTable_name"))

# Extract the unit (last word in the column name, e.g., "g/kg", "mg/kg")
FeedTable_diet$unit <- sapply(as.character(FeedTable_diet$variable), function(x) {
  parts <- strsplit(x, " ")[[1]]
  parts[length(parts)]
})

# Extract the metabolite name (everything except the trailing unit)
FeedTable_diet$met_name <- mapply(function(var, unit) {
  gsub(paste0(" ", unit), "", var)
}, FeedTable_diet$variable, FeedTable_diet$unit)


###############################################################################
#  MAP GENERIC NUTRIENTS TO SPECIFIC METABOLITES
###############################################################################

# --- Fibre -> Cellulose ---
# "Fibre" in feed tables refers predominantly to cellulose (the primary
# structural polysaccharide and major fibre component in cereal grains).
# We represent it as cellulose for metabolic modeling.
FeedTable_diet$met_name[FeedTable_diet$met_name == "Fibre"] <- "cellulose"

# --- Starch -> two structural variants (equal split) ---
# Feed starch is a mixture of amylose and amylopectin with varying chain lengths.
# gapseq represents starch as two chain-length variants:
#   - starch (n=19, 3xalpha1-6, 15xalpha1-4): shorter chains
#   - starch (n=27, 3xalpha1-6, 23xalpha1-4): longer chains
# In the absence of measured chain-length distributions for these feed ingredients,
# we split the total starch equally between the two variants (50:50).
starch19_subset <- FeedTable_diet[FeedTable_diet$met_name == "Starch", ]
starch19_subset$met_name <- "starch (n=19, 3xalpha1-6, 15xalpha1-4)"
starch19_subset$value <- starch19_subset$value / 2

FeedTable_diet$value[FeedTable_diet$met_name == "Starch"] <-
  FeedTable_diet$value[FeedTable_diet$met_name == "Starch"] / 2
FeedTable_diet$met_name[FeedTable_diet$met_name == "Starch"] <-
  "starch (n=27, 3xalpha1-6, 23xalpha1-4)"

FeedTable_diet <- rbind(FeedTable_diet, starch19_subset)

# --- Fe -> Fe2+ and Fe3+ (equal split) ---
# Dietary iron exists in both ferrous (Fe2+) and ferric (Fe3+) forms.
# Since the exact Fe2+/Fe3+ ratio depends on feed processing and is not
# typically reported in standard feed tables, we assume an equal 50:50 split.
fe3_subset <- FeedTable_diet[FeedTable_diet$met_name == "Fe", ]
fe3_subset$met_name <- "Fe3+"
fe3_subset$value <- fe3_subset$value / 2

FeedTable_diet$value[FeedTable_diet$met_name == "Fe"] <-
  FeedTable_diet$value[FeedTable_diet$met_name == "Fe"] / 2
FeedTable_diet$met_name[FeedTable_diet$met_name == "Fe"] <- "Fe2+"

FeedTable_diet <- rbind(FeedTable_diet, fe3_subset)


###############################################################################
#  ASSIGN ModelSEED METABOLITE IDs
###############################################################################

# ---  matching by metabolite name ---
FeedTable_diet$met_id <- NA
for (met in unique(FeedTable_diet$met_name)) {
  if (met %in% seed_metabolites$name) {
    FeedTable_diet$met_id[FeedTable_diet$met_name == met] <-
      seed_metabolites$id[seed_metabolites$name == met]
  }
}

# --- Correct IDs for metabolites where the default match is not the
#     biologically relevant form used in bacterial models ---
# Fe2+: cpd00021 (generic Fe) -> cpd10515 (Fe2+ ion, used in bacterial models)
FeedTable_diet$met_id[FeedTable_diet$met_id == "cpd00021"] <- "cpd10515"
# Calcium: cpd20429 (Ca, generic) -> cpd00149 (Ca2+, used in bacterial models)
FeedTable_diet$met_id[FeedTable_diet$met_id == "cpd20429"] <- "cpd00149"
# Cobalt: cpd29674 (generic) -> cpd00063 (Co2+, required cofactor form)
FeedTable_diet$met_id[FeedTable_diet$met_id == "cpd29674"] <- "cpd00063"

# --- Manual ID assignments for nutrients not found by name matching ---
# (different naming conventions between feed tables and ModelSEED)
manual_met_ids <- c(
  "Palmitic acid"    = "cpd00214",
  "Stearic acid"     = "cpd01080",
  "Oleic acid"       = "cpd15269",
  "Linoleic acid"    = "cpd01122",
  "Linolenic acid"   = "cpd03850",
  "Potassium"        = "cpd00205",
  "Chloride"         = "cpd00099",
  "Linoleic acid"    = "cpd01122",   
  "Palmitoleic acid" = "cpd15237",
  "Magnesium"        = "cpd00254",
  "Phosphorus"       = "cpd00009",
  "Zinc"             = "cpd00034",
  "Iodine"           = "cpd00994",
  "Manganese"        = "cpd00030",
  "Copper"           = "cpd00058"
)

for (met_name in names(manual_met_ids)) {
  FeedTable_diet$met_id[FeedTable_diet$met_name == met_name] <- manual_met_ids[met_name]
}

# --- D-Methionine -> L-Methionine ---
# D-Methionine is not utilized by the bacterial models in our set. In poultry
# nutrition, supplemental DL-Methionine is converted to L-Methionine in vivo.
# We reassign D-Methionine amounts to L-Methionine for modeling purposes.
FeedTable_diet$met_id[FeedTable_diet$met_name == "D-Methionine"] <- "cpd00060"
FeedTable_diet$variable[FeedTable_diet$variable == "D-Methionine g/kg"] <- "L-Methionine g/kg"
FeedTable_diet$met_name[FeedTable_diet$met_name == "D-Methionine"] <- "L-Methionine"


###############################################################################
#   ADD INDIVIDUAL SUGAR COMPOSITIONS
###############################################################################
# Standard feed tables report only total sugars for each ingredient.
# We replace total sugars with individual sugar species based on published
# compositional analyses of each ingredient:
#
# --- CORN (total sugars: 17 g/kg) ---
# Sugar proportions from:
#   - Mabelebele et al. (2018) S Afr J Plant Soil 35:237-243
#     (fructose:glucose:sucrose ~1.5:1.5:1)
#   - Nieto-Calvache et al. (2011) Ciencia y Tecnologia Alimentaria
#     (fructose:glucose:sucrose ~1:1:1.5)
# We use an equal 1:1:1 split (5.67 g/kg each) as a simplification,
# given the variability across sources and corn varieties.
#
# --- WHEAT (total sugars: 26 g/kg) ---
# Sugar proportions from Ziegler et al. (2021) Eur Food Res Technol 247:213
# (Table 6: whole meal wheat flour), rescaled to 26 g/kg total:
#   Glucose 9.36, Fructose 2.8, Sucrose 2.16, Raffinose 2.6, Maltose 9.07 g/kg
#
# --- SOYBEAN MEAL (total sugars: 87 g/kg, rescaled from 137 g/kg) ---
# Individual sugars from Giannoccaro et al. (2006) J Agric Food Chem 51:7684-7691,
# rescaled proportionally to the 87 g/kg total in our feed table:
#   Galactose 2.35, Glucose 1.97, Fructose 1.91, Sucrose 30.5,
#   Raffinose 3.94, Stachyose 24.4, Verbascose 1.02, Uronic acid 21 g/kg
#
# --- CORN DISTILLERS GRAINS AND SOLUBLES (DDGS) (total sugars: 8 g/kg) ---
# Sugar proportions from Belyea et al. (1989) J Agric Food Chem 37:1024,
# ratio glycerol:arabinose:xylose:mannose:glucose:galactose ~ 8:6:8:2:12:2:
#   Glycerol 0.22, Arabinose 0.17, Xylose 0.22, Mannose 0.055,
#   Galactose 0.055, Glucose 0.33 g/kg
#
# --- CANOLA MEAL (total sugars: 95 g/kg) ---
# Sugar proportions from Maison (2013) PhD thesis, U Illinois,
# rescaled to 95 g/kg total:
#   Fructose 2.17, Glucose 0.41, Sucrose 71.1, Raffinose 5.28,
#   Stachyose 16.04 g/kg

Sugars_table <- as.data.frame(
  readxl::read_xlsx("Diet_FeedTables_ingredients.xlsx", sheet = "Sugars_content"))
Sugars_table <- Sugars_table[, -3]  # remove total sugars column
Sugars_table <- melt(Sugars_table, id.vars = c("DietTable_name", "FeedTable_name"))

Sugars_table$unit <- sapply(as.character(Sugars_table$variable), function(x) {
  parts <- strsplit(x, " ")[[1]]
  parts[length(parts)]
})

Sugars_table$met_name <- mapply(function(var, unit) {
  gsub(paste0(" ", unit), "", var)
}, Sugars_table$variable, Sugars_table$unit)

# Glucuronic acid -> Glucuronate (name used in ModelSEED)
Sugars_table$met_name[Sugars_table$met_name == "beta-D-Glucuronic acid"] <- "Glucuronate"

# Assign ModelSEED IDs to sugars
Sugars_table$met_id <- NA
for (met in unique(Sugars_table$met_name)) {
  if (met %in% seed_metabolites$name) {
    # If multiple IDs match, take the first (lowest cpd number) as it is
    # most likely the non-deprecated entry present in bacterial models
    ids <- seed_metabolites$id[seed_metabolites$name == met]
    Sugars_table$met_id[Sugars_table$met_name == met] <- ids[1]
  }
}

# Replace aggregated "Sugars" row with individual sugar species
FeedTable_diet <- FeedTable_diet[FeedTable_diet$met_name != "Sugars", ]
FeedTable_diet_complete <- rbind(FeedTable_diet, Sugars_table)


###############################################################################
#  HELPER FUNCTIONS: compile diet and convert units
###############################################################################

#' Compile a diet from ingredient percentages and the nutrient composition table
#' @param diet_components  data.frame with columns: Ingredient, <diet_column>
#' @param diet_column      name of the column containing ingredient percentages
#' @param nutrient_table   the complete FeedTable with per-ingredient nutrient data
#' @return data.frame with per-component nutrient contributions
compile_diet_nutrients <- function(diet_components, diet_column, nutrient_table) {
  diet_nutrients <- NULL
  for (component in diet_components$Ingredient) {
    pct <- diet_components[[diet_column]][diet_components$Ingredient == component] / 100
    if (pct != 0) {
      component_nutrients <- nutrient_table[nutrient_table$DietTable_name == component, ]
      component_nutrients$value <- pct * component_nutrients$value
      diet_nutrients <- rbind(diet_nutrients, component_nutrients)
    }
  }
  return(diet_nutrients)
}

#' Aggregate nutrient contributions across ingredients, convert to g/100g and mM
#' @param diet_nutrients     output of compile_diet_nutrients()
#' @param seed_metabolites   ModelSEED metabolite table (with molar mass)
#' @return aggregated data.frame with met_mM column
aggregate_and_convert <- function(diet_nutrients, seed_metabolites) {
  # Sum nutrient contributions across all ingredients
  diet_agg <- aggregate(value ~ variable + unit + met_name + met_id,
                        diet_nutrients, sum)

  # Convert all units to g/100g of feed
  diet_agg$value[diet_agg$unit == "mg/kg"]     <- diet_agg$value[diet_agg$unit == "mg/kg"] / 1000
  diet_agg$value[diet_agg$unit == "microg/kg"] <- diet_agg$value[diet_agg$unit == "microg/kg"] / 1e6
  diet_agg$value <- diet_agg$value / 10  # g/kg -> g/100g
  diet_agg$unit  <- "g_per100g"

  # Remove metabolites with zero concentration
  diet_agg <- diet_agg[diet_agg$value != 0, ]

  # Convert g/100g to millimolar (mM) using molar mass from ModelSEED
  diet_agg$met_mM <- NA
  for (id in diet_agg$met_id) {
    molar_mass <- as.numeric(seed_metabolites$mass[seed_metabolites$id == id])
    if (length(molar_mass) > 0 && !is.na(molar_mass)) {
      diet_agg$met_mM[diet_agg$met_id == id] <-
        diet_agg$value[diet_agg$met_id == id] * 1000 / molar_mass
    }
  }

  # Format metabolite IDs as exchange reaction IDs (BacArena convention)
  diet_agg$met_id <- paste0("EX_", diet_agg$met_id, "_e0")

  return(diet_agg)
}

#' Assign nutrient categories (used by absorption functions in the simulation)
#' @param diet_agg  aggregated diet table, sorted by variable name
#' @param carb_idx  row indices for carbohydrates
#' @param prot_idx  row indices for proteins
#' @param fat_idx   row indices for fats
#' @return diet_agg with category column
assign_categories <- function(diet_agg, carb_idx, prot_idx, fat_idx) {
  diet_agg <- diet_agg[order(diet_agg$variable), ]
  diet_agg$category <- "Vitamins and minerals"  # default
  diet_agg$category[carb_idx] <- "Carbohydrates"
  diet_agg$category[prot_idx] <- "Proteins"
  diet_agg$category[fat_idx]  <- "Fats"
  return(diet_agg)
}

#' Increase cobalt concentration by a scaling factor.
#' Cobalt is a required cofactor (vitamin B12 precursor) and the primary
#' growth-limiting "shadow metabolite" for most gapseq-reconstructed models.
#' Feed-table cobalt levels are below the threshold needed for in silico
#' growth of many species. We scale cobalt by 12x to ensure it does not
#' artifactually limit growth in simulations.
#' @param diet_agg        aggregated diet table
#' @param scale_factor    multiplicative factor (default: 12)
adjust_cobalt <- function(diet_agg, scale_factor = 12) {
  cobalt_rows <- diet_agg$met_name == "Cobalt"
  diet_agg$met_mM[cobalt_rows]  <- scale_factor * diet_agg$met_mM[cobalt_rows]
  diet_agg$value[cobalt_rows]   <- scale_factor * diet_agg$value[cobalt_rows]
  return(diet_agg)
}


###############################################################################
#   COMPILE AND EXPORT DIETS
###############################################################################

# ---- STARTER DIETS ----
Diet_components <- Diet_components_starter

# Corn-based starter
corn_starter_nutrients <- compile_diet_nutrients(Diet_components, "corn_diet", FeedTable_diet_complete)
write.table(corn_starter_nutrients, "corn_diet_starter_allcomponents_allnutrients_2023.txt",
            sep = '\t', quote = FALSE)

corn_starter_agg <- aggregate_and_convert(corn_starter_nutrients, seed_metabolites)
corn_starter_agg <- assign_categories(corn_starter_agg,
                                      carb_idx = c(1:3, 58:69),
                                      prot_idx = 18:36,
                                      fat_idx  = 37:45)
corn_starter_agg <- adjust_cobalt(corn_starter_agg)
write.table(corn_starter_agg, "corn_diet_starter_aggregated_2023.txt",
            sep = '\t', quote = FALSE)

# Wheat-based starter
wheat_starter_nutrients <- compile_diet_nutrients(Diet_components, "wheat_diet", FeedTable_diet_complete)
write.table(wheat_starter_nutrients, "wheat_diet_starter_allcomponents_allnutrients_2023.txt",
            sep = '\t', quote = FALSE)

wheat_starter_agg <- aggregate_and_convert(wheat_starter_nutrients, seed_metabolites)
wheat_starter_agg <- assign_categories(wheat_starter_agg,
                                       carb_idx = c(1:3, 56:64),
                                       prot_idx = 18:36,
                                       fat_idx  = 37:45)
wheat_starter_agg <- adjust_cobalt(wheat_starter_agg)
write.table(wheat_starter_agg, "wheat_diet_starter_aggregated_2023.txt",
            sep = '\t', quote = FALSE)


# ---- GROWER DIETS ----
Diet_components <- Diet_components_grower

# Corn-based grower
corn_grower_nutrients <- compile_diet_nutrients(Diet_components, "corn_diet", FeedTable_diet_complete)
write.table(corn_grower_nutrients, "corn_diet_grower_allcomponents_allnutrients_2023.txt",
            sep = '\t', quote = FALSE)

corn_grower_agg <- aggregate_and_convert(corn_grower_nutrients, seed_metabolites)
corn_grower_agg <- assign_categories(corn_grower_agg,
                                     carb_idx = c(1:3, 58:69),
                                     prot_idx = 18:36,
                                     fat_idx  = 37:45)
corn_grower_agg <- adjust_cobalt(corn_grower_agg)
write.table(corn_grower_agg, "corn_diet_grower_aggregated_2023.txt",
            sep = '\t', quote = FALSE)

# Wheat-based grower
wheat_grower_nutrients <- compile_diet_nutrients(Diet_components, "wheat_diet", FeedTable_diet_complete)
write.table(wheat_grower_nutrients, "wheat_diet_grower_allcomponents_allnutrients_2023.txt",
            sep = '\t', quote = FALSE)

wheat_grower_agg <- aggregate_and_convert(wheat_grower_nutrients, seed_metabolites)
wheat_grower_agg <- assign_categories(wheat_grower_agg,
                                      carb_idx = c(1:3, 56:64),
                                      prot_idx = 18:36,
                                      fat_idx  = 37:45)
wheat_grower_agg <- adjust_cobalt(wheat_grower_agg)
write.table(wheat_grower_agg, "wheat_diet_grower_aggregated_2023.txt",
            sep = '\t', quote = FALSE)
