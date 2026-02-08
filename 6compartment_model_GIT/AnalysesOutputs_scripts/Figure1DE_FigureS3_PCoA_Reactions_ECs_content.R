# Figure 1D/E and Figure S3: PCoA of reaction and EC content across samples/compartments.
# Loads 16S data, builds abundance-weighted reaction/EC matrices, performs PCoA,
# analyzes reaction loadings, runs PERMANOVA and dispersion tests, and generates
# density plots for reaction/EC counts per compartment.

library(data.table)
library(magrittr)
library(BacArena)
library(dplyr)

Get16SAbundance_table <- function(OTU_table, Compartment, Day, With_Antibiotics = F) {
  p <- regexpr(Compartment, colnames(OTU_table))
  OTU_table_Comp <- cbind(OTU_table[, c(1:8)], OTU_table[, p != -1])
  q <- regexpr(Day, colnames(OTU_table_Comp))
  OTU_table_Comp_Day <- cbind(OTU_table_Comp[, c(1:8)], OTU_table_Comp[, q != -1])

  OTU_table_Comp_Day <- stats::aggregate(.~taxonomy, data = OTU_table_Comp_Day[, c(1, 9:ncol(OTU_table_Comp_Day))], FUN = sum)
  OTU_table_Comp_Day <- reshape2::melt(OTU_table_Comp_Day)

  if (With_Antibiotics == F) {
    metadata <- subset(metadata, metadata$Antibiotics %in% "None")
  }
  metadata$SampleID <- gsub('-', '.', metadata$SampleID)

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

  OTU_table_Comp_Day <- OTU_table_Comp_Day %>% dplyr::group_by(variable) %>% dplyr::mutate(value = value / sum(value))
  OTU_table_Comp_Day$taxonomy <- factor(OTU_table_Comp_Day$taxonomy,
                                        levels = c("Lachnospiraceae", "Clostridiaceae", "Ruminococcaceae",
                                                   "Eubacteriaceae", "Peptostreptococcaceae", "Other Clostridiales",
                                                   "Lactobacillaceae", "Streptococcaceae", "Other Bacilli",
                                                   "Erysipelotrichaceae", "Other Firmicutes", "Enterobacteriaceae",
                                                   "Alphaproteobacteria", "Other Proteobacteria", "Chloroplast",
                                                   "Other Bacteria"))
  OTU_table_Comp_Day$Diet <- stringr::str_remove(gsub('[0-9]', '', OTU_table_Comp_Day$variable), '-')

  return(OTU_table_Comp_Day)
}


OTU_table <- read.csv("OTUtable.csv")
metadata <- read.table("ModelGeneration_Files/metadata.txt", header = T)
Compartments_models_noAGP_D10_gapseq <- readxl::read_xlsx("/Users/irina/Desktop/Study/UofT/ParkinsonLab/Modeling_chickens/ModelGeneration_Files/SpecieModels_present_in_compartments_D10_noAGP_2022upd.xlsx",
                                                          sheet = "Gapseq_ModelList_upd")
Compartments_models_noAGP_D10_gapseq <- as.data.frame(Compartments_models_noAGP_D10_gapseq)

Compartments_list <- list("Gizzard", "Duodenum", "Jejunum", "Ileum", "Ceca", "Colon")

# Load names of organisms in each compartment community model and corresponding tax category
file_path = ""
organism_files <- list.files(file.path(file_path, "ModelGeneration_Files"), pattern = "Chicken_.*_organisms_tax-category_table_gapseq_MAGs_gapseqCobraAdapt_Aug2023.rds", full.names = TRUE)
# rearrange in the GIT order
organism_files <- organism_files[c(4, 3, 6, 5, 1, 2)]
organism_lists <- lapply(organism_files, readRDS)
names(organism_lists) = unlist(Compartments_list)

# get family-level abundances in each sample-compartment and add species from each compartment:
abund_table_10samples_allcomp <- c()
for (comp in Compartments_list) {
  abund_table_comp <- Get16SAbundance_table(OTU_table, Compartment = comp, Day = "D10")
  colnames(abund_table_comp)[1:3] <- c("TaxCategory", "Sample", "RelAbundance")
  abund_table_comp$Compartment = comp
  abund_table_comp_species <- merge(abund_table_comp, organism_lists[[comp]], by = "TaxCategory")
  abund_table_10samples_allcomp <- rbind(abund_table_10samples_allcomp,
                                         abund_table_comp_species)
}

# divide relative abundance (currently for family level) by the number of species within a family -> specie-level abundance
abund_table_10samples_allcomp <-
  abund_table_10samples_allcomp %>%
  group_by(Compartment, Sample, TaxCategory, Diet) %>%
  mutate(RelAbundance_specie = RelAbundance / n())

write.table(abund_table_10samples_allcomp,
            paste0(file_path, "/ModelsAnalyses_outputs/RelativeAbundanceTable_allSamples_bothDiets_D10_noAGP.txt"),
            sep = '\t', quote = F, row.names = F)

### now construct sample_compartment_specie_reacs and sample_compartment_reacs lists
# load all sbml models
model_list_all_compartments <-
  readRDS(paste0(file_path, "/ModelGeneration_Files/ChickenAllSpecies_noAGP_D10_gapseq_ModelList_MAGs_gapseqCobraAdapt_Aug2023.rds"))

# generate abundance table with model descriptions
specie_abundance_table <- abund_table_10samples_allcomp
specie_abundance_table$Sample_Comp <- paste0(specie_abundance_table$Sample, "_",
                                             specie_abundance_table$Compartment)

# make a list of reactions in each model:
models_reacs_list <- list()
for (model in model_list_all_compartments) {
  model_name <- model@mod_name
  models_reacs_list[[model_name]] <- model@react_id
  # remove biomass reacs:
  models_reacs_list[[model_name]] <- models_reacs_list[[model_name]][-grep("bio1|EX_cpd11416_c0",
                                                                           models_reacs_list[[model_name]])]
}

sample_compartment_models_reactions_list <- list()
sample_compartment_reactions_list <- list()
for (sample in unique(abund_table_10samples_allcomp$Sample)) {
  for (compartment in unique(abund_table_10samples_allcomp$Compartment)) {
    sample_compartment_models_reactions_list[[paste0(sample, "_", compartment)]] <- list()
    sample_compartment_reactions_list[[paste0(sample, "_", compartment)]] <- list()
    for (specie in unique(abund_table_10samples_allcomp$Specie[
      abund_table_10samples_allcomp$Sample == sample & abund_table_10samples_allcomp$Compartment == compartment
    ])) {
      sample_compartment_models_reactions_list[[paste0(sample, "_", compartment)]][[specie]] <- models_reacs_list[[specie]]
    }
    # make sample_compartment_reactions_list by unlisting and dereplicating the reactions:
    sample_compartment_reactions_list[[paste0(sample, "_", compartment)]] <-
      unique(unlist(sample_compartment_models_reactions_list[[paste0(sample, "_", compartment)]]))
  }
}

saveRDS(models_reacs_list, "ModelsAnalyses_outputs/AllModels_reacs_list.RDS")
saveRDS(sample_compartment_reactions_list, "ModelsAnalyses_outputs/SampleCompartment_reacs_list.RDS")

######### construct similar list but for ECs:
ECs_models <- read.delim("ModelsAnalyses_outputs/Table_ECnumbers_inEachModel_gapseqAdapted_cobrapyAdapted_July2023.txt", sep = '\t', header = T)

ECs_sample_compartment_models <- list()
for (sample in sort(unique(abund_table_10samples_allcomp$Sample))) {
  for (compartment in unique(abund_table_10samples_allcomp$Compartment)) {
    ECs_sample_compartment_models[[paste0(sample, "_", compartment)]] <- list()
    models_sample_comp <- abund_table_10samples_allcomp$Specie[abund_table_10samples_allcomp$Sample == sample &
                                                                 abund_table_10samples_allcomp$Compartment == compartment]
    for (specie in models_sample_comp) {
      ECs <- unlist(strsplit(ECs_models$EC_numbers[ECs_models$SpecieName == specie], ", ")[[1]])
      ECs_sample_compartment_models[[paste0(sample, "_", compartment)]][[specie]] <- ECs
    }
  }
}


############################################################################################
###########      PCOA:  REACTIONS IN SAMPLES_COMPARTMENTS             ######################
########      SUMMED (BY REL.ABUNDANCES OF SPECIES) -  BRAY           ######################
############################################################################################
library(vegan)

get_weightedByRelAbund_reacs_list <- function(sample_compartment_models_reactions_list) {
  weighted_reactions <- list()
  for (sample in names(sample_compartment_models_reactions_list)) {
    species_reactions <- sample_compartment_models_reactions_list[[sample]]
    weighted_reactions_sample <- c()

    for (specie in names(species_reactions)) {
      rel_abundance <- specie_abundance_table %>%
        filter(Specie == specie, Sample_Comp == sample) %>%
        pull(RelAbundance_specie)

      reactions <- species_reactions[[specie]]

      for (reaction in reactions) {
        if (reaction %in% names(weighted_reactions_sample)) {
          weighted_reactions_sample[reaction] <- weighted_reactions_sample[reaction] + rel_abundance
        } else {
          weighted_reactions_sample[reaction] <- rel_abundance
        }
      }
    }

    weighted_reactions[[sample]] <- weighted_reactions_sample
  }
  return(weighted_reactions)
}


#### FOR REACTIONS
weighted_reactions <- get_weightedByRelAbund_reacs_list(sample_compartment_models_reactions_list)
all_reactions <- unique(unlist(lapply(weighted_reactions, names)))
reaction_matrix <- matrix(0, nrow = length(weighted_reactions), ncol = length(all_reactions),
                          dimnames = list(names(weighted_reactions), all_reactions))

# Fill in the matrix with the weighted reaction counts
for (sample in names(weighted_reactions)) {
  for (reaction in names(weighted_reactions[[sample]])) {
    reaction_matrix[sample, reaction] <- weighted_reactions[[sample]][reaction]
  }
}

# Calculate dissimilarity using Bray-Curtis index
dissimilarity_matrix <- vegdist(reaction_matrix, method = "bray")


# Perform PCoA
pcoa_results <- cmdscale(dissimilarity_matrix, eig = TRUE, k = 2)
eigenvalues <- pcoa_results$eig
# Calculate the percentages of variance explained by each axis
percent_variance_explained <- eigenvalues / sum(eigenvalues) * 100

# Convert to a dataframe for easier handling with ggplot2
pcoa_df <- as.data.frame(pcoa_results$points)
names(pcoa_df) <- c("PC1", "PC2")
pcoa_df$Sample_Compartment <- rownames(pcoa_df)
pcoa_df$Sample <- unlist(lapply(pcoa_df$Sample_Compartment, function(x) strsplit(x, "_")[[1]][1]))
pcoa_df$Compartment <- unlist(lapply(pcoa_df$Sample_Compartment, function(x) strsplit(x, "_")[[1]][2]))
pcoa_df$Diet <- gsub('[0-9]|-', '', pcoa_df$Sample)

color_compartments = c(Gizzard = wesanderson::wes_palette(name = "Darjeeling2")[1],
                       Duodenum = wesanderson::wes_palette(name = "Darjeeling2")[3],
                       Jejunum = wesanderson::wes_palette(name = "Darjeeling1")[4],
                       Ileum = wesanderson::wes_palette(name = "Darjeeling1")[1],
                       Ceca = wesanderson::wes_palette(name = "Darjeeling1")[5],
                       Cecum = wesanderson::wes_palette(name = "Darjeeling1")[5],
                       Colon = wesanderson::wes_palette(name = "Darjeeling2")[2])

pcoa_plot <- ggplot(pcoa_df, aes(x = PC1, y = PC2, fill = Compartment)) +
  geom_point(aes(shape = Diet), size = 4, alpha = 0.7) +
  scale_fill_manual(values = color_compartments) +
  scale_shape_manual(values = c(21, 25)) +
  theme_bw() +
  labs(title = "PCoA of samples by reaction content",
       x = paste("PCoA 1 (", format(percent_variance_explained[1], digits = 2), "%)", sep = ""),
       y = paste("PCoA 2 (", format(percent_variance_explained[2], digits = 2), "%)", sep = "")) +
  theme(legend.title = element_blank(),
        title = element_text(size = 18),
        panel.grid.major = element_blank(),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 16),
        axis.title = element_text(size = 17)) +
  stat_ellipse(aes(group = Compartment, color = Compartment),
               type = "t", linetype = "dashed", size = 0.8, level = 0.95) +
  scale_color_manual(values = color_compartments)

pdf("ModelsAnalyses_outputs/PCoA_10samples_compartments_byReactionContent_weightedRelAbund_BrayCurtis_colorByComp.pdf",
    width = 11, height = 9)
print(pcoa_plot)
dev.off()


#### #### #### #### #### #### #### #### #### ####
#### GET REACTIONS DRIVING THE SEPARATION OF COMPARTMENTS ON PCOA
#### #### #### #### #### #### #### #### #### ####
# Correlation-based approach for reaction loadings
scores <- pcoa_results$points[, 1:2]
PC1 <- scores[, 1]
PC2 <- scores[, 2]

# Spearman correlation of each reaction with PC1/PC2
cor_PC1 <- apply(reaction_matrix, 2, function(x) cor(x, PC1, method = "spearman"))
cor_PC2 <- apply(reaction_matrix, 2, function(x) cor(x, PC2, method = "spearman"))

# Get top loadings
top_PC1 <- sort(cor_PC1, decreasing = TRUE)[1:30]
top_PC1_negative <- sort(cor_PC1, decreasing = F)[1:20]
top_PC2 <- sort(cor_PC2, decreasing = TRUE)[1:20]
top_PC2_negative <- sort(cor_PC2, decreasing = F)[1:20]

reac_names_PC1_neg <- unlist(lapply(names(top_PC1_negative),
                                    function(x) {seed_reactions$name[seed_reactions$id == gsub("_c0", "", x)]}))

reac_names_PC1 <- unlist(lapply(names(top_PC1), function(x) {seed_reactions$name[seed_reactions$id == gsub("_c0", "", x)]}))

reac_names_PC2_neg <- unlist(lapply(names(top_PC2_negative),
                                    function(x) {seed_reactions$name[seed_reactions$id == gsub("_c0", "", x)]}))

reac_names_PC2 <- unlist(lapply(names(top_PC2), function(x) {seed_reactions$name[seed_reactions$id == gsub("_c0", "", x)]}))

# SCFA production: #19 (3-hydroxybutanoyl-CoA oxidoreductase - butyrate pathway), #24-26 (propanoyl-CoA acyltransferases - propionate/butyrate pathways)
# Amino acid metabolism: Heavy presence (#1, 5-6, 11-12, 14, 16-18) - glycine, proline, aspartate
# Electron transport: #3-4, 8-9 (NAD/quinone oxidoreductases)
### Upper GIT enrichment = aerobic metabolism + diverse amino acid processing + electron transport

# Butyrate synthesis: #10-11 (butyryl-CoA dehydrogenase, butanoyl-CoA phosphotransferase)
# Fiber degradation: #12, 15 (galacturonan, arabinofuranosidase - polysaccharide breakdown)
# Fatty acid beta-oxidation: #7-9, 13 (various acyl-CoA dehydrogenases)
### Lower GIT enrichment = fermentation (SCFA synthesis) + polysaccharide degradation + specific amino acid metabolism

top_reactions_df <- data.frame(
  reaction = c(names(top_PC1)[1:15], names(top_PC2)[1:15]),
  loading = c(top_PC1[1:15], top_PC2[1:15]),
  PC = c(rep("PC1", 15), rep("PC2", 15)),
  reaction_name = c(reac_names_PC1[1:15], reac_names_PC2[1:15])
)


######
##### Add compartment info to PC scores
pc_scores_with_compartment <- data.frame(
  Sample = rownames(pcoa_results$points),
  PC1 = pcoa_results$points[, 1],
  PC2 = pcoa_results$points[, 2]
)

pc_scores_with_compartment$Compartment <- unlist(lapply(pc_scores_with_compartment$Sample,
                                                        function(x) strsplit(x, "_")[[1]][2]))

# Check mean PC1 by compartment
pc1_by_compartment <- pc_scores_with_compartment %>%
  group_by(Compartment) %>%
  summarise(mean_PC1 = mean(PC1),
            sd_PC1 = sd(PC1),
            mean_PC2 = mean(PC2),
            sd_PC2 = sd(PC2))

print(pc1_by_compartment)

# Visualize
ggplot(pc_scores_with_compartment, aes(x = factor(Compartment,
                                                  levels = c("Gizzard", "Duodenum", "Jejunum", "Ileum", "Ceca", "Colon")),
                                       y = PC1, fill = Compartment)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_fill_manual(values = color_compartments) +
  theme_bw() +
  labs(x = "", y = "PC1 Score", title = "PC1 distribution by compartment") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

### Ceca/Colon have positive PC1, i.e., positive PC1 loadings = lower GIT enrichment.

# Plot PC2 distribution by compartment
ggplot(pc_scores_with_compartment, aes(x = factor(Compartment, levels = c("Gizzard", "Duodenum", "Jejunum", "Ileum", "Ceca", "Colon")),
                                       y = PC2, fill = Compartment)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_fill_manual(values = color_compartments) +
  theme_bw() +
  labs(x = "", y = "PC2 Score", title = "PC2 distribution by compartment") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Also check means
pc_scores_with_compartment %>%
  group_by(Compartment) %>%
  summarise(mean_PC2 = mean(PC2),
            sd_PC2 = sd(PC2)) %>%
  arrange(mean_PC2)


### Verify Reaction Enrichment in Compartments
calculate_reaction_enrichment <- function(reaction_matrix, top_reactions) {
  reaction_df <- as.data.frame(reaction_matrix)
  reaction_df$Sample <- rownames(reaction_df)
  reaction_df$Compartment <- unlist(lapply(reaction_df$Sample, function(x) strsplit(x, "_")[[1]][2]))

  top_reaction_ids <- names(top_reactions)
  enrichment_results <- list()

  for (reaction_id in top_reaction_ids) {
    if (reaction_id %in% colnames(reaction_df)) {
      comp_means <- reaction_df %>%
        group_by(Compartment) %>%
        summarise(mean_abundance = mean(.data[[reaction_id]], na.rm = TRUE)) %>%
        arrange(desc(mean_abundance))

      enrichment_results[[reaction_id]] <- comp_means
    }
  }

  return(enrichment_results)
}

# Analyze top PC1 reactions
pc1_enrichment <- calculate_reaction_enrichment(reaction_matrix, top_PC1)

# Check which compartments have highest abundance for top reactions
for (i in 1:5) {
  reaction_id <- names(top_PC1)[i]
  cat("\n", i, ". ", reac_names_PC1[i], "\n")
  print(pc1_enrichment[[reaction_id]])
}


######################################################################
##########     FIGURE S3 - REACTION LOADINGS PC1/2      #############
######################################################################
library(pheatmap)
library(viridis)
library(dplyr)
library(tidyr)

# Select top reactions from each PC direction
top_n <- 20

selected_reactions <- c(
  names(top_PC1)[1:top_n],
  names(top_PC1_negative)[1:10],
  names(top_PC2)[1:top_n],
  names(top_PC2_negative)[1:10]
)

# Get unique reactions (in case of overlap)
selected_reactions <- unique(selected_reactions)

# Calculate mean abundance by compartment for selected reactions
heatmap_data <- matrix(NA,
                       nrow = 6,
                       ncol = length(selected_reactions),
                       dimnames = list(
                         c("Gizzard", "Duodenum", "Jejunum", "Ileum", "Ceca", "Colon"),
                         selected_reactions
                       ))

for (i in 1:length(selected_reactions)) {
  reaction_id <- selected_reactions[i]

  if (reaction_id %in% colnames(reaction_matrix)) {
    reaction_df <- data.frame(
      Compartment = unlist(lapply(rownames(reaction_matrix),
                                  function(x) strsplit(x, "_")[[1]][2])),
      Abundance = reaction_matrix[, reaction_id]
    )

    comp_means <- reaction_df %>%
      group_by(Compartment) %>%
      summarise(mean_abund = mean(Abundance, na.rm = TRUE))

    for (comp in comp_means$Compartment) {
      heatmap_data[comp, i] <- comp_means$mean_abund[comp_means$Compartment == comp]
    }
  }
}

# Get reaction names for selected reactions
reaction_labels <- sapply(selected_reactions, function(x) {
  idx <- which(names(top_PC1) == x)
  if (length(idx) > 0) return(reac_names_PC1[idx])

  idx <- which(names(top_PC1_negative) == x)
  if (length(idx) > 0) return(reac_names_PC1_neg[idx])

  idx <- which(names(top_PC2) == x)
  if (length(idx) > 0) return(reac_names_PC2[idx])

  idx <- which(names(top_PC2_negative) == x)
  if (length(idx) > 0) return(reac_names_PC2_neg[idx])

  return(x)
})

colnames(heatmap_data) <- reaction_labels

# Create annotation for which PC/direction
annotation_col <- data.frame(
  PC_Direction = c(
    rep("PC1_Lower_GIT", top_n),
    rep("PC1_Upper_GIT", 10),
    rep("PC2_Ceca_Colon", top_n),
    rep("PC2_Gizzard", 10)
  )[1:length(selected_reactions)]
)

# Make reaction labels unique by appending counter for duplicates
reaction_labels_unique <- make.unique(reaction_labels, sep = "_")

colnames(heatmap_data) <- reaction_labels_unique
rownames(annotation_col) <- reaction_labels_unique

# Define colors for annotation
ann_colors <- list(
  PC_Direction = c(
    "PC1_Lower_GIT" = "#542294",
    "PC1_Upper_GIT" = "#E05536",
    "PC2_Ceca_Colon" = "#4EA0AE",
    "PC2_Gizzard" = "#8B4513"
  )
)

# Plot heatmap
pdf("ModelsAnalyses_outputs/FigureS3_PCoA_reaction_loadings_top20PC1_heatmap_angle90.pdf", width = 16, height = 12)
pheatmap(heatmap_data,
         scale = "column",
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         color = viridis(100, option = "D"),
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         show_colnames = TRUE,
         fontsize_col = 10,
         fontsize_row = 13,
         angle_col = 90,
         main = "Metabolic reaction enrichment across gut compartments\n(Top loadings from PCoA)")
dev.off()


##############################################################################
####################   PERMANOVA - REACTION ABUNDANCES      ##################
##############################################################################

reactions_matrix_abund_byComp <- as.data.frame(reaction_matrix)
reactions_matrix_abund_byComp$Diet <- unlist(lapply(rownames(reactions_matrix_abund_byComp), function(x)
  gsub('[0-9]|-', '', strsplit(x, "_")[[1]][1])))
reactions_matrix_abund_byComp$Compartment <- unlist(lapply(rownames(reactions_matrix_abund_byComp), function(x)
  strsplit(x, "_")[[1]][2]))
reactions_matrix_abund_byComp_dist <-
  vegan::vegdist(as.matrix(reactions_matrix_abund_byComp[, -which(colnames(reactions_matrix_abund_byComp) %in% c("Diet", "Compartment"))]))

reactions_matrix_abund_byComp_permanova <- adonis2(reactions_matrix_abund_byComp_dist ~ Compartment,
                                                    data = reactions_matrix_abund_byComp, permutations = 999, method = "bray")
reactions_matrix_abund_byComp_permanova

reactions_matrix_abund_byCompDiet_permanova <- adonis2(reactions_matrix_abund_byComp_dist ~ Compartment*Diet,
                                                       data = reactions_matrix_abund_byComp, permutations = 999, method = "bray")
reactions_matrix_abund_byCompDiet_permanova

dispersion_reacs_byComp <- vegan::betadisper(reactions_matrix_abund_byComp_dist,
                                             group = reactions_matrix_abund_byComp$Compartment)
dispersion_reacs_byCompDiet <- vegan::betadisper(reactions_matrix_abund_byComp_dist,
                                                 group = paste0(reactions_matrix_abund_byComp$Compartment, "_", reactions_matrix_abund_byComp$Diet))
dispersion_reacs_byDiet <- vegan::betadisper(reactions_matrix_abund_byComp_dist,
                                             group = reactions_matrix_abund_byComp$Diet)

TukeyHSD(dispersion_reacs_byComp)
TukeyHSD(dispersion_reacs_byCompDiet)
TukeyHSD(dispersion_reacs_byDiet)

boxplot(dispersion_reacs_byComp)
boxplot(dispersion_reacs_byCompDiet)
# Permutation test for homogeneity of multivariate dispersions
test_reacsAbund_compartments <- vegan::permutest(dispersion_reacs_byComp)
test_reacsAbund_compartments_diet <- vegan::permutest(dispersion_reacs_byCompDiet)

## plot:
dispersion_reacs_abund_dist_df <- as.data.frame(dispersion_reacs_byComp$distances)
dispersion_reacs_abund_dist_df$SampleName <- rownames(dispersion_reacs_abund_dist_df)
dispersion_reacs_abund_dist_df$Compartment <- reactions_matrix_abund_byComp$Compartment

dispersion_reacs_abund_plot_byComp <- ggplot(dispersion_reacs_abund_dist_df, aes(x = factor(Compartment,
                                                                                            levels = c("Gizzard", "Duodenum", "Jejunum",
                                                                                                       "Ileum", "Ceca", "Colon")),
                                                                                 y = dispersion_reacs_byComp$distances)) +
  geom_boxplot(aes(fill = factor(Compartment,
                                 levels = c("Gizzard", "Duodenum", "Jejunum",
                                            "Ileum", "Ceca", "Colon")))) +
  labs(y = "Distance to centroid", x = "", title = "Beta-dispersion by reaction abundance") +
  scale_fill_manual(values = color_compartments) +
  theme_classic() +
  theme(
    legend.text = element_text(size = 15),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 15),
    axis.text.x = element_blank(),
    legend.title = element_blank(),
    axis.title.x = element_blank())


#### check by reaction presence
reactions_matrix_presence_byComp <- reactions_matrix_abund_byComp
reactions_matrix_presence_byComp[, 1:3699] <- sapply(reactions_matrix_presence_byComp[, 1:3699],
                                                      function(x) { ifelse(x == 0, 0, 1) })
reactions_matrix_presence_byComp_dist <-
  vegan::vegdist(as.matrix(reactions_matrix_presence_byComp[, -which(colnames(reactions_matrix_presence_byComp) %in% c("Diet", "Compartment"))]))

reactions_matrix_presence_byComp_permanova <- adonis2(reactions_matrix_presence_byComp_dist ~ Compartment,
                                                       data = reactions_matrix_presence_byComp, permutations = 999, method = "bray")
reactions_matrix_presence_byComp_permanova

dispersion_reacs_pres_byComp <- vegan::betadisper(reactions_matrix_presence_byComp_dist,
                                                   group = reactions_matrix_presence_byComp$Compartment)

TukeyHSD(dispersion_reacs_pres_byComp)

# Permutation test for homogeneity of multivariate dispersions
test_reacsPres_compartments <- vegan::permutest(dispersion_reacs_pres_byComp)

## plot:
dispersion_reacs_pres_dist_df <- as.data.frame(dispersion_reacs_pres_byComp$distances)
dispersion_reacs_pres_dist_df$SampleName <- rownames(dispersion_reacs_pres_dist_df)
dispersion_reacs_pres_dist_df$Compartment <- reactions_matrix_presence_byComp$Compartment

dispersion_reacs_presence_plot_byComp <- ggplot(dispersion_reacs_pres_dist_df, aes(x = factor(Compartment,
                                                                                              levels = c("Gizzard", "Duodenum", "Jejunum",
                                                                                                         "Ileum", "Ceca", "Colon")),
                                                                                   y = dispersion_reacs_pres_byComp$distances)) +
  geom_boxplot(aes(fill = factor(Compartment,
                                 levels = c("Gizzard", "Duodenum", "Jejunum",
                                            "Ileum", "Ceca", "Colon")))) +
  labs(y = "Distance to centroid", x = "", title = "Beta-dispersion by reaction presence") +
  scale_fill_manual(values = color_compartments) +
  theme_classic() +
  theme(
    legend.text = element_text(size = 15),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 15),
    axis.text.x = element_blank(),
    legend.title = element_blank(),
    axis.title.x = element_blank())


