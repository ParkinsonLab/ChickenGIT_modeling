# Figure 1C: Species characteristics for 29 GEM models.
# PCoA by EC content, EC sharing analysis, phylogenetic tree, compartment presence
# dot plots, bar plots for reactions/ECs/pathways, KEGG pathway mapping, and
# gapseq subsystems pathway diversity analysis.

library(taxize)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ape)
library(BacArena)

##########################################################################################
################      PCOA FOR REACS / ECs FOR EACH MODEL            #####################
##########################################################################################
ECs_models <- read.delim("ModelsAnalyses_outputs/Table_ECnumbers_inEachModel_gapseqAdapted_cobrapyAdapted_July2023.txt", sep='\t', header=T)
models_reacs_list <- readRDS("ModelsAnalyses_outputs/AllModels_reacs_list.RDS")

CustomColors_species <- c("Blautia_sp" = "#542294", "Mediterraneibacter_sp002161355" = "#412566", "Lachnoclostridium_phocaeensis" = "#2b1e3d",
                          "Clostridium_saudiense" = colpal3[12], "Faecalibacterium_sp002160895" = "#d95dc4",
                          "Agathobaculum_sp" = "#9c3b8b", "Flavonifractor_plautii" = "#8f1d7b",
                          "Oscillospiraceae_genusUBA9475" = "#73015f",
                          "Romboutsia_sp" = colpal3[9], "Borkfalkia_sp" = "#FFD4EB", "Lactobacillus_crispatus" = "#f72548",
                          "Limosilactobacillus_reuteri" = "#bd112e", "Lactobacillus_johnsonii" = "#99021b",
                          "Streptococcus_equinus" = colpal3[19], "Staphylococcus_saprophyticus" = "#f58971",
                          "Enterococcus_cecorum" = "#e05536", "Enterococcus_faecalis" = "#ab4d38",
                          "Erysipelatoclostridium_sp002160495" = "#BA5D00", "Massiliomicrobiota_sp002160815" = "#8a4601",
                          "Megasphaera_stantonii" = colpal3[15], "Escherichia_flexneri" = "#0569ab", "Klebsiella_pneumoniae" = "#044d7d",
                          "Sphingomonas_sp" = colpal3[4],
                          "Acinetobacter_lwoffii" = "#e6ba7a", "Pseudomonas_aeruginosa" = "#b08b56",
                          "Bifidobacterium_gallinarum" = "#bfb5ae",
                          "Corynebacterium_amycolatum" = "#1b7a3c", "Corynebacterium_kroppenstedtii" = "#1d4f2f",
                          "Phocaeicola_vulgatus" = "#595653")

ECs_models_list <- list()
for (specie in ECs_models$SpecieName) {
  ECs <- unlist(strsplit(ECs_models$EC_numbers[ECs_models$SpecieName == specie], ", ")[[1]])
  ECs_models_list[[specie]] <- ECs
}

all_ECs <- unique(unlist(ECs_models_list))

EC_models_matrix <- matrix(0, nrow = length(ECs_models_list), ncol = length(all_ECs),
                           dimnames = list(names(ECs_models_list), all_ECs))

# Fill in the matrix with the reac presence/absence
for (specie in names(ECs_models_list)) {
  for (EC in ECs_models_list[[specie]]) {
    EC_models_matrix[specie, EC] <- 1
  }
}

dissimilarity_matrix <- vegdist(EC_models_matrix, method = "bray")

# Perform PCoA
pcoa_results <- cmdscale(dissimilarity_matrix, eig = TRUE, k = 2)
eigenvalues <- pcoa_results$eig

# Calculate the percentages of variance explained by each axis
percent_variance_explained <- eigenvalues / sum(eigenvalues) * 100
# Convert to a dataframe for easier handling with ggplot2
pcoa_df <- as.data.frame(pcoa_results$points)
names(pcoa_df) <- c("PC1", "PC2")
pcoa_df$taxon <- rownames(pcoa_df)

pcoa_plot <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = taxon)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = CustomColors_species) +
  theme_bw() +
  labs(title = "PCoA of models by EC content",
       x = paste("PCoA 1 (", format(percent_variance_explained[1], digits = 2), "%)", sep = ""),
       y = paste("PCoA 2 (", format(percent_variance_explained[2], digits = 2), "%)", sep = "")) +
  theme(legend.title = element_blank(),
        title = element_text(size = 18),
        legend.position = "none",
        panel.grid.major = element_blank(),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 16),
        axis.title = element_text(size = 17)) +
  ggrepel::geom_text_repel(aes(label = taxon),
                           size = 5, seed = 43,
                           max.overlaps = 10, force = 15)

pdf("ModelsAnalyses_outputs/Plots/PCoA_models_byECcontent_BrayCurtis_colorBySpecie_Dec2024.pdf",
    width = 9, height = 8)
print(pcoa_plot)
dev.off()



#######################################################################################################
#######    LIST OF MODELS WITH PHYLOGENETIC TREE AND PRESENCE/ABSENCE IN EACH GUT COMPARTMENT  #########
########################################################################################################
Compartments_models_noAGP_D10_gapseq <- readxl::read_xlsx("ModelGeneration_Files/SpecieModels_present_in_compartments_D10_noAGP_2022upd.xlsx",
                                                          sheet = "Gapseq_ModelList_upd")

############ MAKE PHYLOGENETIC TREE
classifications_specie <- classification(Compartments_models_noAGP_D10_gapseq$Specie_ncbi_name, db = "ncbi")

tree_result <- class2tree(unique(classifications_specie))
tree <- tree_result$phylo
plot.phylo(tree, cex = 0.7, no.margin = TRUE)

# save order of the species for further plots that will be panels on the right
tree <- ladderize(tree)
plot(tree)
# get the plotting info from 'ape'
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
# Determine the order from last_plot.phylo
N <- length(tree$tip.label)
# Order tips by their vertical coordinates
ordered_indices <- order(lastPP$yy[1:N], decreasing = TRUE)

# re-plot the tree with colored species and labeled as their names:
tree$tip.label <- Compartments_models_noAGP_D10_gapseq$Specie
# Get the tip labels (species) in the plot order (top to bottom)
specie_order <- tree$tip.label[ordered_indices]
print(specie_order)

plot_color_vector <- CustomColors_species[tree$tip.label]
pdf("ModelsAnalyses_outputs/Plots/PhylogeneticTree_29species_coloredBySpecie.pdf",
    height = 6, width = 4)
plot(tree, tip.color = plot_color_vector)
dev.off()

############ MAKE DOT PLOT WITH PRESENCE/ABSENCE IN COMPARTMENTS
specie_presence_comp <-
  Compartments_models_noAGP_D10_gapseq[, c("Specie", "Gizzard", "Duodenum", "Jejunum", "Ileum", "Ceca", "Colon")] %>%
  melt(id = "Specie") %>%
  rename("Compartment" = "variable")

ggplot(specie_presence_comp, aes(x = factor(Compartment),
                                 y = factor(Specie, levels = rev(specie_order)))) +
  geom_point(aes(fill = ifelse(value == "yes", "darkblue", "white")),
             shape = 21,
             size = 5,
             stroke = 1) +
  scale_fill_manual(values = c("darkblue" = "black", "white" = "white")) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    strip.text = element_text(size = 12),
    legend.position = "none"
  ) +
  labs(y = "", x = "")

############# DOT PLOT COLORED BY AVERAGE RELATIVE ABUNDANCE IN COMPARTMENTS (GRADIENT)
abund_table_10samples_allcomp <- read.table("ModelsAnalyses_outputs/RelativeAbundanceTable_allSamples_bothDiets_D10_noAGP.txt",
                                            sep = '\t', header = T)
relabund_avg_specie_allcomp <- abund_table_10samples_allcomp %>%
  group_by(Specie, Compartment, Diet) %>%
  summarise(RelAbundance_specie_mean = mean(RelAbundance_specie))

specie_presence_comp_withMeanAbund <- merge(specie_presence_comp, relabund_avg_specie_allcomp, by = c("Specie", "Compartment"), all.y = T)

# dots split per diet
library(ggnewscale)
corn_data <- specie_presence_comp_withMeanAbund %>% filter(Diet == "Corn")
wheat_data <- specie_presence_comp_withMeanAbund %>% filter(Diet == "Wheat")

dotplot_2diets <- ggplot() +
  # Corn diet layer (darkblue gradient) - positioned slightly to the left
  geom_point(data = corn_data,
             aes(x = factor(Compartment),
                 y = factor(Specie, levels = rev(specie_order)),
                 fill = RelAbundance_specie_mean),
             shape = 21,
             size = 5,
             stroke = 1,
             position = position_nudge(x = -0.1)) +
  scale_fill_gradient(
    low = "white",
    high = "darkblue",
    na.value = "white",
    trans = "sqrt",
    breaks = c(0, 1e-03, 1e-02, 0.05, 0.10, 0.15, 0.20),
    labels = scales::scientific,
    name = "Corn\nRelative\nAbundance"
  ) +
  new_scale_fill() +
  # Wheat diet layer (brown gradient) - positioned slightly to the right
  geom_point(data = wheat_data,
             aes(x = factor(Compartment),
                 y = factor(Specie, levels = rev(specie_order)),
                 fill = RelAbundance_specie_mean),
             shape = 21,
             size = 5,
             stroke = 1,
             position = position_nudge(x = 0.1)) +
  scale_fill_gradient(
    low = "white",
    high = "#6e2e00",  # brown
    na.value = "white",
    trans = "sqrt",
    breaks = c(0, 1e-03, 1e-02, 0.05, 0.10, 0.2, 0.30),
    labels = scales::scientific,
    name = "Wheat\nRelative\nAbundance"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    strip.text = element_text(size = 12),
    legend.position = "bottom",
    legend.text = element_text(angle = 45, hjust = 1, vjust = 1.2),
    legend.title.position = "top"
  ) +
  labs(y = "", x = "")

pdf("ModelsAnalyses_outputs/Plots/DotPlot_29species_coloredByRelAbund_perDiet.pdf",
    height = 8, width = 6)
plot(dotplot_2diets)
dev.off()


##########################################################################################
#######     BARPLOTS FOR NUMBERS OF REACS, ECs, PATHWAYS IN EACH MODEL      ##############
##########################################################################################

### plot number of pathways
SubSystems_models <- readRDS("ModelsAnalyses_outputs/SubsystemsList_eachModel_D10_noAGP_adapted_Jul2023_min3reacs.RDS")
path_counts <- as.data.frame(sapply(SubSystems_models, length))
path_counts$specie <- rownames(path_counts)
path_counts$specie <- gsub('_rect|_genomic|_ileum', '', path_counts$specie)
path_counts$specie <- gsub('_sp_cec', '_sp', path_counts$specie)
colnames(path_counts)[1] <- "path_count"

pathways_n <- ggplot(path_counts, aes(x = factor(specie, levels = rev(specie_order)),
                                      y = path_count)) +
  geom_bar(position = "stack", stat = "identity", fill = "#4EA0AE") +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 16),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size = 14),
        legend.position = "bottom",
        legend.title = element_blank()) +
  scale_y_continuous(breaks = c(0, 500)) +
  coord_flip()

### plot number of ECs
EC_counts <- as.data.frame(sapply(ECs_models_list, length))
EC_counts$specie <- rownames(EC_counts)
colnames(EC_counts)[1] <- "EC_count"

ECs_n <- ggplot(EC_counts, aes(x = factor(specie, levels = rev(specie_order)),
                               y = EC_count)) +
  geom_bar(position = "stack", stat = "identity", fill = "#765396") +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 16),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size = 14),
        legend.position = "bottom",
        legend.title = element_blank()) +
  scale_y_continuous(breaks = c(0, 1000)) +
  coord_flip()

### plot number of reactions
reacs_counts <- as.data.frame(sapply(models_reacs_list, length))
reacs_counts$specie <- rownames(reacs_counts)
colnames(reacs_counts)[1] <- "reacs_count"

reacs_n <- ggplot(reacs_counts, aes(x = factor(specie, levels = rev(specie_order)),
                                    y = reacs_count)) +
  geom_bar(position = "stack", stat = "identity", fill = "#107478") +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 16),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size = 14),
        legend.position = "bottom",
        legend.title = element_blank()) +
  scale_y_continuous(breaks = c(0, 2000)) +
  coord_flip()

### combine them all together:
pdf("ModelsAnalyses_outputs/Plots/Figure1C_29models_NofECreacsPathways_vDec2024.pdf",
    width = 5, height = 6)
cowplot::plot_grid(ECs_n, reacs_n, pathways_n, rel_widths = c(1, 1, 1), nrow = 1)
dev.off()


