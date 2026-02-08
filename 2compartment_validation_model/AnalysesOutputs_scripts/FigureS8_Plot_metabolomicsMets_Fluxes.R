# Figure S8: stacked bar plots of metabolomics-measured metabolite fluxes
# by family, across Ileum and Cecum compartments and treatment groups.

library(ggplot2)
library(scales)
library(tidyverse)
library(dplyr)
library(tidyr)
library(readr)
library(cowplot)

otu_mag_metadata <- as.data.frame(read_tsv("Chicken_IleumCecum_2024trial_modeling/ModelGeneration_Files/RepresentativeOTUs_matchedMAGs_with_metadata_paired_birds_corrected.tsv"))

path_new_predictions <- "Chicken_IleumCecum_2024trial_modeling/simulation_outputs/"

flux_specie_48h_G1 <- read.table(paste0(path_new_predictions, "Tables/flux_table_sum_metabolomics_bySpecie_meanByRep_G1_updmods_lessIlealBac_sim20250801.txt"),
                                 sep = '\t', header = T)
flux_specie_48h_treatments <- read.table(paste0(path_new_predictions, "Tables/flux_table_sum_metabolomics_bySpecie_meanByRep_treatment_updmods_lessIlealBac_sim20250802.txt"),
                                         sep = '\t', header = T)

# Summarize at the family level
otu_mag_metadata_sub <- otu_mag_metadata %>%
  rename(Sample = bird_id, Specie = species_ID, Compartment = site) %>%
  dplyr::select(Sample, Specie, genus, family, Compartment, group_short)

flux_specie_48h <- bind_rows(flux_specie_48h_G1, flux_specie_48h_treatments)
flux_specie_48h$Specie <- gsub("_ileum_str|_cecum_str", "", flux_specie_48h$Specie)

flux_specie_48h <- flux_specie_48h %>%
  left_join(otu_mag_metadata_sub, by = c("Sample", "Specie")) %>%
  select(-Compartment.y) %>%
  rename(Compartment = Compartment.x)

flux_family_48h <- flux_specie_48h %>%
  group_by(Compartment, Sample, Metabolite_Name, group_short, family) %>%
  summarize(Flux.sum = sum(Flux.mean))

flux_family_48h$Compartment <- factor(flux_family_48h$Compartment, levels = c("Ileum", "Cecum"))

flux_family_48h_treatMean <- flux_family_48h %>%
  group_by(Compartment, Metabolite_Name, group_short, family) %>%
  summarize(Flux.sum_mean = mean(Flux.sum))

CustomColors_families <- c(
  "[Eubacterium]_coprostanoligenes_group" = "#dbd50f",
  "Anaerovoracaceae" = "#ff8c7f", "Aerococcaceae" = "#ba6d4e", "Bacillaceae" = BacArena::colpal3[15],
  "Bacteroidaceae" = BacArena::colpal3[7], "Butyricicoccaceae" = "#a53e76",
  "Christensenellaceae" = "#fff897", "Clostridiaceae" = BacArena::colpal3[12],
  "Clostridia_vadinBB60_group" = "#8c76ab", "Clostridia_UCG-014" = "#a463ff",
  "Comamonadaceae" = BacArena::colpal3[4], "Defluviitaleaceae" = "#FFD4EB",
  "Eggerthellaceae" = "#2dbabd", "Enterobacteriaceae" = BacArena::colpal3[10],
  "Enterococcaceae" = "#cc9616",
  "Erysipelotrichaceae" = "#BA5D00",
  "Erysipelatoclostridiaceae" = "#8a2f0e", "Lachnospiraceae" = "#2A114A",
  "Lactobacillaceae" = BacArena::colpal3[5], "Muribaculaceae" = "#99ad72",
  "Moraxellaceae" = "#dec39b", "Monoglobaceae" = "#196163",
  "Oscillospiraceae" = "#bd7eb2", "Peptostreptococcaceae" = BacArena::colpal3[9],
  "Pseudomonadaceae" = BacArena::colpal3[6], "RF39" = "#700f6b", "Rhizobiaceae" = "#a5abed",
  "Rhizobiaceae_A" = "#a5abed",
  "Pseudonocardiaceae" = "#e6e640",
  "Ruminococcaceae" = BacArena::colpal3[2],
  "Streptococcaceae" = BacArena::colpal3[19], "Staphylococcaceae" = "#3fbf7d",
  "Synergistaceae" = "#f15613",
  "Borkfalkiaceae" = "#3edbde", "UBA1242" = "#489799", "CAG-314" = "#436f70",
  "Cellulomonadaceae" = "#899c90",
  "Helicobacteraceae" = "#e6b681",
  "UBA660" = "#821f7d",
  "Acutalibacteraceae" = "#ab5ab8",
  "Anaerotignaceae" = "#410973",
  "CAG-508" = "#926294",
  "CAG-239" = "#7e87ed",
  "CAG-465" = "#926294",
  "UBA1390" = "#6b18b5",
  "CAG-74" = "#4798d1",
  "Mycobacteriaceae" = "#708775",
  "UBA932" = "#a19087",
  "Burkholderiaceae" = "#0257bf",
  "UBA3700" = "#9fdee0",
  "Planococcaceae" = "#9c2f0e",
  "Bifidobacteriaceae" = "#7e8280",
  "Other" = "black"
)

flux_barplots_scfa <- list()
for (met in c("Acetic acid", "Propionic acid", "Butyric acid",
              "Lactic acid", "Fumaric acid", "Succinic acid",
              "alpha-Ketoglutaric acid")) {
  flux_barplots_scfa[[met]] <-
    ggplot(flux_family_48h_treatMean[which(flux_family_48h_treatMean$Metabolite_Name == met), ],
           aes(x = factor(Compartment, levels = c("Cecum", "Ileum")), y = Flux.sum_mean, fill = family)) +
    geom_bar(position = "stack", stat = "identity", alpha = 0.9) +
    scale_fill_manual(values = CustomColors_families) +
    facet_grid(group_short + Metabolite_Name ~ ., scales = "free") +
    theme_bw() +
    coord_flip() +
    theme(legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(size = 10),
          title = element_text(size = 11),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 13),
          strip.text.x = element_text(size = 16),
          strip.text.y = element_text(size = 10),
          panel.spacing = unit(x = 0.1, units = "lines"),
          axis.ticks = element_blank())
}

flux_plots_bars <- ggpubr::ggarrange(flux_barplots_scfa$`Acetic acid`,
                                     flux_barplots_scfa$`Propionic acid`,
                                     flux_barplots_scfa$`Butyric acid`,
                                     flux_barplots_scfa$`Lactic acid`,
                                     flux_barplots_scfa$`Fumaric acid`,
                                     flux_barplots_scfa$`Succinic acid`,
                                     flux_barplots_scfa$`alpha-Ketoglutaric acid`,
                                     ncol = 1,
                                     common.legend = TRUE, legend = "bottom")
flux_plots_bars

pdf(paste0(path_new_predictions,
           "Plots/FermentationMetsFluxByFamily_BarPlots_coordFlip_separateScales_sumMeanWithinGroup.pdf"),
    height = 30, width = 10)
flux_plots_bars
dev.off()
