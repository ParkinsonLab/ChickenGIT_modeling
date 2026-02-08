
# SCFA flux fold-change analysis for 96h compound simulations vs control.
# Generates boxplots for Supplementary Figure 5

library(cowplot)
library(ggplot2)
library(dplyr)
library(reshape2)
library(magrittr)
library(ggpp)
library(tidyr)
library(ggrepel)
      
## load flux tables with wilcoxon test pvalues:
setwd("WithPrebiotics/2024/")

path = "Tables/"

pattern_arg_compounds = "diet_grower_96h_econc1_propGizzard_corn180_wheat130_CobaltX12_rmAnaero_less02_mucUrea_wCoA_moreNightFlow"
pattern_arg_control = "diet_grower_96h_econc1_propGizzard_corn180_wheat130_CobaltX12_rmAnaero_less02_mucUreaCecaColon_wCoA_moreNightFlow_fix"
sim_dates_compounds = "corn20240721_wheat20240721"
sim_dates_control = "corn20240719_wheat20240719"

Compartments_list <- list(c("Gizzard","Duodenum","Jejunum","Ileum","Cecum", "Colon"))

# load compound fluxes
flux_table_allReps_96h <- read.delim(paste0(path,"FluxTable_ToPlot_PairedWilcoxPvalsBH_FoldChange_allCompounds_",
                                            sim_dates_compounds,"_vs_control_", sim_dates_control, ".txt"),
                                              sep ='\t', header=T)

# subset to SCFA only:
flux_table_allReps_96h_scfa <- flux_table_allReps_96h[which(flux_table_allReps_96h$Metabolite %in% c("Acetate","Propionate","Butyrate")),]

# remove L_Threonine_3 and 4:
flux_table_allReps_96h_scfa <- flux_table_allReps_96h_scfa[which(!(flux_table_allReps_96h_scfa$compound_amount %in% c("L-Threonine_3", "L-Threonine_4"))),]

scfa_colors <- setNames(c("#ab3613","#576991","#702963"), c("Acetate","Propionate","Butyrate"))

# plot fold changes of selected compounds for each compound_amount for each replicate, with dots colored by changes in commsize
plotselectedMets_foldChange_boxplots <- function(flux_table,
                                                 plot_name, logScale = F,
                                                 selectedMets = c("Acetate","Butyrate","Propionate"),
                                                 sim_date, orientation, wilcox_pval = FALSE,
                                                 adjust_asterisk = FALSE,
                                                 lowerGIT_only = FALSE,
                                                 compartment=NULL,
                                                 height_def=13, width_def=7) {

  flux_table$pval_sign <- c("***", "**", "*", "")[findInterval(flux_table$pval_adj, c(0.001, 0.01, 0.05)) + 1]

  selectedMets_plots <- list()
  for (met in selectedMets) {

    df = flux_table[which(flux_table$Metabolite == met),]

    if (lowerGIT_only == TRUE) {
      df <- df[which(df$Compartment %in% c("Ileum","Cecum","Colon")),]
      df$Compartment <- factor(df$Compartment, levels = c("Ileum","Cecum","Colon"))

    } else if (!(is.null(compartment))) {
      df <- df[which(df$Compartment == compartment),]
    }

    if (orientation == "horizontal") { col_number = min(6,length(unique(df$Compartment)))  }
    if (orientation == "vertical") { col_number = 1 }

    if (wilcox_pval == TRUE & adjust_asterisk == TRUE) {
      # to plot only asterisk on top of the maximum FC point, add another column:
      df$pval_sign_max <- ""
      max_FoldChanges_subdf <- df %>% group_by(Metabolite, compound_amount, Compartment) %>% top_n(1, FoldChange)
      # now check which rows from bigger df exist in subset with max foldChanges:
      df$is.max <- do.call(paste0, df) %in% do.call(paste0, max_FoldChanges_subdf)
      # assign pval_sign_max to those with TRUE in "is.max" column:
      df$pval_sign_max[df$is.max == TRUE] <-
        df$pval_sign[df$is.max == TRUE]
    }

    plot <- ggplot(data = df,
                   aes(x = reorder(compound_amount, FoldChange, median),
                       y = FoldChange, colour = Metabolite )) +
      geom_boxplot(
        position = position_dodge(width=0.9),
        outlier.alpha=0.7,
        outlier.size = 0.9,
        lwd=0.8,
        fatten=2.5) +
      theme_bw() +
      geom_hline(yintercept=1, linetype="dashed", color = "darkred", size = 1.5) +
      ylab(label = paste0(met, " production FoldChange"))+
      facet_wrap( ~ Compartment , scales = "free", ncol=col_number) +
      theme(legend.background=element_rect(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            axis.text.x=element_text(angle=0, size=25),
            axis.text.y=element_text(angle=0, size=21),
            axis.title=element_text(size=14),
            axis.title.y = element_blank(),
            axis.title.x = element_text(size = 30),
            strip.text = element_text(size = 28),
            legend.position="none",
            legend.key=element_blank(),
            legend.title=element_blank(),
            legend.text = element_text(size=14)) +
      coord_flip() +
      scale_color_manual(values = scfa_colors)

    if (wilcox_pval == TRUE & adjust_asterisk == TRUE) {
      plot <- plot + geom_text(aes(label = pval_sign_max),
                               position = position_dodge(width = .9), angle=270, vjust = 0.3, size = 9, inherit.aes = T)
    }
    if (wilcox_pval == TRUE & adjust_asterisk == FALSE) {
      plot <- plot + geom_text(aes(label = pval_sign),
                               position = position_dodge(width = .9), vjust = -0.05, size = 8, inherit.aes = T)
    }

    selectedMets_plots[[met]] <- plot

    plot_name_met <-  paste0(met, plot_name, as.character(orientation),"_", sim_date, ".pdf")

    pdf(paste0(path, "../Plots/",
               plot_name_met), width = width_def, height=height_def)
    plot(selectedMets_plots[[met]])
    dev.off()

  }

}

plot_name_allCompounds_SCFA_lowerGIT_allReps <- "FoldChange_BoxPlot_allCompounds_allReps_lowerGIT_withWilcoxPvals_"

plotselectedMets_foldChange_boxplots(flux_table = flux_table_allReps_96h_scfa,
                                     plot_name = plot_name_allCompounds_SCFA_lowerGIT_allReps,
                                     sim_date=sim_dates_compounds,
                                     orientation = "horizontal", wilcox_pval = TRUE,
                                     adjust_asterisk = TRUE,
                                     lowerGIT_only = T,
                                     height_def=17, width_def=40)

