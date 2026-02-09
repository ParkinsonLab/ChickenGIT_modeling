# Multi-compartment spatiotemporal metabolic modeling of the chicken gut guides dietary intervention design

Irina Utkina, Mohammadali Alizadeh, Shayan Sharif, John Parkinson

## ABSTRACT
Understanding the interactions between diet and the gut microbiome is critical for identifying dietary interventions that support gut health. This is of particular importance for poultry where the elimination of antibiotic growth promoters has resulted in an alarming rise in enteric infections with significant economic consequences. While the application of computational models capable of dissecting the metabolic interactions supporting gut communities has shown promise, they remain limited, largely ignoring the physiological and geographical considerations of the poultry gastrointestinal tract. To address these limitations, we developed the first multi-compartment, spatiotemporally-resolved metabolic model of the chicken gastrointestinal tract. Our six-compartment framework integrates avian-specific physiological features including bidirectional flow, feeding-fasting cycles, and compartment-specific environmental parameters. The model captured distinct metabolic specialization along the gut, with upper compartments enriched for biosynthetic pathways and lower compartments specialized for fermentation. In silico screening of 34 dietary supplements revealed context-dependent metabolic responses and predicted cellulose, starch, and L-threonine as robust enhancers of short-chain fatty acid production. A controlled feeding trial validated key predictions, particularly for butyrate, with prediction accuracy further improved through integration of trial-specific microbial community data. Our findings demonstrate that community composition is a major driver of metabolic outcomes and underscore the need for context-specific modeling. Our framework provides a mechanistic platform for rational dietary intervention design and is broadly adaptable to other animal or human gastrointestinal systems. 


---

## Overview

This repository contains the code, metabolic models, and supplementary data accompanying the manuscript. The project implements a multi-compartment, community-level metabolic modeling framework for the chicken gastrointestinal tract (GIT), simulating microbial metabolism across six GIT compartments (gizzard, duodenum, jejunum, ileum, ceca, colon) using genome-scale metabolic models reconstructed from metagenome-assembled genomes.

The framework couples compartment-specific microbial community simulations via BacArena (an agent-based modeling platform built on flux balance analysis) with compartment-to-compartment metabolite transfer, dietary input, and host absorption. 
A two-compartment (ileum–cecum) model is also provided for validation against independent 16S rRNA survey and metabolomics data from the validation chicken feeding trials.

## Repository Structure

```
.
├── 6compartment_model_GIT/          # Full 6-compartment GIT model
│   ├── ModelGeneration_scripts/     # Scripts to build simulation arenas; includes functions used in the main simulation scripts 
|   ├── ModelGeneration_Files/       # Input files and parameters used in the simulation scripts
│   ├── Simulations_scripts/         # Main simulation scripts (96 h, with/without supplements)
│   ├── AnalysesOutputs_scripts/     # Figure generation, statistical analyses, output extraction
│   └── Input_diets_chicken/         # Diet composition files and conversion script
│
├── 2compartment_validation_model/   # 2-compartment (ileum–cecum) validation model
│   ├── ModelGeneration_scripts/     # Scripts to build simulation arenas; includes functions used in the main simulation scripts 
|   ├── ModelGeneration_Files/       # Input files and parameters used in the simulation scripts
│   ├── AnalysesOutputs_scripts/     # Validation analyses (metabolomics correlation, PCoA)
│   └── Simulation_Ileum_Cecum_*.R   # Simulation script
│
├── Metabolic_models/                # SBML metabolic models
│   ├── final_models_29taxa_6comp/   # 29 species models for the 6-compartment model
│   └── final_models_ileum_cecum/    # 169 MAG-derived models for the 2-compartment model


## Scripts description

### 1. Six-Compartment GIT Model (`6compartment_model_GIT/`)

#### Model Generation (`ModelGeneration_scripts/`)

| Script | Description |
|--------|-------------|
| `Generate_SpecieModelLists_FromTable.R` | Loads gapseq-reconstructed SBML models and assembles per-compartment species lists |
| `Function_ConstructArena_universal.R` | Constructs BacArena simulation arenas with diet, spatial structure, and mucus gradients |
| `Functions_AddBacMetabolites_upd.R` | Adds bacterially produced metabolites to the shared environment |
| `Functions_RemoveBacMetabolites_TypeSpecieAssoc.R` | Removes metabolites based on consumer–species associations |
| `Functions_AbsorptionConversionMetabolites_gapseq_upd.R` | Simulates host absorption and conversion of metabolites between compartments |
| `Function_AddReactionConstraints_toModelList.R` | Applies reaction-level flux constraints to models |
| `Function_RenameModelReactions.R` | Standardizes reaction naming across models |

#### Simulations (`Simulations_scripts/`)

| Script | Description |
|--------|-------------|
| `Simulation_coupling_6comps_args_universal_growerdiet_add24hMets.R` | Control simulation: couples 6 compartments over 96 h with 24 h metabolite cycling |
| `Simulation_coupling_6comps_args_universal_growerdiet_add24hMets_addCompounds.R` | Supplement simulation: same as above with dietary supplement addition |
| `Simulation_coupling_6comps_args_universal_growerdiet_continueAfter48h.R` | Continuation run for controls (extends from 48 h checkpoint) |
| `Simulation_coupling_6comps_args_universal_growerdiet_continueAfter48h_addCompounds.R` | Continuation run for supplement simulations |

#### Analysis and Figures (`AnalysesOutputs_scripts/`)

| Script | Description |
|--------|-------------|
| `Extract_OUTfiles_from_simulations_Controls.R` | Extracts abundance, flux, and concentration data from control simulations |
| `Extract_OUTfiles_from_simulations_withCompounds.R` | Same extraction for supplement simulations |
| `CompileTables_fromOUTfiles_Controls.R` | Aggregates extracted control data into summary tables |
| `CompileTables_fromOUTfiles_withCompounds.R` | Aggregates extracted supplement data |
| `Figure1C_29SpeciesCharacteristics.R` | Characterizes the 29 modeled taxa (reaction/metabolite counts) |
| `Figure1DE_FigureS3_PCoA_Reactions_ECs_content.R` | PCoA of reaction and EC content across species |
| `Figure2_FigureS2_paper_conc_fluxes_PCoA_heatmap_barplots.R` | Metabolite concentration dynamics, flux profiles, PCoA, heatmaps |
| `Figure3_lowerGIT_SCFA_FC_bubblePlot_quantitativeBarChart.R` | SCFA fold-change analysis with dietary supplements |
| `FigureS5_FluxFoldChange_BoxlPlots_96h_compounds_vs_control.R` | Flux fold-change box plots for supplement vs. control |
| `Run_PairedWilcoxTest_FoldChanges_TotalFluxesAllMets_Conc_h96allMets_compound_vs_control.R` | Paired Wilcoxon tests for statistical comparisons |
| `FindAndWrite_uniqueECs_inXMLmodels.py` | Extracts unique EC numbers from SBML model XML files |

#### Diet Input (`Input_diets_chicken/`)

| File | Description |
|------|-------------|
| `Diet_feedtables_design_2023.R` | Converts feed ingredient tables to ModelSEED metabolite concentrations (mM) |
| `Diet_FeedTables_ingredients.xlsx` | Raw ingredient composition of corn- and wheat-based diets |
| `all_seed_metabolites_edited.tsv` | ModelSEED metabolite reference table |

### 2. Two-Compartment Validation Model (`2compartment_validation_model/`)

#### Model Generation (`ModelGeneration_scripts/`)

| Script | Description |
|--------|-------------|
| `Generate_SpecieModelLists_fromMAGs_singleBird_niagara.R` | Builds bird-specific model lists from MAG-to-ASV mapping |
| `Function_ConstructArena_IleumCecum.R` | Constructs the 2-compartment (ileum → cecum) arena |
| `Function_ctFBA_BacArena.R` | Implements L2-regularized context-specific FBA within BacArena |
| `Compile_input_diet_for_ileum.R` | Prepares ileal diet input from upstream compartment outputs |
| `Functions_AbsorptionMetabolites_gapseq.R` | Host absorption functions for the 2-compartment model |
| `Functions_RemoveBacMetabolites_TypeSpecieAssoc.R` | Metabolite removal by consumer associations |
| `Function_AddBacMetabolites_fromUpstream.R` | Adds metabolites transferred from the ileum to the cecum |

#### Simulation

| Script | Description |
|--------|-------------|
| `Simulation_Ileum_Cecum_universal_compounds_L2reg_ctFBA.R` | Main simulation script for the 2-compartment model with L2-regularized ctFBA implementation |

#### Analysis and Figures (`AnalysesOutputs_scripts/`)

| Script | Description |
|--------|-------------|
| `Extract_2comp_SimulationData_upd.R` | Extracts simulation output data from 2-compartment runs |
| `CombineOutputTables_2compartments_batches.R` | Combines outputs across batched simulation runs |
| `Figure4A_CorrelationAnalyses_metabolomics_wet_NORM_vs_predicted_medians.R` | Correlates predicted metabolite concentrations with experimental metabolomics |
| `Figure4B_PCoA_composition_initial16S_simRelAbund_trial16S_ileum_cecum.R` | PCoA comparing predicted community composition to observed 16S profiles |
| `FigureS8_Plot_metabolomicsMets_Fluxes.R` | Supplementary metabolomics vs. predicted flux comparisons |


## Usage

### Running 6-compartment simulations

Simulations are designed to run on an HPC cluster. Each script accepts command-line arguments, e.g.:

```bash
Rscript Simulation_coupling_6comps_args_universal_growerdiet_add24hMets.R \
  <sample_id> <diet_type> <simulation_date>
```

For supplement simulations, additional arguments specify the compound and concentration, e.g.:

```bash
Rscript Simulation_coupling_6comps_args_universal_growerdiet_add24hMets_addCompounds.R \
  <sample_id> <diet_type> <simulation_date> <compound_name> <compound_concentration>
```

### Running 2-compartment validation simulations

```bash
Rscript Simulation_Ileum_Cecum_universal_compounds_L2reg_ctFBA.R \
  <econc_proportion> <diet_proportion> <bird_id> <timesteps> <sim_date> [compound] [compound_conc]
```


