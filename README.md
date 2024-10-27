DOSE-L1000
================
Junmin Wang
10/29/2023

## Latest Update (10/27/2024)

I have created a Shiny app that allows users to **visualize** the compound-induced 
gene expression changes and efficacy and potency of 20 million compound-gene 
pairs from the DOSE-L1000 database. The app can be accessed at 
https://www.dosel1000.com

## Introduction

The LINCS L1000 project has collected gene expression profiles for
thousands of compounds across a wide array of concentrations, cell
lines, and time points. However, conventional analysis methods often
fall short in capturing the rich information encapsulated within the
L1000 transcriptional dose-response data.

We present DOSE-L1000, a database that unravels the intricate landscape
of compound-induced transcriptional changes and compound-gene
interactions. By fitting over 140 million generalized additive models
and robust linear models encompassing the entire LINCS L1000 database,
our work provides quantitative insights into differential gene
expression and the potency and efficacy of compound-gene pairs across
diverse cellular contexts.

Creation of the DOSE-L1000 database comprises a series of steps: model
fitting, differential expression analysis, and potency/efficacy
calculation. This repository contains the scripts used to generate the
database as well as the plots and analysis in the paper.

## Prerequisite

To execute the scripts used to generate the DOSE-L1000 database, make
sure that the original LINCS L1000 data are downloaded. The complete
list of files that need to be downloaded is provided in Supplementary
Table 1. To execute the scripts used to generate the plots, ensure that
the DOSE-L1000 database is downloaded from Zenodo.org. All R packages
that need to be installed are listed in the header of each script.

## Creation of DOSE-L1000

Scripts inside the “model” folder are used to generate the DOSE-L1000
database.

- main_GSE92742.R conducts differential expression analysis for all
  combinations of compounds, cell lines, and genes in the LINCS L1000
  phase 1 data.
- main_GSE70138.R conducts differential expression analysis for all
  combinations of compounds, cell lines, and genes in the LINCS L1000
  phase 2 data.
- main_lx50_GSE92742.R characterizes the efficacy and potency of
  compound-gene interactions for all combinations of compounds, cell
  lines, and genes in the LINCS L1000 phase 1 data.
- main_lx50_GSE70138.R characterizes the efficacy and potency of
  compound-gene interactions for all combinations of compounds, cell
  lines, and genes in the LINCS L1000 phase 2 data.

Scripts inside the “funs” folder are required to execute the scripts
inside other folders.

- get.gam.tests.R conducts t-tests to compare the mean responses for
  generalized additive models.
- get.rlm.tests.R conducts t-tests to compare the mean responses for
  robust linear models.
- get.gam.stats.R calculates the estimates and standard errors of
  efficacy and potency.
- get.lm.tests.R conducts t-tests to compare the mean responses for
  linear models.
- get.rlmer.log2fc.R estimates log2 fold change for robust linear mixed
  models.
- get.gamer.log2fc.R estimates log2 fold change for generalized additive
  mixed models.
- get.gam.stats.boot.R calculates the estimates and standard errors of
  efficacy and potency via bootstrapping.

## Plotting

Scripts inside the “plot” folder are used to generate the plots in the
paper.

- plot_volcano.R makes the volcano plots in Figure 3a.
- plot_potency_vs_efficacy.R makes the scatter plots in Figures 3b and
  3c.
- plot_heatmap.R makes the the heat map in Figure 4.

## Assessment of Normality Assumption

Scripts inside the “normal_dist” folder are used to evaluate the
normality assumption.

- check_norm_assump_GSE92742.R calculates the Pearson correlation
  coefficients from the Q-Q plots for all combinations of compounds,
  cell lines, and time in the phase 1 data.
- check_norm_assump_GSE70138.R calculates the Pearson correlation
  coefficients from the Q-Q plots for all combinations of compounds,
  cell lines, and time in the phase 2 data.

## Examination of Inter- and Intra-Batch Correlations

Scripts inside the “robust” folder are used to evaluate the robustness
of the models.

- find_common_cpd_cell_pairs.R finds repeated compound-cell line-project
  combinations.
- fit_all_models_GSE92742.R fits GAMs, RLMs, and LMs and runs
  differential expression analysis for all repeated combinations of
  compounds and cell lines in the phase 1 data.
- fit_all_models_GSE70138.R fits GAMs, RLMs, and LMs and runs
  differential expression analysis for all repeated combinations of
  compounds and cell lines in the phase 2 data.
- calc_sd_log2fc.R calculates standard deviation of log2 fold changes
  for all repeated perturbation conditions.
- calc_inter_batch_corr.R calculates pairwise inter-batch Pearson
  correlation coefficients.
- calc_intra_batch_corr.R calculates pairwise intra-batch Pearson
  correlation coefficients.

## Examination of Plate Effects

Scripts inside the “plate” folder are used to evaluate plate effects.

- fit_mixed_model_GSE92742.R fits mixed models and estimates log2 fold
  change for all combinations of compounds and cell lines in the phase 1
  data.
- fit_mixed_model_GSE70138.R fits mixed models and estimates log2 fold
  change for all combinations of compounds and cell lines in the phase 2
  data.
- compare_fixed_vs_mixed.R calculates the Pearson correlation between
  the results of fixed models and the results of mixed models.

## Comparison between Delta Method and Bootstrapping

Scripts inside the “delta_vs_bootstrap” folder are used to compare the
delta method against bootstrapping.

- compare_delta_vs_boot.R compares the delta method against
  bootstrapping.

## Predictions of Drug-Target Interactions Based on Gene Expression Similarities

Scripts inside the “LOOCV” folder are used to predict drug-target
interactions based on expression similarities.

- find_HDACi_DTC.R retrieves HDAC inhibitors from Drug Target Commons.
- LOOCV_HDACi.R builds a classifier for HDAC inhibitors and runs
  leave-one-out cross-validation (LOOCV) on the classifier.
