DOSE-L1000
================
Junmin Wang
8/14/2023

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
database as well as the plots in the paper.

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
inside the “model” folder.

- get.gam.tests.R conducts t-tests to compare the mean responses for
  generalized additive models.
- get.rlm.tests.R conducts t-tests to compare the mean responses for
  robust linear models.
- get.gam.stats.R calculates the estimates and standard errors of
  efficacy and potency.

## Plotting

Scripts inside the “plot” folder are used to generate the plots in the
paper.

- plot_volcano.R makes the volcano plots in Figure 3a.
- plot_potency_vs_efficacy.R makes the scatter plots in Figures 3b and
  3c.
- plot_heatmap.R makes the the heat map in Figure 4.
