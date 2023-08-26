
## This script makes the volcano plots in Figure 3a.
## Make sure to update the paths to test.rds, model.rds, combination.rds, and GSE70138_Broad_LINCS_pert_info_2017-03-06.txt.gz

## load libraries
library(ggplot2)
library(dplyr)
library(latex2exp)

## define the labeling function
pert_dose_labeller <- function(x) paste(x, "uM vs DMSO")

## read data
test <- readRDS("path/to/test.rds")
model <- readRDS("path/to/model.rds")
combination <- readRDS("path/to/combination.rds")

## read metadata
pert_info_phase2 <- read.delim("path/to/GSE70138_Broad_LINCS_pert_info_2017-03-06.txt.gz",
                               sep = "\t")

## extract phase 2 combinations; add compound names
combination_phase2_ann <- combination %>%
  filter(phase == 2) %>%
  inner_join(pert_info_phase2, by = "pert_id")

## extract data for Tamoxifen, MCF7; adjust p-values
test_tamoxifen_mcf7 <- test %>%
  filter(comb_index %in% (combination_phase2_ann %>%
            filter(pert_iname == "tamoxifen",
                   cell_id == "MCF7") %>%
            pull(comb_index)),
         pert_time == 24) %>%
  mutate(pert_dose = round(pert_dose, 2)) %>%
  group_by(pert_dose) %>%
  mutate(adj.pval = p.adjust(pval, method = "fdr"))

## make volcano plot
ggplot(data = test_tamoxifen_mcf7) +
  geom_point(size = 0.2,
             aes(x = ifelse(abs(Diff) > 3,
                            3 * sign(Diff), Diff),
                 y = ifelse(-log10(adj.pval) > 30, 
                            30, -log10(adj.pval)), 
                 col = ifelse(adj.pval <= 0.05, 
                            "p < 0.05", 
                            "p > 0.05"))) +
  facet_wrap(. ~ pert_dose, ncol = 2,
             labeller = as_labeller(pert_dose_labeller)) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text = element_text(color = "black"),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(limits = c(-3.2, 3.2)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = TeX("Log$_2$ (Fold Change)"),
       y = TeX("-Log$_{10}$ (Adj. P-Value)"),
       color = "") +
  scale_color_manual(breaks = c( "p < 0.05", 
                                 "p > 0.05"),
                     values = c("red", "gray")) +
  geom_hline(yintercept = -log10(0.05), 
             color = "black", linetype = "dashed",
             size = 0.5) +
  ggtitle("Tamoxifen, MCF-7")

## extract data for Tamoxifen, A375; adjust p-values
test_tamoxifen_a375 <- test %>%
  filter(comb_index %in% (combination_phase2_ann %>%
                            filter(pert_iname == "tamoxifen",
                                   cell_id == "A375") %>%
                            pull(comb_index)),
         pert_time == 24) %>%
  mutate(pert_dose = round(pert_dose, 2)) %>%
  group_by(pert_dose) %>%
  mutate(adj.pval = p.adjust(pval, method = "fdr"))

## make volcano plot
ggplot(data = test_tamoxifen_a375) +
  geom_point(size = 0.2,
             aes(x = ifelse(abs(Diff) > 3,
                            3 * sign(Diff), Diff),
                 y = ifelse(-log10(adj.pval) > 30, 
                            30, -log10(adj.pval)), 
                 col = ifelse(adj.pval <= 0.05, 
                            "p < 0.05", 
                            "p > 0.05"))) +
  facet_wrap(. ~ pert_dose, ncol = 2,
             labeller = as_labeller(pert_dose_labeller)) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text = element_text(color = "black"),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(limits = c(-3.2, 3.2)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = TeX("Log$_2$ (Fold Change)"),
       y = TeX("-Log$_{10}$ (Adj. P-Value)"),
       color = "") +
  scale_color_manual(breaks = c( "p < 0.05", 
                                 "p > 0.05"),
                     values = c("red", "gray")) +
  geom_hline(yintercept = -log10(0.05), 
             color = "black", linetype = "dashed",
             size = 0.5) +
  ggtitle("Tamoxifen, A375")
