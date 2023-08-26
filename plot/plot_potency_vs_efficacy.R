
## This script makes the scatter plots in Figures 3b and 3c.
## Make sure to update the paths to interaction.rds, combination.rds, GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz, and GSE70138_Broad_LINCS_pert_info_2017-03-06.txt.gz

## load library
library(ggplot2)
library(ggrepel)
library(dplyr)
library(latex2exp)

## read data
interaction <- readRDS('path/to/interaction.rds')
combination <- readRDS("path/to/combination.rds")

## read metadata
gene_info_phase2 <- read.delim("path/to/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz",
                               sep = "\t")
pert_info_phase2 <- read.delim("path/to/GSE70138_Broad_LINCS_pert_info_2017-03-06.txt.gz",
                               sep = "\t")

## extract phase 2 combinations; add compound names
combination_phase2_ann <- combination %>%
  filter(phase == 2) %>%
  inner_join(pert_info_phase2, by = "pert_id")

## extract data (MCF-7, CTSD)
interaction_mcf7_ctsd <- interaction %>% 
  inner_join(combination_phase2_ann,
             by = "comb_index") %>%
  filter(gene == 1509,  ## CTSD
         cell_id == "MCF7",
         pert_iname %in% c("tamoxifen",
                           "raloxifene",
                           "fulvestrant",
                           "bazedoxifene",
                           "toremifene",
                           "clomifene")) %>%
  mutate(efficacy = 2^lefficacy, 
         se_efficacy = 2^se_lefficacy,
         potency = 10^lpotency - pseudo_conc,
         se_potency = 10^se_lpotency)

## plot potency vs efficacy
p <- ggplot(interaction_mcf7_ctsd,
            aes(x = potency, 
                y = efficacy * 100,
                color = pert_iname,
                fill = pert_iname,
                label = pert_iname)) +
  geom_point() +
  geom_errorbarh(aes(xmin = potency / se_potency,
                     xmax = potency * se_potency),
                 alpha = 0.5) +
  geom_errorbar(aes(ymin = efficacy / se_efficacy * 100,
                    ymax = efficacy * se_efficacy * 100),
                alpha = 0.5) +
  geom_rect(aes(xmin = potency / se_potency,
                xmax = potency * se_potency,
                ymin = efficacy / se_efficacy * 100,
                ymax = efficacy * se_efficacy * 100),
            alpha = 0.2,
            color = NA) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  labs(x="IC50 (Unit: uM)",
       y="% Maximum Inhibition of CTSD") +
  geom_text_repel(size = 3.5) +
  ggtitle("MCF-7")

## extract data (MCF-7, tamoxifen)
interaction_tamoxifen_mcf7 <- interaction %>% 
  inner_join(combination_phase2_ann,
             by = "comb_index") %>%
  filter(gene %in% c(1509, 595, 4609,
                     332, 5111),
         cell_id == "MCF7",
         pert_iname == "tamoxifen") %>%
  inner_join(gene_info_phase2 %>%
               mutate(pr_gene_id = as.character(pr_gene_id)), 
             by = c("gene" = "pr_gene_id")) %>%
  mutate(efficacy = 2^lefficacy, 
         se_efficacy = 2^se_lefficacy,
         potency = 10^lpotency - pseudo_conc,
         se_potency = 10^se_lpotency)  

## plot potency vs efficacy
p <- ggplot(interaction_tamoxifen_mcf7,
            aes(x = potency, 
                y = efficacy * 100,
                color = pr_gene_symbol,
                fill = pr_gene_symbol,
                label = pr_gene_symbol)) +
  geom_point() +
  geom_errorbarh(aes(xmin = potency / se_potency,
                     xmax = potency * se_potency),
                 alpha = 0.5) +
  geom_errorbar(aes(ymin = efficacy / se_efficacy * 100,
                    ymax = efficacy * se_efficacy * 100),
                alpha = 0.5) +
  geom_rect(aes(xmin = potency / se_potency,
                xmax = potency * se_potency,
                ymin = efficacy / se_efficacy * 100,
                ymax = efficacy * se_efficacy * 100),
            alpha = 0.2,
            color = NA) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  labs(x="IC50 (Unit: uM)",
       y="% Maximum Inhibition in MCF-7 Cells") +
  geom_text_repel(size = 3.5) +
  ggtitle("Tamoxifen")
