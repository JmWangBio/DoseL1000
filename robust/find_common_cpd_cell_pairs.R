
## This script finds repeated compound-cell line-project combinations.
## Make sure to update the paths to combination.rds and model.rds.

## load library
library(dplyr)

## read database
combination <- readRDS("path/to/combination.rds")
model <- readRDS("path/to/model.rds")

## find all unique combinations
combination_gam <- combination %>% 
  inner_join(model, by = "comb_index") %>%
  filter(model == "gam") %>%
  dplyr::select(pert_id, 
                group_id, 
                cell_id,
                phase) %>% 
  distinct() 

## find repeated compound-cell line pairs
repeated_cpd_cell_pairs <- combination_gam %>%
  group_by(pert_id, cell_id) %>%
  summarise(num = n()) %>%
  filter(num > 1) %>%
  dplyr::select(c("pert_id", "cell_id"))

## keep repeated pairs only
combination_gam_nonuniq <- combination_gam %>%
  inner_join(repeated_cpd_cell_pairs, by = c("pert_id", 
                                             "cell_id"))

## separate into phase 1 and phase 2
combination_nonuniq_GSE92742 <- combination_gam_nonuniq %>%
  filter(phase == 1)
combination_nonuniq_GSE70138 <- combination_gam_nonuniq %>%
  filter(phase == 2)
