
## This script calculates the Pearson correlation between the results of fixed models and the results of mixed models.
## It assumes "fit_mixed_model_GSE92742.R" and "fit_mixed_model_GSE70138.R" have been run. The output "all_tests_df" should be saved as "all_tests_df_%s.rds", where "%s" represents the index of the gene (1, 2, ..., 978).
## Make sure to update the paths to test.rds, combination.rds, and all_tests_df_%s.rds.

## load library
library(dplyr)

## read data
test <- readRDS(file = "path/to/test.rds")
combination <- readRDS(file = "path/to/combination.rds")

## read results from mixed models
all_log2fc_mixed_GSE92742_total <- do.call('rbind', 
                                         lapply(1:978,
                                                function(i) {
                                                  print(i)
                                                  all_log2fc_df <- readRDS(sprintf("path/to/all_tests_df_%s.rds", i))
                                                }
                                         ))
all_log2fc_mixed_GSE92742_total$phase <- 1
all_log2fc_mixed_GSE92742_total <- all_log2fc_mixed_GSE92742_total %>%
  rename(det_plate = rna_plate)

all_log2fc_mixed_GSE70138_total <- do.call('rbind', 
                                         lapply(1:978,
                                                function(i) {
                                                  print(i)
                                                  all_log2fc_df <- readRDS(sprintf("path/to/all_tests_df_%s.rds", i))
                                                }
                                         ))
all_log2fc_mixed_GSE70138_total$phase <- 2

## combine phase1 and phase2 results
all_log2fc_mixed_total <- rbind(all_log2fc_mixed_GSE92742_total, 
                                all_log2fc_mixed_GSE70138_total)

## create data frame combining results from fixed and mixed models
all_log2fc_mixed_total_refmt <- all_log2fc_mixed_total %>%
  dplyr::select(pert_time, pert_dose, gene, Diff,
                pert_id, group_id, cell_id, phase, model) %>%
  rename(Diff2 = Diff) %>%
  inner_join(combination, ., by = c("pert_id", "group_id", 
                                    "cell_id", "phase")) %>%
  dplyr::select(comb_index, pert_time, pert_dose, gene, Diff2, model)

all_log2fc_fixed_mixed_total <- test %>%
  dplyr::select(comb_index, pert_dose, pert_time,
                gene, Diff) %>%
  inner_join(all_log2fc_mixed_total_refmt,
             by = c("comb_index", "pert_time", 
                    "pert_dose", "gene"))

## calculate correlation within each perturbed condition
all_log2fc_fixed_mixed_corr_df <- all_log2fc_fixed_mixed_total %>%
  group_by(comb_index, pert_dose, pert_time, model) %>%
  dplyr::summarize(corP = cor(Diff, Diff2, method = "pearson"))
