
## This script assumes "fit_all_models_GSE92742.R" and "fit_all_models_GSE70138.R" have been run. The output "all_tests_df" should be saved as "all_tests_df_%s.rds", where "%s" represents the index of the gene (1, 2, ..., 978).
## It fulfills two objectives:
## 1. calculate standard deviation of log2 fold changes for all repeated perturbation conditions
## 2. identify all "complete" perturbation conditions with results from all three models for all genes; this is required for calculating inter- and intra-batch correlation.
## Make sure to update the paths to "all_tests_df_%s.rds" and "complete_pairs.rds".

## load library
library(dplyr)

## load data; calculate standard deviation on the fly
all_Diff_sd_df <- do.call('rbind', 
                          lapply(1:978, 
                                 function(i) {
                                   print(i)
                                   ## load results from fit_all_models_GSE92742.R
                                   all_tests_df_GSE92742 <- readRDS(sprintf("path/to/all_tests_df_%s.rds", i))
                                   ## load results from fit_all_models_GSE70138.R
                                   all_tests_df_GSE70138 <- readRDS(sprintf("path/to/all_tests_df_%s.rds", i))
                                   ## combine results
                                   all_tests_df <- rbind(all_tests_df_GSE92742, all_tests_df_GSE70138)
                                   ## calculate sd of Diff within each combination of compound, cell line, time, logConc, and model
                                   all_tests_df_stats <- all_tests_df %>%
                                     dplyr::select(pert_id, cell_id, pert_time, logConc, model, Diff) %>%
                                     group_by(pert_id, cell_id, pert_time, logConc, model) %>%
                                     summarize(Diff.sd = sd(Diff)) %>%
                                     ungroup() %>%
                                     mutate(gene = unique(all_tests_df$gene))
                                   return(all_tests_df_stats)
                                 }
                          ))

## remove NA values
all_Diff_sd_df_nona <- all_Diff_sd_df %>% filter(!is.na(Diff.sd))

## find "complete" perturbation conditions with results from all three models for all genes
complete_pairs <- all_Diff_sd_df_nona %>% 
  group_by(pert_id,
           cell_id,
           pert_time,
           logConc) %>%
  summarise(n = n()) %>%
  filter(n == 3 * 978) %>%
  dplyr::select(pert_id, cell_id, pert_time, logConc)

## keep data from "complete" perturbation conditions; ready for violin plot
all_Diff_sd_df_complete <- complete_pairs %>%
  inner_join(all_Diff_sd_df_nona, ., 
             by = c("pert_id", "cell_id",
                    "pert_time", "logConc"))

## calculate medians
all_Diff_sd_df_complete %>%
  group_by(model) %>%
  summarise(Diff.sd.median = median(Diff.sd))

## save complete pairs (required for calculating inter- and intra-batch correlation)
saveRDS(complete_pairs, file = "path/to/complete_pairs.rds")
