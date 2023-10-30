
## This script assumes "fit_all_models_GSE92742.R", "fit_all_models_GSE70138.R", and "calc_sd_log2fc.R" have been run. The output "all_tests_df" should be saved as "all_tests_df_%s.rds", where "%s" represents the index of the gene (1, 2, ..., 978).
## It calculates pairwise intra-batch Pearson correlation coefficients of log2 fold changes for perturbation conditions within the same batches.
## Make sure to update the paths to "all_tests_df_%s.rds" and "complete_pairs.rds".

## load library
library(dplyr)

## load data; keep only "complete" perturbation conditions
complete_pairs <- readRDS("path/to/complete_pairs.rds")
all_Diff_df_complete <- do.call('rbind', 
                                lapply(1:978, 
                                       function(i) {
                                         print(i)
                                         ## load results from fit_all_models_GSE92742.R
                                         all_tests_df_GSE92742 <- readRDS(sprintf("path/to/all_tests_df_%s.rds", i))
                                         ## load results from fit_all_models_GSE70138.R
                                         all_tests_df_GSE70138 <- readRDS(sprintf("path/to/all_tests_df_%s.rds", i))
                                         ## combine results
                                         all_tests_df <- rbind(all_tests_df_GSE92742, all_tests_df_GSE70138)
                                         ## keep only data from "complete" perturbation conditions
                                         all_tests_df_complete <- all_tests_df %>%
                                           inner_join(complete_pairs, 
                                                      by = c("pert_id", "cell_id",
                                                             "pert_time", "logConc"))
                                       }
                                ))

## create a new column "condition"
all_Diff_df_complete_new <- all_Diff_df_complete %>%
  mutate(condition = paste0(pert_id, "_", logConc, "_",
                            pert_time, "_", cell_id)) %>%
  dplyr::select(condition, group_id, gene, model, Diff)

## create list to store result
pw_cor_by_grp_df_lst <- list()
k <- 1

## loop through all groups, extract two cpds, calculate intra-batch correlation
for (group in unique(all_Diff_df_complete_new$group_id)) {
  tmp_all_cond_gene_pairs <- all_Diff_df_complete_new %>%
    filter(group_id == group) %>%
    dplyr::select(condition, gene) %>%
    distinct()
  tmp_all_cond_gene_combs <- tmp_all_cond_gene_pairs %>%
    group_by(gene) %>%
    do(as_data_frame(t(combn(.$condition, m = 2)))) %>%
    ungroup() 
  for (tmp_model in c("gam", "rlm", "lm")) {
    tmp_all_Diff_df_complete_new <- all_Diff_df_complete_new %>%
      filter(group_id == group, model == tmp_model) %>%
      dplyr::select(condition, gene, Diff)
    tmp_pw_cor_by_grp_df <- tmp_all_cond_gene_combs %>%
      inner_join(tmp_all_Diff_df_complete_new, 
                 by = c("gene" = "gene", "V1" = "condition")) %>% 
      inner_join(tmp_all_Diff_df_complete_new, 
                 by = c("gene" = "gene", "V2" = "condition")) %>%
      group_by(V1, V2) %>% 
      summarize(corP = cor(Diff.x, Diff.y, method = "pearson")) %>%
      mutate(model = tmp_model,
             group_id = group)
    pw_cor_by_grp_df_lst[[k]] <- tmp_pw_cor_by_grp_df
    print(k)
    k <- k + 1
  }
}

## combine list into df
total_pw_cor_by_grp_df <- do.call('rbind', pw_cor_by_grp_df_lst)

## calculate medians
total_pw_cor_by_grp_df %>% 
  group_by(model) %>% 
  summarise(corP.median = median(corP, na.rm=TRUE))
