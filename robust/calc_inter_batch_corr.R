
## This script assumes "fit_all_models_GSE92742.R", "fit_all_models_GSE70138.R", and "calc_sd_log2fc.R" have been run. The output "all_tests_df" should be saved as "all_tests_df_%s.rds", where "%s" represents the index of the gene (1, 2, ..., 978).
## It calculates pairwise inter-batch Pearson correlation coefficients of log2 fold changes for all repeated perturbation conditions.
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

## create data frame to store result
pw_cor_by_cpd_df <- data.frame()
k <- 1

## loop through all conditions, extract two groups, calculate inter-batch correlation
for (cond in unique(all_Diff_df_complete_new$condition)) {
  tmp_all_Diff_df_complete_new <- all_Diff_df_complete_new %>%
    filter(condition == cond)
  all_tmp_groups <- unique(tmp_all_Diff_df_complete_new$group_id)
  for (tmp_model in c("gam", "rlm", "lm")) {
    for (i in 1:(length(all_tmp_groups)-1)) {
      for (j in (i+1):length(all_tmp_groups)) {
        group_1 <- all_tmp_groups[i]
        group_2 <- all_tmp_groups[j]
        tmp_Diff_complete_grp1_grp2 <- tmp_all_Diff_df_complete_new %>% 
          filter(group_id == group_1,
                 model == tmp_model) %>%
          dplyr::select(gene, Diff) %>%
          rename(Diff1 = Diff) %>%
          inner_join(tmp_all_Diff_df_complete_new %>% 
                       filter(group_id == group_2,
                              model == tmp_model) %>%
                       dplyr::select(gene, Diff) %>%
                       rename(Diff2 = Diff),
                     by = "gene")
        tmp_corP <- cor(tmp_Diff_complete_grp1_grp2$Diff1, 
                        tmp_Diff_complete_grp1_grp2$Diff2, method = "pearson")
        tmp_cor_df <- data.frame(condition = cond,
                                 model = tmp_model,
                                 grp1 = group_1,
                                 grp2 = group_2,
                                 corP = tmp_corP)
        pw_cor_by_cpd_df <- rbind(pw_cor_by_cpd_df, tmp_cor_df)
      }
    }
  }
  print(k)
  k <- k + 1
}

## calculate medians
pw_cor_by_cpd_df %>% 
  group_by(model) %>% 
  summarise(corP.median = median(corP, na.rm=TRUE))
