
## This script builds a binary classifier for HDAC inhibitors and runs leave-one-out cross-validation (LOOCV) on the classifier.
## It requires output from "find_HDACi_DTC.R", i.e., DTC_HDACi_df.
## On line 28, change cell_id to "MCF7" or "A549" to run the script for other cell lines.
## Make sure to update the paths to interaction.rds, combination.rds, test.rds, and DTC_HDACi_df.rds.

## load libraries
library(dplyr)
library(tibble)
library(fgsea)
library(yardstick)

## load data
interaction <- readRDS("path/to/interaction.rds")
combination <- readRDS("path/to/combination.rds")
test <- readRDS("path/to/test.rds")
HDACi_cpd_df <- readRDS("path/to/DTC_HDACi_df.rds")

## keep 6 hr only; extract degrees of freedom from test
degf <- test %>%
  filter(pert_time == 6) %>%
  dplyr::select(comb_index, gene, degf) %>%
  distinct()

## extract efficacy data within each cell line
interaction_sig_by_eff <- interaction %>%
  inner_join(combination, by = "comb_index") %>%
  filter(cell_id == "PC3", pert_time == 6) %>%
  inner_join(degf, by = c("comb_index", "gene")) %>%
  ## calculate p-value
  mutate(tstat = lefficacy / se_lefficacy,
         pval = 2 * pt(-abs(tstat), degf)) %>%
  ## calculate adjusted p-value
  group_by(comb_index) %>%
  mutate(adj.pval = p.adjust(pval, method = "fdr")) %>%
  ungroup() %>%
  ## order by comb_index
  arrange(comb_index) %>%
  ## create column called condition
  mutate(cond = paste0(pert_id, "_", group_id)) %>%
  dplyr::select(comb_index, cond, pert_id, group_id,
                gene, adj.pval, lefficacy)

## get conditions of all cpds
all_cond_lst <- interaction_sig_by_eff %>%
  pull(cond) %>%
  unique(.)

## get conditions of all cpds targeting HDAC
HDACi_cond_lst <- interaction_sig_by_eff %>%
  filter(pert_id %in% HDACi_cpd_df$pert_id) %>%
  pull(cond) %>%
  unique(.)

## get conditions of all cpds not targeting HDAC
nonHDACi_cond_lst <- setdiff(all_cond_lst, HDACi_cond_lst)

#################################
############ efficacy ###########
#################################
## convert to list with lefficacy
interaction_sig_by_eff_ranks_lst <- lapply(interaction_sig_by_eff %>%
                                                  group_split(comb_index), function(x) {
                                                    x %>% 
                                                      arrange(desc(lefficacy)) %>%
                                                      dplyr::select(gene, lefficacy) %>%
                                                      deframe()
                                                  })
names(interaction_sig_by_eff_ranks_lst) <- all_cond_lst

## convert to list with gene IDs only
interaction_sig_by_eff_genesets_lst <- lapply(interaction_sig_by_eff %>%
                                                     group_split(comb_index), function(x) {
                                                       x %>% 
                                                         filter(adj.pval < 0.05) %>%
                                                         arrange(desc(lefficacy)) %>%
                                                         pull(gene)
                                                     })
names(interaction_sig_by_eff_genesets_lst) <- all_cond_lst

####################################
### run fgsea (background drugs) ###
####################################
fgseaRes_bg_drug_lst_1 <- lapply(1:length(interaction_sig_by_eff_ranks_lst[nonHDACi_cond_lst]),
                                 function(i) {
                                   print(i)
                                   x <- interaction_sig_by_eff_ranks_lst[nonHDACi_cond_lst][[i]]
                                   fgseaRes <- fgsea(pathways = interaction_sig_by_eff_genesets_lst[HDACi_cond_lst],
                                                     stats = x,
                                                     minSize = 15,
                                                     maxSize = 500)
                                   fgseaRes <- fgseaRes %>% 
                                     dplyr::select(pathway, ES) %>%
                                     rename(pathway.B = pathway) %>%
                                     mutate(pathway.A = nonHDACi_cond_lst[i]) %>%
                                     relocate(pathway.A, .before = pathway.B)
                                   return(fgseaRes)
                                 })

fgseaRes_bg_drug_lst_2 <- lapply(1:length(interaction_sig_by_eff_ranks_lst[HDACi_cond_lst]),
                                 function(i) {
                                   print(i)
                                   x <- interaction_sig_by_eff_ranks_lst[HDACi_cond_lst][[i]]
                                   fgseaRes <- fgsea(pathways = interaction_sig_by_eff_genesets_lst[nonHDACi_cond_lst],
                                                     stats = x,
                                                     minSize = 15,
                                                     maxSize = 500)
                                   fgseaRes <- fgseaRes %>%
                                     dplyr::select(pathway, ES) %>%
                                     rename(pathway.A = pathway) %>%
                                     mutate(pathway.B = HDACi_cond_lst[i]) %>% 
                                     relocate(pathway.A, .before = pathway.B)
                                   return(fgseaRes)
                                 })

fgseaRes_bg_drug_df <- do.call('rbind', c(fgseaRes_bg_drug_lst_1, fgseaRes_bg_drug_lst_2))

## calculate BAES score
BAES_bg_drug_df <- fgseaRes_bg_drug_df %>%
  group_by(pathway.A, pathway.B) %>%
  summarize(BAES = sum(ES) / 2) %>%
  ungroup()

################################
### run fgsea (HDAC ligands) ###
################################
fgseaRes_HDAC_ligand_lst_1 <- lapply(1:length(interaction_sig_by_eff_ranks_lst[HDACi_cond_lst]),
                                 function(i) {
                                   print(i)
                                   x <- interaction_sig_by_eff_ranks_lst[HDACi_cond_lst][[i]]
                                   fgseaRes <- fgsea(pathways = interaction_sig_by_eff_genesets_lst[HDACi_cond_lst][-i],
                                                     stats = x,
                                                     minSize = 15,
                                                     maxSize = 500)
                                   fgseaRes <- fgseaRes %>% 
                                     dplyr::select(pathway, ES) %>%
                                     rename(pathway.B = pathway) %>%
                                     mutate(pathway.A = HDACi_cond_lst[i]) %>%
                                     relocate(pathway.A, .before = pathway.B)
                                   return(fgseaRes)
                                 })

fgseaRes_HDAC_ligand_df_1 <- do.call('rbind', fgseaRes_HDAC_ligand_lst_1)
fgseaRes_HDAC_ligand_df_2 <- fgseaRes_HDAC_ligand_df_1 %>%
  `colnames<-`(c("pathway.B", "pathway.A", "ES"))
fgseaRes_HDAC_ligand_df <- rbind(fgseaRes_HDAC_ligand_df_1,
                               fgseaRes_HDAC_ligand_df_2)

## calculate BAES score
BAES_HDAC_ligand_df <- fgseaRes_HDAC_ligand_df %>%
  group_by(pathway.A, pathway.B) %>%
  summarize(BAES = sum(ES) / 2) %>%
  ungroup()

########################
### cross-validation ###
########################
## combine pos and neg data
BAES_total_df <- rbind(BAES_bg_drug_df %>% 
                         mutate(class = "neg"),
                       BAES_HDAC_ligand_df %>%
                         mutate(class = "pos")) %>%
  mutate(class = factor(class, 
                        levels = c("pos", "neg")))
  
## leave-one-out cross-validation
cv_res_lst <- list()
k <- 1
for (cond in all_cond_lst) {
  ## separate data into training and validation
  BAES_total_df_train <- BAES_total_df %>%
    filter(pathway.A != cond)
  BAES_total_df_val <- BAES_total_df %>%
    filter(pathway.A == cond)
  ## calculate LOI
  LOI_total_df_train <- BAES_total_df_train %>%
    group_by(pathway.A) %>%
    summarize(LOI = sum(BAES) / n(),
              class = first(class)) %>%
    ungroup()
  LOI_total_df_val <- BAES_total_df_val %>%
    group_by(pathway.A) %>%
    summarize(LOI = sum(BAES) / n(),
              class = first(class)) %>%
    ungroup()
  ## fit binary classifier
  roc_curve_df <- yardstick::roc_curve(LOI_total_df_train, class, LOI)
  roc_curve_df <- roc_curve_df %>%
    mutate(bal_accuracy = (specificity + sensitivity) / 2)
  ## determine threshold
  pred_class <- "neg"
  LOI_threshold <- roc_curve_df %>%
    filter(bal_accuracy == max(bal_accuracy)) %>%
    pull(.threshold)
  ## make prediction
  if (LOI_total_df_val$LOI > LOI_threshold) {
    pred_class <- "pos"
  }
  LOI_total_df_val <- LOI_total_df_val %>%
    mutate(pred = pred_class,
           threshold = LOI_threshold)
  cv_res_lst[[k]] <- LOI_total_df_val
  k <- k + 1
  print(k)
}

cv_res_df <- do.call('rbind', cv_res_lst)
cv_res_df <- cv_res_df %>%
  mutate(pred = factor(pred, 
                       levels = c("pos", "neg")))

## sensitivity and specificity
sens(cv_res_df, truth = class, estimate = pred)
spec(cv_res_df, truth = class, estimate = pred)
