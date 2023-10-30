
## This function calculates the Pearson correlation coefficients from the Q-Q plots for all combinations of compounds, cell lines, and time.
## Note:
## 1. Change line 42 to run the analysis for different genes (i = 1, 2, ..., 978).
## 2. We need the parse_gctx() function in the cmapR package to read gctx files. One could either install the cmapR package or download the source code from GitHub and load only relevant functions, as shown in this script.
## 3. Make sure to update the paths to the cmapR package, funs folder, GSE92742_Broad_LINCS_gene_info.txt.gz, GSE92742_Broad_LINCS_inst_info.txt.gz,
## GSE92742_Broad_LINCS_pert_info.txt.gz, and GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx. GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx.gz needs to be unzipped.

## load libraries
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(mgcv)
library(MASS)
library(rhdf5)

## Load custom functions
sapply(c("path/to/cmapR/R/GCT.R",
         "path/to/cmapR/R/io.R",
         "path/to/cmapR/R/utils.R"), 
       source)

sapply(list.files("path/to/funs", full.names=TRUE), 
       source)

## load metadata
gene_info <- read.delim("path/to/GSE92742_Broad_LINCS_gene_info.txt.gz",
                        sep = "\t")
inst_info <- read.delim("path/to/GSE92742_Broad_LINCS_inst_info.txt.gz",
                        sep = "\t")
pert_info <- read.delim("path/to/GSE92742_Broad_LINCS_pert_info.txt.gz",
                        sep = "\t")

## retrieve all landmark genes
lm_rid_lst <- gene_info %>% 
  filter(pr_is_lm == 1) %>%
  pull(pr_gene_id) %>%
  as.character(.)

## specify index (which landmark gene?)
i <- 1
cat(paste("Starting gene #", i, "\n"))

## read data
gct_lv3 <- parse_gctx(fname="path/to/GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx",
                      rid = lm_rid_lst[i])

## build gct data with perturbation information integrated; convert concentration units to um
gct_df <- gct_lv3@mat %>%
  as.data.frame(.) %>%
  rownames_to_column(var="gene") %>%
  pivot_longer(!gene, names_to = "inst_id", values_to = "abundance") %>%
  inner_join(., inst_info, by = "inst_id") %>%
  mutate(group_id = as.character(sapply(inst_id, function(x) str_split(x, "_")[[1]][1]))) %>%
  filter((pert_dose_unit %in% c("mm", "um") & pert_id != "DMSO") | pert_id == "DMSO") %>%
  mutate(pert_dose = ifelse(pert_dose_unit == "mm", pert_dose * 1000, pert_dose)) %>%
  mutate(pert_dose_unit = ifelse(pert_dose_unit == "mm", "um", pert_dose_unit))

## obtain unique gene ID
gene_id <- unique(gct_df$gene)

## retrieve list of all perturbation IDs
all_pert_ids <- pert_info %>% 
  filter(pert_type == "trt_cp") %>% 
  pull(pert_id) %>%
  unique(.)

## create list to store result
all_corrs_lst <- list()
k <- 1

## loop through all perturbations
for (tmp_pert_id in all_pert_ids) {
  print(tmp_pert_id)
  
  ## keep selected perturbation only
  gct_df_pert_half <- gct_df %>% 
    filter(pert_id == tmp_pert_id)
  
  ## retrieve all unique perturbagen group IDs (batch IDs)
  all_tmp_group_ids <- gct_df_pert_half %>%
    pull(group_id) %>%
    unique(.)
  
  ## loop through all group IDs (batch IDs)
  for (tmp_group_id in all_tmp_group_ids) {
    gct_df_pert <- gct_df_pert_half %>%
      filter(group_id == tmp_group_id)
    
    ## retrieve all unique cell lines
    all_tmp_cell_ids <- gct_df_pert %>% 
      pull(cell_id) %>%
      unique(.)
    
    ## loop through all cell lines
    for (tmp_cell_id in all_tmp_cell_ids) {
      ## extract all relevant DMSO controls (on the same plates as selected perturbations)
      tmp_DMSO_inst_ids <- inst_info %>% 
        filter(rna_plate %in% gct_df_pert$rna_plate,
               pert_id == "DMSO", 
               cell_id == tmp_cell_id) %>%
        pull(inst_id)
      
      ## select relevant DMSO data
      tmp_gct_df_DMSO <- gct_df %>%
        filter(inst_id %in% tmp_DMSO_inst_ids)
      
      ## loo through all times
      for (tmp_time in unique(tmp_gct_df_DMSO$pert_time)) {
        ## extract DMSO controls belonging to a time point
        tmp_gct_df_DMSO_time <- tmp_gct_df_DMSO %>%
          filter(pert_time == tmp_time)
        tmp_qqp <- qqnorm(tmp_gct_df_DMSO_time$abundance, plot.it = FALSE)
        tmp_qq_corP <- cor(tmp_qqp$x, tmp_qqp$y, method = "pearson")
        all.corrs <- data.frame(pert_id = tmp_pert_id,
                                gene = gene_id,
                                pert_time = as.numeric(as.character(tmp_time)),
                                cell_id = tmp_cell_id,
                                group_id = tmp_group_id,
                                corP = tmp_qq_corP)
        ## store result to list
        all_corrs_lst[[k]] <- all.corrs
        
        k <- k + 1
        if (k %% 100 == 0) {
          print(k)
        }
      }
    }
  }
}

## concatenate list into df
all_corrs_df <- do.call('rbind', all_corrs_lst)
