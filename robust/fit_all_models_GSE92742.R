
## This function fits GAMs, RLMs, and LMs and runs differential expression analysis for all repeated combinations of compounds and cell lines.
## It requires output from 'find_common_cpd_cell_pairs.R': combination_nonuniq_GSE92742.
## Note:
## 1. Change line 43 to run the analysis for different genes (i = 1, 2, ..., 978).
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

## read common cpd-cell line-project combinations
repeated_cpd_cell_proj_combs <- readRDS("path/to/combination_nonuniq_GSE92742.rds")

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
all_tests_lst <- list()
k <- 1

## loop through all perturbations
for (tmp_pert_id in all_pert_ids[1:300]) {
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
      ## check if a pair of perturbation id, cell id, and group id are repeated
      if (repeated_cpd_cell_proj_combs %>% 
          filter(pert_id == tmp_pert_id,
                 cell_id == tmp_cell_id,
                 group_id == tmp_group_id) %>%
          nrow() == 0) {
        next
      }
      
      ## extract all relevant DMSO controls (on the same plates as selected perturbations)
      tmp_DMSO_inst_ids <- inst_info %>% 
        filter(rna_plate %in% gct_df_pert$rna_plate,
               pert_id == "DMSO", 
               cell_id == tmp_cell_id) %>%
        pull(inst_id)
      
      ## combine selected perturbation and relevant DMSO
      tmp_gct_df_pert <- gct_df_pert %>%
        filter(cell_id == tmp_cell_id) %>%
        rbind(., gct_df %>%
                filter(inst_id %in% tmp_DMSO_inst_ids)) %>%
        dplyr::select(gene, inst_id, abundance, cell_id, 
                      pert_dose, pert_time, rna_plate, pert_id) %>%
        mutate(pert_dose = ifelse(pert_id == "DMSO", 0, pert_dose),
               pert_time = factor(pert_time, 
                                  levels = unique(pert_time)))
      
      ## calculate pseudo-concentration
      min_conc <- 0
      if (length(unique(tmp_gct_df_pert$pert_dose)) >= 5) {
        min_conc <- sort(unique(gct_df_pert$pert_dose), decreasing = FALSE)[1]^3 / 
          sort(unique(gct_df_pert$pert_dose), decreasing = FALSE)[2]^2
      }
      
      ## ln-transform concentration
      tmp_gct_df_pert <- tmp_gct_df_pert %>% 
        mutate(logConc = log(pert_dose + min_conc + 1 - 1))
      
      ## make a data frame covering all combinations of time and concentration of interest
      d.grid <- expand.grid(pert_time = unique(tmp_gct_df_pert$pert_time), 
                            logConc = unique(tmp_gct_df_pert$logConc)) 
      d.grid <- d.grid[order(d.grid$pert_time, 
                             d.grid$logConc), ]
      
      #######################
      #### fit GAM model ####
      #######################
      ## check whether there is more than one time point
      if (length(unique(tmp_gct_df_pert$pert_time)) > 1) {
        all.fits.gam <- tryCatch({gam(abundance ~ pert_time + s(logConc, k = 4,
                                                                by = pert_time),
                                      data = tmp_gct_df_pert)}, error = function(e) {
                                        e
                                      })
      } else {
        all.fits.gam <- tryCatch({gam(abundance ~ s(logConc, k = 4,
                                                    by = pert_time),
                                      data = tmp_gct_df_pert)}, error = function(e) {
                                        e
                                      })
      }
      
      ## differential expression analysis
      all.tests.gam <- tryCatch({
        get.gam.tests(all.fits = all.fits.gam, 
                      d.grid = d.grid)
      }, error = function(e) {
        e
      })
      
      #################
      #### fit RLM ####
      #################
      ## find all complete pairs of log conc. and time
      complete_logConc_pert_time_pairs <- tmp_gct_df_pert %>%
        dplyr::select(logConc, pert_time, abundance) %>%
        group_by(logConc, pert_time) %>%
        summarise(abundance = mean(abundance)) %>%
        pivot_wider(id_cols = "logConc", 
                    names_from = "pert_time",
                    values_from = "abundance") %>%
        na.omit() %>%
        pivot_longer(!logConc,
                     names_to = "pert_time",
                     values_to = "abundance") %>%
        dplyr::select(logConc, pert_time)
      
      ## keep data with complete pairs 
      tmp_gct_df_pert <- tmp_gct_df_pert %>%
        inner_join(., complete_logConc_pert_time_pairs,
                   by = c("logConc", "pert_time")) %>%
        mutate(logConc = factor(logConc, levels = sort(unique(logConc),
                                                       decreasing = FALSE)),
               pert_time = factor(pert_time, levels = unique(pert_time)))
      
      ## keep grids with complete pairs  
      d.grid <- d.grid %>%
        inner_join(., 
                   complete_logConc_pert_time_pairs,
                   by = c("pert_time", "logConc")) %>%
        mutate(pert_time = factor(pert_time, 
                                  levels = levels(tmp_gct_df_pert$pert_time)),
               logConc = factor(logConc, 
                                levels = levels(tmp_gct_df_pert$logConc)))
      
      ## check whether there is more than one time point
      if (length(unique(tmp_gct_df_pert$pert_time)) >= 2) {
        all.fits.rlm <- tryCatch({rlm(abundance ~ logConc * pert_time, 
                                  data = tmp_gct_df_pert)}, error = function(e) {
                                    e
                                  })
      } else {
        all.fits.rlm <- tryCatch({rlm(abundance ~ logConc, 
                                  data = tmp_gct_df_pert)}, error = function(e) {
                                    e
                                  })
      }
      
      ## differential expression analysis
      all.tests.rlm <- tryCatch({
        get.rlm.tests(all.fits = all.fits.rlm, 
                      d.grid = d.grid)
      }, error = function(e) {
        e
      })
      
      ################
      #### fit LM ####
      ################
      ## check whether there is more than one time point
      if (length(unique(tmp_gct_df_pert$pert_time)) >= 2) {
        all.fits.lm <- tryCatch({lm(abundance ~ logConc * pert_time, 
                                      data = tmp_gct_df_pert)}, error = function(e) {
                                        e
                                      })
      } else {
        all.fits.lm <- tryCatch({lm(abundance ~ logConc, 
                                      data = tmp_gct_df_pert)}, error = function(e) {
                                        e
                                      })
      }
      
      ## differential expression analysis
      all.tests.lm <- tryCatch({
        get.lm.tests(all.fits = all.fits.lm, 
                     d.grid = d.grid)
      }, error = function(e) {
        e
      })
      
      ## if error, skip
      if ((inherits(all.tests.gam, "error")) | (inherits(all.tests.rlm, "error")) | (inherits(all.tests.lm, "error"))) {
        next
      }
      
      ## add model info and combine the results
      all.tests.gam$model <- "gam"
      all.tests.rlm$model <- "rlm"
      all.tests.lm$model <- "lm"
      all.tests <- rbind(all.tests.gam, all.tests.rlm, all.tests.lm)
      
      ## annotate
      all.tests <- all.tests %>%
        remove_rownames() %>%
        mutate(pert_id = tmp_pert_id,
               gene = gene_id,
               pert_time = as.numeric(as.character(pert_time)),
               logConc = as.numeric(as.character(logConc)),
               cell_id = tmp_cell_id,
               group_id = tmp_group_id) %>%
        mutate(pert_dose = exp(logConc) - min_conc)
      
      ## store result to list
      all_tests_lst[[k]] <- all.tests
      
      k <- k + 1
      if (k %% 100 == 0) {
        print(k)
      }
    }
  }
}

## concatenate list into df
all_tests_df <- do.call('rbind', all_tests_lst)
