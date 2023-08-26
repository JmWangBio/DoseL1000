
## This script makes the heat map in Figure 4.
## Make sure to update the paths to interaction.rds, combination.rds, and GSE92742_Broad_LINCS_pert_info.txt.gz.

## load libraries
library(dplyr)
library(tidyr)
library(tibble)
library(latex2exp)
library(stringr)
library(ComplexHeatmap)
library(grid)

## load data
interaction <- readRDS(file = "path/to/interaction.rds")
combination <- readRDS(file = "path/to/combination.rds")

## extract phase 1 combinations
combination_phase1 <- combination %>%
  filter(phase == 1)

## read metadata
pert_info_phase1 <- read.delim(file = "path/to/GSE92742_Broad_LINCS_pert_info.txt.gz",
                               sep = "\t")
pert_info_phase1 <- pert_info_phase1 %>%
  dplyr::select(pert_id, pert_iname)

## keep HEPG2 interaction data only
interaction_hepg2 <- interaction %>%
  filter(comb_index %in% (combination_phase1 %>% 
                            filter(cell_id == "HEPG2") %>% 
                            pull(comb_index)))

## add compound names
interaction_hepg2_annot <- interaction_hepg2 %>%
  inner_join(combination_phase1, by = "comb_index") %>%
  inner_join(pert_info_phase1, by = "pert_id") %>%
  distinct()

## extract 24hr data; keep compounds with valid generic names
interaction_hepg2_24h_annot_wide <- interaction_hepg2_annot %>% 
  filter(pert_time == 24,
         grepl("^[a-z]", pert_iname)) %>%
  mutate(pert_comb = paste0(pert_iname, "_",
                            group_id)) %>%
  pivot_wider(id_cols = "gene",
              names_from = "pert_comb",
              values_from = "lefficacy") %>%
  column_to_rownames(var = "gene") %>%
  as.matrix(.)

## make heat map
col_vec <- ifelse(grepl("tamoxifen|toremifene|clomifene|raloxifene|fulvestrant", 
                        colnames(interaction_hepg2_24h_annot_wide)), 
                  "black", ifelse(grepl("apicidin|trichostatin-a|vorinostat|panobinostat|curcumin|entinostat|tubastatin-a",
                                        colnames(interaction_hepg2_24h_annot_wide)),
                                  "purple", "darkgray"))

group_vec <- as.character(sapply(colnames(interaction_hepg2_24h_annot_wide), 
                                 function(x) {
                                   str_split(x, "_")[[1]][2]
                                 }))

colnames(interaction_hepg2_24h_annot_wide) <- gsub(
  "\\.1", "", gsub("_", "-", make.names(
    gsub("-", "_", as.character(sapply(colnames(interaction_hepg2_24h_annot_wide), 
                                       function(x) {
                                         str_split(x, "_")[[1]][1]
                                       }))), unique = TRUE
  ))
)

htmp <- Heatmap(interaction_hepg2_24h_annot_wide, 
                show_row_names = FALSE,
                heatmap_legend_param = list(
                  title = expression(atop(log[2]*'['*E[max]*']' , 
                                          bgroup('(', log[2]*'['*I[max]*']', ')') )),
                  legend_direction = "vertical"),
                show_column_dend = TRUE, 
                show_row_dend = FALSE,
                column_names_gp = grid::gpar(col = col_vec,
                                             fontsize = 12),
                clustering_method_columns = "ward.D2",
                top_annotation = HeatmapAnnotation(batch = group_vec,
                                                   col = list(batch = c("PCLB001" = "lightblue",
                                                                        "PCLB002" = "lavender",
                                                                        "PCLB003" = "lightgreen"))))
draw(htmp, heatmap_legend_side = "right", 
     annotation_legend_side = "right",
     merge_legend = TRUE)
