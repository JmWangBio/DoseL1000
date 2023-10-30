
## This script retrieves the complete list of HDAC inhibitors shared by DOSE-L1000 and Drug Target Commons.
## "DTC_data.csv" can be downloaded from https://drugtargetcommons.fimm.fi/

## load library
library(tidyverse)

## read data
DTC_data <- read.delim("path/to/DTC_data.csv", sep = ",")

## get all HDAC inhibitors
DTC_HDACi_data <- DTC_data %>%
  filter(grepl("^HDAC", gene_names),
         grepl("^inhibition$", standard_type, ignore.case = TRUE))

## load perturbation info
pert_info_2 <- read.delim("path/to/GSE70138_Broad_LINCS_pert_info_2017-03-06.txt.gz",
                        sep = "\t")
pert_info_1 <- read.delim("path/to/GSE92742_Broad_LINCS_pert_info.txt.gz",
                          sep = "\t")
pert_info <- rbind(pert_info_1 %>% 
                     dplyr::select(pert_id, pert_iname, inchi_key) %>% 
                     mutate(phase = 1),
                   pert_info_2 %>% 
                     dplyr::select(pert_id, pert_iname, inchi_key) %>%
                     mutate(phase = 2))

## calculate # of instances of each compound
DTC_HDACi_freq_df <- DTC_HDACi_data %>%
  inner_join(pert_info, 
             by = c("standard_inchi_key" = "inchi_key")) %>%
  group_by(standard_inchi_key) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  filter(n > 5)

## obtain a list of all HDAC inhibitors
DTC_HDACi_df <- DTC_HDACi_data %>%
  inner_join(pert_info, 
             by = c("standard_inchi_key" = "inchi_key")) %>%
  inner_join(DTC_HDACi_freq_df, by = "standard_inchi_key") %>%
  distinct(pert_id, pert_iname)
