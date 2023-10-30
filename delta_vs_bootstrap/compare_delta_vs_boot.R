
## This script compares delta method against bootstrapping.
## Make sure to update the path to the funs folder.

## load libraries
library(dplyr)
library(mgcv)

## load functions
sapply(list.files("path/to/funs", full.names=TRUE), 
       source)

## simulate data
set.seed(1)
sim_data <- data.frame(pert_time = factor(24), logConc = seq(-3, 3, 0.5))
sim_data$abundance <- 5 + 10^sim_data$logConc / (10^sim_data$logConc + 1) + 
  rnorm(13, mean = 0, sd = 0.1)
all.fits <- gam(abundance ~ s(logConc, k = 4,
                              by = pert_time),
                data = sim_data)

####################
### delta method ###
####################
start <- Sys.time()
mfc_lx50_delta <- get.gam.stats(all.fits, B = 10001)
print(paste0("Delta method takes ", Sys.time() - start, " seconds."))

#####################
### bootstrapping ###
#####################
start <- Sys.time()
mfc_lx50_boot_1000 <- get.gam.stats.boot(all.fits, B = 10001, 
                                         bootn = 1000)
print(paste0("Bootstrapping takes ", Sys.time() - start, " seconds."))
