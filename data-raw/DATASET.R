## code to prepare `DATASET` dataset goes here


dat_uni <- read.csv("data-raw/unified_tissue_comp.csv") # unified tissue composition
dat_PT <- read.csv("data-raw/tissue_comp_P&T.csv") # data reported by Poulin and Theil
dat_Berez <- read.csv("data-raw/PKSim_tissue_comp_PT_Berez.csv") # data used by PK-Sim for PT and Berez methods
dat_RR <- read.csv("data-raw/tissue_comp_R&R.csv") # data reported by Rodgers and Rowland
dat_Schmitt <- read.csv("data-raw/PKSim_tissue_comp_Schmitt.csv") # data used by PK-Sim for Schmitt method
dat_pksim <- read.csv("data-raw/PKSim_tissue_comp_pksim.csv") # data used by PK-Sim for PK-Sim method

usethis::use_data(dat_uni, dat_PT, dat_Berez, dat_RR, dat_Schmitt, dat_pksim, overwrite = TRUE, internal = TRUE)

