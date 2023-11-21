library(devtools)
remotes::install_github("sueinchoi/KpCoeff")
library(KpCoeff)
library(tidyverse)
logP <- 2
pKa <- 6
fup <- 0.1
BP <- 1
type <- 'neutral'
type_num <- switch(type,
                   'neutral' = 1, 
                   'acid' = 2,
                   'base' = 3,
                   'diprotic acid' = 4,
                   'diprotic base' = 5,
                   'zwitterion(ampholyte)' = 6
)

predname <- 'Poulin & Theil'
pred <- switch(predname, 
               'Poulin & Theil' = 'P&T',
               'Rowland & Rodgers' = 'R&R',
               'Berezhesky' = 'Berez',
               'PK-sim' = 'pksim',
               'Schmitt' = 'Schmitt'
)


Kpresult <- Kpcoeff(logP, pKa, fup, BP, type_num, pred) %>%
    map_df(~as.data.frame(.x), .id = 'names') %>%
    rename('value' = 2)

Kpcoeff_round <- Kpresult %>%
    mutate_at(2, round, 4)
Kpcoeff_round

