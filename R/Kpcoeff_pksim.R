
#' Prediction of tissue-prediction partition coefficient using PK-sim method
#'
#' @param logP Partition coefficient
#' @param fup Plasma unbound fraction
#' @param dattype Human physiology dataset
#' @importFrom rlang .data
#' @return A list of tissue partition coefficient in each organ
#' @export
#'
#' @examples
#' Kpcoeff_pksim(2.7,0.9,0)
#'
Kpcoeff_pksim <- function(logP, fup, dattype=0){
  if(dattype == 0){
    dat <- dat_pksim
  } else {
    dat <- dat_uni
  }
  #logMA is the log of membrane affinity = phosphatidylcholin:water (neutral phospholipid:water) partition coefficient;
  #we can use the available measurement of lipophilicity instead (logP or logD); from Schmitt, Walter (2008)

  dat_all <- dat[which(dat$tissue != "Plasma" & dat$tissue != "RBCs"), , drop = FALSE]



  logMA <- logP  #in case we don't have a direct logMA
  K_n_pl <- 10^logMA    #neutral phospholipids:water partition coefficient
  #K_protein <- 0.163 + 0.0221*K_n_pl    #protein:water partition; Schmitt, Walter (2008)
  K_protein <- ((0.81 + 0.11 * K_n_pl)/24.92)*5 # From PK-Sim (very similar value to the other method)


  kp <- (dat_all$f_water + (K_n_pl*dat_all$f_lipids) + (K_protein*dat_all$f_proteins))*fup

  #denom <- 0.945 + (10^logMA*0.00575) + (0.93*fup)  #plasma fractions
  #kp <- kp/denom  #according to Willmann et al. (2005)
  dat2 <- data.frame(tissue=dat_all$tissue, Kp=kp)
  name <- dat2$tissue %>% substr(1,2) %>% tolower()
  name <- paste("Kp", name, sep="")
  Kp <- as.list(dat2$Kp)
  names(Kp) <- name

  return(Kp)
}
