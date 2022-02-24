
#' Prediction of tissue-prediction partition coefficient using R&R method
#'
#' @param logP Partition coefficient
#' @param pKa Acid dissociation constant (for ampholyte or zwitterion, enter Acid pKa first)
#' @param fup Plasma unbound fraction
#' @param BP Blood-plasma ratio
#' @param type Neutral(1)/acid(2)/base(3)/diprotic acid(4)/diprotic base(5)/zwitterion(ampholyte)(6)
#' @param dattype Human physiology dataset(0 - original dataset, 1 - unified dataset)
#' @importFrom rlang .data
#'
#' @return A list of tissue partition coefficient in each organ
#' @export
#'
#' @examples
#' Kpcoeff_RR(2.7, 6, 0.9, 1, 1, 0)
#'
Kpcoeff_RR <- function(logP, pKa=0, fup, BP=1, type=1, dattype=0){
  if(dattype == 0){
    dat <- dat_RR
  } else {
    dat <- dat_uni
  }

  dat_all <- dat %>% dplyr::filter(!.data$tissue %in% c("RBCs", "Adipose", "Plasma"))  #df for all tissues except for adipose, RBCs, and plasma
  dat_ad <- dat %>% dplyr::filter(.data$tissue == "Adipose")  #df for adipose
  dat_rbc <- dat %>% dplyr::filter(.data$tissue == "RBCs") #df for RBCs
  dat_plas <- dat %>% dplyr::filter(.data$tissue == "Plasma") #df for aplasma
  pH_IW <- 7       #pH of intracellular tissue water
  pH_P <- 7.4      #pH of plasma
  pH_RBC <- 7.22    #pH of blood cells
  P <- 10^(logP)   # octonal:water partition coeff
  logP_OW <- 1.115*logP - 1.35 #oil:water partition coeff
  Ka <- 10^(-pKa)
  HCT <- 0.45 #hematocrit

  #Calculate Kp values
  Kpu_bc <- (HCT - 1 + BP)/(HCT*fup)

  X <- switch(type,
              #1-neutral
              0,
              #2-monoprotic acid
              10^(pH_IW-pKa),
              #3-monoprotic base
              10^(pKa-pH_IW),
              #4-diprotic acid
              10^(pH_IW-pKa[1])+10^(2*pH_IW-pKa[1]-pKa[2]),
              #5-diprotic base
              10^(pKa[2]-pH_IW)+10^(pKa[1]+pKa[2]-2*pH_IW),
              #6-monoprotic acid monoprotic base (acid comes first)
              10^(pH_IW-pKa[1])+10^(pKa[2]-pH_IW),
              #7-triprotic acid
              10^(pH_IW-pKa[1])+10^(2*pH_IW-pKa[1]-pKa[2])+10^(3*pH_IW-pKa[1]-pKa[2]-pKa[3]),
              #8-triprotic base
              10^(pKa[3]-pH_IW)+10^(pKa[3]+pKa[2]-2*pH_IW)+10^(pKa[1]+pKa[2]+pKa[3]-3*pH_IW),
              #9-diprotic acid monoprotic base (first two are acid)
              10^(pKa[3]-pH_IW)+10^(pH_IW-pKa[1])+10^(2*pH_IW-pKa[1]-pKa[2]),
              #10-diprotic base monoprotic acid (first one is acid)
              10^(pH_IW-pKa[1])+10^(pKa[3]-pH_IW)+10^(pKa[2]+pKa[3]-2*pH_IW))

  Y <- switch(type,
              #1-neutral
              0,
              #2-monoprotic acid
              10^(pH_P-pKa),
              #3-monoprotic base
              10^(pKa-pH_P),
              #4-diprotic acid - Which one comes pKa1???
              10^(pH_P-pKa[1])+10^(2*pH_P-pKa[1]-pKa[2]),
              #5-diprotic base - Which one come pKa2?? bigger? smallar?
              10^(pKa[2]-pH_P)+10^(pKa[1]+pKa[2]-2*pH_P),
              #6-monoprotic acid monoprotic base (acid comes first)
              10^(pH_P-pKa[1]) + 10^(pKa[2]-pH_P),
              #7-triprotic acid
              10^(pH_P-pKa[1])+10^(2*pH_P-pKa[1]-pKa[2])+10^(3*pH_P-pKa[1]-pKa[2]-pKa[3]),
              #8-triprotic base
              10^(pKa[3]-pH_P)+10^(pKa[3]+pKa[2]-2*pH_P)+10^(pKa[1]+pKa[2]+pKa[3]-3*pH_P),
              #9-diprotic acid monoprotic base (first two are acid)
              10^(pKa[3]-pH_P)+10^(pH_P-pKa[1])+10^(2*pH_P-pKa[1]-pKa[2]),
              #10-diprotic base monoprotic acid (first one is acid)
              10^(pH_P-pKa[1])+10^(pKa[3]-pH_P)+10^(pKa[2]+pKa[3]-2*pH_P))

  Z <- switch(type,
              #1-neutral
              0,
              #2-monoprotic acid
              10^(pH_RBC - pKa),
              #3-monoprotic base
              10^(pKa-pH_RBC),
              #4-diprotic acid
              10^(pH_RBC - pKa[1]) + 10^(2*pH_RBC - pKa[1] - pKa[2]),
              #5-diprotic base
              10^(pKa[2]-pH_RBC)+10^(pKa[1]+pKa[2]-2*pH_RBC),
              #6-monoprotic acid monoprotic base (acid comes first)
              10^(pH_RBC-pKa[1]) + 10^(pKa[2]-pH_RBC),
              #7-triprotic acid
              10^(pH_RBC-pKa[1])+10^(2*pH_RBC-pKa[1]-pKa[2])+10^(3*pH_RBC-pKa[1]-pKa[2]-pKa[3]),
              #8-triprotic base
              10^(pKa[3]-pH_RBC)+10^(pKa[3]+pKa[2]-2*pH_RBC)+10^(pKa[1]+pKa[2]+pKa[3]-3*pH_RBC),
              #9-diprotic acid monoprotic base (first two are acid)
              10^(pKa[3]-pH_RBC)+10^(pH_RBC-pKa[1])+10^(2*pH_RBC-pKa[1]-pKa[2]),
              #10-diprotic base monoprotic acid (first one is acid)
              10^(pH_RBC-pKa[1])+10^(pKa[3]-pH_RBC)+10^(pKa[2]+pKa[3]-2*pH_RBC))
  logD_OW <- logP_OW + log10(1/(Y+1)) # Parition coefficient to distribution coefficient (for SIMCYP only - original reference used partition coefficient instead)
  P_OW <- 10^(logD_OW)

  Ka_PR <- (1/fup - 1 - (P*dat_plas$f_n_l + (0.3*P + 0.7)*dat_plas$f_n_pl)/(1+Y))
  Ka_AP <- (Kpu_bc - (1 + Z)/(1 + Y)*dat_rbc$f_iw - (P*dat_rbc$f_n_l + (0.3*P + 0.7)*dat_rbc$f_n_pl)/(1 + Y)) * (1 + Y)/(dat_rbc$f_a_pl*Z)


  # Assign the moderate to strong bases type_calc=1 and everything else type_calc=2
  type_calc <- ifelse((type==3 & pKa[1]>7) | (type==5 & pKa[1] >7) | (type==8 & pKa[1] > 7) | (type==9 & pKa[3]>7) | (type==10 & pKa[2]>7), 1,2)

  # Re-assign the neutrals type_calc=3
  if(type==1){type_calc=3}  #neutrals
  if((type==6 & pKa[2] > 7)){type_calc=4}

  if(type_calc==1){  #moderate to strong bases
    Kp_all <- (dat_all$f_ew + ((1 + X)/(1 + Y))*dat_all$f_iw + (P*dat_all$f_n_l + (0.3*P + 0.7)*dat_all$f_n_pl)/(1 + Y) + (Ka_AP*dat_all$f_a_pl*X)/(1 + Y))*fup  #non lipid
    Kp_ad <- (dat_ad$f_ew + ((1 + X) / (1 + Y)) * dat_ad$f_iw + (P_OW * dat_ad$f_n_l + (0.3 * P_OW + 0.7) * dat_ad$f_n_pl) / (1 + Y) + (Ka_AP * dat_ad$f_a_pl * X) / (1 + Y)) * fup # lipid
  }else if(type_calc==2){   #acidic and zwitterions
    Kp_all <- (dat_all$f_ew + ((1 + X) / (1 + Y)) * dat_all$f_iw + (P * dat_all$f_n_l + (0.3 * P + 0.7) * dat_all$f_n_pl) / (1 + Y) + Ka_PR * dat_all$AR) * fup # non lipid
    Kp_ad <- (dat_ad$f_ew + ((1 + X) / (1 + Y)) * dat_ad$f_iw + (P_OW * dat_ad$f_n_l + (0.3 * P_OW + 0.7) * dat_ad$f_n_pl) / (1 + Y) + Ka_PR * dat_ad$AR) * fup # lipid
  }else if(type_calc==3){  #neutrals
    Kp_all <- (dat_all$f_ew + ((1 + X) / (1 + Y)) * dat_all$f_iw + (P * dat_all$f_n_l + (0.3 * P + 0.7) * dat_all$f_n_pl) / (1 + Y) + Ka_PR * dat_all$LR) * fup # non lipid
    Kp_ad <- (dat_ad$f_ew + ((1 + X) / (1 + Y)) * dat_ad$f_iw + (P_OW * dat_ad$f_n_l + (0.3 * P_OW + 0.7) * dat_ad$f_n_pl) / (1 + Y) + Ka_PR * dat_ad$LR) * fup # lipid
  }else {
    Kp_all <- (dat_all$f_ew + ((1 + X) / (1 + Y)) * dat_all$f_iw + (P * dat_all$f_n_l + (0.3 * P + 0.7) * dat_all$f_n_pl) / (1 + Y) + (Ka_AP * dat_all$f_a_pl * 10^(pKa[2] - pH_IW) + 10^(pH_IW - pKa[1])) / (1 + Y)) * fup # non lipid
    Kp_ad <- (dat_ad$f_ew + ((1 + X) / (1 + Y)) * dat_ad$f_iw + (P_OW * dat_ad$f_n_l + (0.3 * P_OW + 0.7) * dat_ad$f_n_pl) / (1 + Y) + (Ka_AP * dat_ad$f_a_pl * 10^(pKa[2] - pH_IW) + 10^(pH_IW - pKa[1])) / (1 + Y)) * fup # lipid
  }


  nms_all <- dat_all$tissue %>% substr(1,2) %>% tolower()
  nms_all <- paste("Kp", nms_all, sep="")
  nms <- c("Kpad",nms_all)
  Kp <- as.list(c(Kp_ad,Kp_all))
  names(Kp) <- nms


  return(Kp)
}

