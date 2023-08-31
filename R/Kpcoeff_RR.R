
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
#' 


library(tidyverse)

Kpcoeff_RR <- function(logP, logD, pKa=0, fup, BP=1, type=1, dattype=0){
  dat_tissue <- read_csv('data-raw/Tissue_comp_R&R.csv')
  
  P_W <- 10^(logP)   # octonal:water partition coeff
  # logP_OW <- 1.115*logP - 1.35 #oil:water partition coeff
  # P_OW <- 10^(logP_OW)
  P_OW <- 10^(logD)
  Ka <- 10^(-pKa)
 

  # PH setting

  pH_IW <- 7       #pH of intracellular tissue water
  pH_P <- 7.4      #pH of plasma
  pH_RBC <- 7.22    #pH of blood cells

  # Plasma data

  fNLp <- 0.0023
  fNPp <- 0.0013
  
  # RBC data

  HCT <- 0.43 #hematocrit
  
  fNLb <- 0.0017
  fNPb <- 0.0029
  fIWb <- 0.603
  fAPb <- 0.5

  Kpu_bc <- (HCT - 1 + BP)/(HCT*fup)

  # Dissociation constant with tissue and AP
  # Ka_AP <- (Kpu_bc - Z/Y*fIWb - (P*fNLb + (0.3*P + 0.7)*fNPb)/Y) * Y/fAPb/(Z - 1)


  X <- switch(type,
              #1-neutral
              1,
              #2-monoprotic acid
              1 + 10^(pH_IW - pKa),
              #3-monoprotic base
              1 + 10^(pKa - pH_IW),
              #4-diprotic acid
              1 + 10^(pH_IW - min(pKa[1], pKa[2])) + 10^(2*pH_IW - pKa[1] - pKa[2]),
              #5-diprotic base
              1 + 10^(max(pKa[1], pKa[2]) - pH_IW) + 10^(pKa[1] + pKa[2] - 2*pH_IW),
              #6-monoprotic acid monoprotic base (acid comes first, pKa[2] = base pKa > 7)
              1 + 10^(pKa[2] - pH_IW) + 10^(pH_IW - pKa[1]),
              #7-all other zwitterions
              1 + 10^(pKa[2] - pH_IW) + 10^(pH_IW - pKa[1])
  )

  Y <- switch(type,
              #1-neutral
              1,
              #2-monoprotic acid
              1 + 10^(pH_P - pKa),
              #3-monoprotic base
              1 + 10^(pKa - pH_P),
              #4-diprotic acid
              1 + 10^(pH_P - min(pKa[1], pKa[2])) + 10^(2*pH_P - pKa[1] - pKa[2]),
              #5-diprotic base
              1 + 10^(max(pKa[1], pKa[2]) - pH_P) + 10^(pKa[1] + pKa[2] - 2*pH_P),
              #6-monoprotic acid monoprotic base (acid comes first)
              1 + 10^(pKa[2] - pH_P) + 10^(pH_P - pKa[1]),
              #7-all other zwitterions
              1 + 10^(pKa[2] - pH_P)+10^(pH_P - pKa[1])
  )
  Z <- switch(type,
              #1-neutral
              1,
              #2-monoprotic acid
              1,
              #3-monoprotic base
              1 + 10^(pKa - pH_RBC),
              #4-diprotic acid
              1,
              #5-diprotic base
              1 + 10^(max(pKa[1], pKa[2])-pH_RBC) + 10^(pKa[1] + pKa[2] - 2*pH_RBC),
              #6-monoprotic acid monoprotic base (acid comes first = pKa[1])
              1 + 10^(pKa[2] - pH_RBC) + 10^(pH_RBC - pKa[1]),
              #7-all other zwitterions
              1
  )



  dat_tissue <- dat_tissue %>%
    mutate(P = ifelse(tissue == "Adipose", P_OW, P_W))
    
  dat_tissue <- dat_tissue %>%
    mutate(
    )

  Ka_AP = (Kpu_bc - Z/Y*fIWb - (P_W*fNLb + (0.3*P_W + 0.7)*fNPb)/Y) * Y/fAPb/(Z - 1)
  Ka_PR <- (1/fup - 1 - (P_W*fNLp + (0.3*P_W + 0.7)*fNPp)/Y)

           
  # Multiply by fup to get Kp rather than Kpu
  if(type == 1) {  # Neutral
    Kpu_tissue <- dat_tissue$f_ew + X/Y*dat_tissue$f_iw + ((dat_tissue$P*dat_tissue$f_n_l + (0.3*dat_tissue$P + 0.7)*dat_tissue$f_n_pl)/Y) + Ka_PR*dat_tissue$LR

  } else if(type %in% c(2, 4)){  #acidic / diprotic acids
    
    Kpu_tissue <- dat_tissue$f_ew + X/Y*dat_tissue$f_iw + ((dat_tissue$P*dat_tissue$f_n_l + (0.3*dat_tissue$P + 0.7)*dat_tissue$f_n_pl)/Y) + Ka_PR*dat_tissue$AR  #non lipid

  }else if(type %in% c(3, 5) & max(pKa) > 7){   # strong base with pKa > 7
    Kpu_tissue <- dat_tissue$f_ew + X/Y*dat_tissue$f_iw + (dat_tissue$P*dat_tissue$f_n_l + (0.3*dat_tissue$P + 0.7)*dat_tissue$f_n_pl)/Y + (dat_tissue$Ka_AP*dat_tissue$f_a_pl*(X - 1))/Y 

  }else if(type %in% c(3, 5) & max(pKa) <= 7){   # weak base
    Kpu_tissue <- dat_tissue$f_ew + X/Y*dat_tissue$f_iw + ((dat_tissue$P*dat_tissue$f_n_l + (0.3*P + 0.7)*dat_tissue$f_n_pl)/Y) + Ka_PR*dat_tissue$AR  #non lipid
  }  else if(type == 6){    # Zwitterion wit pKa[2] > 7
    Kpu_tissue <- dat_tissue$f_ew + X/Y*dat_tissue$f_iw + ((dat_tissue$P*dat_tissue$f_n_l + (0.3*dat_tissue$P + 0.7)*dat_tissue$f_n_pl)/Y) + ((dat_tissue$Ka_AP*dat_tissue$f_a_pl*10^(pKa[2] - pH_IW)) + 10^(pH_IW - pKa[1]))/Y 
  } else{    # Zwitterion wit pKa[2] <= 7
    Kpu_tissue <- dat_tissue$f_ew + X/Y*dat_tissue$f_iw + ((dat_tissue$P*dat_tissue$f_n_l + (0.3*dat_tissue$P + 0.7)*dat_tissue$f_n_pl)/Y) + Ka_PR*dat_tissue$AR  #non lipid
  }


  # Kp <- c(Kp_ad, Kp_all)
  # Kp_rest <- mean(Kp)
  # Kp <- c(Kp, Kp_rest)
  name <- dat_tissue$tissue %>% substr(1,2) %>% tolower()
  name <- paste("Kp", name, sep="")

  result <- data.frame(name = name, Kp = Kpu_tissue)
  return(result)
  # nms_all <- dat_all$tissue %>% substr(1,2) %>% tolower()
  # nms_all <- paste("Kp", nms_all, sep="")
  # nms <- c("Kpad",nms_all, "Kprest")
  # Kp <- as.list(c(Kp_ad,Kp_all, mean(Kp_all, Kp_ad)))
  # names(Kp) <- nms

  # vols <- select(dat, FV) %>% pull()
  # prod <- vols[1:11]*Kp[1:11]

  # FV_rest <- 0.00726
  # Vss <- sum(prod[1:11]) + 0.0347*(BP - (1-0.45))/0.45 + vols[13] + FV_rest*Kp[12]

  # return(Vss)
}

Kpcoeff_RR(logP = 2.59, logD = 1.54, pKa = 9.4, fu = 0.53, BP = 2.06, type = 3)
