
Kpcoeff_RR <- function(logP, pKa=0, fup, BP=1, type=1, dattype=0){

  dat_tissue <- dat_RR

  Volume_pl <- 3.15

  P_W <- 10^(logP)   # octonal:water partition coeff
  logD<- 1.115*logP - 1.35 #oil:water partition coeff

  P_OW <- 10^(logD)

  # PH setting

  pH_IW <- 7       #pH of intracellular tissue water
  pH_P <- 7.4      #pH of plasma
  pH_RBC <- 7.22    #pH of blood cells

  # Plasma data

  fNLp <- 0.0023
  fNPp <- 0.0013
  
  # RBC data

  HCT <- 0.44 #hematocrit
  
  fNLb <- 0.0017
  fNPb <- 0.0029
  fIWb <- 0.603
  fAPb <- 0.5

  Kpu_bc <- (HCT - 1 + BP)/(HCT*fup)

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
              #6-zwitterion
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
              #6-zwitterion
              1 + 10^(pKa[2] - pH_P) + 10^(pH_P - pKa[1])
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
              #6-zwitterion
              1 + 10^(pKa[2] - pH_RBC) + 10^(pH_RBC - pKa[1])
  )

  dat_tissue <- dat_tissue %>%
    mutate(P = ifelse(tissue == "Adipose", P_OW, P_W))

  Ka_AP = (Kpu_bc - Z/Y*fIWb - (P_W*fNLb + (0.3*P_W + 0.7)*fNPb)/Y) * Y/fAPb/(Z - 1)
  Ka_PR <- (1/fup - 1 - (P_W*fNLp + (0.3*P_W + 0.7)*fNPp)/Y)
  
  Ka_AP <- ifelse(Ka_AP<0, 0, Ka_AP)
  Ka_PR <- ifelse(Ka_PR<0, 0, Ka_PR)
           
  # Multiply by fup to get Kp rather than Kpu
  if(type == 1) {  # Neutral
    Kpu_tissue <- dat_tissue$f_ew + X/Y*dat_tissue$f_iw + ((dat_tissue$P*dat_tissue$f_n_l + (0.3*dat_tissue$P + 0.7)*dat_tissue$f_n_pl)/Y) + Ka_PR*dat_tissue$LR

  } else if(type %in% c(2, 4)){  #acidic / diprotic acids
    
    Kpu_tissue <- dat_tissue$f_ew + X/Y*dat_tissue$f_iw + ((dat_tissue$P*dat_tissue$f_n_l + (0.3*dat_tissue$P + 0.7)*dat_tissue$f_n_pl)/Y) + Ka_PR*dat_tissue$AR  #non lipid

  }else if(type %in% c(3, 5) & max(pKa) >= 7){   # strong base with pKa > 7
    Kpu_tissue <- dat_tissue$f_ew + X/Y*dat_tissue$f_iw + ((dat_tissue$P*dat_tissue$f_n_l + (0.3*dat_tissue$P + 0.7)*dat_tissue$f_n_pl)/Y) + (Ka_AP*dat_tissue$f_a_pl*(X - 1))/Y 

  }else if(type %in% c(3, 5) & max(pKa) < 7){   # weak base
    Kpu_tissue <- dat_tissue$f_ew + X/Y*dat_tissue$f_iw + ((dat_tissue$P*dat_tissue$f_n_l + (0.3*dat_tissue$P + 0.7)*dat_tissue$f_n_pl)/Y) + Ka_PR*dat_tissue$AR  #non lipid
  }  else if(type == 6 & pKa[2] >= 7){    # Zwitterion wit pKa[2] > 7
    Kpu_tissue <- dat_tissue$f_ew + X/Y*dat_tissue$f_iw + ((dat_tissue$P*dat_tissue$f_n_l + (0.3*dat_tissue$P + 0.7)*dat_tissue$f_n_pl)/Y) + ((Ka_AP*dat_tissue$f_a_pl*10^(pKa[2] - pH_IW)) + 10^(pH_IW - pKa[1]))/Y 
  } else{    # Zwitterion wit pKa[2] <= 7
    Kpu_tissue <- dat_tissue$f_ew + X/Y*dat_tissue$f_iw + ((dat_tissue$P*dat_tissue$f_n_l + (0.3*dat_tissue$P + 0.7)*dat_tissue$f_n_pl)/Y) + Ka_PR*dat_tissue$AR  #non lipid
  }

  Kp_tissue <- as.list(Kpu_tissue * fup)

  name <- dat_tissue$tissue %>% substr(1,2) %>% tolower()
  name <- paste("Kp", name, sep="")

  names(Kp_tissue) <- name
  return(Kp_tissue)
}


