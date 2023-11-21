
pcoeffs <- function(logP, pKa, fup, BP, type = 3, pred = "P&T", dat) {
  if (pred == "P&T") {
    source("Kp_calculator/code/CalcKp_P&T.R")
    pcoeff <- calcKp_PT(logP = logP, pKa = pKa, fup = fup, BP = BP, type = type)
  } else if (pred == "Berez") { # Berezhkovskiy
    source("Kp_calculator/code/CalcKp_Berez.R")
    pcoeff <- calcKp_Berez(logP = logP, pKa = pKa, fup = fup, BP = BP, type = type)
  } else if (pred == "Schmitt") { # Schmitt, Walter 2008
    source("Kp_calculator/code/CalcKp_Schmitt.R")
    pcoeff <- calcKp_Schmitt(logP = logP, pKa = pKa, fup = fup, type = type)
  } else {
    source("Kp_calculator/code/CalcKp_R&R.R") # Rodgers and Rowland 2006
    pcoeff <- calcKp_RR(logP = logP, pKa = pKa, fup = fup, BP = BP, type = type)
  }
  
  return(pcoeff)
}