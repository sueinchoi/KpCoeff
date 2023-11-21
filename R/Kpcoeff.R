
KpCoeff <- function(logP, pKa, fup, BP, type = 3, pred = "P&T", dat) {
  if (pred == "P&T") {
    pcoeff <- Kpcoeff_PT(logP = logP, pKa = pKa, fup = fup, BP = BP, type = type)
  } else if (pred == "Berez") { 
    pcoeff <- Kpcoeff_Berez(logP = logP, pKa = pKa, fup = fup, BP = BP, type = type)
  } else if (pred == "Schmitt") { 
    pcoeff <- Kpcoeff_Schmitt(logP = logP, pKa = pKa, fup = fup, type = type)
  } else {
    pcoeff <- Kpcoeff_RR(logP = logP, pKa = pKa, fup = fup, BP = BP, type = type)
  }
  
  return(pcoeff)
}