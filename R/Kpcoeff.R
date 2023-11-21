
KpCoeff <- function(logP, pKa, fup, BP, type = 3, pred = "P&T", dat) {
  if (pred == "P&T") {
    pcoeff <- calcKp_PT(logP = logP, pKa = pKa, fup = fup, BP = BP, type = type)
  } else if (pred == "Berez") { 
    pcoeff <- calcKp_Berez(logP = logP, pKa = pKa, fup = fup, BP = BP, type = type)
  } else if (pred == "Schmitt") { 
    pcoeff <- calcKp_Schmitt(logP = logP, pKa = pKa, fup = fup, type = type)
  } else {
    pcoeff <- calcKp_RR(logP = logP, pKa = pKa, fup = fup, BP = BP, type = type)
  }
  
  return(pcoeff)
}