
#' Calculation of tissue-plasma partition coefficient
#'
#' It calculates tissue-plasma partition coefficient using the specified prediction method and the human physiology dataset
#'
#' @param logP partition coefficient (octanol-to-water)
#' @param pKa Acid dissociation constant (for ampholyte or zwitterion, enter Acid pKa first)
#' @param fup Plasma unbound fraction
#' @param BP Blood-to-plasma ratio
#' @param type Neutral(1)/acid(2)/base(3)/diprotic acid(4)/diprotic base(5)/zwitterion(ampholyte)(6)
#' @param pred Prediction method (Poulin&Theil - "P&T", Rowland&Rowdgers - "R&R", Berezhesky - "Berez", PKsim - "pksim", Schmitt - "Schmitt")
#' @param dattype Human physiology dataset(0 - original dataset, 1 - unified dataset)
#'
#' @return A list of partition coefficient in each organ
#' @export
#'
#' @examples
#' Kpcoeff(2.6, 2, 0.9, 1.5, 2, "P&T", 0)
#' Kpcoeff(2.6, c(1.1, 8), 0.5, 1, 6, "R&R", 1)
Kpcoeff <- function(logP, pKa, fup, BP, type = 3, pred = "P&T", dattype = 0) {
    if (pred == "P&T") {
      pcoeff <- Kpcoeff_PT(logP = logP, pKa = pKa, fup = fup, BP = BP, type = type, dattype = dattype)
    } else if (pred == "Berez") { # Berezhkovskiy
      pcoeff <- Kpcoeff_Berez(logP = logP, pKa = pKa, fup = fup, BP = BP, type = type, dattype = dattype)
    } else if (pred == "pksim") { # standard PK-Sim, Willmann et al. 2008
      pcoeff <- Kpcoeff_pksim(logP = logP, fup = fup, dattype = dattype)
    } else if (pred == "Schmitt") { # Schmitt, Walter 2008
      pcoeff <- Kpcoeff_Schmitt(logP = logP, pKa = pKa, fup = fup, type = type, dattype = dattype)
    } else {
      pcoeff <- Kpcoeff_RR(logP = logP, pKa = pKa, fup = fup, BP = BP, type = type, dattype = dattype)
  }

  return(pcoeff)
}
