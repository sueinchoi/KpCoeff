
# mrgsolve <img align="right" src = "man/figures/mrgsolve_sticker_812418_1.png" width="135px">

KpCoeff is an R package for prediction of tissue-plasma partition coefficients (Kp) using published method with published/unified dataset. KpCoeff is free and open-source software.

## Resources

Unified dataset is from the reference 

[![Reference](Reference/903.full.pdf)

## Installation

Install the current development version

``` r
library(devtools)
remotes::install_github("sueinchoi/KpCoeff")

```


## Interaction

We welcome **questions** about anything KpCoeff: installation, dataset, how it works, reference, understanding better how KpCoeff works. We also
welcome **suggestions** for how to make KpCoeff more useful to you and
to the PBPK community.


## Some examples

### A simple simulation

``` r
library(KpCoeff)
```

Predict the tissue-plasma partition coefficient using P&T method with reference dataset

``` r
KpCoeff(1.9, 6.4, 0.01, 1, 3, "P&T", 0)
```

Predict the tissue-plasma partition coefficient using R&R method with unified dataset


``` r
KpCoeff(1.9, 6.4, 0.01, 1, 3, "R&R", 1)

```

The result can be direcly incorporated to the simulation model from mrgsolve package

Example of using KpCoeff result with mrgsolve package will be soon updated