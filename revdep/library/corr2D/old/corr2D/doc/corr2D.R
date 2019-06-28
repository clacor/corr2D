## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)

## ------------------------------------------------------------------------
library("corr2D")
data(FuranMale, package = "corr2D")

## ------------------------------------------------------------------------
twod <- corr2d(FuranMale, Ref1 = FuranMale[1, ], corenumber = 1)

## ---- eval=FALSE, fig.show='hold', fig.cap = "Your figure caption."------
#  plot_corr2d(twod, xlab = expression(paste("relative Wavenumber" / cm^-1)),
#              ylab = expression(paste("relative Wavenumber" / cm^-1)))

