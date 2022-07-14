#' Two-trace two-dimensional (2T2D) correlation spectroscopy
#' 
#' \code{corr2t2d} compares a pair of spectra in the form of a cross
#'      correlation analysis.
#'
#' \code{corr2t2d} implements the Two-trace two-dimensional (2T2D) approach
#'     as described by I. Noda (2018) <DOI:10.1016/j.molstruc.2018.01.091>.
#'     The idea is to compare two spectra in a 2D correlation-like
#'     approach which was previously not possible as 2D correlation analysis
#'     usually needs at least three spectra.
#'     
#' @param Sam Numeric vector containing the sample spectrum to be correlated.
#'     Can contain the spectral variable of the sample and reference spectrum
#'     as \code{names}.
#' @param Ref Numeric vector containing the sample spectrum to be correlated.
#'     Can contain the spectral variable of the sample and reference spectrum
#'     as \code{names}.
#' @param Wave Numeric vector containing the spectral variable. Needs to be
#'     specified if names of \code{Sam} and \code{Ref} are undefined.
#' @param preview Logical: Should a 3D preview of the asynchronous codistribution
#'     spectrum be drawn at the end? Uses \code{\link[rgl]{persp3d}} from \pkg{rgl}
#'     package.
#'     
#' @return \code{corr2t2d} returns a list of class "corr2d" containing the
#'     complex correlation matrix (\code{$FT}), the correlation and
#'     disrelation coefficient as a complex matrix ($coef), the sample
#'     \code{$Ref1} and reference spectrum \code{$Ref2} as well as the 
#'     spectral variable \code{$Wave1} and \code{$Wave2}.
#' 
#' @references
#'       I. Noda (2018) <DOI:10.1016/j.molstruc.2018.01.091>
#'       
#' @seealso For plotting of the resulting list containing the 2D correlation
#'     spectra or correlation coefficient see \code{\link{plot_corr2d}} and
#'     \code{\link{plot_corr2din3d}}.
#' 
#' @examples
#'     testdata <- sim2ddata()
#'     
#'     twodtest <- corr2t2d(testdata[4, ], testdata[5, ])
#'     
#'     plot_corr2d(twodtest, Im(twodtest$FT))
#' 
#' @export
corr2t2d <- function(Sam, Ref, Wave = NULL, preview = FALSE)
    {
    if (!identical(length(Sam), length(Ref))) {
      stop("length(Sam) and length(Ref) must be equal")
    }

    if (is.null(names(Sam)) && is.null(names(Ref)) &&
        is.null(Wave)) {
      stop("Spectral variable must be specified at Wave,
           names(Sam) or names(Ref)")
    }
  
    if(is.null(Wave) && is.null(names(Sam))) {
      Wave1 <- as.numeric(names(Ref))
      Wave2 <- as.numeric(names(Ref))
    } else if(is.null(Wave)) {
      Wave1 <- as.numeric(names(Sam))
      Wave2 <- as.numeric(names(Sam))
    } else {
      Wave1 <- Wave
      Wave2 <- Wave
    }
  
    # Calculate 2T2D spectra
    syn2t2d <- Sam %o% Sam + Ref %o% Ref
    asyn2t2d <- Sam %o% Ref - Ref %o% Sam
    
    # Caluclate 2T2D correlation and disrelation coefficient
    corrcoef <- syn2t2d / sqrt(diag(syn2t2d) %o% diag(syn2t2d))
    disrcoef <- asyn2t2d / sqrt(diag(syn2t2d) %o% diag(syn2t2d))
    
    Obj <- list(FT = syn2t2d + asyn2t2d*1i, coef = corrcoef + disrcoef*1i,
                Ref1 = Sam, Ref2 = Ref, Wave1 = Wave1, Wave2 = Wave2)
  
    # 3d preview of the asynchronous 2T2D spectrum ----------------------------
    if (preview == TRUE) {
      if (dim(asyn2t2d)[1] > 700) {
        tmp1 <- round(seq(1, dim(asyn2t2d)[1], length = 700))
      } else {
        tmp1 <- seq(1, dim(asyn2t2d)[1], 1)
      }
      if (dim(asyn2t2d)[2] > 700) {
        tmp2 <- round(seq(1, dim(asyn2t2d)[2], length = 700))
      } else {
        tmp2 <- seq(1, dim(asyn2t2d)[2], 1)
      }
      rgl::persp3d(Wave1[tmp1], Wave2[tmp2], asyn2t2d[tmp1, tmp2], col = "grey")
    }
  
    class(Obj) <- "corr2d"
    return(Obj)
}
