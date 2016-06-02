#' Two-dimensional correlation analysis.
#'
#' \code{corr2d} calculates the synchronous and asynchronous correlation
#'     spectra between \code{Mat1} and \code{Mat1} (homo correlation)
#'     or between \code{Mat1} and \code{Mat2} (hetero correlation).
#'
#' \code{corr2d} uses a parallel fast Fourier transformation
#'     (\code{\link[stats]{fft}}) to calculate the complex correlation matrix.
#'     For parallelization the \code{\link[foreach]{foreach}} function is used.
#'     Large input matrices (> 4000 columns) can lead to long calculation times
#'     depending on the number of cores used. Also note that the resulting
#'     matrix can become very large, adjust the RAM limit with
#'     \code{\link[utils]{memory.limit}} accordingly. For a detailed description
#'     of the underlying math see references.
#'
#' @param Mat1,Mat2 Numeric matrix containing the data which will be correlated;
#'     '\emph{Spectral variable}' by columns and '\emph{perturbation}' by rows. For hetero
#'     correlations \code{Mat1} and \code{Mat2} must have the same number of rows.
#' @param Ref1,Ref2 Numeric vector containg a single spectrum, which will be
#'     substracted from \code{Mat1} (or \code{Mat2}, respectivly) to generate dynamic spectra
#'     for 2D correlation analysis. Defaults to \code{NULL} in which case the \code{colMeans()}
#'     of \code{Mat1} (or \code{Mat2}, respectivly) is used as reference. The length of \code{Ref1}
#'     (or \code{Ref2}) needs to be equal to the number of columns in \code{Mat1} (or \code{Mat2}).
#' @param Wave1,Wave2 Numeric vector containing the spectral variable. Needs to be
#'     specified if column names of \code{Mat1} (or \code{Mat2}) are undefinied.
#' @param Time Numeric vector containing the perturbation. If specified, \code{Mat1}
#'     (and \code{Mat2} if given) will be interpolated to \code{N} equally spaced perturbation
#'     values using \code{Int} to speed up the fft algorithm.
#' @param Int Function specifing how the dataset will be interpolated to give
#'     \code{N} equally spaced perturbation values.
#' @param N Positive, non-zero integer specifing how many equally spaced
#'     perturbation values should be interpolated using \code{Int}. \code{corr2d}
#'     is fastest if \code{N} is a power of 2.
#' @param Norm Numeric vector specifing how the correlation matrix should be
#'     normalized. Needs to have length \code{(NCOL(Mat1) - 1)}.
#' @param scaling Positive real number used as exponent when scaling the dataset
#'     with its standard deviation. Defaults to 0 meaning no scaling. 0.5
#'     (\emph{Pareto scaling}) and 1 (\emph{Pearson scaling}) are commonly used to enhance
#'     weak correlations relative to strong correlations.
#' @param corenumber Positive, non-zero integer specifing how many CPU cores
#'     should be used for parallel fft computation.
#' @param preview Logical scalar, should a 3D preview of the synchronous
#'     correlation spectrum be drawn at the end? Uses \code{\link[rgl]{persp3d}} from \pkg{rgl}
#'     package.
#'     
#' @return \code{corr2D} returns a list of class "corr2d" containing the complex
#'     correlation matrix (\code{$FT}), the used reference spectra (\code{$Ref1},
#'     \code{$Ref2}), the spectral variables (\code{$Wave1}, \code{$Wave2}), the
#'     (interpolated) perturbation (\code{$Time}) and logical scalar (\code{$Het})
#'     indicating if homo (\code{FALSE}) or hetero (\code{TRUE}) correlation was done.
#' 
#' @references 
#'     I. Noda, \emph{Appl. Spectrosc.}, 1993, \strong{47}, 1329-1336.\cr
#'     I. Noda, \emph{Vib. Spectrosc.}, 2012, \strong{60}, 146-153.
#'     
#' 
#' @seealso For plotting of the resulting list containing the 2D correlation
#'     spectra see \code{\link{plot_corr2d}} and \code{\link{plot_corr2din3d}}.
#' 
#' @examples
#'     data(FuranMale, package = "corr2D")
#'     twod <- corr2d(FuranMale, Ref1 = FuranMale[1, ], corenumber = 1)
#'     
#'     plot_corr2d(twod, xlab = expression(paste("relative Wavenumber" / cm^-1)),
#'                       ylab = expression(paste("relative Wavenumber" / cm^-1)))
#'
#' @aliases corr2d
#' 
#' @export
#' @importFrom foreach %dopar%
#' @importFrom stats fft sd
corr2d <-
    function(Mat1, Mat2 = NULL, Ref1 = NULL, Ref2 = NULL,
             Wave1 = NULL, Wave2 = NULL, Time = NULL, Int = stats::splinefun,
             N = 2^ceiling(log2(NROW(Mat1))), Norm = 1/(pi * (NCOL(Mat1) - 1)), 
             scaling = 0, corenumber = parallel::detectCores(), preview = FALSE)
    { 
        # create an R session for every detected CPU core for parallel computing
        cl <- parallel::makeCluster(corenumber)
        doParallel::registerDoParallel(cl)
        
        # spectral variable must be incresing, also switch Ref if needed ------
        if (is.unsorted(as.numeric(colnames(Mat1)))) {
            Mat1 <- Mat1[, NCOL(Mat1):1]
            if (!is.null(Ref1)) {
                Ref1 <- Ref1[length(Ref1):1]
            }
        }
        if (!is.null(Mat2) && is.unsorted(as.numeric(colnames(Mat2)))) {
            Mat2 <- Mat2[, NCOL(Mat2):1]
            if (!is.null(Ref2)) {
                Ref2 <- Ref2[length(Ref2):1]
            }
        }
        
        # define Wave1/2, D1/2, Het -------------------------------------------
        if (is.null(Wave1)) {
            Wave1 <- as.numeric(colnames(Mat1))
        }
        if (is.null(Mat2) || identical(Mat1, Mat2)) {
            cat(c("HOMO-Correlation:", corenumber, "cores used for calculation\n"))
            Het <- FALSE
            D1 <- D2 <- dim(Mat1)
            Wave2 <- Wave1
            Mat2 <- Mat1
            Ref2 <- Ref1
        } else {
            cat(c("HETERO-Correlation:", corenumber, "cores used for calculation\n"))
            Het <- TRUE
            D1 <- dim(Mat1)
            D2 <- dim(Mat2)
            if (is.null(Wave2)) {
                Wave2 <- as.numeric(colnames(Mat2))
            }
        }
        
        # Interpolate values of perturbation variable for equidistant ---------
        # datapoints for fft --------------------------------------------------
        if ((is.null(Time) == FALSE) & (length(Time) == D1[1])) {
            cat(c(format(Sys.time(), "%X -"), "Interpolate Time from",
                  min(Time), "to", max(Time), "\n", "to obtain", N,
                  "equidistant datapoints for FFT", "\n"))
            MAT <- c()
            TIME <- seq(min(Time), max(Time), length.out = N)
            for (i in 1:NCOL(Mat1)) {
                tmp <- Int(x = Time, y = Mat1[, i])
                MAT <- cbind(MAT, tmp(TIME))
            }
            Mat1 <- MAT
            if (Het == FALSE) {
                Mat2 <- Mat1
            } else {
                MAT <- c()
                for (i in 1:NCOL(Mat2)) {
                    tmp <- Int(x = Time, y = Mat2[, i])
                    MAT <- cbind(MAT, tmp(TIME))
                }
                Mat2 <- MAT
            }
        }
        # Get perturbation variable from rownames
        if (is.null(Time) == TRUE && is.null(rownames(Mat1)) == FALSE) {
            Time <- as.numeric(rownames(Mat1))
        } else {
            Time <- NULL
        }
        
        # Substract reference -------------------------------------------------
        if (is.null(Ref1)) {
            cat(c(format(Sys.time(), "%X -"), "using mean values as reference\n"))
            Ref1 <- colMeans(Mat1)
            Ref2 <- colMeans(Mat2)
        }
        Mat1<-sweep(Mat1, 2, Ref1, "-")
        if (Het == FALSE) {
            Mat2 <- Mat1
        } else {
            Mat2 <- sweep(Mat2, 2, Ref2, "-")
        }
        
        # Apply scaling -------------------------------------------------------
        if (scaling > 0) {
            cat(c(format(Sys.time(), "%X -"), "apply scaling\n"))
            sd1 <- apply(Mat1, 1, sd)
            sd2 <- apply(Mat2, 1, sd)
            Mat1 <- Mat1 / (sd1^scaling)
            Mat2 <- Mat2 / (sd2^scaling)
        }
        
        # Do fft for every wavenumber in both Mattices to obtain frequency ----
        # dependent vectors for every wavenumber and scalar-multiply them -----
        # for every point in the 2D-correlation spectrum ----------------------
        FT <- NULL
        cat(c(format(Sys.time(), "%X -"), "Fast Fourier Transformation and Multiplication \n",
              "to obtain a",D1[2], "x", D2[2], "Correlation-Matrix","\n"))
        
        FT <- matrix(NA, NCOL(Mat1), NCOL(Mat2))
        ft1 <- foreach::foreach(i = 1:NCOL(Mat1), .combine = 'cbind') %dopar% {
            fft(Mat1[, i])[1:NROW(Mat1) %/% 2]
        }
        if (Het == FALSE) {
            ft2 <- ft1
        } else {
            ft2 <- foreach::foreach(i = 1:NCOL(Mat2), .combine = 'cbind') %dopar% {
                fft(Mat2[, i])[1:NROW(Mat2) %/% 2]
            }
        }
        
        FT<-matrix(Norm * parallel::parCapply(cl, ft1, get("%*%"), Conj(ft2)), NCOL(ft1), NCOL(ft2), byrow = T)
        cat(c(format(Sys.time(), "%X -"), "Done\n"))
        
        Obj<-list(FT = FT, Ref1 = Ref1, Ref2 = Ref2,
                  Wave1 = Wave1, Wave2 = Wave2,
                  Time = Time, Het = Het)
        
        parallel::stopCluster(cl)
        #closeAllConnections()
        
        # 3d preview of the synchronous correlation spectrum ------------------
        if (preview == TRUE) {
            if (dim(FT)[1] > 700) {
                tmp1 <- round(seq(1, dim(FT)[1], length = 700))
            } else {
                tmp1 <- seq(1, dim(FT)[1], 1)
            }
            if (dim(FT)[2] > 700) {
                tmp2 <- round(seq(1, dim(FT)[2], length = 700))
            } else {
                tmp2 <- seq(1, dim(FT)[2], 1)
            }
            rgl::persp3d(Wave1[tmp1], Wave2[tmp2], Re(FT)[tmp1, tmp2], col = "grey")
        }
        
        class(Obj) <- "corr2d"
        return(Obj)
    }

#' @export
summary.corr2d <- function(object, ...)
    {
        cat(c(NROW(object$FT), "x", NCOL(object$FT), if (object$Het == TRUE){"Hetero"} else {"Homo"}, "correlation spectra\n"))
        if (object$Het == TRUE) {
            cat(c("Spectral variable 1:", min(object$Wave1), "-", max(object$Wave1), "\n"))
            cat(c("Spectral variable 2:", min(object$Wave2), "-", max(object$Wave2), "\n"))
        } else {
            cat(c("Spectral variable:", min(object$Wave1), "-", max(object$Wave1), "\n"))
        }
        if (!is.null(object$Time)) {
            cat(c("Perturbation variable:", length(object$Time), "values,", min(object$Time), "-", max(object$Time), "\n"))
        }
    
}

#' Check for object class "corr2d"
#'
#' The function checks if an object is of class "corr2d".
#' 
#' The function uses the \code{\link{inherits}} function.
#' 
#' @param x An object which should be check if it is of class "corr2d".
#' 
#' @return A logical scalar
#'
#' @examples 
#'     data(FuranMale, package = "corr2D")
#'     twod <- corr2d(FuranMale, Ref1 = FuranMale[1, ], corenumber = 1)
#'     
#'     # TRUE
#'     is.corr2d(twod) 
#'     # FALSE
#'     is.corr2d(2) 
#'
#' @export
is.corr2d <- function(x)
    {
    inherits(x, "corr2d", which = FALSE)
    }