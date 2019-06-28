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
#'     '\emph{spectral variable}' by columns and '\emph{perturbation variables}'
#'     by rows. For hetero correlations \code{Mat1} and \code{Mat2} must have
#'     the same number of rows.
#' @param Ref1,Ref2 Numeric vector containing a single spectrum, which will be
#'     subtracted from \code{Mat1} (or \code{Mat2}, respectively) to generate dynamic spectra
#'     for 2D correlation analysis. Default is \code{NULL} in which case the \code{colMeans()}
#'     of \code{Mat1} (or \code{Mat2}, respectively) is used as reference. The length of \code{Ref1}
#'     (or \code{Ref2}) needs to be equal to the number of columns in \code{Mat1} (or \code{Mat2}).
#' @param Wave1,Wave2 Numeric vector containing the spectral variable. Needs to be
#'     specified if column names of \code{Mat1} (or \code{Mat2}) are undefined.
#' @param Time Numeric vector containing the perturbation variables. If specified, \code{Mat1}
#'     (and \code{Mat2} if given) will be interpolated to \code{N} equally spaced perturbation
#'     varibales using \code{Int} to speed up the fft algorithm.
#' @param Int Function specifying how the dataset will be interpolated to give
#'     \code{N} equally spaced perturbation variables. \code{\link[stats]{splinefun}}
#'     (default) or \code{\link[stats]{approxfun}} can for example be used.
#' @param N Positive, non-zero integer specifying how many equally spaced
#'     perturbation variables should be interpolated using \code{Int}. \code{N}
#'     should be higher than 4.\code{corr2d} is fastest if \code{N} is a power
#'     of 2.
#' @param Norm A number specifying how the correlation matrix should be
#'     normalized.
#' @param scaling Positive real number used as exponent when scaling the dataset
#'     with its standard deviation. Defaults to 0 meaning no scaling. 0.5
#'     (\emph{Pareto scaling}) and 1 (\emph{Pearson scaling}) are commonly used to enhance
#'     weak correlations relative to strong correlations.
#' @param corenumber Positive, non-zero integer specifying how many CPU cores
#'     should be used for parallel fft computation.
#' @param preview Logical: Should a 3D preview of the synchronous correlation
#'     spectrum be drawn at the end? Uses \code{\link[rgl]{persp3d}} from \pkg{rgl}
#'     package.
#'     
#' @return \code{corr2D} returns a list of class "corr2d" containing the complex
#'     correlation matrix (\code{$FT}), the used reference spectra (\code{$Ref1},
#'     \code{$Ref2}), the spectral variables (\code{$Wave1}, \code{$Wave2}), the
#'     Fourier transformed data (\code{$ft1}, \code{$ft2}), the (interpolated)
#'     perturbation variables (\code{$Time}) and logical variable (\code{$Het})
#'     indicating if homo (\code{FALSE}) or hetero (\code{TRUE}) correlation was done.
#' 
#' @references 
#'     I. Noda (1993) <DOI:10.1366/0003702934067694>\cr
#'     I. Noda (2012) <DOI:10.1016/j.vibspec.2012.01.006>
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
             N = 2^ceiling(log2(NROW(Mat1))), Norm = 1/(pi * (NROW(Mat1) - 1)), 
             scaling = 0, corenumber = parallel::detectCores(), preview = FALSE)
    {
        # check user input for errors -----------------------------------------
        if (!is.null(Ref1)) {
            if (!identical(length(Ref1), NCOL(Mat1))) {
                stop("length(Ref1) must be equal to NCOL(Mat1)")
            }
        }
        if (!is.null(Mat2) && !is.null(Ref2)) {
            if (!identical(length(Ref2), NCOL(Mat2))) {
                stop("length(Ref2) must be equal to NCOL(Mat2)")
            }
        }
        
        if (is.null(colnames(Mat1)) &&
            is.null(Wave1)) {
            stop("Spectral variable 1 must be specified at Wave1 or in
                 colnames(Mat1)")
        }
        if (!is.null(Mat2)) {
            if (is.null(colnames(Mat2)) &&
                is.null(Wave2)) {
                stop("Spectral variable 2 must be specified at Wave2 or in
                     colnames(Mat2)")
            }
        }
        
        if (N <= 0 || N%%1 != 0) {
            stop("N must be a positive, non-zero integer")
        }
        if (corenumber <= 0 || corenumber%%1 != 0) {
            stop("corenumber must be a positive, non-zero integer")
        }
        
        if (!is.logical(preview)) {
            stop("preview needs to be logical")
        }
        
        # create an R session for every detected CPU core for parallel computing
        cl <- parallel::makeCluster(corenumber)
        doParallel::registerDoParallel(cl)
        
        # define Wave1/2 from colnames if needed ------------------------------
        if (is.null(Wave1)) {
            Wave1 <- as.numeric(colnames(Mat1))
        }
        if (!is.null(Mat2) && !identical(Mat1, Mat2) && is.null(Wave2)) {
            Wave2 <- as.numeric(colnames(Mat2))
        }
        
        # spectral variable must be increasing, also switch Ref if needed ------
        if (is.unsorted(Wave1)) {
            Wave1 <- Wave1[length(Wave1):1]
            Mat1 <- Mat1[, NCOL(Mat1):1]
            if (!is.null(Ref1)) {
                Ref1 <- Ref1[length(Ref1):1]
            }
        }
        if (!is.null(Mat2) && !identical(Mat1, Mat2) && is.unsorted(Wave2)) {
            Wave2 <- Wave2[length(Wave2):1]
            Mat2 <- Mat2[, NCOL(Mat2):1]
            if (!is.null(Ref2)) {
                Ref2 <- Ref2[length(Ref2):1]
            }
        }
        
        # define Wave2 (if needed), D1/2, Het ---------------------------------
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
        }
        
        # Interpolate values of perturbation variable for equidistant ---------
        # datapoints for fft --------------------------------------------------
        if ((is.null(Time) == FALSE) & (length(Time) == D1[1])) {
            cat(c(format(Sys.time(), "%X -"), "Interpolate Time from",
                  min(Time), "to", max(Time), "\n", "to obtain", N,
                  "equidistant data points for FFT", "\n"))
            TIME <- seq(min(Time), max(Time), length.out = N)
            
            tmp <- apply(Mat1, 2, function(y) Int(x = Time, y = y))
            Mat1 <- sapply(tmp, function(x) x(TIME))
            
            if (Het == FALSE) {
                Mat2 <- Mat1
            } else {
                tmp <- apply(Mat2, 2, function(y) Int(x = Time, y = y))
                Mat2 <- sapply(tmp, function(x) x(TIME))
            }
            Time <- TIME
        }
        # Get perturbation variables from rownames ----------------------------
        if (is.null(Time) && !is.null(rownames(Mat1))) {
            Time <- as.numeric(rownames(Mat1))
        }
        
        # Subtract reference -------------------------------------------------
        if (is.null(Ref1)) {
            cat(c(format(Sys.time(), "%X -"), "using mean values as reference\n"))
            Ref1 <- colMeans(Mat1)
        }
        Mat1 <- sweep(Mat1, 2, Ref1, "-")
        if (Het == FALSE) {
            Mat2 <- Mat1
        } else {
          if (is.null(Ref2)) {
            Ref2 <- colMeans(Mat2)
          }
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
        
        # Do fft for every wavenumber in both matrices to obtain frequency ----
        # dependent vectors for every wavenumber and scalar-multiply them -----
        # for every point in the 2D-correlation spectrum ----------------------
        i <- NULL # Dummy to trick R CMD check 
        FT <- NULL
        cat(c(format(Sys.time(), "%X -"), "Fast Fourier Transformation and multiplication \n",
              "to obtain a", D1[2], "x", D2[2], "correlation matrix", "\n"))
        
        FT <- matrix(NA, NCOL(Mat1), NCOL(Mat2))
        # Selection of fft elements
        ftele1 <- if(NROW(Mat1) %% 2 == 0) {
            1:(NROW(Mat1) - 1) %/% 2 + 1
          } else {
            1:(NROW(Mat1)) %/% 2 + 1
          }
        
        ft1 <- foreach::foreach(i = 1:NCOL(Mat1), .combine = 'cbind') %dopar% {
            fft(Mat1[, i])[ftele1]
        }
        if (Het == FALSE) {
            ft2 <- ft1
        } else {
          # Selection of fft elements
          ftele2 <- if(NROW(Mat2) %% 2 == 0) {
            1:(NROW(Mat2) - 1) %/% 2 + 1
          } else {
            1:(NROW(Mat2)) %/% 2 + 1
          }
          
          ft2 <- foreach::foreach(i = 1:NCOL(Mat2), .combine = 'cbind') %dopar% {
              fft(Mat2[, i])[ftele2]
          }
        }
        
        FT <- matrix(Norm * parallel::parCapply(cl, ft1, get("%*%"), Conj(ft2)), NCOL(ft1), NCOL(ft2), byrow = TRUE)
        cat(c(format(Sys.time(), "%X -"), "Done\n"))
        
        Obj <- list(FT = FT, Ref1 = Ref1, Ref2 = Ref2,
                  Wave1 = Wave1, Wave2 = Wave2, ft1 = ft1, ft2 = ft2,
                  Time = Time, Het = Het)
        
        parallel::stopCluster(cl)
        
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

#' @method summary corr2d
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

#' @method print corr2d
#' @export
print.corr2d <- function(x, ...)
{
    print.default(x, ...)
  
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

