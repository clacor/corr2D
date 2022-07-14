#' Two-dimensional codistribution spectroscopy.
#'
#' \code{codis2d} calculates the synchronous and asynchronous codistribution
#'     spectra.
#'
#' \code{codis2d} calculates the the synchronous 2D correlation spectrum and
#'     uses the 2D spectrum to calculate the synchronous and asynchronous
#'     codistribution spectra. For parallelization the
#'     \code{\link[parallel]{parCapply}} function is used. Large input matrices
#'     (> 4000 columns) can lead to long calculation times depending on the
#'     number of cores used. Also note that the resulting matrix can become
#'     very large, adjust the RAM limit with \code{\link[utils]{memory.limit}}
#'     accordingly. For a detailed description of the underlying math see
#'     references.
#'
#' @param Mat Numeric matrix containing the data which will be correlated;
#'     '\emph{spectral variable}' by columns and '\emph{perturbation variables}'
#'     by rows.
#' @param Ref Numeric vector containing a single spectrum, which will be
#'     subtracted from \code{Mat} to generate dynamic spectra for 2D
#'     correlation analysis. Default is \code{NULL} in which case the
#'     \code{colMeans()} of \code{Mat} is used as reference. The length of
#'     \code{Ref} needs to be equal to the number of columns in \code{Mat}.
#'     2D codistribution spectroscopy is only strictly defined using the
#'     perturbation-mean spectrum as reference spectrum. Thus, any deviation
#'     from this definition can lead to unexpected results.
#' @param Wave Numeric vector containing the spectral variable. Needs to be
#'     specified if column names of \code{Mat} are undefined.
#' @param Time Numeric vector containing the perturbation variables. If
#'     specified, \code{Mat} will be interpolated to \code{N} equally spaced
#'     perturbation variables using \code{Int}.
#' @param Int Function specifying how the dataset will be interpolated to give
#'     \code{N} equally spaced perturbation variables. \code{\link[stats]{splinefun}}
#'     (default) or \code{\link[stats]{approxfun}} can for example be used.
#' @param N Positive, non-zero integer specifying how many equally spaced
#'     perturbation variables should be interpolated using \code{Int}. \code{N}
#'     should be higher than 4. \code{corr2d} is fastest if \code{N} is a power
#'     of 2.
#' @param Norm A number specifying how the correlation matrix should be
#'     normalized.
#' @param scaling Positive real number used as exponent when scaling the dataset
#'     with its standard deviation. Defaults to 0 meaning no scaling. 0.5
#'     (\emph{Pareto scaling}) and 1 (\emph{Pearson scaling}) are commonly used
#'     to enhance weak correlations relative to strong correlations. 2D
#'     codistribution spectroscopy is only strictly defined without the usage
#'     of any scaling techniques. Thus, any deviation from this definition can
#'     lead to unexpected results.
#' @param corenumber Positive, non-zero integer specifying how many CPU cores
#'     should be used for parallel fft computation.
#' @param preview Logical: Should a 3D preview of the asynchronous codistribution
#'     spectrum be drawn at the end? Uses \code{\link[rgl]{persp3d}} from \pkg{rgl}
#'     package.
#'     
#' @return \code{codis2D} returns a list of class "corr2d" containing the complex
#'     codistribution matrix (\code{$FT}), the synchronous correlation spectrum
#'     ($corr), the used reference spectrum \code{$Ref1} and \code{$Ref2}, the 
#'     spectral variables \code{$Wave1} and \code{$Wave2} as well as the
#'     (interpolated) perturbation variables (\code{$Time}).
#' 
#' @references 
#'     I. Noda (2014) <DOI:10.1016/j.molstruc.2014.01.024>
#'     
#' @seealso For plotting of the resulting list containing the 2D codistribution
#'     spectra see \code{\link{plot_corr2d}} and \code{\link{plot_corr2din3d}}.
#' 
#' @examples
#'     testdata <- sim2ddata(C = NULL, Camp = NULL)
#'     codis <- codis2d(testdata, corenumber = 1)
#'     
#'     plot_corr2d(codis, Im(codis$FT),
#'                 xlab = expression(paste("Wavenumber" / cm^-1)),
#'                 ylab = expression(paste("Wavenumber" / cm^-1)))
#'
#' @aliases codis2d
#' 
#' @export
#' @importFrom stats sd
codis2d <-
    function(Mat, Ref = NULL, Wave = NULL, Time = NULL, Int = stats::splinefun,
             N = 2^ceiling(log2(NROW(Mat))), Norm = 1/(NROW(Mat) - 1), 
             scaling = 0, corenumber = parallel::detectCores(), preview = FALSE)
    {
        # check user input for errors -----------------------------------------
        if (!is.null(Ref)) {
            if (!identical(length(Ref), NCOL(Mat))) {
                stop("length(Ref) must be equal to NCOL(Mat)")
            }
        }
        
        if (is.null(colnames(Mat)) &&
            is.null(Wave)) {
            stop("Spectral variable must be specified at Wave or in
                 colnames(Mat)")
        }
        
        if (N <= 0 || N%%1 != 0) {
            stop("N must be a positive, non-zero integer")
        }
        if (corenumber <= 0 || corenumber%%1 != 0) {
            stop("corenumber must be a positive, non-zero integer")
        }
        
        if (!is.logical(preview)) {
            stop("preview needs to logical")
        }
        
        # create an R session for every detected CPU core for parallel computing
        cl <- parallel::makeCluster(corenumber)
        doParallel::registerDoParallel(cl)
        
        # define Wave from colnames if needed ---------------------------------
        if (is.null(Wave)) {
            Wave <- as.numeric(colnames(Mat))
        }
        
        # spectral variable must be increasing, also switch Ref if needed ------
        if (is.unsorted(Wave)) {
            Wave <- Wave[length(Wave):1]
            Mat <- Mat[, NCOL(Mat):1]
            if (!is.null(Ref)) {
                Ref <- Ref[length(Ref):1]
            }
        }
        
        D1 <- dim(Mat)

        # Interpolate values of perturbation variable for equidistant ---------
        # datapoints for fft --------------------------------------------------
        if ((is.null(Time) == FALSE) & (length(Time) == D1[1])) {
            cat(c(format(Sys.time(), "%X -"), "Interpolate Time from",
                  min(Time), "to", max(Time), "\n", "to obtain", N,
                  "equidistant datapoints for FFT", "\n"))
            TIME <- seq(min(Time), max(Time), length.out = N)
            
            tmp <- apply(Mat, 2, function(y) Int(x = Time, y = y))
            Mat <- sapply(tmp, function(x) x(TIME))

            Time <- TIME
        }
        # Get perturbation variables from rownames ----------------------------
        if (is.null(Time) && !is.null(rownames(Mat))) {
            Time <- as.numeric(rownames(Mat))
        }
        
        # Subtract reference -------------------------------------------------
        if (is.null(Ref)) {
            cat(c(format(Sys.time(), "%X -"), "using mean values as reference\n"))
            Ref <- tavspec <- colMeans(Mat)
        } else {
            tavspec <- colMeans(Mat)
        }
        Mat <- sweep(Mat, 2, Ref, "-")
        
        # Apply scaling -------------------------------------------------------
        if (scaling > 0) {
            cat(c(format(Sys.time(), "%X -"), "apply scaling\n"))
            sd1 <- apply(Mat, 1, sd)
            Mat <- Mat / (sd1^scaling)
        }
        
        Syncorr <- matrix(Norm * parallel::parCapply(cl, Mat, get("%*%"), Mat),
                           NCOL(Mat), NCOL(Mat), byrow = T)

        cat(c(format(Sys.time(), "%X -"), "Done calculating synchronous 2D spectrum\n"))
        
        parallel::stopCluster(cl)
        
        # build the sum weigthed sum spectrum
        su <- colSums(1:NROW(Mat) * Mat)
        
        # calculate characteristic perturbation index k and -------------------
        # characteristic perturbation of the distribution t -------------------
        k <- 1 / (NROW(Mat) * tavspec) * su  + (NROW(Mat) + 1) / 2
        dist <- (max(Time) - min(Time)) * (k - 1) / (NROW(Mat) - 1) + min(Time)
        
        # calculate total joint variance tjv from synchronous correlation spectrum
        tjv <- sqrt(outer(diag(Syncorr), diag(Syncorr)))
        tmat <- -1 * outer(dist, dist, "-")
        
        # calculate asynchronous and synchronous codistribution spectra -------
        # from total joint variance tjv and characteristic perturbation of the 
        # distribution t, respectivly tmat ------------------------------------
        Asyncds <- (tmat * tjv) / (max(Time) - min(Time))
        Syncds <- sqrt(tjv^2 - Asyncds^2)
        FT <- Syncds + Asyncds * 1i
        cat(c(format(Sys.time(), "%X -"), "Finished\n"))
        
        Obj <- list(FT = FT, corr = Syncorr, Ref1 = Ref, Ref2 = Ref,
                    Wave1 = Wave, Wave2 = Wave, Time = Time)
        
        # 3d preview of the asynchronous codistribution spectrum --------------
        if (preview == TRUE) {
            if (dim(Asyncds)[1] > 700) {
                tmp1 <- round(seq(1, dim(Asyncds)[1], length = 700))
            } else {
                tmp1 <- seq(1, dim(Asyncds)[1], 1)
            }
            rgl::persp3d(Wave[tmp1], Wave[tmp1], Asyncds[tmp1, tmp1],
                          col = "grey")
        }
        
        class(Obj) <- "corr2d"
        return(Obj)
    }