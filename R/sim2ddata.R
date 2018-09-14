#' Simulate kinetic data from two-step sequential first-order reactions
#' 
#' \code{sim2ddata} simulates kinetic data for the sequential reaction
#'     A -> B -> C with the time constants k1 and k2.
#'
#' The simulation assumes 2 spectral signals for each of the 3 species A, B
#'     and C. The sequential reaction is defined by 2 time constants k1 and k2.
#'     The spectral information can be sampled at every point during the
#'     reaction to get an arbitrary profile of the kinetic data. The signals of
#'     the three species are modeled by a normal distribution. In addition the
#'     spectral variable is assumed to be equidistant and the number of spectral
#'     variables can also be chosen arbitrary.
#'     
#' @param L Positive, non-zero integer specifying how many spectral variables
#'     should be used to describe the kinetic dataset.
#' @param t Numeric vector containing non-negative real numbers describing at
#'     which reaction times the kinetic data should be sampled.
#' @param k1,k2 Positive, non-zero real numbers describing the time constants
#'     used to simulate the reactions A -> B (\code{k1}) and B -> C (\code{k2}).
#' @param X Numeric vector with two values specifying the range of the simulated
#'     spectral variables.
#' @param A,B,C Numeric vector with two real values specifying the two signal
#'     positions of species A, B and C, respectively. It's the \code{mean} used
#'     in \code{\link[stats]{dnorm}} to simulate the signal. C and Camp may be
#'     NULL in which case only the reaction A -> B is simulated and sampled.
#' @param Aamp,Bamp,Camp Numeric vector with two values specifying the signal
#'     width of species A, B and C, respectively. It's the standard deviation
#'     (\code{sd}) used in \code{\link[stats]{dnorm}} to simulate the signal. 
#'     C and Camp may be NULL in which case only the reaction A -> B is
#'     simulated and sampled.
#' 
#' @return \code{sim2ddata} returns a matrix containing the kinetic data. The
#'     matrix contains the sampled reaction times by rows and the spectral
#'     variables by columns. The reaction times are the row names while the
#'     spectral variables are saved as the column names. The matrix has the
#'     ideal format to be analyzed by \code{\link[corr2D]{corr2d}}.
#' 
#' @references The default values are inspired by:
#'     I. Noda (2014) <DOI:10.1016/j.molstruc.2014.01.024>
#' 
#' @examples
#'     testdata <- sim2ddata()
#'     
#'     twodtest <- corr2d(testdata, corenumber = 1)
#'     
#'     plot_corr2d(twodtest)
#' 
#' @export
sim2ddata <- function(L = 400, t = 0:10, k1 = 0.2, k2 = 0.8, X = c(1000, 1400),
                      A = c(1080, 1320), Aamp = c(3, 8),
                      B = c(1120, 1280), Bamp = c(5, 15),
                      C = c(1160, 1240), Camp = c(4, 9))
    {
    # check user input for errors -----------------------------------------
    if (L <= 0 || is.complex(L) || L%%1 != 0) {
        stop("L must be a positive, non-zero integer")
    }
    
    if (k1 <= 0 || is.complex(k1)) {
        stop("k1 must be a positive, non-zero real number")
    }
    if (k2 <= 0 || is.complex(k2)) {
        stop("k2 must be a positive, non-zero real number")
    }
    
    if (length(A) != 2 || is.complex(A)) {
        stop("A must have exactly 2 real values")
    }
    if (length(B) != 2 || is.complex(B)) {
        stop("B must have exactly 2 real values")
    }
    if (!is.null(C)) {
        if (is.null(Camp)) {
            stop("If C is not NULL, than Camp also needs to be defined")
        }
        if (length(C) != 2 || is.complex(C)) {
            stop("C must have exactly 2 real values")
        }
    }
    if (length(X) != 2 || is.complex(X)) {
        stop("X must have exactly 2 real values")
    }
    
    if (length(Aamp) != 2 || is.complex(Aamp)) {
        stop("Aamp must have exactly 2 real values")
    }
    if (length(Bamp) != 2 || is.complex(Bamp)) {
        stop("Bamp must have exactly 2 real values")
    }
    if (!is.null(Camp)) {
        if (is.null(C)) {
            stop("If Camp is not NULL, than C also needs to be defined")
        }
        if (length(Camp) != 2 || is.complex(Camp)) {
            stop("Camp must have exactly 2 real values")
        }
    }
    
    
    X <- seq(X[1], X[2], length.out = L)
    
    A <- stats::dnorm(X, A[1], Aamp[1]) + stats::dnorm(X, A[2], Aamp[2])
    B <- stats::dnorm(X, B[1], Bamp[1]) + stats::dnorm(X, B[2], Bamp[2])
    
    if (is.null(C)) {
        At <- exp(-k1 * t) %o% A
        
        Bt <- (1 - exp(-k1 * t)) %o% B
        
        Ct <- 0
    }
    
    if (!is.null(C)) {
    C <- stats::dnorm(X, C[1], Camp[1]) + stats::dnorm(X, C[2], Camp[2])
    
    tstar <- log(k2 / k1)
    
    At <- exp(-k1 * t) %o% A
    
    Bt <- ((exp(-k2 * t) - exp(-k1 * t)) / (exp(-k2 * tstar) - exp(-k1 * tstar))) %o% B
    
    Ct <- (1 - (k2 * exp(-k1 * t) - k1 * exp(-k2 * t)) / (k2 - k1)) %o% C
    }
    
    Gt <- At + Bt + Ct
    
    colnames(Gt) <- X
    rownames(Gt) <- t
    
    return(Gt)
}
