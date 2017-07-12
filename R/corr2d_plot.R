#' Plot two-dimensional correlation spectra.
#'
#' \code{plot_corr2d} plots two-dimensional correlation spectra either
#'     as an image or a contour plot. Red color indicates positive
#'     correlations, while blue color shows negative ones.
#'
#' For the synchronous correlation spectrum the real component (\code{Re})
#'     of the complex correlation matrix must be plotted. The asynchronous
#'     spectrum is the respective imaginary component (\code{Im}).
#'     \code{Cutout} can be used to leave out smaller (noise) contributions,
#'     but should be used with care as it can be used to create misleading
#'     2D correlation plots. See references for interpretation rules (so
#'     called Noda rules).
#' 
#' @param Obj List from \code{corr2d} containing the 2D correlation data.
#' @param what Real numeric matrix containing the z-values that should be plotted.
#' @param specx,specy Numeric vector containing the data that should be plotted
#'     on top (\code{specx}) and/or on the left (\code{specy}) of
#'     the 2D spectrum. \code{Mat}, \code{specx} and/or \code{specy} should
#'     have the same dimensions, respectively. If \code{NULL} nothing will
#'     be plotted.
#' @param xlim,ylim Numeric vector with two values indicating the borders
#'     of the 2D plot. Also truncates \code{specx} and/or \code{specy} to
#'     match the new plot range.
#' @param xlab,ylab Character or expression containing the text that will
#'     be plotted on the bottom (\code{xlab}) and/or to the right
#'     (\code{ylab}) of the 2D plot. Labels can be suppressed with \code{NA}.
#' @param Contour Logical: Should a contour (\code{TRUE}) or image
#'     (\code{FALSE}) be drawn?
#' @param axes Integer ranging from 0 to 3. Should the axis of the 2D plot
#'     be drawn? "0" means no axes, "1" only bottom axis, "2" only right axis and
#'     "3" both axes are drawn.
#' @param Legend Logical: Should a color legend be plotted in the top
#'     right corner?
#' @param N Positive, non-zero integer indicating how many contour or image
#'     levels should be plotted.
#' @param zlim Numeric vector with two values defining the z-range of the 2D
#'     plot.
#' @param Cutout Numeric vector with two values defining which z-values should
#'     not be plotted. Use with care, because this can generate misleading
#'     2D plots.
#' @param col A specification for the plotting color of the reference spectra
#'     (top and left), axes, axes ticks and the central plot surrounding box.
#'     See \code{\link[graphics]{par}} and \code{\link[graphics]{contour}} for
#'     additional information.
#' @param lwd A numeric value which sets the line width in the contour plot. See
#'     \code{\link[graphics]{par}} and \code{\link[graphics]{contour}} for
#'     additional information.
#' @param lwd.axis A numeric value which sets the line width for axes and the
#'     central plot surrounding box. See \code{\link[graphics]{par}} and
#'     \code{\link[graphics]{axis}} for additional information.
#' @param lwd.spec A numeric value which sets the line width in the reference
#'     spectra on top and to the left. See \code{\link[graphics]{par}} and 
#'     \code{\link[graphics]{plot.default}} for additional information.
#' @param cex.leg A numerical value giving the amount by which numbers at the
#'     legend should be magnified. See \code{\link[graphics]{par}} and 
#'     \code{\link[fields]{image.plot}} for additional information.
#' @param at.xaxs,at.yaxs The points at which tick-marks are to be drawn at the
#'     x- and y-axis, respectively. See \code{\link[graphics]{axis}} for
#'     additional information.
#' @param label.xaxs,label.yaxs This can either be a logical value specifying
#'     whether (numerical) annotations are to be made at the tickmarks of the
#'     x- and y-axis, or a character or expression vector of labels to be
#'     placed at the tickpoints of the x- and y-axis. See
#'     \code{\link[graphics]{axis}} for additional information.
#' @param line.xlab,line.ylab Numeric value on which MARgin line the x- and
#'     y-label is plotted, respectively, starting at 0 counting outwards. See
#'     \code{\link[graphics]{mtext}} for additional information.
#' @param ... Additional arguments either passed to
#'     \code{\link[graphics]{image}} or \code{\link[graphics]{contour}}. Can
#'     include graphics parameters \code{\link[graphics]{par}} which are in
#'     part also used by other functions. This includes \code{cex.axis}
#'     (influences axes and thier label magnification), \code{cex.lab}
#'     (influences label magnification), \code{col.axis} (influences axes
#'     label color), \code{col.lab} (influences label color),
#'     \code{font.axis} (influences axes label font), \code{font.lab}
#'     (influences label font) and \code{lty} (influences line type for contour
#'     plot).
#'
#' @references For interpretation rules see:
#'     I. Noda (2006) <DOI:10.1016/j.molstruc.2005.12.060>
#' 
#' @seealso See \code{\link{plot_corr2din3d}} for 3D plots.
#'
#' @aliases plot_corr2d
#' 
#' @examples
#'     data(FuranMale, package = "corr2D")
#'     twod <- corr2d(FuranMale, Ref1 = FuranMale[1, ], corenumber = 1)
#'     
#'     plot_corr2d(twod, xlab = expression(paste("relative Wavenumber" / cm^-1)),
#'                       ylab = expression(paste("relative Wavenumber" / cm^-1)))
#'                       
#'     plot_corr2d(twod, at.xaxs = c(1560, 1585, 1610),
#'                 label.xaxs = c(1560, 1585, 1610),
#'                 col = 2, lwd = 3, col.axis = 3, col.lab = 4, Legend = FALSE,
#'                 cex.lab = 3, xlab = "Large x label", ylab = "Large y label",
#'                 line.xlab = 5, line.ylab = 5)
#'
#' @export
#' @importFrom grDevices rgb
#' @importFrom graphics abline axis box close.screen mtext par plot.default screen split.screen
#' @importFrom stats quantile
#' @importFrom utils modifyList
plot_corr2d <-
    function(Obj, what = Re(Obj$FT), specx = Obj$Ref1, specy = Obj$Ref2,
             xlim = NULL, ylim = NULL,
             xlab = expression(nu[1]), ylab = expression(nu[2]),
             Contour = TRUE, axes = 3, Legend = TRUE, N = 20,
             zlim = NULL, Cutout = NULL, col = par("col"), lwd = par("lwd"), 
             lwd.axis = NULL, lwd.spec = NULL, cex.leg = NULL,
             at.xaxs = NULL, label.xaxs = TRUE,
             at.yaxs = NULL, label.yaxs = TRUE,
             line.xlab = 3.5, line.ylab = 3.5, ...)
    {
        # check user input for errors -----------------------------------------
        if (!is.null(xlim)) {
            if (length(xlim) != 2 || is.complex(xlim) || xlim[1] > xlim[2]) {
                stop("xlim must have exactly 2 real values or must be NULL.
                     The first value needs to be smaller than the second one.")
            }
        }
        if (!is.null(ylim)) {
            if (length(ylim) != 2 || is.complex(ylim) || ylim[1] > ylim[2]) {
                stop("ylim must have exactly 2 real values or must be NULL.
                     The first value needs to be smaller than the second one.")
            }
        }
        if (!is.null(zlim)) {
            if (length(zlim) != 2 || is.complex(zlim) || zlim[1] > zlim[2]) {
                stop("zlim must have exactly 2 real values or must be NULL.
                     The first value needs to be smaller than the second one.")
            }
        }
        
        if (!is.logical(Contour)) {
            stop("Contour needs to be logical")
        }
        
        if (axes %in% c(0, 1, 2, 3) == FALSE) {
            stop("axes my only be 0 (no axes), 1 (only bottom axis),
                 2 (only right axis) or 3 (both axes)")
        }
        
        if (!is.logical(Legend)) {
            stop("Legend needs to be logical")
        }
        
        if (N <= 0 || N%%1 != 0) {
            stop("N must be a positive, non-zero integer")
        }
        
        if (!is.null(Cutout)) {
            if (length(Cutout) != 2 || is.complex(Cutout) || Cutout[1] > Cutout[2]) {
                stop("Cutout must have exactly 2 real values or must be NULL.
                     The first value needs to be smaller than the second one.")
            }
        }
        
        par_old <- par(no.readonly = TRUE)
        # avoid "invalid screen(1)" error in RStudio --------------------------
        close.screen(all.screens = TRUE)
        graphics::plot.new()
        
        # get graphics parameters from "..."
        getparm <- list(..., lwd = lwd, col = col)
        graphparm <- utils::modifyList(par(), getparm)
        if(is.null(lwd.axis)) {lwd.axis <- graphparm$lwd + 1}
        if(is.null(lwd.spec)) {lwd.spec <- graphparm$lwd}
        if(is.null(cex.leg)) {cex.leg <- graphparm$cex.axis}
        
        # calculate x- and y-window range -------------------------------------
        if (is.null(xlim)) {
            Which1 <- 1:NROW(what)
        } else {
            Which1 <- which(xlim[1] < Obj$Wave1 & Obj$Wave1 < xlim[2])
        }
        
        if (is.null(ylim)) {
            Which2 <- 1:NCOL(what)
        } else {
            Which2 <- which(ylim[1] < Obj$Wave2 & Obj$Wave2 < ylim[2])
        }
        
        # create splitscreen for plotting -------------------------------------
        OFF <- 0.05
        split.screen(rbind(c(0, 0.15 + OFF, 0.15 + OFF, 0.85 - OFF),          # Spectrum left
                           c(0.15 + OFF, 0.85 - OFF, 0.85 - OFF, 1),          # Spectrum top
                           c(0.15 + OFF, 0.85 - OFF, 0.15 + OFF, 0.85 - OFF), # Main
                           c(0.85 - OFF, 1, 0.15 + OFF, 0.85 - OFF),          # right
                           c(0.15 + OFF, 0.85 - OFF, 0, 0.15 + OFF),          # bottom
                           c(0, 0.15 + OFF, 0.85 - OFF, 1),                   # top left
                           c(0.85 - OFF, 1, 0.85 - OFF, 1)                    # Legend top right
                           )
        )
        
        # plot one dimensional spectra top and left ---------------------------
        if (!is.null(specy)) {
            # Spec left -------------------------------------------------------
            screen(1)
            par(xaxt = "n", yaxt = "n", mar = c(0, 0, 0, 0), bty = "n", yaxs = "i")
            plot.default(x = max(specy[Which2]) - specy[Which2],
                         y = Obj$Wave2[Which2], col = graphparm$col, 
                         type = "l", lwd = lwd.spec, ann = FALSE)
        }
        
        if (!is.null(specx)) {
            # Spec top -------------------------------------------------------
            screen(2)
            par(xaxt = "n", yaxt = "n", mar = c(0, 0, 0, 0), bty = "n", xaxs = "i")
            plot.default(x = Obj$Wave1[Which1],
                         y = specx[Which1], col = graphparm$col,
                         type = "l", lwd = lwd.spec, ann = FALSE)
        }
        
        # main Part -----------------------------------------------------------
        screen(3)
        if (is.null(zlim)) {
            zlim <- range(what[Which1, Which2])
        }
        # Number of levels is always odd --------------------------------------
        if (N%%2 == 0){
            N <- N + 1
        }
        # Symmetric distribution of color code --------------------------------
        Where <- seq(-max(abs(zlim)), max(abs(zlim)), length.out = N) 
        
        if (is.null(Cutout)) {
            OM <- which(Where < 0)
            OP <- which(Where > 0)
        } else {
            OM <- which(Where <= Cutout[1])
            OP <- which(Where >= Cutout[2])
        }
        
        COL <- rep("transparent", length(Where))
        COL[OM] <- fields::designer.colors(col = c("darkblue", "cyan"), n = length(OM))
        COL[OP] <- fields::designer.colors(col = c("yellow", "red", "darkred"), n = length(OP))
        COL[(N + 1)/2] <- "transparent"
        COL <- COL[which(zlim[1] < Where & Where < zlim[2])]
        Where <- seq(zlim[1], zlim[2], length.out = length(COL))
        
        par(xaxt = "n", yaxt = "n", mar = c(0, 0, 0, 0), bty = "n", xaxs = "i", yaxs = "i")
        if (Contour == TRUE){
            graphics::contour(x = Obj$Wave1[Which1], y = Obj$Wave2[Which2], z = what[Which1, Which2],
                              col = COL, levels = Where, zlim = zlim, drawlabels = FALSE,
                              lwd = graphparm$lwd, ...)
        } else {
            graphics::image(x = Obj$Wave1[Which1], y = Obj$Wave2[Which2], z = what[Which1, Which2],
                            col = COL, xlab = "", ylab = "", zlim = zlim,
                            lwd = graphparm$lwd, ...)
        }
        
        abline(a = 0, b = 1, col = rgb(red = 1, green = 1, blue = 1, alpha = 0.5),
            lwd = graphparm$lwd)
        par(xpd = NA, xaxt = "s", yaxt = "s", xaxs = "i", yaxs = "i", cex = graphparm$cex,
            mar = c(0, 0, 0, 0))
        box(which = "figure", lwd = lwd.axis, col = graphparm$col)
        if ((axes == 1) | (axes == 3)){
            axis(side = 1, at = at.xaxs, labels = label.xaxs, lwd = lwd.axis,
                 col = graphparm$col, col.ticks = graphparm$col,
                 cex.axis = graphparm$cex.axis, col.axis = graphparm$col.axis,
                 font.axis = graphparm$font.axis)
        }
        if ((axes == 2) | (axes == 3)){
            axis(side = 4, las = 2, at = at.yaxs, labels = label.yaxs,
                 lwd = lwd.axis, col = graphparm$col, col.ticks = graphparm$col,
                 cex.axis = graphparm$cex.axis, col.axis = graphparm$col.axis,
                 font.axis = graphparm$font.axis)
        }
        
        mtext(side = 1, xlab, line = line.xlab, cex = graphparm$cex.lab * 1.3,
              col = graphparm$col.lab, font = graphparm$font.lab)
        mtext(side = 4, ylab, line = line.ylab, cex = graphparm$cex.lab * 1.3,
              col = graphparm$col.lab, font = graphparm$font.lab)
        
        if(Legend == TRUE){
            # top right -------------------------------------------------------
            screen(7)
            # avoid par(par.old) error from image.plot() by setting par(pin) value positive
            par(pin = abs(par()$pin))
            
            if (Contour == TRUE){
                fields::image.plot(z = what[Which1,Which2], legend.only = TRUE,
                    smallplot = c(0.15, 0.3, 0.2, 0.8), col = COL,
                    axis.args = list(at = quantile(Where, prob = c(0.1, 0.9)),
                        labels = format(x = quantile(Where, prob = c(0.1, 0.9)),
                        digit = 2, scientific = TRUE), cex.axis = cex.leg),
                    zlim = zlim, graphics.reset = TRUE)
            } else {
                fields::image.plot(z = what[Which1, Which2],legend.only = TRUE,
                    smallplot = c(0.15, 0.3, 0.2, 0.8), col = COL,
                    axis.args = list(at = range(what[Which1, Which2]),
                        labels = format(x = range(what[Which1, Which2]),
                        digits = 2, scientific = TRUE), cex.axis = cex.leg),
                    graphics.reset = TRUE)
            }
        
        }
        
        screen(3, new = FALSE)
        close.screen(c(1,2,4,5,6,7))
        on.exit(options(par_old), add = TRUE)
    }

#' 3D plot of two-dimensional correlation spectra.
#'
#' \code{plot_corr2din3d} plots two-dimensional correlation spectra as an 3D surface.
#'
#' For the synchronous correlation spectrum the real component (\code{Re})
#'     of the complex correlation matrix must be plotted. The asynchronous
#'     spectrum is the respective imaginary component (\code{Im}).
#' 
#' @param Mat Real numeric matrix containing the z-values to plot.
#' @param specx,specy Numeric vector containing the data, that will be
#'     plotted at the x and y axis. Can be any data and does not need to have
#'     the same dimensions as \code{Mat}.
#' @param scalex,scaley A real number which describes how \code{specx}
#'     (or \code{specy}) get scaled. Positive numbers lead to a spectrum
#'     plotted inside the box, while negative numbers lead to a spectrum
#'     plotted outside the box.
#' @param Col Vector containing colors used to plot the 3D plot and the
#'     respective projection. 
#' @param reduce Non-zero rational number describing how to
#'     \code{\link[mmand]{resample}} the data values. Can reduce the 
#'     computational demand and can be used for fast previews.
#' @param zlim Numeric vector with two values indicating the z-range of
#'     the 3D plot.
#' @param projection Logical: Should a 2D projection of the 3D surface
#'     be plotted a the bottom of the box?
#' @param ... Additional arguments passed to \code{\link[fields]{drape.plot}}.
#' 
#' @seealso See \code{\link{plot_corr2d}} for 2D plots.
#'     See \code{\link[fields]{drape.plot}} for information on the plot function.
#' 
#' @aliases plot_corr2din3d
#' 
#' @examples
#'    data(FuranMale, package = "corr2D")
#'    twod <- corr2d(FuranMale, Ref1 = FuranMale[1, ], corenumber = 1)
#'    
#'    plot_corr2din3d(Mat = Re(twod$FT), specx = twod$Ref1,
#'        specy = twod$Ref1, reduce = 2, scalex = -175, scaley = -175,
#'        zlim = c(-1.5, 2.2)*10^-3, projection = FALSE,
#'        border = gray(0.2), theta = 25, phi = 15, add.legend = FALSE,
#'        Col = fields::tim.colors(64))
#'    
#' @export
#' @importFrom graphics par polygon lines
#' @importFrom stats median
#' @importFrom colorspace diverge_hcl
plot_corr2din3d <-
    function(Mat, specx = NULL, specy = NULL,
             scalex = NULL, scaley = NULL,
             Col = colorspace::diverge_hcl(64, h = c(240, 0),
                       c = 100, l = c(20, 100), power = 0.4),
             reduce = NULL, zlim = NULL, projection = FALSE, ...)
    {
        # check user input for errors -----------------------------------------
        if (!is.null(zlim)) {
            if (length(zlim) != 2 || is.complex(zlim) || zlim[1] > zlim[2]) {
                stop("zlim must have exactly 2 real values or must be NULL.
                     The first value needs to be smaller than the second one.")
            }
        }
        
        if (!is.null(scalex)) {
            if (length(scalex) != 1 || is.complex(scalex)) {
                stop("scalex must be a real number")
            }
        }
        if (!is.null(scaley)) {
            if (length(scaley) != 1 || is.complex(scaley)) {
                stop("scaley must be a real number")
            }
        }
        
        if (!is.logical(projection)) {
            stop("Contour needs to logical")
        }
        
        # generate x- and y-axis ----------------------------------------------
        par_old <- par(no.readonly = TRUE)
        x <- 1:NROW(Mat)
        y <- 1:NCOL(Mat)
        
        # resample input matrix -----------------------------------------------
        if (!is.null(reduce)) {
            Which.x <- (1:length(x))[which(1:length(x)%%reduce == 0)]
            Which.y <- (1:length(y))[which(1:length(y)%%reduce == 0)]
            Mat <- mmand::resample(x = Mat, points = list(x = x[Which.x],
                                                          y = y[Which.y]),
                                   kernel = mmand::boxKernel())
            x <- x[Which.x]
            y <- y[Which.y]
        }

        # Calculate breaks
        N <- length(Col)
        Zero <- 0
        Max <- max(Mat)
        Min <- min(Mat)
        
        if (N%%2 == 0){
            Breaks <- c(seq(Min, Zero, length.out = round(N / 2, 0) + 1),
                        seq(Zero, Max, length.out = round(N / 2, 0) + 1)[2:(round(N / 2, 0) + 1)])
        } else {
            Breaks <- c(seq(Min, Zero, length.out = round(N / 2, 0) + 2)[1:(round(N / 2, 0) + 1)],
                        seq(Zero, Max, length.out = round(N / 2, 0) + 2)[2:(round(N / 2, 0) + 2)])
        }

        if (is.null(zlim)){
            zlim <- range(Mat, na.rm = TRUE)
        }
        
        # plot 3D surface and maybe 2D projection of it
        if (projection == TRUE){
            WW <- fields::drape.plot(x = x, y = y, z = Mat, col = Col, breaks = Breaks, zlim = zlim, ...)
            COL <- fields::drape.color(z = Mat, col = Col, zlim = zlim, breaks = Breaks)$color.index
            for (i in 2:NROW(Mat)) {
                for (j in 2:NCOL(Mat)) {
                    Points <- grDevices::trans3d(y = y[c(j - 1, j, j, j - 1, j - 1)],
                                                 x = x[c(i - 1, i - 1, i, i, i - 1)],
                                                 z = rep(zlim[1], length(5)), pmat = WW)
                    polygon(Points$x, Points$y,
                            border = NA, col = COL[i - 1, j - 1])
                }
            }
            par(new=T)
            fields::drape.plot(x = x, y = y, z = Mat, col = Col, breaks = Breaks, zlim = zlim, ...)
        } else {
            WW <- fields::drape.plot(x = x, y = y, z = Mat, col = Col, breaks = Breaks, zlim = zlim, ...)
        }
        
        # add spectra to x- and/or y-axis
        if (!is.null(specx)) {
            if (is.null(scalex)) {
                scalex <- 1
            }
            X <- seq(min(x), max(x), length.out = length(specx))
            Points.x <- grDevices::trans3d(x = X, y = min(y) + scalex * specx,
                                           z = rep(zlim[1], length(X)), pmat = WW)
            lines(x = Points.x$x, Points.x$y, lwd = 2)
        }
        
        if (!is.null(specy)) {
            if (is.null(scaley)) {
                scaley <- 1
            }
            Y <- seq(min(y), max(y), length.out = length(specy))
            Points.y <- grDevices::trans3d(y = Y, x = max(x) - scaley * specy,
                                           z = rep(zlim[1], length(Y)), pmat = WW)
            lines(x = Points.y$x, Points.y$y, lwd = 2)
        }
        
        on.exit(options(par(par_old)), add = TRUE)
    }




#' @seealso See \code{plot\_corr2d()} for further information.
#' 
#' @method plot corr2d
#' @export
plot.corr2d <- function(x, ...)
{
    plot_corr2d(x, ...)
}
