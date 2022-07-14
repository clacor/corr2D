# corr2D 1.0.3
- Moved package profr and xtable to Suggests. They are optionally used to check
    corr2ds performance and are thus not essential to the
    package.

# corr2D 1.0.1
- Moved package rgl to Suggests as the future of OpenGL is uncertain. rgl is
    used as an optional preview in corr2d and is thus not essential to the
    package.

# corr2D 1.0.0
- Added publication in Journal of Statistical Software as preferred
    citation. The publication is a good start for newcomers to the field of
    2D correlation spectroscopy (in R)!

# corr2D 0.4.0
- Added internal 'testthat' tests to corr2D to easier check correlation results
    during package development.
- Added a default print() method for corr2d objects (plot.default() for now).
- Changed plot_corr2d() to use plot() instead of plot.default() to enable
    method dispatch.
- Fixed a small error with the calculation of reference spectrum 2 (Ref2) if
    reference spectrum 1 (Ref1) is not present (Thanks to Bettina Gruen).

# corr2D 0.3.0
- Added 2T2D correlation analysis to the package. Function uses the approach
    described by I. Noda (2018)
    <DOI: https://doi.org/10.1016/j.molstruc.2018.01.091>.
- A freshly published paper at The Journal of Statistical Software was added
    as a vignette. This will be updated in the future to also include
    descriptions of new functions.
- Added citation referencing the published paper at The Journal of Statistical
    Software.

# corr2D 0.2.0
- Added 2D codistibution analysis to the packages functions. Function uses
    the approach described by I. Noda (2014)
    <DOI: https://doi.org/10.1016/j.molstruc.2014.01.024>.
- Rewrote the function plot_corr2d() for better control about the plot
    appearence. Introduced the graphical parameters "lwd", "lwd.axis",
    "lwd.spec", "col", "col.axis", "col.lab", "cex.axis", "cex.lab", "cex.leg",
    "font.axis" and "font.lab" which are/are derived from "par".
    "at.xaxs"/"at.yaxs" and "label.xaxs"/"label.yaxs" allow control over the
    axes ticks and thier labels. "line.xlab"/"line.ylab" control the position
    of the axes label.

# corr2D 0.1.12
- Rewrote the function plot_corr2d(): The graphical parameters specified at
    "..." are now partly transferred to all parts of the plot, not just the
    main screen.
- Cleaned up the image.plot() code inside plot_corr2d().
- Fixed a bug which prevented the use of an individual color palette in
    plot_corr2din3d().
- Changed the default color palette in plot_corr2din3d() from
    fields::timcolors() to colorspace::diverge_hcl(). This change should help
    to improve the default graphics quality.

# corr2D 0.1.11
- Corrected the default value for the normalization factor. It's now the number
    of perturbation variables which is NROW() of Mat1.
- Added error massages to all functions.
- Reworte the function sim2ddata() to allow to set C and Camp NULL. In that
    case only the first order reaction A -> B will be simulated and sampled.

# corr2D 0.1.10
- Added the function sim2ddata to simulate artificial data.
- Fixed the plotting of the reference spectra on the 2D correlation spectra.
    Before the reference spectra used an arbitrary x-axis. Now they use their
    correct spectral variables for plotting. 2D correlation peaks should now
    align with their 1D counterparts in the reference spectra.
- Rewrote the interpolation of the perturbation variables to no longer use the
    ineffective cbind in a for loop. Instead apply and sapply are now used
    instead.
- The interpolated perturbation variables are now saved when the perturbation
    values get interpolated.
- Cleaned up some left/right confusions in the documentation of plot_corr2d.

# corr2D 0.1.9
- Package is released on CRAN

## Known Bugs
- A new plot after a call of plot_corr2d() is plotted inside the main part.
    It's a consequence of keeping the split.screen of the main part active
    for interactive data readout and manipulation.