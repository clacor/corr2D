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