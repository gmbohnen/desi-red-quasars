import matplotlib.pyplot as plt
from astropy.convolution import convolve, Gaussian1DKernel
from .my_colors import pink, palette_dark


settings = {
        'font.size':16,
        'xtick.direction':'in', 
        'xtick.minor.visible':True,
        'xtick.top':True,
        'ytick.direction':'in', 
        'ytick.minor.visible':True,
        'ytick.right':True,
        'xtick.major.size':6.0,
        'xtick.minor.size':4.0,
        'ytick.major.size':6.0,
        'ytick.minor.size':4.0,
        'legend.frameon':False
    }


def plot_spectrum(lam, flux, title="", flux_label="raw", add_smoothed=True, line_lambda=None, line_label=None, additional_lam=None, additional_flux=None, additional_flux_label=None, additional_flux_color=palette_dark, line_width_factor=1):
    '''
    Plot spectrum, optionally smoothed spectrum on top of it.

    Parameters
    ----------
    lam : list of float
        Wavelength data.
    flux : list of float
        Flux data.
    title : str, optional
        Specify plot title.
    flux_label : str, optional
        Specify label of primary flux.
    add_smoothed : bool, optional
        Specifies whether smoothed spectrum is added on top of primary spectrum.
    line_lambda : list of float or None, optional
        Passing values adds vertical lines.
    line_label : list of float or None, optional
        Specify labels of lines.
    additional_lam : list of list of float or None, optional
        Passing values adds a secondary spectrum on top of the other. Cannot be used with add_smoothed=True. Can be None if it is equivalent to primary lam.
    additional_flux : list of list of float or None, optional
        Passing values adds a secondary spectrum on top of the other. Cannot be used with add_smoothed=True. Uses lam if additional_lam is not specified.
    additional_flux_label : list of str or None, optional
        Specify label of additional flux.
    additional_flux_color : list of str, optional
        Specify colors of additional flux.
    color : str, optional
        Specify color of primary spectrum.
    line_width_factor : int, optional
        Scales line widths by that factor.    
    '''

    assert additional_flux is None or not add_smoothed, "Can only have one of the two: additional_flux, add_smoothed"

    plt.rcParams.update(**settings)  # update settings to apply styles

    plt.subplots(figsize=(15,8))  

    # plot primary flux
    plt.plot(lam, flux, linewidth=0.4*line_width_factor, color=pink, label=flux_label, alpha=0.7)
    
    # plot smoothed version of primary flux if specified
    if add_smoothed:
        plt.plot(lam, convolve(flux, Gaussian1DKernel(5)), linewidth=0.8*line_width_factor, color="black", label="smoothed")   # 5 is value from datalab tutorial

    # plot additional fluxes if specified
    if additional_flux is not None:
        if additional_lam is not None:
            for i in range(len(additional_flux)):
                plt.plot(additional_lam[i], additional_flux[i], linewidth=0.8*line_width_factor, color=additional_flux_color[i], label=additional_flux_label[i])
        else:
            for i in range(len(additional_flux)):
                plt.plot(lam, additional_flux[i], linewidth=0.8*line_width_factor, color=additional_flux_color[i], label=additional_flux_label[i])

    # plot vertical lines if specified
    if line_lambda:
        for i, elem in enumerate(line_lambda):
            plt.axvline(line_lambda[i], label=line_label[i], color=palette_dark[i+1])


    plt.title(title)
    plt.xlabel(r"$\lambda$ [$\AA$]")
    plt.ylabel(r"$F_{\lambda}~[10^{-17}~ergs~s^{-1}~cm^{-2}~{\AA}^{-1}]$")
    leg = plt.legend()

    # increase and align line width in legend
    for legobj in leg.legend_handles:
        legobj.set_linewidth(3)

    plt.show()

    plt.rcdefaults()