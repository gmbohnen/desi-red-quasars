import matplotlib.pyplot as plt
from astropy.convolution import convolve, Gaussian1DKernel
from .my_colors import pink


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
        'ytick.minor.size':4.0
    }


def plot_spectrum(lam, flux, title="DESI-DR1 galaxy spectrum", flux_label="raw", line_lambda=None, line_label=None, add_smoothed=True, color=pink, additional_flux=None, additional_flux_label=None, line_width_factor=1):
    '''
    Plot spectrum, optionally smoothed spectrum on top of it.

    Parameters
    ----------
    lam : list of float
        Wavelength data
    flux : list of float
        Flux data
    title : str or None, optional
        Specify plot title
    add_smoothed : bool, optional
        Specifies whether smoothed spectrum is added on top of raw spectrum.
    '''

    assert additional_flux is None or not add_smoothed, "Can only have one of the two: additional_flux, add_smoothed"

    plt.rcParams.update(**settings)

    plt.subplots(figsize=(15,8))

    plt.plot(lam, flux, linewidth=0.4*line_width_factor, color=color, label=flux_label, alpha=0.7)
    
    if add_smoothed:
        plt.plot(lam, convolve(flux, Gaussian1DKernel(5)), linewidth=0.8*line_width_factor, color="black", label="smoothed")   # 5 is value from datalab tutorial

    if additional_flux is not None:
        plt.plot(lam, additional_flux, linewidth=0.8*line_width_factor, color="black", label=additional_flux_label)

    if line_lambda:
        plt.axvline(line_lambda, label=line_label)

    plt.title(title)
    plt.xlabel(r"$\lambda$ [$\AA$]")
    plt.ylabel(r"$F_{\lambda}~[10^{-17}~ergs~s^{-1}~cm^{-2}~{\AA}^{-1}]$")
    plt.legend(frameon=False)

    plt.show()

    plt.rcdefaults()