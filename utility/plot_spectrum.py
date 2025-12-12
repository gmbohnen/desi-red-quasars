import matplotlib.pyplot as plt
from astropy.convolution import convolve, Gaussian1DKernel


settings = {
        'font.size':16,
        'xtick.direction':'in', 
        'xtick.minor.visible':True,
        'xtick.top':True,
        'ytick.direction':'in', 
        'ytick.minor.visible':True,
        'ytick.right':True
    }

color = "#f86588"


def plot_spectrum(lam, flux, title="DESI-DR1 galaxy spectrum", add_smoothed=True):
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

    plt.rcParams.update(**settings)

    plt.subplots(figsize=(15,8))

    plt.plot(lam, flux, linewidth=0.4, color=color, label="raw", alpha=0.7)
    
    if add_smoothed:
        plt.plot(lam, convolve(flux, Gaussian1DKernel(5)), linewidth=0.8, color="black", label="smoothed")

    plt.title(title)
    plt.xlabel(r"$\lambda$ [nm]")
    plt.ylabel(r"$F_{\lambda}~[10^{-17}~ergs~s^{-1}~cm^{-2}~{\AA}^{-1}]$")
    plt.legend(loc="lower right")

    plt.show()