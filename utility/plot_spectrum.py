import matplotlib.pyplot as plt
from astropy.convolution import convolve, Gaussian1DKernel
from colors import pink


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


def plot_spectrum(lam, flux, title="DESI-DR1 galaxy spectrum", add_smoothed=True, color=pink):
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
    plt.xlabel(r"$\lambda$ [$\AA$]")
    plt.ylabel(r"$F_{\lambda}~[10^{-17}~ergs~s^{-1}~cm^{-2}~{\AA}^{-1}]$")
    plt.legend(frameon=False)

    plt.show()