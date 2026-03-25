from utility.adjust_redshift import adjust_redshift

from processing.remove_spikes import remove_spikes
from processing.fit_continuum import fit_continuum
from processing.quality_cuts import quality_cuts


def pipeline(z, lam, flux, ivar=None):
    if ivar is not None:
        lam, flux, ivar = adjust_redshift(z, lam, flux, ivar=ivar)
    else:
        lam, flux = adjust_redshift(z, lam, flux)

    flux_clean, spike_mask = remove_spikes(lam,flux,ivar=ivar)

    result, continuum, window_mask = fit_continuum(lam, flux_clean, ivar=ivar)

    rejected, reason, flags = quality_cuts(lam, flux_clean, continuum, result)