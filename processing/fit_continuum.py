from .parameters import CONTINUUM_WINDOWS
from lmfit import Model
import numpy as np


def _powerlaw(x, amplitude, alpha):
    '''F(λ) = amplitude * (λ / λ_ref) ^ alpha'''
    return amplitude * (x ** alpha)


def fit_continuum(lam, flux, ivar=None,
                  windows=CONTINUUM_WINDOWS,
                  lambda_ref=1700.0):  # TODO whatever is a reasonable value for that
    '''
    Fit a power-law continuum using median fluxes in emission-free windows.

    Parameters
    ----------
    lam : array_like
        Rest-frame wavelength array (Å).
    flux : array_like
        Flux array (spike-cleaned).
    ivar : array_like or None
        Inverse variance. Used to weight the window medians if provided.
    windows : list of (float, float)
        Rest-frame wavelength windows free of emission lines.
    lambda_ref : float
        Reference wavelength for the power law (normalisation pivot).

    Returns
    -------
    result : lmfit ModelResult
        Full lmfit fit result (contains best-fit params, uncertainties, etc.)
    continuum : ndarray
        Best-fit continuum evaluated at every wavelength in `lam`.
    window_mask : ndarray of bool
        True for pixels used in the fit.
    '''
    lam = np.array(lam, dtype=float)
    flux = np.array(flux, dtype=float)

    # Build window mask and compute per-window median flux
    window_lams = []
    window_fluxes = []
    window_weights = []

    for w0, w1 in windows:
        mask = (lam >= w0) & (lam <= w1) & np.isfinite(flux)
        if mask.sum() < 3:
            continue
        wc = 0.5 * (w0 + w1)          # window centre wavelength
        med_flux = np.median(flux[mask])
        if ivar is not None:
            med_ivar = np.median(ivar[mask])
            weight = np.sqrt(max(med_ivar, 0))
        else:
            weight = 1.0
        window_lams.append(wc)
        window_fluxes.append(med_flux)
        window_weights.append(weight)

    if len(window_lams) < 2:
        raise ValueError("Fewer than 2 continuum windows have data; cannot fit continuum.")

    xdata = np.array(window_lams) / lambda_ref
    ydata = np.array(window_fluxes)
    weights = np.array(window_weights)
    if weights.sum() == 0:
        weights = np.ones_like(weights)

    # lmfit power-law model
    plaw_model = Model(_powerlaw)
    params = plaw_model.make_params(
        amplitude=dict(value=np.median(ydata), min=0),
        alpha=dict(value=-1.5, min=-10, max=10),
    )

    result = plaw_model.fit(ydata, params, x=xdata, weights=weights)

    # Evaluate continuum across the full spectrum
    continuum = result.best_values['amplitude'] * ((lam / lambda_ref) ** result.best_values['alpha'])

    # Full-spectrum window mask (for plotting)
    window_mask = np.zeros(len(lam), dtype=bool)
    for w0, w1 in windows:
        window_mask |= (lam >= w0) & (lam <= w1)

    return result, continuum, window_mask   