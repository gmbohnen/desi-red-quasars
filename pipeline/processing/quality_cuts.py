from .parameters import CIV_CONTINUUM_WINDOWS, CIV_REGION, MAX_POWERLAW_SLOPE, BAL_SIGMA_THRESHOLD, MIN_EMISSION_SIGMA
import numpy as np


def _estimate_noise(residual, lam, window=CIV_CONTINUUM_WINDOWS):
    '''Estimate noise from the scatter in continuum window residuals.'''

    mask = np.zeros(len(lam), dtype=bool)

    for w0, w1 in window:
        mask |= (np.array(lam) >= w0) & (np.array(lam) <= w1)

    if mask.sum() < 5:
        return np.std(residual)

    return 1.4826 * np.median(np.abs(residual[mask] - np.median(residual[mask])))  # MAD as sigma estimator


def quality_cuts(lam, flux, continuum, fit_result,
                 civ_region=CIV_REGION,
                 max_slope=MAX_POWERLAW_SLOPE,
                 bal_sigma=BAL_SIGMA_THRESHOLD,
                 min_emission_sigma=MIN_EMISSION_SIGMA):
    '''
    Apply the three quality cuts from Hamann.

    Parameters
    ----------
    lam : array_like
        Rest-frame wavelength array.
    flux : array_like
        Spike-cleaned flux.
    continuum : array_like
        Best-fit continuum flux.
    fit_result : lmfit ModelResult
        Result from fit_continuum.
    civ_region : (float, float)
        Wavelength range of the CIV emission line.
    max_slope : float
        Maximum allowed |alpha| before the fit is considered unrealistic.
    bal_sigma : float
        Sigma threshold for broad absorption detection (negative value).
    min_emission_sigma : float
        Minimum sigma for a significant emission detection.

    Returns
    -------
    rejected : bool
        True if the spectrum fails any quality cut.
    reason : str or None
        Human-readable reason for rejection, or None if accepted.
    flags : dict
        Dict of individual flag values for inspection.
    '''
    
    alpha = fit_result.best_values['alpha']
    residual = np.array(flux) - np.array(continuum)
    noise_scalar = _estimate_noise(residual, lam)

    civ_mask = (np.array(lam) >= civ_region[0]) & \
               (np.array(lam) <= civ_region[1])

    flags = {}

    # Cut 1: Unrealistically steep continuum
    flags['steep_continuum'] = abs(alpha) > max_slope
    if flags['steep_continuum']:
        return True, f"Continuum too steep (alpha={alpha:.2f})", flags

    # Cut 2: Broad absorption at CIV wavelengths
    if civ_mask.sum() > 0 and noise_scalar > 0:
        med_residual_civ = np.median(residual[civ_mask])
        bal_detection = med_residual_civ / noise_scalar
        flags['bal_detected'] = bal_detection < bal_sigma
        if flags['bal_detected']:
            return True, \
                f"Broad absorption detected at CIV (residual S/N={bal_detection:.1f})", flags
    else:
        flags['bal_detected'] = False

    # Cut 3: No significant emission above continuum
    if civ_mask.sum() > 0 and noise_scalar > 0:
        peak_emission = np.max(residual[civ_mask]) / noise_scalar
        flags['no_emission'] = peak_emission < min_emission_sigma
        if flags['no_emission']:
            return True, \
                f"No significant CIV emission (peak S/N={peak_emission:.1f})", flags
    else:
        flags['no_emission'] = True
        return True, "CIV region not covered by spectrum", flags

    return False, None, flags