from scipy.interpolate import interp1d
import numpy as np
import warnings


def remove_spikes(lam, flux, ivar=None, sigma_thresh=5.0, window=11):
    '''
    Detect and remove narrow spikes (cosmic rays / noise anomalies).

    A pixel is flagged as a spike if its flux deviates from the local median
    by more than `sigma_thresh` * local MAD (= median absolute deviation). Flagged pixels are replaced by
    linear interpolation from their neighbours.

    Parameters
    ----------
    lam : array_like
        Wavelength array (Å).
    flux : array_like
        Flux array (arbitrary units).
    ivar : array_like or None
        Inverse variance array. If provided, zero-ivar pixels are also masked.
    sigma_thresh : float
        Number of sigma above/below the local median to flag a spike.
    window : int
        Half-width (in pixels) of the local median filter window.

    Returns
    -------
    flux_clean : ndarray
        Flux array with spikes replaced by interpolated values.
    spike_mask : ndarray of bool
        True where a spike was detected and replaced.
    '''

    flux = np.array(flux, dtype=float)
    lam = np.array(lam, dtype=float)
    n = len(flux)
    spike_mask = np.zeros(n, dtype=bool)

    # Flag zero/negative ivar pixels
    if ivar is not None:
        spike_mask |= (np.array(ivar) <= 0)

    # Rolling median and MAD-based sigma clipping
    for i in range(n):
        lo = max(0, i - window)
        hi = min(n, i + window + 1)
        neighbours = np.concatenate([flux[lo:i], flux[i+1:hi]])

        if len(neighbours) < 3:
            continue

        med = np.median(neighbours)
        mad = np.median(np.abs(neighbours - med))

        sigma = 1.4826 * mad  # MAD*k is an estimator for sigma, where k is a scaling factor depending on the dstribution, 1.4826 for gaussian

        if sigma > 0 and np.abs(flux[i] - med) > sigma_thresh * sigma:
            spike_mask[i] = True

    # Interpolate over flagged pixels using clean neighbours
    flux_clean = flux.copy()

    if spike_mask.any():
        good = ~spike_mask

        if good.sum() < 2:  # TODO: 2 is only the technical limit, might wanna raise this lol
            warnings.warn("Too many spikes — cannot interpolate reliably.")
            return flux_clean, spike_mask

        interp_fn = interp1d(lam[good], flux[good], kind='linear', bounds_error=False, fill_value='extrapolate')
        
        flux_clean[spike_mask] = interp_fn(lam[spike_mask])

    return flux_clean, spike_mask