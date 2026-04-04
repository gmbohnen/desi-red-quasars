from .parameters import BAL_SIGMA_MASK, BAL_MIN_WIDTH, FIT_WINDOW
import numpy as np

def make_bal_mask(wave, flux_sub, ivar=None,
                  sigma_thresh=BAL_SIGMA_MASK,
                  min_width_aa=BAL_MIN_WIDTH,
                  window=FIT_WINDOW):
    win      = (wave >= window[0]) & (wave <= window[1])
    bal_mask = np.zeros(len(wave), dtype=bool)

    if ivar is not None and np.any(ivar[win] > 0):
        noise = np.where(ivar > 0, 1.0 / np.sqrt(np.maximum(ivar, 1e-30)), np.inf)
    else:
        noise_est = np.percentile(np.abs(flux_sub[win]), 10)
        noise     = np.full(len(wave), max(noise_est, 1e-30))

    below   = flux_sub < -sigma_thresh * noise
    in_run  = False
    run_start = 0

    for i in range(len(wave)):
        if below[i] and not in_run:
            in_run    = True
            run_start = i
        elif not below[i] and in_run:
            in_run = False
            if wave[i - 1] - wave[run_start] >= min_width_aa:
                bal_mask[run_start:i] = True
    if in_run and wave[-1] - wave[run_start] >= min_width_aa:
        bal_mask[run_start:] = True
 
    return bal_mask