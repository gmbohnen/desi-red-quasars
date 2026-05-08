from .parameters import BAL_SIGMA_MASK, BAL_MIN_WIDTH, CIV_FIT_WINDOW
import numpy as np


def make_bal_mask(lam, flux_sub, ivar=None, sigma_thresh=BAL_SIGMA_MASK, min_width_aa=BAL_MIN_WIDTH, window=CIV_FIT_WINDOW):
    '''Takes continuum subtracted flux, checks for broad absorption lines (BALs), and returns a mask of BAL detected regions.'''

    win = (lam >= window[0]) & (lam <= window[1])
    bal_mask = np.zeros(len(lam), dtype=bool)

    if ivar is not None and np.any(ivar[win] > 0):
        noise = np.where(ivar > 0, 1.0 / np.sqrt(np.maximum(ivar, 1e-30)), np.inf)  # basically gets the variance back out of the invariance and use that as noise estimate
    else:
        # TODO: check if that is a good way to go
        noise_est = np.percentile(np.abs(flux_sub[win]), 20)
        noise = np.full(len(lam), max(noise_est, 1e-30))

    below   = flux_sub < -sigma_thresh * noise
    in_run  = False
    run_start = 0

    # check if the index is below the threshold (i.e. absorption happens), if it is, start a run (i.e. an absorption section)
    # continue until next index that is not below, then check if section is long (broad) enough to be a BAL, if yes add to bal_mask, of not, continue
    for i in range(len(lam)):
        if below[i] and not in_run:
            in_run = True
            run_start = i
        elif not below[i] and in_run:
            in_run = False
            if lam[i - 1] - lam[run_start] >= min_width_aa:
                bal_mask[run_start:i] = True

    if in_run and lam[-1] - lam[run_start] >= min_width_aa:
        bal_mask[run_start:] = True
 
    return bal_mask