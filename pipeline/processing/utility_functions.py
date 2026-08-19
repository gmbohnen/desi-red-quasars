import numpy as np


def snr_around_lam(lam,flux,noise,ref_lam=1700):
    '''Calculates the signal-to-noise ratio in a +-5A window around the specified ref_lam as the median of the absolute pixel wise ratios.'''
    mask = (lam > ref_lam-5) & (lam < ref_lam+5)
    snr = np.median(np.abs(flux[mask] / noise[mask]))

    return snr


def pixel_used_in_window(lam,mask,lower_lim,upper_lim):
    '''Calculates how many pixels within the interval [lower_lim,upper_lim] are not masked away by the mask.'''
    util_mask = (lam >= lower_lim) & (lam <= upper_lim)
    pixel_used = (util_mask & ~mask).sum() / util_mask.sum()
    
    return pixel_used