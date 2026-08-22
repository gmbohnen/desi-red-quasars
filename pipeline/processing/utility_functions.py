from .parameters import CIV_AIR, NV_AIR, LYA_AIR

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


def make_scaling_mask(lam,profile_civ,lims=False):
    '''Creates mask that masks the left side of LYA to enable proper scaling of the CIV template to NV.
    Mask is +-60% of the midpoint between NV and LYA, adjusted for possible redshift inaccuracies by adding the difference between the theoretical CIV line and the peak of the actual fit.'''
    diff_civ = lam[np.argmax(profile_civ)] - CIV_AIR
    m = (NV_AIR-LYA_AIR)  # distance between NV and LYA
    c = LYA_AIR + diff_civ  # center of the window, adjusted for redshift inaccuracies using the difference between CIV_AIR and the actual profile peak
    lims = (c-0.6*m,c+0.6*m)
    lya_mask = (lam > lims[0]) & (lam < lims[1])

    if not lims:
        return lya_mask
    else:
        return lya_mask, lims