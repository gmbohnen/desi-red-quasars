from .parameters import FIT_WINDOW, MAX_WING_SIGMA, MIN_CORE_SIGMA, MAX_CORE_SIGMA, CIV_VACUUM, MIN_WIDTH_RATIO, MAX_WIDTH_RATIO, MAX_WING_SIGMA, FTEST_PVALUE
import numpy as np
from lmfit import Model, Parameters
from scipy.stats import f as f_dist

## single gaussian

def _gaussian(x, amplitude, center, sigma):
    return amplitude * np.exp(-0.5 * ((x - center) / sigma) ** 2)


def fit_single_gaussian(wave, flux_sub, ivar=None, mask=None,
                         window=FIT_WINDOW):
    win = (wave >= window[0]) & (wave <= window[1])
    if mask is not None:
        win &= ~mask
    if win.sum() < 5:
        raise ValueError("Too few unmasked pixels in CIV window.")

    xd      = wave[win]
    yd      = flux_sub[win]
    weights = np.sqrt(ivar[win]) if ivar is not None else np.ones(win.sum())
    weights = np.where(np.isfinite(weights) & (weights > 0), weights, 1.0)

    gmodel = Model(_gaussian)
    params = gmodel.make_params(
        amplitude=dict(value=max(yd.max(), 1e-30), min=0),
        center   =dict(value=xd[np.argmax(yd)], min=window[0], max=window[1]),
        sigma    =dict(value=20.0, min=MIN_CORE_SIGMA, max=MAX_CORE_SIGMA),
    )
    result  = gmodel.fit(yd, params, x=xd, weights=weights, method='least_squares')
    profile = _gaussian(wave, result.best_values['amplitude'],
                        result.best_values['center'],
                        result.best_values['sigma'])
    return result, profile


## two gaussians

def _two_gaussian(x, amp_c, cen_c, sigma_c, amp_w, cen_w, sigma_w):
    core = _gaussian(x, amp_c, cen_c, sigma_c)
    wing = _gaussian(x, amp_w, cen_w, sigma_w)
    return core + wing


def fit_two_gaussians(wave, flux_sub, ivar=None, mask=None,
                      single_result=None, window=FIT_WINDOW):
    win = (wave >= window[0]) & (wave <= window[1])
    if mask is not None:
        win &= ~mask
    if win.sum() < 8:
        raise ValueError("Too few pixels for two-Gaussian fit.")

    xd      = wave[win]
    yd      = flux_sub[win]
    weights = np.sqrt(ivar[win]) if ivar is not None else np.ones(win.sum())
    weights = np.where(np.isfinite(weights) & (weights > 0), weights, 1.0)

    if single_result is not None:
        bv   = single_result.best_values
        amp0 = bv['amplitude']
        cen0 = bv['center']
        sig0 = bv['sigma']
    else:
        amp0 = max(yd.max(), 1e-30)
        cen0 = CIV_VACUUM
        sig0 = 20.0

    sig_w0 = np.clip(sig0 * 2.0, sig0 * MIN_WIDTH_RATIO, MAX_WING_SIGMA)
    hw     = 1.1775 * sig_w0   # half-FWHM of wing for core centroid bound

    tmodel = Model(_two_gaussian)
    params = Parameters()

    # Wing
    params.add('amp_w',      value=amp0 * 0.5, min=0)
    params.add('cen_w',      value=cen0, min=window[0], max=window[1])
    params.add('sigma_c',    value=sig0,  min=MIN_CORE_SIGMA, max=MAX_CORE_SIGMA)
    params.add('width_ratio',value=2.0,   min=MIN_WIDTH_RATIO, max=MAX_WIDTH_RATIO)
    params.add('sigma_w',    expr='width_ratio * sigma_c')

    # Core
    params.add('amp_c',      value=amp0 * 0.7, min=0)
    params.add('delta_cen',  value=0.0, min=-hw, max=hw)
    params.add('cen_c',      expr='cen_w + delta_cen')

    result  = tmodel.fit(yd, params, x=xd, weights=weights, method='least_squares')
    bv      = result.best_values
    profile = _two_gaussian(wave,
                             bv['amp_c'], bv['cen_c'], bv['sigma_c'],
                             bv['amp_w'], bv['cen_w'], bv['sigma_w'])
    return result, profile


def ftest_which_gaussian(result1, result2, alpha=FTEST_PVALUE):
    '''F-test whether to choose the single or double gaussian line fit.'''
    rss1 = result1.residual @ result1.residual
    rss2 = result2.residual @ result2.residual
    n    = len(result1.residual)
    p1   = result1.nvarys
    p2   = result2.nvarys

    if rss2 <= 0 or rss1 <= rss2:
        return False, 1.0

    f_stat = ((rss1 - rss2) / (p2 - p1)) / (rss2 / (n - p2))
    p_val  = 1.0 - f_dist.cdf(f_stat, dfn=(p2 - p1), dfd=(n - p2))
    return p_val < alpha, p_val