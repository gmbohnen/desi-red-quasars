from .parameters import C_KMS, CIV_AIR, CIV_FIT_WINDOW
import numpy as np


def _profile_width(lam, profile, level):
    above     = profile >= level
    crossings = np.where(np.diff(above.astype(int)))[0]
    if len(crossings) < 2:
        return np.nan

    i         = crossings[0]
    lam_left  = np.interp(level, [profile[i], profile[i+1]], [lam[i], lam[i+1]])
    j         = crossings[-1]
    lam_right = np.interp(level, [profile[j+1], profile[j]], [lam[j+1], lam[j]])

    return lam_right - lam_left


def _kt80(lam, profile):
    '''
    Kurtosis index:
    kt80 = Dv(80%) / Dv(20%)
    where Dv(X%) is the velocity width of the profile at X% of peak height.
    Equivalent to the ratio of line widths at 80% and 20% of the peak flux.
    '''

    peak = profile.max()
    if peak <= 0:
        return np.nan

    dv_80 = _profile_width(lam, profile, 0.80 * peak)
    dv_20 = _profile_width(lam, profile, 0.20 * peak)

    # convert to km/s
    dv_80_kms = C_KMS * dv_80 / CIV_AIR
    dv_20_kms = C_KMS * dv_20 / CIV_AIR

    if np.isnan(dv_80_kms) or np.isnan(dv_20_kms) or dv_80_kms <= 0:
        return np.nan

    return dv_80_kms / dv_20_kms


def lam_to_vel(lam, lam0=CIV_AIR):
    return C_KMS * (lam - lam0) / lam0

def measure_line(lam, profile, flux_sub, continuum_flux, ci, window=CIV_FIT_WINDOW, get_window=False):
    # integration window constraints
    # 1. in the window (1450,1650)
    # 2. where the line fit profile is larger than 10E-4
    # 3. inside the given confidence interval
    win = (lam >= window[0]) & (lam <= window[1]) & (profile > 10E-4) & (lam >= ci[0]) & (lam <= ci[1])
    if win.sum() < 5:
        return {k: np.nan for k in ['fwhm_kms','ew_aa','line_flux','centroid_kms','kt80']}

    w  = lam[win]
    p  = profile[win]
    fs = flux_sub[win]
    fc = continuum_flux[win] if continuum_flux is not None else np.ones(win.sum())

    fwhm_aa      = _profile_width(w, p, 0.5 * p.max())
    fwhm_kms     = C_KMS * fwhm_aa / CIV_AIR

    centroid_aa  = np.trapz(w * p, w) / np.trapz(p, w)
    centroid_kms = lam_to_vel(centroid_aa)

    line_win  = win & (profile > 0.01 * profile[win].max())
    line_flux = np.trapz(flux_sub[line_win], lam[line_win])

    with np.errstate(divide='ignore', invalid='ignore'):
        ew_integrand = np.where(fc > 0, fs / fc, 0.0)
        fit_integrand = np.where(p > 0, p / fc, 0.0)
    ew_aa = np.trapz(ew_integrand, w)
    fit_rew = np.trapz(fit_integrand, w)

    kt80 = _kt80(w, p)

    if get_window:
        return dict(fwhm_kms=fwhm_kms, ew_aa=ew_aa, fit_rew=fit_rew, line_flux=line_flux, centroid_kms=centroid_kms, kt80=kt80), win
    else:
        return dict(fwhm_kms=fwhm_kms, ew_aa=ew_aa, fit_rew=fit_rew, line_flux=line_flux, centroid_kms=centroid_kms, kt80=kt80)