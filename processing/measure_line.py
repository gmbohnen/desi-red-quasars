from .parameters import C_KMS, CIV_VACUUM, FIT_WINDOW
import numpy as np


def _profile_width(wave, profile, level):
    above     = profile >= level
    crossings = np.where(np.diff(above.astype(int)))[0]
    if len(crossings) < 2:
        return np.nan
    i         = crossings[0]
    lam_left  = np.interp(level, [profile[i], profile[i+1]], [wave[i], wave[i+1]])
    j         = crossings[-1]
    lam_right = np.interp(level, [profile[j+1], profile[j]], [wave[j+1], wave[j]])
    return lam_right - lam_left


def _kt80(wave, profile):
    '''
    Kurtosis index:
    kt80 = Dv(80%) / Dv(20%)
    where Dv(X%) is the velocity width of the profile at X% of peak height.
    Equivalent to the ratio of line widths at 80% and 20% of the peak flux.
    '''
    peak = profile.max()
    if peak <= 0:
        return np.nan

    dv_80 = _profile_width(wave, profile, 0.80 * peak)
    dv_20 = _profile_width(wave, profile, 0.20 * peak)

    # Convert widths from Angstrom to km/s
    dv_80_kms = C_KMS * dv_80 / CIV_VACUUM
    dv_20_kms = C_KMS * dv_20 / CIV_VACUUM

    if np.isnan(dv_80_kms) or np.isnan(dv_20_kms) or dv_80_kms <= 0:
        return np.nan

    return dv_80_kms / dv_20_kms


def wave_to_vel(wave, wave0=CIV_VACUUM):
    return C_KMS * (wave - wave0) / wave0


def measure_line(wave, profile, flux_sub, continuum_flux, window=FIT_WINDOW):
    win = (wave >= window[0]) & (wave <= window[1]) & (profile > 0)
    if win.sum() < 5:
        return {k: np.nan for k in
                ['fwhm_kms','ew_aa','line_flux','centroid_kms','kt80']}

    w  = wave[win]
    p  = profile[win]
    fs = flux_sub[win]
    fc = continuum_flux[win] if continuum_flux is not None else np.ones(win.sum())

    fwhm_aa      = _profile_width(w, p, 0.5 * p.max())
    fwhm_kms     = C_KMS * fwhm_aa / CIV_VACUUM
    centroid_aa  = np.trapz(w * p, w) / np.trapz(p, w)
    centroid_kms = wave_to_vel(centroid_aa)

    line_win  = win & (profile > 0.01 * profile[win].max())
    line_flux = np.trapz(flux_sub[line_win], wave[line_win])

    with np.errstate(divide='ignore', invalid='ignore'):
        ew_integrand = np.where(fc > 0, fs / fc, 0.0)
    ew_aa = np.trapz(ew_integrand, w)

    kt80 = _kt80(w, p)

    return dict(fwhm_kms=fwhm_kms, ew_aa=ew_aa, line_flux=line_flux,
                centroid_kms=centroid_kms, kt80=kt80)