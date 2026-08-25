from .parameters import LYA_NV_FIT_WINDOW, NV_AIR, CIV_AIR, C_KMS
from .fit_gaussians import _gaussian, _two_gaussian

from lmfit import Model, Parameters
import numpy as np


def _civ_template(x,scale,shift,civ_params,two_gaussians):
    # idea: keep CIV parameters as they are, i.e. dont shift them, but multiply the flux values around NV by the lam_shift to see what their values would be in the CIV model
    # lam_shift: ~0.8, i.e. NV = 0.8*CIV, so a feature that is a certain point in the CIV line, is at 0.8 times the wavelength within NV
    # actual parameter that will be fit is an additional shift on top of the inherent one due to the lambda difference
    lam_shift = NV_AIR / CIV_AIR 
    dv_factor = 1 + (shift / C_KMS)

    # rescaled wavelenghts
    w = x / (lam_shift * dv_factor)

    if two_gaussians:
        gaussian = _two_gaussian(w,civ_params["amp_c"],civ_params["cen_c"],civ_params["sigma_c"],civ_params["amp_w"],civ_params["cen_w"],civ_params["sigma_w"])
    else:
        gaussian = _gaussian(w,civ_params["amplitude"],civ_params["center"],civ_params["sigma"])
    
    return scale * gaussian


def _make_nv_model(civ_params,two_gaussians):
    def nv_model(x, scale, shift):
        return _civ_template(x, scale, shift, civ_params, two_gaussians)
    return nv_model


def scale_civ_to_nv(lam,flux_sub,params_dict,accept_two,ivar=None,mask=None,fit_window=LYA_NV_FIT_WINDOW,verbose=True):
        win = (lam > fit_window[0]) & (lam < fit_window[1])

        if mask is not None:
            win &= ~mask
        if win.sum() < 8:
            if verbose:
                raise ValueError("Too few pixels for scaling.")
            else:
                return None, None

        xd = lam[win]
        yd = flux_sub[win]
        weights = np.sqrt(ivar[win]) if ivar is not None else np.ones(win.sum())


        model = Model(_make_nv_model(params_dict,accept_two),independent_vars=["x"])

        params = Parameters()
        params.add('scale',min=10E-2)
        params.add('shift',min=-500,max=500)

        result = model.fit(yd,params,x=xd,weights=weights)

        bv = result.best_values

        profile = _civ_template(lam,bv["scale"],bv["shift"],params_dict,accept_two)

        return result, profile