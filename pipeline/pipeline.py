from .processing.remove_spikes import remove_spikes
from .processing.fit_continuum import fit_continuum
from .processing.get_absorption_intervals import get_absorption_intervals
from .processing.fit_gaussians import fit_single_gaussian, fit_two_gaussians, ftest_which_gaussian
from .processing.measure_line import measure_line
from .processing.scale_civ_to_nv import scale_civ_to_nv
from .processing.utility_functions import *
from .processing.parameters import CIV_CONTINUUM_WINDOWS, CIV_FIT_WINDOW, LYA_NV_FIT_WINDOW, NV_LYA_CONTINUUM_WINDOWS, CIV_AIR, NV_AIR, LYA_AIR

import numpy as np
import h5py
import warnings


spectra_hdf5_path = "data/spectra.h5"  # TODO move to a config file?


def pipeline(idx,verbose=False):

    if not verbose:
        warnings.filterwarnings('ignore')


    with h5py.File(spectra_hdf5_path, "r") as f:
        data = f[f"{idx}"][()]
        
    lam = data[0,:]
    flux = data[1,:]
    ivar = data[2,:]

    # Note: preprocessing (redshift removal and bad pixel masking) already done in retrieval pipeline

    # create noise
    noise = np.abs(1/np.sqrt(ivar))

    # remove spikes
    flux_clean, _ = remove_spikes(lam=lam,flux=flux,ivar=ivar)


    #------------------------------------------------------------------------------------

    ## FIT CONTINUUM BELOW CIV
    result_civ_cont, continuum_civ, _ = fit_continuum(lam=lam, flux=flux_clean, ivar=ivar, windows=CIV_CONTINUUM_WINDOWS, lambda_ref=1700.0)
      

    ## FIT CIV LINE

    # subtract CIV continuum
    flux_sub_civ = flux_clean - continuum_civ

    # get absorption mask
    civ_absorption_intervals, civ_absorption_mask, civ_blue_absorption_detected, civ_center_absorption_detected = get_absorption_intervals(lam=lam, flux=flux_clean, flux_sub=flux_sub_civ, continuum=continuum_civ, noise=noise, window=CIV_FIT_WINDOW, verbose=verbose)

    # fit single gaussian, use result as parameter guesses for two gaussian fit
    result1_civ, profile1_civ = fit_single_gaussian(lam=lam, flux_sub=flux_sub_civ, ivar=ivar, mask=civ_absorption_mask, window=CIV_FIT_WINDOW)
    result2_civ, profile2_civ = fit_two_gaussians(lam=lam, flux_sub=flux_sub_civ, ivar=ivar, mask=civ_absorption_mask, single_result=result1_civ, window=CIV_FIT_WINDOW)

    # F-test which fit is better
    accept_two_civ, p_val_civ = ftest_which_gaussian(result1=result1_civ, result2=result2_civ)

    # depending on F-test result, assign respective profile, result and best parameter dictionary to general variables
    profile_civ = profile2_civ if accept_two_civ else profile1_civ
    result_civ = result2_civ if accept_two_civ else result1_civ
    params_dict_civ = result_civ.params.valuesdict()

    #------------------------------------------------------------------------------------

    ## CONTINUUM BELOW NV AND LYA

    result_lya_nv_cont, continuum_nv_lya, _ = fit_continuum(lam=lam, flux=flux_clean, ivar=ivar, windows=NV_LYA_CONTINUUM_WINDOWS+CIV_CONTINUUM_WINDOWS, lambda_ref=1290.0)  # fit continuum
    flux_sub_nv_lya = flux_clean - continuum_nv_lya  # subtract continuum
    

    ## NV PROFILE

    # make LYA mask
    lya_mask, lya_mask_lims = make_scaling_mask(lam=lam, profile_civ=profile_civ, lims=True)

    # fit NV by scaling CIV template
    result_nv, profile_nv = scale_civ_to_nv(lam=lam, flux_sub=flux_sub_nv_lya, params_dict=params_dict_civ, accept_two=accept_two_civ, ivar=ivar, mask=lya_mask)


    ## LYA PROFILE

    # subtract NV profile
    flux_sub_lya = flux_sub_nv_lya - profile_nv

    # get absorption mask
    lya_absorption_intervals, lya_absorption_mask, lya_blue_absorption_detected, lya_center_absorption_detected = get_absorption_intervals(lam=lam, flux=flux_clean, flux_sub=flux_sub_lya, continuum=continuum_nv_lya, noise=noise, window=LYA_NV_FIT_WINDOW, verbose=verbose)

    # fit single gaussian, use result as parameter guesses for two gaussian fit
    result1_lya, profile1_lya = fit_single_gaussian(lam=lam, flux_sub=flux_sub_lya, ivar=ivar, mask=lya_absorption_mask, window=LYA_NV_FIT_WINDOW)
    result2_lya, profile2_lya = fit_two_gaussians(lam=lam, flux_sub=flux_sub_lya, ivar=ivar, mask=lya_absorption_mask, single_result=result1_lya, window=LYA_NV_FIT_WINDOW)

    # F-test which fit is better, assign to general variable
    accept_two_lya, p_val_lya = ftest_which_gaussian(result1=result1_lya, result2=result2_lya)

    # depending on F-test result, assign respective profile and result to general variables
    profile_lya = profile2_lya if accept_two_lya else profile1_lya
    result_lya = result2_lya if accept_two_lya else result1_lya

    #------------------------------------------------------------------------------------

    ## LINE MEASUREMENT

    # CIV
    line_stats_civ, civ_integration_win = measure_line(lam=lam, profile=profile_civ, flux_sub=flux_sub_civ, continuum=continuum_civ, window=CIV_FIT_WINDOW,get_window=True, verbose=verbose) 

    # NV
    flux_sub_nv = flux_sub_nv_lya - profile_lya  # subtract LYA profile
    line_stats_nv, nv_integration_win = measure_line(lam=lam, profile=profile_nv, flux_sub=flux_sub_nv, continuum=continuum_nv_lya, window=LYA_NV_FIT_WINDOW, get_window=True, verbose=verbose)

    # LYA
    line_stats_lya, lya_integration_win = measure_line(lam=lam, profile=profile_lya, flux_sub=flux_sub_lya, continuum=continuum_nv_lya, window=LYA_NV_FIT_WINDOW, get_window=True, verbose=verbose)

    #------------------------------------------------------------------------------------

    ## ADDITIONAL METRICS

    # signal-to-noise ratio at each line + 1700A
    snr_1700 = snr_around_lam(lam=lam, flux=flux_sub_civ, noise=noise, ref_lam=1700)
    snr_civ = snr_around_lam(lam=lam, flux=flux_sub_civ, noise=noise, ref_lam=CIV_AIR)
    snr_nv = snr_around_lam(lam=lam, flux=flux_sub_nv_lya, noise= noise, ref_lam=NV_AIR)
    snr_lya = snr_around_lam(lam=lam, flux=flux_sub_nv_lya, noise=noise, ref_lam=LYA_AIR)

    # ratio of pixels used in fits for validation
    civ_pixel_used = pixel_used_in_window(lam=lam, mask=civ_absorption_mask, lower_lim=CIV_FIT_WINDOW[0], upper_lim=CIV_FIT_WINDOW[1])
    lya_blue_end_pixel_used = pixel_used_in_window(lam=lam, mask=lya_absorption_mask, lower_lim=LYA_NV_FIT_WINDOW[0], upper_lim=NV_AIR)  # only consider the left part of the fit window, up to NV

    #------------------------------------------------------------------------------------

    results = {
        "targetid":                 idx,

        # REW
        "rew_civ":                  line_stats_civ["ew_aa"],
        "rew_nv":                   line_stats_nv["ew_aa"],
        "rew_lya":                  line_stats_lya["ew_aa"],

        # SNR
        "snr_1700A":                snr_1700,
        "snr_civ":                  snr_civ,
        "snr_nv":                   snr_nv,
        "snr_lya":                  snr_lya,

        # ratio of pixels used in fits
        "pixel_used_civ":           civ_pixel_used,
        "pixel_used_blue_end_lya":  lya_blue_end_pixel_used,

        # continuum reduced chi square
        "continuum_redchi":         result_lya_nv_cont.redchi,

        
        # other line measurements
        "measurements_civ": {
            "integration_window":   civ_integration_win,
            "further_measurements": line_stats_civ
            },
        "measurements_nv": {
            "integration_window":   nv_integration_win,
            "further_measurements": line_stats_nv
            },
        "measurements_lya": {
            "integration_window":   lya_integration_win,
            "further_measurements": line_stats_lya
            },
            
        # CIV continuum fit details
        "continuum_fit_civ": {
            "redchi":               result_civ_cont.redchi,
            "fit_params":           result_civ_cont.params.valuesdict()
            },

        # CIV absorption details
        "absorption_civ": {
            "intervals":            civ_absorption_intervals,
            "blue_detected":        civ_blue_absorption_detected,
            "center_detected":      civ_center_absorption_detected
            },
        
        # CIV line fit details
        "line_fit_civ": {
            "accept_two":           accept_two_civ,
            "p_val":                p_val_civ,
            "redchi":               result_civ.redchi,
            "fit_params":           params_dict_civ
            },

        # LYA & NV continuum fit details
        "continuum_fit_lya_nv": {
            "redchi":               result_lya_nv_cont.redchi,
            "fit_params":           result_lya_nv_cont.params.valuesdict()
            },

        # NV line fit details
        "line_fit_nv": {
            "lya_mask":             lya_mask_lims,
            "redchi":               result_nv.redchi,
            "fit_params":           result_nv.params.valuesdict()
            },

        # LYA absorption details
        "absorption_lya": {
            "intervals":            lya_absorption_intervals,
            "blue_detected":        lya_blue_absorption_detected,
            "center_detected":      lya_center_absorption_detected
            },  

        # LYA line fit details
        "line_fit_lya": {
            "accept_two":           accept_two_lya,
            "p_val":                p_val_lya,
            "redchi":               result_lya.redchi,
            "fit_params":           result_lya.params.valuesdict()
            },    

        }

    return results