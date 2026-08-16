from processing import measure_line
from utility.adjust_redshift import adjust_redshift

from processing.remove_spikes import remove_spikes
from processing.fit_continuum import fit_continuum
from processing.quality_cuts import quality_cuts
from processing.get_absorption_intervals import get_absorption_intervals
from processing.fit_gaussians import fit_single_gaussian, fit_two_gaussians, ftest_which_gaussian
from processing.measure_line import measure_line
from processing.scale_civ_to_nv import scale_civ_to_nv
from processing.parameters import CIV_CONTINUUM_WINDOWS, CIV_FIT_WINDOW, LYA_NV_FIT_WINDOW, NV_LYA_CONTINUUM_WINDOWS, CIV_AIR, NV_AIR, LYA_AIR

import numpy as np
import h5py


spectra_hdf5_path = "data/batches/spectra.h5"  # TODO move to a config file?


def pipeline(idx):

    with h5py.File(spectra_hdf5_path, "r") as f:
        data = f[f"{idx}"][()]
        
    lam = data[0,:]
    flux = data[1,:]
    ivar = data[2,:]

    # preprocessing (redshift removal and bad pixel masking) already done in retrieval pipeline

    # create noise
    noise = np.abs(1/np.sqrt(ivar))

    # remove spikes
    flux_clean, _ = remove_spikes(lam,flux,ivar=ivar)


    ## CHECK SNR @ 1700A 
    # TODO


    ## FIT CONTINUUM BELOW CIV
    result_civ_cont, continuum_civ, _ = fit_continuum(lam, flux_clean, ivar=ivar, windows=CIV_CONTINUUM_WINDOWS, lambda_ref=1700.0)


    ## QUALITY CUTS (based on CIV region)

    # TODO: fix all this, like i think it doesnt do shit currently and you need to get the values fixed
    rejected, reason, flags = quality_cuts(lam, flux_clean, continuum_civ, result_civ_cont, max_slope=10.0)

    if rejected:
        print(f'Rejected because: {reason}')
        
    else:
        ## FIT CIV LINE

        # subtract CIV continuum
        flux_sub_civ = flux_clean - continuum_civ if not rejected else None

        # get absorption mask
        _, civ_absorption_mask = get_absorption_intervals(lam,flux_sub_civ,noise,continuum_civ, window=CIV_FIT_WINDOW)
        # civ_absorption_mask = (flux_sub_civ < -noise)

        # fit single gaussian, use result as parameter guesses for two gaussian fit
        result1_civ, profile1_civ = fit_single_gaussian(lam, flux_sub_civ, ivar=ivar, mask=civ_absorption_mask, window=CIV_FIT_WINDOW)
        result2_civ, profile2_civ = fit_two_gaussians(lam, flux_sub_civ, ivar=ivar, mask=civ_absorption_mask, single_result=result1_civ, window=CIV_FIT_WINDOW)

        # F-test which fit is better
        accept_two_civ, p_val_civ = ftest_which_gaussian(result1_civ, result2_civ)

        # depending on F-test result, assign respective profile and best parameter dictionary to general variables
        profile_civ = profile2_civ if accept_two_civ else profile1_civ
        result_civ = result2_civ if accept_two_civ else result1_civ
        params_dict_civ = result_civ.params.valuesdict()


        ## CONTINUUM BELOW NV AND LYA
        _, continuum_nv_lya, _ = fit_continuum(lam, flux_clean, ivar=ivar, windows=NV_LYA_CONTINUUM_WINDOWS+CIV_CONTINUUM_WINDOWS, lambda_ref=1290.0)  # fit continuum
        flux_sub_nv_lya = flux_clean - continuum_nv_lya  # subtract continuum


        ## NV PROFILE

        # make LYA mask
        # idea: mask left and right by a certain fixed window (left: 60% of lam(NV)-lam(LYA), right 40%),
        # then shift by the difference between CIV fit peak and air wavelengths to make it more robust in case of bad redshift values
        # TODO move to different file
        diff_civ = lam[np.argmax(profile_civ)] - CIV_AIR
        m = (NV_AIR-LYA_AIR)
        c = LYA_AIR + diff_civ
        lims = (c-0.6*m,c+0.6*m)
        lya_mask = (lam > lims[0]) & (lam < lims[1])
        
        # fit NV by scaling CIV template
        result_nv, profile_nv = scale_civ_to_nv(lam,flux_sub_nv_lya,params_dict_civ,accept_two_civ,ivar=ivar,mask=lya_mask)

        ## LYA PROFILE

        # subtract NV profile
        flux_sub_lya = flux_sub_nv_lya - profile_nv

        # get absorption mask
        _, lya_absorption_mask = get_absorption_intervals(lam,flux_sub_lya,continuum_nv_lya,noise,window=LYA_NV_FIT_WINDOW)  # first is intervals for plotting

        # fit single gaussian, use result as parameter guesses for two gaussian fit
        result1_lya, profile1_lya = fit_single_gaussian(lam, flux_sub_lya, ivar=ivar, mask=lya_absorption_mask, window=LYA_NV_FIT_WINDOW)
        result2_lya, profile2_lya = fit_two_gaussians(lam, flux_sub_lya, ivar=ivar, mask=lya_absorption_mask, single_result=result1_lya, window=LYA_NV_FIT_WINDOW)

        # F-test which fit is better, assign to general variable
        accept_two_lya, p_val_lya = ftest_which_gaussian(result1_lya, result2_lya)
        profile_lya = profile2_lya if accept_two_lya else profile1_lya
        result_lya = result2_lya if accept_two_lya else result1_lya

        

        ## LINE MEASUREMENT

        # CIV
        line_stats_civ, civ_integration_win = measure_line(lam, profile_civ, flux_sub_civ, continuum_civ, window=CIV_FIT_WINDOW,get_window=True) 

        # NV
        flux_sub_nv = flux_sub_nv_lya - profile_lya  # subtract Lya profile
        line_stats_nv, nv_integration_win = measure_line(lam,profile_nv,flux_sub_nv,continuum_nv_lya,window=LYA_NV_FIT_WINDOW,get_window=True)

        # Lya
        line_stats_lya, lya_integration_win = measure_line(lam,profile_lya,flux_sub_lya,continuum_nv_lya,window=LYA_NV_FIT_WINDOW,get_window=True)
        

        ## RESULTS
        
        rew_dic = {
            "rew_civ" : line_stats_civ["ew_aa"],
            "rew_nv"  : line_stats_nv["ew_aa"],
            "rew_lya" : line_stats_lya["ew_aa"]
        }

        return rew_dic