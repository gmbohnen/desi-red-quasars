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


    ## FIT CONTINUUM BELOW CIV
    result_civ_cont, continuum_civ, _ = fit_continuum(lam, flux_clean, ivar=ivar, windows=CIV_CONTINUUM_WINDOWS, lambda_ref=1700.0)


    ## QUALITY CUTS similar to Hamann (based on CIV region)

    # TODO: fix all this, like i think it doesnt do shit currently and you need to get the values fixed
    # rejected, reason, flags = quality_cuts(lam, flux_clean, continuum_civ, result_civ_cont, max_slope=10.0)

    # if rejected:
    #     pass
        

    ## FIT CIV LINE

    # subtract CIV continuum
    flux_sub_civ = flux_clean - continuum_civ

    # get absorption mask
    civ_absorption_intervals, civ_absorption_mask, civ_blue_absorption_found, civ_center_absorption_found  = get_absorption_intervals(lam,flux_clean,flux_sub_civ,continuum_civ,noise, window=CIV_FIT_WINDOW)

    ## ADDITIONAL METRICS

    # signal-to-noise ratio at each line + 1700A
    snr_1700 = snr_around_lam(lam,flux_sub_civ,noise)
    snr_civ = snr_around_lam(lam,flux_sub_civ,noise,ref_lam=CIV_AIR)
    

    results = {
        "targetid":         idx,

        # REW
        "rew_civ":          ,
        
        # SNR
        "snr_1700A":        snr_1700,
        "snr_civ":          snr_civ,

        # CIV absorption
        "absorption_civ":   {
            "intervals":        civ_absorption_intervals,
            "blue_found":       civ_blue_absorption_found,
            "center_found":     civ_center_absorption_found
            }
        }

    return rew_dic