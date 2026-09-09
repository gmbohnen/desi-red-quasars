from .plot_spectrum import plot_spectrum
from .my_colors import pink, palette_dark
from pipeline.processing.fit_gaussians import _gaussian, _two_gaussian
from pipeline.processing.scale_civ_to_nv import _civ_template
from pipeline.processing.parameters import CIV_FIT_WINDOW, LYA_NV_FIT_WINDOW

import pandas as pd
from astropy.convolution import convolve, Box1DKernel

import h5py
import ast

SPECTRA_PATH = "/home/leya/Code/Uni/desi-red-quasars/data/spectra.h5"
DETAILS_PATH = "/home/leya/Code/Uni/desi-red-quasars/data/result_details.csv"


def plot_from_hdf5(targetid, df=None, CIV=False, NV=False, LYA=False, smoothing_strength=3, **kwargs):
    targetid = str(targetid)

    with h5py.File(SPECTRA_PATH, "r") as f:
        data = f[f"{targetid}"][()]
    
    lam = data[0,:]
    flux = data[1,:]


    # load dataframe if it not given
    if (CIV or NV or LYA) and (df is None):
        df = pd.read_csv(DETAILS_PATH,index_col="targetid",dtype={"targetid":"str"})

    #---------------------------------------------------------------------------------------

    if CIV:
        # load civ params
        civ_params = ast.literal_eval(df.loc[targetid,"line_fit_civ"])
        civ_integration_window = ast.literal_eval(df.loc[targetid,"measurements_civ"])["integration_window"]

        # create profile
        if civ_params["accept_two"]:
            civ_profile = _two_gaussian(lam,civ_params["fit_params"]["amp_c"],civ_params["fit_params"]["cen_c"],civ_params["fit_params"]["sigma_c"],civ_params["fit_params"]["amp_w"],civ_params["fit_params"]["cen_w"],civ_params["fit_params"]["sigma_w"])
        else:
            civ_profile = _gaussian(lam, civ_params["fit_params"]["amplitude"],civ_params["fit_params"]["center"],civ_params["fit_params"]["sigma"])

        # create continuum profile
        civ_cont_params = ast.literal_eval(df.loc[targetid,"continuum_fit_civ"])
        civ_cont_profile = civ_cont_params["fit_params"]["amplitude"] * ((lam / 1700.0) ** civ_cont_params["fit_params"]["alpha"])

        # add continuum to line profile
        civ_profile += civ_cont_profile

        # create arrays to plot later
        civ_mask = (lam > civ_integration_window[0]) & (lam < civ_integration_window[1])
        civ_lam = lam[civ_mask]
        civ_profile = civ_profile[civ_mask]

    else:
        civ_profile = None
        civ_lam = None

    #---------------------------------------------------------------------------------------

    # create continuum below NV and LYA
    if NV or LYA:
        cont_params = ast.literal_eval(df.loc[targetid,"continuum_fit_lya_nv"])
        
        cont_profile = cont_params["fit_params"]["amplitude"] * ((lam / 1290.0) ** cont_params["fit_params"]["alpha"])

    # else:
    #     lya_lam = None
    #     nv_lam = None

    #---------------------------------------------------------------------------------------

    if NV:
        nv_params = ast.literal_eval(df.loc[targetid,"line_fit_nv"])
        nv_integration_window = ast.literal_eval(df.loc[targetid,"measurements_nv"])["integration_window"]

        if "civ_params" not in locals():  # get CIV values if it does not yet exist
            civ_params = ast.literal_eval(df.loc[targetid,"line_fit_civ"])

        # create profile
        nv_profile = _civ_template(lam, nv_params["fit_params"]["scale"], nv_params["fit_params"]["shift"], civ_params["fit_params"], civ_params["accept_two"])

        nv_profile += cont_profile  # add continuum to line profile
        
        # create arrays to plot later
        nv_mask = (lam > nv_integration_window[0]) & (lam < nv_integration_window[1])
        nv_lam = lam[nv_mask]
        nv_profile = nv_profile[nv_mask]

    else:
        nv_profile = None
        nv_lam = None

    #---------------------------------------------------------------------------------------

    if LYA:
        lya_params = ast.literal_eval(df.loc[targetid,"line_fit_lya"])
        lya_integration_window = ast.literal_eval(df.loc[targetid,"measurements_lya"])["integration_window"]

        # create profile
        if lya_params["accept_two"]:
            lya_profile = _two_gaussian(lam,lya_params["fit_params"]["amp_c"],lya_params["fit_params"]["cen_c"],lya_params["fit_params"]["sigma_c"],lya_params["fit_params"]["amp_w"],lya_params["fit_params"]["cen_w"],lya_params["fit_params"]["sigma_w"])
        else:
            lya_profile = _gaussian(lam, lya_params["fit_params"]["amplitude"],lya_params["fit_params"]["center"],lya_params["fit_params"]["sigma"])

        lya_profile += cont_profile  # add continuum to line profile

        # create arrays to plot later
        lya_mask = (lam > lya_integration_window[0]) & (lam < lya_integration_window[1])
        lya_lam = lam[lya_mask]
        lya_profile = lya_profile[lya_mask]


    else:
        lya_profile = None
        lya_lam = None

    #---------------------------------------------------------------------------------------

    plot_spectrum(lam,convolve(flux,Box1DKernel(smoothing_strength)),
        color="black",
        additional_lam=[lya_lam, nv_lam, civ_lam],
        additional_flux=[lya_profile, nv_profile, civ_profile],
        additional_flux_color=[palette_dark[0],palette_dark[2],palette_dark[1]],
        line_width_factor=2,
        **kwargs
        )