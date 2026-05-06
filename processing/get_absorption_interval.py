import numpy as np
import pandas as pd
from astropy.convolution import convolve, Gaussian1DKernel
from scipy.signal import argrelmin, argrelmax

def get_absorption_interval(lam, flux, ivar):
    mask = (lam > 1450) & (lam < 1650)

    lam, flux, ivar = lam[mask], flux[mask], ivar[mask]
    flux_conv = convolve(flux, Gaussian1DKernel(3))

    noise = np.where(ivar > 0, 1.0 / np.sqrt(np.maximum(ivar, 1e-30)), np.inf)

    mins = argrelmin(flux_conv)
    maxs = argrelmax(flux_conv)
    extrema_idxs = sorted(list(mins[0])+list((maxs[0])))
    mins_list = [(lam[i],flux_conv[i],noise[i],0) for i in mins[0]]
    maxs_list = [(lam[i],flux_conv[i],noise[i],1) for i in maxs[0]]
    extrema = pd.DataFrame(sorted(mins_list + maxs_list),columns=["lam","flux_conv","noise","max_bool"])
    extrema["flux_diff_backward"] = extrema["flux_conv"].diff().abs()
    extrema["flux_diff_forward"] = extrema["flux_conv"].diff(-1).abs()
    extrema["more_than_noise"] = (extrema["flux_diff_backward"] > 1.5*extrema["noise"]) & (extrema["flux_diff_forward"] > 1.5*extrema["noise"])
    extrema["overall_max_bool"] = extrema["flux_conv"] == extrema["flux_conv"].max()

    query = extrema.query("more_than_noise == 1 and max_bool == 0")

    if len(query) > 0:
        idx = query.index
        try:
            if extrema["overall_max_bool"][idx-1].squeeze() == True or extrema["overall_max_bool"][idx+1].squeeze() == True:

                lower = extrema["lam"][idx-1].squeeze()
                upper = extrema["lam"][idx+1].squeeze()
                five_percent = (upper - lower) * 0.05  # um nen puffer zu haben quasi

                absorption_interval = (lower+five_percent, upper-five_percent)

                return absorption_interval

        except ValueError:
            print("ValueError in get_absorption_interval()")

    return None