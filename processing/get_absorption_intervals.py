import sys
sys.path.append('../utility')

from utility.get_intercepts_near_line import get_intercepts_near_line
from .parameters import BAL_MIN_WIDTH

import numpy as np
import pandas as pd
from astropy.convolution import convolve, Gaussian1DKernel, Box1DKernel
from scipy.signal import argrelmin, argrelmax


def get_absorption_intervals(lam, flux, ivar, continuum, window=(1450,1650), min_width_aa=BAL_MIN_WIDTH):
    mask = (lam > window[0]) & (lam < window[1])

    lam, flux, ivar = lam[mask], flux[mask], ivar[mask]

    noise = np.where(ivar > 0, 1.0 / np.sqrt(np.maximum(ivar, 1e-30)), np.inf)
    continuum = continuum[mask]

    center_absorption, blue_absorption = None, None
    absorption_intervals = []

    ## BLUE ABSORPTION
    # use rolling mean because we want to smooth but stay as much to the original values as possible
    flux_conv_mean = convolve(flux,Box1DKernel(5))

    # get indeces of "global" minimum and maximum
    max_idx = np.argmax(flux_conv_mean)
    min_idx = np.argmin(flux_conv_mean[:max_idx])

    if flux_conv_mean[min_idx] < continuum[min_idx] - noise[min_idx]:
        try:
            blue_absorption = get_intercepts_near_line(lam,lower_flux=continuum,upper_flux=flux_conv_mean,x_val=lam[min_idx])
            if blue_absorption[1] - blue_absorption[0] >= min_width_aa:
                absorption_intervals.append(blue_absorption)
        except IndexError:
            print("IndexError while trying to find blue absorption boundaries")
        

    ## CENTER ABSORPTION
    # use gaussian smoothing since it strongly removes noise
    flux_conv_gauss = convolve(flux, Gaussian1DKernel(3))

    # extract indeces of local minima and maxima
    local_mins_idx = argrelmin(flux_conv_gauss)
    local_maxs_idx = argrelmax(flux_conv_gauss)
    # make mins and max lists with lam, flux and noise for each extremum
    mins_list = [(lam[i],flux_conv_gauss[i],noise[i],0) for i in local_mins_idx[0]]
    maxs_list = [(lam[i],flux_conv_gauss[i],noise[i],1) for i in local_maxs_idx[0]]
    # put into dataframe
    extrema = pd.DataFrame(sorted(mins_list + maxs_list),columns=["lam","flux_conv","noise","max_bool"])

    # TODO: do all this with numpy and not pandas

    # calculate forward and backward flux differences
    extrema["flux_diff_backward"] = extrema["flux_conv"].diff().abs()
    extrema["flux_diff_forward"] = extrema["flux_conv"].diff(-1).abs()

    # for each extremum, check if the absolute flux difference to both neighbors is more than 1.5*noise
    extrema["more_than_noise"] = (extrema["flux_diff_backward"] > 1.5*extrema["noise"]) & (extrema["flux_diff_forward"] > 1.5*extrema["noise"])
    # boolean column where the maximum value is true
    extrema["overall_max_bool"] = (extrema["flux_conv"] == extrema["flux_conv"].max())

    # select minima (max_bool==0, i.e. False) that are more than noise
    query = extrema.query("more_than_noise == 1 and max_bool == 0")

    if len(query) > 0:
        indeces = query.index  # get indeces of selected minima
        for idx in indeces:
            # for each minima, check if one of its neighbors the overall maximum
            if extrema["overall_max_bool"][idx-1].squeeze() == True or extrema["overall_max_bool"][idx+1].squeeze() == True:
                
                # get the lambdas of the neighbors
                lower = extrema["lam"][idx-1].squeeze()
                upper = extrema["lam"][idx+1].squeeze()
                five_percent = (upper - lower) * 0.05  # add puffer zone

                center_absorption = (lower+five_percent, upper-five_percent)
                absorption_intervals.append(center_absorption)


    # if no absorption, return nothing
    if len(absorption_intervals) == 0:
        return None, None

    # if absorption, make mask from intervals
    elif len(absorption_intervals) == 1:
        absorption_mask = (lam > absorption_intervals[0][0]) & (lam < absorption_intervals[0][1])
    else:
        absorption_mask = ((lam > absorption_intervals[0][0]) & (lam < absorption_intervals[0][1])) | ((lam > absorption_intervals[1][0]) & (lam < absorption_intervals[1][1]))
    
    return absorption_intervals, absorption_mask
