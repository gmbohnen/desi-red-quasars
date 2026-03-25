import numpy as np


def get_intercepts_near_line(lam, upper_flux, lower_flux, line_lambda=1549):
    '''Find two intercepting points of two fluxes that are closest to a certain x value.'''
    diffs = upper_flux - lower_flux  # calculate difference
    sign_changes = np.where(np.diff(np.sign(diffs)))  # find where sign changes to get all intercepts

    vals_at_sign_change = lam[sign_changes]  # get lambda values at intercepts

    # separate into arrays for values below and over the line
    before_line = [x for x in vals_at_sign_change if x < line_lambda]
    behind_line = [x for x in vals_at_sign_change if x > line_lambda]

    # sort lambda values by their distance to the line, then select closest each
    lower_intercept = sorted(before_line, key=lambda x: abs(x - line_lambda))[0]  
    upper_intercept = sorted(behind_line, key=lambda x: abs(x - line_lambda))[0]

    return (lower_intercept, upper_intercept)