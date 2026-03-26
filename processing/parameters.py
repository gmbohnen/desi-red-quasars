# TODO check all these
continuum_windows = [
    (1425, 1470),   # blue side of CIV, relatively clean
    (1680, 1710),   # red side of CIV, between HeII and CIII]
    (1760, 1810),   # far red continuum anchor
]
 
# CIV emission line region (used for significance check and BAL detection)
civ_region = (1500, 1600)   # Å, rest-frame
 
# Quality cut thresholds
max_powerlaw_slope = 3.0    # reject if |alpha| > this (unrealistically steep)
bal_sigma_threshold = -2.0  # reject if median flux in CIV window is this many sigma BELOW the continuum (broad absorption)
min_emission_sigma = 2.0    # reject if peak emission above continuum is less than this many sigma (no significant emission)

