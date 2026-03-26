# TODO check all these
CONTINUUM_WINDOWS = [
    (1425, 1470),   # blue side of CIV, relatively clean
    (1680, 1710),   # red side of CIV, between HeII and CIII]
    (1760, 1810),   # far red continuum anchor
]
 
# CIV emission line region (used for significance check and BAL detection)
CIV_REGION = (1500, 1600)   # Å, rest-frame
 
# Quality cut thresholds
MAX_POWERLAW_SLOPE = 3.0    # reject if |alpha| > this (unrealistically steep)
BAL_SIGMA_THRESHOLD = -2.0  # reject if median flux in CIV window is this many sigma BELOW the continuum (broad absorption)
MIN_EMISSION_SIGMA = 2.0    # reject if peak emission above continuum is less than this many sigma (no significant emission)

