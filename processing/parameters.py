# TODO check all these

CONTINUUM_WINDOWS = [
    (1425, 1470),             # blue side of CIV, relatively clean
    (1680, 1710),             # red side of CIV, between HeII and CIII]
    (1760, 1810),             # far red continuum anchor
]

# CIV emission line region (used for significance check and BAL detection)
CIV_REGION = (1500, 1600)     # Angstrom, rest-frame

# Quality cut thresholds
MAX_POWERLAW_SLOPE = 3.0      # reject if |alpha| > this (unrealistically steep)
BAL_SIGMA_THRESHOLD = -2.0    # reject if median flux in CIV window is this many sigma BELOW the continuum (broad absorption)
MIN_EMISSION_SIGMA = 2.0      # reject if peak emission above continuum is less than this many sigma (no significant emission)

C_KMS        = 2.99792458e5   # speed of light (km/s)
CIV_VACUUM   = 1549.48        # CIV rest wavelength (Angstrom), doublet ignored
FIT_WINDOW   = (1450, 1650)   # fitting window (rest-frame Angstrom)
 
# Two-Gaussian constraints
MIN_WIDTH_RATIO = 1.2         # wing sigma >= this * core sigma
MAX_WIDTH_RATIO = 4.0         # wing sigma <= this * core sigma
MIN_CORE_SIGMA  = 3.0         # Angstrom
MAX_CORE_SIGMA  = 80.0        # Angstrom (~15000 km/s FWHM)
MAX_WING_SIGMA  = 200.0       # Angstrom
 
# F-test significance threshold
FTEST_PVALUE    = 0.05
 
# BAL masking thresholds
BAL_SIGMA_MASK  = 2.5
BAL_MIN_WIDTH   = 10          # Angstrom