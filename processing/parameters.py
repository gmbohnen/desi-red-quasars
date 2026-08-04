# TODO check all these

CIV_CONTINUUM_WINDOWS = [ 
    (1425, 1470),
    (1680, 1710),
    (1760, 1810)
]

NV_LYA_CONTINUUM_WINDOWS = [   # Kramer2009
    (1155, 1165),
    (1280,1293),
    (1315,1325),
    (1340,1360)
]


# CIV emission line region (used for significance check and BAL detection) 
# i think this isnt used at this point?
CIV_REGION = (1500, 1600)      # Angstrom, rest-frame


# Quality cut thresholds
MAX_POWERLAW_SLOPE  = 3.0      # reject if |alpha| > this (unrealistically steep)
BAL_SIGMA_THRESHOLD = -2.0     # reject if median flux in CIV window is this many sigma BELOW the continuum (broad absorption)
MIN_EMISSION_SIGMA  = 2.0      # reject if peak emission above continuum is less than this many sigma (no significant emission)

# Fit quality: reduced chi square threshold
MAX_REDCHI = 3                 # Meija-Restrepo2018

# line rest wavelength (Angstrom) using vacuum wavelengths from Harris2016, then Edlén formula applied to convert to air
LYA_AIR        = 1215.34     
NV_AIR         = 1239.81
CIV_AIR        = 1549.06       # doublet ignored; (CIV-1 + CIV-2) / 2


C_KMS          = 2.99792458e5  # speed of light (km/s) 

# Fit regions
CIV_FIT_WINDOW = (1500,1600)   # Shen2011
LYA_NV_FIT_WINDOW = (1150,1290)# from Shen2019


# BAL masking thresholds
BAL_SIGMA_MASK  = 2.5
BAL_MIN_WIDTH   = 5            # Angstrom


# Two-Gaussian constraints
MIN_WIDTH_RATIO = 1.2          # wing sigma >= this * core sigma
MAX_WIDTH_RATIO = 3.0#4.0      # wing sigma <= this * core sigma
MIN_CORE_SIGMA  = 3.0          # Angstrom
MAX_CORE_SIGMA  = 40.0#80.0    # Angstrom (~15000 km/s FWHM)
MAX_WING_SIGMA  = 80.0#200.0   # Angstrom


# F-test significance threshold
FTEST_PVALUE    = 0.05
