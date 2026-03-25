def adjust_redshift(z, wavelength, flux, ivar=None):
    wavelength = wavelength/(1+z)
    flux = flux*(1+z)
    
    if ivar is not None:
        ivar = ivar/((1+z)**2)
        return wavelength, flux, ivar

    return wavelength, flux