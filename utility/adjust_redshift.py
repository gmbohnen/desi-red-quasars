def adjust_redshift(wavelength, flux, z):
    wavelength = wavelength/(1+z)
    flux = flux*(1+z)
    
    return wavelength, flux