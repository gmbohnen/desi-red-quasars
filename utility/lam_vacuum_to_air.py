def lam_vacuum_to_air(lam_vac):
    '''
        Convert vacuum wavelengths to air wavelenghts, using the Edlén formula.
        
        Context: Wavelenghts at which emission/absorption is observed differs in vacuum and air. Wavelenghts below 2000 Angstrom are usually given in the vacuum.
    '''

    n = 1.00027
    lam_air = lam_vac / n

    return lam_air


if __name__ == "__main__":
    import sys

    lam_vac = float(sys.argv[1])
    lam_air = lam_vacuum_to_air(lam_vac)
    print(round(lam_air,4))
