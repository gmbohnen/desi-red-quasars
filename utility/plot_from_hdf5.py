from matplotlib.pyplot import plot
from .plot_spectrum import plot_spectrum

import h5py

spectra_hdf5_path = #TODO

def plot_from_hdf5(targetid):
    with h5py.File(spectra_hdf5_path, "r") as f:
        data = f[f"{idx}"][()]
        
    lam = data[0,:]
    flux = data[1,:]

    plot_spectrum(lam,flux,add_smoothed=True)