target_properties.csv:

        File containing the list of DESI objects with 'targetid' and relevant properties

        Query used to retrieve it from the DataLab QueryClient:

        query = """
                SELECT zp.targetid, zp.z, zp.zerr, agn.civ_1549_flux, agn.civ_1549_flux_ivar, agn.civ_1549_sigma, agn.flux_w3, agn.flux_ivar_w3, zp.desiname
                FROM desi_dr1.zpix AS zp
                JOIN desi_dr1.agngal as agn
                ON zp.targetid = agn.targetid
                WHERE (zp.z BETWEEN 2.2 AND 4.38) AND (zp.zwarn=0) AND (agn.civ_1549_flux>0) AND (zp.spectype='QSO') AND zp.zcat_primary
                """

        - zwarn=0 means no issues with the pipeline
        - zcat_primary means it is the primary spectrum for the target


duplicate_targets.txt:

        File containing the 'targetid's of all objects that had more than one spectrum in SPARCL


spectra.h5:

        File containing the spectra, i.e. wavelength, flux and ivar arrays, for each object listed in 'target_properties.csv' that had only one spectrum associated with its 'targetid' in SPARCL.

        Each spectrum is a separate dataset, the arrays can be accessed as follows:

                '''
                with h5py.File(spectra_hdf5_path, "r") as f:
                        data = f[targetid][()]
                
                lam = data[0,:]
                flux = data[1,:]
                ivar = data[2,:]
                '''