File containing the list of DESI objects with targetid and relevant properties: target_properties.csv

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