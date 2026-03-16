from sparcl.client import SparclClient
from dl import queryClient as qc
from desitarget import targetmask

from polars import DataFrame as pl_df
from pandas import DataFrame as pd_df
from astropy.table import Table


def retrieve_desi_spectra(query, output_fields, target_mask=None, format="polars"):
    '''
    Retrieve spectra from DESI-DR1.
    
    Parameters
    ----------
    query : str
        SQL query, must include 'targetid' or 'specid' in the SELECT statement.
    output_fields : list of str
        List of fields to be included in the output, must include 'flux' and 'wavelength'.
    target_mask : list of str or None, optional
        List of targetmask categories to be filtered for.
        See `desitarget.targetmask.desi_mask` for valid categories.
    format : {'polars', 'pandas', 'dict'}, optional
        Format in which the output is provided.

    Returns
    -------
    table
        Table in the specified format containing the columns specified in output_fields.
    '''

    assert ("targetid" in query) or ("specid" in query) or ("SELECT *" in query) or ("select *" in query), "'targetid' or 'specid' must be in the SELECT statement of the query to retrieve spectra."

    assert ("flux" in output_fields) and ("wavelength" in output_fields), "'flux' and 'wavelength' must be in output_fields to retrieve spectra."

    assert (format in ['polars', 'pandas', 'dict']), "'format must be one of {'polars', 'pandas', 'dict'}."

    if target_mask:
        assert ("desi_target" in query) or ("SELECT *" in query) or ("select *" in query), "To filter for targetmask categories, 'desi_target' must be in the SELECT statement."

        if len(target_mask) > 1:
            raise NotImplementedError("Filtering for more than one targetmask category is currently not implemented.")


    # retrieve data matching the query
    sql_results = qc.query(sql=query, fmt="table")
    
    # filter for specific targetmask categories if wanted
    if target_mask:
        desi_mask = targetmask.desi_mask  # retrieve the mask
        desi_target = sql_results["desi_target"]  # get the column that contains target info

        target_filter = (desi_target & desi_mask[target_mask[0]] != 0)  # create filter using binary &

        sql_results = sql_results[target_filter]
    

    ## RETRIEVE SPECTRA
    client = SparclClient()

    idxs = [int(x) for x in sql_results["targetid"]]  # extract IDs of desired objects, must be python list of python integers

    spectra_results = client.retrieve_by_specid(specid_list=idxs, include=output_fields, dataset_list=["DESI-DR1"])


    # output in desired format
    if format == "polars":
        out = pl_df(spectra_results.data[1:])
    elif format == "pandas":
        out = pd_df(spectra_results.data[1:])
    else:
        out = Table(spectra_results)
    

    return out