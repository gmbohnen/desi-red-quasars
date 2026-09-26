results.csv:

        File containing main line measurements for all processed spectra


relevant_results.csv:

        Combination results.csv and relevant further measurements from result_details.csv, filtered for good fits ('warning'==False)


result_details.csv:

        File containing fitting parameters and minor line measurements for all processed spectra


fit_params.csv:

        File containing just the fitting parameters from result_details.csv


civ_blueshift_vals.csv:

        File containing CIV blueshift values.


subset_desi_spectra.h5

        File containing subset of spectra for various testing purposes


target_properties.csv / target_properties_v2.csv:

        Files containing prefiltered spectra properties; result from SQL query


duplicate_targets.txt:

        File containing the 'targetid's of all objects that had more than one spectrum in SPARCL
