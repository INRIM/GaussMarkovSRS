# Usage

## Input data

gauss_markov_srs requires 3 types of inputs: reference values, input data and correlations between input data.

### Reference values
Reference frequency values should be specified in a file similar to:

    #Id	Atom    	nu0
    0	133Cs   	1.
    1	115In+  	1267402452901041.3
    2	1H      	1233030706593514.
    3	199Hg   	1128575290808154.32
    ...

where columns are an integer id, a string labelling the reference and a string with an approximate absolute frequency value. The first line specifies the vale for caesium, whose value is not used. This file can be loaded as a numpy structured array with:
    
    reference = gm.load_txt("reference.txt")

### Input values
Input values should be specified in a file similar to:

    #Id	Ref         	Atom1   	Atom2	Sup	            Value	            Unc 	Note
    1	vonZanthier2000	115In+  	133Cs		            1267402452899920.0	690.0	Uncertx3
    2	Wang2007	    115In+  	133Cs		            1267402452901265.0	256.0	
    3	Ohtsubo2017	    115In+  	133Cs		            1267402452901049.8	7.5	
    4	Parthey2011	    1H  	    133Cs		            1233030706593517.5	5.0	
    ...
    84	Riedel2020  	171Yb+E3	87Sr	NPL/SYRTE-TW	1.495991618544900976	4.94e-16
    85	Riedel2020  	171Yb+E3	87Sr	NPL/SYRTE-PPP	1.495991618544901113	4.04e-16
    86	Riedel2020  	171Yb+E3	87Sr	NPL/PTB-TW  	1.495991618544900644	5.24e-16
    ...


whit columns:
* `Id` for integer ids
* `Ref` for strings specifying a reference
* `Atom1` and `Atom2` for numerator and denominator of the measured frequency ratios
* `Sup` for strings with supplementary information in case multiple ratios reported in a single reference (can be empty)
* `Value` with a string representing the measured value
* `Unc` with a string representing the uncertainty of each measurements (absolute, not relative)
* `Note` for comments, this field is unused by gauss_markov_srs

This file can be loaded as a numpy structured array from text or excel files, for example:

    input_data = gm.load_excel("data.xlsx")

### Input correlations
Correlations between input data should be specified in a file similar to:

    #Id1	Id2	Corr
    3   	7	0.001
    3   	25	0.004
    3      	43	0.001
    3   	48	0.011
    3   	49	0.010
    3   	50	0.002
    3   	56	0.005
    ...

whit columns:
* `Id1` and `Id2` for the integer ids of the measurements in the input value file.
* `Corr` for the correlation value 

This file can be loaded as a numpy structured array, for example:

    input_corr = gm.load_txt("correlations.txt")

## Preprocessing
Input data, loaded as numpy structured arrays, can be converted to an array of fractional frequency values `y` and their uncertainties `yu` with:

    y, yu = gm.get_y(input_data, reference)

The model matrix can be obtained with:

    M = gm.get_model_matrix(input_data, reference)

Finally the corrlation and covariance matrices of the input data can be calculated with:

    cyy = gm.get_corr_matrix(input_data, input_corr)
    Cyy = gm.corr2cov(cyy, yu)

## Running the fit
The fit proper can be run as:

    fit_results = gm.gauss_markov_fit(M, y, Cyy)

Extra stats can be calculated with:    
    
    stats = gm.fit_stats(M, y, Cyy, fit_results)

## Outputs
Results can be printed with

    gm.pretty_print_fit_result(fit_results, reference)
    gm.pretty_print_stats(stats)

saved in a directory with:
    
    gm.save(dir, label, reference, fit_results, stats)

and plotted to a directory with
    
    gm.plot_residuals(figdir, label, input_data, reference, fit_results, stats)


