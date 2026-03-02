the notebook 3Dplots.nb reads in the Fortran background data for 512 modes and produces interpolating functions which can then be plotted. 

How this works:

- copy the relevant .out files together with the shoot_inner.par file into this directory
- run bg_to_file.exe
- in notebook: check that the values for Nt, xcut0 and xcut1 are correct
- run notebook