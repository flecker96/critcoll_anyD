This repository contains the code used to compute critical solutions of spherically symmetric D-dimensional general relativity coupled to a massless scalar field as well as aspects of linear perturbations around them. 

It is a reworked version of the program by Jose M. Martin-Garcia and Carsten Gundlach used in the work [Phys.Rev. D68 (2003) 024011](https://arxiv.org/abs/gr-qc/0304070). In short, the main new features of the current program are 

- Generalized to arbitrary real values of the dimension D.
- Included Taylor expansions up to order 5 around the center and order 3 around the self-similarity horizon. 
- Implemented adaptive stepsize algorithm for generating optimal grids in the x-direction. (as an alternative to a fixed logarithmic grid)
- Perturbations can be computed at lower t-resolution than the background (see parameterfile `tstep`)
- Implemented bisection + Brent algorithm for faster root finding in the perturbative computation



.) background/

contains the routine for computing the background solution

.) back_to_pert/

routine for converting the output of the background routine into a format for
the input of the perturbation routine

.) perturbations/

code for extracting the exponent of the single unstable linear perturbation mode

.) Mathematica/

contains general derivations of the EMKG system and preparation of the EOM for the 
Fortran code.
