# Computing critical spacetimes for scalar field collapse
This repository contains a `Fortran` implementation used to compute critical solutions of spherically symmetric D-dimensional general relativity coupled to a massless scalar field as well as aspects of linear perturbations around them. It was used to generate data for the publications 

- Angle of null energy condition lines in critical spacetimes ([Phys. Rev. D112 (2025) L021505](https://journals.aps.org/prd/abstract/10.1103/vpyn-b2fn))
- Critical spacetime crystals in continuous dimensions ([SciPost submission](https://arxiv.org/abs/2602.10185))

The code is a reworked version of the program by Jose M. Martin-Garcia and Carsten Gundlach used in [Phys.Rev. D68 (2003) 024011](https://arxiv.org/abs/gr-qc/0304070). In short, the main new features are 

- Generalized to arbitrary real values of the dimension D.
- Included Taylor expansions up to order 5 around the center and order 3 around the self-similarity horizon. 
- Implemented adaptive stepsize algorithm for generating optimal grids in the x-direction. (as an alternative to a fixed logarithmic grid)
- Perturbations can be computed at lower t-resolution than the background (see `tstep` in parameterfile)
- Implemented bisection + Brent algorithm for faster root finding in the perturbative computation

For a detailed description of the numerical procedure please refer to the [SciPost submission](https://arxiv.org/abs/2602.10185). 
A performance-optimized and parallelized `C++` implementation may be found in [tobjec/parallel-critical-collapse](https://github.com/tobjec/parallel-critical-collapse).

## Background and echoing period
First, build the executable by running 


 `background/source/` 

## Linear Perturbations and critical exponent

## Plot scripts and postprocessing
