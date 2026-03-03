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

For a detailed description of the numerical procedure please refer to section 4 in the [SciPost submission](https://arxiv.org/abs/2602.10185). 
 A performance-optimized and parallelized `C++` implementation may be found in [tobjec/parallel-critical-collapse](https://github.com/tobjec/parallel-critical-collapse).

## Background and echoing period
First, build the executable by running 

```make shoot_inner.exe```

in background/source/ (need to install gfortran first). The maximal values for nx and ntau (here called ny) are alredy fixed at compile time.
nxmax can be set in inner/shootmain_inner.f and nymax is set in `nymax.inc`.  

Now, copy the executable to wherever your initial data guess is saved and run it there. 
Required files for correct execution: 

+ fc.dat
+ psic.dat
+ Up.dat
+ Delta.dat
+ shoot_inner.par

The first three contain an initial guess for the boundary data at the center and the SSH and Delta.dat is the initial guess for the echoing period. 
The parameter file contains the following options:

+ `ny` number of points in the tau direction. Should be a power of two for the FFT to work properly. 
+ `dim` spacetime dimension
+ `xleft` position of left cutoff surface (around the center)
+ `xmid` position of matching surface
+ `xright` position of left cutoff surface (around the SSH)
+ `eps` finite difference approximation for the Jacobian
+ `errmax` precision goal for mismatch Newton algorithm
+ `verbose` additional command line output info at execution
+ `slowerr` damping for Newton algorithm (for better stability)
+ `outevery` stepsize in x for writing output to files (only needed for plotting)
+ `useloggrid` if T, uses fixed logarithmid grid in x. if F, uses adaptive stepsize algorithm to generate grid data (ignores values for xmid, nleft, nright). The matching surface in the latter case is put at the value of x where Max(Psi) = Max(Psic)/2
+  `nleft` number of x-points left to xmid
+  `nright` number of x-points right to xmid
+  `tol` local error tolerance of adaptive stepsize (ignored if useloggrid=T)
+  `tstep` can be 1 or 2. Only relevant for reading out fields to generate perturbations (see later)
+  `crit` precision of IRK2 solver in x-direction
+  `debug` additional output of data for debugging

After convergence, the boundary data corresponding to the critical solution is saved to `fc.out`, `psic.out`, `Up.out` and `Delta.out`. 

## Linear Perturbations and critical exponent
For computing the critical exponent \lambda one needs two steps. First, the background fields need to be written out (over x and tau) and then linear perturbations are computed to find the single unstable mode. 

### 1. Read out background

### 2. Run perturbation code


## Plot scripts and postprocessing
