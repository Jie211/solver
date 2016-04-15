Iterative Krylov-Subspace Solvers
====
*-last update 16/4/15-*
## News
* 2016/4/12 Now Working on GMRES Method.
* 2016/3/20 Now Working on CUDA version Kskip-CG.

## Branches
* **cuda branche** is for parallel with CUDA.(now work at here)
* **makefile branche** is for make this solvers as a linux package by using GNUAutomake.(**last update 2015/12**)
* **cmdline branche** is a dev branche.(now not use)
* **master branche** is **not using** now.


## Description
Iterative Krylov-Subspace Solvers for large systems of linear equations. 
All the code is write be pure C.
This repository is part of my project.

## About Kskip method
NOLTA2015(English) [keynote](https://www.dropbox.com/s/ni0gt1m93izdhem/NOLTA2015_12_3.key?dl=0)
## About Variable Preconditioned method
Abe, Kuniyoshi, and Shao-Liang Zhang. "A variable preconditioning using the SOR method for GCR-like methods." Int. J. Numer. Anal. Model 2.2 (2005): 147-161.
Ikuno, Soichiro, et al. "Iterative solver for linear system obtained by edge element: variable preconditioned method with mixed precision on GPU." Magnetics, IEEE Transactions on 48.2 (2012): 467-470.

## Now Support method
* Conjugate Gradient method(CG).
* Conjugate Residual method(CR).
* Generalized Conjugate Residual method(GCR).
* Generalized Minimal REsidual method(GMRES).(Anyway, implementation done, meybe rewrite itlater.)
* Kskip method
	- Kskip-CG.(implementation done, but algorithm is not stability)
	- Kskip-CR.(implementation done, stil working...)
* Variable Preconditioned method(a.k.a VP)
	- VP Conjugate Gradient method(VPCG).
	- VP Conjugate Residual method(VPCR).
	- VP Generalized Conjugate Residual method(VPGCR).
  - VP Generalized Minimal REsidual method(VPGMRES).(Anyway, implementation done, meybe rewrite itlater.)

## Requirement
Hard
* CPU(of course...)
* NVIDIA GPU(capability version >= 2.x, architecture >= Kepler like NVIDIA GTX Titan)

Soft
* GCC compiler in CentOS need >= 4.4.7.(see CUDA Toolkit Documentation 1.1 System Requirements)
* CUDA Ver>6.0.
* OpenMP.

## Licence
The MIT License (MIT)

Copyright (c) 2015-2016 Jie211
## Author
[Jie211@github](https://github.com/Jie211)
[blog](https://www.jie211.me)
