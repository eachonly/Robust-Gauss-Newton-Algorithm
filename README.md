# Robust Gauss-Newton Algorithm
This repository contains the Robust Gauss-Newton (RGN) algorithm developed by Youwei Qin, Dmitri Kavetski and George Kuczera. 

When using RGN please cite the following articles:

Qin Y, Kavetski D, Kuczera G (2018) A robust Gauss-Newton algorithm for the optimization of hydrological models: From standard Gauss-Newton to robust Gauss-Newton. Water Resources Research, 54. https://doi.org/10.1029/2017WR022488

Qin Y, Kavetski D, Kuczera G (2018) A robust Gauss-Newton algorithm for the optimization of hydrological models: Benchmarking against industry-standard algorithms. Water Resources Research, 54. https://doi.org/10.1029/2017WR022489

# Robust Gauss-Newton Algorithm Description
The Robust Gauss-Newton (RGN) algorithm is designed for solving optimization problems with a sum of least squares objective function. The RGN algorithm introduces three heuristics to enhance its performance: (i) the large sampling scale scheme to capture large-scale features of the objective function, (ii) the best-sampling point scheme to take advantage of free information, and (iii) the null-space jump scheme to escape near-flat regions.

This repository includes two examples to illustrate the application of the RGN algorithm: optimisation of a 2D Rosenbrock function and calibration of the 5 parameter hydrological model HYMOD. The following folders are included:

  - SCR_RGN: the RGN algorithm (rgn.f90) and an auxiliary module (constantsMod.f90)
  - SCR_DEMO\rosenbrock: driver code (rgnMain_Rosenbrock.f90)
  - PROJ\rosenbrock:vfproj files, batches files, makefiles, input/output files for rosenbrock example
  - SCR_DEMO\hymod: driver code (rgnMain_Hymod.f90) and HYMOD model code
  - PROJ\hymod:vfproj files, batches files, makefiles, input/output files for hymod example
  - SLN: the sln files from Visual Studio

This repository contains the RGN algorithm implementation in Fortran-95, which has been tested using the GNU gfortran and Intel Fortran compilers, on Windows and Linux.
