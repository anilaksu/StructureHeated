# StructureHeated
This repository will include Fortran 90 codes developed for the structural deformation under heating 
In the end, I hope it will solve a heated block with finite deformation which will be performed with spectral methods

1) Heat Transfer problem is solved, it will be coupled with linear elasticity problem.
2) Ways to include (Spectral Exterior Calculus) SPEX formulation, which is explained in Rfat Zhelil Phd thesis, 
will be tried to be incorporated into the system to capture the finite deformation. 
3) It will start in a relatively small domain which requires less subgrids, but it will be extended to larger grids. As a result, 
parallelized linear solvers will be developed.
