# me573_final_cfd_solver
Finite Difference Solver for Lid Driven (Shear Driven) Cavity Flow

University of Wisconsin Madison ME 573 - Computational Fluid Dynamics
12/13/2019
Author: Piotr Liebersbach

Final Projecct for UW Madison Mechanical Engineering Course on Computational Fluid Dynamics.  CFD solver simulates a lid driven cavity flow at a Reynolds Number Re = 100.
Solution computed using transient explicit finite difference method on a staggared grid.  Navier Stokes equations are discretized using operator splitting method, and each time step is updated as follows:
  1. Intermediate Velocity is calculated
  2. Presure is computed by solving the Pressure Poisson Equation using Intermediate Velocity
  3. Intedmediate velocity is corrected using newly computed pressure

Solution validated against published data from:
  U. GHIA, K. N. GHIA, AND C. T. SHIN (1982) "High-Re Solutions for Incompressible Flow Using the Navier-Stokes Equations and a Multigrid Method,"
  JOURNAL OF COMPUTATIONAL PHYSICS 48, 387-411
  https://doi.org/10.1016/0021-9991(82)90058-4

Runs in Matlab (and Octave although slower) and has no external dependencies.
