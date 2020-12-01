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

Included are figures of the solution:
  1. ke.png - Plots the kinetic energy of the fluid as a function of time.  Shows the flow reaching a steady state
  2. pressure.png - Plots the pressure contour of the flow at final time step
  3. u_velocity.png - Plots a contour of the U component of the velocity at the final time step
  4. u_velocity_profile.png - Plots U component velocity against Y location at a constant x = 0.5 as a line, and plots the validation data as points to compare against
  5. v_velocity.png - Plots a contour of the V component of the velocity at the final time step
  6. v_velocity.png - Plots V component velocity against X location at a constant y = 0.5 as a line, and plots the validation data as points to compare against
  7. velocity_magnitude.png - Plots contour of the velocity magnitude with the velocity vectors at each grid points, at the final time step

Runs in Matlab (and Octave although slower) and has no external dependencies.
