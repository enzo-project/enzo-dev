.. _gravity:


Gravity
======================================

The current implementation of self-gravity in Enzo uses a fast Fourier
technique (Hockney & Eastwood 1988) to solve Poisson’s equation on the
root grid on each timestep. The advantage of using this method is that
it is fast, accurate, and naturally allows both periodic and isolated
boundary conditions for the gravity, choices which are very common in
astrophysics and cosmology. On subgrids, we interpolate the boundary
conditions from the parent grid (either the root grid or some other
subgrid). The Poisson equation is then solved on every timestep using
a multigrid technique on one subgrid at a time. Aside from
self-consistently calculating the gravitational potential arising from
the baryon fields and particles in the simulation, there are also a
number of options for specifying static gravitational fields.  Enzo
parameters relating to gravity can be found in
:ref:`gravity-parameters` .
