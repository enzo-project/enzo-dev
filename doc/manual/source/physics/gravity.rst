.. _gravity:


Gravity
=======


The Multigrid solver
--------------------

The default implementation of self-gravity in Enzo uses a fast Fourier
technique (`Hockney & Eastwood 1988 <http://adsabs.harvard.edu/abs/1988csup.book.....H>`_)
to solve Poisson’s equation on the
root grid on each timestep. The advantage of using this method is that
it is fast, accurate, and naturally allows both periodic and isolated
boundary conditions for the gravity, choices which are very common in
astrophysics and cosmology (with isolated boundary conditions on the
root grid being implemented with the `James (1977) method
<https://doi.org/10.1016/0021-9991(77)90013-4>`_).
On subgrids, we interpolate the boundary
conditions from the parent grid (either the root grid or some other
subgrid). The Poisson equation is then solved on every timestep using
a multigrid technique on one subgrid at a time. Aside from
self-consistently calculating the gravitational potential arising from
the baryon fields and particles in the simulation, there are also a
number of options for specifying static gravitational fields
(including, for example, gravitational acceleration from NFW halos,
galactic disks, and point sources).  Enzo
parameters relating to gravity can be found in
:ref:`gravity_parameters`, and a brief description .

The APM solver
--------------

Self-gravity can also be solved the Adaptive Particle-Mesh (APM) technique from
`Passy & Bryan 2014 <https://ui.adsabs.harvard.edu/abs/2014ApJS..215....8P/abstract>`_.
The general idea is to split the gravitational force between a long-range component
and one or more short-range components that are non-zero only for a narrow range of wavenumbers.
More details on the algorithm can be found in the paper above.
Enzo parameters related to the APM solver are listed and briefly described in :ref:`gravity_parameters`.