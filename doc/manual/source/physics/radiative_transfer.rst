.. _radiative_transfer:

Radiative Transfer
==================
.. sectionauthor:: John Wise <jwise@physics.gatech.edu>
.. versionadded:: 2.0

For relevant parameters, please also see :ref:`radiative_transfer_ray_tracing` and :ref:`radiative_transfer_fld`.


Adaptive Ray Tracing
--------------------

Solving the radiative transfer equation can be computed with adaptive
ray tracing that is fully coupled with the hydrodynamics and energy /
rate solvers.  The adaptive ray tracing uses the algorithm of Abel &
Wandelt (2002) that is based on the HEALPix framework.

For the time being, a detailed description and test suite can be found
in the paper Wise & Abel (2011, MNRAS 414, 3458).

Flux Limited Diffusion
----------------------

More details can be found in the paper Reynolds et al. (2009, Journal
of Computational Physics 228, 6833).
