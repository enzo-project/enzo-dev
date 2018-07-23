.. _radiative_transfer:

Radiative Transfer
==================

Enzo has two options for radiation transport: an adaptive ray-tracing
method, and an implicit flux-limited diffusion solver that is coupled
to Enzo's internal chemistry and cooling solvers.  Both are described
in more detail below.


Adaptive Ray Tracing
--------------------

Enzo includes a photon-conserving radiative transfer algorithm that
is based on an adaptive ray-tracing method utilizing the HEALPix
pixelization of a sphere (`Abel & Wandelt 2002 <http://adsabs.harvard.edu/abs/2002MNRAS.330L..53A>`_). Photons are integrated
outward from sources using an adaptive timestepping scheme that
preserves accuracy in ionization fronts even in the optically-thin
limit. This has been coupled to the chemistry and cooling network to
provide ionization and heating rates on a cell-by-cell basis, and
has the ability to follow multiple radiation groups, as well as
capturing H-minus and H2-photodissociating radiation as well as
hydrogen and helium-ionizing radiation. The
method is described in detail in `Wise & Abel (2011)
<http://adsabs.harvard.edu/abs/2011MNRAS.414.3458W>`_, and a listing
of parameters can be found at :ref:`radiative_transfer_ray_tracing`.
      

Flux-Limited Diffusion
----------------------

A second option for radiative transfer is a moment-based method that
adds an additional field tracking the radiation energy density. This
field is evolved using the flux-limited diffusion method, which
transitions smoothly between streaming (optically thin) and opaque
limits and is coupled to an ionization network of either purely
hydrogen, or both hydrogen and helium. The resulting set of linear
equations is solved using the parallel `HYPRE framework <https://computation.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods>`_. Full details
on the Enzo implementation of this method can be found in
`Reynolds et al. (2009)
<http://adsabs.harvard.edu/abs/2009JCoPh.228.6833R>`_, and a listing
of parameters can be found at :ref:`radiative_transfer_fld`.


A Practical Comparison of Methods
---------------------------------

Both the adaptive ray-tracing and flux-limited diffusion methods work
in both unigrid and adaptive mesh simulations.  In general, the
adaptive ray-tracing method provides a more accurate solution for
point-based radiation sources (i.e., it captures radiation shadowing
more accurately), but the computational cost scales roughly with the
number of sources.  The cost of the flux-limited diffusion solver, on
the other hand, has a cost that is independent of the number of
sources, which can make it more efficient for large-volume calculations.
