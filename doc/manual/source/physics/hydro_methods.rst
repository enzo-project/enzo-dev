.. _hydro_methods:

Hydro and MHD Methods
=====================

There are four available methods in Enzo for calculating the evolution
of the gas with and without magnetic fields. Below is a brief
description of each method, including the parameters associated with
each one and a link to further reading. 
For relevant parameters please also see :ref:`hydrodynamics_parameters`.

Additionally, there are two MHD methods, which are described in detail in :ref:`mhd_methods`

Method 0: Piecewise Parabolic Method (PPM)
------------------------------------------
*Source:  Grid_SolvePPM_DE.C*

The PPM scheme uses a parabolic function to estimate the left and
right states of the Godunov problem. This method has a third-order
accurate piecewise parabolic monotonic interpolation and a nonlinear
Riemann solver for shock capturing. It does an excellent job of
capturing strong shocks across a few cells. This more accurately
represents both smooth gradients and discontinuities over linear
interpolation, i.e. PLM (piecewise linear method). See
:ref:`hydrodynamics_parameters` for more details about parameters.

Parameters
^^^^^^^^^^

Main call: ``HydroMethod = 0``

``RiemannSolver``: specifies the type of solver, where the following
only works with the PPM solver.

1. HLL (Harten-Lax-van Leer) a two-wave, three-state solver with no
   resolution of contact waves.  This is the most diffusive of the
   available three solvers in PPM.  *New for version 2.1*

4. HLLC (Harten-Lax-van Leer with Contact) a three-wave, four-state
   solver with better resolution of contacts.  The most resilient to
   rarefaction waves (e.g. blastwave interiors). *New for version 2.1*

5. **Default** Two-shock approximation.  Iterative solver.

``RiemannSolverFallback``: allows for the Riemann solver to "fallback"
to the more diffusive HLL solver when negative energies or densities
are computed.  Only applicable when using the HLLC and Two-shock
solvers.  The fluxes in the failing cell are recomputed and used in
the Euler update of the gas quantities. *New for version 2.1*

``ConservativeReconstruction``: When interpolating (PPM) to the left
and right states, interpolation occurs in the conserved variables
(density, momentum, and energy) instead of the primitive variables
(density, velocity, and pressure).  This results in more accurate
results in unigrid simulations but can cause errors with AMR.  See
Section 4.2.2 (steps 1-5) and Appendices A1 and B1 in Stone et
al. (2008, ApJS 178, 137).  *New for version 2.1*

``DualEnergyFormalism``: allows the total and thermal energy to be
followed seperately during the simulation. Helpful when the velocities
are high such that E\ :sub:`total`\ >> E\ :sub:`thermal`.

``PPMFlatteningParameter``

``PPMSteepeningParameter``

Links
^^^^^

\ P. R. Woodward and P. Colella. "A piecewise parabolic method for gas
dynamical simulations," *J. Comp. Phys*, 54:174, 1984 `link
<https://seesar.lbl.gov/anag/publications/colella/A_1_4_1984.pdf>`__


Method 2: ZEUS
--------------
*Source: ZeusSource.C, Zeus_xTransport.C, Zeus_yTransport.C,
Zeus_zTransport.C, Grid_ZeusSolver.C, ZeusUtilities.C*

ZEUS is a finite-difference method of solving hyperbolic PDEs instead
of solving the Godunov problem. This method uses hydrodynamical
algorithm originally used in ZEUS, but the MHD and radiation
hydrodynamics schemes are not implemented from ZEUS. This method is
formally second-accurate in space but first-order accurate in time. It
is a very robust but relatively diffusive scheme.

Parameters
^^^^^^^^^^

Main call: ``HydroMethod = 2``

``ZEUSQuadraticArtificialViscosity``

``ZEUSLinearArtificialViscosity`` 


Links
^^^^^

\ J. M. Stone and M. L. Norman. "Zeus-2D: A radiation
magnetohydrodynamics code for astrophysical flows in two space
dimensions. I. The hydrodynamics algorithms and tests."  *The
Astrophysical Journal Supplement*, 80:753, 1992 `link
<http://adsabs.harvard.edu/abs/1992ApJS...80..753S>`__

\ J. M. Stone and M. L. Norman. "Zeus-2D: A radiation
magnetohydrodynamics code for astrophysical flows in two space
dimensions. II. The magnetohydrodynamic algorithms and tests." *The
Astrophysical Journal Supplement*, 80:791, 1992 `link
<http://adsabs.harvard.edu/abs/1992ApJS...80..791S>`__

Method 3: MUSCL
---------------

.. versionadded:: 2.0

The MUSCL [#f1]_ scheme is a second-order accurate extensive of Godunov's
method for solving the hydrodynamics in one dimension.  The
implementation in Enzo uses second-order Runge-Kutta time
integration.  In principle, it can use any number of Riemann solvers
and interpolation schemes.  Here we list the compatible ones that are
currently implemented.

Parameters
^^^^^^^^^^
Parameter file call: ``HydroMethod = 3``

``RiemannSolver``: specifies the type of solver, where the following
only works with the MUSCL solver.

1. HLL (Harten-Lax-van Leer): a two-wave, three-state solver with no
   resolution of contact waves.

3. LLF (Local Lax-Friedrichs) is based on central differences instead
   of a Riemann problem.  It requires no characteristic information.
   This is the most diffusive of the available three solvers in
   MUSCL.

4. HLLC (Harten-Lax-van Leer with Contact): a three-wave, four-state
   solver with better resolution of contacts.  The most resilient to
   rarefaction waves (e.g. blastwave interiors).

If negative energies or densities are computed, the solution is
corrected using a more diffusive solver, where the order in decreasing
accuracy is HLLC -> HLL -> LLF.

``ReconstructionMethod``: specifies the type of interpolation scheme
used for the left and right states in the Riemann problem.

0. PLM: **default**

Method 4: MHD with Hyperbolic Cleaning (Dedner)
-----------------------------------------------

The two MHD methods in Enzo differ primarily in the mechanism for maintaining
:math:`\nabla \cdot B = 0`.  
These are described in more detail in :ref:`mhd_methods`.

Parameters
^^^^^^^^^^

``HydroMethod = 4`` uses the hyperbolic cleaning method of Dedner et
al. (2002, JCP 175, 645).  The basic integration is the MUSCL 2nd
order Runga Kutta method described above. This class of solvers has
been ported to nVidia's CUDA framework.  As ``HydroMethod = 3``, there
are three Riemann solver options, though instead of HLLC, HLLD is
available

1. HLL (Harten-Lax-van Leer): a two-wave, three-state solver with no
   resolution of contact waves.

3. LLF (Local Lax-Friedrichs) is based on central differences instead
   of a Riemann problem.  It requires no characteristic information.
   This is the most diffusive of the available three solvers in
   MUSCL.

6. HLLD (Harten-Lax-van Leer with Discontinuities): a 5-wave, six-state
   solver.  HLLD includes two fast waves, two Alfven waves, and one contact
   discontinuity.  

``ReconstructionMethod``: specifies the type of interpolation scheme
used for the left and right states in the Riemann problem.

0. PLM: **default**

``UsePoissonDivergenceCleaning``:
Enables additional divergence cleaning by solving a Poisson equation.
This works on top of the standard mixed hyperbolic/parabolic divergence
leaning
and is in most cases not required.
Works on indiviual grids, i.e., it's *not* a global divergence purge.
Use with care as this feature is not extensively tested.

Default: 0 (off)

Please see  for all relevant parameters, see :ref:`mhd_dender_parameters`.


Links
^^^^^

\ Dedner et al. "Hyperbolic Divergence Cleaning for the MHD
Equations,"
*Journal of Computational Physics*, 175, 645, 2002 `link
<https://https://ui.adsabs.harvard.edu/#abs/2010ApJS..186..308C/abstract>`__

Method 6: MHD with Constrained Transport (CT)
---------------------------------------------

``HydroMethod = 6`` uses the CT method, which computes an electric field from
the Riemann solver, then uses that electric field to update the magnetic field.
This MHD method is second-order in space and timee, and preserves
the divergence constraint, ∇ · B = 0, to machine precision through
the Constrained Transport (CT) method (Collins et al. 2010)

Links
^^^^^

\ Collins et al. "Cosmological Adaptive Mesh Refinement
Magnetohydrodynamics with Enzo,"
*The Astrophysical Journal Supplement*, 186:308, 2010 `link
<https://https://ui.adsabs.harvard.edu/#abs/2010ApJS..186..308C/abstract>`__

Parameters
^^^^^^^^^^
Parameter file call: ``HydroMethod = 6``

Method 5: No Hydro
------------------

.. versionadded:: 2.0

For testing non-hydro machinery in Enzo, one can turn hydro off.

Parameters
^^^^^^^^^^
Parameter file call: ``HydroMethod = 5``

Notes
-----

``HydroMethod = 1`` was an experimental implementation that is now
obsolete, which is why it is skipped in the above notes.

.. rubric:: Footnotes

.. [#f1] Monotone Upstream-centered Schemes for Conservation Laws
