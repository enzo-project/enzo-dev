.. _cosmic_rays:

Cosmic Ray Two-Fluid Model
==========================
.. sectionauthor:: Munier Salem <msalem@astro.columbia.edu>
.. versionadded:: 2.2


For relevant parameters, please also see :ref:`cosmic_ray_parameters`.


*Source: Grid_ZeusSolver.C*


Models dynamical role of cosmic rays via a set of two-fluid hydro equations
(see `Jun et. al. 1994
<http://adsabs.harvard.edu/abs/1994ApJ...429..748J>`_ ). Central to the effort
is a new baryon field, CREnergyDensity, which is in units of ergs/cm^3, and is
advected along with the gas. Gradients in the CR field result in a pressure
felt by the gas. The CR gas is also diffusive and rays can be produced during
star formation. See :ref:`cosmic_ray_parameters` for information on how to control all
these options. But most important:



  - ``CRModel`` - Switches on the CR physics (0 = off, 1 = on)

  - ``CRgamma`` - For polytropic equation of state. 4/3 = relativistic, adiabatic gas (default)

  - ``CRDiffusion`` - turns on diffusion of CREnergyDensity field

  - ``CRkappa`` - Diffusion coefficient (currently constant, isotropic)

  - ``CRFeedback`` - Controls production of rays in star forming regions


For this model to run properly you *must be running the Zeus Hydro 
Solver:* ``HydroMethod = 2``. The model has not been implemented for
higher order solvers.


If you plan on including cosmic rays, definitely first verify the solver is working by running
the Cosmic Ray Shocktube problem, which ought to match the analytic solution described in
`Pfrommer 2006 <http://adsabs.harvard.edu/abs/2006MNRAS.367..113P>`_ . See the Test Problem 
Parameter list for more information.


Cosmic Rays have also been implemented in the isolated galaxy simulation. They initialize with
a profile identical to the thermal gas, multiplied by a constant, ``GalaxySimulationCR``, typically
set to 0.1.
