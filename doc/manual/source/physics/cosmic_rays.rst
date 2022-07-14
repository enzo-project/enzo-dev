.. _cosmic_rays:

Cosmic Ray Two-Fluid Model
==========================


This section documents the two-fluid cosmic ray model implemented in
Enzo, which was first used (and is described in detail) in `Salem &
Bryan 2014 <http://adsabs.harvard.edu/abs/2014MNRAS.437.3312S>`_ . 
The implementation of anisotropic cosmic ray physics is described and tested
in `Butsky and Quinn 2018 <https://ui.adsabs.harvard.edu/abs/2018ApJ...868..108B/abstract>`_.
For relevant parameters, please also see
:ref:`cosmic_ray_two_fluid_model_parameters`.  The bulk of the advection code itself can be found in *Grid_ZeusSolver.C* and *hydro_rk/Riemann_LLF_MHD.C*. Cosmic ray transport is solved in *Grid_ComputeCRDiffusion.C*, *Grid_ComputeAnisotropicCRDiffusion.C*, and *Grid_ComputeCRStreaming.C*. 


This module models the dynamical role of cosmic rays via a set of two-fluid hydro equations
(see `Jun et. al. 1994
<http://adsabs.harvard.edu/abs/1994ApJ...429..748J>`_ ). Central to the effort
is a new baryon field, CREnergyDensity, which is in units of ergs/cm^3, and is
advected along with the gas. Gradients in the CR field result in a pressure
felt by the gas. The CR gas is also diffusive and rays can be produced during
star formation. See :ref:`cosmic_ray_two_fluid_model_parameters` for information on how to control all
these options. But most important are:


  - ``CRModel`` - Switches on the CR physics (0 = off, 1 = on). Compatible with HydroMethod = 2 or 4.

  - ``CRgamma`` - For polytropic equation of state. 4/3 = relativistic, adiabatic gas (default)

  - ``CRDiffusion`` - Turns on diffusion of CREnergyDensity field (0 = off, 1 = isotropic, 2 = anisotropic). CRDiffusion = 2 is only Compatible with HydroMethod = 4. 

  - ``CRkappa`` - Diffusion coefficient (currently constant, isotropic)

  - ``CRStreaming`` - Turns on CR streaming along mangetic field lines (0 = off, 1 = on). Compatible with HydroMethod = 4. 

  - ``CRFeedback`` - Controls production of rays in star forming regions (0 = off, 1 = on) 


For this model to run properly you *must be running either the Zeus Hydro 
Solver!* (``HydroMethod = 2``) or the Dedner MHD Solver (``HydroMethod = 4``). The model has not yet been implemented for
any of the other fluid solvers in Enzo.

If you plan on including cosmic rays, definitely first verify the solver is working by running
the Cosmic Ray Shocktube problem, which ought to match the analytic solution described in
`Pfrommer 2006 <http://adsabs.harvard.edu/abs/2006MNRAS.367..113P>`_
. See :ref:`cr_shocktube_param` for more detailed information on this
test problem.

Cosmic Rays have also been implemented in the isolated galaxy simulation. They initialize with
a profile equal to the density of the thermal gas, multiplied by a constant, ``GalaxySimulationCR``, typically
set to 0.1 (all in code units).
