.. _hydro_methods:

Hydro and MHD Methods
=============================

There are five possible routines in Enzo for calculating the evolution of the gas with and without magnetic fields. Below is a brief description of each method, including the parameters associated with each one and a link to further reading. 

Method 0: Piecewise Parabolic Method (PPM)
------------------------------------
*Source:  euler.F, euler_sweep.h, feuler_sweep.F*

Brief description here.

Parameters
^^^^^^^^^^

Main call: ``HydroMethod = 0``

``RiemannSolver = 0``: specifies the type of solver xxxxx

``DuelEnergyFormalism``: allows the total and thermal energy to be followed seperately during the simulation. Helpful when the velocities are high such that E\ :sub:`total`\ >> E\ :sub:`thermal`. 

``PPMFlatteningParameter``

``PPMSteepeningParamter``


Links
^^^^^^

\ P. R. Woodward and P. Colella. "A piecewise parabolic method for gas dynamical simulations," *J. Comp. Phys*, 54:174, 1984 `link <https://seesar.lbl.gov/anag/publications/colella/A_1_4_1984.pdf>`_


Method 2: ZEUS
---------------
*Source: ZeusSource.C, Zeus_xTransport.C, Zeus_yTransport.C, Zeus_zTransport.C, Grid_ZeusSolver.C, ZeusUtilities.C*

Description here.

Parameters
^^^^^^^^^^

Main call: ``HydroMethod = 2``

``ZEUSQuadraticArtificialViscosity``

``ZEUSLinearArtificialViscosity`` 


Links
^^^^^^

\  J. M. Stone and M. L. Norman. "Zeus-2D: A radiation magnetohydrodynamics code for astrophysical flows in two space dimensions. I. The hydrodynamics algorithms and tests."  *The Astrophysical Journal Supplement*, 80:753, 1992 `link <http://adsabs.harvard.edu/abs/1992ApJS...80..753S>`_

\ J. M. Stone and M. L. Norman. "Zeus-2D: A radiation magnetohydrodynamics code for astrophysical flows in two space dimensions. II. The magnetohydrodynamic algorithms and tests." *The Astrophysical Journal Supplement*, 80:791, 1992 `link <http://adsabs.harvard.edu/abs/1992ApJS...80..791S>`_

Method 3: MUSCL
---------------

Parameter file call: ``HydroMethod = 3``

Method 4: ???
---------------

Parameter file call: ``HydroMethod = 4``

.. raw:: html
   
   <font color="red">Only available in the unstable release</font>

Method 5: ???
---------------

Parameter file call: ``HydroMethod = 5``

.. raw:: html
   
   <font color="red">Only available in the unstable release</font>


Notes
------

``HydroMethod = 1`` was an experimental implementation that is now obsolute, which is why it is skipped in the above notes.
