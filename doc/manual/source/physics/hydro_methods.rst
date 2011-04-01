.. _hydro_methods:

Hydro Methods
=============================

There are five possible routines in Enzo for calculating the evolution of the gas. 

Method 0: Piecewise Parabolic Method (PPM)
------------------------------------

Parameter file call: ``HydroMethod = 0``

Associated parameters: 

``RiemannSolver = 0``; specifies the type of solver xxxxx

``DuelEnergyFormalism``; allows the total and thermal energy to be followed seperately during the simulation. Helpful when the velocities are high such that E\ :sub:`total`\ >> E\ :sub:`thermal`. 

``PPMFlatteningParameter``;

``PPMSteepeningParamter``;

*Reference: Colella & Woodward (1992)*


fsdfjsisuruserusirlfsjdf;
Method 2: ZEUS
---------------

Parameter file call: ``HydroMethod = 2``

Associated parameters:

``ZEUSQuadraticArtificialViscosity``; 

``ZEUSLinearArtificialViscosity``;

*Reference: Stone & Norman*

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
