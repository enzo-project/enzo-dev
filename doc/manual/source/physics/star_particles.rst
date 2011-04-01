Active Particles: Stars, BH, and Sinks
======================================

There are many different subgrid models of star formation and feedback
in the astrophysical literature, and we have included several of them
in Enzo.  There are also methods include routines for black hole,
sink, and Pop III stellar tracer formation.  Here we give the details
of each implementation and the parameters that control them.

Method 1: Cen & Ostriker
------------------------
*Source: star_maker2.F*

This routine uses the algorithm from Cen & Ostriker (1992, ApJL 399,
113) that creates star particles when the following six criteria are
met

#. The gas density is greater than the threshold set in the parameter
   ``StarMakerOverDensityThreshold``.  This parameter is in code units
   (i.e. overdensity with respect to the mean matter density)

#. The divergence is negative

#. The dynamical time is less than the cooling time or the temperature
   is less than 11,000 K.

#. The cell is Jeans unstable.

#. The star particle mass is greater than ``StarMakerMinimumMass``,
   which is in units of solar masses.

#. The cell does not have finer refinement underneath it.

These particles add thermal and momentum feedback to the grid cell
that contains it until 12 dynamical times after its creation.  In each
timestep,

.. math::
   
   M_{\rm form} &= M_0 [ (1+x_1) \exp(-x_1) - (1+x_2) \exp(-x_2) ]\\
   x_1 &= (t - t_0) / t_{\rm dyn}\\
   x_2 &= (t + dt - t_0) / t_{\rm dyn}

of stars are formed, where M\ :sub:`0` and t\ :sub:`0` are the initial
star particle mass and creation time, respectively.  

* M\ :sub:`ej` = M\ :sub:`form` * ``StarMakerEjectionFraction`` of gas
  are returned to the grid and removed from the particle.

* M\ :sub:`ej` * v\ :sub:`particle` of momentum are added to the cell.

* M\ :sub:`form` * c\ :sup:`2` * ``StarMakerEnergyToThermalFeedback``
  of energy is deposited into the cell.

* M\ :sub:`form` * ((1 - Z\ :sub:`star`) * ``StarMetalYield`` + M\
  :sub:`ej` * Z\ :sub:`star`) of metals are added to the cell.  This
  formulation accounts for gas recycling back into the stars.

Method 2: Cen & Ostriker with Stochastic Star Formation
-------------------------------------------------------
*Source: star_maker3.F*

Method 3: Global Schmidt Law
----------------------------
*Source: star_maker4.F*

Method 4: Population III Stars
------------------------------
*Source: pop3_maker.F*

Method 5: Sink particles
------------------------
*Source: sink_maker.C*

Method 6: Radiative Stellar Clusters
------------------------------------
*Source: cluster_maker.F*

Method 7: Cen & Ostriker with no delay in formation
---------------------------------------------------
*Source: star_maker7.F*

Method 8: Springel & Hernquist
------------------------------
*Source: star_maker5.F*

Method 9: Massive Black Holes
-----------------------------
*Source: mbh_maker.C*

Method 10: Population III stellar tracers
-----------------------------------------
*Source: pop3_color_maker.F*

Notes
------------------------

The routines included in ``star_maker1.F`` are obsolete and not
compiled into the executable.  For a more stable version of the
algorithm, use Method 1.

