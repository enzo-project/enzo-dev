Active Particles: Stars, BH, and Sinks
======================================

There are many different subgrid models of star formation and feedback
in the astrophysical literature, and we have included several of them
in Enzo.  Here we give the details of each implementation and the
parameters that control them.

Method 1: Cen & Ostriker
------------------------
*Source: star_maker2.F*

Method 2: Cen & Ostriker with Stochastic Star Formation
-------------------------------------------------------
*Source: star_maker3.F*

Method 3: Global Schmidt Low
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

