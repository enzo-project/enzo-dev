.. _shock_finding:

Shock Finding
==================
.. sectionauthor:: Sam Skillman <samskillman@gmail.com>
.. versionadded:: 2.1

For relevant parameters, please also see :ref:`shock_finding_parameters`.

*Source: Grid_FindShocks.C*

Shock finding is accomplished using one of four methods.  The primary
method uses a coordinate-unsplit temperature jump (method 1), as described in
`Skillman et. al. 2008
<http://adsabs.harvard.edu/abs/2008ApJ...689.1063S>`_ with the
exception that instead of searching across multiple grids for the pre-
and post-shock cells, we terminate the search at the edge of the ghost
zones within each grid.  

Shock finding is controlled by the ``ShockMethod`` parameter, which
can take the following values:

0 - Off

1 - Unsplit Temperature Jumps

2 - Dimensionally Split Temperature Jumps

3 - Unsplit Velocity Jumps

4 - Dimensionally Split Velocity Jumps

When ``ShockMethod`` nonzero, this will create a "Mach" field in the
output files.   

Note: Method 1 has been used the most by the developer, and therefore
is the primary method.  Method 2 has been tested quite a bit, but the
downsides of using a dimensionally split method are outlined in the
above paper.  Methods 3 and 4 are more experimental and will run, but
results may vary.

Additional Shock Finding Parameters:

``ShockTemperatureFloor`` - When calculating the mach number using temperature jumps, set the temperature floor in the calculation to this value.

``StorePreShockFields`` - Optionally store the Pre-shock Density and Temperature during data output.





