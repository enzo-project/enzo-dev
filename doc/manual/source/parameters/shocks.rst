.. _shock_finding_parameters:

Shock Finding Parameters
~~~~~~~~~~~~~~~~~~~~~~~~

For details on shock finding in Enzo see :ref:`shock_finding`.

``ShockMethod`` (external)
    This parameter controls the use and type of shock finding. Default: 0
    
    ::

	0 - Off
	1 - Temperature Dimensionally Unsplit Jumps
	2 - Temperature Dimensionally Split Jumps
	1 - Velocity Dimensionally Unsplit Jumps
	2 - Velocity Dimensionally Split Jumps

``ShockTemperatureFloor`` (external)
    When calculating the mach number using temperature jumps, set the
    temperature floor in the calculation to this value.

``StorePreShockFields`` (external)
    Optionally store the Pre-shock Density and Temperature during data output.

``FindShocksOnlyOnOutput`` (external)
    0: Finds shocks during Evolve Level and just before writing out data. 1: Only find shocks just before writing out data.  2: Only find shocks during EvolveLevel. Default: 0
