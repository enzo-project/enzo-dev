.. _cosmic_ray_parameters:

Cosmic Ray Two-Fluid Model Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For details on the cosmic ray solver in Enzo see :ref:`cosmic_rays`.

``CRModel`` (external)
    This parameter turns on the model. Default: 0

    0. Off
    1. On

``CRgamma``
    For CR equation of state. Default: 4.0/3.0 (relativistic, adiabatic gas)

``CRDiffusion`` (external)
    Switches on diffusion of the cosmic ray energy density. Default: 0

    0. Off
    1. On with constant coefficient (``CRkappa``)


``CRkappa`` (external)
    Cosmic ray diffusion coefficient in CGS units (cm^2/s), Default: 0.0. For MW-like galaxies: 1E28.

``CRCourantSafetyNumber`` (external)
    Multiplies CR diffusion timestep, for stability should be <= 0.5. Default: 0.5

``CRFeedback`` (external)
    Specify fraction of star formation feedback energy should be diverted into the cosmic
    ray energy density. implemented ONLY for star_maker3 (feedback method 2). Default: 0.0

``CRdensFloor`` (external)
    Floor in gas density, can be imposed, for speed purposes (default 0.0 = off). Any value
    larger than 0.0 is on with that value as the floor in code units.
