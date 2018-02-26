.. _cosmology-parameters:

Cosmology Parameters
~~~~~~~~~~~~~~~~~~~~

``ComovingCoordinates`` (external)
    Flag (1 - on, 0 - off) that determines if comoving coordinates are
    used or not. In practice this turns on or off the entire cosmology
    machinery. Default: 0
``CosmologyFinalRedshift`` (external)
    This parameter specifies the redshift when the calculation will
    halt. Default: 0.0
``CosmologyOmegaMatterNow`` (external)
    This is the contribution of all non-relativistic matter (including
    HDM) to the energy density at the current epoch (z=0), relative to
    the value required to marginally close the universe. It includes
    dark and baryonic matter. Default: 0.279
``CosmologyOmegaLambdaNow`` (external)
    This is the contribution of the cosmological constant to the energy
    density at the current epoch, in the same units as above. Default:
    0.721
``CosmologyOmegaRadiationNow`` (external)
    This is the contribution of all relativistic matter to the energy
    density at the current epoch (z=0), in the same units as above.
    Default: 0.0.
``CosmologyHubbleConstantNow`` (external)
    The Hubble constant at z=0, in units of 100 km/s/Mpc. Default:
    0.701
``CosmologyComovingBoxSize`` (external)
    The size of the volume to be simulated in Mpc/h (at z=0). Default:
    64.0
``CosmologyInitialRedshift`` (external)
    The redshift for which the initial conditions are to be generated.
    Default: 20.0
``CosmologyMaxExpansionRate`` (external)
    This float controls the timestep so that cosmological terms are
    accurate followed. The timestep is constrained so that the relative
    change in the expansion factor in a step is less than this value.
    Default: 0.01
``CosmologyTableNumberOfBins`` (external)
    Conversions between time and redshift are computed by interpolating
    from a numerically integrated table of log(scale factor) vs. time.
    This parameter sets the number of bins in the table. Default: 1000.
``CosmologyTableLogaInitial`` (external)
    This sets the lower bound of the table used to convert between time
    and redshift. This is log10 of the lowest value of the scale factor.
    This value will be automatically adjusted if
    ``CosmologyInitialRedshift`` is set to an earlier time.
    Default: -6.0, (i.e., z = 999,999.)
``CosmologyTableLogaFinal`` (external)
    This sets the upper bound of the table used to convert between time
    and redshift. This is log10 of the highest value of the scale factor.
    This value will be automatically adjusted if
    ``CosmologyFinalRedshift`` is set to a later time.
    Default: 0.0, (i.e., z = 0.)
``CosmologyCurrentRedshift`` (information only)
    This is not strictly speaking a parameter since it is never
    interpreted and is only meant to provide information to the user.
    Default: n/a

