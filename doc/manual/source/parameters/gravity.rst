Gravity Parameters
~~~~~~~~~~~~~~~~~~

General
^^^^^^^

``TopGridGravityBoundary`` (external)
    A single integer which specified the type of gravitational boundary
    conditions for the top grid. Possible values are 0 for periodic and
    1 for isolated (for all dimensions). The isolated boundary
    conditions have not been tested recently, so caveat emptor.
    Default: 0
``SelfGravity`` (external)
    This flag (1 - on, 0 - off) indicates if the baryons and particles
    undergo self-gravity.
``SelfGravityGasOff`` (external)
    This parameter is used in conjuction with SelfGravity so that only particles contribute to potential, not gas. Default = False (i.e. gas does contribute)
``GravitationalConstant`` (external)
    This is the gravitational constant to be used in code units. For cgs units it
    should be 4\*pi\*G. For cosmology, this value must be 1 for the
    standard units to hold. A more detailed decription can be found at :ref:`EnzoInternalUnits`. Default: 4\*pi.
``PotentialIterations`` (external)
    Number of iterations to solve the potential on the subgrids. Values
    less than 4 sometimes will result in slight overdensities on grid
    boundaries. Default: 4.
``MaximumGravityRefinementLevel`` (external)
    This is the lowest (most refined) depth that a gravitational
    acceleration field is computed. More refined levels interpolate
    from this level, provided a mechanism for instituting a minimum
    gravitational smoothing length. Default: ``MaximumRefinementLevel``
    (unless ``HydroMethod`` is ZEUS and radiative cooling is on, in which
    case it is ``MaximumRefinementLevel`` - 3).
``MaximumParticleRefinementLevel`` (external)
    This is the level at which the dark matter particle contribution to
    the gravity is smoothed. This works in an inefficient way (it
    actually smoothes the particle density onto the grid), and so is
    only intended for highly refined regions which are nearly
    completely baryon dominated. It is used to remove the discreteness
    effects of the few remaining dark matter particles. Not used if set
    to a value less than 0. Default: -1
``ParticleSubgridDepositMode`` (external)
    This parameter controls how particles stored in subgrid are deposited
    into the current grid.  Options are:
      0 (CIC_DEPOSIT) - This is a second-order, cloud-in-cell deposition
         method in which the cloud size is equal to the cell size in
         the target grid (particles are in source grid, deposited into
         target grid).  This method preserves the correct center-of-mass
         for a single particle but smears out boundaries and can result
         in small artifacts for smooth particle distributions (e.g.
         nested cosmological simulations with low perturbations).
      1 (CIC_DEPOSIT_SMALL) - This is also a CIC method, but the cloud
         size is taken to be the cell size in the source grid, so for
         subgrids, the cloud is smaller than the grid size.  This
         is an attempt to compromise between the other two methods.
      2 (NGP_DEPOSIT) - This uses a first order, nearest-grid-point
        method to deposit particle mass.  It does not preserve center-
        of mass position and so for single particle results in noisy
        accelerations.  However, it does correctly treat nested
        cosmology simulations with low initial perturbations.
     Default: 1
``BaryonSelfGravityApproximation`` (external)
    This flag indicates if baryon density is derived in a strange,
    expensive but self-consistent way (0 - off), or by a completely
    reasonable and much faster approximation (1 - on). This is an
    experiment gone wrong; leave on. Well, actually, it's important for
    very dense structures as when radiative cooling is turned on, so
    set to 0 if using many levels and radiative cooling is on [ignored
    in current version]. Default: 1

External Gravity Source
^^^^^^^^^^^^^^^^^^^^^^^

These parameters set up an external static background gravity source that is
added to the acceleration field for the baryons and particles.

``PointSourceGravity`` (external)
    This parameter indicates that there is to be a
    (constant) gravitational field with a point source profile (``PointSourceGravity`` =
    1) or NFW profile (``PointSourceGravity`` = 2). Default: 0
``PointSourceGravityConstant`` (external)
    If ``PointSourceGravity`` = 1, this is the magnitude of the point
    source acceleration at a distance of 1
    length unit (i.e. GM in code units). If ``PointSourceGravity`` =
    2, then it takes the mass of the dark matter halo in CGS
    units. ``ProblemType`` = 31 (galaxy disk simulation) automatically calculates
    values for ``PointSourceGravityConstant`` and
    ``PointSourceGravityCoreRadius``. Default: 1
``PointSourceGravityCoreRadius`` (external)
    For ``PointSourceGravity`` = 1, this is the radius inside which
    the acceleration field is smoothed in code units. With ``PointSourceGravity`` =
    2, it is the scale radius, rs, in CGS units (see Navarro, Frank & White,
    1997). Default: 0
``PointSourceGravityPosition`` (external)
    If the ``PointSourceGravity`` flag is turned on, this parameter
    specifies the center of the point-source gravitational field in
    code units. Default: 0 0 0
``ExternalGravity`` (external)
   This fulfills the same purpose as ``PointSourceGravity`` but is
   more aptly named. ``ExternalGravity = 1`` turns on an alternative
   implementation of the NFW profile with properties
   defined via the parameters ``HaloCentralDensity``, ``HaloConcentration`` and ``HaloVirialRadius``. Boxsize is assumed to be 1.0 in this case. ``ExternalGravity = 10`` gives a gravitational field defined by the logarithmic potential in Binney & Tremaine, corresponding to a disk with constant circular velocity.  Default: 0 
``ExternalGravityConstant`` (external)
    If ``ExternalGravity = 10``, this is the circular velocity of the disk in code units. Default: 0.0
``ExternalGravityDensity`` 
   Reserved for future use.
``ExternalGravityPosition`` (external)
    If ``ExternalGravity = 10``, this parameter specifies the center of the gravitational field in code units. Default: 0 0 0
``ExternalGravityOrientation`` (external)
    For ``ExternalGravity = 10``, this is the unit vector of the disk's angular momentum (e.g. a disk whose face-on view is oriented in the x-y plane would have ``ExternalGravityOrientation = 0 0 1``). Default: 0 0 0 
``ExternalGravityRadius`` (external)
   If ``ExternalGravity = 10``, this marks the inner radius of the disk in code units within which the velocity drops to zero. Default: 0.0
``UniformGravity`` (external)
    This flag (1 - on, 0 - off) indicates if there is to be a uniform
    gravitational field. Default: 0
``UniformGravityDirection`` (external)
    This integer is the direction of the uniform gravitational field: 0
    - along the x axis, 1 - y axis, 2 - z axis. Default: 0
``UniformGravityConstant`` (external)
    Magnitude (and sign) of the uniform gravitational acceleration.
    Default: 1
``DiskGravity`` (external)
    This flag (1 - on, 0 - off) indicates if there is to be a
    disk-like gravity field (Berkert 1995; Mori & Burkert 2000).  Default: 0
``DiskGravityPosition`` (external)
    This indicates the position of the center of the disk gravity.
    Default: 0 0 0
``DiskGravityAngularMomentum`` (external)
    Specifies the unit vector of the disk angular momentum.
    Default: 0 0 1
``DiskGravityStellarDiskMass`` (external)
    Total mass of stellar disk (in solar masses)
    Default: 1e11
``DiskGravityDiskScaleHeightR`` (external)
    Disk scale length in radius (in Mpc)
    Default: 4.0e-3
``DiskGravityDiskScaleHeightz`` (external)
    Disk scale height in z (in Mpc)
    Default: 2.5e-4
``DiskGravityStellarBulgeMass`` (external)
    Disk stellar bulge mass (in solar masses)
    Default: 1.0e10
``DiskGravityStellarBulgeR`` (external)
    Disk stellar bulge scalue radius (in Mpc)
    Default: 1.0e-4
``DiskGravityDarkMatterR`` (external)
    Dark matter halo scale radius (in Mpc)
    Default: 2.3e-2
``DiskGravityDarkMatterDensity`` (external)
    Dark matter effective density (in cgs)
    Default: 3.81323e-25
