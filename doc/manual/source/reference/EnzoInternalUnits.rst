.. _EnzoInternalUnits:

Enzo Internal Unit System
=========================

The units of the physical quantities used in Enzo depend on the problem being
run. For most test problems there is no physical length or time specified, so
the units can be be simply scaled. For cosmology, there are a set of units
designed to make most quantities of order unity so that single precision
floating-point variables can be used. These units are defined in
:ref:`EnzoOutputFormats`.  Additionally, discussion of how particle masses are
stored in Enzo can be found at :ref:`EnzoParticleMass`.  However, with the
broader use of Enzo for non-cosmological astrophysics applications, it has
become necessary to add a new set of units into the code. This page describes
how to set these units.

In order to have a self-consistent set of units, the user has to set
appropriate length, time, and mass OR density scales.  Simulations that include
gravity also need to have a self-consistent gravitational constant that is
scaled to the other variables. The four parameters that the user can set are
``LengthUnits``, ``TimeUnits``, ``DensityUnits``, and ``MassUnits``. Only one of ``DensityUnits``
or ``MassUnits`` needs to be set, since ``MassUnits`` = ``DensityUnits`` * ``LengthUnits``
:sup:`3` . Additionally, if the parameter ``SelfGravity`` is turned on (set to 1),
the parameter ``GravitationalConstant`` must be set to 4\*pi\*G, where G is
Newton's gravitational constant as a dimensionless quantity (that is, with all
units scaled out).

The primary motivation for using a non-arbitrary set of units is to take
advantage of Enzo's various chemistry and cooling algorithms, some of which
have been scaled assuming CGS units. To do this, one chooses physical units
assuming the simulation box size is unity in code units, and that a
density/mass and time value of 1 in code units are something meaningful in CGS.
For example, if one is interested in setting the box size to one parsec, a
density of 1 in code units equivalent to to 10 hydrogen atoms per cubic
centimeter (in CGS units), and the time unit to one million years, the
appropriate settings of the parameters would be as follows:

::

    DensityUnits = 1.67e-23    # 10 hydrogen atoms per cc in CGS (c/cm^3)
    LengthUnits = 3.0857e+18   # one parsec in cm
    TimeUnits = 3.1557e+13     # one megayear in seconds

If we then wish to use gravity, the gravitational constant must be set
explicitly to 4\*pi\*G expressed in a unitless fashion. Since the gravitational
constant in CGS has units of cm\ :sup:`3`\ /(g\*s\ :sup:`2`\ ), this means that
the value should be 4\*pi\*G\ :sub:`cgs` \* ``DensityUnits`` * ``TimeUnits`` :sup:`2`\ . So,
in the units expressed above, that means the gravity parameters must be set as
follows:

::

    SelfGravity                = 1
    GravitationalConstant      = 0.0139394         # 4*pi*G_{cgs}*DensityUnits*TimeUnits^2

Note that if gravity is turned on, the parameter ``TopGridGravityBoundary`` must
also be set to either 0 (periodic) or 1 (isolated).

If you set only ``LengthUnits`` and ``DensityUnits`` but not ``TimeUnits`` the code will
calculate it for you using the actual gravitational constant. You see it
printed out in the terminal when the code starts up and you can also find it
towards the end of the parameter file of any output.

