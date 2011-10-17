Header files in Enzo
====================

Here is a complete list of the Enzo 2.0 header files and a brief
description of what they do.

``src/enzo/CoolData.h`` 

  Contains parameters for cooling tables and radiation fields.  Most
  importantly this struct has the pointers to the tabulated cooling
  functions that are used in ``cool1d_multi.src``.  This type is used
  for the global variable CoolData.

``src/enzo/CosmologyParameters.h`` 

  Defines the global variables that are used in cosmology
  simulations, e.g. cosmological parameters, initial redshift,
  redshift outputs.

``src/enzo/ealFloat.h`` 

  Class for floating-point arrays that supports array arithmetic.
  Mainly used by the Enzo Analysis class.

``src/enzo/ealInt.h`` 

  Same as ``ealFloat.h`` but for integers.

``src/enzo/EnzoArray.h`` 

  Templated class that is a container for grid and particle quantities
  in the Enzo Analysis class.

``src/enzo/enzo_unit_tests.h`` 
  
  Framework for simple tests on Enzo.  Not used in typical
  simulations.

``src/enzo/ExternalBoundary.h`` 

  The ExternalBoundary class definition.

``src/enzo/FastSiblingLocator.h`` 

  Structure definitions for the chaining mesh and sibling lists.

``src/enzo/flowdefs.h`` 

  Function prototypes and variables for FLOW_TRACE define.  Currently
  not used.

``src/enzo/Fluxes.h`` 

  The fluxes structure, used to contain the Coarse and Refined fluxes
  for each parent/subgrid pair.

``src/enzo/global_data.h`` 

  This houses all global parameters for Enzo, which is most of them.
  Variables defined here are defined as extern in all routines but
  ``src/enzo/enzo.C`` (see the ``DEFINE_STORAGE`` #define there) and
  are initialized with ``src/enzo/SetDefaultGlobalValues.C``.

``src/enzo/Grid.h`` 

  This defines the primary God Class, ``grid``.

``src/enzo/GridList.h`` 

  Structure for a linked list of grids.  Used when identifying new
  subgrids, ``Grid_IdentifyNewSubgrids.C`` and
  ``Grid_IdentifyNewSubgridsSmall.C``.

``src/enzo/Hierarchy.h`` 

  Defines the HierarchyEntry linked list structure. More can be found
  about this in :doc:`LinkedLists`.

``src/enzo/ImplosionGlobalData.h`` 

  Contains global variables that have store the parameters in the
  Implosion problem type.

``src/enzo/LevelHierarchy.h`` 

  Defines the ``LevelHierarchyEntry`` linked list structure. More can
  be found about this in :doc:`LinkedLists`.

``src/enzo/ListOfParticles.h`` 

  Structure for a linked list of particle lists.  Used in
  ``OutputAsParticleData.C``.

``src/enzo/macros_and_parameters.h`` 

  This is the home for all preprocessor directives, and is responsible
  for overloading floating point precision keywords.

``src/enzo/message.h`` 

  Defines to handle error, warning, and debug messages.

``src/enzo/MTLPARAM.h`` 

  Common variables for the Cen's metal cooling routines,
  ``mcooling.src``

``src/enzo/performance.h`` 

  Defines for the interface between Enzo and LCAperf.

``src/enzo/phys_constants.h`` 

  Defines for physical constants

``src/enzo/ProtoSubgrid.h`` 

  Defines the ProtoSubgrid class, used in ``src/enzo/FindSubgrids.C``.

``src/enzo/RadiationFieldData.h`` 

  Structure that contains the parameters and variables that describe
  the background radiation field.  Only used for the global variable
  RadiationData in ``global_data.h``.

``src/enzo/RateData.h`` 

  Structure that holds all of the parameters and arrays of the rate
  equations for the non-equilibrium chemistry.  Only used for the
  global variable RateData.

``src/enzo/region.h`` 

  Structures that describe a region when computing the parallel FFT.

``src/enzo/SedovBlastGlobalData.h`` 

  Contains global variables that have store the parameters in the
  Sedov blast problem type.

``src/enzo/ShockPoolGlobalData.h`` 

  Contains global variables that have store the parameters in the
  shock pool problem type.

``src/enzo/SphericalInfall.h`` 

  Contains global variables that have store the parameters in the
  spherical infall problem type.

``src/enzo/StarParticleData.h`` 

  Global variables that store parameters about the star formation
  routines.  It also has variables that keep track of the number of
  stars.

``src/enzo/TestGravitySphereGlobalData.h`` 

  Contains global variables that have store the parameters in the test
  gravity sphere problem type.

``src/enzo/TestProblemData.h`` 

  Structure that stores parameters that describe a problem
  initialization.

``src/enzo/TopGridData.h`` 

  Defines the TopGrid structure, which houses the global parameters of
  the simulation.

``src/enzo/typedefs.h`` 

  Has all the enumerate lists used to give words to
  parameters. Defines types for field (density, etc), interpolation
  method, hydro method, boundary type, gravity boundary type.

``src/enzo/units.h`` 

  Global variables that store the units in CGS.  Used when
  ComovingCoordinates is *off*.

``src/enzo/WavePoolGlobalData.h`` 

  Contains global variables that have store the parameters in the wave
  pool problem type.


