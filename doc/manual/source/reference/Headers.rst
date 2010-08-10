Header files in Enzo
====================

Here is a complete list of the Enzo 2.0 header files and a brief
description of what they do.

- ``src/enzo/AnalysisBaseClass.h`` *Needs to be filled in*

- ``src/enzo/AnalyzeClusters.h`` *Needs to be filled in*

- ``src/enzo/CoolData.h`` *Needs to be filled in*

- ``src/enzo/CosmologyParameters.h`` *Needs to be filled in*

- ``src/enzo/ealFloat.h`` *Needs to be filled in*

- ``src/enzo/ealInt.h`` *Needs to be filled in*

- ``src/enzo/EnzoArray.h`` *Needs to be filled in*

- ``src/enzo/enzo\_unit\_tests.h`` *Needs to be filled in*

- ``src/enzo/error.h`` Houses one macro to
check and deal with MPI errors

- ``src/enzo/ExternalBoundary.h`` The ExternalBoundary class definition.

- ``src/enzo/FastSiblingLocator.h`` Structure definitions for the chaining mesh
and sibling lists

- ``src/enzo/flowdefs.h`` *Needs to be filled in*

- ``src/enzo/Fluxes.h`` The fluxes
structure, used to contain the Coarse and Refined fluxes for each
parent/subgrid pair.

- ``src/enzo/global\_data.h`` This
houses all global parameters for Enzo, which is most of them.
Variables defined here are defined as extern in all routines but
``src/enzo/enzo.C`` (see the
``DEFINE\_STORAGE`` #define there) and are initialized with
``src/enzo/SetDefaultGlobalValues.C``.

- ``src/enzo/Grid\_AnalyzeClusters.h``  *Needs to be filled in*

- ``src/enzo/Grid.h`` This defines the
primary God Class, grid

- ``src/enzo/GridList.h`` *Needs to be filled in*

- ``src/enzo/Hierarchy.h`` Defines the
HierarchyEntry linked list structure. More can be found about this
in the `Linked List page? </wiki/Tutorials/LinkedLists>`_

- ``src/enzo/ImplosionGlobalData.h`` *Needs to be filled in*

- ``src/enzo/LevelHierarchy.h`` Defines the ``LevelHierarchyEntry``
linked list structure. More can be
found about this in the
`Linked List page? </wiki/Tutorials/LinkedLists>`_

- ``src/enzo/ListOfParticles.h`` *Needs to be filled in*

- ``src/enzo/macros\_and\_parameters.h`` This is the home for all preprocessor
directives, and is responsible for overloading floating point
precision keywords.

- ``src/enzo/message.h`` *Needs to be filled in*

- ``src/enzo/MTLPARAM.h`` *Needs to be filled in*

- ``src/enzo/performance.h`` *Needs to be filled in*

- ``src/enzo/phys\_constants.h`` *Needs to be filled in*

- ``src/enzo/ProtoSubgrid.h`` Defines the ProtoSubgrid class, used in
``src/enzo/FindSubgrids.C``.

- ``src/enzo/RadiationFieldData.h`` *Needs to be filled in*

- ``src/enzo/RateData.h`` *Needs to be filled in*

- ``src/enzo/region.h`` *Needs to be filled in*

- ``src/enzo/SedovBlastGlobalData.h`` Problem specific data.

- ``src/enzo/ShockPoolGlobalData.h`` Problem specific data.

- ``src/enzo/SphericalInfall.h`` Problem specific data.

- ``src/enzo/StarParticleData.h`` Problem specific data.

- ``src/enzo/STD\_typedefs.h`` *Needs to be filled in*

- ``src/enzo/TestGravitySphereGlobalData.h`` Problem specific data.

- ``src/enzo/TestProblemData.h`` *Needs to be filled in*

- ``src/enzo/TopGridData.h`` Defines
the TopGrid structure, which houses the global parameters of the
simulation.

- ``src/enzo/typedefs.h`` Has all the
enumerate lists used to give words to parameters. Defines types for
field (density, etc), interpolation method, hydro method, boundary
type, gravity boundary type.

- ``src/enzo/units.h`` *Needs to be filled in*

- ``src/enzo/WavePoolGlobalData.h`` *Needs to be filled in*


