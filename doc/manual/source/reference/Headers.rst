Header files in Enzo
====================

Here is a complete list of the Enzo 1.5 header files and a brief
description of what they do. I'm not 100% sure that this list is
actually going to be useful. I'm not filling it in completely, only
the things I need to be referencing in the tutorials.

[source:/public/trunk/src/enzo/AnalysisBaseClass.h
AnalysisBaseClass.h] *Needs to be filled in*

[source:/public/trunk/src/enzo/AnalyzeClusters.h AnalyzeClusters.h]
*Needs to be filled in*

[source:/public/trunk/src/enzo/CoolData.h CoolData.h]
*Needs to be filled in*

[source:/public/trunk/src/enzo/CosmologyParameters.h
CosmologyParameters.h] *Needs to be filled in*

[source:/public/trunk/src/enzo/ealFloat.h ealFloat.h]
*Needs to be filled in*

[source:/public/trunk/src/enzo/ealInt.h ealInt.h]
*Needs to be filled in*

[source:/public/trunk/src/enzo/EnzoArray.h EnzoArray.h]
*Needs to be filled in*

[source:/public/trunk/src/enzo/enzo\_unit\_tests.h
enzo\_unit\_tests.h] *Needs to be filled in*

[source:/public/trunk/src/enzo/error.h error.h] Houses one macro to
check and deal with MPI errors

[source:/public/trunk/src/enzo/ExternalBoundary.h
ExternalBoundary.h] The ExternalBoundary class definition.

[source:/public/trunk/src/enzo/FastSiblingLocator.h
FastSiblingLocator.h] Structure definitions for the chaining mesh
and sibling lists

[source:/public/trunk/src/enzo/flowdefs.h flowdefs.h]
*Needs to be filled in*

[source:/public/trunk/src/enzo/Fluxes.h Fluxes.h] The fluxes
structure, used to contain the Coarse and Refined fluxes for each
parent/subgrid pair.

[source:/public/trunk/src/enzo/global\_data.h global\_data.h] This
houses all global parameters for enzo, which is most of them.
Variables defined here are defined as extern in all routines but
[source:/public/trunk/src/enzo/enzo.C enzo.C] (see the
DEFINE\_STORAGE define there) and are initialized with
[source:/public/trunk/src/enzo/SetDefaultGlobalValues.C
SetDefaultGlobalValues.C]

[source:/public/trunk/src/enzo/Grid\_AnalyzeClusters.h
Grid\_AnalyzeClusters.h] *Needs to be filled in*

[source:/public/trunk/src/enzo/Grid.h Grid.h] This defines the
primary God Class, grid

[source:/public/trunk/src/enzo/GridList.h GridList.h]
*Needs to be filled in*

[source:/public/trunk/src/enzo/Hierarchy.h Hierarchy.h] Defines the
HierarchyEntry linked list structure. More can be found about this
in the `Linked List page? </wiki/Tutorials/LinkedLists>`_

[source:/public/trunk/src/enzo/ImplosionGlobalData.h
ImplosionGlobalData.h] *Needs to be filled in*

[source:/public/trunk/src/enzo/LevelHierarchy.h LevelHierarchy.h]
Defines the LevelHierarchyEntry linked list structure. More can be
found about this in the
`Linked List page? </wiki/Tutorials/LinkedLists>`_

[source:/public/trunk/src/enzo/ListOfParticles.h ListOfParticles.h]
*Needs to be filled in*

[source:/public/trunk/src/enzo/macros\_and\_parameters.h
macros\_and\_parameters.h] This is the home for all preprocessor
directives, and is responsible for overloading floating point
precision keywords.

[source:/public/trunk/src/enzo/message.h message.h]
*Needs to be filled in*

[source:/public/trunk/src/enzo/MTLPARAM.h MTLPARAM.h]
*Needs to be filled in*

[source:/public/trunk/src/enzo/performance.h performance.h]
*Needs to be filled in*

[source:/public/trunk/src/enzo/phys\_constants.h phys\_constants.h]
*Needs to be filled in*

[source:/public/trunk/src/enzo/ProtoSubgrid.h ProtoSubgrid.h]
Defines the ProtoSubgrid class, used in
[source:/public/trunk/src/enzo/FindSubgrids.C FindSubgrids.C]

[source:/public/trunk/src/enzo/RadiationFieldData.h
RadiationFieldData.h] *Needs to be filled in*

[source:/public/trunk/src/enzo/RateData.h RateData.h]
*Needs to be filled in*

[source:/public/trunk/src/enzo/region.h region.h]
*Needs to be filled in*

[source:/public/trunk/src/enzo/SedovBlastGlobalData.h
SedovBlastGlobalData.h] Problem specific data.

[source:/public/trunk/src/enzo/ShockPoolGlobalData.h
ShockPoolGlobalData.h] Problem specific data.

[source:/public/trunk/src/enzo/SphericalInfall.h SphericalInfall.h]
Problem specific data.

[source:/public/trunk/src/enzo/StarParticleData.h
StarParticleData.h] Problem specific data.

[source:/public/trunk/src/enzo/STD\_typedefs.h STD\_typedefs.h]
*Needs to be filled in*

[source:/public/trunk/src/enzo/TestGravitySphereGlobalData.h
TestGravitySphereGlobalData.h] Problem specific data.

[source:/public/trunk/src/enzo/TestProblemData.h TestProblemData.h]
*Needs to be filled in*

[source:/public/trunk/src/enzo/TopGridData.h TopGridData.h] Defines
the TopGrid structure, which houses the global parameters of the
simulation.

[source:/public/trunk/src/enzo/typedefs.h typedefs.h] Has all the
enumerate lists used to give words to parameters. Defines types for
field (density, etc), interpolation method, hydro method, boundary
type, gravity boundary type.

[source:/public/trunk/src/enzo/units.h units.h]
*Needs to be filled in*

[source:/public/trunk/src/enzo/WavePoolGlobalData.h
WavePoolGlobalData.h] *Needs to be filled in*


