.. _cosmology:


Cosmology
=========

Enzo's fluid and gravity solvers are solved in a comoving coordinate
system, which allows computation in an expanding universe (i.e.,
cosmological expansion).  Other physics modules (such as cooling,
chemistry, and radiation) interface with this via Enzo's physical
units infrastructure, which converts comoving units to the proper (e.g.,
physical) units that are needed by those modules. 
This functionality is turned on by setting
``ComovingCoordinates = 1``, and further documentation can be found in
:ref:`cosmology_parameters`.
The source code for computing the expansion factor and its rate of change is in
*CosmologyComputeExpansionFactor.C*, and the code for computing
physical units can be found in *CosmologyGetUnits.C*.

At present, Enzo fully supports both flat (:math:`\Omega_k = 0`) and non-flat cosmologies with a
cosmological constant.  While it does not fully support models with
variable dark energy equations of state (e.g., quintessence models) or
other non-standard cosmologies, *CosmologyComputeExpansionFactor.C*
and related files can be modified to do so.  Note that Enzo's internal
units system for the expansion parameter (*a*) differs from the modern
standard; it is
set to 1 at the initial redshift of the simulation, rather than at the
present day (z=0).
