.. _tabulated_feedback:

Using Tabulated Feedback & Source Tracking
=================================================

By default, all stellar feedback routines use analytic prescriptions
for mass, metal, and energy yields. Tabular yields computed
from stellar population synthesis models such as the
`SYGMA <https://nugrid.github.io/NuPyCEE/overview.html>`_
chemical evolution code can be used instead.
These tables make it easy to track metals from specific feedback sources 
such as different types of supernovae or even pre-supernova feedback.
This feature is called "process tracking."

A subset of methods support the use of tabulated yields. 
See :doc:`../physics/star_particles`
for a current list of feedback methods and problem types
that support tabulated yields and source tracking respectively.

Feedback Table Structure
------------------------

The feedback tables are stored in HDF5 files with the following structure::

  /
    indexer/
      attributes:
        type2_index
        type1a_index
        agb_index
        massive_index
      datasets:
        initial_metal_fraction
        population_age
    sygma_models/
      attributes:
        *SYGMA model parameters*
      datasets:
        ejecta_mass
        ejecta_metal_mass
        sne_event_rate

Datasets in the ``indexer`` group act as keys for the first two dimensions 
of each table in the ``sygma_models`` group. The ``initial_metal_fraction``
key describes the first dimension and ``population_age`` the second.

The third dimension of the ``sygma_models`` datasets are the sources.
These are keyed by the ``indexer`` group attributes. For example, the
attribute ``type2_index`` should have a value of 0 and ``agb_index`` a value of 2.
**Enzo makes assumptions about the structure of these HDF5 tables so the
order of sources should not be changed.** Rather, these attributes
can be used as a reminder of how to index the third dimension of the yield tables.

The ``ejecta_mass`` and ``ejecta_metal_mass`` tables have total sizes
``N(initial_metal_fraction) * N(population_age) * 4``. 
The ``sne_event_rate`` table has size ``N(initial_metal_fraction) * N(population_age) * 2``.

The various parameters used to generate the table are stored as attributes in
the ``sygma_models`` group for posterity.

Adding Tabular Feedback in a New Feedback Routine
-------------------------------------------------

To add tabulated feedback to your desired feedback routine,
you'll need to modify both ``Grid_StarParticleHandler.C`` and
the file containing your feedback scheme. 

You'll want to check the state of the parameter ``StarFeedbackUseTabularYields``
either in ``Grid_StarParticleHandler.C`` or inside your feedback routine.
If this parameter is false, you should not be using tabular feedback!
Checking the state of this parameter in ``Grid_StarParticleHandler.C`` is
recommended as you won't need to pass an additional parameter
to your feedback routine. See the logic under ``if (STARFEED_METHOD(NORMAL_STAR))`` 
in ``Grid_StarParticleHandler.C`` for an example.

Changes to ``Grid_StarParticleHandler.C``
+++++++++++++++++++++++++++++++++++++++++

You'll need to pass the following pieces of grid-accessible data
to your feedback routine for the feedback tables to be integrated
properly:

* ``ParticleInitialMass`` (float pointer): array of initial particle masses
* ``FBTable.n_met`` (int): length of the table's initial metallicity dimension
* ``FBTable.n_age`` (int): length of the table's population age dimension
* ``FBTable.ini_met`` (float pointer): key for the table's initial metallicity dimension
* ``FBTable.pop_age`` (float pointer): key for the table's population age dimension
* ``FBTable.mass_yield`` (float pointer): table of mass yields
* ``FBTable.metal_yield`` (float pointer): table of metal yields
* ``FBTable.event_rate`` (float pointer): table of supernova event rates
* ``StarFeedbackTabularSNIIEnergy`` (float): parameter for the energy yield of a single Type II SN
* ``StarFeedbackTabularSNIaEnergy`` (float): parameter for the energy yield of a single Type Ia SN

Changes to Feedback Routine
+++++++++++++++++++++++++++

Your feedback routine will need to accept the data listed in the section above.

The Fortran subroutines for extracting and integrating yields from the
feedback tables are contained in a file called ``tabular_feedback.F``.
The subroutines contained within can be used as a module in any Fortran-based
feedback subroutine (e.g. ``subroutine star_feedback2_tab`` in ``star_feedback2_tab.F``)
or (theoretically) wrapped into a C(++)-based feedback function.
To use the ``tabular_feedback.F`` subroutines in your Fortran feedback subroutine,
add the line ``use tabular_feedback`` just under your subroutine declaration.
See ``star_feedback2_tab.F`` for an example.

The subroutines available for use in your feedback routine are as follows:

* ``sne_mass``: total mass yield for Type II and Ia supernovae
* ``sne_metal``: total metal yield for Type II and Ia supernovae
* ``sne_II_metal``: metal yield for only Type II supernovae
* ``sne_Ia_metal``: metal yield for only Typa Ia supernova
* ``sne_energy``: total energy yield for Type II and Ia supernovae
* ``AGB_mass``: mass yield for AGB winds
* ``AGB_metal``: metal yield for AGB winds

The ``mass`` and ``metal`` subroutines share the same call signature::

  routine(yield, 
          initial_particle_mass, particle_metal_frac, particle_age, 
          timestep_dt, time_units, 
          yield_table_pointer, table_key_pointer_metal, table_key_age_pointer,
          table_key_size_metal, table_key_size_age)

The ``energy`` subroutines have the following call signature::

  routine(yield, 
          initial_particle_mass, particle_metal_frac, particle_age, 
          timestep_dt, time_units, 
          yield_table_pointer, table_key_pointer_metal, table_key_age_pointer,
          table_key_size_metal, table_key_size_age)

The ``initial_particle_mass`` should be pulled from the ``ParticleInitialMass`` array.
Pointers and sizes related to the ``yield_table`` and ``table_key`` parameters should be
the ``FBTable`` members that were passed into your feedback routine from
``Grid_StarParticleHandler.C``.


Adding Source Tracking to a New Problem Type
--------------------------------------------

If you would like to use source tracking in a problem type not listed in 
:doc:`../physics/star_particles`, you'll need to follow steps 3 and 4 in
:doc:`HowToAddNewBaryonField` to add the following code
blocks::

  char *MetalIIName = "MetalSNII_Density";
  char *MetalIaName = "MetalSNIa_Density";
  char *MetalAGBName = "MetalAGB_Density";
  char *MetalNSMName = "MetalNSM_Density"

and ::

  if (StarMakerTypeIaSNe || StarFeedbackTrackMetalSources)
    DataLabel[count++] = MetalIaName;
  if (StarFeedbackTrackMetalSources) {
    DataLabel[count++] = MetalIIName;
    DataLabel[count++] = MetalAGBName;
    DataLabel[count++] = MetalNSMName;
  }

Make sure to follow the order in which these fields were added to ``Grid_InitializeUniformGrid.C``!