How to add a new baryon field
=============================

If you wish to add a new field to Enzo, there are a few files that you
will need to modify:

1. Add the field to ``typedefs.h`` in the ``field_type`` data structure. To do this, look for the last field in the list before ``FieldUndefined``. Give your problem a new number and increment ``FieldUndefined``. For example, let's say you wanted to add a new field called ``SodiumDensity``. If the last field in field_types is ``FieldUndefined  = 96``, you would add your field as

.. code-block:: c

   SodiumDensity = 96,
   FieldUndefined = 97;

2. Next, you need to add your field to ``Grid_InitializeUniformGrid.C``. At the top of the file you need to declare an ``int`` to hold the number which is used to reference the field, for example ``NaNum`` for the ``SodiumDensity`` example. Further down, you need to add your field to the FieldType list. After the other fields have been added, add the line

.. code-block:: c

   FieldType[NaNum    = NumberOfBaryonFields++] = SodiumDensity;

In theory, you could add the fields and allocate the fields in your problem initializer code (as some test problems currently in Enzo do), but it is cleaner and simpler to have your problem initializer call ``InitializeUniformGrid`` before doing the setup for your problem type. For more details, see 
:ref:`AddingANewTestProblem`.

3. Finally, you need to modify the initializer of problem types using your new field to make sure that the field is written out. Add the lines

.. code-block:: c

   char* SodiumName = "Sodium_Density";
   ...
   DataLabel[i++] = SodiumName;

after the other DataLabels are set. Note that you need to set the Data Labels in the same order that the fields were added in Grid_InitializeUniformGrid.C or the fields will be written incorrectly.

4. You can access the field in your problem initializer or elsewhere using the ``FindField`` function. To get the field number, you would use

.. code-block:: c

   int NaNum;
   NaNum = FindField(SodiumDensity, FieldType, NumberOfBaryonFields);

Now you can access the field as ``BaryonField[NaNum]``. For example, to set the value to 0 everywhere,

.. code-block:: c

   for (int i = 0; i < size; i++)
      BaryonField[NaNum][i] = 0.0;

For a more detailed discussion of how data in BaryonFields is accessed, see
:ref:`BaryonFieldAccess`.

Conservation
------------

For the purpose of advection and interpolation, Enzo assumes that all fields are densities unless told otherwise. If your field is not a density field, you will need to make some adjustments to make sure that the field is properly conserved. To do this, you can modify the macros in ``typedefs.h`` under ``FieldTypeIsDensity``. Non-density fields will be multiplied by density prior to flux correction and converted back afterwards. This process will make the field be conserved in the same way as density fields. To see how Enzo decides whether a field needs to be multiplied by density, take a look at the file ``MakeFieldConservative.C``. The actual manipulation is done in the flux correction and interpolation routines, and should not need to be modified.
