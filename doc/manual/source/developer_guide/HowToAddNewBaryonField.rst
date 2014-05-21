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

3. You can access the field in your problem initializer or elsewhere using the ``FindField`` function. To get the field number, you would use

.. code-block:: c

   int NaNum;
   NaNum = FindField(SodiumDensity, FieldType, NumberOfBaryonFields);

Now you can access the field as ``BaryonField[NaNum]``. For example, to set the value to 0 everywhere,

.. code-block:: c

   for (int i = 0; i < size; i++)
      BaryonField[NaNum][i] = 0.0;

Conservation
------------

By default, Enzo assumes that fields represent densities for purposes of advection and interpolation. If your field is not a density field, you will
need to make some adjustments to make sure that the field is properly conserved. To do this, you can modify the macros in ``typedefs.h`` under ``FieldTypeIsDensity``.
