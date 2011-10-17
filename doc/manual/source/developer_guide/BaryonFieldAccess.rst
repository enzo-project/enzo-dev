.. _BaryonFieldAccess:

Accessing Data in BaryonField
=============================

For performance reasons, Enzo uses Fortran source to do all the
important work. Because of this, it doesn't use the standard C/C++
data structure for the 3D ``BaryonField`` array, which stores all the
Eulerian data.

``BaryonField`` is stored as a one dimensional array. Typically C/C++
data is stored in row major order. **ENZO DATA IS STORED IN COLUMN
MAJOR ORDER** because of its Fortran underpinnings.

To map between one and three dimensions, in column major order, use
the following:

.. code-block:: c

    OneDindex = i + nx*j + nx*ny*k

in Enzo grid member functions, this can be done like this:

.. code-block:: c

    index = i + GridDimension[0]*(j + GridDimension[1]*k);

It should also be mentioned that it is always important to access
data in 'stride 1' order. That means accessing data in the order it
is stored in memory. So to set all ``BaryonFields`` to the number
12345.6:

.. code-block:: c

    int index;
    for(int field=0;field<NumberOfBaryonFields;field++)
    for(int k=0;k<GridDimension[2]; k++)
      for(int j=0;j<GridDimension[1]; j++)
        for(int i=0;i<GridDimension[0]; i++){
          index = i + GridDimension[0]*(j + GridDimension[1]*k);
          BaryonField[field][index] = 12345.6;
        }

This loops over the ghost zones as well as the active zones. To
loop over only active zones, use ``GridStartIndex`` and ``GridEndIndex``.
Note that this loop must include ``GridEndIndex``

.. code-block:: c

    int index;
    for(int field=0;field<NumberOfBaryonFields;field++)
    for(int k=GridStartIndex[2];k<=GridEndIndex[2]; k++)
      for(int j=GridStartIndex[1];j<=GridEndIndex[1]; j++)
        for(int i=GridStartIndex[0];i<=GridEndIndex[0]; i++){
          index = i + GridDimension[0]*(j + GridDimension[1]*k);
          BaryonField[field][index] = 12345.6;
        }


