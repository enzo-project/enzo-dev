Grid Field Arrays
=================

Field arrays are convenient ways (within code linked against the
Enzo code base -- including within Enzo itself!) to access grid
data such as the baryon fields, or particle lists. They can also be
used to get pre-defined derived fields, such as temperature. They
are intended to be used by solvers, initializers, and analysis
routines. The hope is provide a clean way for classes other than
the grid to get to grid data, and to help make the current code
base more modular.

Class Description
-----------------

The array class is pretty simple: just enough to represent an
N-dimensional grid, without any spatial information. Here is the
heart of it, from ``EnzoArray.h``:

.. code-block:: c

    template<typename T>
    class EnzoArray
    {
    public:
    
      EnzoArray(int rank, int *dims, int *start, int *end,
            FLOAT *cell_size=NULL, int derived=FALSE){
    
    ...
    
      int Rank;                        // number of dimensions
      int Dimension[MAX_DIMENSION];    // total dimensions of all grids
      int StartIndex[MAX_DIMENSION];   // starting index of the active region
                                       //   (zero based)
      int EndIndex[MAX_DIMENSION];     // stoping index of the active region
                                       //   (zero based)
      FLOAT CellWidth[MAX_DIMENSION];
      
      T *Array;
    
      // used for velocities and positions
      T *Vector[MAX_NUMBER_OF_PARTICLE_ATTRIBUTES];
    
    ...
    };
    
    #define EnzoArrayFLOAT EnzoArray<FLOAT>
    #define EnzoArrayFloat EnzoArray<float>
    #define EnzoArrayInt   EnzoArray<int>

The array classes are really a single template, but the macros at
the bottom of the header file will hide that from you.

Array vs. Vector
~~~~~~~~~~~~~~~~

In the above code block, you'll notice two pointers: ``T \*Array``; and
``T \*Vector``. Here are the rules that these attributes follow: Only
one of these will be used, and which one is used depends on the
type of data you try to access. Namely, field data, such as
density, will be pointed to by ``Array``, and vector data, such as
velocities or particle positions, will be pointed to by ``Vector``.

Destructor (What Gets Deleted)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When the destructor is called, ``Array`` and ``Vector`` get deleted
*only if* ``derived`` is TRUE. This is to keep the usage (declare and
delete) similar for both derived and underived data. We really
don't want to delete the density field on accident.

Access Methods
--------------

There are six accessor methods declared in
``Grid.h``, two per data type
(``float``, ``int``, and ``FLOAT``).

.. code-block:: c

       EnzoArrayInt *CreateFieldArrayInt(field_type field);
       EnzoArrayInt *CreateFieldArrayInt(char *field_name);
      
       EnzoArrayFloat *CreateFieldArrayFloat(field_type field);
       EnzoArrayFloat *CreateFieldArrayFloat(char *field_name);
      
       EnzoArrayFLOAT *CreateFieldArrayFLOAT(field_type field);
       EnzoArrayFLOAT *CreateFieldArrayFLOAT(char *field_name);

These methods are defined in
``Grid_CreateFieldArray.C``.
Basically, they allocate a new
``EnzoArray``, fill in the dimensions, attach the relevant pointers,
and hand it back to. All you need to do is delete the return
object.

Field Numbers and Names
-----------------------

The arguments to are either a field number, defined in
``typedefs.h``, or the
string version of the same. The string versions are defined in a
``long`` array, named ``field_map`` in
``Grid_CreateFieldArray.C``.
This means you can access something as

.. code-block:: c

       EnzoArrayFloat *density_array = mygrid->CreateFieldArrayFloat(Density);

or

.. code-block:: c

       EnzoArrayFloat *density_array = mygrid->CreateFieldArrayFloat("Density");

There are some fields which have names that are the same as grid
attributes, like ``ParticlePosition``. Rather than have a huge
namespace conflict, these have field numbers prefixed with a "g",
e.g., ``gParticlePosition``. The string called is still just
"ParticlePosition", like

.. code-block:: c

       EnzoArrayFloat *ppos = mygrid->CreateFieldArrayFloat(gParticlePosition);

or

.. code-block:: c

       EnzoArrayFloat *ppos = mygrid->CreateFieldArrayFloat("ParticlePosition");

The important part of the map is that it knows the data type of the
fields, which you need to know, so you can call the right method.
This is really pretty simple, since just about everything returned
is a ``float``. For a complete list of the (hopefully current) fields,
see the section Field_List_Reference_. For the best reference,
check in ``typedefs.h``,
and ``Grid_CreateFieldArray.C``.

Using the Methods
-----------------

Here's a somewhat long-winded example of how to use the arrays.
First, here's function to create a non-uniform grid

.. code-block:: c

    grid *Linear3DGrid(){
      // Create a new 3D grid                                                                                                        
      float dens = M_PI, total_energy = 0.5, internal_energy = 0.0;
      float vel[3];
      int dims[3];
      FLOAT left[3], right[3];
    
      grid *lineargrid = new grid;
      int i, j, k, rank = 3;
      int index;
    
      for (i = 0; i < rank; i++) {
        dims[i] = 134;
        left[i] = 0.0;
        right[i] = 1.0;
        vel[i] = (i+1) * 0.125;
      }
    
      NumberOfParticleAttributes = 0;
      lineargrid->PrepareGrid(3, dims,
                              left, right, 2);
    
      int result = lineargrid->InitializeUniformGrid(dens, total_energy, internal_energy, vel);
      assert(result != FAIL);
    
      EnzoArrayFloat *dens_field = lineargrid->CreateFieldArrayFloat("Density");
    
      for (k = 3; k <= 130; k++) {
        for (j = 3; j <= 130; j++) {
          index =  k*(134)*(134) +
            j*(134) + 3;
          for (i = 3; i <= 130; i++, index++) {
            dens_field->Array[index] = (float)(i + 1000*j + 1000000*k);
          }
        }
      }
    
      delete dens_field;
    
      return lineargrid;
    }

Notice how this function uses ``CreateFieldArrayFloat`` to set the
values of the density array.

Now, here's a program that creates a uniform grid, and looks at
some of the attributes:

.. code-block:: c

    Eint32 main(Eint32 argc, char *argv[]) {
    
      CommunicationInitialize(&argc, &argv);
    
      grid *agrid = Linear3DGrid();
    
      EnzoArrayFloat *dens = agrid->CreateFieldArrayFloat(Density);
    
      Eint32 index = 7 + 8*134 + 9*134*134;
    
      printf("density rank = %"ISYM"\n", dens->Rank);
      printf("density dim[0]  = %"ISYM"\n", dens->Dimension[0]);
      printf("density start[0]  = %"ISYM"\n", dens->StartIndex[0]);
      printf("density end[0]  = %"ISYM"\n", dens->EndIndex[0], 130);
      printf("density field[7 + 8*134 + 9*134*134] = %"FSYM"\n", dens->Array[index]);
    
      delete dens;
      delete agrid;
    
      // End the overall test suite                                                                                                  
      CommunicationFinalize();
    
      return 0;
    }

This is a complete program,
``field_array_example.C``;
what this snippet lacks is the fairly
long list of header files that need to be included. You can compile
this by calling ``make field_array_example.exe`` in source directory.

.. _Field_List_Reference:

Field List Reference
--------------------

The following table is a partial list of the fields in Enzo.  The **Field Type ID** is defined in the ``typedef.h`` file.

=============  ======================  ======================  ==========  ===============
Field Type ID  Field Number            Field Name              Data Type   Array or Vector
=============  ======================  ======================  ==========  ===============
0              Density                 "Density"               float       Array
1              TotalEnergy             "TotalEnergy"           float       Array
2              InternalEnergy          "InternalEnergy"        float       Array
3              Pressure                "Pressure"              float       Array
4              Velocity1               "Velocity1"             float       Array
5              Velocity2               "Velocity2"             float       Array
6              Velocity3               "Velocity3"             float       Array
7              ElectronDensity         "ElectronDensity"       float       Array
8              HIDensity               "HIDensity"             float       Array
9              HIIDensity              "HIIDensity"            float       Array
10             HeIDensity              "HeIDensity"            float       Array
11             HeIIDensity             "HeIIDensity"           float       Array
12             HeIIIDensity            "HeIIIDensity"          float       Array
13             HMDensity               "HMDensity"             float       Array
14             H2IDensity              "H2IDensity"            float       Array
15             H2IIDensity             "H2IIDensity"           float       Array
16             DIDensity               "DIDensity"             float       Array
17             DIIDensity              "DIIDensity"            float       Array
18             HDIDensity              "HDIDensity"            float       Array
19             SNColour
20             Metallicity             "Metallicity"           float       Array
21             ExtraType0              "ExtraType0"            float       Array
22             ExtraType1              "ExtraType1"            float       Array
30             GravPotential           "GravPotential"         float       Array
31             Acceleration0           "Acceleration0"         float       Array
32             Acceleration1           "Acceleration1"         float       Array
33             Acceleration2           "Acceleration2"         float       Array
37             gParticlePosition       "ParticlePosition"      FLOAT       Vector
38             gParticleVelocity       "ParticleVelocity"      float       Vector
39             gParticleMass           "ParticleMass"          float       Array
40             gParticleAcceleration   "ParticleAcceleration"  float       Vector
41             gParticleNumber         "ParticleNumber"        int         Array
42             gParticleType           "ParticleType"          int         Array
43             gParticleAttribute      "ParticleAttribute"     float       Vector
44             gPotentialField         "PotentialField"        float       Array
45             gAccelerationField      "AccelerationField"     float       Vector
46             gGravitatingMassField   "GravitatingMassField"  float       Array
47             gFlaggingField          "FlaggingField"         int         Array
48             gVelocity               "Velocity"              float       Vector
=============  ======================  ======================  ==========  ===============

