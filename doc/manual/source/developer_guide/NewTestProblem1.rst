.. _AddingANewTestProblem:

Adding a new Test Problem.
==========================

This is the best place to start your Enzo Development Career. Even
if you're not interested in actually writing a new problem
generator, in this page I'll discuss the basic Enzo datastructures
and programming patterns.

One deficiency in this tutorial is the lack of Particles. This is
not an oversight, but due to the fact that the author of the
article doesn't really use particles, as he's not a cosmologist.
These will be added in the future, but particles are really not
that big of a deal when it comes to the general Enzo data
structures. All the information herein is still essential.

Overview
--------

Essentially, you need two write files: ``MyProblemInitialize.C`` and
``Grid_MyProblemInitializeGrid.C``. We'll be discussing these two
files. ``MyProblemInitialize`` is the basic setup code that sets up
parameters and the hierarchy, and ``MyProblemInitializeGrid`` is a
member function of the grid class, and actually allocates and
assigns data. There are several pitfalls to setting up these files,
so read these pages carefully.

We strongly recommend reading everything that proceeds this page on
the :doc:`../tutorials/index` page and the page about version control
and regression testing, :doc:`ModificationIntro`.

Lastly, please give your problem a reasonable name. I'll be using
``MyProblem`` throughout this tutorial. Please change this to something
that reflects the problem you're installing.

Setup and Installation
----------------------

Please follow the general Enzo naming convention and call your
routines ``MyProblemInitialize`` and store it in ``MyProblemInitialize.C``,
and ``MyProblemInitializeGrid`` and store it in
``Grid_MyProblemInitializeGrid.C``

You'll need to install your code in three places.

.. highlight:: none

#. ``Make.config.objects``
    is the file that lists all the source file objects needed to build Enzo. Put

    ::

       MyProblemInitialize.o\ 
       Grid_MyProblemInitializeGrid.o\

    somewhere in the list of objects. If you want to make things really
    clean, you can add your own variable to the Makefile and have it
    driven by a command line switch, but this isn't necessary.

#. ``Grid.h``. You'll need to put
    ``MyProblemInitializeGrid`` in this the grid class definition. Put it
    with the rest of the ``*InitializeGrid`` routines.

#. ``InitializeNew.C``. Put
    ``MyProblemInitialize`` in ``InitializeNew``. At the end of the large block
    of ``*Initialize``, take the next unused ProblemType number and
    install your code. It should look something like this:

    .. code-block:: c

          // 61) Protostellar Collapse                                                                                 
          if (ProblemType == 61)
            ret = ProtostellarCollapseInitialize(fptr, Outfptr, TopGrid, MetaData);
        
          // 62) My New Problem 
          if ( ProblemType == 62 )
            ret = MyProblemInitialize(fptr, Outfptr, TopGrid, MetaData);
        
          // Insert new problem intializer here...                                                                     
        
          if (ret == INT_UNDEFINED) {
            fprintf(stderr, "Problem Type %"ISYM" undefined.\n", ProblemType);
            return FAIL;
          }

To call your problem generator, make sure ProblemType = 62 is in
your parameter file. (Or, if 62 is taken, whatever the next unused
value is.)

The return value ``ret`` is used to check for errors and invalid values
of ``ProblemType``. The function signature will be discussed in the
next section.

Also, don't forget to put the proto type at the top:

.. code-block:: c

    int MyProblemInitialize(FILE *fptr, FILE *Outfptr,
                                       HierarchyEntry &TopGrid,
                                       TopGridData &MetaData);

We will revisit ``InitializeNew`` at the end. For almost all problems,
this will be all you do for these three files.

MyProblemInitialize
-------------------

The primary drive routine is called ``MyProblemInitialize``. It
basically sets up some global values, problem specific values, and
the hierarchy before calling ``MyProblemInitializeGrid``.

Function Signature
~~~~~~~~~~~~~~~~~~

The function signature of ``MyProblemInitialize`` is fairly rigid. It
should look exactly like the prototype you installed in
``InitializeNew``. There are 4 arguments that you'll almost certainly
need, and one additional argument that only rare problems will
need. You won't likely have any need to add any other arguments. In
order, they are:

#. ``FILE *fptr`` This is the pointer to the parameter file argument to
    Enzo. It's opened and closed in InitializeNew You can read
    parameters if you like, see below.

#. ``FILE *Outfptr`` This is the output pointer, a file called "amr.out."
    This file contains the derived details of your problem setup for
    your record. There is no necessary output for this, it's for the
    users convenience.

#. ``HierarchyEntry &TopGrid`` This is the pointer to the top of the
    Hierarchy Linked List. For details of the linked list,
    :doc:`../reference/LinkedLists`. For most problem types, it
    points to the undivided root grid, which is a grid the full size
    of the top grid, where you will be initializing your data. For
    problems that are too large for the entire root grid to be
    allocated, we use the ParallelRootGridIO functionality, to be
    discussed later. (Please read everything between here and there.)

#. ``TopGridData &MetaData`` This is the structure that contains the meta
    data describing the Top Grid. Things like boundary condition,
    problem domain size, rank, and dimension are stored here.
    See ``TopGridData.h`` for a complete list of the contents.

If you want to write a problem with Dirichlet boundary conditions,
for instance jet inflow, you will need to add a fifth argument to
the function (and, of course, it's called in ``InitializeNew``). This is
the external boundary, ``ExternalBoundary &Exterior``. This is the
External Boundary object, which you will need to deal with.  We will
not be discussing this here. If you need to be
doing a problem with boundary conditions other than the big 3
(periodic, reflecting, outflow) then we recommend you read the
entirety of this tutorial, then follow what's done with the
DoubleMach problem, which is problem type 4. You will also need to
examine ``Grid_SetExternalBoundaryValues.C``

Necessary Headers
~~~~~~~~~~~~~~~~~

The essential header files for ``MyProblemInitialize`` are the
following:

.. code-block:: c

    #include <stdio.h>
    #include <string.h>
    #include "macros_and_parameters.h"
    #include "typedefs.h"
    #include "global_data.h"
    #include "Fluxes.h"
    #include "GridList.h"
    #include "ExternalBoundary.h"
    #include "Grid.h"
    #include "Hierarchy.h"
    #include "TopGridData.h"

These should be in this order, to ensure proper definitions across
different header files. You should be familiar with the two
standard headers <stdio.h> and <string.h>

In brief, these are:

- ``macros_and_parameters.h`` The standard set of macros. This takes
    care of the float promotion so its inclusion is
    **ABSOLUTELY ESSENTIAL**

- ``typedefs.h`` This takes
    care of enumerates for parameters like the hydro method.

- ``global_data.h`` There
    are a lot of global parameters in Enzo. This houses them.

- ``Fluxes.h`` Definition of the
    flux object. Not necessary for your objects, but I think its
    necessary for the later

- ``GridList.h`` I don't think
    this is necessary, but it's usually included.

- ``ExternalBoundary.h`` This defines the external boundary object. Even
    if you're not including the external boundary, it's
    necessary for the following headers.

- ``Grid.h`` This defines the grid
    class, which you'll definitely need.

- ``Hierarchy.h`` This defines the Hierarchy Entry linked list.

- ``TopGridData.h`` This defines the meta data object.

More information can be found in :doc:`../reference/Headers`.

Necessary Assignments
~~~~~~~~~~~~~~~~~~~~~

There are two arrays that need to be filled in ``MyProblemInitialize``.
One of them is **ABSOLUTELY ESSENTIAL** for the functioning of the
code. These are ``DataLabel`` and ``DataUnits``. Both of these are arrays
of strings that will be used to label the HDF5 output files. Each
element of the array corresponds to an element of the BaryonField
array (more on this later) and MUST be defined in the same order.
There is not a mechanism to ensure that you do this right, so don't
screw it up.

DataLabel
^^^^^^^^^

This is the actual name of the field in the HDF5 file. Messing this
up is asking for trouble. If you're not using chemistry, you'll
want something that looks like this. If you change the actual
names, you guarantee that an analysis tool somewhere will break, so
don't do it. See
``CosmologySimulationInitialize.C`` for
a more complete list, including extra chemical species.

.. code-block:: c

      char *DensName = "Density";
      char *TEName   = "TotalEnergy";
      char *GEName   = "GasEnergy";
      char *Vel1Name = "x-velocity";
      char *Vel2Name = "y-velocity";
      char *Vel3Name = "z-velocity";
      i = 0;
      DataLabel[i++] = DensName;
      DataLabel[i++] = TEName;
      if (DualEnergyFormalism)
        DataLabel[i++] = GEName;
      DataLabel[i++] = Vel1Name;
      DataLabel[i++] = Vel2Name;
      DataLabel[i++] = Vel3Name;

DataUnits
^^^^^^^^^

The units really don't matter very much. They're usually set to
NULL

Reading from the Parameter File
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You may want to read in problem specific parameters. PLEASE do not
put problem specific parameters in the main parameter file reader.

The usual pattern reads each line of the parameter file, and tries
to match each line with a parameter. This allows the parameter file
to be independent of of order. The typical pattern looks like
this:

.. code-block:: c

      float MyVelocity, MyDensity;
      char line[MAX_LINE_LENGTH];
      while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
       ret = 0;
    
        /* read parameters */
    
        ret += sscanf(line, "MyProblemVelocity      = %"FSYM,
                      &MyVelocity);
        ret += sscanf(line, "MyProblemDensity      = %"FSYM,
                      &MyDensity);
        if (ret == 0 && strstr(line, "=") && strstr(line, "MyProblem") &&
            line[0] != '#' && MyProcessorNumber == ROOT_PROCESSOR)
          fprintf(stderr,
             "warning: the following parameter line was not interpreted:\n%s\n",
                  line);
      }

If you're not familiar with these functions,
`here is a good list of standard C functions <http://www.cppreference.com/all_c_functions.html>`_.

The last line checks for errors in parameters that start with
``MyProblem``. Everything involving this routine should be prepended
with ``MyProblem``. In the file ``ReadParameterFile.C``, the parameter file
is read and any lines not recognized are thrown as errors; this is
the section identified with

.. code-block:: c

    \* check to see if the line belongs to one of the test problems \*/.
    
You must add your prefix (in this
case, ``MyProblem``) to the list of test problem prefixes considered in
this section:

.. code-block:: c

        if (strstr(line, "MyProblem")           ) ret++;

or else it will register as an error.

.. _UnigridInitialize:

Calling the Grid Initializer: Unigrid
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For a small, unigrid problem, the problem initializer is called
using the standard Enzo function call procedure.

.. code-block:: c

    if( TopGrid.GridData->MyProblemInitializeGrid(MyVelocity, MyDensity) == FAIL ){
      fprintf(stderr,"MyProblemInitialize: Error in MyProblemInitializeGrid\n");
      return FAIL;

``TopGrid`` is the ``HierarchyEntry`` that starts the hierarchy linked
list. It's member ``GridData`` is a pointer to the actual grid object
that you will be modifying.

We will be discussing AMR problems, and large problems that require
parallel startup later.

.. _InitializeGrid:

MyProblemInitializeGrid
-----------------------

``MyProblemInitializeGrid`` is the member function of the grid class.
As a member function, it can access the private data, most
importantly ``BaryonField``. ``BaryonField`` is an array of pointers that
stores the actual data that the simulator is interested in.

.. code-block:: c

    float *BaryonField[MAX_NUMBER_OF_BARYON_FIELDS];

Necessary Actions
~~~~~~~~~~~~~~~~~

There are four things that this routine ABSOLUTELY MUST do, and
they MUST BE DONE IN THIS ORDER.

#. Set up the FieldType array and define NumberOfBaryonFields.

    The FieldType array is an array of type field_type, a type defined
    in ``src/enzo/typedefs.h``. It is used to relate physics to the
    actual BaryonField element.

    ``NumberOfBaryonFields`` is the number of valid, allocated fields. This
    can be as little as 5 for pure fluid dynamics, or as many as you
    have chemistry to deal with.

    A typical pattern looks like this:

    .. code-block:: c

      NumberOfBaryonFields = 0;
      FieldType[NumberOfBaryonFields++] = Density;
      FieldType[NumberOfBaryonFields++] = TotalEnergy;
      if (DualEnergyFormalism)
        FieldType[NumberOfBaryonFields++] = InternalEnergy;
      FieldType[NumberOfBaryonFields++] = Velocity1;
      vel = NumberOfBaryonFields - 1;
      if (GridRank > 1)
        FieldType[NumberOfBaryonFields++] = Velocity2;
      if (GridRank > 2)
        FieldType[NumberOfBaryonFields++] = Velocity3;

    All the right hand side of those assigns can be found in
    ``typedefs.h``.

    Note that all processors must have this information defined for all
    grids, so this MUST come before step 2.

#. Exit for remote grids.

    Generally, grid member functions have two modes: things that all
    processors can do, and things that only processors that own the
    data can do. Usually, the routine simply exits if the processor
    doesn't own the data:

    .. code-block:: c

      if (ProcessorNumber != MyProcessorNumber)
        return SUCCESS;

    ``ProcessorNumber`` is a grid member that stores which processor
    actually has the data, and ``MyProcessorNumber`` is the global number
    of the processor.
    
    Processors that don't get this data need to not execute the rest of the code.

#. Allocate the BaryonFields

    .. code-block:: c
    
        this->AllocateGrids()
    
    does the trick. More details on this in the Parallel section below.

#. Assign values to the ``BaryonField``. See the page on Baryon Field Access for details.
    Note that for problems that are perturbations on a homogenous background,
    the routine ``InitializeUniformGrid`` has been provided. See that routine for its function signature.

Initializing AMR problems
~~~~~~~~~~~~~~~~~~~~~~~~~

For problems that you want to initialize in an AMR fashion, all the previous
steps apply. However, instead of simply calling the problem initializer on the
Top Grid, one must now initialize a ``HierarchyEntry`` linked list (of which ``TopGrid``
is the head) and call the problem initializer on each subgrid. There are several
ways to do this, depending on the complexity of the code. One first needs to
understand the ``HierarchyEntry`` linked list. This Page gives a tutorial on the
linked lists, and links to examples in the code.

Using ParallelRootGridIO
~~~~~~~~~~~~~~~~~~~~~~~~

Main article: :doc:`NewTestProblem3`

``ParallelRootGridIO`` is a fairly complex piece of code. If you absolutely
must do this in the code, it is recommended that you read the description
of the inner workings of ``ParallelRootGridIO`` and then cloning what's done
for the ``CosmologyInitialize`` routines.
