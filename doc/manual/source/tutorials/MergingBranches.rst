Merging Branches
================

This is an example of a development process using Subversion to
help branch and merge from the trunk and back into it. This
particular example happens to be the result of my
(` Rick's <http://lca.ucsd.edu/projects/rpwagner>`_) comment about
`GridAccessors? </wiki/GridAccessors>`_.

Phase I: Branching the Trunk
----------------------------

Create a New Branch
~~~~~~~~~~~~~~~~~~~

I didn't want to make this change in [browser:branches/rpwagner my
branch], since I have a lot of crap in there. Instead, I wanted to
make a change to the trunk, then merge that down to my branch. This
example is the same as the one in
`EnzoBranches? </wiki/EnzoBranches>`_.

::

    cable:~ rpwagner$ svn mkdir -m "creating a temporary working branch" \
              http://lca.ucsd.edu/svn/Enzo/branches/fieldarray
    cable:~ rpwagner$ svn copy -m "copying the trunk to the new branch" \
              http://lca.ucsd.edu/svn/Enzo/trunk/devel/Enzo http://lca.ucsd.edu/svn/Enzo/branches/fieldarray

Check Out the Branch
~~~~~~~~~~~~~~~~~~~~

Next, I check out the branch to get to work. This step creates a
new local directory NewFields where I will hack on this copy of the
trunk.

::

    cable:~ rpwagner$ svn co http://lca.ucsd.edu/svn/Enzo/branches/fieldarray/Enzo NewFields

Time Passes
~~~~~~~~~~~

About 48 hours go by, and I hack a bunch of stuff. In particular, I
added a [browser:trunk/devel/Enzo/amr\_mpi/src/EnzoArray.h new
array class] that can be used to
[browser:trunk/devel/Enzo/amr\_mpi/src/Grid\_CreateFieldArray.C
access grid data]. There are examples of using this in the
[browser:trunk/devel/Enzo/amr\_mpi/src/gridtests.C test code].

When I'm done, I commit my changes to the branch. Remember: at this
point, I am only committing changes into the branch I just made.
Moving it to the trunk comes next.

::

    cable:~ rpwagner$ svn commit -m "added new field arrays"

Phase II: Testing
-----------------

Before I merge anything back to the trunk, I go and find
` James <http://jbpc.ucsd.edu/~jbordner/>`_ (he's the one in the
middle). James glances at the code and runs
` lcatest <http://lca.ucsd.eud/projects/lcatest>`_ to make sure
that I haven't broken anything.

The punch line to this phase is that I'm checking with the trunk
owner before I try to move my stuff into the trunk.

Phase III: Merging Back
-----------------------

Check Out the Trunk
~~~~~~~~~~~~~~~~~~~

We're going to need a local copy of the trunk to merge into. This
is a much better idea than just merging to URL's.

::

    cable:~/tmp/merge rpwagner$ svn co http://lca.ucsd.edu/svn/Enzo/trunk/devel/Enzo  

Merge From the Branch
~~~~~~~~~~~~~~~~~~~~~

Now we merge our changes into the local working copy. The revision
range I'm using is based on the point where I copied the trunk
over, since this is the last time the two branches were synced.

::

    cable:~/tmp/merge/Enzo rpwagner$ svn merge -r1157:HEAD http://lca.ucsd.edu/svn/Enzo/branches/fieldarray/Enzo .
    A    amr_mpi/src/Grid_CreateFieldArray.C
    U    amr_mpi/src/Make.config.objects
    U    amr_mpi/src/Make.config.assemble
    A    amr_mpi/src/enzo_unit_tests.h
    U    amr_mpi/src/Grid.h
    U    amr_mpi/src/Make.mach.rpwagner-cable
    A    amr_mpi/src/EnzoArray.h
    U    amr_mpi/src/typedefs.h
    U    amr_mpi/src/Makefile.config
    A    amr_mpi/src/gridtests.C
    cable:~/tmp/merge/Enzo rpwagner$ 

This applies the changes from revision 1157 to the latest, that
have been done to the fieldarry branch, to the local copy of the
trunk. From this, you can see the list of files that were updated,
added, or deleted. If something was really bad, there could even be
conflicts.

Check Updated Files
~~~~~~~~~~~~~~~~~~~

If you wonder what changes are going to be applied, you should use
svn diff and svn status. svn status will show you which files in
the local copy have been added or modified after the merge; these
are the ones that will get committed back to the server.

::

    cable:~/tmp/merge/Enzo/amr_mpi/src rpwagner$ svn status
    A  +   Grid_CreateFieldArray.C
    M      Make.config.objects
    M      Make.config.assemble
    A  +   enzo_unit_tests.h
    M      Grid.h
    M      Make.mach.rpwagner-cable
    A  +   EnzoArray.h
    M      typedefs.h
    M      Makefile.config
    A  +   gridtests.C

svn diff can give you the details on a particular file, so you can
double check that there won't be any surprises. We'll look at
Grid.h, since that has the potential to cause the most havoc.

::

    cable:~/tmp/merge/Enzo/amr_mpi/src rpwagner$ svn diff Grid.h
    Index: Grid.h
    ===================================================================
    --- Grid.h      (revision 1174)
    +++ Grid.h      (working copy)
    @@ -13,7 +13,6 @@
     #ifndef GRID_DEFINED__
     #define GRID_DEFINED__
     
    -
     #include "ProtoSubgrid.h"
     #include "ListOfParticles.h"
     #include "region.h"
    @@ -35,6 +34,8 @@
     
     struct LevelHierarchyEntry;
     
    +#include "EnzoArray.h"
    +
     class grid
     {
      private:
    @@ -138,8 +139,20 @@
     
     /* Read grid data from a file (returns: success/failure) */
     
    -   int ReadGrid(FILE *main_file_pointer, int GridID);
    +  int ReadGrid(FILE *main_file_pointer, int GridID);
     
    +/* Get field or particle data based on name or integer 
    +   defined in typedefs.h. Details are in Grid_CreateFieldArray.C. */
    +
    +  EnzoArrayInt *CreateFieldArrayInt(field_type field);
    +  EnzoArrayInt *CreateFieldArrayInt(char *field_name);
    +  
    +  EnzoArrayFloat *CreateFieldArrayFloat(field_type field);
    +  EnzoArrayFloat *CreateFieldArrayFloat(char *field_name);
    +  
    +  EnzoArrayFLOAT *CreateFieldArrayFLOAT(field_type field);
    +  EnzoArrayFLOAT *CreateFieldArrayFLOAT(char *field_name);
    +
     /* Write unigrid cubes to a file (returns: success/failure) */
     
        int WriteCube(char *base_name, int grid_id, int TGdims[]);

OK, I added a header file and new functions, no surprises.

Build and Test
~~~~~~~~~~~~~~

Nothing new here, just compile the code and make sure it works. For
the trunk, it should pass
` lcatest <http://lca.ucsd.eud/projects/lcatest>`_, and whatever
other tests are defined. I happened to have used this as a chance
to get unit tests building in the trunk, so that's the test I'm
going show.

::

    cable:~/tmp/merge/Enzo/amr_mpi/src rpwagner$ make gridtests.exe
    Compiling gridtests.C
    Linking
    Success!
    cable:~/tmp/merge/Enzo/amr_mpi/src rpwagner$ ./gridtests.exe 
    > get density float array...
    > get density float array by name...
    > get particle type array...
    > get particle type array by name...
    > check for NULL - float, int...
    > check for NULL - int, float...
    > get velocity array...
    > get velocity array by name...
    
    --- Results ---
    Tests run:    8
    Passes:      48
    Failures:     0

This was done after the main executable was built, so the needed
object files were already there.

Commit
~~~~~~

The merge is basically done at this point, we just need to commit
this local copy to have it show up in the trunk.

::

    cable:~/tmp/merge/Enzo/amr_mpi/src rpwagner$ svn commit -m "committing new array fields to the trunk"
    Adding         src/EnzoArray.h
    Sending        src/Grid.h
    Adding         src/Grid_CreateFieldArray.C
    Sending        src/Make.config.assemble
    Sending        src/Make.config.objects
    Sending        src/Make.mach.rpwagner-cable
    Sending        src/Makefile.config
    Adding         src/enzo_unit_tests.h
    Adding         src/gridtests.C
    Sending        src/typedefs.h
    Transmitting file data ......
    Committed revision 1175


