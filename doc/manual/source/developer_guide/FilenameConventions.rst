File naming conventions
-----------------------

The large number of source files can be intimidating even to the
experienced Enzo developer, and this page describes some of the naming
conventions.  Familiarity with ``grep`` or ``ack`` and pipes like ``ls
-1 | grep`` are essential.  Also routines with a similar functionality
are grouped together with a common name.

Here are some file naming rules that are used.

1. Internal capitalization is used for C files, all lowercase with
   underscores for fortran files and header files. All Fortran files
   end with .F.

2. With very few exceptions, Enzo has a one function per file layout, with the
   file name being the function name. 

3. Object methods have the object name prepended to the beginning,
   such as the member of the grid class ``SolveHydroEquations`` lives
   in the file ``Grid_SolveHydroEquations.C``.


