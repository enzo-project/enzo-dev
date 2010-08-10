.. _AddingNewParameters:

Adding a new parameter to Enzo
==============================

If your parameter is only used for a problem initialization, this
page is not relevant for you. You should just read it in during
``ProblemInitialize.C`` where ``Problem`` is replaced by the name of
the problem type.

If you're extending Enzo for any reason, you'll probably need to
add a new switch or parameter to the code. Currently, this page
describes the simplest, most brute force method. There are four
files you'll need to edit to make this happen.


-  ``global_data.h`` holds all the global data. It's included in
   almost all Enzo files. Your parameter should be added like this:
   ::

       EXTERN int MyInt;
       EXTERN float MyFloat;

   EXTERN is a macro that either maps to extern if ``USE\_STORAGE`` is
   defined, or nothing if ``USE\_STORAGE`` is not defined. ``USE\_STORAGE`` is
   defined in ``enzo.C`` before the inclusion of ``global\_data.h``, and
   undefined after.


-  ``SetDefaultGlobalValues.C`` sets the default global values. Set
   your value here.


-  ``ReadParameterFile.C`` reads the parameter file. In this routine,
   each line is read from the file and is compared to the given
   parameters with ``sscanf()``. Your line should look like this:
   ::

        ret += sscanf(line, "MyFloat      = %"FSYM, &MyFloat);
        ret += sscanf(line, "MyInt        = %"ISYM, &MyInteger);

   and should be inserted somewhere in the loop where line is
   relevant. Note that ``ISYM`` and ``FSYM`` are the generalized integer and
   float I/O macro, which exist to take care of the dynamic hijacking
   of 'float'.
   `See this page for more information? </wiki/Tutorials/FloatIsDouble>`_.
   The ``ret +=`` controls whether the line has been read, or if Enzo
   should issue a warning about the line. Note also that ``sscanf()`` is
   neutral to the amount of consecutive whitespace in the format
   string argument.


-  ``WriteParameterFile.C`` writes the restart parameter file.
   Somewhere before the end of the routine, you should add something
   that looks like
   ::

         fprintf(fptr, "MyFloat       = %"GSYM"\n", MyFloat);
         fprintf(fptr, "MyInt         = %"ISYM"\n", MyInt);

   Note the use of quotes here and in the previous code snippet. This
   is correct.

For testing purposes you can verify that your new parameter is being correctly read in by 
adding a line like this at the bottom of ``ReadParameterFile.C``:

::

      fprintf(stdout, "MyFloat %f MyInt %d\n", MyFloat, MyInt);
      return SUCCESS;
    }

