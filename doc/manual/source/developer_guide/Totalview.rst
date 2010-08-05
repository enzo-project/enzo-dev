Debugging Enzo .src files with Totalview
========================================

The totalview debugger is a widely popular parallel debugging tool,
and makes things like finding segfaults extremely easy.

This is not intended to be a tutorial on totalview-- your super
computer center and the internet will provide you with better
tutorials that I will. This is just a few caveats when using enzo.


#. Compile with -g and -O1 or-O1. Missing the -g will cause
   totalview to not be able to find your code. Higher optimizations do
   really weird things, like displaying nonsensical values and
   executing things out of order. (Unsurprising behavior, if you've
   optimized the code.)
#. Totalview doesn't know anything about the preprocessor that the
   fortran .src (and fortran 90 .src90) files go through. So while you
   may see the C/C++ routines, it won't initially see the fortran
   source because it's looking for .f or .f90 files. You can get
   around this by building the actual .f or .f90 files with the
   following script, before launching totalview:
   ::

       foreach i (`ls *.src`)
               make `basename $i src`f
       end
       
       foreach i (`ls *.src90`)
               make `basename $i src90`f90
       end

   Note that you'll now be looking at the pre-processed version-- line
   numbers will not match the X.src file you're actually interested
   in. Don't try to edit these .f files, because they will be deleted
   when the code next compiles. This is for debugging purposes ONLY.

David Collins, Nov. 2008.


