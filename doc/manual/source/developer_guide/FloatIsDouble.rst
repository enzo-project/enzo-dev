.. _FloatIsDouble:
.. _VariablePrecisionInEnzo:

Variable precision in Enzo
==========================

In order to provide some global control over variable precision,
Enzo uses a set of macros that control how the code treats integer
and floating-point precision by overriding the float and int data
types, and by introducing the ``FLOAT`` macro. This is a major sticking
point for new users to the code, and this page is an attempt to
clarify the issue as much as possible.

Floating-point precision
------------------------

There are two different kinds of floating-point quantities in Enzo, those that
explicitly deal with positional information (grid edges/locations, cell sizes,
particle positions, and so on), and those that deal with non-position
information (baryon density, temperature, velocity, etc.) Any variables that
deal with position information should be declared as the ``FLOAT`` data type. For
example:

.. code-block:: c

    FLOAT xpos, ypos, zpos;

A quantity that deals with non-positional information would be
declared using the ``float`` data type:

.. code-block:: c

    float cell_HI_density, cell_H2I_density, cell_temperature;

The actual precision of ``float`` and ``FLOAT`` are controlled by the
Makefile system (see :ref:`obtaining_and_building_enzo`.) To set the
non-positional precision to 64-bit (double), you would issue this
command:

::

    make precision-64

before compiling the code. Similarly, to set the positional
precision to 64-bit (double), you would issue this command:

::

    make particles-64

The allowable values for non-positional precision are 32 and 64
bits, and for positional precision are 32, 64, and 128 bits. It is
not recommended that you use particles-128 unless you need more
than 30 or so levels of AMR, since ``long double`` arithmetic generally
requires software libraries and can be very slow. Also note that
the 128-bit precision code is not terribly stable, and only works
on some systems (and with some sets of compilers). Use this with
**extreme caution**.

**Mixing ``float`` and ``FLOAT``**: One can mix the ``float`` and ``FLOAT`` data
types, but some care is required since the two are not necessarily
the same precision. Compilers will generally promote the variables
to the higher precision of the two, but this is not always true.
The Enzo developers have chosen to make the assumption that the
precision of ``FLOAT`` is always the same as, or greater than, the
precision of ``float``. So, when precision is critical or when mixing
``float`` and ``FLOAT``, we recommend that you always promote all variables
to ``FLOAT``. Regardless, it is a good idea to check that your code is
producing sensible results.

Integer precision
-----------------

There is only one commonly-used type of integer in Enzo, which is
``int``. This is controlled by the the integers- makefile command. For
example,

::

    make integers-64

would force all ints to be 64-bit integers (``long int``). The
allowable integer values are 32 and 64 bit. In general, the only
time one would need 64-bit ints is if you are using more than
2\ :sup:`31`\  particles, since signed integers are used for the
particle index numbers, and chaos will ensue if you have duplicate
(or, worse, negative) particle indices.

Precision macros and printf/scanf
---------------------------------

In order to keep the printf family of commands happy, Enzo uses
several macros. ``ISYM`` is used for integers, ``FSYM`` and ``ESYM`` for ``float``, and
``PSYM`` and ``GSYM`` for ``FLOAT`` (the latter of each pair outputs floats in
exponential notation). Additionally, when writing ``FLOAT`` data to a
file that will be read back in by Enzo (such as to the parameter or
hierarchy file), it is wise to use ``GOUTSYM``. In a printf or scanf
statement, this macro will be replaced with the actual string
literal statement.

An example of this usage macro in a printf statement to write out a
float is:

.. code-block:: c

    printf("Hello there, your float value is %"FSYM".\n", some_float);

and to read in a set of three position coordinates using scanf out
of a string named line:

.. code-block:: c

    sscanf(line,"PartPos  = %"PSYM" %"PSYM" %"PSYM, &XPOS, &YPOS, &ZPOS);

Note the somewhat counterintuitive use of quotation marks after the
3rd ``PSYM``. For a large number of examples of how to use these
macros, please refer to the files ``ReadParameterFile.C`` and
``WriteParameterFile.C`` in the Enzo source code.

The Fortran-C++ interface
-------------------------

It is critical to make sure that if you are interfacing Fortran
and C/C++ code, the variable precision agrees between the two
languages. Compilers do not attempt to ensure that calls from C/C++
to Fortran make any sense, so the user is manifestly on their own.
To this end, when writing Fortran code you must ensure that your
variables are declared with the correct type.  Unlike Enzo's C/C++
routines that overwrite the default ``float`` and ``int``
types with their single/double precision equivalents, Enzo's Fortran
routines do not overwrite the basic data types.  Hence, we have
created unique type identifiers for the Fortran routines that map to
Enzo's ``float``, ``FLOAT`` and ``int`` types, as specified below:

==================  ==============
**Enzo C/C++**      **Enzo F/F90**
``float``           ``R_PREC``
``int``             ``INTG_PREC``
``FLOAT``           ``P_PREC``
==================  ==============

In addition, Fortran allows additional data types for both ``logical``
and ``complex`` variables.  In Enzo, the precision of these variables
may be chosen to match Enzo's ``int`` and ``float`` values from C/C++
using the F/F90 types ``LOGIC_PREC`` and ``CMPLX_PREC`` respectively.

Moreover, unlike C/C++, hard-coded constants in Fortran routines
default to single-precision values.  This can be especially
troublesome when calling a Fortran subroutine or function with
constants as their inputs, or when writing complicated formulas using
constants that must be of higher precision.  To this end, we have
defined four type-modifier Fortran suffixes, that can be used to
declare constants of differing precision:

===================  ==========
**Variable Type**    **Suffix**
``R_PREC``           ``RKIND``
``INTG_PREC``        ``IKIND``
``P_PREC``           ``PKIND``
``LOGIC_PREC``       ``LKIND``
===================  ==========

Note: since a complex number in Fortran is defined through a pair of
real numbers, to create a complex constant of type ``CMPLX_PREC`` you
would use the ``RKIND`` suffix on both components.

For example, the type specifiers and constant suffixes could be used
in the following ways: 

.. code-block:: c

    c     Declarations
          R_PREC     third, tenth
          INTG_PREC  one
          P_PREC     fifth
          CMPLX_PREC two_i
          LOGIC_PREC test

    c     Calculations
          third = 1._RKIND / 3._RKIND
	  tenth = 1.e-1_RKIND
          one   = 1_IKIND
          fifth = real(1, PKIND) / 5._PKIND
          two_i = (0._RKIND, 2._RKIND)
          test  = .true._LKIND


All of these type definitions are supplied in the file
``fortran_types.def`` and should be included within a Fortran routine
within the scope of the function, after any ``implicit none``
declaration, and before declaring any variables, e.g.

.. code-block:: c

          subroutine foo(a)
             implicit none
    #include "fortran_types.def"
             R_PREC a

The Enzo build system will preprocess this file to include
``fortran_types.def`` at the specified location in the file, prior to
compilation.  Moreover, the spacing in this file is usable using
either fixed-source-form or free-source-form Fortran files.

**A word of warning:** mismatching the data types between C/C++ and
Fortran codes can cause misalignment in the data, and will often
result in nonsense values that will break Enzo elsewhere in the
code. This can be particularly tricky to debug if the values are
not used immediately after they are modified!

If you need more details…
-------------------------

If you need more detailed information on this particular subject,
there is no substitute for looking at the source code. All of these
macros are defined in the Enzo source code file
``macros_and_parameters.h``. Just look for this comment:

.. code-block:: c

    /* Precision-dependent definitions */

There are many examples of using the IO macros in
``ReadParameterFile.C`` and ``WriteParameterFile.C``.

Also, please note that this set of macros may be replaced with a
more robust set of macros in future versions.
