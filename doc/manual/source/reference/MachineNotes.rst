.. _MachineNotes:

Machine Specific Notes
======================

Here we will mention some miscellaneous notes on specific machines.
This is merely a list of pitfalls or things we have found useful,
and by no means a replacement to the documentation.

NICS: Kraken
------------

`http://www.nics.tennessee.edu/computing-resources/kraken <http://www.nics.tennessee.edu/computing-resources/kraken>`_

Important
~~~~~~~~~

Serious errors have been found with a few Enzo routines when using
-O2 and the PGI compilers on Kraken. Use with caution.

Trace Trap Flags
~~~~~~~~~~~~~~~~

Useful for debugging, but slows the code down. You can find this
info in the pgCC man page. (Not all compilers have decent trace
trapping, so it deserves a mention here.)

::

     -Ktrap=[option,[option]...]
            Controls the behavior of the processor when
            exceptions occur.  Possible options include
            -Ktrap=divz  Trap on divide by zero.
            -Ktrap=fp  Trap on floating point exceptions.          
            -Ktrap=align Trap on memory alignment errors, currently ignored
            -Ktrap=denorm Trap on denormalized operands.
            -Ktrap=inexact Trap on inexact result.
            -Ktrap=inv Trap on invalid operands.
            -Ktrap=none (default)   Disable all traps.
            -Ktrap=ovf Trap on floating point overflow.
            -Ktrap=unf Trap on floating point underflow.
                          


