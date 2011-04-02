Embedded Python
===============

Python can now be embedded inside Enzo, for inline analysis as well as
interaction. This comes with several shortcomings, but some compelling strong
points.

How To Compile
--------------

The configure option that controls compilation of the Python code
can be toggled with

::

    make python-yes

or to turn it off,

::

    make python-no

This will look for the following variables in the machine-specific Makefile:

::

    MACH_INCLUDES_PYTHON
    MACH_LIBS_PYTHON

for an example of how to define these variables, see
Make.mach.orange in the source repository.

How it Works
------------

On Enzo startup, the Python interface will be initialized. This constitutes the
creation of an interpreter within the memory-space of each Enzo process, as
well as import and construct the `NumPy <http://numpy.scipy.org/>`_ function
table. Several Enzo-global data objects for storing grid parameters and
simulation parameters will be initialized and the Enzo module will be created
and filled with those data objects.

Once the Python interface and interpreter have finished initializing, the
module user_script will be imported -- typically this means that a script named
``user_script.py`` in the current directory will be imported, but it will
search the entire import path as well. Every ``PythonSubcycleSkip`` subcycles,
at the bottom of the hierarchy in ``EvolveLevel.C`` the entire grid hierarchy
and the current set of parameters will be exported to the Enzo module and then
user_script.main() will be called.

How to Run
----------

By constructing a script inside ``user_script.py``, the Enzo hierarchy can be
accessed and modified. The analysis toolkit `yt <http://yt.enzotools.org/>`_
has functionality that can abstract much of the data-access and handling.
Currently several different plotting methods -- profiles, phase plots, slices
and cutting planes -- along with all derived quantities can be accessed and
calculated. Projections cannot yet be made, but halo finding can be performed
with Parallel HOP only. The following script is an example of a script that
will save a slice as well as print some information about the simulation. Note
that, other than the instantiation of ``lagos.EnzoStaticOutputInMemory``, this
script is identical to one that would be run on an output located on disk.

Recipes and convenience functions are being created to make every aspect of
this simpler.

::

    from yt.mods import *
    
    def main():
         pf = lagos.EnzoStaticOutputInMemory()
         pc = PlotCollection(pf)
         pc.add_slice("Density", 0)
         pc.save("%s" % pf)
         v, c = pf.h.find_max("Density")
         sp = pf.h.sphere(c, 1.0/pf['mpc'])
         totals = sp.quantities["TotalQuantity"](["CellMassMsun","Ones"], lazy_reader=True)
         print "Total mass within 1 mpc: %0.3e total cells: %0.3e" % (totals[0], totals[1])

Which Operations Work
---------------------

The following operations in yt work:

  * Derived quantities
  * Slices
  * Cutting planes
  * Fixed Resolution Projections (i.e., non-adaptive)
  * 1-, 2-, 3-D Profiles

This should enable substantial analysis to be conducted in-line.  Unfortunate
adaptive projections require a domain decomposition as they currently stand (as
of yt-1.7) but this will be eliminated with a quad-tree projection method
slated to come online in yt-2.0.  In future versions of yt the volume rendering
approach will be parallelized using kD-tree decomposition and it will also
become available for inline processing.

Please drop a line to the yt or Enzo mailing lists for help with any of this!

Things Not Yet Done
-------------------

-  Adaptive Projections do not work.
-  Particles are not yet exported correctly
-  Speed could be improved, but should be extremely efficient for a small
   number of grids.  Future versions will utilize intercommunicators in MPI to
   allow for asynchronous analysis.

