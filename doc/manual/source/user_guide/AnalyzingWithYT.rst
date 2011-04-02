.. _analyzing_with_yt:

Analyzing With YT
=================

What is YT?
-----------

YT is a python-based tool designed for analyzing and visualizing
Adaptive Mesh Refinement data, specifically as output from Enzo. YT
is completely free and open source, with an active and expanding
development community, and it presents to the user both high-level
and low-level APIs. The
`documentation <http://yt.enzotools.org/doc/>`_ contains a
tutorial as well as an API reference, but here we will step through
some simple steps toward creating script to make simple plots of a
cosmological simulation.

This brief tutorial presupposes that you have run the installation
script and are comfortable launching python.  (The install script will
tell you how!) It's also encouraged to launch the special YT-enhanced
`IPython <http://ipython.scipy.org/>`_ shell via the command ``iyt``,
which (thanks to IPython!) features filesystem navigation and tab
completion, along with interactive plotting capabilities.

Making Slices
-------------

Here is a sample script that will make a set of slices centered on
the maximum density location, with a width of 100 kpc.

.. code-block:: python

   from yt.mods import *
   pf = EnzoStaticOutput("RedshiftOutput0035.dir/RedshiftOutput0035")
   
   pc = raven.PlotCollection(pf)
   pc.add_slice("Density",0)
   pc.add_slice("Density",1)
   pc.add_slice("Density",2)
   pc.set_width(100.0,'kpc')
   pc.save("z35_100kpc")

If you put this into a file called ``my_script.py``, you can execute
it with ``python2.5 my_script.py`` and it will save out a set of
images prefixed with ``z35_100kpc`` in PNG format.

Making Simple Radial Profiles
-----------------------------

If you want to make radial profiles, you can generate and plot them
very easily with YT. Here is a sample script to do so.

.. code-block:: python

    from yt.mods import *
    pf = EnzoStaticOutput("RedshiftOutput0035.dir/RedshiftOutput0035")
    
    pc = PlotCollection(pf)
    
    pc.add_profile_sphere(100.0, 'kpc', ["Density", "Temperature"])
    pc.save("z35_100kpc")
    
    pc.switch_z("VelocityMagnitude")
    pc.save("z35_100kpc")

To show the mass distribution in the Density-Temperature plane, we
would make a phase diagram.

.. code-block:: python

    from yt.mods import *
    pf = EnzoStaticOutput("RedshiftOutput0035.dir/RedshiftOutput0035")
    
    pc = PlotCollection(pf)
    
    pc.add_phase_sphere(100.0, 'kpc', ["Density", "Temperature", "CellMassMsun"], weight=None)
    pc.save("z35_100kpc")

More Information
----------------

For more information on yt, see the `yt website <http://yt.enzotools.org>`_,
where you will find mailing lists, documentation, API documentation, a cookbook
and even a gallery of images.
