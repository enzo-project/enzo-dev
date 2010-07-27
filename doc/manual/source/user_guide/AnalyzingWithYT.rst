Analyzing With YT
=================

What is YT?
-----------

YT is a python-based tool designed for analyzing and visualizing
Adaptive Mesh Refinement data, specifically as output from Enzo. YT
is completely free and open source, with an active and expanding
development community, and it presents to the user both high-level
and low-level APIs. The
` documentation <http://yt.enzotools.org/doc/>`_ contains a
tutorial as well as an API reference, but here we will step through
some simple steps toward creating script to make simple plots of a
cosmological simulation.

This brief tutorial presupposes that you have run the
`installation script? </wiki/Devel/UserGuide/BuildingEnzo#YT>`_ and
are comfortable launching python. (The install script will tell you
how!) It's also encouraged to launch the special YT-enhanced
` IPython <http://ipython.scipy.org/>`_ shell via the command iyt,
which (thanks to IPython!) features filesystem navigation and tab
completion, along with interactive plotting capabilities.

Making Slices
-------------

Here is a sample script that will make a set of slices centered on
the maximum density location, with a width of 100 kpc.

::

    from yt.mods import *
    pf = EnzoStaticOutput("RedshiftOutput0035.dir/RedshiftOutput0035")
    
    pc = raven.PlotCollection(pf)
    pc.add_slice("Density",0)
    pc.add_slice("Density",1)
    pc.add_slice("Density",2)
    pc.set_width(100.0,'kpc')
    pc.save("z35_100kpc")

If you put this into a file called my\_script.py, you can execute
it with python2.5 my\_script.py and it will save out a set of
images prefixed with z35\_100kpc in PNG format.

Making Simple Radial Profiles
-----------------------------

If you want to make radial profiles, you can generate and plot them
very easily with YT. Here is a sample script to do so.

::

    from yt.mods import *
    pf = EnzoStaticOutput("RedshiftOutput0035.dir/RedshiftOutput0035")
    
    pc = PlotCollection(pf)
    
    pc.add_profile_sphere(100.0, 'kpc', ["Density", "Temperature"])
    pc.save("z35_100kpc")
    
    pc.switch_z("VelocityMagnitude")
    pc.save("z35_100kpc")

To show the mass distribution in the Density-Temperature plane, we
would make a phase diagram.

::

    from yt.mods import *
    pf = EnzoStaticOutput("RedshiftOutput0035.dir/RedshiftOutput0035")
    
    pc = PlotCollection(pf)
    
    pc.add_phase_sphere(100.0, 'kpc', ["Density", "Temperature", "CellMassMsun"], weight=None)
    pc.save("z35_100kpc")

More Analysis Tasks
-------------------

YT has
` many more features <http://yt.enzotools.org/doc/intro.html#what-functionality-does-yt-offer>`_
than just these, including halo finding, projections, oblique
slices, a wide variety of derived fields, many derived quantities
and so on and so forth. On the documentation site you will find a
` quick guide <http://yt.enzotools.org/doc/quick_guide/index.html>`_,
a
` longer tutorial <http://yt.enzotools.org/doc/tutorial/index.html>`_,
and an
` API reference <http://yt.enzotools.org/doc/modules/index.html>`_.

YT is under active development with a release forthcoming featuring
fully-parallel analysis capabilities, along with lightweight
frontends to 3D visualization toolkits like VTK and S2Plot. Please
feel free to visit ` the website <http://yt.enzotools.org/>`_,
check out the ` documentation <http://yt.enzotools.org/doc/>`_ and
ask questions on the
a class="ext-link"
href="http://lists.spacepope.org/listinfo.cgiide/RunningEnzo,
Devel/UserGuide/EnzoTestSuite, Devel/UserGuide/RunningInits,
Devel/UserGuide/RunningMPgrafic, Devel/UserGuide/EnzoOutputFormat,
Devel/UserGuide/AnalyzingWithYT, Devel/UserGuide/HierarchyFile,
Devel/UserGuide/EnzoInternalUnits,
Devel/UserGuide/EnzoParticleMass, Devel/UserGuide/GridFieldArrays,
Devel/UserGuide/FlowChart)?

