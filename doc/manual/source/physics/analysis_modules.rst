.. _analysis_modules:

Analysis Modules
================



Tracer Particles
----------------
Enzo has extensive analysis methods to track tracer particles. Details about this will be updated soon.



Inline Halo Finder
------------------

*Source:  FOF.C*

Enzo has in-built Friends-of-Friends halo finder to identify dark matter halos, originally written by Volker Springel. Dark matter halos are identified by linking "friend" dark matter particles that lie within a specified linking length. All the output files are written in the directory FOF/.  See :ref:`inline_analysis` for more details on parameters. 



Embedded Python
---------------

*Source:  InitializePythonInterface.C*

Python can now be embedded inside Enzo, for inline analysis as well as interaction. See :ref:`embedded-python` for the details about enabling it and using it (out of date).



