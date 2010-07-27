Enzo Output Formats
===================

Although there are a number of ways of specifying when (and how
often) enzo outputs information, there is only one type of output
'dump' (well, not quite -- there are now movie dumps, see below),
which can also be used to restart the simulation. The output format
uses the following files, each of which begins with the output
name, here we use the example base\_name, and are then followed by
the output number, ranging from 0000 to 9999 (if more than 10000
grids are generated then the number goes to 10000, etc.). When
restarting, or other times when an output filename needs to be
specified, use the name without any extension (e.g. enzo -r
base\_name0000).

Summary of Files
----------------

**base\_name0000**
    This ascii file contains a complete listing of all the parameter
    settings, both those specified in the initial parameter file, as
    well as all those for which default values were assumed. The
    parameters (see the
    `page on parameters? </wiki/Devel/UserGuide/EnzoParameters>`_ for a
    discussion) are in the same format as that used in the input file:
    parameter\_name = value. This file is modifiable if you would like
    to restart from a certain point with different parameter values.
**base\_name0000.hierarchy**
    This ascii file specifies the hierarchy structure as well as the
    names of the grid files, their sizes, and what they contain. It
    should not be modified and was not intended for general human
    consumption.
**base\_name0000.cpu00001**
    The field information for each cpu (padded with zeros) is contained
    in separate files with a root 'Node' for each grid, padded with
    zeros to be eight digits. The format is the Hierarchy Data Format
    (HDF) version 5, a self-describing machine-independent data format
    developed and supported by the National Center for Supercomputing
    Applications (NCSA). More information can be found on their
    ` home page <http://www.hdfgroup.org>`_. Most scientific
    visualization packages support this format. Each field is stored as
    it's own one-, two- or three-dimensional Scientific Data Set (SDS),
    and is named for identification. Particles (if any) are included
    with a set of one-dimensional datasets under the top 'grid' node.
**base\_name0000.boundary**
    An ascii file which specifies boundary information. It is not
    generally useful to modify.
**base\_name0000.boundary.hdf**
    Contains field-specific boundary information, in HDF format.
**base\_name0000.radiation**
    This ascii file is only generated if using the self-consistent
    radiation field.

Output Units
------------

The units of the physical quantities in the grid SDS's are depend
on the problem being run. For most test problems there is no
physical length or time specified, so they can be be simply scaled.
For cosmology there are a set of units designed to make most
quantities of order unity (so single precision variables can be
used). These units are defined below (rho0 =
3\*OmegaMatterNow\*(100\*HubbleConstantNow
km/s/Mpc)\ :sup:`2`\ /(8\*Pi\*G)).


-  length: ComovingBoxSize/HubbleConstantNow \* Mpc / (1+z)
-  density: rho0 \* (1+z)\ :sup:`3`\ 
-  time: 1/sqrt(4\*Pi\*G\*rho0\*(1+InitialRedshift)\ :sup:`3`\ )
-  temperature: K
-  velocity: (length/time)\*(1+z)/(1+InitialRedshift) (this is z
   independent)

The conversion factor is also given in the ascii output file
(base\_name0000): search for DataCGSConversionFactor. Each field
has its own conversation factor, which converts that field to cgs
units. Users can also set completely arbitrary internal units, as
long as they are self-consistent: to see how to do this, go to
a class="missing wiki" href="/wiki/Devel/UserGuide/EnzoIntern

