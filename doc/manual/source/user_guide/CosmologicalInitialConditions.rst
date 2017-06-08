.. _CosmologicalInitialConditions:

Creating Cosmological Initial Conditions
========================================

Enzo does not, by itself, generate cosmological initial conditions.
Instead, it relies on other packages to do so.  Here we describe a
variety of external packages, beginning with ``MUSIC``, a widely used
and flexible method for generating "uniform" or "zoomed" initial
condition files that can be read by Enzo.  We also describe the
original mechanism, ``inits``, has long been distributed with Enzo.
It is exclusively serial.  We also now distribute ``mpgrafic`` with
modifications to support Enzo data formats.

.. _using_music:

Using MUSIC
-----------
The MUSIC package is described in `Hahn & Abel (2011)
<http://arxiv.org/abs/1103.6031>`_ and is available from a bit bucket
repository `here <https://bitbucket.org/ohahn/music>`_.   After you
download and compiled it (note that MUSIC has its own depenencies
including FFTW and GSL), you should have the ``MUSIC`` executable.  

To use this to generate an enzo initial condition file, you can either
modify the MUSIC example parameter file, or use the sample one here:
provided :download:`music_example.conf
<./samples/music_example.conf>`  which can be used to set a small
(16^3) unigrid run.  To generate the initial conditions file, use:

::

   ./MUSIC music_example.conf

which should generate a directory called ``ic.enzo``, containing a
number of files.  MUSIC will generate the baryon densities, velocities
as well as dark matter particle displacements and velocities.
Finally, it will also generate a rudimentary enzo parameter file
``parameter_file.txt``.  You can modify this (text) file to change or
add any enzo parameters (although it does contain all the basic
settings to set the grid and cosmology parameters correctly.  To run
the enzo simulation, you will need to copy ``enzo.exe`` into that
directory and run enzo:

::
   
   [mpirun] enzo.exe parameter_file.txt

A more complicated MUSIC initial parameter file example for a zoomed
simulation (one in which one or more static refined regions are
created from the beginning, designed to model a small region at high
resolution) can be found here: :download:`music_example_zoom.conf
<./samples/music_example_zoom.conf>`

.. _using_inits:

Using inits
-----------
Inits can also generate 2D cosmological initial conditions.
The inits program uses one or more ASCII input files to set
parameters, including the details of the power spectrum, the grid
size, and output file names. Each line of the parameter file is
interpreted independently and can contain only a single parameter.
Parameters are specified in the form:

::

     ParameterName = VALUE

Spaces are ignored, and a parameter statement must be contained on
a single line. Lines which begin with the pound symbol (#) are
assumed to be comments and ignored.

First, set the parameters in the file. There are a large number of
parameters, but many don't need to be set since reasonable default
values are provided. Modifying a provided example (see the sample
files below) is probably the easiest route, but for reference there is
a list of the parameters, their meanings, and their default values.

Generating a single grid initialization (for simple Enzo runs) is
relatively straightforward. Generating a multi-grid initialization for
Enzo is somewhat more complicated, and we only sketch the full
procedure here.

Single Grid Initialization
++++++++++++++++++++++++++

To run a single grid initialization, you must set at least the
following parameters: Rank, GridDims, ParticleDims, as well as the
appropriate Cosmology and Power Spectrum parameters. A sample
parameter file is available, which sets up a small, single grid
cosmology simulation (that is, single grid for the initial
conditions, once Enzo is used, additional grids will be created).
It is available here: :download:`gas_plus_dm.inits
<./samples/gas_plus_dm.inits>` and the corresponding enzo parameter
file that can be used to run a simulation with those initial
conditions is here: 

After creating or modifying a parameter file, and compiling inits,
run the code with: :download:`gas_plus_dm_amr_adia.enzo
<./samples/gas_plus_dm_amr_adia.enzo>` 

::

     inits [-d] parameter_file

Where parameter_file is the name of your modified parameter file
(the -d turns on a debug option). This will produce a number of HDF
files containing the initial grids and particles, which are in the
correct units for use in Enzo.

For example, to use the provided sample files, you would use the
following commands:

::

   inits gas_plus_dm.inits
   [mpirun] enzo gas_plus_dm_amr_adia.enzo


Multiple-grid Initialization
++++++++++++++++++++++++++++
.. versionadded:: 2.1

The multi-grid (or nested) initialization can be used to refine in a
specific region, such as the Lagrangian sphere of a halo.  We assume
that you have first run a single-grid simulation and identified a
region out of which a halo will form and can put this in the form of
the left and right corners of a box which describes the region.  Then
you add the following parameters to the single-grid initialization
code:

::

     MaximumInitialRefinementLevel = 2
     RefineRegionLeftEdge          = 0.15523 0.14551 0.30074
     RefineRegionRightEdge         = 0.38523 0.37551 0.53074
     NewCenterFloat                = 0.270230055 0.260508984 0.415739357
     AutomaticSubgridBuffer        = 4

MaximumInitialRefinementLevel indicates how many extra levels you want
to generate (in this case two additional levels, or 3 in total,
including the root grid).  The next two parameters
(RefineRegionLeftEdge and RefineRegionRightEdge) describe the region
to be refined.  The fourth (optional) parameter re-centers the grid on
the halo to be resimulated.  The fifth parameter (AutomaticSubgridBuffer)
indicates how many course cells should be added around each refined
region.

Once you have added these parameters, run inits once on the new
parameter file in the standard way:  

::

     inits [-d] MultiGridParameterFile

It will give you a progress report as it runs (note
that if MaximumInitialRefinementLevel is large, this can take a long
time), and generate all of the necessary files (e.g.  GridDensity.0,
GridDensity.1, etc.).

It will also generate a file called EnzoMultigridParameters which you
can then copy directly into the enzo parameter file, and it specifies
the positions of the new grids.  You will still need to set a few
other parameters in the enzo parameter file, including
RefineRegionLeftEdge and RefineRegionRightEdge so that it only refines
in the specified region (typically this should match the most refined
initial grid).  Also set the MaximumRefinementLevel parameter and the
parameter controlling the density to be refined
(MinimumOverDensityForRefinement -- this also applies to the root
grid, so it needs to be divided by 8^l where l is the value of
MaximumInitialRefinementLevel).


Note that it is also possible to generate each level of initial
conditions manually.  This should not really be necessary, but a rough
guideline is given here.  To do this, prepare multiple parameter file
describing the individual parameter regions, and then top grid can be
generated with:

::

     inits [-d] -s SubGridParameterFile TopGridParameterFile

The -s flag provides the name of the sub-grid parameter file, which
is required by inits so that the particles are not replicated in
the sub-grid region. The sub-grids are made with the usual command
line:

::

     inits [-d] SubGridParameterFile

!Subgrids with MaxDims of 512 or larger will take some time and
require a fair amount of memory since the entire region is
generated and then the desired section extracted.

Inits Parameter List
++++++++++++++++++++

Cosmology Parameters
~~~~~~~~~~~~~~~~~~~~

**CosmologyOmegaMatterNow**
    This is the contribution of all non-relativistic matter (including
    HDM) to the energy density at the current epoch (z=0), relative to
    the value required to marginally close the universe. It includes
    dark and baryonic matter. Default: 1.0
**CosmologyOmegaLambdaNow**
    This is the contribution of the cosmological constant to the energy
    density at the current epoch, in the same units as above. Default:
    0.0
**CosmologyOmegaWDMNow**
    This is the contribution due to warm dark matter alone. Ignored
    unless PowerSpectrumType = 13 or 14. Default: 0.0
**CosmologyOmegaHDMNow**
    This is the contribution due to hot dark matter alone. Default: 0.0
**CosmologyOmegaBaryonNow**
    The baryonic contribution alone. Default: 0.06
**CosmologyComovingBoxSize**
    The size of the volume to be simulated in Mpc/h (at z=0). Default:
    64.0
**CosmologyHubbleConstantNow**
    The Hubble constant at z=0, in units of 100 km/s/Mpc. Default: 0.5
**CosmologyInitialRedshift**
    The redshift for which the initial conditions are to be generated.
    Default: 20.0

Power Spectrum Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~

**PowerSpectrumType**
    This integer parameter indicates the routine to be used for
    generating the power spectrum. Default: 1 The following are
    currently available:
    
    -  1 - CDM approximation from BBKS (Bardeen et al 1986) as modified
       by Peacock and Dodds (1994), to include, very roughly, the effect
       of baryons. This should not be used for high baryon universes or
       for simulations in which precision in the PS is important.
    -  2 - CHDM approximate PS from Ma (1996). Roughly good for hot
       fractions from 0.05 to 0.3.
    -  3 - Power-law (scale-free) spectra.
    -  4 - Reads in a power-spectrum from a file (not working).
    -  5 - CHDM approximate PS from Ma (1996), modified for 2 equal
       mass neutrinos.
    -  6 - A CDM-like Power spectrum with a shape parameter (Gamma),
       that is specified by the parameter PowerSpectrumGamma.
    -  11 - The Eisenstein and Hu fitting functions for low and
       moderate baryon fraction, including the case of one massive
       neutrino.
    -  12 - The Eisenstein and Hu fitting functions for low and
       moderate baryon fraction, for the case of two massive neutrinos.
    -  13 - A Warm Dark Matter (WDM) power spectrum based on the
       formulae of Bode et al. (2001 ApJ 556, 93). The WDM equivalent of
       the Eisenstein & Hu fitting function with one massive neutrino (so
       a WDM version of #11).
    -  14 - A Warm Dark Matter (WDM) power spectrum based on the
       formulae of Bode et al. (2001 ApJ 556, 93). The WDM equivalent of
       the CDM BBKS approximation of Bardeen et al 1986 (the WDM version
       of #1).
    -  20 - A transfer function from CMBFast is input for this option,
       based on the filenames described below.


**PowerSpectrumSigma8**
    The amplitude of the linear power spectrum at z=0 as specified by
    the rms amplitude of mass-fluctuations in a top-hat sphere of
    radius 8 Mpc/h. Default: 0.6
**PowerSpectrumPrimordialIndex**
    This is the index of the mass power spectrum before modification by
    the transfer function. A value of 1 corresponds to the scale-free
    primordial spectrum. Default: 1.0.
**PowerSpectrumRandomSeed**
    This is the initial seed for all random number generation, which
    should be negative. The random number generator (Numerical Recipes
    RAN3) is machine-independent, so the same seed will produce the
    same results (with other parameters unchanged). Note also that
    because the spectrum is sampled strictly in order of increasing
    k-amplitude, the large-scale power will be the same even if you
    increase or decrease the grid size. Default: -123456789
**PowerSpectrumkcutoff**
    The spectrum is set to zero above this wavenumber (i.e. smaller
    scales are set to zero), which is in units of 1/Mpc. It only works
    for power spectrum types 1-6. A value of 0 means no cutoff.
    Default: 0.0
**PowerSpectrumkmin/kmax**
    These two parameters control the range of the internal lookup table
    in wavenumber (units 1/Mpc). Reasonably sized grids will not
    require changes in these parameters. Defaults: kmin = 1e-3, kmax =
    1e+4.
**PowerSpectrumNumberOfkPoints**
    This sets the number of points in the PS look-up table that is
    generated for efficiency purposes. It should not require changing.
    Default: 10000.
**PowerSpectrumFileNameRedshiftZero**
    For input power spectra, such as those from CMBFAST, two transfer
    functions are required: one at z=0 to fix the amplitude (via
    Sigma8) and the other at the initial redshift to give the shape and
    amplitude relative to z=0. No default.
**PowerSpectrumFileNameInitialRedshift**
    see above.
**PowerSpectrumGamma**
    The shape parameter (Omega\*h); ignored unless PowerSpectrumType =
    6.
**PowerSpectrumWDMParticleMass**
    The mass of the dark matter particle in KeV for the Bode et al.
    warm dark matter (WDM) case. Ignored unless PowerSpectrumType = 13
    or 14. Default: 1.0.
**PowerSpectrumWDMDegreesOfFreedom**
    The number of degrees of freedom of the warm dark matter particles
    for the Bode et al. warm dark matter model. Ignored unless
    PowerSpectrumType = 13 or 14. Default: 1.5.
**PowerSpectrumGamma**
    The shape parameter (Omega\*h); ignored unless PowerSpectrumType =
    6.

Grid Parameters: Basic
~~~~~~~~~~~~~~~~~~~~~~

**Rank**
    Dimensionality of the problem, 1 to 3 (warning: not recently tested
    for Rank !=2). Default: 3
**GridDims**
    This sets the actual dimensions of the baryon grid that is to be
    created (and so it may be smaller than MaxDims in some cases).
    Example: 64 64 64 No default.
**ParticleDims**
    Dimensions of the particle grid that is to be created. No default.
**InitializeGrids**
    Flag indicating if the baryon grids should be produced (set to 0 if
    inits is being run to generate particles only). Default: 1
**InitializeParticles**
    Flag indicating if the particles should be produced (set to 0 if
    inits is being run to generate baryons only). Default: 1
**ParticlePositionName**
    This is the name of the particle position output file. This HDF
    file contains one to three Scientific Data Sets (SDS), one for
    dimensional component. Default: ParticlePositions
**ParticleVelocityName**
    The particle velocity file name, which must(!) be different from
    the one above, otherwise the order of the SDS's will be incorrect.
    Default: ParticleVelocities
**ParticleMassName**
    This is the name of the particle mass file, which is generally not
    needed (enzo generates its own masses if not provided). Default:
    None
**GridDensityName**
    The name of the HDF file which contains the grid density SDS. Default: 
    GridDensity
**GridVelocityName**
    The name of the HDF file which contains the SDS's for the baryonic
    velocity (may be the same as GridDensityName). Default:
    GridVelocity

Grid Parameters: Advanced
~~~~~~~~~~~~~~~~~~~~~~~~~

**MaximumInitialRefinementLevel**
    Used for multi-grid (nested) initial code generation.  This
    parameter speciesi the level (0-based) that the initial conditions
    should be generated to.  So, for example, setting it to 1
    generates the top grid and one additional level of refinement.
    Note that the additional levels are nested, keeping at least one
    coarse cell between the edge of a coarse grid and its refined grid.
    Default: 0
**RefineRegionLeftEdge, RefineRegionRightEdge**
    Species the left and right corners of the region that should be
    refined using the AutomaticSubgridGeneration method (see above
    parameter).  Default: 0 0 0 - 1 1 1
**NewCenterFloat**
    Indicates that the final grid should be recenter so that this point
    is the new center (0.5 0.5 0.5) of the grid.
**AutomaticSubgridBuffer**
    For multi-grid (nested) initial code generation (with the above
    parameters).  This parameter controls how many coarse cells are
    added around each refined region as buffer zones.  The value
    of 1 is probably ok, but larger values (4?) are probably safer.
    Default: 1
**MaxDims**
    All dimensions are specified as one to three numbers deliminated by
    spaces (and for those familiar with the KRONOS or ZEUS method of
    specifying dimensions, the ones here do not include ghost zones).
    An example is: 64 64 64. MaxDims are the dimensions of the
    conceptual high-resolution grid that covers the entire
    computational domain. For a single-grid initialization this is just
    the dimension of the grid (or of the particle grid if there are
    more particles than grid points). For multi-grid initializations,
    this is the dimensions of the grid that would cover the region at
    the highest resolution that will be used. It must be identical
    across all parameter files (for multi-grid initializations). The
    default is the maximum of GridDims or ParticleDims, whichever is
    larger (in other words unless you are using a multi-grid
    initialization, this parameter does not need to be set). Confused
    yet?
**GridRefinement**
    This integer is the sampling, for the baryon grid, in each
    dimension, relative to MaxDims. For single-grid initializations,
    this is generally 1. For multi-grids, it is the refinement factor
    relative to the finest level. In other words, if the grid covered
    the entire computational region, then each value in MaxDims would
    equal GridDims times the GridRefinement factor. Default: 1
**ParticleRefinement**
    Similar function as above, but for the particles. Note that it can
    also be used to generate fewer particles than grids (i.e. the
    GridRefinement and ParticleRefinement factors do not have to be the
    same). Default: 1
**StartIndex**
    For single-grid initializations, this should be the zero vector.
    For multi-grid initializations it specifies the index (a triplet of
    integers in 3D) of the left-hand corner of the grid to be
    generated. It is specified in terms of the finest conceptual grid
    and so ranges from 0 to MaxDims-1. Note also that for AMR, the
    start and end of a sub-grid must lie on the cell-boundary of it's
    parent. That means that this number must be divisible by the
    Refinement factor. The end of the sub-grid will be at index:
    StartIndex + GridRefinement\*GridDims. The co-ordinate system used
    by this parameter is always the unshifted one (i.e. it does not
    change if NewCenter is set).

Using mpgrafic
--------------
.. versionadded:: 2.0

This version of mpgrafic is a modified version of the public version of
mpgrafic, found at

` http://www2.iap.fr/users/pichon/mpgrafic.html <http://www2.iap.fr/users/pichon/mpgrafic.html>`_

to produce files readable by Enzo. It has been modified to write HDF5 files in
parallel.

Dependencies
++++++++++++


-  HDF5 with parallel and FORTRAN support (flags --enable-parallel
   --enable-fortran)
-  FFTW v2 with MPI support and different single and double
   precision versions. It must be compiled once for single precision
   and another time for double precision. For the former, use the
   flags --enable-mpi --enable-type-prefix --enable-float. For double
   precision, use --enable-mpi --enable-type-prefix.

Approach
++++++++

Non-nested initial conditions are created only using mpgrafic.  However if the
user wants nested initial conditions, a full-resolution grid (e.g. 256\
:sup:`3`\  grid for a 64\ :sup:`3`\  top grid with 2 nested grids) must be
created first and then post-processed with degraf to create a degraded
top-level grid and cropped (and degraded if not the finest level)
grids for the nested grids.

As with the original inits Enzo package, the baryon density and velocities are
written in a 3 dimensional array. The original inits writes the particle data
in 1-d arrays. In mpgrafic, only the particle velocities are written in a 3-d
array. Enzo has been modified to create the particle positions from the
Zel'dovich approximation from these velocities, so it is not needed to write
the positions anymore. Also it does not create particles that are represented
by a finer grid at the same position.

One big benefit of writing the particle velocities in a 3-d array is avoiding
the use of the RingIO tool because each processor knows which subvolume to read
within the velocity data.

As of HDF5 version 1.8.2, there exists a bug that creates corrupted datasets
when writing very large (e.g. >2048\ :sup:`3`\ ) datasets with multiple
components (4-d arrays). The HDF5 I/O in mpgrafic works around this bug by
creating one file per velocity component for both the baryons and particles.

How to run
++++++++++

First the user needs to compile both mpgrafic and degraf. The
configure / make systems are set up similarly.

**Configure flags:**

--enable-enzo
    turns on I/O for Enzo
--enable-double
    creates files in double precision
--enable-onedim
    creates one file per velocity component
--with-hdf=HDF5_DIR
    sets directory for parallel HDF5

If FFTW is not present in the user's library path, the following
variables must be also set

::

      CFLAGS="-I ${FFTW_DIR}/include"
      FCFLAGS="-I ${FFTW_DIR}/include"
      LDFLAGS="-L ${FFTW_DIR}/lib"

To run in parallel, you can use FC=mpif90 and LD=h5pfc, which the
compiler wrapper for parallel HDF5.

**Example configure (for Mac OSX):**

::

    ./configure LD="-bind_at_load" FC=mpif90 CC=mpicc --enable-enzo \
    --enable-double --enable-onedim --with-hdf=/usr/local/hdf5/1.8.2p

Example configure scripts can be found in mpgrafic/mpgrafic-0.2/conf.\*. After
a successful configure, you can make mpgrafic or degraf by typing 'make'.

After the programs are compiled, you make the initial conditions by using a
python script, make_ic.py, in the top directory that simplifies the user input
into mpgrafic and degraf and the moving of files.

make_ic.py parameters
~~~~~~~~~~~~~~~~~~~~~~

nprocs
    number of processors
boxsize
    box size in comoving Mpc (not Mpc/h)
resolution
    top-level grid resolution
n_levels
    level of the finest nested grid
inner_width
    width of the finest nested grid
buffer_cells
    number of cells separating nested grids
seed
    random seed (must be 9 digits)
name
    name of the data directory (saved in mpgrafic/data/name/)
center
    how much to shift the data in order to center on a particular
    region.
LargeScaleCorrection
    whether to use a noise file from a lower-resolution run
LargeScaleFile
    noise file from that lower-resolution run
OneDimPerFile
    whether we're using one file per velocity component
omega_m
    Omega matter
omega_v
    Omega lambda
omega_b
    Omega baryon
h0
    Hubble constant in units of [km/s/Mpc]
sigma8
    sigma_8
n_plawslope
    slope of power spectrum

After you set your parameters, run this script with

::

    python make_ic.py 

and it will re-compile mpgrafic and (for nested grids) degraf. Then it will run
mpgrafic for the full-resolution box. If the user wants nested grids, it will
copy the data files to mpgrafic/degraf and create the set of nested grid files.

The user cannot specify the initial redshift because mpgrafic determines it
from the parameter sigstart that is the maximum initial density fluctuation.
From this, mpgrafic calculates the initial redshift. This file is overwritten
by the python script, so if you want to change this parameter, change it in the
python script (routine write_grafic1inc).

The noise file is always kept in mpgrafic/mpgrafic-0.2/src and is named
$seed_$resolution.dat, where $resolution is the top-level grid resolution. It
can be re-used with LargeScaleFile if the user wants to re-simulate the volume
at a higher resolution.

The data files are moved to mpgrafic/data/$name. If nested grids were created,
degraf writes a set of parameters in enzo.params for copy-pasting into an Enzo
parameter file. Now you can move the files to the simulation directory and
start your Enzo cosmology simulation!


