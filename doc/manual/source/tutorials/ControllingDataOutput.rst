Controlling Enzo data output
============================

How and when Enzo outputs data is described below.
There are five ways to control when data is output,
two output formats, and two pitfalls when
determining how to output data from your Enzo simulation.

Data Formats and Files
----------------------

There are two output formats for Enzo data. In both cases, each
data dump gets its own directory.

Each data dump writes several key files. NNNN denotes the dump
number (i.e. 0001) and basename is something like RedshiftOutput or
data or DD}.

All output files are also restart files. It's not necessarily wise
to write in 32 bit format if you're computing in 64, though, as
you'll lose all the extra precision when you restart. (These are
makefile flags.)

basenameNNNN::

    The parameter file. This contains general simulation parameters,
    dump time, cycle, and all the parameters defined here. It's worth
    your time to be familiar with what's in this file.


basenameNNNN::

    The hierarchy file. Contains a description of the hierarchy. One
    entry for each grid, including information like the Grid Size, the
    position in the volume, it's position in the hierarchy.


basenameNNNN.boundary::

    A description of the boundary (plain text.) Basically a meta
    description and filename for the next file


basenameNNNN.boundary.hdf5::

    Actually contains the boundary information.


Other versions may output other files. It would be nice to have a
description of them.

Packed AMR
~~~~~~~~~~

This is the default output format. Each processor outputs all the grids it owns.
This is done to avoid the hassle that comes with a 500,000 grid,
512 processor AMR sim. 512 files are much easier to deal with than
500,000

In addition to the parameter, hierarchy, and boundary files which
may or may not be described elsewhere, data is output in one
basenameNNNN.taskmapCCCC} file for each processor, which contains a
map between grid number and hdf5 file, and one basenameNNNN.cpuCCCC
for each processor NNNN and CCCC are the dump number and cpu
number, respectively.

basenameNNNN.cpuCCCC is an hdf5 file which contains an hdf5 group
for each grid. Each grid in turn contains a dataset for each of the
fields in the simulation.

::

    ~/DD0100>h5ls data0100.cpu0003 
    Grid00000002             Group
    Grid00000026             Group
    ~/DD0100>h5ls data0100.cpu0003/Grid00000002
    Density                  Dataset {16, 16, 32}
    z-velocity               Dataset {16, 16, 32}

Not Packed AMR
~~~~~~~~~~~~~~

The second output format is not packed. Not packed amr is included
for legacy reasons only. It writes one hdf5 file per grid. Use is
strongly discouraged. It will likely be removed in subsequent
versions.

Pitfall - Hard Coded Pathnames
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pathnames are hard coded into the parameter, boundary, and
hierarchy files. If you need to move a simulation, put the data
where you want it, and run the following script in the directory
containing the data.

::

    #!/bin/tcsh                                                                                                                        
    if ( `ls -1 *.hierarchy |wc -l ` != 1 ) exit
    set paramfile = `basename *.hierarchy .hierarchy`           # Find the parameter file.                                             
    set dir1 = `grep GlobalDir  $paramfile | awk '{print $3}'`  # Find the hard-coded directory                                        
    set dir2 = `cd ..; pwd`                                     # The working directory                                                
    foreach i ($paramfile $paramfile.hierarchy $paramfile.boundary) # Find the files to alter                                          
        sed -e "s:"$dir1":"$dir2":g" $i > tmp; mv tmp $i        # Replace the directory. sed is rad.                                   
    end

This script is included in the Enzo source distribution, as
[browser:public/trunk/bin/update\_path bin/update\_path]. Assuming
you have this in your path, here's how to loop over the all of the
subdirectories with a particular prefix, and update them. You can
modify the prefix (i.e., DD), to match other prefixes, such as RD.
Please note, this happens to be bash syntax, so you may need to
adjust it for other shells.

::

    for d in `find . -type 'd' -name 'DD*'`
    do
        cd $d
        update_path
        cd ../
    done

*If sed scares you, here's a snippet of Python to iterate over a file, replace a string and write to a new file. It could be a bit more terse, but hopefully it is clear.*

::

    hierarchy = 'RD0033.hierarchy'
    old, new = '/dsgpfs/harkness/NewL7/Dumps/RD0033/','/gpfs/ux455215/L7/RD0033/'
    new_lines = (line.replace(old, new) for line in open(hierarchy))
    open("%s.new" % hierarchy).writelines(new_lines)

Timing Methods
--------------

There are 6 ways to trigger output from enzo.

Cycle Based Output
~~~~~~~~~~~~~~~~~~

::

    CycleSkipDataDump = N
    CycleLastDataDump = W
    DataDumpName = data

One can trigger output every N cycles starting with cycle W using
CycleSkipDataDump and CycleLastDataDump. Outputs are put in the
directory DD0000 (or DD0001, etc.) and the basename is determined
by DataDumpName.

CycleSkipDataDump <= 0 means cycle based output is skipped. The
default is 0.

Pitfall 2: CycleLastDataDump defaults to zero and is incremented by
CycleSkipDataDump every time output is done. If you change the
value of CycleSkipDataDump and neglect to change CycleLastDataDump,
Enzo will dump as long as CycleNumber >= CycleSkipDataDump +
CycleLastDataDump. (So if you change CycleSkipDataDump from 0 to 10
from a Redshift dump at n=70, you'll get an output every timestep
for 7 timesteps.)

Time Based Output
~~~~~~~~~~~~~~~~~

::

    TimeLastDataDump = V
    dtDataDump = W

Exactly like Cycle based output, but triggered whenever time >=
TimeLastDataDump + dtDataDump. The same pitfall applies.

Redshift Based Output
~~~~~~~~~~~~~~~~~~~~~

::

    CosmologyOutputRedshift[ 0 ] = 12
    CosmologyOutputRedshiftName[ 0 ] = Redshift12
    RedshiftDumpName             = RedshiftOutput

Outputs at the specified redshift. Any number of these can be
specified.

CosmologyOutputRedshift[ i ] is the only necessary parameter, and
is the ith redshift to output.

Any outputs with CosmologyOutputRedshiftName[ i ] specified has
that name used for the output, and no number is appended. (so if
CosmologyOutputRedshiftName[ 6 ] = BaconHat, the outputs will be
BaconHat, BaconHat.hierarchy, etc.)

If CosmologyOutputRedshiftName[ i ] is omitted, RedshiftDumpName is
used for the basename, and the output number is taken from the
array index. (So CosmologyOutputRedshift[19] = 2.34 and
RedshiftDumpName = MonkeyOnFire, at dump will be made at z=2.34
with files called MonkeyOnFire0019.hierarchy, etc.)

Force Output Now
~~~~~~~~~~~~~~~~

The following two options are run time driven. These are especially
useful for very deep simulations that spend the majority of their
time on lower levels.

To force an output as soon as the simulation finished the next step
on the finest resolution, make a file called outputNow:

::

    touch outputNow

This will remove the file as soon as the output has finished.

Sub Cycle Based Output
~~~~~~~~~~~~~~~~~~~~~~

To get the simulation to output every 10 subsycles (again at the
finest level of resolution) put the number of subcycles to skip in
a file called subcycleCount:

::

    echo 10 > subcycleCount

Time Based Interpolated Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Even when you are running simulations with a long dtDataDump, sometimes you may
want to see or analyze the interim datadumps.  Using dtInterpolatedDataDump,
you can control enzo to check if it should start outputting interpolated data
based on the time passed (dtInterpolatedDataDump < dtDataDump).

::

    dtDataDump = 1e-4
    dtInterpolatedDataDump = 1e-5

This is mostly for making movies or looking at the interim data where the
TopGrid dt is too long, and in principle, this output shouldn't be used for
restart.

Friendly Note on Data Output
----------------------------

Enzo is content to output enough data to fill up a hard drive --
for instance, your home directory. This should be noted before
output parameters are set, particularly the Sub Cycle outputs, as
Enzo has no prohibition against causing problems with quotas and
file system size.


