lcatest Tutorial
================

This is a somewhat redundant tutorial about how to set up lcatest
on a new machine, and run a test. Now, I say redundant, because
most of this is in the README at the top level of the lcatest
source directory. Of course, if you don't know where to check out
the code, then you can't read the README...

Get lcatest
-----------

To get lcatest, head over to the
` http://lca.ucsd.edu/projects/lcatest lcatest homepage <http://lca.ucsd.edu/projects/lcatest%20lcatest%20homepage>`_,
and download the latest version.

Configure lcatest
-----------------

Now we need to setup lcatest for the machine we're on. In this
case, we have to creating a new configuration. Most of this was
done by James, and don't be surprised if you need his help to get
started.

We'll do this the classic way: copy over a config file for a
different machine and modify it.

::

    [rpwagner@co-login1 rpwagner]$ cd machines/
    [rpwagner@co-login1 rpwagner]$ mkdir ncsa-cobalt
    [rpwagner@co-login1 rpwagner]$ cp sdsc-datastar/????* ncsa-cobalt/
    [rpwagner@co-login1 rpwagner]$ cd ncsa-cobalt/
    [rpwagner@co-login1 ncsa-cobalt]$ ls
    batch.template  interactive.template  sdsc-datastar.config
    [rpwagner@co-login1 ncsa-cobalt]$ mv sdsc-datastar.config ncsa-cobalt.config 

The next step is modify the configuration and batch templates to
match this machine. Here's the current version, still showing the
`DataStar? </wiki/DataStar>`_ setup.

lcatest Configuration File
~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    [rpwagner@co-login1 ncsa-cobalt]$ cat ncsa-cobalt.config
    cpp_path /usr/lib/cpp -P 
    hdf5_inc -I/usr/local/apps64/hdf5/include
    hdf5_lib -L/usr/local/apps64/hdf5/lib -lhdf5
    hdf4_lib -L/usr/local/apps64/hdf4/lib -lz -lsz
    cxx xlC
    num_nodes 272
    num_procs 2176
    num_cores 2176
    num_procs_per_node 8
    num_cores_per_proc 1
    
    
    [rpwagner@co-login1 ncsa-cobalt]$ 

And now for a really fun question: should you set up your test
environment first, and watch the configure and build tests fail;
or, should you try to configure your code, so that you'll have some
settings to configure lcatest? For this one, we happen to have an
answer, since James configured enzo for Cobalt before I started
testing.

So, we were able to look at the machine configuration file
$ENZO\_ROOT/amr\_mpi/src/Make.mach.ncsa-cobalt to find some of the
settings.

::

    [rpwagner@co-login1 ncsa-cobalt]$ emacs -nw ncsa-cobalt.config 
    [rpwagner@co-login1 ncsa-cobalt]$ cat ncsa-cobalt.config 
    cpp_path /lib/cpp -P 
    hdf5_inc -I/usr/apps/hdf5/hdf5-1.6.5/include
    hdf5_lib -L/usr/apps/hdf5/hdf5-1.6.5/lib -lhdf5
    hdf4_lib -L/usr/apps/hdf4/HDF4.2r1/lib -lz -lsz
    cxx icc
    num_nodes 512
    num_procs 1024
    num_cores 1024
    num_procs_per_node 2
    num_cores_per_proc 1
    
    
    [rpwagner@co-login1 ncsa-cobalt]$ 

lcatest Job Templates
---------------------

lcatest uses two templates to run:

interactive.template
    Used to run tests on an interactive node.

batch.template
    Used to run tests via a batch job. (Not used currently.)

interactive.template
~~~~~~~~~~~~~~~~~~~~

Again, here's the current contents of the template, with the
settings from `DataStar? </wiki/DataStar>`_. As you might guess,
the things that look like [FOO]/GOO are variables set by lcatest
when submitting the job.

::

    [rpwagner@co-login1 ncsa-cobalt]$ cat interactive.template 
    poe [EXE_DIR]/EXE ARGS -nodes [NODES] -tasks_per_node [PROCS_PER_NODE] -rmpool 1

And now for our first crack at a template for a
a class

