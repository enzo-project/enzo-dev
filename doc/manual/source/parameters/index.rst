.. _parameters:

Enzo Parameter List
===================

The following is a largely complete list of the parameters that Enzo
understands, and a brief description of what they mean. They are grouped
roughly by meaning; an alphabetical list is also available. Parameters for
individual test problems are also listed here.

This parameter list has two purposes. The first is to describe and explain the
parameters that can be put into the initial parameter file that begins a run.
The second is to provide a comprehensive list of all parameters that the code
uses, including those that go into an output file (which contains a complete
list of all parameters), so that users can better understand these output
files.

The parameters fall into a number of categories:

**external**
    These are user parameters in the sense that they can be set in the
    parameter file, and provide the primary means of communication
    between Enzo and the user.
**internal**
    These are mostly not set in the parameter file (although strictly
    speaking they can be) and are generally used for program to
    communicate with itself (via the restart of output files).
**obsolete**
    No longer used.
**reserved**
    To be used later.

Generally the external parameters are the only ones that are modified or set,
but the internal parameters can provide useful information and can sometimes be
modified so I list them here as well. Some parameters are true/false or on/off
boolean flags.  Eventually, these may be parsed, but in the meantime, we use the
common convention of 0 meaning false or off and 1 for true or on.

This list includes parameters for the Enzo 2.3 release.

.. toctree::
   :maxdepth: 2

   initialization.rst
   io.rst       
   hierarchy.rst
   gravity.rst  
   hydro.rst    
   cooling.rst  
   particles.rst
   starform.rst 
   radiation.rst
   cosmology.rst
   bhform.rst   
   shocks.rst   
   conduction.rst
   analysis.rst 
   other.rst    
   problemtypes.rst
