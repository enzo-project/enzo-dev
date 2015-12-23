.. _mhd_methods:

MHD Methods
===========

Dedner cleans :math:`\nabla \cdot B` with a wave-like hyperbolic cleaner, and is
easy to use.  

CT preserves :math:`\nabla \cdot B` to machine precision, but is slightly harder to use.


    ====================== ==================== ===============
    quantity               MHD-RK               MHD-CT
    ====================== ==================== ===============
    Reconstruction         PLM                  PLM
    Splitting              Unsplit, Runge Kutta Split (strang)
    :math:`\nabla \cdot B` Few Percent          Machine Noise
    Difficulty             Easy                 Two pitfalls 
    ====================== ==================== ===============



Use of Dedner
============= 

The Dedner method (``HydroMethod = 4``) is relatively straight forward.
The three magnetic components are stored in ``BaryonField``, with the relevant
indices found through ``IdentifyFieldQuantities``.  Since they hyperbolic
cleaner will largely eliminate any divergence, initialization an injection of
magnetic energy is straight forward.

AMR is done in the same manner as other fluid quantities in Enzo.

The method is described in `Dedner et al. 2002 JCP 175, 2, 645-673
<http://adsabs.harvard.edu/abs/2002JCoPh.175..645D>`_

The implementation and test problems can be found in `Wang & Abel 2009, ApJ 696 96 <http://adsabs.harvard.edu/abs/2009ApJ...696...96W>`_.


Use of MHD-CT
=============

Use of MHD-CT (``HydroMethod = 6``) is somewhat complicated by the staggered nature of the magnetic field.  This allows the
field to be updated by the curl of an electric field, thus preserving
:math:`\nabla \cdot B = 0` to machine precision, but requires some additional
machinery to ensure consistency of the data structure.

The magnetic field is represented by two data structures.  

- The cell centered magnetic field is accessed the same manner as done in
  ``HydroMethod==4``.  However, this is a *read only quantitiy*.  Periodically
  it is replaced by the spatial average of the face centered magnetic field.
  
- The face centered magnetic field is stored in ``MagneticField``.  This is a
  staggered structure, and is divergence free.  If a modification need to be made
  to the magnetic field when using MHDCT (such as injection by a source or
  initialization), the modification should be done in a divergence free manner.
  More details below.

The primary references are:

CT algorithms: 
`Balsara & Spicer 1999, JCP, 149, 2, 270-292
<http://adsabs.harvard.edu/abs/1999JCoPh.149..270B>`_

`Gardiner & Stone 2005, JCP, 205, 2, 509-539
<http://adsabs.harvard.edu/abs/2005JCoPh.205..509G>`_

AMR Algorithm:
`Balsara 2001 JCP, 174, 2, 614-648
<http://adsabs.harvard.edu/abs/2001JCoPh.174..614B>`_

Implementation and test problems:
`Collins, Xu, Norman, Li & Li, ApJS 186, 2, 308
<http://adsabs.harvard.edu/abs/2010ApJS..186..308C>`_.

More Details
============ 

Within the code, there are several flags to control use of magnetic fields.

``UseMHD`` is true when *either* MHD module is on.  

``UseMHDCT`` is true when ``HydroMethod=6``, and controls use of the face- and
edge- centered fields.

``HydroMethod==4`` or ``HydroMethod==MHD_RK`` controls things that only pertain
to that hydro method.  This is either control of the solver itself, or any
setting of the centered magnetic field (e.g. ``BaryonField[B1Num]``).

``HydroMethod==6`` or ``HydroMethod==MHD_Li`` typically only deals with control
of the ``MHD_Li`` solver.  Things dealing with the data structures should be
controlled with ``UseMHDCT``

Implementation details for MHDCT can be found in :ref:`mhdct_details`
