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
    Difficulty             Easy                 Less Easy
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

Use of MHD-CT (``HydroMethod = 6``) is complicated by the staggered nature of the magnetic field.  This allows the
field to be updated by the curl of an electric field, thus preserving
:math:`\nabla \cdot B = 0` to machine precision, but requires some additional
machinery to ensure consistency of the data structure.

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

Enzo uses two representations of the magnetic field.
Uses two filds: staggered field, ``MagneticField`` and Centered field,
``CenteredB``.  Also uses an edge staggered field, ``ElectricField``.  ``MagneticField``, being stored on the faces of the
zones, has one additional point along each component.  For instance, if a ``grid`` had dimensions :math:`n_x, n_y, n_z` then
 :math:`B_x` will have dimensions :math:`n_x+1, n_y, n_z`.  ``ElectricField`` has additional points transverse to the direction
of the component, so :math:`E_x` has dimensions :math:`n_x, n_y+1, n_z+1`.
There are several helper variables, such as ``MagneticDims[3][3]``,
``ElectricDims[3][3]``, ``MagneticSize[3]``, and ``ElectricSize[3]`` to describe
these variables.

Dealing with the Magnetic Field
-------------------------------

``CenteredB`` should be considered a read-only quantity-- it is
replaced with a centered spatial average of ``MagneticField`` as necessary.  
``MagneticField`` should only be modified in a manner that is definitely
divergence free.  For more general initialization, one can use the function ``MHD_Curl``
for fields that can be represented by a vector potential.

Interpolation
------------- 

Interpolation must be done, obviously, in a divergence-free manner.  Balsara
2001 describes this method.  Interpolation is done on all three components of
``MagneticField`` at once.  This method only allows ``RefineBy = 2``.  

One challenge of this method is that newly interpolated regions require
knowledge of any fine-grid data at the same level that may share a face.  Thus
instead of simply interpolating from parent grids, then copying from old fine
grids, MHDCT must use the magnetic information from the old fine grids.  This is
done by first computing interpolation derivatives (done in ``Grid_MHD_CID.C``
and stored in ``DyBx``, etc) then communicating this information to the relevant
parent grids (done in ``Grid_SendOldFineGrids.C``)  This makes MHD-CT
interpolation a 3 grid interaction (Parent, Child, Old Child) rather than a 2
body interaction (Parent and Child) as all other fields.

Projection and Flux Correction
------------------------------

As with other quantities, magnetic fields need to be projected to parents, then
coarse zones next to projected zones need to be corrected to ensure
conservation.  As described by Balsara 2001, this involves area weighted
projection of face centered field on the fine grid, then a correction using the
electric field.  In order to simplify the logic and machinery, Enzo MHD-CT
actually projects the ``ElectricField``, then takes the curl over the new
magnetic field.  This is formally equivalent to projection plus flux correction,
but doesn't have as many cases to check and grid interactions to worry about.
This is done in ``EvolveLevel`` by the routine ``Grid_MHD_UpdateMagneticField``

``MHDCTUseSpecificEnergy``
--------------------------

Historically, Enzo MHDCT used conserved energy throughout, rather than specific
energy as done in the rest of Enzo.  Upon porting to Enzo2.3, this was changed,
but due to some unforeseen issues, this changes the answer.  This is provided to
ensure compatibility with old answers, and because the author is suspicious
about that which changes the answer.  This will be removed in the future.


Future Work (or, "Glitches")
----------------------------

Most neighbor searching throughout Enzo is done with the Fast Neighbor Locator,
which uses a chaining mesh to identify neighbors.  This is not done for the
communication done in ``SendOldFineGrids,`` but should be.

Additionally, both ``SendOldFineGrids`` and the electric field projection need
to be updated to work with the 3 phase non-blocking communication

In principle, the CT machinery can be used in conjunction with the MHD-RK
machinery.  Interested students can contact dcollins for further instruction.

Presently MHD-CT needs additional layers of ghost zones over the base hydro.  I
believe that I can reduce this by communicating the electric field, which will
improve memory overhead.  Again, interested parties can contact me for details.

Multi-species needs to be tested.

Presently, the centered magnetic field is stored in ``CenteredB`` rather than an
element of the ``BaryonField`` array.  This was done in order to signify that
the array is "Read Only."  This should probably be changed.

The mhd interpolation routine, ``mhd_interpolate.F``, is an abomination.  I
apologize.  I'll fix it some day.
