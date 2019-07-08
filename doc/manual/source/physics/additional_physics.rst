.. _additional_physics:


Additional Physics
==================



Thermal Conduction
------------------

Enzo supports both isotropic and anisotropic thermal conduction within the gas.
The isotropic conduction was use in `Smith et al. (2013)
<http://adsabs.harvard.edu/abs/2013ApJ...778..152S>`_. 
The anisotropic conduction is directionally split.
The methods are adapted from `Parrish and Stone, 2005
<http://adsabs.harvard.edu/abs/2005ApJ...633..334P>`_ and the parameters can be
found under :ref:`conduction`.

Subgrid Turbulence Models
-------------------------

This will be updated soon.


Fuzzy Dark Matter Model
-----------------------

This module was written by Xinyu Li and is described in `Li, Hui &
Bryan (2019) <https://ui.adsabs.harvard.edu/abs/2019PhRvD..99f3509L/abstract>`_.
There are a few 1D tests (LightBoson) and a facility to generate and
read in a Gaussian density with wavefunction (FDMCollapse).  There is
also a framework to read in cosmology fields if one had a mechanism
to generate them, but no cosmological field generator is provided.

Turbulence Driving Modules
--------------------------

This will be updated soon.


