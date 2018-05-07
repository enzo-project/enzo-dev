Subgrid-scale (SGS) turbulence model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following parameter allow the use of an SGS turbulence model in
Enzo, see test problem 
``run/MHD/3D/StochasticForcing/StochasticForcing.enzo``.

It is recommended to not arbitrarily mix model terms, but either 
stick to one model family (nonlinear, dissipative, or scale-similarity)
or conduct additional *a priori* test first.

Best fit model coefficients based on *a priori* testing of compressible
MHD turbulence for a wide range of sonic Mach numbers (0.2 to 20) can
be found in Table II in Grete et al. (2016) Physics of Plasmas 23 062317, 
where all models are presented in more detail.

Overall, the nonlinear model (effectively parameter free) 
with an explicit 3-point stencil showed the
best performance in decaying MHD test problem, see 
Grete et al. (2017) Phys. Rev. E. 95 033206.


``UseSGSModel`` (external)
    This parameter generally turns the SGS machinery on (even though
    no SGS term is added by default as every term needs a coefficient, 
    see below). 
    1: Turn on. Default: 0

Explicit filtering
^^^^^^^^^^^^^^^^^^

All SGS models rely on the notion that they are calculated based on
filtered/resolved quantities.
The spatial discretization itself acts as one filter.
However, in shock-capturing schemes it is questionable how "resolved"
quantities are at the grid-scale.
The following three variables enable the explicit filtering of the grid-scale
quantites as they are used in the SGS terms.

See Table 1 in Grete et al. (2017) Phys. Rev. E. 95 033206 for coefficients 
of a discrete box filter.
The recommended values correspond to a discrete representation of a box filter
on a 3-point stencil.


``SGSFilterWidth`` (external)
    Width (in units of cell widths) of the discrete filter. 
    Default: 0; 
    Recommended: 2.711;

``SGSFilterStencil`` (external)
    Discrete width of filter stencil in numbers of cells. 
    Default: 0;
    Recommended: 3;
    Maximum: 7;

``SGSFilterWeights`` (external)
    Symmetic filter weights that are used in the stencil. List of four floats.
    First number corresponds to weight of central point X_i, 
    second number corresponds to weight of points X_i+1 and X_i-1, and so on.
    Default: 0. 0. 0. 0.;
    Recommended: 0.40150 0.29925 0.00000 0.0;

Nonlinear model
^^^^^^^^^^^^^^^

``SGScoeffNLu`` (external)
    Coefficient for nonlinear Reynolds stress model.
    Default: 0;
    Recommended: 1;
    
``SGScoeffNLb`` (external)
    Coefficient for nonlinear Maxwell stress (only MHD).
    Default: 0;
    Recommended: 1;

``SGScoeffNLemfCompr`` (external)
    Coefficient for nonlinear compressive EMF model (only MHD).
    Default: 0;
    Recommended: 1;

Dissipative model
^^^^^^^^^^^^^^^^^

``SGScoeffEVStarEnS2Star`` (external)
    Coefficient for traceless eddy-viscosity Reynolds stress model
    (scaled by realizability condition in the kinetic SGS energy).
    Default: 0;
    Recommended: 0.01;

``SGScoeffEnS2StarTrace`` (external)
    Coefficient for the trace of the eddy-viscosity Reynolds stress model,
    i.e., the kinetic SGS energy (derived from realizability condition).
    Default: 0;
    Recommended: 0.08;

``SGScoeffERS2M2Star`` (external)
    Coefficient for eddy-resistivity EMF model (only MHD; scaled by
    realizable SGS energies)
    Default: 0;
    Recommended: 0.012;

``SGScoeffERS2J2`` (external)
    Coefficient for eddy-resistivity EMF model (only MHD; scaled by
    Smagorinsky SGS energies)
    Default: 0;

Scale-similarity model
^^^^^^^^^^^^^^^^^^^^^^

``SGScoeffSSu`` (external)
    Coefficient for scale-similarity Reynolds stress model.
    Default: 0;
    Recommended: 0.67;
    
``SGScoeffSSb`` (external)
    Coefficient for scale-similarity Maxwell stress (only MHD).
    Default: 0;
    Recommended: 0.9;

``SGScoeffNLemfCompr`` (external)
    Coefficient for scale-similarity EMF model (only MHD).
    Default: 0;
    Recommended: 0.89;

