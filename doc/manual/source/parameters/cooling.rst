.. _cooling_parameters:

Cooling Parameters
~~~~~~~~~~~~~~~~~~

Simple Cooling Options
^^^^^^^^^^^^^^^^^^^^^^

``RadiativeCooling`` (external)
    This flag (1 - on, 0 - off) controls whether or not a radiative
    cooling module is called for each grid. There are currently several
    possibilities, controlled by the value of another flag. See :ref:`cooling` 
    for more information on the various cooling methods.  Default: 0
    
    -  If the ``MultiSpecies`` flag is off, then equilibrium cooling is
       assumed and one of the following two will happen. If the parameter
       ``GadgetCooling`` is set to 1, the primordial equilibrium code is
       called (see below). If ``GadgetCooling`` is set to 0, a file called
       ``cool_rates.in`` is read to set a cooling curve. This file consists
       of a set of temperature and the associated cgs cooling rate; a
       sample compute with a metallicity Z=0.3 Raymond-Smith code is
       provided in ``input/cool_rates.in``. This has a cutoff at 10000 K
       (Sarazin & White 1987). Another choice will be
       ``input/cool_rates.in_300K`` which goes further down to 300 K (Rosen
       & Bregman 1995).
    -  If the ``MultiSpecies`` flag is on, then the cooling rate is
       computed directly by the species abundances. This routine (which
       uses a backward differenced multi-step algorithm) is borrowed
       from the Hercules code written by Peter Anninos and Yu Zhang,
       featuring rates from Tom Abel. Other varieties of cooling are
       controlled by the ``MetalCooling`` parameter, as discused below.
``RadiativeCoolingModel`` (external)
    This switches between the tabular look up cooling that is standard (RadiativeCoolingModel=1) and an analytic fit to the Wolfire et al 2003, ApJ, 587, 278 made by Koyama and Inutsuka 2006 (RadiativeCoolingModel = 3, arXiv:astro-ph/0605528).  Default: 1
``GadgetCooling`` (external)
    This flag (1 - on, 0 - off) turns on (when set to 1) a set of
    routines that calculate cooling rates based on the assumption of a
    six-species primordial gas (H, He, no H2 or D) in equilibrium, and
    is valid for temperatures greater than 10,000 K. This requires the
    file ``TREECOOL`` to execute. Default: 0
``GadgetEquilibriumCooling`` (external)
    An implementation of the ionization equilibrium cooling code used
    in the GADGET code which includes both radiative cooling and a
    uniform metagalactic UV background specified by the ``TREECOOL`` file
    (in the ``amr_mpi/exe`` directory). When this parameter is turned on,
    ``MultiSpecies`` and ``RadiationFieldType`` are forced to 0 and
    ``RadiativeCooling`` is forced to 1.
    [Not in public release version]
``MetalCooling`` (external)
    This flag (0 - off, 1 - metal cooling from Glover & Jappsen 2007,
    2 - Cen et al (1995), 3 - Cloudy cooling from Smith, Sigurdsson, &
    Abel 2008) turns on metal cooling for runs that track
    metallicity. Option 1 is valid for temperatures between 100 K and
    10\ :sup:`8`\ K because it considers fine-structure line emission
    from carbon, oxygen, and silicon and includes the additional metal
    cooling rates from Sutherland & Dopita (1993). Option 2 is only
    valid for temperatures above 10\ :sup:`4`\ K. Option 3 uses
    multi-dimensional tables of heating/cooling values created with
    Cloudy and optionally coupled to the ``MultiSpecies``
    chemistry/cooling solver. This method is valid from 10 K to 10\
    :sup:`8`\ K. See the Cloudy Cooling parameters below.  Default: 0.
``MetalCoolingTable`` (internal)
    This field contains the metal cooling table required for
    ``MetalCooling`` option 1. In the top level directory input/, there are
    two files ``metal_cool.dat`` and ``metal_cool_pop3.dat`` that consider
    metal cooling for solar abundance and abundances from
    pair-instability supernovae, respectively. In the same directory,
    one can find an IDL routine (``make_Zcool_table.pro``) that generates
    these tables. Default: ``metal_cool.dat``
``MultiSpecies`` (external)
    If this flag (1, 2, 3- on, 0 - off) is on, then the code follows
    not just the total density, but also the ionization states of
    Hydrogen and Helium. If set to 2, then a nine-species model
    (including H2, H2+ and H-) will be computed, otherwise only six
    species are followed (H, H+, He, He+, He++, e-). If set to 3, then
    a 12 species model is followed, including D, D+ and HD. This
    routine, like the last one, is based on work done by Abel, Zhang
    and Anninos. Default: 0
``MultiMetals`` (external)
    This was added so that the user could turn on or off additional
    metal fields - currently there is the standard metallicity field
    (Metal_Density) and two additional metal fields (Z_Field1 and
    Z_Field2). Acceptable values are 1 or 0, Default: 0 (off).
``ThreeBodyRate`` (external)
    Which Three Body rate should be used for H2 formation?: 0 = Abel, Bryan, Norman 2002, 1 = PSS83, 2= CW83, 3 = FH07, 4= G08.  (Turk et al 2011 covers these)
``CIECooling`` (external)
    Should CIE (Ripamonti & Abel 2004) cooling be included at high densities?
``H2OpticalDepthApproximation`` (external)
    Should the H2 cooling be attenuated (RA04)?
``H2FormationOnDust`` (external)
    Turns on H2 formation on dust grains and gas-grain heat transfer following Omukai (2000). Default: 0 (OFF)
``NumberOfDustTemperatureBins`` (external)
    Number of dust temperature bins for the dust cooling and H2 formation rates.  Default: 250
``DustTemperatureStart`` (external)
    Minimum dust temperature for dust rates.  Default: 1.0
``DustTemperatureEnd`` (external)
    Maximum dust temperature for dust rates.  Default: 1500
``OutputDustTemperature`` (external)
    Flag to write out the dust temperature field.  Default: 0
``PhotoelectricHeating`` (external)
    If set to be 1, the following parameter will be added uniformly
    to the gas without any shielding (Tasker & Bryan 2008). Default: 0
``PhotoelectricHeatingRate`` (external)
    This is the parameter used as Gamma_pe for uniform photoelectric heating.
    Default: 8.5e-26 erg s^-1 cm^-3

.. _cloudy_cooling:

Cloudy Cooling
^^^^^^^^^^^^^^

Cloudy cooling from Smith, Sigurdsson, & Abel (2008) interpolates
over tables of precomputed cooling data. Cloudy cooling is turned
on by setting ``MetalCooling`` to 3. ``RadiativeCooling`` must also be set
to 1. Depending on the cooling data used, it can be coupled with
``MultiSpecies`` = 1, 2, or 3 so that the metal-free cooling comes from
the ``MultiSpecies`` machinery and the Cloudy tables provide only the
metal cooling. Datasets range in dimension from 1 to 5. Dim 1:
interpolate over temperature. Dim 2: density and temperature. Dim
3: density, metallicity, and temperature. Dim 4: density,
metallicity, electron fraction, and temperature. Dim 5: density,
metallicity, electron fraction, spectral strength, and temperature.
See Smith, Sigurdsson, & Abel (2008) for more information on
creating Cloudy datasets.

``CloudyCoolingGridFile`` (external)
    A string specifying the path to the Cloudy cooling dataset.
``IncludeCloudyHeating`` (external)
    An integer (0 or 1) specifying whether the heating rates are to be
    included in the calculation of the cooling. Some Cloudy datasets
    are made with the intention that only the cooling rates are to be
    used. Default: 0 (off).
``CMBTemperatureFloor`` (external)
    An integer (0 or 1) specifying whether a temperature floor is
    created at the temperature of the cosmic microwave background
    (T\ :sub:`CMB`\  = 2.72 (1 + z) K). This is accomplished in the
    code by subtracting the cooling rate at T\ :sub:`CMB`\  such that
    Cooling = Cooling(T) - Cooling(T\ :sub:`CMB`\ ). Default: 1 (on).
``CloudyElectronFractionFactor`` (external)
    A float value to account for additional electrons contributed by
    metals. This is only used with Cloudy datasets with dimension
    greater than or equal to 4. The value of this factor is calculated
    as the sum of (A\ :sub:`i`\  \* i) over all elements i heavier than
    He, where A\ :sub:`i`\  is the solar number abundance relative to
    H. For the solar abundance pattern from the latest version of
    Cloudy, using all metals through Zn, this value is 9.153959e-3.
    Default: 9.153959e-3.

.. _grackle_pars:

The Grackle
^^^^^^^^^^^

The Grackle is an external chemistry and cooling library originally derived from 
Enzo's MultiSpecies chemistry and Cloudy cooling modules.  See :ref:`here <Grackle>` 
for a full description, including why you might use this over Enzo's internal 
chemistry and cooling.  For more information on Grackle parameter, see also the 
`Grackle documentation <https://grackle.readthedocs.org/>`_.  Note, some Grackle 
parameters have been mapped to Enzo parameters for simplicity.

``use_grackle`` (int)
    Flag to use the Grackle machinery (1 - on, 0 - off). Default: 0.

``with_radiative_cooling`` (int)
    Flag to include radiative cooling and actually update the thermal energy during the chemistry solver.  If off, the chemistry species will still be updated.  The most common reason to set this to off is to iterate the chemistry network to an equilibrium state (1 - on, 0 - off).  Default: 1.

``MultiSpecies`` (int) [mapped to Grackle parameter ``primordial_chemistry``]
    Flag to control which primordial chemistry network is used.  Default: 0.

    - 0: no chemistry network.  Radiative cooling for primordial species is solved by interpolating from lookup tables calculated with Cloudy.
    - 1: 6-species atomic H and He.  Active species: H, H\ :sup:`+`, He, He\ :sup:`+`, \ :sup:`++`, e\ :sup:`-`.
    - 2: 9-species network including atomic species above and species for molecular hydrogen formation.  This network includes formation from the H\ :sup:`-` and H\ :sub:`2`\ :sup:`+` channels, three-body formation (H+H+H and H+H+H\ :sub:`2`), H\ :sub:`2` rotational transitions, chemical heating, and collision-induced emission (optional).  Active species: above + H\ :sup:`-`, H\ :sub:`2`, H\ :sub:`2`\ :sup:`+`.
    - 3: 12-species network include all above plus HD rotation cooling.  Active species: above plus D, D\ :sup:`+`, HD.

``H2FormationOnDust`` (int) [mapped to Grackle parameter ``h2_on_dust``]
    See Enzo equivalent above.  Default: 0.

``MetalCooling`` (int) [mapped to Grackle parameter ``metal_cooling``]
    Flag to enable metal cooling using the Cloudy tables.  If enabled, the cooling table to be used must be specified with the ``grackle_data_file`` parameter (1 - on, 0 - off).  Default: 0.

``CMBTemperatureFloor`` (int) [mapped to Grackle parameter ``cmb_temperature_floor``]
    See Enzo equivalent above.  Default: 1.

``UVbackground`` (int)
    Flag to enable a UV background.  If enabled, the cooling table to be used must be specified with the ``grackle_data_file`` parameter (1 - on, 0 - off).  Default: 0.

``grackle_data_file`` (string)
    Path to the data file containing the metal cooling and UV background tables.  Default: "".

``Gamma`` (float)
    See Enzo equivalent above.  Default:  5/3.

``ThreeBodyRate`` (int) [mapped to Grackle parameter ``three_body_rate``]
    See Enzo equivalent above.  Default: 0.

``CIECooling`` (int) [mapped to Grackle parameter ``cie_cooling``]
    See Enzo equivalent above.  Default: 0.

``H2OpticalDepthApproximation`` (int) [mapped to Grackle parameter ``h2_optical_depth_approximation``]
    See Enzo equivalent above.  Default: 0.

``PhotoelectricHeating`` (int) [mapped to Grackle parameter ``photoelectric_heating``]
    See Enzo equivalent above.  Default: 0.

``PhotoelectricHeatingRate`` (float) [mapped to Grackle parameter ``photoelectric_heating_rate``]
    See Enzo equivalent above.  Default: 8.5e-26.

``Compton_xray_heating`` (int)
   Flag to enable Compton heating from an X-ray background following `Madau & Efstathiou (1999) <http://adsabs.harvard.edu/abs/1999ApJ...517L...9M>`_.  Default: 0.

``LWbackground_intensity`` (float)
    Intensity of a constant Lyman-Werner H\ :sub:`2` photo-dissociating radiation field in units of 10\ :sup:`-21` erg s\ :sup:`-1` cm\ :sup:`-2` Hz\ :sup:`-1` sr\ :sup:`-1`.  Default: 0.

``LWbackground_sawtooth_suppression`` (int)
    Flag to enable suppression of Lyman-Werner flux due to Lyman-series absorption (giving a sawtooth pattern), taken from `Haiman & Abel, & Rees (2000) <http://adsabs.harvard.edu/abs/2000ApJ...534...11H>`_.  Default: 0.
