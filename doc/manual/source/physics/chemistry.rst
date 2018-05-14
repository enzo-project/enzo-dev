.. _cooling:

Gas Chemistry, Cooling, and Heating
==========================

There are two ways that Enzo can perform radiative cooling.
Natively, Enzo provides a number of different methods ranging from
simple tabulated, analytical approximations to very sophisticated 
non-equilibrium primoridal chemistry. These are described below. Alternatively, Enzo supports the external `Grackle library <http://grackle.readthedocs.io>`_ for tracking chemistry and performing radiative cooling. 
See :ref:`here <Grackle>` for a full description and reasons to use Grackle.
Both Grackle and the native Enzo methods support metal cooling via Cloudy tables.

All of these methods require the parameter ``RadiativeCooling`` be set to 1.
Other parameters are required for using the various native methods, including
Grackle, which are described below.
For all relevant parameters, please also see :ref:`cooling_parameters`.

MultiSpecies = 0: Sarazin & White
---------------------------------
*Source: solve_cool.F, cool1d.F*

``RadiativeCooling`` = 1

``MultiSpecies`` = 0

This method uses an analytical approximation from `Sarazin & White
(1987, ApJ, 320, 32) <http://adsabs.harvard.edu/abs/1987ApJ...320...32S>`_
for a fully ionized gas with metallicity of 0.5
solar.  This cooling curve is valid over the temperature range from
10,000 K to 10\ :sup:`9`\  K.  Since this assumes a fully ionized gas, the
cooling rate is effectively zero below 10,000 K.

*Note: In order use this cooling method, you must copy the file,
cool_rates.in, from the input directory into your simulation directory.*

MultiSpecies = 1, 2, or 3: Primordial Chemistry and Cooling
-----------------------------------------------------------
*Source: multi_cool.F, cool1d_multi.F*

This method follows the nonequilibrium evolution of primordial
(metal-free) gas.  The chemical rate equations are solved using a
semi-implicit backward differencing scheme described by 
`Abel et al. (1997, New Astronomy, 181) 
<http://adsabs.harvard.edu/abs/1997NewA....2..181A>`_
and `Anninos et al. (1997, New Astronomy, 209)
<http://adsabs.harvard.edu/abs/1997NewA....2..209A>`_.
Heating and cooling processes include atomic line
excitation, recombination, collisional excitation, free-free
transitions, Compton scattering of the cosmic microwave background and
photoionization from a variety of metagalactic UV backgrounds.  For 
``MultiSpecies`` > 1, molecular cooling is also included and UV
backgrounds that include photodissociation may also be used.
Numerous chemistry and cooling rates have been added or updated.  For
the exact reference for any given rate, users are encouraged to
consult *calc_rates.F*.

#. Atomic

   ``RadiativeCooling`` = 1

   ``MultiSpecies`` = 1

   Only atomic species, H, H\ :sup:`+`\, He, He\ :sup:`+`\, He\
   :sup:`++`\, and e\ :sup:`-`\  are followed.  Since 
   molecular species are not treated, the cooling is effectively zero for
   temperatures below roughly 10,000 K.

#. Molecular Hydrogen

   ``RadiativeCooling`` = 1

   ``MultiSpecies`` = 2

   Along with the six species above, H\ :sub:`2`\, H\
   :sub:`2`:sup:`+`\, and H\ :sup:`-`\  are also followed.
   In addition to the rates described in `Abel et al. (1997, New Astronomy, 181) 
   <http://adsabs.harvard.edu/abs/1997NewA....2..181A>`_
   and `Anninos et al. (1997, New Astronomy, 209)
   <http://adsabs.harvard.edu/abs/1997NewA....2..209A>`_, 
   H2 formation via three-body reactions as described by
   `Abel, Bryan, and Norman (2002, Science, 295, 93)
   <http://adsabs.harvard.edu/abs/2002Sci...295...93A>`_ is also included.
   This method is valid in the temperature range of 1 K to 10\
   :sup:`8`\  K and up to number densities of roughly 10\ :sup:`9`\  cm\ :sup:`-3`\.
   Additionally, three-body heating (4.48eV per molecule formed or dissociated)
   is added as appropriate.

#. Deuterium

   ``RadiativeCooling`` = 1

   ``MultiSpecies`` = 3

   In addition to the nine species solved with ``MultiSpecies`` = 2,
   D, D\ :sup:`+`\, and HD are also followed.  The range of validity
   is the same as for ``MultiSpecies`` = 2.

Metal Cooling
-------------

Three distinct methods to calculate the cooling from elements heavier
than He exist.  These are selected by setting the ``MetalCooling``
parameter to 1, 2, or 3.

#. John Wise's metal cooling.

   ``RadiativeCooling`` = 1

   ``MetalCooling`` = 1

#. `Cen et al (1995) <http://adsabs.harvard.edu/abs/1995ApJ...451..436C>`_ 
   cooling. This uses output from a Raymond-Smith
   code to determine cooling rates from T > 10\ :sup:`4`\ K.  No ionizing
   background is used in computing cooling rates.  This method has
   not been extensively tested in the context of Enzo.

   ``RadiativeCooling`` = 1

   ``MetalCooling`` = 2

#. Cloudy cooling.

   *Source: cool1d_cloudy.F*

   ``RadiativeCooling`` = 1

   ``MetalCooling`` = 3

   ``MultiSpeces`` = 1, 2, or 3

   Cloudy cooling operates in conjunction with the primordial
   chemistry and cooling from ``MultiSpecies`` set to 1, 2, or 3.
   As described in `Smith, Sigurdsson, & Abel (2008)
   <http://adsabs.harvard.edu/abs/2008MNRAS.385.1443S>`_, Cloudy cooling
   interpolates over tables of precomputed cooling data using the
   Cloudy photoionization software (`Ferland et al. 1998, PASP, 110,
   761
   <http://adsabs.harvard.edu/abs/1998PASP..110..761F>`_, 
   `<http://nublado.org>`_).  The cooling datasets can be from one to
   five dimensional.  The range of validity will depends on the
   dataset used.

   #. Temperature
   #. Density and temperature.
   #. Density, metallicity, and temperature.
   #. Density, metallicity, electron fraction, and temperature.
   #. Density, metallicity, electron fraction, redshift of UV
      background, and temperature.

   See :ref:`cloudy_cooling` for additional parameters that control
   the behavior of the Cloudy cooling.  For more information on
   obtaining or creating Cloudy cooling datasets, contact Britton
   Smith (brittonsmith@gmail.com).

Using Grackle
-------------

To use the Grackle library for radiative cooling:

    ``use_grackle`` = 1

    ``with_radiative_cooling`` = 1 

The ``MultiSpecies`` and ``MetalCooling`` Enzo parameters are mapped to their
respective Grackle parameters, with ``MultiSpecies`` behaving the same way as
described above. 

When using Grackle, ``MetalCooling`` can only be 0 or 1 for off
and on respectively. If ``MetalCooling`` is on, a Cloudy cooling table must be
specified using the ``grackle_data_file`` parameter.

See :ref:`grackle_pars` for more parameters that can be used with Grackle.

UV Meta-galactic Backgrounds
----------------------------
*Source: RadiationFieldCalculateRates.C*

Enzo natively includes a variety of spatially uniform photoionizing 
and photodissociating backgrounds, mainly by setting the parameter
``RadiationFieldType``.  These radiation backgrounds are redshift
dependent and work by setting the photoionization and photoheating
coeffiecients for H, He, and He\ :sup:`+`\.  See
:ref:`radiation_backgrounds` for the additional parameters that
control the UV backgrounds.

**If using Grackle**, the UV background is instead enabled via the ``UVbackground``
parameter, and specified with the ``grackle_data_file`` parameter.
