/************************************************************************
/
/   RETURN VALUES AFTER READING (BLACK-BODY) SPECTRUM TABLE
/
/   written by: Ji-hoon Kim
/   date: February, 2010
/   modified1:
/
/   PURPOSE: For a given black body spectrum read in ReadSpectrumTable,
/            return the fraction of absorbed photons or the mean energy
/            of the spectrum at this column density
/
/   INPUTS:  type: type of absorber
/            ColumnDensity: column density so far
/            dColumnDensity: column density added in this cell
/
/   RETURNS: for mode = 0,1,2: the fraction of absorbed photons 
/                              by the absorber species HI, HeI, HeII
/            for mode = 3    : the mean energy of the spectrum
/
/************************************************************************/

#ifdef TRANSFER

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "ExternalBoundary.h"
#include "Fluxes.h"
#include "GridList.h"
#include "Grid.h"
#include "CosmologyParameters.h"

FLOAT FindCrossSection(int type, float energy);

float ReturnValuesFromSpectrumTable(float ColumnDensity, float dColumnDensity, 
				    int mode)
{

  if (mode < 0 || mode > 3) {
    ENZO_FAIL("ReturnValuesFromSpectrumTable: mode unrecognized\n");
  }
    
  int index_in, index_out;
  float frac_in = 0.0, frac_out = 0.0, photon_fraction = 0.0, mean_energy;
  float change, tau;
  float logC_start, logC_end, logC_step, logC_in, logC_out;
  FLOAT pseudo_CrossSection;

  int nbins  = RadiativeTransferSpectrumTable.NumberOfColumnDensityBins;
  logC_start = log(RadiativeTransferSpectrumTable.columndensity_table[0]);
  logC_end   = log(RadiativeTransferSpectrumTable.columndensity_table[nbins-1]);
  logC_step  = (logC_end - logC_start) / (nbins-1);

  /* calculate indices for ColumnDensity and ColumnDensity+dColumnDensity */

  logC_in   = log(ColumnDensity);
  logC_out  = log(ColumnDensity + dColumnDensity);

  index_in  = min(nbins-1, max(1, int((logC_in  - logC_start)/logC_step)+1));
  index_out = min(nbins-1, max(1, int((logC_out - logC_start)/logC_step)+1));

  /* find mean energy */
  
  mean_energy = RadiativeTransferSpectrumTable.meanenergy_table[index_in];

  /* return values */

  if (mode == 3) {

    return mean_energy;

  } else if (mode >= 0 && mode <= 2) {

    /* find photons left for ColumnDensity and ColumnDensity+dColumnDensity */

    frac_in  = RadiativeTransferSpectrumTable.fractionphotons_table[mode][index_in];
    frac_out = RadiativeTransferSpectrumTable.fractionphotons_table[mode][index_out];

    /* find fraction of photons absorbed by species w.r.t. the incoming number of photons;
       using the same logic for tau in Grid_WalkPhotonPackage, if "change" is smaller 
       than float precision, try linear interpolation to get absorbed photon fraction */

    /* tau = dColumnDensity * FindCrossSection(0, mean_energy), but we cannot use this!
       instead, we find the proportionality constant near ColumnDensity using
       frac_in = exp(-pseudo_CrossSection * ColumnDensity) */
    
    pseudo_CrossSection = -log(frac_in) / 
      max(ColumnDensity, RadiativeTransferSpectrumTable.columndensity_table[0]);

    /* expf(-tau) = frac_out / frac_in, below is to avoid cases such as frac_in = 0.0 */

    change = frac_out / frac_in;
    tau = (isnan(change)) ? 
      dColumnDensity * pseudo_CrossSection : -log(frac_out / frac_in);   

    /* return photon_fraction */

    if (tau > 2.e1) 
      photon_fraction = (1.0+BFLOAT_EPSILON);
    else if (tau > 1.e-4) 
      photon_fraction = min(1 - frac_out / frac_in, 1.0);
    else
      photon_fraction = min(dColumnDensity * pseudo_CrossSection, 1.0);  

//    fprintf(stderr, "RVFST: id_in = %d, id_out = %d, f_in =%f, f_out = %f, tau = %f, photon_f = %f\n", 
//	    index_in, index_out, frac_in, frac_out, tau, photon_fraction); 

    return photon_fraction;

  } else {
    ENZO_FAIL("Unreconized Mode!\n");
  }





#ifdef TEST_WITH_FORMULA  // this was for test

  mean_energy = 2000;  

  if (mode == 3) {

    return mean_energy;

  } else if (mode >= 0 && mode <= 2) {

    // optical depth of ray segment (by HI, note that dColumnDensity already has LengthUnits in it)

    float tau = dColumnDensity * FindCrossSection(0, mean_energy);
  
    // calculate the fraction of photons absorbed at this column density
    if (tau > 2.e1) 
      photon_fraction = (1.0+BFLOAT_EPSILON);
    else if (tau > 1.e-4) 
      photon_fraction = min((1-expf(-tau)), 1.0);
    else
      photon_fraction = min(tau, 1.0);  

    return photon_fraction;

  } else {

    ENZO_FAIL("ReturnValuesFromSpectrumTable: mode unrecognized\n");


  }

#endif


}

#endif
