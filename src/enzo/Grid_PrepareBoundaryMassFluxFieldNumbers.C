#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

int FindField(int f, int farray[], int n);


int grid::PrepareBoundaryMassFluxFieldNumbers(){

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  if (!StoreDomainBoundaryMassFlux) return SUCCESS;

  if (BoundaryMassFluxFieldNumbers[0] >= 0) return SUCCESS; // we've already initialized

  /* obtain baryon field indexes */
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, CRNum, B1Num, B2Num, B3Num;
  this->DebugCheck("StarParticleHandler");
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                       Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }

  /* identify species fields if they exist */
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
  if ( MultiSpecies ){
    IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
                          HMNum, H2INum, H2IINum, DINum, DIINum, HDINum);
  }

  /* get metallicity tracer field number */
  int MetalNum;
  MetalNum   = FindField(Metallicity, this->FieldType, this->NumberOfBaryonFields);

  /* pre-compute field numbers so we only need to do this once */
  int count = 0;
  BoundaryMassFluxFieldNumbers[count++] = DensNum;

  /* H, He and ionization states */
  if (MultiSpecies > 0){
    BoundaryMassFluxFieldNumbers[count++] = HINum;
    BoundaryMassFluxFieldNumbers[count++] = HIINum;
    BoundaryMassFluxFieldNumbers[count++] = HeINum;
    BoundaryMassFluxFieldNumbers[count++] = HeIINum;
    BoundaryMassFluxFieldNumbers[count++] = HeIIINum;
  }

  /* H-, H2, H2- */
  if (MultiSpecies > 1){
    BoundaryMassFluxFieldNumbers[count++] = HMNum;
    BoundaryMassFluxFieldNumbers[count++] = H2INum;
    BoundaryMassFluxFieldNumbers[count++] = H2IINum;
  }

  if (MultiSpecies > 2){
    BoundaryMassFluxFieldNumbers[count++] = DINum;
    BoundaryMassFluxFieldNumbers[count++] = DIINum;
    BoundaryMassFluxFieldNumbers[count++] = HDINum;
  }

  /* metal density */
  if (MultiMetals > 0){
    BoundaryMassFluxFieldNumbers[count++] = MetalNum;
  }

  /* stellar yields tracers */
  if (MultiMetals > 1 && STARMAKE_METHOD(INDIVIDUAL_STAR)){
    for (int yield = 0; yield < StellarYieldsNumberOfSpecies; yield++){
      int anum = StellarYieldsAtomicNumbers[yield];
      if (anum == 1 || anum == 2) continue; // handled above

      int field_num;
      this->IdentifyChemicalTracerSpeciesFieldsByNumber(field_num, anum);
      BoundaryMassFluxFieldNumbers[count++] = field_num;
    }
  }

  return SUCCESS;
}
