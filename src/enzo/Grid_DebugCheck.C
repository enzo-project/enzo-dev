/***********************************************************************
/
/  GRID CLASS (CHECKS FOR NANS)
/
/  written by: Greg Bryan
/  date:       1999
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/
************************************************************************/
 
 
#include <stdlib.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "fortran.def"
#include "Grid.h"
 
// Turn the first on to check the gas and dark matter data for nan's.
// Turn the second on just to print the message (i.e. which routine
//     has called it.
 
#define DEBUG_CHECK_ON
#define TRACE_OFF
 
int grid::DebugCheck(char *message)
{
 
#ifdef TRACE_ON
 
  if (ProcessorNumber == MyProcessorNumber)
    fprintf(stderr, "P(%"ISYM"): %s\n", MyProcessorNumber, message);
 
#endif /* TRACE_ON */
 
#ifdef DEBUG_CHECK_ON
 
  // Return if this doesn't concern us
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  // Set this to zero so the compiler doesn't optimize everything away
 
  int ThisIsZero = GridStartIndex[0] - NumberOfGhostZones, size = 1,
      dim, k1, k2;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];



  int failflag = 0; 
  for(k1=0; k1<NumberOfBaryonFields; ++k1) {
    int fi = this->getField(k1);
    if (  fi==TotalEnergy ) {
      for (k2=0; k2<size; ++k2) {
        double ke=0.0;
        for(int i=0; i<NumberOfBaryonFields; ++i) {
          int ii = this->getField(i);
          if( ii == Velocity1 || ii==Velocity2 || ii==Velocity3) {
            ke += 0.5 *BaryonField[i][k2]*BaryonField[i][k2];
          }
        }
        double de = 0.0;
        if( this->getField(0) == Density )
          de = BaryonField[ 0 ][k2]; 
        double inte = 0.0;
        if( this->getField(2) == InternalEnergy)
          inte = BaryonField[2][k2];
        double tote = 0.0;
        if( this->getField(1) == TotalEnergy)
          tote = BaryonField[1][k2];
        //if( tote > EnergyCeiling || ke>EnergyCeiling || inte > EnergyCeiling ||  (de < 25.0 && de>1.0)) {
        if( inte > EnergyCeiling ||  (de < 15.0 && de>1.0)) {
            fprintf(stderr, "DebugCheck[%s](Proc %"ISYM"): %"ISYM", %"GSYM" %"GSYM" %"GSYM" %"GSYM"\n", message,
            	MyProcessorNumber, k2, inte, tote, ke, float(CellWidth[0][0]));
	  for(int i=0; i<NumberOfBaryonFields; ++i) {
            int ii = this->getField(i);
	    fprintf(stderr, "DebugCheckDump %"ISYM" %"ISYM" %"GSYM" \n",i,ii,BaryonField[i][k2]);
          }
          failflag = 1;
        }
      }
    } 
  }
  if(failflag==1) {
    return FAIL;
  }

//  for (k1 = 0; k1 < NumberOfBaryonFields; k1++)
//    for (k2 = 0; k2 < size; k2++)
//      if (BaryonField[k1][k2+ThisIsZero] != BaryonField[k1][k2]) {
//	fprintf(stderr, "DebugCheck[%s](Proc %"ISYM"): gas %"ISYM" %"ISYM" %"GSYM"\n", message,
//		MyProcessorNumber, k1, k2, BaryonField[k1][k2]);
//	exit(EXIT_FAILURE);
//      }
// 
//  for (k1 = 0; k1 < GridRank; k1++)
//    for (k2 = 0; k2 < NumberOfParticles; k2++)
//      if (ParticlePosition[k1][k2] != ParticlePosition[k1][k2+ThisIsZero] ||
//	  ParticleVelocity[k1][k2] != ParticleVelocity[k1][k2+ThisIsZero]  ) {
//	fprintf(stderr, "DebugCheck[%s](Proc %"ISYM"): dm %"ISYM" (%"ISYM"/%"ISYM") %"GSYM" %"GSYM"\n",
//		message, MyProcessorNumber, k1, k2, NumberOfParticles,
//		ParticlePosition[k1][k2], ParticleVelocity[k1][k2]);
//	exit(EXIT_FAILURE);
//      }
 
#if 0		
  if (NumberOfBaryonFields > 0)
    for (k1 = 0; k1 < 2+DualEnergyFormalism; k1++)
      for (k2 = 0; k2 < size; k2++)
	if (BaryonField[k1][k2+ThisIsZero] <= 0) {
	  fprintf(stderr, "DebugCheck[%s](Proc %"ISYM"): <0 %"ISYM" %"ISYM" %"GSYM"\n", message,
		  MyProcessorNumber, k1, k2, BaryonField[k1][k2]);
	  exit(EXIT_FAILURE);
	}
#endif
 
#endif /* DEBUG_CHECK_ON */
		
  return SUCCESS;
}
