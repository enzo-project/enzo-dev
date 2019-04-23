/*------------------------------------------------------------------------
  READ EVOLVING REFINE REGION FILE
  By John Wise

  File format: 
    (time or redshift) x_left y_left z_left x_right y_right z_right

  History:
     03 May 2005 : JHW -- Created
------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

int ReadEvolveRefineFile(void)
{

  FILE *fptr;

  if((RefineRegionTimeType == 0) || (RefineRegionTimeType == 1)){
    
    if ((fptr = fopen(RefineRegionFile, "r")) == NULL) {
      fprintf(stderr, "Error opening refine region file %s.\n", RefineRegionFile);
      return FAIL;
    }

    char line[MAX_LINE_LENGTH];
    int nret, i=0;
    EvolveRefineRegionNtimes=0;
    while ((fgets(line, MAX_LINE_LENGTH, fptr) != NULL)){
      nret = sscanf(line, "%"FSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM,
		    &(EvolveRefineRegionTime[i]),
		    &(EvolveRefineRegionLeftEdge[i][0]),
		    &(EvolveRefineRegionLeftEdge[i][1]),
		    &(EvolveRefineRegionLeftEdge[i][2]),
		    &(EvolveRefineRegionRightEdge[i][0]),
		    &(EvolveRefineRegionRightEdge[i][1]),
		    &(EvolveRefineRegionRightEdge[i][2])); 
      if( nret != 7 ){
	fprintf(stderr,"WARNING: ReadEvolveRefineFile cannot interpret line %s",line);
	continue;
      }
      i++;
      EvolveRefineRegionNtimes++;
    }
    
    fclose(fptr);

    // Error check - are there too many input times?
    if(EvolveRefineRegionNtimes > MAX_REFINE_REGIONS){
      fprintf(stderr, "Too many EvolveRefineRegion times in your file!\nIncrease MAX_REFINE_REGIONS in macros_and_parameters.h!\n");
      return FAIL;
    }
    
  }

  /* READ IN MustRefineRegion information */
  if((MustRefineRegionTimeType == 0) || (MustRefineRegionTimeType == 1)){

    fprintf(stderr,"ReadEvolveRefineFile: I am reading in a MustRefineRegion time evolution file!\n");
    
    if ((fptr = fopen(MustRefineRegionFile, "r")) == NULL) {
      fprintf(stderr, "Error opening MustRefine region file %s.\n", MustRefineRegionFile);
      return FAIL;
    }

    char line[MAX_LINE_LENGTH];
    int nret, i=0;
    EvolveMustRefineRegionNtimes=0;
    while ((fgets(line, MAX_LINE_LENGTH, fptr) != NULL)){
      nret = sscanf(line, "%"FSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"ISYM,
		    &(EvolveMustRefineRegionTime[i]),
		    &(EvolveMustRefineRegionLeftEdge[i][0]),
		    &(EvolveMustRefineRegionLeftEdge[i][1]),
		    &(EvolveMustRefineRegionLeftEdge[i][2]),
		    &(EvolveMustRefineRegionRightEdge[i][0]),
		    &(EvolveMustRefineRegionRightEdge[i][1]),
		    &(EvolveMustRefineRegionRightEdge[i][2]),
      		    &(EvolveMustRefineRegionMinLevel[i]));
      if(debug && MyProcessorNumber == ROOT_PROCESSOR){
         fprintf(stderr,"Here is the line (MustRefineRegion): %s \n",line);
         fprintf(stderr,". . . and here is the value (MustRefineRegion): %i \n",EvolveMustRefineRegionMinLevel[i]);
         } 
      if( nret != 8 ){
	fprintf(stderr,"WARNING: ReadEvolveRefineFile (MustRefineRegion) cannot interpret line %s",line);
	continue;
      }
      i++;
      EvolveMustRefineRegionNtimes++;
    }

    fclose(fptr);

    // Error check - are there too many input times?
    if(EvolveMustRefineRegionNtimes > MAX_REFINE_REGIONS){
      fprintf(stderr, "Too many EvolveMustRefineRegion times in your file!\nIncrease MAX_REFINE_REGIONS in macros_and_parameters.h!\n");
      return FAIL;
    }

    
    if(debug && MyProcessorNumber == ROOT_PROCESSOR){

      printf("ReadEvolveRefineFile: I have a MustRefineRegion with TimeType %"ISYM" \n",
	     MustRefineRegionTimeType);
      
      printf("ReadEvolveRefineFile: And here is what I think my times, edges, and minimum levels are:\n");

      for(int i=0; i<EvolveMustRefineRegionNtimes; i++){

	printf("ReadEvolveRefineFile (MustRefineRegion): %"FSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"ISYM"\n",
	       EvolveMustRefineRegionTime[i],
	       EvolveMustRefineRegionLeftEdge[i][0],
	       EvolveMustRefineRegionLeftEdge[i][1],
	       EvolveMustRefineRegionLeftEdge[i][2],
	       EvolveMustRefineRegionRightEdge[i][0],
	       EvolveMustRefineRegionRightEdge[i][1],
	       EvolveMustRefineRegionRightEdge[i][2],
	       EvolveMustRefineRegionMinLevel[i]);

      } // for loop

      fflush(stdout);

    } // if(debug && MyProcessorNumber == ROOT_PROCESSOR)

  } // if((MustRefineRegionTimeType == 0) || (MustRefineRegionTimeType == 1))



  /* READ IN CoolingRefineRegion information.  Note that this requires a file that is 
     EXACTLY like the MustRefineRegion file, which includes a level as the last entry in 
     each row of the file.  This is NOT USED but must be there.  */
  if((CoolingRefineRegionTimeType == 0) || (CoolingRefineRegionTimeType == 1)){

    fprintf(stderr,"ReadEvolveRefineFile: I am reading in a CoolingRefineRegion time evolution file!\n");
    
    if ((fptr = fopen(CoolingRefineRegionFile, "r")) == NULL) {
      fprintf(stderr, "Error opening CoolingRefine region file %s.\n", CoolingRefineRegionFile);
      return FAIL;
    }

    char line[MAX_LINE_LENGTH];
    int nret, i=0, dummy;
    EvolveCoolingRefineRegionNtimes=0;
    while ((fgets(line, MAX_LINE_LENGTH, fptr) != NULL)){
      nret = sscanf(line, "%"FSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"ISYM,
		    &(EvolveCoolingRefineRegionTime[i]),
		    &(EvolveCoolingRefineRegionLeftEdge[i][0]),
		    &(EvolveCoolingRefineRegionLeftEdge[i][1]),
		    &(EvolveCoolingRefineRegionLeftEdge[i][2]),
		    &(EvolveCoolingRefineRegionRightEdge[i][0]),
		    &(EvolveCoolingRefineRegionRightEdge[i][1]),
		    &(EvolveCoolingRefineRegionRightEdge[i][2]),
      		    &dummy);
      if(debug && MyProcessorNumber == ROOT_PROCESSOR){
         fprintf(stderr,"Here is the line (CoolingRefineRegion): %s \n",line);
         fprintf(stderr,". . . and here is the value (CoolingRefineRegion): %i \n",dummy);
         } 
      if( nret != 8 ){
	fprintf(stderr,"WARNING: ReadEvolveRefineFile cannot interpret line %s",line);
	continue;
      }
      i++;
      EvolveCoolingRefineRegionNtimes++;
    }

    fclose(fptr);

    // Error check - are there too many input times?
    if(EvolveCoolingRefineRegionNtimes > MAX_REFINE_REGIONS){
      fprintf(stderr, "Too many EvolveCoolingRefineRegion times in your file!\nIncrease MAX_REFINE_REGIONS in macros_and_parameters.h!\n");
      return FAIL;
    }

    
    if(debug && MyProcessorNumber == ROOT_PROCESSOR){

      printf("ReadEvolveRefineFile: I have a CoolingRefineRegion with TimeType %"ISYM" \n",
	     CoolingRefineRegionTimeType);
      
      printf("ReadEvolveRefineFile: And here is what I think my times, edges, and minimum levels are:\n");

      for(int i=0; i<EvolveCoolingRefineRegionNtimes; i++){

	printf("ReadEvolveRefineFile (CoolingRefineRegion): %"FSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM" %"PSYM"\n",
	       EvolveCoolingRefineRegionTime[i],
	       EvolveCoolingRefineRegionLeftEdge[i][0],
	       EvolveCoolingRefineRegionLeftEdge[i][1],
	       EvolveCoolingRefineRegionLeftEdge[i][2],
	       EvolveCoolingRefineRegionRightEdge[i][0],
	       EvolveCoolingRefineRegionRightEdge[i][1],
	       EvolveCoolingRefineRegionRightEdge[i][2]);

      } // for loop

      fflush(stdout);

    } // if(debug && MyProcessorNumber == ROOT_PROCESSOR)

  } // if((CoolingRefineRegionTimeType == 0) || (CoolingRefineRegionTimeType == 1))

  
  
  
  return SUCCESS;

}
