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
  int i = 0;

  if ((fptr = fopen(RefineRegionFile, "r")) == NULL) {
    fprintf(stderr, "Error opening refine region file %s.\n", RefineRegionFile);
    return FAIL;
  }

  char line[MAX_LINE_LENGTH];
  int nret;
  EvolveRefineRegionNtimes=0;
  while ((fgets(line, MAX_LINE_LENGTH, fptr) != NULL)){
    nret = sscanf(line, "%"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM,
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

  return SUCCESS;

}
