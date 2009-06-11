/***********************************************************************
/
/  GRID CLASS (HANDLE CALLING AND SOLVING COOLING/CHEMISTRY)
/
/  written by: Matthew Turk
/  date:       June, 2009
/  modified1:
/
/  PURPOSE: Move logic for chemistry/cooling module selection here
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/
 
int grid::MultiSpeciesHandler()
{
  if ((!MultiSpecies) && (!RadiativeCooling)) return SUCCESS; 

  if (MultiSpecies && RadiativeCooling) {
	  Grids[grid1]->GridData->SolveRateAndCoolEquations();
  } else {
    if (MultiSpecies)
      Grids[grid1]->GridData->SolveRateEquations();
    if (RadiativeCooling)
      Grids[grid1]->GridData->SolveRadiativeCooling();
  }

  return SUCCESS;
}
