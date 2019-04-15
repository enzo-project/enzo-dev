/***********************************************************************
/
/  GALAXY PARTICLE TYPE
/
/  written by: Stephen Skory
/  date:       August, 2012
/
/  PURPOSE:
/
************************************************************************/

#include "ActiveParticle_GalaxyParticle.h"

/* We need to make sure that we can operate on the grid, so this dance is
 * necessary to make sure that grid is 'friend' to this particle type. */

class GalaxyParticleGrid : private grid {
  friend class ActiveParticleType_GalaxyParticle;
};

FLOAT calc_dist2(FLOAT x1, FLOAT y1, FLOAT z1,
    FLOAT x2, FLOAT y2, FLOAT z2, FLOAT period[])
{
    int dim;
    float part_dist2, xdist, ydist, zdist;
    // Periodicity
    xdist = fabs(x1 - x2);
    xdist = min(xdist, period[0] - xdist);
    ydist = fabs(y1 - y2);
    ydist = min(ydist, period[1] - ydist);
    zdist = fabs(z1 - z2);
    zdist = min(zdist, period[2] - zdist);
    part_dist2 = xdist * xdist + ydist * ydist + zdist * zdist;
    return part_dist2;
}

/* Note that we only refer to GalaxyParticleGrid here. 
 * Given a grid object, we static cast to get this:
 *
 *    GalaxyParticleGrid *thisgrid =
 *      static_cast<GalaxyParticleGrid *>(thisgrid_orig); */

int ActiveParticleType_GalaxyParticle::InitializeParticleType(void)
{

  // Need to turn on particle mass flagging if it isn't already turned on.

  bool TurnOnParticleMassRefinement = true;
  int method;
  for (method = 0; method < MAX_FLAGGING_METHODS; method++)
    if (CellFlaggingMethod[method] == 8 || CellFlaggingMethod[method] == 4) {
      TurnOnParticleMassRefinement = false;
      break;
    }
	
  if (TurnOnParticleMassRefinement) {
    method = 0;
    while(CellFlaggingMethod[method] != INT_UNDEFINED)
      method++;
    CellFlaggingMethod[method] = 4;
  }

  /* Add on the Particle Array Handlers */
  typedef ActiveParticleType_GalaxyParticle ap;
  AttributeVector &ah = ap::AttributeHandlers;
  ActiveParticleType::SetupBaseParticleAttributes(ah);

  ah.push_back(new Handler<ap, float, &ap::Radius>("Radius"));
  ah.push_back(new Handler<ap, int, &ap::initialized>("initialized"));

  return SUCCESS;
}

int ActiveParticleType_GalaxyParticle::EvaluateFormation
(grid *thisgrid_orig, ActiveParticleFormationData &data)
{
  // We actually allow insertion where there isn't maximum refinement,
  // although it's probably the case already, because we care about halos
  // not gas properties.
  //if (data.level != MaximumRefinementLevel) {
  //  return SUCCESS;
  //}
  //fprintf(stderr, "Checking formation of galaxy particles.\n");
  // Let's read in some galaxy particle data. This will be replaced by
  // the Python interface someday.
  FILE *fptr;
  char line[MAX_LINE_LENGTH];
  char fname[20];
  sprintf(fname, "gp%04"ISYM".txt", LevelCycleCount[0]);

  if ((fptr = fopen(fname, "r")) == NULL) {
    //fprintf(stderr, "galaxy particle text file not found for this cycle.\n");
    return SUCCESS;
  }

  float x, y, z, dens, radius, dist, temp, vx, vy, vz;
  float xx, yy, zz;
  int i, j, k, index, mark;

  GalaxyParticleGrid *thisGrid =
    static_cast<GalaxyParticleGrid *>(thisgrid_orig);


  int GridDimension[3] = {thisGrid->GridDimension[0],
                          thisGrid->GridDimension[1],
                          thisGrid->GridDimension[2]};

  float *density = thisGrid->BaryonField[data.DensNum];
  float *velx = thisGrid->BaryonField[data.Vel1Num];
  float *vely = thisGrid->BaryonField[data.Vel2Num];
  float *velz = thisGrid->BaryonField[data.Vel3Num];
  //float *tvel;

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
    sscanf(line, "%"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM,
      &x, &y, &z, &vx, &vy, &vz, &dens, &radius);

    // If this particle is outside this grid,
    // we don't do anything.
    if (x < thisGrid->GridLeftEdge[0] || x >= thisGrid->GridRightEdge[0] ||
        y < thisGrid->GridLeftEdge[1] || y >= thisGrid->GridRightEdge[1] ||
        z < thisGrid->GridLeftEdge[2] || z >= thisGrid->GridRightEdge[2]) {
      continue;
    }

    fprintf(stderr,"%d %d inserting particle %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM"\n",
        MyProcessorNumber, data.GridID, x, y, z, dens, radius);

	// If no more room for particles, throw an ENZO_FAIL
	if (data.NumberOfNewParticles >= data.MaxNumberOfNewParticles)
	  return FAIL;

    ActiveParticleType_GalaxyParticle *np = new ActiveParticleType_GalaxyParticle();
    data.NumberOfNewParticles++;
    data.NewParticles.insert(*np);
	
	np->type = np->GetEnabledParticleID();
	np->BirthTime = thisGrid->ReturnTime();
	
	np->level = data.level;
	np->GridID = data.GridID;
	np->CurrentGrid = thisGrid;
	
	np->pos[0] = x;
	np->pos[1] = y;
	np->pos[2] = z;
    
    // This particle is born with zero mass (density), we'll fix this
    // when we apply 'feedback' in EvaluateFeedback.
	np->Mass = 0.;
	// For now, give it the velocity of the cell it lives in, but in the
	// future we will either want this calculated in concert with the mass
	// above, or have it be supplied by the halo finder.
    // Find the closest cell to x, y, z.
//     dist = huge_number;
//     for (k = thisGrid->GridStartIndex[2]; k <= thisGrid->GridEndIndex[2]; k++) {
//       for (j = thisGrid->GridStartIndex[1]; j <= thisGrid->GridEndIndex[1]; j++) {
//         for (i = thisGrid->GridStartIndex[0]; i <= thisGrid->GridEndIndex[0]; i++) {
//           index = GRIDINDEX_NOGHOST(i, j, k);
//    	      xx = thisGrid->CellLeftEdge[0][i] + 0.5*thisGrid->CellWidth[0][i];
//           yy = thisGrid->CellLeftEdge[1][j] + 0.5*thisGrid->CellWidth[1][j];
//           zz = thisGrid->CellLeftEdge[2][k] + 0.5*thisGrid->CellWidth[2][k];
//           temp = (xx - x) * (xx - x) + (yy - y) * (yy - y) + (zz - z) * (zz - z);
//           if (temp < dist) {
//             dist = temp;
//              mark = index;
//           }
//         }
//       }
//     }
// 	tvel = thisGrid->AveragedVelocityAtCell(mark ,data.DensNum,data.Vel1Num);
	
	np->vel[0] = vx;
	np->vel[1] = vy;
	np->vel[2] = vz;
	
	np->Radius = radius;
	np->Metallicity = 0.0;
	np->initialized = 0;
	
	// Clean up
	//delete [] tvel;
  } // end while fgets...

  return SUCCESS;
}

int ActiveParticleType_GalaxyParticle::EvaluateFeedback
(grid *thisgrid_orig, ActiveParticleFormationData &data)
{
  // Feedback takes place in grid::ApplyGalaxyParticleFeedback.

  return SUCCESS;
}

void ActiveParticleType_GalaxyParticle::DescribeSupplementalData
(ActiveParticleFormationDataFlags &flags)
{
  flags.DarkMatterDensity = true;
  flags.Temperature = true;
  flags.MetalField = true;
}

int ActiveParticleType_GalaxyParticle::SetFlaggingField
(LevelHierarchyEntry *LevelArray[], int level,
							int TopGridDims[], int GalaxyParticleID)
{
  /* Generate a list of all galaxy particles in the simulation box */
  int i, nParticles;
  FLOAT *pos = NULL, dx=0, rad;
  ActiveParticleList<ActiveParticleType> GalaxyParticleList;
  LevelHierarchyEntry *Temp = NULL;
  
  ActiveParticleFindAll(LevelArray, &nParticles, GalaxyParticleID, 
      GalaxyParticleList);
  
  /* Calculate CellWidth on maximum refinement level */
  
  // this will fail for noncubic boxes or simulations with MinimimMassForRefinementLevelExponent
  dx = (DomainRightEdge[0] - DomainLeftEdge[0]) /
    (TopGridDims[0]*POW(FLOAT(RefineBy),FLOAT(MaximumRefinementLevel)));
  
  for (i=0 ; i<nParticles; i++){
    pos = GalaxyParticleList[i]->ReturnPosition();
    rad = static_cast<ActiveParticleType_GalaxyParticle*>(
        GalaxyParticleList[i])->Radius;
    for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel)
      if (Temp->GridData->DepositRefinementZone(level,pos,rad) == FAIL) {
	ENZO_FAIL("Error in grid->DepositRefinementZone.\n")
	  }
  }

  return SUCCESS;
}

grid* ConstructFeedbackZone(ActiveParticleType* ThisParticle, int FeedbackRadius, 
			    FLOAT dx, HierarchyEntry** Grids, int NumberOfGrids,
			    int SendField);

int DistributeFeedbackZone(grid* FeedbackZones, HierarchyEntry** Grids, 
			   int NumberOfGrids, int SendField);

int ActiveParticleType_GalaxyParticle::SubtractMassFromGrid(int nParticles,
    ActiveParticleList<ActiveParticleType>& ParticleList, LevelHierarchyEntry *LevelArray[],
    FLOAT dx, int ThisLevel)
{

  // This function is a work in progress.
  return SUCCESS;

  /* Skip subtraction if we're not on the maximum refinement level. 
     This should only ever happen right after creation and then
     only in pathological cases where sink creation is happening at 
     the edges of two regions at the maximum refinement level.
     At any rate, this means that the subtraction may be delayed by a cycle
     at most.
   */

//   if (ThisLevel < MaximumRefinementLevel)
//     return SUCCESS;
// 
//   /* For each particle, loop over all of the grids and subtract
//      if the grid overlaps with the accretion zone. */
//   
//   int i, NumberOfGrids;
//   int *FeedbackRadius = NULL;
//   HierarchyEntry **Grids = NULL;
//   grid *particleGrid = NULL;
// 
//   bool ParticleIsOnThisProc, ParticleIsOnThisGrid;
//   
//   float SubtractedMass, SubtractedMomentum[3] = {};
//   
//   NumberOfGrids = GenerateGridArray(LevelArray, ThisLevel, &Grids);
//   
//   // FeedbackRadius is different for all galaxy particles, and we need to
//   // convert it to cells.
//   FeedbackRadius = new int[nParticles];
//   for (i = 0; i < nParticles; i++) {
//     FeedbackRadius[i] = int(static_cast<ActiveParticleType_GalaxyParticle*>(ParticleList[i])->Radius / dx);
//   }
//   grid** FeedbackZones = ConstructFeedbackZones(ParticleList, nParticles,
//     FeedbackRadius, dx, Grids, NumberOfGrids);
// 
//   for (i = 0; i < nParticles; i++) {
//     grid* FeedbackZone = FeedbackZones[i];
//     if (MyProcessorNumber == FeedbackZone->ReturnProcessorNumber()) {
//     
//         
//       if (FeedbackZone->AccreteOntoAccretingParticle(&ParticleList[i],FeedbackRadius[i]*dx,&AccretionRate) == FAIL)
// 	return FAIL;
//   
//       // No need to communicate the accretion rate to the other CPUs since this particle is already local.
//       static_cast<ActiveParticleType_AccretingParticle*>(ParticleList[i])->AccretionRate = AccretionRate;
//     }
//   }
//   
//   DistributeFeedbackZones(FeedbackZones, nParticles, Grids, NumberOfGrids);
// 
//   for (i = 0; i < nParticles; i++) {
//     delete FeedbackZones[i];    
//   }
// 
//   delete [] FeedbackZones;
// 
//   if (AssignActiveParticlesToGrids(ParticleList, nParticles, LevelArray) == FAIL)
//     return FAIL;
// 
//   delete [] Grids;
//   return SUCCESS;
}

int ActiveParticleType_GalaxyParticle::GalaxyParticleFeedback(int nParticles,
    ActiveParticleList<ActiveParticleType>& ParticleList, FLOAT dx, 
	LevelHierarchyEntry *LevelArray[], int ThisLevel, FLOAT period[MAX_DIMENSION])
{
  
  /* Skip if we're not on the maximum refinement level. 
     This should only ever happen right after creation and then
     only in pathological cases where creation is happening at 
     the edges of two regions at the maximum refinement level */

  if (ThisLevel < MaximumRefinementLevel)
    return SUCCESS;

  /* For each particle, loop over all of the grids and do feedback 
     if the grid overlaps with the feedback zone                   */
  
  int i, NumberOfGrids;
  int *FeedbackRadius = NULL;
  HierarchyEntry **Grids = NULL;
  
  NumberOfGrids = GenerateGridArray(LevelArray, ThisLevel, &Grids);
  
  FeedbackRadius = new int[nParticles];
  for (i = 0; i < nParticles; i++) {
    FeedbackRadius[i] = nint(static_cast<ActiveParticleType_GalaxyParticle*>(
            ParticleList[i])->Radius / dx);
  }
  
  for (i = 0; i < nParticles; i++) {
    grid* FeedbackZone = ConstructFeedbackZone(ParticleList[i], FeedbackRadius[i], dx, 
					       Grids, NumberOfGrids, ALL_FIELDS);

    if (MyProcessorNumber == FeedbackZone->ReturnProcessorNumber()) {
        
      if (FeedbackZone->ApplyGalaxyParticleFeedback(&ParticleList[i]) == FAIL)
	return FAIL;

    }

    DistributeFeedbackZone(FeedbackZone, Grids, NumberOfGrids, ALL_FIELDS);

    delete FeedbackZone;
  }

  delete [] FeedbackRadius;
  
  if (AssignActiveParticlesToGrids(ParticleList, nParticles, LevelArray) == FAIL)
    return FAIL;

  delete [] Grids;
  return SUCCESS;
}

int ActiveParticleType_GalaxyParticle::GalaxyParticleGravity(int nParticles,
    ActiveParticleList<ActiveParticleType>& ParticleList, FLOAT dx, 
	LevelHierarchyEntry *LevelArray[], int ThisLevel, FLOAT period[MAX_DIMENSION])
{
  
  /* Skip if we're not on the maximum refinement level. 
     This should only ever happen right after creation and then
     only in pathological cases where creation is happening at 
     the edges of two regions at the maximum refinement level */

  if (ThisLevel < MaximumRefinementLevel)
    return SUCCESS;

  /* For each particle, loop over all of the grids and do gravity 
     if the grid overlaps with the feedback zone                   */
  
  int i, NumberOfGrids;
  int *FeedbackRadius = NULL;
  HierarchyEntry **Grids = NULL;
  
  NumberOfGrids = GenerateGridArray(LevelArray, ThisLevel, &Grids);
  
  FeedbackRadius = new int[nParticles];
  for (i = 0; i < nParticles; i++) {
    FeedbackRadius[i] = nint(static_cast<ActiveParticleType_GalaxyParticle*>(
            ParticleList[i])->Radius / dx);
  }
  
  for (i = 0; i < nParticles; i++) {
    grid* FeedbackZone = ConstructFeedbackZone(ParticleList[i], FeedbackRadius[i], 
				     dx, Grids, NumberOfGrids, GRAVITATING_MASS_FIELD);

    if (MyProcessorNumber == FeedbackZone->ReturnProcessorNumber()) {
        
      if (FeedbackZone->ApplyGalaxyParticleGravity(&ParticleList[i]) == FAIL)
	return FAIL;

    }

    DistributeFeedbackZone(FeedbackZone, Grids, NumberOfGrids, GRAVITATING_MASS_FIELD);
    
    delete FeedbackZone;
  }

  delete [] FeedbackRadius;
  
  delete [] Grids;
  return SUCCESS;
}


/* 
 * This function can be used to reset the particle acceleration if required.
 * For example if a massless particle needs to be fixed in space. 
 * See ActiveParticle_RadiationParticle.C for details. 
 */
int ActiveParticleType_GalaxyParticle::ResetAcceleration(float *ActiveParticleAcceleration)
{
  return SUCCESS;
}

/*
 * For brute force creation of a particle. Useful for converting from 
 * star objects to active particles.
 */
int ActiveParticleType_GalaxyParticle::CreateParticle(grid *thisgrid_orig,
						      ActiveParticleFormationData &supp_data,
						      int particle_index)
{
  return SUCCESS;
} 

namespace {
  ActiveParticleType_info *GalaxyParticleInfo = 
    register_ptype <ActiveParticleType_GalaxyParticle> 
    ("GalaxyParticle");
}

std::vector<ParticleAttributeHandler*>
  ActiveParticleType_GalaxyParticle::AttributeHandlers;
