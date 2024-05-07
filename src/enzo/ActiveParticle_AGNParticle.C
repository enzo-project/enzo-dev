/***********************************************************************
/
/  A "AGN" active particle that compiles.
/
************************************************************************/

#include "ActiveParticle_AGNParticle.h"
#include "phys_constants.h"
#include "ActiveParticle.h"
#include "communication.h"
#include "TopGridData.h"
#include "macros_and_parameters.h"

#include <stdio.h>
/*
 * A 'friend' class can can be used to access private and protected member
 * variables of an object.  Here, we declare AGNParticleGrid, which is a trivial
 * subclass of the grid class and is a friend of the AGNParticle active particle.
 *
 * If we did not do this, we would need to modify grid.h to make it 'aware' of
 * the AGNParticle particle.  This way we can create new particle types without
 * modifying the grid class at all.
 */

class AGNParticleGrid : private grid {
  friend class ActiveParticleType_AGNParticle;
};

int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);

/* Note that we only refer to AGNParticleGrid here.
 * Given a grid object, we static cast to get this:
 * 
 * AGNParticleGrid *thisgrid =
 * static_cast<AGNParticleGrid *>(thisgrid_orig);
 * */

float ActiveParticleType_AGNParticle::static_cooling_radius= FLOAT_UNDEFINED;
float ActiveParticleType_AGNParticle::static_feedback_radius= FLOAT_UNDEFINED;


/* Initialize Particle Type
 *  * Called at startup
 *   */
int ActiveParticleType_AGNParticle::InitializeParticleType() {

  // if (debug)
  //    printf("Entering InitializeParticleType [%"ISYM"]\n", MyProcessorNumber);

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

      while(CellFlaggingMethod[method] != INT_UNDEFINED){
        method++;
      }
         

      CellFlaggingMethod[method] = 4; 
      } 

   TotalAGNParticlesCreated = 0;


  // This sets up the particle attributes that are defined on the base
  // ActiveParticleType class.  This includes thing like mass, position,
  // velocity, unique ID, etc.
  AttributeVector &ah = ActiveParticleType_AGNParticle::AttributeHandlers;
  ActiveParticleType::SetupBaseParticleAttributes(ah);

  // Register instance member variables here.  Warning: if you do not do this,
  // this data will not be saved to disk or communicated over MPI.
  ah.push_back(new Handler<ActiveParticleType_AGNParticle, float,
      &ActiveParticleType_AGNParticle::JetPhi>("JetPhi"));
  ah.push_back(new Handler<ActiveParticleType_AGNParticle, float,
      &ActiveParticleType_AGNParticle::JetAngle>("JetAngle"));
  ah.push_back(new Handler<ActiveParticleType_AGNParticle, float,
      &ActiveParticleType_AGNParticle::JetTheta>("JetTheta"));
  ah.push_back(new Handler<ActiveParticleType_AGNParticle, float,
      &ActiveParticleType_AGNParticle::Edot>("Edot"));
  ah.push_back(new Handler<ActiveParticleType_AGNParticle, float,
      &ActiveParticleType_AGNParticle::StoredEnergy>("StoredEnergy"));  
  ah.push_back(new Handler<ActiveParticleType_AGNParticle, float,
      &ActiveParticleType_AGNParticle::StoredMass>("StoredMass"));
  ah.push_back(new Handler<ActiveParticleType_AGNParticle, float,
      &ActiveParticleType_AGNParticle::CoolingRadius>("CoolingRadius"));
  ah.push_back(new Handler<ActiveParticleType_AGNParticle, float,
      &ActiveParticleType_AGNParticle::FeedbackRadius>("FeedbackRadius"));
  ah.push_back(new Handler<ActiveParticleType_AGNParticle, float,
      &ActiveParticleType_AGNParticle::CondensationFraction>("CondensationFraction"));
  ah.push_back(new Handler<ActiveParticleType_AGNParticle, float,
      &ActiveParticleType_AGNParticle::FeedbackEfficiency>("FeedbackEfficiency"));
  ah.push_back(new Handler<ActiveParticleType_AGNParticle, float,
      &ActiveParticleType_AGNParticle::TimescaleThreshold>("TimescaleThreshold"));
  ah.push_back(new Handler<ActiveParticleType_AGNParticle, float,
      &ActiveParticleType_AGNParticle::KineticFraction>("KineticFraction"));
  ah.push_back(new Handler<ActiveParticleType_AGNParticle, int,
      &ActiveParticleType_AGNParticle::FixedJetIsOn>("FixedJetIsOn"));
  ah.push_back(new Handler<ActiveParticleType_AGNParticle, float,
      &ActiveParticleType_AGNParticle::TimeOfLastShock>("TimeOfLastShock"));

  ah.push_back(new ArrayHandler<ActiveParticleType_AGNParticle, float, 256,
      &ActiveParticleType_AGNParticle::AccretionHistoryTime>("AccretionHistoryTime"));
  ah.push_back(new ArrayHandler<ActiveParticleType_AGNParticle, float, 256, 
      &ActiveParticleType_AGNParticle::AccretionHistoryRate>("AccretionHistoryRate"));


  if (debug)
      fprintf(stderr,"Leaving InitializeParticleType [%"ISYM"]\n", MyProcessorNumber);
  return SUCCESS;
}

bool ActiveParticleType_AGNParticle::AGNInsideGrid(grid* grid){
  float gridLeftX = grid->GridLeftEdge[0];
  float gridLeftY = grid->GridLeftEdge[1];
  float gridLeftZ = grid->GridLeftEdge[2];

  float gridRightX = grid->GridRightEdge[0];
  float gridRightY = grid->GridRightEdge[1];
  float gridRightZ = grid->GridRightEdge[2];

  bool in_grid_x = (gridLeftX <= AGNParticleInsert_x) && (AGNParticleInsert_x < gridRightX);
  bool in_grid_y = (gridLeftY <= AGNParticleInsert_y) && (AGNParticleInsert_y < gridRightY);
  bool in_grid_z = (gridLeftZ <= AGNParticleInsert_z) && (AGNParticleInsert_z < gridRightZ);

  bool inside_grid = in_grid_x && in_grid_y && in_grid_z;
  return inside_grid;
}

// Added by CJL
// For directly inserting the particle in the most highly refined subgrid
// that matches the desired position and time coordinates of the AGN
int ActiveParticleType_AGNParticle::InsertAGN
(grid *thisGrid, HierarchyEntry* SubgridPointer, ActiveParticleFormationData &data){

  /****************************************************************************/
  /*         Don't create more than 1 particle on a grid
  /****************************************************************************/
  // This should never actually happen but better safe than sorry I guess
  // This function should never be used for anything but the direct creation of 
  // a single AGN particle at a particular place and time in the simulation
  // For creating multiple AGN particles based on physical criteria, please
  // use/modify EvaluateFormation instead
  if (thisGrid->NumberOfActiveParticles >= 1)
  {
    return SUCCESS;
  }

  bool create_particle = false;
  float cell_x, cell_y, cell_z;
  float dx, dy, dz;
  bool in_x, in_y, in_z;
  float current_time = thisGrid->ReturnTime();
  float dt = thisGrid->ReturnTimeStep();
  

  // don't use this, it will give you the wrong value sometimes for some reason
  //int grid_level = thisGrid->GetLevel();
  //int grid_id = thisGrid->GetGridID();

  // use this instead
  int grid_level = data.level;
  int grid_id = data.GridID;
  

  /****************************************************************************/
  //                      Checking timestep   
  // If the insert time is outside of the acceptable formation time range, exit
  /****************************************************************************/
  
  bool in_timestep = (current_time <= AGNParticleInsert_Time) && 
                     (AGNParticleInsert_Time < (current_time+dt));

  if (!in_timestep){
    return SUCCESS;
  }
  

  /****************************************************************************/
  /*                      Checking grid boundary                    
  // If we're not inside the right grid, exit
  /****************************************************************************/
  // If the particle will not appear in this grid, exit
  
  
  if (!AGNInsideGrid(thisGrid)){
    return SUCCESS;
  }

  /****************************************************************************/
  /*                      Checking subgrid boundaries                    
  // If we're in the right place, try to figure out if the particle
  // should be placed on a subgrid instead
  /****************************************************************************/
  HierarchyEntry *Subgrid;
  for (Subgrid = SubgridPointer; Subgrid; Subgrid = Subgrid->NextGridThisLevel){

    if(AGNInsideGrid(Subgrid->GridData)){
      return SUCCESS;
    }
  }
    
  /****************************************************************************/
  /* At this point we have confirmed that the Insertion Time is now
  // The particle belongs on this grid and not on any of its subgrids
  // Now all we have left is to create the particle
  /****************************************************************************/

  CreateAGN(thisGrid, data);
  

  return SUCCESS;
}

int ActiveParticleType_AGNParticle::CreateAGN(grid *thisGrid, ActiveParticleFormationData &data){
  
  // we're gonna need these later
  int i, j, k, cell_index, method, MassRefinementMethod;
  float cell_x, cell_y, cell_z, dx, dy, dz;
  bool in_x, in_y, in_z;
  float *density = thisGrid->BaryonField[data.DensNum];
  float *velx = thisGrid->BaryonField[data.Vel1Num];
  float *vely = thisGrid->BaryonField[data.Vel2Num];
  float *velz = thisGrid->BaryonField[data.Vel3Num];
  float cell_mass;

  int GridDimension[3] = {thisGrid->GridDimension[0],
                          thisGrid->GridDimension[1],
                          thisGrid->GridDimension[2]};
                          
  /****************************************************************************/
  /*                      Find the right cell to place the particle in                    
  /****************************************************************************/
  bool create_particle_here = false;
  for (k = thisGrid->GridStartIndex[2]; k <= thisGrid->GridEndIndex[2]; k++)
  {
    for (j = thisGrid->GridStartIndex[1]; j <= thisGrid->GridEndIndex[1]; j++)
    {
      for (i = thisGrid->GridStartIndex[0]; i <= thisGrid->GridEndIndex[0]; i++)
      {

        

        cell_x = thisGrid->CellLeftEdge[0][i];
        cell_y = thisGrid->CellLeftEdge[1][j];
        cell_z = thisGrid->CellLeftEdge[2][k];

        dx = thisGrid->CellWidth[0][i];
        dy = thisGrid->CellWidth[1][j];
        dz = thisGrid->CellWidth[2][k];

        in_x = cell_x <= AGNParticleInsert_x && AGNParticleInsert_x < cell_x + dx;
        in_y = cell_y <= AGNParticleInsert_y && AGNParticleInsert_y < cell_y + dy;
        in_z = cell_z <= AGNParticleInsert_z && AGNParticleInsert_z < cell_z + dz;

        create_particle_here = in_x && in_y && in_z;
        if(create_particle_here){
          cell_index = GRIDINDEX_NOGHOST(i, j, k);
          break;
        }
      }
      if(create_particle_here) break;
    }
    if(create_particle_here) break;
  }

  /****************************************************************************/
  /*              We found the right cell, now create the particle                    
  /****************************************************************************/

  const int offset[] = {1, GridDimension[0], GridDimension[0]*GridDimension[1]};
  float DensUnits, LengthUnits, TempUnits, TimeUnits, VelUnits, Time;
  GetUnits(&DensUnits, &LengthUnits, &TempUnits, &TimeUnits, &VelUnits, AGNParticleInsert_Time);

  ActiveParticleType_AGNParticle *np = new ActiveParticleType_AGNParticle();
  data.NumberOfNewParticles++;

  // ok for some weird reason, the code is set up to add new particles to this 
  // data object, and then they get copied to the grid right after this function returns
  // why on earth would somebody throw the particles into this data object when they're 
  // just going to get added to the grid anyway is completely beyond me
  // my best guess is that the data object gets destroyed at the end and therefore
  // so does the original particle
  // is this some kind of defensive programming to make sure memory gets freed? 
  // why do particles need to be inserted as a copy?
  data.NewParticles.insert(*np);

  np->level = data.level;

  cell_mass = density[cell_index] * pow(thisGrid->CellWidth[0][0], 3.0);
  np->Mass = cell_mass; //TODO: Change this to use seed mass value from parameter file
  np->type = np->GetEnabledParticleID();
  np->BirthTime = thisGrid->ReturnTime();

  np->pos[0] = thisGrid->CellLeftEdge[0][i] + 0.5*thisGrid->CellWidth[0][i];
  np->pos[1] = thisGrid->CellLeftEdge[1][j] + 0.5*thisGrid->CellWidth[1][j];
  np->pos[2] = thisGrid->CellLeftEdge[2][k] + 0.5*thisGrid->CellWidth[2][k];

  np->CoolingRadius = AGNParticleAccretionRadiusKpc * 1.0e-3*Mpc_cm / LengthUnits;
  np->FeedbackRadius = AGNParticleFeedbackRadiusKpc * 1.0e-3*Mpc_cm / LengthUnits;
  np->CondensationFraction = AGNParticleCondensationFraction;
  np->FeedbackEfficiency = AGNParticleFeedbackEfficiency;
  np->TimescaleThreshold = 10.0;
  np->KineticFraction = AGNParticleKineticFraction;

  // These are supposed to control jet precession, not the orientation of the jet
  // Currently they don't actually do anything
  np->JetAngle = (20.0/360.0) * 2.0 * M_PI;
  np->JetPhi = 0.0*(15.0/360.0) * 2.0 * M_PI;
  np->JetTheta = (20.0/360.0)*2.0 * M_PI;

  np->FixedJetIsOn = 1; // 0 = Off, 1 = On
  np->TimeOfLastShock = 0.0;
  // AccretionHistory Bins
  for (int i = 0; i < 256; i++) {
    np -> AccretionHistoryTime[i] = np -> BirthTime - 
          (float)(255 - i) * 2.0 * AGNParticleAccretionDelayTime / 255.0;
    np -> AccretionHistoryRate[i] = 0.0;
  }

  static_cooling_radius = np -> CoolingRadius;
  static_feedback_radius = np -> FeedbackRadius;

  fprintf(stderr,"Creating new AGN Particle at [%"GSYM" %"GSYM" %"GSYM"] [%"ISYM"]\n", np->ReturnPosition()[0], 
          np->ReturnPosition()[1], np->ReturnPosition()[2], MyProcessorNumber);
  // modified by Deovrat Prasad, AGNParticle given initial velocity.
  if (HydroMethod == PPM_DirectEuler)
  {
    np->vel[0] = velx[cell_index];
    np->vel[1] = vely[cell_index];
    np->vel[2] = velz[cell_index];
  }
  else if (HydroMethod == Zeus_Hydro)
  {
    np->vel[0] = 0.5 * (velx[cell_index] + velx[cell_index+offset[0]]);
    np->vel[1] = 0.5 * (vely[cell_index] + vely[cell_index+offset[1]]);
    np->vel[2] = 0.5 * (velz[cell_index] + velz[cell_index+offset[2]]);
  } else {
    ENZO_FAIL("AGNParticle does not support RK Hydro or RK MHD");
  }
  np->Metallicity = 1.0;

  return SUCCESS;
}

int ActiveParticleType_AGNParticle::EvaluateFormation
(grid *thisgrid_orig, TopGridData *MetaData, ActiveParticleFormationData &data)
{       

  if(debug){
    fprintf(stderr,"Entering EvaluateFormation [%"ISYM"]\n", MyProcessorNumber);
  }

   // No need to do the rest if we're not on the maximum refinement level.
  if (data.level != MaximumRefinementLevel){
    return SUCCESS;
  }

  

  // This function assumes an algorithm exists that considers the local physics
  // of the grid before applying that algorithm to create an AGN particle
  // No such algorithm currently exists.
  // Instead consider using the InsertAGN function to directly place an AGN particle
  // where and when you want it to be. 
  //
  // In the future, when an algorithm exists to physically determine when/where an AGN
  // will spawn, consider implementing that algorithm within this function
  // Until then, this function should do nothing. 
  // -- CJL 2023
  bool in_development = true;
  if(in_development){
    return SUCCESS;
  }
  
  // Create a 'friend' grid alias that we can use to access private grid data.
  // AGNParticleGrid *thisGrid =
  //   static_cast<AGNParticleGrid *>(thisgrid_orig);

  // int i, j, k, index, method, MassRefinementMethod;
  // float *density = thisGrid->BaryonField[data.DensNum];
  // float *velx = thisGrid->BaryonField[data.Vel1Num];
  // float *vely = thisGrid->BaryonField[data.Vel2Num];
  // float *velz = thisGrid->BaryonField[data.Vel3Num];
  // float cell_mass;
  // int GridDimension[3] = {thisGrid->GridDimension[0],
  //                         thisGrid->GridDimension[1],
  //                         thisGrid->GridDimension[2]};

  // const int offset[] = {1, GridDimension[0], GridDimension[0]*GridDimension[1]};
  // float DensUnits, LengthUnits, TempUnits, TimeUnits, VelUnits; //added by DP

  // /****************************************************************************/
  // /*                      Writing AGN data out to file                        */
  // /****************************************************************************/

  // FILE * AGN_data_file = fopen("AGN_data.out", "a");
  
 
  // // Search for the cell nearest the origin
  // for (k = thisGrid->GridStartIndex[2]; k <= thisGrid->GridEndIndex[2]; k++)
  // {
  //   for (j = thisGrid->GridStartIndex[1]; j <= thisGrid->GridEndIndex[1]; j++)
  //   {
  //     for (i = thisGrid->GridStartIndex[0]; i <= thisGrid->GridEndIndex[0]; i++)
  //     {

  //       index = GRIDINDEX_NOGHOST(i, j, k);
  //       index = i + j * GridDimension[0] + k * GridDimension[1] * GridDimension[0];

        
          
  //       cell_x = thisGrid->CellLeftEdge[0][i];
  //       cell_y = thisGrid->CellLeftEdge[1][j];
  //       cell_z = thisGrid->CellLeftEdge[2][k];

  //       dx = thisGrid->CellWidth[0][i];
  //       dy = thisGrid->CellWidth[1][j];
  //       dz = thisGrid->CellWidth[2][k];

  //       in_x = cell_x <= AGNParticleInsert_x && AGNParticleInsert_x < cell_x + dx;
  //       in_y = cell_y <= AGNParticleInsert_y && AGNParticleInsert_y < cell_y + dy;
  //       in_z = cell_z <= AGNParticleInsert_z && AGNParticleInsert_z < cell_z + dz;

  //       create_particle = in_x && in_y && in_z;

  //       //Passed creation tests, create AGN particle
  //       if (create_particle) {
  //         // Get the units
  //         float DensUnits, LengthUnits, TempUnits, TimeUnits, VelUnits, Time;
  //         GetUnits(&DensUnits, &LengthUnits, &TempUnits, &TimeUnits, &VelUnits, AGNParticleInsert_Time);

  //         ActiveParticleType_AGNParticle *np = new ActiveParticleType_AGNParticle();
  //         data.NumberOfNewParticles++;
  //         data.NewParticles.insert(*np);

  //         np->level = data.level;
  //         np->GridID = data.GridID;
  //         np->CurrentGrid = thisGrid;

  //         cell_mass = density[index] * pow(thisGrid->CellWidth[0][0], 3.0);
  //         np->Mass = cell_mass;
	//   //np->Mass = 1.0e6*SolarMass/(DensUnits*LengthUnits*LengthUnits*LengthUnits); // Introduced by DP
  //         np->type = np->GetEnabledParticleID();
  //         np->BirthTime = thisGrid->ReturnTime();

  //         np->pos[0] = thisGrid->CellLeftEdge[0][i] + 0.5*thisGrid->CellWidth[0][i];
  //         np->pos[1] = thisGrid->CellLeftEdge[1][j] + 0.5*thisGrid->CellWidth[1][j];
  //         np->pos[2] = thisGrid->CellLeftEdge[2][k] + 0.5*thisGrid->CellWidth[2][k];

  //         np->CoolingRadius = AGNParticleAccretionRadiusKpc * 1.0e-3*Mpc_cm / LengthUnits;
  //         np->FeedbackRadius = AGNParticleFeedbackRadiusKpc * 1.0e-3*Mpc_cm / LengthUnits;
  //         np->CondensationFraction = AGNParticleCondensationFraction;
  //         np->FeedbackEfficiency = AGNParticleFeedbackEfficiency;
  //         np->TimescaleThreshold = 10.0;
  //         np->KineticFraction = AGNParticleKineticFraction;
  //         np->JetAngle = (20.0/360.0) * 2.0 * M_PI;
  //         np->JetPhi = 0.0*(15.0/360.0) * 2.0 * M_PI;
  //         np->JetTheta = (20.0/360.0)*2.0 * M_PI;
  //         np->FixedJetIsOn = 1; // 0 = Off, 1 = On
  //         np->TimeOfLastShock = 0.0;
  //         // AccretionHistory Bins
  //         for (int i = 0; i < 256; i++) {
  //                 np -> AccretionHistoryTime[i] = np -> BirthTime - 
  //                       (float)(255 - i) * 2.0 * AGNParticleAccretionDelayTime / 255.0;
  //                 np -> AccretionHistoryRate[i] = 0.0;
  //                 }

  //         static_cooling_radius = np -> CoolingRadius;
  //         static_feedback_radius = np -> FeedbackRadius;

  //         fprintf(AGN_data_file,"Creating new AGN Particle at [%"GSYM" %"GSYM" %"GSYM"] [%"ISYM"]\n", np->ReturnPosition()[0], 
  //                np->ReturnPosition()[1], np->ReturnPosition()[2], MyProcessorNumber);
  //         // modified by Deovrat Prasad, AGNParticle given initial velocity.
  //         if (HydroMethod == PPM_DirectEuler)
  //         {
  //           np->vel[0] = velx[index];
  //           np->vel[1] = vely[index];
  //           np->vel[2] = velz[index];
  //         }
  //         else if (HydroMethod == Zeus_Hydro)
  //         {
  //           np->vel[0] = 0.5 * (velx[index] + velx[index+offset[0]]);
  //           np->vel[1] = 0.5 * (vely[index] + vely[index+offset[1]]);
  //           np->vel[2] = 0.5 * (velz[index] + velz[index+offset[2]]);
  //         } else {
  //           ENZO_FAIL("AGNParticle does not support RK Hydro or RK MHD");
  //         }
  //         np->Metallicity = 1.0;

  //         thisGrid->AddActiveParticle(np);
  //         return SUCCESS;
  //       }
  //     } //i
  //   } //j
  // } // k

  // if (debug)
  //   fprintf(stderr, "Leaving EvaluateFormation [%"ISYM"]\n",
  // 	    MyProcessorNumber);

  // fclose(AGN_data_file);
  return SUCCESS;
}


/* Currently all feedback is done using HandleAGNFeedback which is 
 * called in BeforeEvolveLevel routine. 
 * This function does nothing.*/
int ActiveParticleType_AGNParticle::EvaluateFeedback
(grid *thisGrid_orig, ActiveParticleFormationData &data)
{
   //if (debug)
   //   printf("Entering EvaluateFeedback [%"ISYM"]\n", MyProcessorNumber);

   //if (debug)
   //   printf("Leaving EvaluateFeedback [%"ISYM"]\n", MyProcessorNumber);

  return SUCCESS;
}

void ActiveParticleType_AGNParticle::DescribeSupplementalData(ActiveParticleFormationDataFlags &flags)
{

  /*
   * For every entry in the ActiveparticleFormationData struct, we have a bool
   * here. If a flag is set to true, some derived data is calculated and attached
   * to the ActiveParticleFormationData struct.
   *
   * DarkMatterDensity, CoolingTime, Temperature, MetalField, H2Fraction, and
   * CoolingRate return an array containing the quantity. Since it is sometimes expensive
   * to cache these fields, they should only be turned on if your formation or feedback
   * algorithm require it.
   */

  flags.DarkMatterDensity = true;
  flags.CoolingTime = true;
  flags.Temperature = true;
  flags.MetalField = true;
  flags.UnitConversions = true;
  flags.DataFieldNumbers = true;
}


int ActiveParticleType_AGNParticle::SetFlaggingField(LevelHierarchyEntry *LevelArray[],int level,
    int TopGridDims[], int AGNParticleID)
{
  /* Generate a list of all AGNParticles in the simulation box */
   int i, nParticles;
   FLOAT* pos = NULL, dx = 0;
   ActiveParticleList<ActiveParticleType> AGNParticleList ;
   LevelHierarchyEntry* Temp = NULL;

   ActiveParticleFindAll(LevelArray, &nParticles, AGNParticleID, AGNParticleList);
  
   /* Calculate CellWidth on maximum refinement level */
   dx = (DomainRightEdge[0] - DomainLeftEdge[0]) /
        (TopGridDims[0] * POW(FLOAT(RefineBy), FLOAT(MaximumRefinementLevel)));

   float DensUnits, LengthUnits, TempUnits, TimeUnits, VelUnits, Time_loc;
   for (i = 0 ; i < nParticles; i++) {
      grid* tg = AGNParticleList[i] -> ReturnCurrentGrid();
      Time_loc = tg -> ReturnTime();

      if (GetUnits(&DensUnits, &LengthUnits, &TempUnits, &TimeUnits, &VelUnits, Time_loc) == FAIL) {
        ENZO_FAIL("Error in GetUnits.\n");
        }

      float ref_radius;
      ref_radius = AGNParticleFeedbackRadiusKpc * 1.0e-3*Mpc_cm / LengthUnits;
      //printf("ref_radius, cells: %"GSYM", %"GSYM"\n", ref_radius, ref_radius/dx);
            //exit(1);

      pos = AGNParticleList[i]->ReturnPosition();

      for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel)
         if (Temp->GridData->DepositRefinementZone(level, pos, ref_radius) == FAIL) {
            ENZO_FAIL("Error in grid->DepositRefinementZone.\n")
            } 
   }
   
   if (NumberOfProcessors > 1)
      for (i = 0; i < nParticles; i++) {
         delete AGNParticleList[i];
         AGNParticleList[i] = NULL;
         } 
   
   return SUCCESS;
}

/* Declaring function for creation fake grids for AGNParticle and distributing 
 * AGNparticle information back to real grids*/
grid* ConstructFeedbackZone(ActiveParticleType* ThisParticle, int FeedbackRadius,
                            FLOAT dx, HierarchyEntry** Grids, int NumberOfGrids,
                            int SendField);

int DistributeFeedbackZone(grid* FeedbackZone, HierarchyEntry** Grids,
                           int NumberOfGrids, int SendField);

/* This function, called from BeforeEvolveLevel, actually does the AGN feedback.
 * This function constructs a feedback grid and calls the DoAGNFeedback method
 *  on the grid. The result is then interpolated back onto the real grids.
 */   
int ActiveParticleType_AGNParticle::Handle_AGN_Feedback(int nParticles, ActiveParticleList<ActiveParticleType>& ParticleList,
      FLOAT dx, LevelHierarchyEntry* LevelArray[], int ThisLevel) {

   if (debug)
      printf("Entering Handle_AGN_Feedback [%"ISYM"]\n", MyProcessorNumber);

   if (ThisLevel < MaximumRefinementLevel)
      return SUCCESS;

   int i, NumberOfGrids;
   HierarchyEntry** Grids = NULL;

   int n_cells, n_cells_cooling, n_cells_feedback, n_disk_cells;

   NumberOfGrids = GenerateGridArray(LevelArray, ThisLevel, &Grids);

   // loop over particles and do feedback
   for (i=0; i < nParticles; i++) {
      float DensUnits, LengthUnits, TempUnits, TimeUnits, VelUnits, Time;
      
      grid* tg = ParticleList[i] -> ReturnCurrentGrid();
      Time = tg -> ReturnTime();

      // if (MyProcessorNumber != tg->ReturnProcessorNumber()) {
      //   continue;
      // }
      
      GetUnits(&DensUnits, &LengthUnits, &TempUnits, &TimeUnits, &VelUnits, Time);

      float code_racc, code_rfeed, code_diskr, code_diskd;
      code_racc = AGNParticleAccretionRadiusKpc * kpc_cm / LengthUnits;
      code_rfeed = AGNParticleFeedbackRadiusKpc * kpc_cm / LengthUnits;
      code_diskr = AGNParticleDiskRadius * kpc_cm / LengthUnits;
      code_diskd = AGNParticleDiskDistance * kpc_cm / LengthUnits;

      // calculating size of the feedback zone
      n_cells_cooling = floor(code_racc/ dx) + 1;
      n_cells_feedback = floor(code_rfeed/ dx) + 1;
      n_disk_cells = floor((code_diskr + code_diskd) / dx) + 1;
      
      n_cells = max(n_cells_cooling, n_cells_feedback);

      // Contructing feedback zone
      grid* FeedbackZone = ConstructFeedbackZone(ParticleList[i],
                           n_cells, dx, Grids, NumberOfGrids, ALL_FIELDS);

      // Calling DoAGNFeedback function which takes care of feedback physics
       
      // if (MyProcessorNumber == FeedbackZone->ReturnProcessorNumber()) {
         if( FeedbackZone -> DoAGNFeedback(ParticleList[i]) == FAIL)
           return FAIL;
      //} 
 
      // Copy the feedback grid to the real grids
      DistributeFeedbackZone(FeedbackZone, Grids, NumberOfGrids, ALL_FIELDS);

      delete FeedbackZone;


   } // end of loop over particles
   if (AssignActiveParticlesToGrids(ParticleList, nParticles, LevelArray) == FAIL)
      return FAIL;

   delete [] Grids;

   // All done
   if (debug)
      printf("Leaving Handle_AGN_Feedback [%"ISYM"]\n", MyProcessorNumber);

   return SUCCESS;
}

 /* 
 * This function can be used to reset the particle acceleration if required.
 * For example if a massless particle needs to be fixed in space. 
 * See ActiveParticle_RadiationParticle.C for details. 
 */
int ActiveParticleType_AGNParticle::ResetAcceleration(float *ActiveParticleAcceleration)
{
  ActiveParticleAcceleration[0] = 0.0;
  ActiveParticleAcceleration[1] = 0.0;
  ActiveParticleAcceleration[2] = 0.0;

  return SUCCESS;
}

/*
 * For brute force creation of a particle. Useful for converting from 
 * star objects to active particles.
 */
int ActiveParticleType_AGNParticle::CreateParticle(grid *thisgrid_orig,
						ActiveParticleFormationData &supp_data,
						int particle_index)
{
  return SUCCESS;
} 

namespace
{

  /*
   * This creates the ActiveParticleType_info singleton for the "AGNParticle"
   * particle type.  This object will be used elsewhere to interface with the
   * active particle API in a manner that is agnostic to the underlying particle
   * type. Without this line, the particle name will not be recognized when used
   * in a parameter file.
   */

  ActiveParticleType_info *AGNParticleInfo =
    register_ptype <ActiveParticleType_AGNParticle> ("AGNParticle");
}

// Instantiate the AttributeHandler singleton for this particle type.

std::vector<ParticleAttributeHandler*>
  ActiveParticleType_AGNParticle::AttributeHandlers;
