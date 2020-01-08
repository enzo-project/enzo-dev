/***********************************************************************
/
/  ADD FEEDBACK TO RADIAL PROFILE OVER MULTIPLE GRIDS
/
/  written by: John Wise
/  date:       September, 2005
/  modified1: Ji-hoon Kim
/             October, 2009
/
/ PURPOSE: To apply feedback effects, we must consider multiple grids
/          since sometimes the feedback radius often exceeds the grid
/          boundaries.  This routine makes sure that all of the grids
/          have the same code time to ensure consistency.
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "performance.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "list.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"

#include "phys_constants.h"

#define MAX_TEMPERATURE 1e8

int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);
int RecalibrateMBHFeedbackThermalRadius(FLOAT star_pos[], LevelHierarchyEntry *LevelArray[],
                                        int level, float &Radius,
                                        double &EjectaDensity, double &EjectaMetalDensity,
                                        double &EjectaThermalEnergy);
int RemoveParticles(LevelHierarchyEntry *LevelArray[], int level, int ID);
FLOAT FindCrossSection(int type, float energy);

int StarParticleAddFeedback(TopGridData *MetaData,
                            LevelHierarchyEntry *LevelArray[], int level,
                            Star *&AllStars, bool *&AddedFeedback)

{

    Star *cstar;
    bool MarkedSubgrids = false;
    bool SphereCheck;
    int i, l, dim, temp_int, SkipMassRemoval, SphereContained,
        SphereContainedNextLevel, dummy, count;
    float influenceRadius, RootCellWidth, SNe_dt, dtForThisStar, MassLoss;
    double EjectaThermalEnergy, EjectaDensity, EjectaMetalDensity;
    FLOAT Time;
    LevelHierarchyEntry *Temp;

    if (AllStars == NULL)
        return SUCCESS;

    LCAPERF_START("StarParticleAddFeedback");

    /* Get time and SNe timestep */

    Temp = LevelArray[level];
    Time = Temp->GridData->ReturnTime();
    if (LastSupernovaTime < 0)
        SNe_dt = 0.0;
    else
        SNe_dt = Time - LastSupernovaTime;
    LastSupernovaTime = Time;
    RootCellWidth = 1.0 / MetaData->TopGridDims[0];

    /* Set the units. */

    float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,
        VelocityUnits;
    GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
             &TimeUnits, &VelocityUnits, Time);

    count = 0;

    for (cstar = AllStars; cstar; cstar = cstar->NextStar, count++)
    {

        AddedFeedback[count] = false;

        /* Special case for "normal" star particles to account for mass
       loss through supernovae. */

        if (cstar->ReturnType() == NormalStar &&
            cstar->ReturnLevel() == level)
        {
            MassLoss = cstar->CalculateMassLoss(SNe_dt);
            cstar->SetAccretionMass(-MassLoss);
        }

        if ((cstar->ReturnFeedbackFlag() != MBH_THERMAL) &&
            (cstar->ReturnFeedbackFlag() != MBH_JETS) &&
            !cstar->ApplyFeedbackTrue(SNe_dt))
            continue;

        dtForThisStar = LevelArray[level]->GridData->ReturnTimeStep();

        /* Compute some parameters */

        cstar->CalculateFeedbackParameters(influenceRadius, RootCellWidth, SNe_dt, EjectaDensity,
                                           EjectaThermalEnergy, EjectaMetalDensity, DensityUnits, LengthUnits,
                                           TemperatureUnits, TimeUnits, VelocityUnits, dtForThisStar,
                                           Time, SphereCheck);

        if (SphereCheck)
        {

            /* Recalibrate MBHFeedbackThermalRadius if requested */

            if (cstar->ReturnFeedbackFlag() == MBH_THERMAL)
                RecalibrateMBHFeedbackThermalRadius(cstar->ReturnPosition(), LevelArray, level, influenceRadius,
                                                    EjectaDensity, EjectaMetalDensity, EjectaThermalEnergy);

            /* Determine if a sphere with enough mass (or equivalently radius
       for SNe) is enclosed within grids on this level */

            LCAPERF_START("star_FindFeedbackSphere");
            cstar->FindFeedbackSphere(LevelArray, level, influenceRadius, EjectaDensity, EjectaThermalEnergy,
                                      SphereContained, SkipMassRemoval, DensityUnits, LengthUnits,
                                      TemperatureUnits, TimeUnits, VelocityUnits, Time, MarkedSubgrids);
            LCAPERF_STOP("star_FindFeedbackSphere");

            /* If the particle already had sufficient mass, we still want to
       mark this particle to activate it. */

            if (SkipMassRemoval == TRUE)
                AddedFeedback[count] = true;

            /* If there's no feedback or something weird happens, don't bother. */

            if (influenceRadius <= tiny_number ||
                SphereContained == FALSE ||
                ((cstar->ReturnFeedbackFlag() == MBH_THERMAL ||
                  cstar->ReturnFeedbackFlag() == MBH_JETS) &&
                 (influenceRadius >= RootCellWidth / 2 ||
                  EjectaThermalEnergy <= tiny_number)))
                continue;

            /* Determine if a sphere is enclosed within the grids on next level
       If that is the case, we perform AddFeedbackSphere not here, 
       but in the EvolveLevel of the next level. */

            SphereContainedNextLevel = FALSE;

            LCAPERF_START("star_FindFeedbackSphere2");
            if ((cstar->ReturnFeedbackFlag() == MBH_THERMAL ||
                 cstar->ReturnFeedbackFlag() == MBH_JETS ||
                 cstar->ReturnFeedbackFlag() == CONT_SUPERNOVA) &&
                LevelArray[level + 1] != NULL)
                cstar->FindFeedbackSphere(LevelArray, level + 1, influenceRadius, EjectaDensity, EjectaThermalEnergy,
                                          SphereContainedNextLevel, dummy, DensityUnits, LengthUnits,
                                          TemperatureUnits, TimeUnits, VelocityUnits, Time, MarkedSubgrids);
            LCAPERF_STOP("star_FindFeedbackSphere2");

            //    if (debug) {
            //      fprintf(stdout, "EjectaDensity=%g, influenceRadius=%g\n", EjectaDensity, influenceRadius);
            //      fprintf(stdout, "SkipMassRemoval=%d, SphereContained=%d, SphereContainedNextLevel=%d\n",
            //	      SkipMassRemoval, SphereContained, SphereContainedNextLevel);
            //    }

            /* Quit this routine when 
       (1) sphere is not contained, or 
       (2) sphere is contained, but the next level can contain the sphere, too. */
            if ((SphereContained == FALSE) ||
                (SphereContained == TRUE && SphereContainedNextLevel == TRUE))
                continue;

        } // ENDIF SphereCheck
        else
        {

            /* When the sphere is completely confined in a grid, only apply
	 feedback at the level at which the star exists. */

            if (level != cstar->ReturnLevel())
                continue;
        }

        /* Now set cells within the radius to their values after feedback.
       While walking through the hierarchy, look for particle to
       change their properties to post-feedback values. */

        int CellsModified = 0;

        if (SkipMassRemoval == FALSE)
        {

            /* Determine the H-ionizing photon luminosity to calculate the
	 photo-ionization and heating rate in the initial Stroemgren
	 sphere. */

            int nbins;
            double Q[MAX_ENERGY_BINS], Q_HI, sigma;
            float energies[MAX_ENERGY_BINS], deltaE;
#ifdef TRANSFER
            if (RadiativeTransfer)
            {
                cstar->ComputePhotonRates(TimeUnits, Time, nbins, energies, Q, dtForThisStar);
                sigma = (double)FindCrossSection(0, energies[0]); // HI (cm^2)
                Q_HI = Q[0];
                deltaE = energies[0] - 13.6; // eV
            }
            else
#endif /* TRANSFER */
            {
                Q_HI = 0.0;
                sigma = 0.0;
                deltaE = 0.0;
            }
            float rho = EjectaDensity;
            float z_rho = EjectaMetalDensity;
            FLOAT AVL0 = 0.0;
            float metalAccreted = 0.0;
            float metalFrac;
            float rescale = 1.0;
            for (l=level; l < MAX_DEPTH_OF_HIERARCHY; l++){ // initially l=level; l<MAX; l++
                if (!LevelArray[l]) continue;

                if (cstar->ReturnFeedbackFlag() == SUPERNOVA || cstar->ReturnFeedbackFlag() == FORMATION)
                {
                                
                    /*
                        Spheres interacting with grids isnt consistent; Do a first pass with no deposition
                        to validate the volume we will deposit into, then rescale the deposition accordingly.
                        --AIW
                    */
                        // only rescale SN if theyre PopIII single supernova--rescaling PopII is too expensive since
                        // it happens at every timestep.
                    bool rescaleSN = cstar->ReturnFeedbackFlag()==SUPERNOVA && (cstar->ReturnType() == PopIII);
                    int nCells = 0;
                    FLOAT vol_modified = 0.0;
                    float mass_dep = 0.0;
                    float metal_dep = 0.0; // track P3 metal
                    float metal2_dep = 0.0; // track p2 metal
                    if (rescaleSN || cstar->ReturnFeedbackFlag() == FORMATION)
                        for (Temp = LevelArray[l]; Temp; Temp = Temp->NextGridThisLevel)
                            Temp->GridData->StarParticleCalculateFeedbackVolume(
                                                        cstar, l, influenceRadius, DensityUnits, 
                                                        LengthUnits,
                                                        VelocityUnits, TemperatureUnits, TimeUnits, nCells,
                                                        mass_dep, metal_dep, metal2_dep, vol_modified);
                    FLOAT AllVol = 0;
                    FLOAT old_vol = influenceRadius * influenceRadius * influenceRadius 
                                        * 4.0 * pi / 3.0;
                            
                            // sum volume across all processes for this level 
                            //      (sphere can be across procs as well!) -AIW

                    if (rescaleSN || cstar->ReturnFeedbackFlag() == FORMATION){
                        MPI_Allreduce(&vol_modified, &AllVol,1,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    if (AVL0 == 0)
                        AVL0 = AllVol;
                    float allMass = 0;
                    float allMetal = 0; // track p3 metal
                    float allMetal2 = 0; // track p2 metal
                    FLOAT vCell = LevelArray[l]->GridData->GetVCell();
                    FLOAT MassUnits = DensityUnits*pow(LengthUnits,3)/SolarMass; //code -> Msun
                    // if forming mass, need to check that mass accreting from grid is consistent
                    if (cstar->ReturnFeedbackFlag() == FORMATION && rho == EjectaDensity && AVL0 > 0){ // set all this on the first valid pass
                        /* sum quantities across tasks */
                        MPI_Allreduce(&vol_modified, &AllVol,1,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                        MPI_Allreduce(&mass_dep, &allMass,1,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                        MPI_Allreduce(&metal_dep, &allMetal,1,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                        MPI_Allreduce(&metal2_dep, &allMetal2,1,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                        // Set the densities to constant value for the interior of Stromgren sphere
                        allMass *= MassUnits;
                        allMetal *= MassUnits;
                        allMetal2 *= MassUnits;
                        metalFrac = allMetal/(allMetal+allMetal2); //fraction of p3/all metals
                        printf("[L = %d : type = %d] AllMass = %g AllZ = %g z_rho = %g New FZ = %g volume = %g Zstar = %g\n",
                                    l, cstar->ReturnType(), allMass, allMetal, allMetal/MassUnits*vCell, 
                                    (allMetal+allMetal2)/allMass, AVL0*pow(LengthUnits,3), cstar->ReturnMetallicity());
                        if (allMass > cstar->ReturnFinalMass() && AVL0 > 0){
                            rho = (allMass - cstar->ReturnFinalMass())/MassUnits/AVL0*vCell;
                            z_rho = max((allMetal+allMetal2-(cstar->ReturnFinalMass()*cstar->ReturnMetallicity()))/MassUnits*vCell/AVL0, 1e-20*rho);
                            printf("New densities rho=%g z_rho=%g M = %g Mz = %g\n",rho,z_rho,
                                   rho*MassUnits/SolarMass*pow(LengthUnits,3)*vCell,
                                   z_rho*MassUnits/SolarMass*pow(LengthUnits,3)*vCell);
                            cstar->SetMetallicity((allMetal+allMetal2)/allMass);
                        
                        }
                        // if (allMass > cstar->ReturnFinalMass() && level < l){
                        //     rho *= 0.578703704;
                        //     z_rho *= 0.578703704;
                        // }
                    }



                        if (vol_modified > 0.0 && rescaleSN && debug)
                            printf("level %d Prior Deposition Masses M = %g Z = %g AllVol = %g rankVol = %g ratio = %f r_mass = %g r_z = %g\n ",
                                l, AllVol*EjectaDensity*pow(LengthUnits,3)*DensityUnits/SolarMass, 
                                AllVol*EjectaMetalDensity*pow(LengthUnits,3)*DensityUnits/SolarMass,
                                AllVol*pow(LengthUnits,3), vol_modified*pow(LengthUnits,3),
                                old_vol/AllVol,
                                vol_modified*EjectaDensity*pow(LengthUnits,3)*DensityUnits/SolarMass, 
                                vol_modified*EjectaMetalDensity*pow(LengthUnits,3)*DensityUnits/SolarMass
                                );
                            
                    }// endif rescale or formation
                    if (rescaleSN){
                        if (l==level || AVL0 > 0){
                                rescale = (old_vol/AVL0);
                                rho = EjectaDensity * rescale;
                                z_rho = EjectaMetalDensity * rescale;
                            }
                        // its a damn dirty hack, but it works. 
                        if (l > level)
                        {
                            rho = EjectaDensity*rescale;//*AllVol/AVL0 ;
                            z_rho = EjectaMetalDensity*rescale;//*AllVol/AVL0;
                        }
                        if (rescale < 1.0 && AllVol > 0)
                            fprintf(stdout, "\n\n[ %d ]Rescaling volume on level %d v = %g/%g  lratio = %f rho = %g/%g z_rho=%g/%g m_d = %g m_z = %g\n\n\n",
                                cstar->ReturnFeedbackFlag(), l, AVL0*pow(LengthUnits,3), 
                                old_vol*pow(LengthUnits,3), AllVol/AVL0, rho * DensityUnits, EjectaDensity*DensityUnits, 
                                z_rho * DensityUnits, EjectaMetalDensity*DensityUnits, 
                                rho*AllVol*DensityUnits*pow(LengthUnits,3)/SolarMass,
                                z_rho*AllVol*DensityUnits*pow(LengthUnits,3)/SolarMass);
                    }
                }
                        /* do the real deposition */

                for (Temp = LevelArray[l]; Temp; Temp = Temp->NextGridThisLevel)
                    Temp->GridData->AddFeedbackSphere(cstar, l, 
                                        influenceRadius,            
                                        DensityUnits, LengthUnits,
                                        VelocityUnits, TemperatureUnits, TimeUnits, rho, z_rho,
                                        EjectaThermalEnergy, Q_HI, sigma, deltaE,
                                        CellsModified, metalFrac);
            }        

        } // ENDIF

        //    fprintf(stdout, "StarParticleAddFeedback[%"ISYM"][%"ISYM"]: "
        //	    "Radius = %e pc_cm, changed %"ISYM" cells.\n",
        //	    cstar->ReturnID(), level, influenceRadius*LengthUnits/pc_cm, CellsModified);

        /* Remove mass from the star that is added to grids. Also, because EjectaDensity 
       is added with zero net momentum, increase the particle's velocity accordingly. 
       Only for MBH_JETS; currently this is done in Grid_AddFeedbackSphere.C */

        /*
    if (EjectaDensity != 0 && CellsModified > 0)
      if (cstar->ReturnFeedbackFlag() == MBH_THERMAL ||
	  cstar->ReturnFeedbackFlag() == MBH_JETS)
	cstar->RemoveMassFromStarAfterFeedback(influenceRadius, EjectaDensity, 
					       DensityUnits, LengthUnits, CellsModified);
    */

        /* Only kill a Pop III star after it has gone SN */

        if (cstar->ReturnFeedbackFlag() == SUPERNOVA)
            cstar->SetFeedbackFlag(DEATH);

        /* We only color the fields once */

        AddedFeedback[count] = true;

#ifdef UNUSED
        temp_int = CellsModified;
        MPI_Reduce(&temp_int, &CellsModified, 1, MPI_INT, MPI_SUM, ROOT_PROCESSOR,
                   MPI_COMM_WORLD);

        if (debug)
        {
            if (cstar->ReturnFeedbackFlag() != FORMATION)
                fprintf(stdout, "StarParticleAddFeedback[%" ISYM "][%" ISYM "]: "
                                "Radius = %" GSYM " pc\n",
                        cstar->ReturnID(), level, influenceRadius * LengthUnits / pc_cm);
            if (cstar->ReturnFeedbackFlag() == DEATH ||

                cstar->ReturnFeedbackFlag() == CONT_SUPERNOVA ||
                cstar->ReturnFeedbackFlag() == MBH_THERMAL ||
                cstar->ReturnFeedbackFlag() == MBH_JETS)
                fprintf(stdout, "StarParticleAddFeedback[%" ISYM "][%" ISYM "]: "
                                "Energy = %" GSYM "  , skip = %" ISYM "\n",
                        cstar->ReturnID(), level, EjectaThermalEnergy, SkipMassRemoval);
            fprintf(stdout, "StarParticleAddFeedback[%" ISYM "][%" ISYM "]: "
                            "changed %" ISYM " cells.  AddedFeedback[%d] = %d\n",
                    cstar->ReturnID(), level, CellsModified,
                    count, AddedFeedback[count]);
        }
#endif

    } // ENDFOR stars

    for (level = MaximumRefinementLevel; level > 0; level--){
            Temp = LevelArray[level];
            while (Temp != NULL) {
                if (Temp->GridData->ProjectSolutionToParentGrid(*Temp->GridHierarchyEntry->ParentGrid->GridData) == FAIL){
                fprintf(stderr, "Error in grid->ProjectSolutionToParentGrid\n");
                return FAIL;
                }
                Temp = Temp->NextGridThisLevel;
            }
    }
    LCAPERF_STOP("StarParticleAddFeedback");
    return SUCCESS;
}
