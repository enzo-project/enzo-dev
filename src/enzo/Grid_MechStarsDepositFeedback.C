#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <string.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "CosmologyParameters.h"
#include "StarParticleData.h"
#include "phys_constants.h"

extern "C" void FORTRAN_NAME(cic_deposit)(float *xPosition, float *yPosition,
                                          float *zPosition, int *gridRank, int *nParticles,
                                          float *DepositQuantity, float *FieldToDepositTo,
                                          float *leftEdge, int *xGridDim, int *yGridDim,
                                          int *zGridDim, float *gridDx, float *cloudsize);
int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, float *MassUnits, float Time);
int transformComovingWithStar(float *Density, float *Metals,
                              float *MetalsSNII, float *MetalsSNIA,
                              float *Vel1, float *Vel2, float *Vel3,
                              float *TE, float *GE,
                              float up, float vp, float wp,
                              int sizeX, int sizeY, int sizeZ, int direction);
int FindField(int field, int farray[], int numfields);
float Window(float xd, float yd, float zd, float width, bool NGP);
int grid::MechStars_DepositFeedback(float ejectaEnergy,
                                    float ejectaMass, float ejectaMetal,
                                    float *totalMetals, float *temperature,
                                    float *up, float *vp, float *wp,
                                    float *xp, float *yp, float *zp,
                                    int ip, int jp, int kp,
                                    int size, float *muField, int winds, int nSNII,
                                    int nSNIA, float starMetal, int isP3)
{

    /*
     This routine will create an isocahedron of coupling particles, where we determine
        the feedback quantities.  The vertices of the isocahedron are ~coupled particles
        and all have radius dx from the source particle. 
        Each vertex particle will then be CIC deposited to the grid!
    */
    //printf("In Feedback deposition\n");
    if (MyProcessorNumber != ProcessorNumber)
        return 0;
    bool debug = true;
    bool criticalDebug = true;
    float min_winds = 1.0;
    bool printout = debug & !winds;
    int index = ip + jp * GridDimension[0] + kp * GridDimension[0] * GridDimension[1];
    if (printout)
        printf("Host index = %d\n", index);
    int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
    /*
        *
        *  flag for NGP deposition
        * 
    */
    bool NGP = false; // true to use NGP deposition instead of CIC
    /*
        *
        *
        * 
    */
    float ntouched = NGP ? (27) : (64); // how many cells get touched by deposition? 
    /* Compute size (in floats) of the current grid. */
    float stretchFactor = 1.0; //1.5/sin(M_PI/10.0);  // How far should cloud particles be from their host
                               // in units of dx. Since the cloud forms a sphere shell, stretchFactor > 1 is not recommended
    size = 1;
    for (int dim = 0; dim < GridRank; dim++)
        size *= GridDimension[dim];
    float DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
        HMNum, H2INum, H2IINum, DINum, DIINum, HDINum;
    /* Find fields: density, total energy, velocity1-3. */

    if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                         Vel3Num, TENum) == FAIL)
    {
        fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
        return FAIL;
    }
    /* Set the units */
    float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1,
          TimeUnits = 1, VelocityUnits = 1, MassUnits = 1;
    if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
                 &TimeUnits, &VelocityUnits,
                 &MassUnits, this->Time) == FAIL)
    {
        fprintf(stderr, "Error in GetUnits.\n");
        return FAIL;
    }
    FLOAT dx = CellWidth[0][0];

    if (printout)
        fprintf(stdout, "depositing quantities: Energy %e, Mass %e, Metals %e\n",
               ejectaEnergy, ejectaMass, ejectaMetal);

    /* 
        get metallicity field and set flag; assumed true thoughout feedback
        since so many quantities are metallicity dependent
     */
    int MetallicityField = FALSE, MetalNum = -1, MetalIaNum = -1, MetalIINum = -1, SNColourNum = -1;
    if ((MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields)) != -1)
        MetallicityField = TRUE;
    else
    {
        fprintf(stdout, "MechStars only functions with metallicity field enabled!");
        ENZO_FAIL("Grid_MechStarsDepositFeedback: 91");
        MetalNum = 0;
    }
    if (StarMakerTypeIaSNe)
        MetalIaNum = FindField(MetalSNIaDensity, FieldType, NumberOfBaryonFields);
    if (StarMakerTypeIISNeMetalField)
        MetalIINum = FindField(MetalSNIIDensity, FieldType, NumberOfBaryonFields);
    if (MechStarsSeedField)
        SNColourNum = FindField(SNColour, FieldType, NumberOfBaryonFields);

    if (MechStarsSeedField && SNColourNum < 0)
        ENZO_FAIL("Cant use Seed Field without SNColour field");
    /* set other units that we need */
    MassUnits = DensityUnits * pow(LengthUnits * dx, 3) / SolarMass;                                 //Msun!
    float EnergyUnits = DensityUnits * pow(LengthUnits * dx, 3) * pow(VelocityUnits, 2.0); //[g cm^2/s^2] -> code_energy
    float MomentaUnits = MassUnits * VelocityUnits / 1e5;

    /* Make copys of pointers fields to work with (theyre just easier to type!). */
    float *density = BaryonField[DensNum];
    float *metals = BaryonField[MetalNum];
    float *metalsII = NULL;
    float *metalsIA = NULL;
    float *metalsIII = NULL;
    if (StarMakerTypeIISNeMetalField)
        metalsII = BaryonField[MetalIINum];
    if (StarMakerTypeIaSNe)
        metalsIA = BaryonField[MetalIaNum];
    if (MechStarsSeedField)
        metalsIII = BaryonField[SNColourNum];
    float *u = BaryonField[Vel1Num];
    float *v = BaryonField[Vel2Num];
    float *w = BaryonField[Vel3Num];
    float *totalEnergy = BaryonField[TENum];
    float *gasEnergy = BaryonField[GENum];

    float phi = (1.0 + sqrt(5)) / 2.0; //Golden Ratio
    float iphi = 1.0 / phi;            // inverse GR
    
    /* Coupling Particle Vectors: A DODECAHEDRON+ISOCAHEDRON */
    int nCouple = 26;
    float A = stretchFactor * dx;
    FLOAT cloudSize = stretchFactor * dx;

    /* 
        Make a cloud of coupling particles; 
        each is on a grid with dx spacing
        from the host particles.
    */

    FLOAT CloudParticleVectorX[] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    // {1, 1, 1, 1, -1, -1, -1, -1, 0, 0, 0, 0,
    // iphi, iphi, -iphi, -iphi, phi, phi, -phi, -phi,
    // 0, 0, 0, 0, 1, 1, -1, -1, phi,-phi, phi,-phi};
    FLOAT CloudParticleVectorY[] = {1, 1, 1, 0, 0, -1, -1, -1, 0, 1, 1, 1, 0, 0, -1, -1, -1, 1, 1, 1, 0, 0, -1, -1, -1, 0};
    // {1,1,-1,-1, 1, 1, -1, -1, iphi, iphi, -iphi, -iphi,
    // phi, -phi, phi,-phi, 0, 0, 0, 0,1, 1, -1, -1,
    // phi, -phi, -phi, phi, 0, 0, 0, 0};
    FLOAT CloudParticleVectorZ[] = { 1, 0, -1, 1, -1, 1, 0, -1, 0, 1, 0, -1, 1, -1, 1, 0, -1, 1, 0, -1, 1, -1, 1, 0, -1, 0};
    // {1,-1, 1,-1, 1,-1, 1,-1, phi,-phi, phi,-phi,
    // 0, 0, 0, 0, iphi, -iphi, iphi, -iphi,
    // phi, -phi, -phi, phi, 0, 0, 0, 0, 1, 1, -1, -1};
    float weightsVector[nCouple];
    /* Set position of feedback cloud particles */

    FLOAT CloudParticlePositionX[nCouple];
    FLOAT CloudParticlePositionY[nCouple];
    FLOAT CloudParticlePositionZ[nCouple];

    /*all possible values of x,y,z with origin particle at x=y=z=0.0 */

    for (int cpInd = 0; cpInd < nCouple; cpInd++)
    {
        float norm = sqrt(CloudParticleVectorX[cpInd] * CloudParticleVectorX[cpInd] + CloudParticleVectorY[cpInd] * CloudParticleVectorY[cpInd] + CloudParticleVectorZ[cpInd] * CloudParticleVectorZ[cpInd]);
        float xbaMag = A * A * norm * norm;
        /* in this cloud, take the coupling particle position as 0.5, 0.5, 0.5 */
        CloudParticlePositionX[cpInd] = *xp + CloudParticleVectorX[cpInd] / norm *
                                                  A;
        CloudParticlePositionY[cpInd] = *yp + CloudParticleVectorY[cpInd] / norm *
                                                  A;
        CloudParticlePositionZ[cpInd] = *zp + CloudParticleVectorZ[cpInd] / norm *
                                                  A;
        weightsVector[cpInd] = 1.0 ; //0.5 * (1. - 1. / (1. + 1. / 4. / M_PI / 26. / xbaMag / M_PI));
        /* turn the vectors back into unit-vectors */
        CloudParticleVectorZ[cpInd] /= norm;
        CloudParticleVectorY[cpInd] /= norm;
        CloudParticleVectorX[cpInd] /= norm;
    }
    float weightsSum = 0.0;
    for (int wind = 0; wind < nCouple; wind++)
    {
        weightsSum += weightsVector[wind];
    }
    for (int wind = 0; wind < nCouple; wind++)
    {
        weightsVector[wind] /= weightsSum;
        if (weightsVector[wind] == 0 || isnan(weightsVector[wind]))
        {
            ENZO_FAIL("NaN weight Vector!")
        }
    }

    /* transform to comoving with the star and take velocities to momenta.
        Take Energy densities to energy.  winds are ngp unless VERY high resolution
     */
    // if (!winds || dx * LengthUnits / pc_cm < min_winds)
        transformComovingWithStar(BaryonField[DensNum], BaryonField[MetalNum],
                                  BaryonField[MetalIINum], BaryonField[MetalIaNum],
                                  BaryonField[Vel1Num], BaryonField[Vel2Num], BaryonField[Vel3Num],
                                  BaryonField[TENum], BaryonField[GENum],
                                  *up, *vp, *wp, GridDimension[0], GridDimension[1],
                                  GridDimension[2], 1);
    /* Use averaged quantities across multiple cells so that deposition is stable.
        vmean is used to determine whether the supernova shell calculation should proceed:
            M_shell > 0 iff v_shell > v_gas */
    float zmean = 0, dmean = 0, nmean = 0, vmean = 0, mu_mean = 0;
    for (int idi = ip-1; idi < ip+2; idi++)
        for (int idj = jp-1;idj < jp+2; idj++)
            for (int idk = kp-1; idk < kp+2; idk++)
                {
                    int ind = idi + idj * GridDimension[0] + idk * GridDimension[0] * GridDimension[1];

                    
                    zmean += metals[ind];
                    if (SNColourNum > 0)
                        zmean += BaryonField[SNColour][ind];
                    mu_mean += muField[ind];
                    dmean += BaryonField[DensNum][ind];
                    // velocities are currently momenta
                    vmean += BaryonField[Vel1Num][ind]/BaryonField[Vel1Num][ind];
                    vmean += BaryonField[Vel2Num][ind]/BaryonField[Vel2Num][ind];
                    vmean += BaryonField[Vel3Num][ind]/BaryonField[Vel3Num][ind];
                }
    vmean = sqrt(vmean) / 27.0 * VelocityUnits; // cm/s!
    zmean = zmean / SolarMetalFractionByMass / (27.0);
    mu_mean /= 27.0;
    dmean = dmean * DensityUnits / (27.);
    // zmean = zmean * DensityUnits / dmean;
    nmean = dmean / (mh/mu_mean);
    //nmean = max(nmean, 1e-1);
    if (printout) printf ("Zmean = %e Dmean = %e (%e) mu_mean = %e ", zmean, dmean, dmean / DensityUnits, mu_mean);
    if (printout) printf ("Nmean = %f vmean = %f\n", nmean, vmean/1e5);
    float zZsun = max(zmean, 1e-8);
    float fz = (zZsun < 0.01) ? (2.0) : (pow(zZsun, -0.14));

    /* Cooling radius as in Hopkins, but as an average over cells */

    float CoolingRadius = 28.4 *
                          pow(max(0.1, nmean), -3.0 / 7.0) * pow(ejectaEnergy / 1.0e51, 2.0 / 7.0) * fz;
    if (printout)
        fprintf(stdout, "cooling radius [pc] = %e\n %f %e %f %e %e \n",
                CoolingRadius, nmean, ejectaEnergy / 1e51, fz, zmean, dmean);

    float coupledEnergy = ejectaEnergy;

    if (printout)
        fprintf(stdout, "Dx [pc] = %f\n", dx * LengthUnits / pc_cm);
    // although the cell width is 'dx', we actually couple to larger radii than that
    // because of CIC, roughly 2.5-3 *dx
    float cellwidth = dx * LengthUnits / pc_cm;

    float dxRatio = cellwidth / CoolingRadius;
    /* We want to couple one of four phases: free expansion, Sedov-taylor, shell formation, or terminal 
    The first three phases are take forms from Kim & Ostriker 2015, the last from Cioffi 1988*/


    // if we resolve free expansion, all energy is thermally coupled

    float p_free = sqrt(ejectaMass*SolarMass*ejectaEnergy)/SolarMass/1e5;//1.73e4*sqrt(ejectaMass*ejectaEnergy/1e51/3.); // free exp. momentum eq 15
    float r_free = 2.75 * pow(ejectaMass / 3 / nmean, 1. / 3.); // free exp radius eq 2
    bool use_free = false; // could just deposit free expansion into host cell, really...
    // assuming r_sedov == dx, solve for t3

    float t3_sedov = pow(max(r_free, cellwidth) * pc_cm / (5.0 * pc_cm * pow(ejectaEnergy / 1e51 / nmean, 1.0 / 5.0)), 5. / 2.);
    float p_sedov = 2.21e4 * pow(ejectaEnergy / 1e51, 4. / 5.) * pow(nmean, 1. / 5.) * pow(t3_sedov, 3. / 5.); // eq 16

    // shell formation radius eq 8
    float r_shellform = 22.6 * pow(ejectaEnergy / 1e51, 0.29) * pow(nmean, -0.42);
    // p_sf = m_sf*v_sf eq 9,11
    float p_shellform = 3.1e5 * pow(ejectaEnergy / 1e51, 0.94) * pow(nmean, -0.13); // p_sf = m_sf*v_sf eq 9,11

    /* termninal momentum */
    float pTerminal = 4.8e5 * pow(nmean, -1.0 / 7.0) * pow(ejectaEnergy / 1e51, 13.0 / 14.0) * fz; // cioffi 1988, as written in Hopkins 2018

    /* fading radius of a supernova, using gas energy of the host cell and ideal gas approximations */
    float T = BaryonField[GENum][index] / BaryonField[DensNum][index] * TemperatureUnits;
    float cSound = sqrt(kboltz * T / mh) / 1e5; // [km/s] 
    float r_fade = max(66.0*pow(ejectaEnergy/1e51, 0.32)*pow(nmean, -0.37)*pow(cSound/10, -2.0/5.0), dxRatio);
    if (r_fade < CoolingRadius)
        r_fade = CoolingRadius * 1.1;
    float fadeRatio = cellwidth/r_fade;
    if (printout) fprintf(stdout, "Fading: T = %e; Cs = %e; R_f = %e; fadeR = %f\n", T, cSound, r_fade, fadeRatio);

    float coupledMomenta = 0.0;
    float eKinetic = 0.0;
    if (printout)
        fprintf(stdout, "RADII: cell = %e, free = %e, shellform = %e, cooling = %e, fade = %e t_3=%e\n", 
                                        cellwidth, r_free, r_shellform, CoolingRadius, r_fade, t3_sedov);

    if (!winds)
    {  // this calculation for SNe only
        float cw_eff = sqrt(2) * cellwidth; // effective cell width couples to farther than just dx
        float dxeff = cw_eff / CoolingRadius;
        float fader = cw_eff / r_fade;
        if (cw_eff < r_free){
            printf("Coupling free phase\n");
            coupledMomenta = p_free * (1 + pow(cellwidth/r_free, 3));
        }
        if (r_free < cw_eff && fader < 1)
            if (p_sedov < pTerminal && dxeff < 1){
                coupledMomenta = p_sedov + pow(dxeff,1)*(pTerminal-p_sedov);
                printf("Coupling Sedov-Terminal phase\n");
            } else {   
                coupledMomenta = pTerminal / sqrt(1+dxeff);
                printf("Coupling Terminal phase\n");
            }
        if (fader > 1){
            printf("Coupling Fading phase\n");
            float red_fact = max(sqrt(3), sqrt(fader));
            if (fader > 4)
                red_fact = 2.0;
            coupledMomenta = pTerminal / red_fact;
        }
        if ((nmean * T > 1e6 && nmean <= 1) || (nmean <=0.01 && T > 1e5)){ // in high-pressure, low nb, p_t doesnt hold since there is essentailly no radiative phase.
                                        // I cannot find a good analytic expression to use, but since there is no snowplough, swept up
                                        // thermal energy dominates the evolution (Tang, 2005, doi 10.1086/430875 )
                                        // for now, couple the free expansion, until a better expression can be found.
            printf("Coupling high-pressure low-n phase (free expansion)\n");
            coupledMomenta = p_free;
        }
    }

    if (winds)
    { // simple divide momenta and thermal energy
        
        coupledMomenta = sqrt(ejectaMass*SolarMass* 0.5 * ejectaEnergy)/SolarMass/1e5;

    }


    if (printout)
        fprintf(stdout, "Calculated p = %e (sq_fact = %e; p_f = %e; p_t = %e; mcell = %e; mcpl = %e)\n", 
                                coupledMomenta, (dmean / DensityUnits * MassUnits) / ejectaMass * ntouched, p_free, pTerminal, dmean / DensityUnits * MassUnits, ejectaMass/27.0);


    //    coupledMomenta = (cellwidth > r_fade)?(coupledMomenta*pow(r_fade/cellwidth,3/2)):(coupledMomenta);
    float shellMass = 0.0, shellVelocity = 0.0;
    if (printout)
        fprintf(stdout, "Coupled momentum: %e\n", coupledMomenta);
    /* 
        If resolution is in a range comparable to Rcool and
        Analytic SNR shell mass is on, adjust the shell mass 
        upper range of applicability for shell mass is determined by
        local average gas velocity (v_shell ~ v_gas = no shell)
    */
    if (cellwidth > r_shellform && coupledEnergy > 0 && AnalyticSNRShellMass)
    {
        shellVelocity = 413.0 * pow(nmean, 1.0 / 7.0) * pow(zZsun, 3.0 / 14.0) * pow(coupledEnergy / EnergyUnits / 1e51, 1.0 / 14.0) * pow(dxRatio, -7.0 / 3.0); //km/s

        if (shellVelocity > vmean / 1e5)
        {
            shellMass = max(8e3, coupledMomenta / shellVelocity); //Msun

            /* cant let host cells evacuate completely!  
                   7.974045e+017.974045e+017.974045e+017.974045e+017.974045e+01 Shell mass will be evacuated from central cells by CIC a negative mass,
                    so have to check that the neighbors can handle it too*/
            for (int ind = -1; ind <= 1; ind++)
            {
                float minD = min(BaryonField[DensNum][index + ind], BaryonField[DensNum][index + GridDimension[0] * ind]);
                minD = min(minD, BaryonField[DensNum][index + GridDimension[0] * GridDimension[1] * ind]);
                if (shellMass >= 0.25 * minD * MassUnits)
                    shellMass = 0.25 * minD * MassUnits;
            }
        }
    }

    if (shellMass < 0.0)
    {
        fprintf(stdout, "Shell mass = %e Velocity= %e P = %e",
                shellMass, shellVelocity, coupledMomenta);
        ENZO_FAIL("SM_deposit: 391");
    }
    float coupledMass = shellMass + ejectaMass;
    // kinetic energy from the momenta, taking mass as ejecta + mass of effected cells (4^3 because of CIC in coupling cloud)
    eKinetic = coupledMomenta * coupledMomenta / (2.0 *(ejectaMass + ntouched * dmean * pow(LengthUnits * CellWidth[0][0], 3) / SolarMass)) * SolarMass * 1e10;
    if (eKinetic > 1e53){
        fprintf(stdout, "Rescaling high kinetic energy %e -> ", eKinetic);
        coupledMomenta = sqrt((2.0 * (ejectaMass*SolarMass + ntouched * dmean * pow(LengthUnits*dx, 3) ) * ejectaEnergy))/SolarMass/1e5;
        eKinetic = coupledMomenta * coupledMomenta / (2.0 *(ejectaMass + ntouched * dmean * pow(LengthUnits * CellWidth[0][0], 3) / SolarMass)) * SolarMass * 1e10;
        
        fprintf(stdout, " %e; new p = %e\n", eKinetic, coupledMomenta);
    }
    
    if (printout)
        fprintf(stdout, "Ekinetic = %e Mass = %e\n",
                eKinetic, dmean * pow(LengthUnits * CellWidth[0][0], 3) / SolarMass);
    if (eKinetic > 1e60 && winds)
    {
        fprintf(stdout, "winds Ekinetic = %e Mass = %e\n",
                eKinetic, dmean * pow(LengthUnits * CellWidth[0][0], 3) / SolarMass);
        ENZO_FAIL("winds Ekinetic > reasonability!\n");
    }
    if (eKinetic > 1e60 && !winds)
    {
        fprintf(stdout, "Ekinetic = %e Mass = %e\n",
                eKinetic, dmean * pow(LengthUnits * CellWidth[0][0], 3) / SolarMass);
        ENZO_FAIL("SNE Ekinetic > reasonability!\n");
    }

    float coupledGasEnergy = max(ejectaEnergy - eKinetic, 0);
    if (printout)
        fprintf(stdout, "Coupled Gas Energy = %e\n", coupledGasEnergy);
    if (dxRatio > 1.0 && !winds) // if we apply this reduction to winds, then there is literally *no* effect, even at Renaissance resolution.
        coupledGasEnergy = (DepositUnresolvedEnergyAsThermal)
                               ? (coupledGasEnergy)
                               : (coupledGasEnergy * pow(dxRatio, -6.5));

    float shellMetals = zZsun * 0.02 * shellMass;
    float coupledMetals = 0.0, SNIAmetals = 0.0, SNIImetals = 0.0, P3metals = 0.0;
    if (winds)
        coupledMetals = ejectaMetal; //+ shellMetals; // winds only couple to metallicity
    if (AnalyticSNRShellMass)
        coupledMetals += shellMetals;
    SNIAmetals = (StarMakerTypeIaSNe) ? nSNIA * 1.4 : 0.0;
    if (!StarMakerTypeIaSNe)
        coupledMetals += nSNIA * 1.4;
    SNIImetals = (StarMakerTypeIISNeMetalField) ? nSNII * (1.91 + 0.0479 * max(starMetal, 1.65)) : 0.0;
    if (!StarMakerTypeIISNeMetalField)
        coupledMetals += nSNII * (1.91 + 0.0479 * max(starMetal, 1.65));
    if (isP3 && MechStarsSeedField)
        P3metals = ejectaMetal;

    if (printout)
        fprintf(stdout, "Coupled Metals: %e %e %e %e %e %e\n", ejectaMetal, SNIAmetals, SNIImetals, shellMetals, P3metals, coupledMetals);

    /* 
        Critical debug compares the pre-feedback and post-feedback field sums.  Note that this doesn't 
        work well on vector quantities outside of an ideal test.
    */
    float preMass = 0, preZ = 0, preP = 0, prePmag = 0, preTE = 0, preGE = 0, preZII = 0, preZIa = 0;
    float dsum = 0.0, zsum = 0.0, psum = 0.0, psqsum = 0.0, tesum = 0.0, gesum = 0.0, kesum = 0.0;
    float postMass = 0, postZ = 0, postP = 0, postPmag = 0, postTE = 0, postGE = 0, postZII = 0, postZIa = 0;
    if (criticalDebug)
    {
        // for (int i = 0; i < size; ++i)
        for (int k = max(0,kp-5); k <= min(kp+5, GridDimension[2]); ++k)
            for (int j = max(0,jp-5); j <= max(jp+5, GridDimension[1]); ++j)
                for (int i = max(0,ip-5); i <= max(ip+5, GridDimension[0]); ++i)

                {
                    int idx = i + j * GridDimension[0] + k * GridDimension[0] * GridDimension[1];
                    preMass += BaryonField[DensNum][idx];
                    preZ += BaryonField[MetalNum][idx];
                    if (StarMakerTypeIISNeMetalField)
                        preZII += BaryonField[MetalIINum][idx];
                    if (StarMakerTypeIaSNe)
                        preZIa += BaryonField[MetalIaNum][idx];
                    if (MechStarsSeedField)
                        preZ += BaryonField[SNColourNum][idx];
                    preP += BaryonField[Vel1Num][idx] + BaryonField[Vel2Num][idx] + BaryonField[Vel3Num][idx];
                    prePmag += pow(BaryonField[Vel1Num][idx], 2) +
                            pow(BaryonField[Vel2Num][idx], 2) + pow(BaryonField[Vel3Num][idx], 2);
                    preTE += BaryonField[TENum][idx];
                    preGE += BaryonField[GENum][idx];
                }
    }
    /* Reduce coupled quantities to per-particle quantity and converMetalNumt to 
        code units.
        Hopkins has complicated weights due to complicated geometry. 
            This implementation is simple since our coupled particles are 
            spherically symmetric about the feedback particle*/

    if (!winds)
        coupledEnergy = min((nSNII + nSNIA) * 1e51, eKinetic);
    if (printout) fprintf(stdout, "Pre: Coupling TE = %e, GE = %e\n", coupledEnergy, coupledGasEnergy);
    
    coupledEnergy = coupledEnergy / EnergyUnits;
    coupledGasEnergy = coupledGasEnergy / EnergyUnits;
    if (printout) fprintf(stdout, "Post: Coupling TE = %e, GE = %e\n", coupledEnergy, coupledGasEnergy);
    if (printout) fprintf(stdout, "Coupled mass: %e Msun", coupledMass);
    coupledMass /= MassUnits;
    if (printout) fprintf(stdout, "Coupled Mass: %e code_mass", coupledMass);
    coupledMetals /= MassUnits;

    // conversion includes km -> cm
    if (printout) fprintf(stdout, "475: Coupling %e Msun km/s\n", coupledMomenta);
    coupledMomenta = coupledMomenta / MomentaUnits;
    if (printout) fprintf(stdout, "477: %e mass * len / time (%e = %e * %e / %e) (vu = %e)\n", coupledMomenta, MomentaUnits, MassUnits, LengthUnits, TimeUnits, VelocityUnits);
    SNIAmetals /= MassUnits;
    SNIImetals /= MassUnits;
    P3metals /= MassUnits;

    /* deposit the particles with their respective quantities */

    FLOAT LeftEdge[3] = {CellLeftEdge[0][0], CellLeftEdge[1][0], CellLeftEdge[2][0]};
    FLOAT pX, pY, pZ;
    float eCouple, geCouple, mCouple, zCouple, zIICouple, zIACouple, p3Couple;
    float expect_momenta = 0;
    float expect_error []= {0,0,0};
    float expect_z = 0;
    int n_pos[] = {0,0,0};
// printf("Coupling for star particle at %e %e \n", *xp, *yp, *zp);
/* LOOP+SETUP FOR CIC DEPOSIT */
        for (int i = -1; i <= 1; i++)
            for (int j = -1; j <= 1; j++)
                for (int k = -1; k <= 1; k++)
                    {
                        
                        // nearest index of coupling particle...
                        
                        float xc, yc, zc;
                        xc = *xp + i * dx;
                        yc = *yp + j * dx;
                        zc = *zp + k * dx;
                        int ic, jc, kc;
                        ic = (xc - CellLeftEdge[0][0] - 0.5 * dx) / dx;
                        jc = (yc - CellLeftEdge[1][0] - 0.5 * dx) / dx;
                        kc = (zc - CellLeftEdge[2][0] - 0.5 * dx) / dx;

                        // printf("Coupling particle near %e %e %e => %d %d %d \n\t\t(%e %e %e), dx = %e\n", xc, yc, zc, 
                        //                                                                 ic, jc, kc,
                        //                                                                 (*xp - CellLeftEdge[0][0] - 0.5 * dx) / dx,
                        //                                                                 (*yp - CellLeftEdge[1][0] - 0.5 * dx) / dx,
                        //                                                                 (*zp - CellLeftEdge[2][0] - 0.5 * dx) / dx,
                        //                                                                 dx);
                        if  ( i == 0 && j == 0 && k == 0) continue;

                        float modi, modj, modk;
                        modi = (float) i;
                        modj = (float) j;
                        modk = (float) k;
                        float uniti, unitj, unitk, normfact;
                        normfact = sqrt(modi*modi + modj*modj + modk*modk);
                        uniti = modi / normfact;
                        unitj = modj / normfact;
                        unitk = modk / normfact;

                        /************************/

                   
                        
                        float sq3 = sqrt((float) nCouple);

                        pX = coupledMomenta * (float) uniti / sq3; //* weightsVector[n];
                        pY = coupledMomenta * (float) unitj / sq3; //* weightsVector[n];
                        pZ = coupledMomenta * (float) unitk / sq3; //* weightsVector[n];

                        eCouple = coupledEnergy / (float) nCouple; //* weightsVector[n];
                        geCouple = coupledGasEnergy / (float) nCouple; //* weightsVector[n];
                        mCouple = coupledMass / (float) nCouple; //* weightsVector[n];
                        zCouple = coupledMetals / (float) nCouple; //* weightsVector[n];
                        if (StarMakerTypeIISNeMetalField)
                            zIICouple = SNIImetals / (float) nCouple; //* weightsVector[n];
                        if (StarMakerTypeIaSNe)
                            zIACouple = SNIAmetals / (float) nCouple; //* weightsVector[n];
                        if (MechStarsSeedField)
                            p3Couple = P3metals / (float) nCouple; //* weightsVector[n];
                        int np = 1;


                        expect_error[0] += pX;
                        expect_error[1] += pY;
                        expect_error[2] += pZ;

             
/*
    THIS SECTION DOES A CIC-DEPOSIT USING THE BUILT IN ROUTINES--IT
    CIC DEPOSITS EACH COUPLING PARTICLE
*/

/* Hand rolled CIC for some stupid reason... */
                    float dep_vol_frac = 0;

                    
                    for (int kk = -1; kk <= 1; kk ++)
                        for (int jj = -1; jj <= 1; jj++)
                            for (int ii = -1; ii <= 1; ii++)
                            {
                                /* position of cloud cell */
                                float xcell, ycell, zcell;
                                xcell = xc + dx * ii;
                                ycell = yc + dx * jj;
                                zcell = zc + dx * kk;

                                /* index of cloud cell */
                                int icell, jcell, kcell;
                                icell = (xcell - CellLeftEdge[0][0] - 0.5 * dx) / dx;
                                jcell = (ycell - CellLeftEdge[1][0] - 0.5 * dx) / dx;
                                kcell = (zcell - CellLeftEdge[2][0] - 0.5 * dx) / dx;

                                /* position of cell to couple to, center value */
                                xcell = CellLeftEdge[0][0] + (0.5 + (float) icell) * dx;
                                ycell = CellLeftEdge[1][0] + (0.5 + (float) jcell) * dx;
                                zcell = CellLeftEdge[2][0] + (0.5 + (float) kcell) * dx;
                                /* flat index ... */
                                int flat = icell + jcell * GridDimension[0] + kcell * GridDimension[0] * GridDimension[1];

                                

                                /* fraction of quantities that go into this cloud cell... */
                                float window = Window(xc - xcell, yc - ycell, zc - zcell, dx, NGP);
                                dep_vol_frac += window;
                                // printf("\t\t\tCloud %e %e %e Coupling at %d %d %d with window %e.  \n\t\t\t\twfactors: %e %e %e (%f %f %f)\n", 
                                //                                             xcell, ycell, zcell, 
                                //                                             icell, jcell, kcell, window,
                                //                                             xc-xcell, yc-ycell, zc-zcell,
                                //                                             (xc-xcell)/dx, (yc-ycell)/dx, (zc-zcell)/dx);
                                /* add quantities to this cell accordingly */

                                BaryonField[DensNum][flat] = BaryonField[DensNum][flat] + mCouple * window;
                                BaryonField[MetalNum][flat] = BaryonField[MetalNum][flat] + zCouple * window;
                                if (eCouple > 0 && DualEnergyFormalism)
                                    {    
                                        BaryonField[TENum][flat] = BaryonField[TENum][flat] + eCouple * window;
                                        BaryonField[TENum][flat] = BaryonField[TENum][flat] + geCouple * window;
                                    }
                                if (geCouple > 0)
                                    BaryonField[GENum][flat] = BaryonField[GENum][flat] + geCouple * window;
                                float prep = BaryonField[Vel1Num][flat]*BaryonField[Vel1Num][flat]
                                                 + BaryonField[Vel2Num][flat]*BaryonField[Vel2Num][flat] 
                                                 + BaryonField[Vel3Num][flat]*BaryonField[Vel1Num][flat];
                                BaryonField[Vel1Num][flat] = BaryonField[Vel1Num][flat] + pX * window;
                                BaryonField[Vel2Num][flat] = BaryonField[Vel2Num][flat] + pY * window;
                                BaryonField[Vel3Num][flat] = BaryonField[Vel3Num][flat] + pZ * window;
                                expect_momenta += BaryonField[Vel1Num][flat]*BaryonField[Vel1Num][flat]
                                                 + BaryonField[Vel2Num][flat]*BaryonField[Vel2Num][flat] 
                                                 + BaryonField[Vel3Num][flat]*BaryonField[Vel1Num][flat]
                                                 - prep;

                                // TODO: add other metal fields for consistency.

                            }
                        // printf("CIC deposited in sum window %e\n", dep_vol_frac);
                    }


            
    if (printout){
        fprintf(stdout, "After deposition, counted %e Msun km/s momenta deposited.  Error = %e..\n", sqrt(expect_momenta) * MomentaUnits, 
                                                                                        (expect_error[0]+expect_error[1]+expect_error[2]) * MomentaUnits);
        fprintf(stdout, "\tTabulated %e Msun deposited metals\n", expect_z * MassUnits);
    }
    // printf("\n");
    /* Deposit one negative mass particle centered on star to account for 
        shell mass leaving host cells .  Same for metals that were evacuated*/
    // int np = 1;
    // shellMass *= -1 / MassUnits;
    // FORTRAN_NAME(cic_deposit)
    //     (xp, yp, zp, &GridRank, &np, &shellMass, density, LeftEdge,
    //     &GridDimension[0], &GridDimension[1], &GridDimension[2], &dx, &cloudSize);
    // shellMetals *= -1 / MassUnits;
    // FORTRAN_NAME(cic_deposit)
    //     (xp, yp, zp, &GridRank, &np, &shellMetals, metals, LeftEdge,
    //     &GridDimension[0], &GridDimension[1], &GridDimension[2], &dx, &cloudSize);


    if (criticalDebug)
    {
        for (int k = max(0,kp-5); k <= min(kp+5, GridDimension[2]); ++k)
            for (int j = max(0,jp-5); j <= max(jp+5, GridDimension[1]); ++j)
                for (int i = max(0,ip-5); i <= max(ip+5, GridDimension[0]); ++i)
                {
                    int idx = i + j * GridDimension[0] + k * GridDimension[0] * GridDimension[1];            
                    postMass += BaryonField[DensNum][idx];
                    postZ += BaryonField[MetalNum][idx];
                    if (StarMakerTypeIISNeMetalField)
                        postZII += BaryonField[MetalIINum][idx];
                    if (StarMakerTypeIaSNe)
                        postZIa += BaryonField[MetalIaNum][idx];
                    if (SNColourNum > 0)
                        postZ += BaryonField[SNColourNum][idx];
                    postP += BaryonField[Vel1Num][idx] + BaryonField[Vel2Num][idx] + BaryonField[Vel3Num][idx];
                    postPmag += pow(BaryonField[Vel1Num][idx], 2) +
                                pow(BaryonField[Vel2Num][idx], 2) + pow(BaryonField[Vel3Num][idx], 2);
                    postTE += BaryonField[TENum][idx];
                    postGE += BaryonField[GENum][idx];
                }

        if (printout)
            fprintf(stderr, "Difference quantities: dxRatio = %f dMass = %e dZ = %e dzII = %e dxIa = %e  P = %e |P| = %e (%e) TE = %e GE = %e coupledGE = %e Ej = %e Mej = %e Zej = %e\n",
                    dxRatio, (postMass - preMass) * MassUnits, (postZ - preZ) * MassUnits,
                    (postZII - preZII) * MassUnits, (postZIa - preZIa) * MassUnits,
                    (postP - preP) * MomentaUnits,
                    (sqrt(postPmag) - sqrt(prePmag))*MomentaUnits, sqrt(expect_momenta) * MomentaUnits,
                    (postTE - preTE) * EnergyUnits, (postGE - preGE) * EnergyUnits,
                    coupledGasEnergy * EnergyUnits, ejectaEnergy,
                    ejectaMass, ejectaMetal);
        if (isnan(postMass) || isnan(postTE) || isnan(postPmag) || isnan(postZ))
        {
            fprintf(stderr, "NAN IN GRID:  %f %e %e %e-%e %e\n", postMass, postTE, postZ, preZ, postP);
            for (float w : weightsVector)
            {
                fprintf(stderr, "%e\t", w);
            }
            fprintf(stderr, "\n");
            exit(3);
            ENZO_FAIL("MechStars_depositFeedback.C: 530\n")
        }
    }

    transformComovingWithStar(BaryonField[DensNum], BaryonField[MetalNum],
                            BaryonField[MetalIINum], BaryonField[MetalIaNum],
                            BaryonField[Vel1Num], BaryonField[Vel2Num],
                            BaryonField[Vel3Num],
                            BaryonField[TENum], BaryonField[GENum],
                            *up, *vp, *wp,
                            GridDimension[0], GridDimension[1],
                            GridDimension[2], -1);

    return SUCCESS;
}
float Window(float xd, float yd, float zd, float width, bool NGP){


    float wx = 0, wy = 0, wz = 0;
    if (!NGP)
    {
        if (fabs(xd) <= width)
                wx = 1.0 - fabs(xd)/width;
        if (fabs(yd) <= width)
                wy = 1.0 - fabs(yd)/width;
        if (fabs(zd) <= width)
                wz = 1.0 - fabs(zd)/width;
    }
    if (NGP) 
    {
        if (fabs(xd) < width/2.0)
            wx = 1.0;
        if (fabs(yd) < width/2.0)
            wy = 1.0;
        if (fabs(zd) < width/2.0)
            wz = 1.0;
    }
    return wx * wy * wz;
}