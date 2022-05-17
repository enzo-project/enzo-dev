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
    //printf("STARSS_FB: In Feedback deposition\n");
    if (MyProcessorNumber != ProcessorNumber)
        return 0;
    bool debug = true;
    bool criticalDebug = false;
    bool printout = debug & !winds;
    int index = ip + jp * GridDimension[0] + kp * GridDimension[0] * GridDimension[1];
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
    bool useFading = MechStarsFadeSNR; // true to fade the SN--it only couples heat after a certain radii that exceeds the expected merging radius of the SN
    float ntouched = NGP ? (26) : (63); // how many cells get touched by deposition? 
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
        fprintf(stdout, "STARSS_FB: depositing quantities: Energy %e, Mass %e, Metals %e\n",
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
        fprintf(stdout, "STARSS_FB: MechStars only functions with metallicity field enabled!");
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
    std::vector <float> prior_ke (125, 0.0);

            for (int i = -2; i <= 2; i++)
            for (int j = -2; j <= 2; j++)
                for (int k = -2; k <= 2; k++)
                    {
                        int flat = ip+i + (jp+j) * GridDimension[0] + (kp+k) * GridDimension[0] * GridDimension[1];
                        int bind = i+2 + (j+2) * 5 + (k+2)*25;

                        prior_ke[bind] += 0.5 * BaryonField[DensNum][flat] * 
                                        (BaryonField[Vel1Num][flat]*BaryonField[Vel1Num][flat]
                                        +BaryonField[Vel2Num][flat]*BaryonField[Vel2Num][flat]
                                        +BaryonField[Vel3Num][flat]*BaryonField[Vel3Num][flat]);

                    }
        transformComovingWithStar(BaryonField[DensNum], BaryonField[MetalNum],
                                  BaryonField[MetalIINum], BaryonField[MetalIaNum],
                                  BaryonField[Vel1Num], BaryonField[Vel2Num], BaryonField[Vel3Num],
                                  BaryonField[TENum], BaryonField[GENum],
                                  *up, *vp, *wp, GridDimension[0], GridDimension[1],
                                  GridDimension[2], 1);
    /* Use averaged quantities across multiple cells so that deposition is stable.
        vmean is used to determine whether the supernova shell calculation should proceed:
            M_shell > 0 iff v_shell > v_gas */
    float zmean = 0, dmean = 0, nmean = 0, mu_mean = 0;
    for (int idi = ip-1; idi < ip+2; idi++)
        for (int idj = jp-1;idj < jp+2; idj++)
            for (int idk = kp-1; idk < kp+2; idk++)
                {
                    int ind = idi + idj * GridDimension[0] + idk * GridDimension[0] * GridDimension[1];

                    float tz = 0;
                    tz  = metals[ind];
                    if (SNColourNum > 0)
                        tz += BaryonField[SNColour][ind];
                    zmean += tz/BaryonField[DensNum][ind];
                    mu_mean += muField[ind];
                    dmean += BaryonField[DensNum][ind];
                }
    zmean = zmean / SolarMetalFractionByMass / (27.0);
    mu_mean /= 27.0;
    dmean = dmean * DensityUnits / (27.);
    // zmean = zmean * DensityUnits / dmean;
    nmean = dmean / (mh/mu_mean);
    // nmean = BaryonField[DensNum][index]*DensityUnits / (mh/muField[index]);
    //nmean = max(nmean, 1e-1);
    if (printout) printf ("STARSS_FB: Zmean = %e Dmean = %e (%e) mu_mean = %e ", zmean, dmean, dmean / DensityUnits, mu_mean);
    if (printout) printf ("Nmean = %f\n", nmean);
    float zZsun = zmean;
    float fz = min(2.0, pow(max(zZsun, 0.01), -0.14));

    /* Cooling radius as in Hopkins, but as an average over cells */

    float CoolingRadius = 28.4 *
                          pow(max(0.01, nmean), -3.0 / 7.0) * pow(ejectaEnergy / 1.0e51, 2.0 / 7.0) * fz;

    float coupledEnergy = ejectaEnergy;

    float cellwidth = dx * LengthUnits / pc_cm;

    float dxRatio = cellwidth / CoolingRadius;
    /* We want to couple one of four phases: free expansion, Sedov-taylor, shell formation, or terminal 
    The first three phases are take forms from Kim & Ostriker 2015, the last from Cioffi 1988*/


    // if we resolve free expansion, all energy is thermally coupled

    float p_free = sqrt(2.0 * ejectaMass*SolarMass*ejectaEnergy)/SolarMass/1e5;//1.73e4*sqrt(ejectaMass*ejectaEnergy/1e51/3.); // free exp. momentum eq 15
    float r_free = 2.75 * pow(ejectaMass / 3 / nmean, 1. / 3.); // free exp radius eq 2
    bool use_free = false; // could just deposit free expansion into host cell, really...
    // assuming r_sedov == dx, solve for t3

    float t3_sedov = pow(max(r_free, max(cellwidth, r_free)) * pc_cm / (5.0 * pc_cm * pow(ejectaEnergy / 1e51 / nmean, 1.0 / 5.0)), 5. / 2.);
    float p_sedov = 2.21e4 * pow(ejectaEnergy / 1e51, 4. / 5.) * pow(nmean, 1. / 5.) * pow(t3_sedov, 3. / 5.); // eq 16

    // shell formation radius eq 8
    float r_shellform = 22.6 * pow(ejectaEnergy / 1e51, 0.29) * pow(nmean, -0.42);
    // p_sf = m_sf*v_sf eq 9,11
    float p_shellform = 3.1e5 * pow(ejectaEnergy / 1e51, 0.94) * pow(nmean, -0.13); // p_sf = m_sf*v_sf eq 9,11

    /* 
        termninal momentum 
            We picked the analytic form from Thornton as it matched high-resolution SN we ran the best.
    */
    float pTerminal;
    if (zZsun > 0.01)
    // Cioffi:
        // pTerminal = 4.8e5 * pow(nmean, -1.0 / 7.0) * pow(ejectaEnergy / 1e51, 13.0 / 14.0) * pow(zZsun, -3.0/14.0); // cioffi 1988, as written in Hopkins 2018
    // Thornton, 1998, M_s * V_s, eqs 22, 23, 33, 34
        pTerminal = 1.6272e5 * pow(ejectaEnergy/1e51, 13./14.) * pow(nmean, -0.25) * pow(zZsun, -0.36);
    // kimm & ostriker 2015 
        // pTerminal = 2.8e5 * pow(ejectaEnergy/1e51, 13./14.) * pow(nmean, -0.17) * pow(zZsun, -0.36);
    else
        pTerminal = 8.3619e5 * pow(ejectaEnergy/1e51, 13./14.) * pow(nmean, -0.25);

    /* fading radius of a supernova, using gas energy of the host cell and ideal gas approximations */
    float T = max(100, BaryonField[GENum][index] * (Gamma-1) * muField[index] * TemperatureUnits);
    if ( BaryonField[GENum][index] * (Gamma-1) * muField[index] * TemperatureUnits < 0 ){
        fprintf(stdout, "Error: Negative temperature encountered. GE = %f; Gamma = %f; mu = %f", 
            BaryonField[GENum][index], Gamma, muField[index]);
    }
    float cSound = sqrt(Gamma * kboltz * T / (mh*mu_mean)) / 1e5; // [km/s] 
    float r_fade = max(66.0*pow(ejectaEnergy/1e51, 0.32)*pow(nmean, -0.37)*pow(cSound/10, -2.0/5.0), CoolingRadius * 1.5);
    float fadeRatio = cellwidth/r_fade;
    if (printout) fprintf(stdout, "STARSS_FB: Fading: T = %e; Cs = %e; R_f = %e; fadeR = %f\n", T, cSound, r_fade, fadeRatio);

    float coupledMomenta = 0.0;
    float eKinetic = 0.0;


    float cw_eff = NGP ? cellwidth : sqrt(1)*cellwidth; // effective cell width couples to farther than just dx.  
                                        // theres a lot of numerical fudge factors here because of that.
                                        // the actual coupling is between 2-3 dx, depending on position within the cell 

    float dxeff = cw_eff / CoolingRadius;
    float fader = cw_eff / r_fade;

    // velocity of shell ala Thornton 1998
    float shellVelocity = 413.0 * pow(nmean, 1.0 / 7.0) * pow(zZsun, 3.0 / 14.0) * pow(coupledEnergy / 1e51, 1.0 / 14.0);
    float ratio_pds = cw_eff/r_shellform;
    shellVelocity *=  ratio_pds > 1 ? pow(dxeff, -7.0 / 3.0) : 1; //km/s
    float beta =  max( 1.0, shellVelocity / max(1, cSound)); 
    float nCritical = 0.0038* (pow(nmean * T / 1e4, 7.0/9.0) * pow(beta, 14.0/9.0))/(pow(ejectaEnergy/1e51, 1.0/9.0) * pow(max(0.0001, zZsun), 1.0/3.0)); // n/cc
    float rmerge = max(151.0 * pow((ejectaEnergy/1e51)/ beta / beta / nmean / T * 1e4, 1./3.), 1.25 * r_fade); // pc
    float merger = cw_eff / rmerge;
    bool faded = fader > 1;
    if (printout)
        fprintf(stdout, "STARSS_FB: RADII: cell = %e, free = %e, shellform = %e, cooling = %e, fade = %e t_3=%e, R_m=%e\n", 
                                        cellwidth, r_free, r_shellform, CoolingRadius, r_fade, t3_sedov, rmerge);
    if (!winds)
    {  // this calculation for SNe only
        // coupledMomenta = p_free * min(sqrt(1+ (nCouple * dmean * pow(cellwidth * pc_cm, 3) / SolarMass)/(ejectaMass)), pTerminal/p_free/pow(1+dxeff));
        if (cw_eff < r_free){
            coupledMomenta = min(p_free * pow(cw_eff/r_free, 3.0), p_sedov);
            printf("STARSS_FB: modifying free phase: p = %e\n", coupledMomenta);
        }
        if (r_free < cw_eff && dxeff <= 1){
                coupledMomenta = min(p_sedov, pTerminal*dxeff);
                printf("STARSS_FB: Coupling Sedov-Terminal phase: p = %e (ps = %e, pt = %e, dxe = %e)\n", coupledMomenta, p_sedov, pTerminal, dxeff);
        }
        if (dxeff > 1){   
                coupledMomenta = pTerminal/ sqrt(min(1.5, dxeff));
               if (printout) printf("STARSS_FB: Coupling Terminal phase: p = %e; dxeff = %e\n", coupledMomenta, dxeff);
            }
        if (fader > 1 && useFading){ // high coupling during the fading regime leads to SNRs on the root-grid in 6-level AMR simulations!
            coupledMomenta = pTerminal * (1.0 - tanh(pow(fader * merger, 2.5)));
           if (printout) printf("STARSS_FB: Coupling Fading phase: p = %e\n", coupledMomenta);
        }
        // critical density to skip snowplough (remnant combines with ISM before radiative phase); eq 4.9 cioffi 1988
        if (printout)
           printf("STARSS_FB: Checking critical density metric... (nmean = %e; N_Crit = %e; factors: %e %e %e; beta = %e/%e == %e; rmerge = %e)\n", 
                                        nmean, nCritical, pow(nmean * T / 1e4, 7.0/9.0), pow(ejectaEnergy/1e51, 1.0/9.0), pow(fz, 1.0/3.0), shellVelocity , cSound, beta, rmerge);
        if (nmean <= 10.0 * nCritical){ // in high-pressure, low nb, p_t doesnt hold since there is essentailly no radiative phase.
                                        // thermal energy dominates the evolution (Tang, 2005, doi 10.1086/430875 )
                                        // We inject 100 % thermal energy to simulate this recombining with the ISM
                                        // and rely on the hydro and the thermal radiation to arrive at the right solution
            coupledMomenta = coupledMomenta * (1.0-tanh(pow(1.45*nCritical/nmean, 6.5)));
            // coupledEnergy = coupledEnergy * (1.0-tanh(pow(4.5*nCritical/nmean, 1.5))); // this is just making a lot of hot gas slowing the sim down... =/
           if (printout) printf("STARSS_FB: Adjusting for high-pressure low-n phase (thermal coupling: Nc = %e): p = %e\n", nCritical, coupledMomenta);
        }
            
        if (T > 1e6 && coupledMomenta > 1e5){
            printf("STARSS_FB: Coupling high momenta to very hot gas!! (p= %e, T= %e, n_c = %e)\n", coupledMomenta, T, nCritical);
        }
    }

    if (winds)
    { // simple divide momenta and thermal energy
        
        coupledMomenta = sqrt(ejectaMass*SolarMass* 0.5 * ejectaEnergy)/SolarMass/1e5;

    }


    if (printout) 
        fprintf(stdout, "STARSS_FB: Calculated p = %e (sq_fact = %e; p_f = %e; p_t = %e; mcell = %e; mcpl = %e)\n", 
                                coupledMomenta, (dmean / DensityUnits * MassUnits) / ejectaMass * ntouched, p_free, pTerminal, dmean / DensityUnits * MassUnits, ejectaMass/27.0);


    //    coupledMomenta = (cellwidth > r_fade)?(coupledMomenta*pow(r_fade/cellwidth,3/2)):(coupledMomenta);
    float shellMass = 0.0;
    /* 
        For various stages, we update the shell mass, ie, the mass that 
        would be swept up by the remnant.  We us the shell-mass expected from 
        Thronton 1998 along with approximations for momenta to derive a shell velocity
        as well. 

    */

    float centralMass = 0;
    float centralMetals = 0;
    float maxEvacFraction = 0.75;
    if (coupledEnergy > 0 && AnalyticSNRShellMass && !winds)
    {
            if (dxeff < 1) // goes like a sphere sweeping up density
                shellMass = 4 * M_PI / 3 * dmean * pow((cw_eff*pc_cm),3) / SolarMass; // min(1e8, coupledMomenta / shellVelocity); //Msun
            //only have that expression of velocity for pds stage.  after, we take the M_s from 
            // thornton 1998, eq 22, 33.
            if (dxeff > 1){ 
                if (zZsun > 0.01)
                    shellMass = 1.41e4 *pow(ejectaEnergy/1e51, 6./7.) * pow(dmean, -0.24) * pow(zZsun, -0.27);
                else{
                    shellMass = 4.89e4* pow(ejectaEnergy/1e51, 6./7.) * pow(dmean, -0.24);
                }
            }
            // else
                // shellMass = ejectaMass; // shell mass increases until comparable to ejecta mass, then we enter snowplough phases
            /* cant let host cells evacuate completely!  
                   Shell mass will be evacuated from central cells by CIC a negative mass,
                    so have to check that the neighbors can handle it too*/
        if (shellMass > 0)
            for (int i = ip-1; i <= ip+1; i++)
                for (int j = jp-1; j <= jp+1; j++)
                   for (int k = kp-1; k <= kp+1; k++)
                        {
                        int flat = i + j * GridDimension[0] + k * GridDimension[0] * GridDimension[1];
                        if (flat < size){
                            // only record if this cell would've been touched by CIC on star particle (estimated as the "blast interior")
                            float xcell, ycell, zcell;
                            /* position of cell to couple to, center value */
                            xcell = CellLeftEdge[0][0] + (0.5 + (float) i) * dx;
                            ycell = CellLeftEdge[1][0] + (0.5 + (float) j) * dx;
                            zcell = CellLeftEdge[2][0] + (0.5 + (float) k) * dx;
                            float window = Window(*xp - xcell, *yp - ycell, *zp - zcell, dx, false); // always use cic to remove the mass
                            if (window > 0){
                                centralMass +=  min(StarMakerMassEfficiency,  window) * BaryonField[DensNum][flat];
                                centralMetals +=  window * BaryonField[MetalNum][flat];
                                if (SNColourNum != -1 && !MechStarsMetallicityFloor)
                                    centralMetals += window * BaryonField[SNColourNum][flat];
                            }
                        }
                    }

        centralMass *= MassUnits;
        if (shellMass > maxEvacFraction * centralMass){
            if (printout) 
                fprintf(stdout, "STARSS: Shell mass too high for host cells: Rescaling %e -> %e\n", shellMass, maxEvacFraction * centralMass);
            shellMass = maxEvacFraction * centralMass; 
            }
        }
    // }

    float shellMetals = min(maxEvacFraction*centralMetals * MassUnits, zZsun * SolarMetalFractionByMass * shellMass);
    if (AnalyticSNRShellMass && printout)
    {
        fprintf(stdout, "STARSS_FB: Shell_m = %e Shell_z = %e shell_V= %e P = %e M_C = %e\n",
                shellMass, shellMetals, shellVelocity, coupledMomenta, centralMass);
        // ENZO_FAIL("SM_deposit: 391");
    }
    float coupledMass = shellMass + ejectaMass;
    eKinetic = coupledMomenta * coupledMomenta / (2.0 *(coupledMass)) * SolarMass * 1e10;
    if (eKinetic > (nSNII+nSNIA) * 1e51 && !winds){
        fprintf(stdout, "STARSS_FB: Rescaling high kinetic energy %e -> ", eKinetic);
        coupledMomenta = sqrt(2.0 * (coupledMass*SolarMass) * ejectaEnergy)/SolarMass/1e5;
        eKinetic = coupledMomenta * coupledMomenta / (2.0 *(coupledMass) * SolarMass) * SolarMass * 1e10;
        
        fprintf(stdout, "STARSS_FB:  %e; new p = %e\n", eKinetic, coupledMomenta);
    }
    
    // if (printout)
        if (printout) fprintf(stdout, "STARSS_FB: Ekinetic = %e Mass = %e\n",
                eKinetic, dmean * pow(LengthUnits * CellWidth[0][0], 3) / SolarMass);
    if (eKinetic > 1e60 && winds)
    {
        fprintf(stdout, "STARSS_FB: winds Ekinetic = %e Mass = %e\n",
                eKinetic, dmean * pow(LengthUnits * CellWidth[0][0], 3) / SolarMass);
        ENZO_FAIL("winds Ekinetic > reasonability!\n");
    }
    if (eKinetic > 1e60 && !winds)
    {
        fprintf(stdout, "STARSS_FB: Ekinetic = %e Mass = %e\n",
                eKinetic, dmean * pow(LengthUnits * CellWidth[0][0], 3) / SolarMass);
        ENZO_FAIL("SNE Ekinetic > reasonability!\n");
    }

    float coupledGasEnergy = max(ejectaEnergy - eKinetic, 0);
    // if (printout)
        if (printout) fprintf(stdout, "STARSS_FB: Coupled Gas Energy = %e\n", coupledGasEnergy);
    if (dxRatio > 1.0 && !winds){ // if we apply this reduction to winds, then there is literally *no* effect, even at Renaissance resolution.
        coupledGasEnergy = (coupledGasEnergy * pow(dxRatio, -6.5));
        if (printout)
            fprintf(stdout, "STARSS_FB: Reducing gas energy... GE = %e\n", coupledGasEnergy);
    }
    float coupledMetals = 0.0, SNIAmetals = 0.0, SNIImetals = 0.0, P3metals = 0.0;
    if (winds)
        coupledMetals = ejectaMetal; //+ shellMetals; // winds only couple to metallicity
    if (AnalyticSNRShellMass && !winds) 
        coupledMetals += shellMetals;
    SNIAmetals = (StarMakerTypeIaSNe) ? nSNIA * 1.4 : 0.0;
    if (!StarMakerTypeIaSNe)
        ejectaMetal += nSNIA * 1.4;
    SNIImetals = (StarMakerTypeIISNeMetalField) ? nSNII * (1.91 + 0.0479 * max(starMetal, 1.65)) : 0.0;
    if (!StarMakerTypeIISNeMetalField)
        ejectaMetal += nSNII * (1.91 + 0.0479 * max(starMetal, 1.65));
    if (isP3 && MechStarsSeedField)
        P3metals = ejectaMetal;
    coupledMetals += ejectaMetal;
    if (printout)
        fprintf(stdout, "STARSS_FB: Coupled Metals: %e %e %e %e %e %e\n", ejectaMetal, SNIAmetals, SNIImetals, shellMetals, P3metals, coupledMetals);

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
        for (int k = max(0,kp-3); k <= min(kp+3, GridDimension[2]); ++k)
            for (int j = max(0,jp-3); j <= max(jp+3, GridDimension[1]); ++j)
                for (int i = max(0,ip-3); i <= max(ip+3, GridDimension[0]); ++i)

                {
                    int idx = i + j * GridDimension[0] + k * GridDimension[0] * GridDimension[1];
                    if (idx > 0 && idx < size){
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
                }}
    }
    /* Reduce coupled quantities to per-particle quantity and converMetalNumt to 
        code units.
        Hopkins has complicated weights due to complicated geometry. 
            This implementation is simple since our coupled particles are 
            spherically symmetric about the feedback particle*/

    if (!winds)
        coupledEnergy = min((nSNII + nSNIA) * 1e51, eKinetic);

    
    /* Subtract shell mass from the central cells */
    float minusRho = 0; 
    float minusZ = 0;
    bool massiveCell;
    float remainMass = shellMass/MassUnits;
    float msubtracted = 0;
    float remainZ = shellMetals/MassUnits;
    if (shellMass > 0 && AnalyticSNRShellMass)
        // do
        {
            float zsubtracted = 0;
            massiveCell=false;
            msubtracted = 0;
                for (int i = ip-1; i <= ip+1; i++)
                    for (int j = jp-1; j <= jp+1; j++)
                        for (int k = kp-1; k <= kp+1; k++)
                            {
                            int flat = i + j * GridDimension[0] + k * GridDimension[0] * GridDimension[1];
                            if (flat > 0 && flat < size){
                                float xcell, ycell, zcell;
                                /* position of cell to couple to, center value */
                                xcell = CellLeftEdge[0][0] + (0.5 + (float) i) * dx;
                                ycell = CellLeftEdge[1][0] + (0.5 + (float) j) * dx;
                                zcell = CellLeftEdge[2][0] + (0.5 + (float) k) * dx;
                                float window = Window(*xp - xcell, *yp - ycell, *zp - zcell, dx, false); // always use cic to remove the mass
                                if (window > 0)
                                {
                                    float dpre = BaryonField[DensNum][flat] ;
                                    float zpre = BaryonField[MetalNum][flat];
                                    float pre_z_frac = zpre / dpre;
                                    if (printout)
                                    fprintf(stdout, "STARSS: Baryon Prior: %e, window = %f; mc = %e, ms = %e; m_z = %e , z = %e\n", BaryonField[DensNum][flat] * MassUnits, window, centralMass, shellMass, shellMetals, pre_z_frac);
                                    BaryonField[DensNum][flat] = max(dpre - remainMass * window, (1.0-maxEvacFraction)* dpre);
                                    minusRho += dpre - BaryonField[DensNum][flat];
                                    msubtracted += dpre - BaryonField[DensNum][flat];
                                    // fprintf(stdout, "STARSS: Baryon Post: %e\n", BaryonField[DensNum][flat] * MassUnits);
                                    BaryonField[MetalNum][flat] =   max(tiny_number, zpre - remainZ*window);
                                    minusZ += zpre - BaryonField[MetalNum][flat];
                                    zsubtracted += zpre - BaryonField[MetalNum][flat];
                                }
                            }
                            // if (BaryonField[DensNum][flat] >= 1.25*remainMass)
                            //     massiveCell=true;
                        }
                    remainMass -= msubtracted;
                    remainZ -= zsubtracted;
        } 
        // while (remainMass > 0 && msubtracted > 0);
    // if there wasnt enough mass to evacuate in host cells, need 
    // to rescale the shell+ejecta that will be deposited below
    // it also means that we couldnt build a shell to host the 
    // real amount of momenta.  Not sure what to do on that one...
    minusRho *= MassUnits; // Msun
    minusZ *= MassUnits;
    if (minusRho != coupledMass - ejectaMass && shellMass > 0 && AnalyticSNRShellMass)
    {
        if (printout) fprintf(stdout, "STARSS_FB: Of %e, only subtracted %e; rescaling the coupling mass\n", shellMass, minusRho);
        float oldcouple = coupledMass;
        coupledMass = minusRho + ejectaMass;
        coupledMetals = minusZ + ejectaMetal;
        // if were in here, we found an inconsistent place where the mass within cannot support the expected shell.  
        // the only real choice, to keep things from exploding in a bad way, is to rescale the momentum accordingly. 
        // If this is triggered, dont expect the terminal momenta relations to hold up.
        if (coupledMomenta * coupledMomenta / (2.0 * coupledMass * SolarMass) * SolarMass * 1e10 > 1e51){
            coupledMomenta = min(coupledMomenta, sqrt(2.0 * (coupledMass*SolarMass) * ejectaEnergy)/SolarMass/1e5);
        if (printout) fprintf(stdout, "STARSS_FB: rescaled momentum to %e (est KE = %e)\n", coupledMomenta, coupledMomenta*coupledMomenta / (2*coupledMass) * SolarMass * 1e10);
        }
    }
    coupledEnergy = coupledEnergy / EnergyUnits;
    coupledGasEnergy = coupledGasEnergy / EnergyUnits;
    coupledMass /= MassUnits;
    coupledMetals /= MassUnits;

    // conversion includes km -> cm
    coupledMomenta = coupledMomenta / MomentaUnits / sqrt((float) nCouple);
    SNIAmetals /= MassUnits;
    SNIImetals /= MassUnits;
    P3metals /= MassUnits;

    /* deposit the particles with their respective quantities */

    FLOAT LeftEdge[3] = {CellLeftEdge[0][0], CellLeftEdge[1][0], CellLeftEdge[2][0]};
    float eCouple, geCouple, mCouple, zCouple, zIICouple, zIACouple, p3Couple;
    float sum_cell_mass = 0;
    float expect_momenta = 0;
    float expect_error []= {0,0,0};
    float expect_z = 0;
    int n_pos[] = {0,0,0};    float ge_deposit = 0;
    float mtouched = 0;
    float ke_deposit = 0;

    std::vector<double> te_buffer(125, 0.0);
    std::vector<double> rho_buffer(125,0.0);
    std::vector<double> ge_buffer(125, 0.0);
    /* get prior mass of cell */
    for (int k=-2; k<3; k++)
        for(int j=-2; j<3; j++)
            for (int i=-2; i<3; i++){
                int flat = ip+i + (jp+j) * GridDimension[0] + (kp+k) * GridDimension[0] * GridDimension[1];
                int buffind = i+2 + (j+2) * 5 + (k+2) * 25;
                rho_buffer[buffind] = BaryonField[DensNum][flat];
            }
    /* Couple the main qtys from cloud particles to grid */
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

                    // printf("STARSS_FB: Coupling particle near %e %e %e => %d %d %d \n\t\t(%e %e %e), dx = %e\n", xc, yc, zc, 
                    //                                                                 ic, jc, kc,
                    //                                                                 (*xp - CellLeftEdge[0][0] - 0.5 * dx) / dx,
                    //                                                                 (*yp - CellLeftEdge[1][0] - 0.5 * dx) / dx,
                    //                                                                 (*zp - CellLeftEdge[2][0] - 0.5 * dx) / dx,
                    //                                                                 dx);
                    if  ( i == 0 && j == 0 && k == 0) continue;
                            // compute coupled momenta here to get the direction right
                            float modi, modj, modk;
                            modi = (float) i;
                            modj = (float) j;
                            modk = (float) k;
                            float uniti, unitj, unitk, normfact;
                            float pX=0, pY=0, pZ = 0;
                            // if (!(modi==0 && modj == 0 && modk==0)){
                                normfact = 1; //sqrt(modi*modi + modj*modj + modk*modk);
                                uniti = modi / normfact;
                                unitj = modj / normfact;
                                unitk = modk / normfact;
                            
                                pX = coupledMomenta * (float) uniti; //* weightsVector[n];
                                pY = coupledMomenta * (float) unitj; //* weightsVector[n];
                                pZ = coupledMomenta * (float) unitk; //* weightsVector[n];
                            // }



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
                    // int np = 1;




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
                            if (flat > 0 && flat < size){
                                int buffind = (ii+i +2) + (jj+j+2) * 5 + (kk+k+2)*25;
                                                            float modi, modj, modk;
                                modi = (float) i;
                                modj = (float) j;
                                modk = (float) k;
                                float uniti, unitj, unitk, normfact;
                                float pX=0, pY=0, pZ = 0;
                                // if (!(modi==0 && modj == 0 && modk==0)){
                                    normfact = sqrt(modi*modi + modj*modj + modk*modk);
                                    uniti = modi / normfact;
                                    unitj = modj / normfact;
                                    unitk = modk / normfact;
                                
                                    pX = coupledMomenta * (float) uniti; //* weightsVector[n];
                                    pY = coupledMomenta * (float) unitj; //* weightsVector[n];
                                    pZ = coupledMomenta * (float) unitk; //* weightsVector[n];
                                // }
                                /* fraction of quantities that go into this cloud cell... */
                                float window = Window(xc - xcell, yc - ycell, zc - zcell, dx, NGP);
                                dep_vol_frac += window;
                                
                                // printf("STARSS_FB: \t\t\tCloud %e %e %e Coupling at %d %d %d with window %e.  \n\t\t\t\twfactors: %e %e %e (%f %f %f)\n", 
                                //                                             xcell, ycell, zcell, 
                                //                                             icell, jcell, kcell, window,
                                //                                             xc-xcell, yc-ycell, zc-zcell,
                                //                                             (xc-xcell)/dx, (yc-ycell)/dx, (zc-zcell)/dx);
                                /* add quantities to this cell accordingly */

                                BaryonField[DensNum][flat] = BaryonField[DensNum][flat] + mCouple * window;
                                sum_cell_mass += window * BaryonField[DensNum][flat];
                                BaryonField[MetalNum][flat] = BaryonField[MetalNum][flat] + zCouple * window;
                                if (DualEnergyFormalism)
                                    {    
                                        // add kinetic energy on a second pass to ensure that all masses have been previously adjusted.
                                        // printf("STARSS_FB: ADDING %e to TE buffer %d (%e %e %e %e)\n",(pX*pX + pY*pY + pZ*pZ) * window * window * MomentaUnits, buffind, pX*MomentaUnits, pY*MomentaUnits, pZ*MomentaUnits, window );
                                        te_buffer[buffind] += (pX*pX + pY*pY + pZ*pZ) * window * window ;
                                        // printf("STARSS_FB: added %e ke\n",(pX*pX + pY*pY + pZ*pZ)/ (2.*coupledMass) * window  * EnergyUnits );
                                        // ke_deposit += (pX*pX + pY*pY + pZ*pZ)/ (2.*coupledMass) * window ;
                                        // BaryonField[TENum][flat] = BaryonField[TENum][flat] + geCouple * window;
                                    }
                                ge_buffer[buffind] += geCouple * window;
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
                            }
                            // TODO: add other metal fields for consistency.

                        }
                    // printf("STARSS_FB: CIC deposited in sum window %e\n", dep_vol_frac);
                } // end coupling main qtys
                // subtract out shell mass before adding in KE


    if (printout) fprintf(stdout, "STARSS_FB: Coupled %e gas energy\n", ge_deposit * EnergyUnits);       
    if (printout){
        fprintf(stdout, "STARSS_FB: After deposition, counted %e Msun km/s momenta deposited.  Error = %e..\n", sqrt(expect_momenta) * MomentaUnits, 
                                                                                        (expect_error[0]+expect_error[1]+expect_error[2]) * MomentaUnits);
        fprintf(stdout, "STARSS_FB: \tTabulated %e Msun deposited metals\n", expect_z * MassUnits);
    }


    if (criticalDebug && !winds)
    {
        for (int k = max(0,kp-3); k <= min(kp+3, GridDimension[2]); ++k)
            for (int j = max(0,jp-3); j <= max(jp+3, GridDimension[1]); ++j)
                for (int i = max(0,ip-3); i <= max(ip+3, GridDimension[0]); ++i)
                {
                    int idx = i + j * GridDimension[0] + k * GridDimension[0] * GridDimension[1];  
                    if (idx > 0 && idx < size){          
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
                }

        if (printout)
            fprintf(stdout, "STARSS_FB: Difference quantities: dxRatio = %f dMass = %e dZ = %e dzII = %e dxIa = %e  P = %e |P| = %e (%e) TE = %e GE = %e coupledGE = %e Ej = %e Mej = %e Zej = %e\n",
                    dxeff, (postMass - preMass) * MassUnits, (postZ - preZ) * MassUnits,
                    (postZII - preZII) * MassUnits, (postZIa - preZIa) * MassUnits,
                    (postP - preP) * MomentaUnits,
                    (sqrt(postPmag) - sqrt(prePmag))*MomentaUnits, sqrt(expect_momenta) * MomentaUnits,
                    (postTE - preTE) * EnergyUnits, (postGE - preGE) * EnergyUnits,
                    coupledGasEnergy * EnergyUnits, ejectaEnergy,
                    ejectaMass, ejectaMetal);
    }

    transformComovingWithStar(BaryonField[DensNum], BaryonField[MetalNum],
                            BaryonField[MetalIINum], BaryonField[MetalIaNum],
                            BaryonField[Vel1Num], BaryonField[Vel2Num],
                            BaryonField[Vel3Num],
                            BaryonField[TENum], BaryonField[GENum],
                            *up, *vp, *wp,
                            GridDimension[0], GridDimension[1],
                            GridDimension[2], -1);
        float affected_mass = 0.0;
        float sum_te = 0.0;
        for (float e: te_buffer)
            sum_te += e;
        if (sum_te > 0){
            for (int i = -2; i <= 2; i++)
                for (int j = -2; j <= 2; j++)
                    for (int k = -2; k <= 2; k++)
                        {
                            int flat = ip+i + (jp+j) * GridDimension[0] + (kp+k) * GridDimension[0] * GridDimension[1];
                            int buffind = i+2 + (j+2) * 5 + (k+2) * 25;
                            if (te_buffer[buffind] > 0)
                                affected_mass += BaryonField[DensNum][index];
                        }
            for (int i = -2; i <= 2; i++)
                for (int j = -2; j <= 2; j++)
                    for (int k = -2; k <= 2; k++)
                        {
                            int flat = ip+i + (jp+j) * GridDimension[0] + (kp+k) * GridDimension[0] * GridDimension[1];
                            int buffind = i+2 + (j+2) * 5 + (k+2) * 25;

                            // printf("STARSS_FB: adding %e from buffer %d to TE\n", te_buffer[buffind] / (2.0 * BaryonField[DensNum][flat]/EnergyUnits), buffind);
                            if (DualEnergyFormalism)
                                {    
                                    // add kinetic energy on a second pass to ensure that all masses have been previously adjusted. TE / old + TEN / (old+new) = TE * (old+new)/old + TEN/ (old+new)
                                    BaryonField[TENum][flat] = BaryonField[TENum][flat] *BaryonField[DensNum][flat] /rho_buffer[buffind]
                                                            +  (te_buffer[buffind] 
                                                            / (2 * coupledMass)) / BaryonField[DensNum][flat];
                                    ke_deposit += (te_buffer[buffind] / (2 * coupledMass));
                                    // BaryonField[TENum][flat] += ge_buffer[buffind];
                                }
                            // BaryonField[GENum][flat] += ge_buffer[buffind];
                                

                                } // end cloud-particle cic
        }
    if (printout) fprintf(stdout, "STARSS_FB: Counted %e Kenergy deposited\n", ke_deposit*EnergyUnits);
    // fill in the remaining energy with thermal deposition
    float remain = 1e51/EnergyUnits - ke_deposit;

    if (dxeff > 1)
        remain = remain * pow(dxeff, -6.5);

    if (printout) fprintf(stdout, "STARSS_FB: Coupling missing energy... remaining = %e\n", remain *EnergyUnits);
    if (remain > 0)
        for (int i = -1; i <= 1; i++)
            for (int j = -1; j <= 1; j++)
                for (int k = -1; k <= 1; k++)
                {
                    float xc, yc, zc;
                    xc = *xp + i * dx;
                    yc = *yp + j * dx;
                    zc = *zp + k * dx;
                    int ic, jc, kc;
                    ic = (xc - CellLeftEdge[0][0] - 0.5 * dx) / dx;
                    jc = (yc - CellLeftEdge[1][0] - 0.5 * dx) / dx;
                    kc = (zc - CellLeftEdge[2][0] - 0.5 * dx) / dx;
                    if  ( i == 0 && j == 0 && k == 0) continue;

                    float cGE = remain / (float) nCouple; //* weightsVector[n];
                    for (int kk = -1; kk <= 1; kk ++)
                        for (int jj = -1; jj <= 1; jj++)
                            for (int ii = -1; ii <= 1; ii++)
                            {
                                int buffind = (i+ii)+2 + ((j+jj)+2) * 5 + ((k+kk)+2) * 25;

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
                                if (flat > 0 && flat < size){
                                    
                                    /* fraction of quantities that go into this cloud cell... */
                                    float window = Window(xc - xcell, yc - ycell, zc - zcell, dx, NGP);
                                    if (DualEnergyFormalism && remain > 0)
                                        {    
                                            BaryonField[TENum][flat] = (BaryonField[TENum][flat]*rho_buffer[buffind] + cGE * window) /BaryonField[DensNum][flat];
                                        }
                                    if (remain > 0)
                                        BaryonField[GENum][flat] = (BaryonField[GENum][flat]*rho_buffer[buffind] + cGE * window) / BaryonField[DensNum][flat];
                                    ge_deposit += cGE * window;
                                }
                            }
                        // printf("STARSS_FB: CIC deposited in sum window %e\n", dep_vol_frac);
                    } // end coupling main qtys   
    if (printout) fprintf(stdout, "STARSS_FB: Counted %e GE deposited\n", ge_deposit*EnergyUnits);

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


    if (printout) printf("STARSS_FB\nSTARSS_FB\nSTARSS_FB\n");
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
