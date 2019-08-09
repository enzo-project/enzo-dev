#include <stdio.h>
#include <math.h>
#include <mpi.h>
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

    extern "C"  void FORTRAN_NAME(cic_deposit)(float* xPosition, float* yPosition,
        float* zPosition, int* gridRank, int* nParticles, 
        float* DepositQuantity, float* FieldToDepositTo,
        float* leftEdge, int* xGridDim, int* yGridDim, 
        int* zGridDim, float* gridDx, float* cloudsize);
    int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, float *MassUnits, float Time);
    int transformComovingWithStar(float* Density, float* Metals, 
        float* Vel1, float* Vel2, float* Vel3, 
        float up, float vp, float wp,
        int sizeX, int sizeY, int sizeZ, int direction);
    int FindField(int field, int farray[], int numfields);


int grid::MechStars_DepositFeedback(float ejectaEnergy, 
                        float ejectaMass, float ejectaMetal,
                        float* up, float* vp, float* wp,
                        float* xp, float* yp, float* zp,
                        int ip, int jp, int kp,
                        int size, float* muField, int winds){
    
    /*
     This routine will create an isocahedron of coupling particles, where we determine
        the feedback quantities.  The vertices of the isocahedron are ~coupled particles
        and all have radius dx from the source particle. 
        Each vertex particle will then be CIC deposited to the grid!
    */
    bool debug = true;
    bool criticalDebug = true;
    int index = ip+jp*GridDimension[0]+kp*GridDimension[0]*GridDimension[1];
    int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
    
    /* Compute size (in floats) of the current grid. */
    float stretchFactor =1.;//1.5/sin(M_PI/10.0);  // How far should cloud particles be from their host
                                // in units of dx
    int usePt = 0;
    size = 1;
    for (int dim = 0; dim < GridRank; dim++)
        size *= GridDimension[dim];
    float DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
                        HMNum, H2INum, H2IINum, DINum, DIINum, HDINum;
    /* Find fields: density, total energy, velocity1-3. */

    if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                        Vel3Num, TENum) == FAIL) {
        fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
        return FAIL;
    }
    /* Set the units */
    float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1,
            TimeUnits = 1, VelocityUnits = 1, MassUnits = 1;
    if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
            &TimeUnits, &VelocityUnits, 
            &MassUnits, this->Time) == FAIL) {
        fprintf(stderr, "Error in GetUnits.\n");
        return FAIL;    
    } 
    FLOAT dx = CellWidth[0][0];        

    if (debug)
        fprintf(stdout, "depositing quantities: Energy %e, Mass %e, Metals %e\n",
            ejectaEnergy, ejectaMass, ejectaMetal);

    /* 
        get metallicity field and set flag; assumed true thoughout feedback
        since so many quantities are metallicity dependent
     */
    int MetallicityField = FALSE, MetalNum;
    if ((MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields))
        != -1)
        MetallicityField = TRUE;
    else
        MetalNum = 0;

    /* set other units that we need */
    MassUnits = DensityUnits*pow(LengthUnits*dx, 3)/SolarMass;
    float EnergyUnits = DensityUnits*pow(LengthUnits*dx, 3) 
                    * VelocityUnits*VelocityUnits;//[g cm^2/s^2] -> code_energy
    float MomentaUnits = VelocityUnits;  

    /* Make copys of fields to work with. These will conatin the added deposition
        of all quantities and be coupled to the grid after the cic deposition. */
    float density [size];
    float metals [size];
    float *u = new float [size];
    float *v = new float[size];
    float *w = new float [size];
    float *totalEnergy =new float [size];
    float *gasEnergy =new float[size];
    for (int i=0; i<size; i++){
        density[i] = 0.0;
        metals[i] = 0.0;
        u[i] = 0.0;
        v[i] =0.0;
        w[i] = 0.0;
        totalEnergy[i] = 0.0;
        gasEnergy[i] = 0.0;
    }
    /* Transform coordinates so that metals is fraction (rho metal/rho baryon)
        u, v, w -> respective momenta.  Use -1 to reverse transform after.*/
    /* these temp arrays are implicitly comoving with the star! */
    /* Array of coordinates of the isocahedron vertices scaled by r=dx */
    
    float phi = (1.0+sqrt(5))/2.0; //Golden Ratio
    float iphi = 1.0/phi; // inverse GR
    /* Particle Vectors */

    
    /* A DODECAHEDRON+ISOCAHEDRON */

        int nCouple = 32;
        float A = stretchFactor*dx;

        /* Dodec points followed by isoca points */
        FLOAT CloudParticleVectorX [] = {1, 1, 1, 1, -1, -1, -1, -1, 0, 0, 0, 0,
                                            iphi, iphi, -iphi, -iphi, phi, phi, -phi, -phi, 
                                            0, 0, 0, 0, 1, 1, -1, -1, phi,-phi, phi,-phi};
        FLOAT CloudParticleVectorY [] = {1,1,-1,-1, 1, 1, -1, -1, iphi, iphi, -iphi, -iphi,
                                            phi, -phi, phi,-phi, 0, 0, 0, 0,1, 1, -1, -1, 
                                            phi, -phi, -phi, phi, 0, 0, 0, 0};
        FLOAT CloudParticleVectorZ []  = {1,-1, 1,-1, 1,-1, 1,-1, phi,-phi, phi,-phi, 
                                            0, 0, 0, 0, iphi, -iphi, iphi, -iphi, 
                                            phi, -phi, -phi, phi, 0, 0, 0, 0, 1, 1, -1, -1};
                            
            /* Set position of feedback cloud particles */

            FLOAT CloudParticlePositionX [nCouple];
            FLOAT CloudParticlePositionY [nCouple];
            FLOAT CloudParticlePositionZ [nCouple];

                /*all possible values of x,y,z with origin particle at x=y=z=0.0 */
            for (int cpInd = 0; cpInd < nCouple; cpInd++){
                FLOAT norm = sqrt(CloudParticleVectorX[cpInd]*CloudParticleVectorX[cpInd]
                    + CloudParticleVectorY[cpInd]*CloudParticleVectorY[cpInd]
                    + CloudParticleVectorZ[cpInd]*CloudParticleVectorZ[cpInd]);
                /* in this cloud, take the coupling particle position as 0.5, 0.5, 0.5 */
                CloudParticlePositionX[cpInd] = *xp-CloudParticleVectorX[cpInd]/norm*
                            A;
                CloudParticlePositionY[cpInd] = *yp-CloudParticleVectorY[cpInd]/norm*
                            A;
                CloudParticlePositionZ[cpInd] = *zp-CloudParticleVectorZ[cpInd]/norm*
                            A;
                /* turn the vectors back into unit-vectors */
                CloudParticleVectorZ[cpInd] /= norm;
                CloudParticleVectorY[cpInd] /= norm;
                CloudParticleVectorX[cpInd] /= norm;  
            }
    /* Each particle gets 1/12 of energy, momenta, mass, and metal.  There
    are no varying vector / scalar weights to worry about.  The momenta coupled
    is simply \hat(r_ba) p/12 for r_ba the vector from source to coupled
    particle.  */

    /* transform to comoving with the star and take velocities to momenta */

    transformComovingWithStar(BaryonField[DensNum], BaryonField[MetalNum], 
                        BaryonField[Vel1Num],BaryonField[Vel2Num],BaryonField[Vel3Num],
                        *up, *vp, *wp, GridDimension[0], GridDimension[1],
                        GridDimension[2], 1);
    /* Use averaged quantities across multiple cells so that deposition is stable.
        vmean is used to determine whether the supernova shell calculation should proceed:
            M_shell > 0 iff v_shell > v_gas */
    float zmean=0, dmean=0, nmean=0, vmean;
    for (int ind = -1; ind <= 1; ind++){
        zmean += BaryonField[MetalNum][index+ind]*BaryonField[DensNum][index+ind];
        zmean += BaryonField[MetalNum][index+GridDimension[0]*ind]*BaryonField[DensNum][index+GridDimension[0]*ind];
        zmean += BaryonField[MetalNum][index+GridDimension[0]*GridDimension[1]*ind]*BaryonField[DensNum][index+GridDimension[0]*GridDimension[1]*ind];

        if (debug) fprintf(stdout, "MuField = %f %f %f\nmax = %d ; %d %d %d\n",
            muField[index+ind], muField[index+GridDimension[0]*ind], muField[index+ind*GridDimension[0]*GridDimension[1]], GridDimension[0]*GridDimension[1]*GridDimension[2],
            index+ind, index+GridDimension[0]*ind, index+ind*GridDimension[0]*GridDimension[1]);
        nmean += BaryonField[DensNum][index+ind]*BaryonField[DensNum][index+ind]*DensityUnits/mh/max(0.6,muField[index+ind]);
        nmean += BaryonField[DensNum][index+GridDimension[0]*ind]*
            BaryonField[DensNum][index+GridDimension[0]*ind]*DensityUnits/mh/max(muField[index+GridDimension[0]*ind],0.6);
        nmean += BaryonField[DensNum][index+GridDimension[0]*GridDimension[1]*ind]
            *BaryonField[DensNum][index+GridDimension[0]*GridDimension[1]*ind]*DensityUnits/mh/max(0.6,muField[index+ind*GridDimension[0]*GridDimension[1]]);

        dmean += BaryonField[DensNum][index+ind];
        dmean += BaryonField[DensNum][index+GridDimension[0]*ind];
        dmean += BaryonField[DensNum][index+GridDimension[0]*GridDimension[1]*ind];
        
        vmean += BaryonField[Vel1Num][index+ind]*BaryonField[Vel1Num][index+ind]*VelocityUnits*VelocityUnits;
        vmean += BaryonField[Vel1Num][index+GridDimension[0]*ind]*BaryonField[Vel1Num][index+GridDimension[0]*ind]*VelocityUnits*VelocityUnits;
        vmean += BaryonField[Vel1Num][index+GridDimension[0]*GridDimension[1]*ind]*BaryonField[Vel1Num][index+GridDimension[0]*GridDimension[1]*ind]*VelocityUnits*VelocityUnits;

        vmean += BaryonField[Vel2Num][index+ind]*BaryonField[Vel2Num][index+ind]*VelocityUnits*VelocityUnits;
        vmean += BaryonField[Vel2Num][index+GridDimension[0]*ind]*BaryonField[Vel2Num][index+GridDimension[0]*ind]*VelocityUnits*VelocityUnits;
        vmean += BaryonField[Vel2Num][index+GridDimension[0]*GridDimension[1]*ind]*BaryonField[Vel2Num][index+GridDimension[0]*GridDimension[1]*ind]*VelocityUnits*VelocityUnits;

        vmean += BaryonField[Vel3Num][index+ind]*BaryonField[Vel3Num][index+ind]*VelocityUnits*VelocityUnits;
        vmean += BaryonField[Vel3Num][index+GridDimension[0]*ind]*BaryonField[Vel3Num][index+GridDimension[0]*ind]*VelocityUnits*VelocityUnits;
        vmean += BaryonField[Vel3Num][index+GridDimension[0]*GridDimension[1]*ind]*BaryonField[Vel3Num][index+GridDimension[0]*GridDimension[1]*ind]*VelocityUnits*VelocityUnits;

    }
    vmean = sqrt(vmean/dmean/dmean);
    zmean /= (dmean*0.02);
    nmean /= (dmean);
    
    dmean /= 9.0;
    nmean = max(nmean, 1e-3);
    float zZsun = max(zmean, 1e-8);
    float fz = (zZsun < 0.01)? (2.0): (pow(zZsun, -0.14));

    /* conversions */
    float CoolingRadius = 28.4 *
        pow(max(0.001,nmean), -3.0/7.0)
        *pow(ejectaEnergy/1.0e51, 2.0/7.0)* fz;
    if (debug)fprintf(stdout, "cooling radius [pc] = %f\n %f %f %f %e %f \n", 
            CoolingRadius, nmean, ejectaEnergy/1e51, fz, zmean, dmean);
    /* Calculate coupled energy scaled by reduction to account for unresolved
    cooling, then use that energy to calculate momenta*/
    float coupledEnergy = ejectaEnergy;

    if (debug)fprintf(stdout, "Dx [pc] = %f\n", dx*LengthUnits/pc_cm);
    float dxRatio = stretchFactor*dx*LengthUnits/pc_cm/CoolingRadius;

    float coupledMomenta = 0.0;
    float eKinetic = 0.0;
    /* Hopkins uses ratio of masses to determine how to couple.
        Radius here is well-known and fixed, so we use that instead */
    if (dxRatio > 1.0){ 
        if (ejectaEnergy < 1e5 || dxRatio > 100){
            coupledEnergy = 0.0;
            coupledMomenta = 0.0;
        }else{
            coupledEnergy = ejectaEnergy*pow(dxRatio, -6.5);
            usePt = 1;
            
        /* Determine coupled momenta if rc < dx 
        else couple p_ej*(1+dx/r_cool)**4 */
            if(debug)fprintf(stdout, "Using P_t with Nb = %f, E= %e",nmean, coupledEnergy/1e51);
            float Efactor = 1.0;
            if (dxRatio > 4) Efactor = coupledEnergy/ejectaEnergy*pow(dxRatio,3);
            coupledMomenta = 4.8e5*pow(nmean, -1.0/7.0)
                * pow(ejectaEnergy/1e51, 13.0/14.0) * fz; //Msun*km/s
        }
    } else {
        if (debug)fprintf(stdout, "Directly calculating momenta using energy = %e and mass = %e ", 
                    ejectaEnergy, ejectaMass);
        coupledMomenta = pow(2.0*ejectaEnergy*(ejectaMass*SolarMass), 0.5) 
                                * pow(1.0+dxRatio, 3.75*pow(nmean, -1./14.))/SolarMass/1e5; //Msun*km/s
        if (debug)fprintf(stdout, "Calculated p = %e ", coupledMomenta);
        if (debug)fprintf(stdout, "Ekinetic = %e\n", coupledMomenta*coupledMomenta
                    /(2.0*ejectaMass)*SolarMass*1e10);

    }
    float shellMass = 0.0, shellVelocity = 0.0;
    /* If resolution is in a range compared to Rcool and
        Analytic SNR shell mass is on, adjust the shell mass 
        Shell is limited on upper end by 1/1000 mass of 
            cell with mean density*/
    float maxShellMass = 1.0*DensityUnits*MassUnits/1000; // large shell mass evacuates too efficently... muField-> 0.0 and NaN ensues!
    if (dxRatio <= 50 && dxRatio >= 0.1 && coupledEnergy > 0
        && AnalyticSNRShellMass){
            shellVelocity = 413.0 *pow(nmean, 1.0/7.0)
                *pow(zZsun, 3.0/14.0)*pow(coupledEnergy/EnergyUnits/1e51, 1.0/14.0)
                *pow(dxRatio, -7.0/3.0);//km/s
            /* Underdense regions can have large momenta with 
                low velocity, leading to shell mass instability.  
                The shell velocity is compared to gas velocity, and 
                can only contribute to the mass if the shell velocity is 
                higher than the gas velocity.*/
            if (shellVelocity > vmean){ 
                shellMass = max(1e3, coupledMomenta/shellVelocity); //Msun
                
                /* cant let host cells evacuate completely!  
                    Shell mass will be evacuated from central cells by CIC a negative mass,
                    so have to check that the neighbors can handle it too*/
                for (int ind = -1; ind<=1; ind++){
                        float minD = min(BaryonField[DensNum][index+ind],BaryonField[DensNum][index+GridDimension[0]*ind]);
                        minD = min(minD,BaryonField[DensNum][index+GridDimension[0]*GridDimension[1]*ind]);
                        if (shellMass >= 0.05*minD*MassUnits)
                        shellMass = 0.05*minD*MassUnits;
                }
                    
                         
            }

    }

    if (shellMass < 0.0){
       fprintf(stdout, "Shell mass = %e Velocity= %e P = %e",
            shellMass, shellVelocity, coupledMomenta);
         ENZO_FAIL("SM_deposit: 252");}
    float coupledMass = shellMass+ejectaMass;
    eKinetic = coupledMomenta*coupledMomenta
                    /(2.0*coupledMass)*SolarMass*1e10;

    /*  rescale momenta if it results in too much energy */
    
    // if (eKinetic > coupledEnergy && !usePt){
    //     float fact = coupledEnergy/eKinetic;
    //     if (debug)fprintf(stdout, "recalculating momenta: e_k > e_cpl: e_k = %e e_cpl = %e factor = %e ",
    //         eKinetic, coupledEnergy, fact);
    //     coupledMomenta = pow(fact*2.0*ejectaEnergy*(coupledMass*SolarMass), 0.5)
    //                             * pow(1.0+dxRatio, 3.75*pow(nmean, -1./14.))/SolarMass/1e5;
    //     eKinetic = coupledMomenta*coupledMomenta
    //                 /(2.0*coupledMass)*SolarMass*1e10; 
    // if (debug)fprintf(stdout, "new e_k = %e p = %e\n",eKinetic, coupledMomenta);
    // }

    // // // /* If p_t gives too much kinetic energy, reduce it
    // // //     to preserve energy conservation */

    // if (eKinetic > coupledEnergy && usePt){
    //     float fact = pow(coupledEnergy/eKinetic,14.0/13.0);
    //     if (debug)fprintf(stdout, "recalculating momenta: e_k > e_cpl e_k = %e e_cpl = %e ",
    //         eKinetic, coupledEnergy);
    //     coupledMomenta = pow(dxRatio, -3)*4.8e5*pow(nmean, -1.0/7.0)
    //             * pow(ejectaEnergy/1e51, 13.0/14.0) * fz;
    //     eKinetic = coupledMomenta*coupledMomenta
    //                 /(2.0*coupledMass)*SolarMass*1e10;
    // if (debug)fprintf(stdout, "new e_k = %e p = %e\n",eKinetic, coupledMomenta);
    // }

    float coupledGasEnergy = max(ejectaEnergy-eKinetic, 0);
    if (debug)fprintf(stdout, "Coupled Gas Energy = %e\n",coupledGasEnergy);
    if (dxRatio > 1.0)
        coupledGasEnergy *= pow(dxRatio, -6.5);
    /* rescale momentum for new shell */
    float shellMetals = zZsun*0.02 * shellMass;
    float coupledMetals = ejectaMetal + shellMetals;



    if (debug) fprintf(stdout, "Coupled Momentum: %e\n", coupledMomenta/float(nCouple));
    /* Reduce coupled quantities to per-particle quantity and convert to 
        code quantities.
        Hopkins has complicated weights due to complicated geometry. 
            This is more simple since our coupled particles are 
            spherically symmetric about the feedback particle*/
    
    coupledEnergy /= float(nCouple);
    coupledGasEnergy /= float(nCouple);
    coupledMass /= float(nCouple);
    coupledMetals /= float(nCouple);
    coupledMomenta /= float(nCouple);
    /* Transform coupled quantities to code units */
    coupledEnergy /= EnergyUnits;
    coupledGasEnergy /= EnergyUnits;
    coupledMass /= MassUnits;
    coupledMetals /= MassUnits;
    coupledMomenta /= MomentaUnits;
    /* CIC deposit the particles with their respective quantities */
    float LeftEdge[3] = {CellLeftEdge[0][0], CellLeftEdge[1][0], CellLeftEdge[2][0]};
    for (int n = 0; n < nCouple; n++){

        FLOAT pX = coupledMomenta*CloudParticleVectorX[n];
        FLOAT pY = coupledMomenta*CloudParticleVectorY[n];
        FLOAT pZ = coupledMomenta*CloudParticleVectorZ[n];
        int np = 1;

        FORTRAN_NAME(cic_deposit)(&CloudParticlePositionX[n], &CloudParticlePositionY[n],
            &CloudParticlePositionZ[n], &GridRank, &np,&coupledMass, &density[0], LeftEdge, 
            &GridDimension[0], &GridDimension[1], &GridDimension[2], &dx, &dx);
        FORTRAN_NAME(cic_deposit)(&CloudParticlePositionX[n], &CloudParticlePositionY[n],
            &CloudParticlePositionZ[n], &GridRank, &np,&coupledMetals, &metals[0], LeftEdge, 
            &GridDimension[0], &GridDimension[1], &GridDimension[2], &dx, &dx);
        if (pX != 0.0){
            FORTRAN_NAME(cic_deposit)(&CloudParticlePositionX[n], &CloudParticlePositionY[n],
                &CloudParticlePositionZ[n], &GridRank, &np,&pX, u, LeftEdge, 
                &GridDimension[0], &GridDimension[1], &GridDimension[2], &dx, &dx);
        }
        if (pY != 0.0)
            FORTRAN_NAME(cic_deposit)(&CloudParticlePositionX[n], &CloudParticlePositionY[n],
                &CloudParticlePositionZ[n], &GridRank,&np,&pY, v, LeftEdge, 
                &GridDimension[0], &GridDimension[1], &GridDimension[2], &dx, &dx);
        if (pZ != 0.0)
            FORTRAN_NAME(cic_deposit)(&CloudParticlePositionX[n], &CloudParticlePositionY[n],
                &CloudParticlePositionZ[n], &GridRank,&np,&pZ, w, LeftEdge, 
                &GridDimension[0], &GridDimension[1], &GridDimension[2], &dx, &dx);
        if (coupledEnergy > 0 && DualEnergyFormalism)
            FORTRAN_NAME(cic_deposit)(&CloudParticlePositionX[n], &CloudParticlePositionY[n],
                &CloudParticlePositionZ[n], &GridRank,&np,&coupledEnergy, totalEnergy, LeftEdge, 
                &GridDimension[0], &GridDimension[1], &GridDimension[2], &dx, &dx);
        if (coupledGasEnergy > 0)
            FORTRAN_NAME(cic_deposit)(&CloudParticlePositionX[n], &CloudParticlePositionY[n],
                &CloudParticlePositionZ[n], &GridRank,&np,&coupledGasEnergy, gasEnergy, LeftEdge, 
                &GridDimension[0], &GridDimension[1], &GridDimension[2], &dx, &dx); 
    }
    /* Deposit one negative mass particle centered on star to account for 
        shell mass leaving host cells */
    int np = 1;
    shellMass *= -1/MassUnits;
    FORTRAN_NAME(cic_deposit)(xp, yp, zp, &GridRank,&np,&shellMass, &density[0], LeftEdge, 
        &GridDimension[0], &GridDimension[1], &GridDimension[2], &dx, &dx);
    shellMetals *= -1/MassUnits;   
    FORTRAN_NAME(cic_deposit)(xp, yp, zp, &GridRank,&np,&shellMetals, &metals[0], LeftEdge, 
        &GridDimension[0], &GridDimension[1], &GridDimension[2], &dx, &dx);
    
    
    /* transform the grid to comoving with star ; wouldnt recommend this on root grid if its too big...*/

    float preMass = 0, preZ = 0, preP = 0, prePmag=0, preTE = 0, preGE = 0;
    float dsum = 0.0, zsum=0.0, psum=0.0, psqsum =0.0, tesum=0.0, gesum=0.0, kesum=0.0;
    float postMass = 0, postZ = 0, postP = 0, postPmag = 0, postTE = 0, postGE = 0;
    if (criticalDebug){
        for (int i=0; i<size; i++){
            preMass += BaryonField[DensNum][i];
            preZ += BaryonField[MetalNum][i];
            preP += BaryonField[Vel1Num][i]+BaryonField[Vel2Num][i]+BaryonField[Vel3Num][i];
            prePmag += pow(BaryonField[Vel1Num][i]*MomentaUnits,2)+
                pow(BaryonField[Vel2Num][i]*MomentaUnits,2)
                +pow(BaryonField[Vel3Num][i]*MomentaUnits,2);
            preTE += BaryonField[TENum][i];
            preGE += BaryonField[GENum][i];
        }
    }
        /* Since wind energy is so low, if we want to couple something
        it will have to thermal at host cell.   */
    if (winds && DepositUnresolvedEnergyAsThermal && coupledEnergy == 0){
        totalEnergy[index] = (double(ejectaEnergy)*1e30)/double(EnergyUnits)/BaryonField[DensNum][index]/1e30;
    }
    for (int i = 0; i < size; i++){ 
                
        float delta = (density[i])
                            /(density[i]+BaryonField[DensNum][i]);
        /* Couple placeholder fields to the grid, account 
            for grids that got initialized to -0.0*/
        BaryonField[DensNum][i] += density[i];

        //Metals transformed back to density in transform routine
        
        BaryonField[MetalNum][i] += metals[i]; 
        BaryonField[TENum][i] += 
                    totalEnergy[i]/BaryonField[DensNum][i];
        
        BaryonField[GENum][i] += 
                    gasEnergy[i]/BaryonField[DensNum][i];
        BaryonField[Vel1Num][i] += u[i];
        BaryonField[Vel2Num][i] += v[i];
        BaryonField[Vel3Num][i] += w[i];
    }


    /* Sum of feedback quantities: */
    if (debug){
        for (int i = 0; i<size; i++){
            dsum += density[i];
            zsum += metals[i];
            psum += u[i]+v[i]+w[i];
            tesum += totalEnergy[i];
            gesum += gasEnergy[i];
            psqsum += (u[i]*u[i]+v[i]*v[i]+w[i]*w[i])*MomentaUnits*MomentaUnits;
            kesum += (density[i] > 0)? 
                        (u[i]*u[i]+v[i]*v[i]+w[i]*w[i])*MomentaUnits*MomentaUnits
                        /(2*density[i]*MassUnits)*SolarMass*1e10
                        : 0;

        }

       fprintf(stdout, "Sum Mass  = %e  ", dsum*MassUnits);
       fprintf(stdout, "Metals = %e ", zsum*MassUnits);
       fprintf(stdout, " momenta magnitude = %e ", sqrt(psqsum));

       fprintf(stdout, " momenta error = %e ", psum*MomentaUnits);
       fprintf(stdout, " KE deposit = %e", kesum);
       fprintf(stdout, " Gas energy = %e ", gesum * EnergyUnits);
       fprintf(stdout, " TE = %e\n", tesum*EnergyUnits);
        /* Break out if something went wrong */
        if (isnan(dsum) || isnan(zsum) || isnan(psqsum)|| isnan(tesum)){
           fprintf(stdout, "MechStars_depositFeedback [370]: Found a nan: %e %f %e %e\n",dsum, zsum, psqsum, tesum);
            ENZO_FAIL("MechStars_depositFeedback NaN in grid field!");
        }
    }
    if (criticalDebug){
        for (int i = 0; i< size ; i++){
            postMass += BaryonField[DensNum][i];
            postZ += BaryonField[MetalNum][i];
            postP += BaryonField[Vel1Num][i]+BaryonField[Vel2Num][i]+BaryonField[Vel3Num][i];
            postPmag += pow(BaryonField[Vel1Num][i]*MomentaUnits,2)+
                pow(BaryonField[Vel2Num][i]*MomentaUnits,2)
                +pow(BaryonField[Vel3Num][i]*MomentaUnits,2);
            postTE += BaryonField[TENum][i];
            postGE += BaryonField[GENum][i];
        }
        fprintf(stderr, "Difference quantities: dxRatio = %f dMass = %e dZ = %e  P = %e |P| = %e TE = %e GE = %e coupledGE = %e Ej = %e\n",
                dxRatio, (postMass-preMass)*MassUnits, (postZ-preZ)*MassUnits, 
                (postP - preP)*MomentaUnits,
                (sqrt(postPmag) - sqrt(prePmag)),
                 (postTE-preTE)*EnergyUnits, (postGE-preGE)*EnergyUnits,
                coupledGasEnergy*EnergyUnits*nCouple, ejectaEnergy);
        if(isnan(postMass) || isnan(postTE) || isnan(postPmag) || isnan(postZ)){
            fprintf(stderr, "NAN IN GRID: %e %e %e %e\n", postMass, postTE, postZ, postP);
            ENZO_FAIL("MechStars_depositFeedback.C: 395\n")
        if (postGE-preGE < 0.0)
            ENZO_FAIL("471");
        }
    }

    /* Transform the grid back */

    transformComovingWithStar(BaryonField[DensNum], BaryonField[MetalNum], 
                        BaryonField[Vel1Num],BaryonField[Vel2Num],
                        BaryonField[Vel3Num],*up, *vp, *wp, 
                        GridDimension[0], GridDimension[1],
                        GridDimension[2], -1);


    delete [] u,v,w, gasEnergy, totalEnergy;
    return SUCCESS;
}