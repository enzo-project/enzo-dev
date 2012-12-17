/***********************************************************************
/
/  GRID CLASS (INITIALIZE HYDRO/MHD/RADHYDRO TURBULENCE SIMULATION)
/
/  written by: Peng Wang
/  date:       September, 2007
/  modified1: Tom Abel, allow for parallel generations of ICs (Sept 2009)
/
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"
#include "EOS.h"
#include "phys_constants.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
void Turbulence_Generator(float **vel, int dim0, int dim1, int dim2, int ind, 
			  float kmin, float kmax, float dk,
			  FLOAT **LeftEdge, FLOAT **CellWidth, int seed);

int grid::TurbulenceInitializeGrid(float CloudDensity, float CloudSoundSpeed, FLOAT CloudRadius, 
				   float CloudMachNumber, float CloudAngularVelocity, float InitialBField,
				   int SetTurbulence, int CloudType, int TurbulenceSeed, int PutSink, int level,
				   int SetBaryonFields)
{

  /* declarations */

  int dim, i, j, k,l, m, n, field, sphere, size, igrid, activesize;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
    DINum, DIINum, HDINum,  kphHINum, gammaNum, kphHeINum,
    kphHeIINum, kdissH2INum, RPresNum1, RPresNum2, RPresNum3;
  int ColourNum;

  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  FieldType[NumberOfBaryonFields++] = Velocity2;
  FieldType[NumberOfBaryonFields++] = Velocity3;
  //  if(HydroMethod == Zeus_Hydro)
  //  FieldType[NumberOfBaryonFields++] = InternalEnergy;
  //else  
  FieldType[NumberOfBaryonFields++] = TotalEnergy;
  

  if (DualEnergyFormalism) {
    FieldType[NumberOfBaryonFields++] = InternalEnergy;
  }
  if (HydroMethod == MHD_RK) {
    FieldType[NumberOfBaryonFields++] = Bfield1;
    FieldType[NumberOfBaryonFields++] = Bfield2;
    FieldType[NumberOfBaryonFields++] = Bfield3;
    FieldType[NumberOfBaryonFields++] = PhiField;
  }
  if(UseDivergenceCleaning)
    FieldType[NumberOfBaryonFields++] = Phi_pField;
  
  if (MultiSpecies) {
    FieldType[DeNum    = NumberOfBaryonFields++] = ElectronDensity;
    FieldType[HINum    = NumberOfBaryonFields++] = HIDensity;
    FieldType[HIINum   = NumberOfBaryonFields++] = HIIDensity;
    FieldType[HeINum   = NumberOfBaryonFields++] = HeIDensity;
    FieldType[HeIINum  = NumberOfBaryonFields++] = HeIIDensity;
    FieldType[HeIIINum = NumberOfBaryonFields++] = HeIIIDensity;
    if (MultiSpecies > 1) {
      FieldType[HMNum    = NumberOfBaryonFields++] = HMDensity;
      FieldType[H2INum   = NumberOfBaryonFields++] = H2IDensity;
      FieldType[H2IINum  = NumberOfBaryonFields++] = H2IIDensity;
    }
    if (MultiSpecies > 2) {
      FieldType[DINum   = NumberOfBaryonFields++] = DIDensity;
      FieldType[DIINum  = NumberOfBaryonFields++] = DIIDensity;
      FieldType[HDINum  = NumberOfBaryonFields++] = HDIDensity;
    }
  }
  //  FieldType[ColourNum = NumberOfBaryonFields++] = Metallicity;

  if (RadiativeTransfer && (MultiSpecies < 1)) {
    fprintf(stderr, "Grid_PhotonTestInitialize: Radiative Transfer but not MultiSpecies set");
    return FAIL;
  }

#ifdef TRANSFER
  if (RadiativeTransfer) {
    if (MultiSpecies) {
      FieldType[kphHINum    = NumberOfBaryonFields++] = kphHI;
      FieldType[gammaNum    = NumberOfBaryonFields++] = PhotoGamma;
      FieldType[kphHeINum   = NumberOfBaryonFields++] = kphHeI;
      FieldType[kphHeIINum  = NumberOfBaryonFields++] = kphHeII;
      if (MultiSpecies > 1) {
        FieldType[kdissH2INum    = NumberOfBaryonFields++] = kdissH2I;
      }
    }
    
    if (RadiationPressure) {
      FieldType[RPresNum1 = NumberOfBaryonFields++] = RadPressure0;
      FieldType[RPresNum2 = NumberOfBaryonFields++] = RadPressure1;
      FieldType[RPresNum3 = NumberOfBaryonFields++] = RadPressure2;
    }
    
    NumberOfPhotonPackages = 0;
    PhotonPackages-> NextPackage= NULL;
    
  }
#endif

  int idrivex, idrivey, idrivez;
  if (UseDrivingField) {
    idrivex = NumberOfBaryonFields;
    idrivey = idrivex + 1;
    idrivez = idrivex + 2;
    FieldType[NumberOfBaryonFields++] = DrivingField1;
    FieldType[NumberOfBaryonFields++] = DrivingField2;
    FieldType[NumberOfBaryonFields++] = DrivingField3;
  }

  if (WritePotential) {
    FieldType[NumberOfBaryonFields++] = GravPotential;
    FieldType[NumberOfBaryonFields++] = AccelerationField1;
    FieldType[NumberOfBaryonFields++] = AccelerationField2;
    FieldType[NumberOfBaryonFields++] = AccelerationField3;
  }


  float DensityUnits = 1.0, LengthUnits = 1.0, TemperatureUnits = 1.0, TimeUnits = 1.0,
    VelocityUnits = 1.0;
  if (UsePhysicalUnit)
    GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, Time);
  double MassUnits = DensityUnits*pow(LengthUnits,3);
  printf("Mass Units = %"GSYM" \n",MassUnits);
  printf("Time Units = %"GSYM" \n",TimeUnits);
  printf("Density Units = %"GSYM" \n",DensityUnits);

  GravitationalConstant = 4.0*pi*GravConst*MassUnits*pow(TimeUnits,2)/pow(LengthUnits,3);

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }

  if (SetBaryonFields == 0) 
    return SUCCESS;


  size = 1;
  for (dim = 0; dim < GridRank; dim++) {
    size *= GridDimension[dim];
  }

  for (field = 0; field < NumberOfBaryonFields; field++) {
    if (BaryonField[field] == NULL) {
      BaryonField[field] = new float[size];
      for (n = 0; n < size; n++) {
	BaryonField[field][n] = 0.0;
      }
    }
  }

  /* Initialize radiation fields */
#ifdef TRANSFER
  if (this->InitializeRadiativeTransferFields() == FAIL) {
    fprintf(stderr, "\nError in InitializeRadiativeTransferFields.\n");
    return FAIL;
  }
#endif
  
  float *TurbulenceVelocity[3], *DrivingField[3];
  float CloudInternalEnergy, CloudPressure;
  activesize = 1;
  for (dim = 0; dim < GridRank; dim++) {
    activesize *= (GridDimension[dim] - 2*NumberOfGhostZones);
  }

  for (dim = 0; dim < GridRank; dim++) {
    TurbulenceVelocity[dim] = new float[size];
    DrivingField[dim] = new float[size];
    for (n = 0; n < activesize; n++) {
      TurbulenceVelocity[dim][n] = 0.0;
      DrivingField[dim][n] = 0.0;
    }
  }


  /* Set density, rotational speed, internal energy and chemical species. */

  /* assume isothermal initially */

  CloudPressure = CloudDensity * pow(CloudSoundSpeed,2);
  //  CloudInternalEnergy = pow(CloudSoundSpeed,2) / (Gamma - 1.0);
  float h, cs, dpdrho, dpde;
  EOS(CloudPressure, CloudDensity, CloudInternalEnergy, h, cs, dpdrho, dpde, EOSType, 1);

  float InitialFractionHII = 1.2e-5;
  float InitialFractionHeII = 1.0e-14;
  float InitialFractionHeIII = 1.0e-17;  

  /* if (CloudType == 3) {
    CloudRadius = 10.0;
    } */

  /* Cloud center is box center. */
  FLOAT xc = 0.5, yc = 0.5, zc = 0.5, xpos, ypos, zpos,
    cosphi, sinphi, x, y, z, r;
  FLOAT costheta = cos(17.5*M_PI/180.0), sintheta = sin(17.5*M_PI/180);
  float Density, eint, Velx, Vely, Velz;
  n = 0;
  for (k = 0; k < GridDimension[2]; k++) {
    for (j = 0; j < GridDimension[1]; j++) {
      for (i = 0; i < GridDimension[0]; i++, n++) {

	x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];
	r = sqrt(pow(x-xc,2) + pow(y-yc,2) + pow(z-zc,2));
	r = max(r, CellWidth[0][0]);

	xpos = x - xc;
	ypos = y - yc;
	zpos = z - zc;

	/* compute the azimuthal angle */
	cosphi = xpos/sqrt(xpos*xpos+ypos*ypos);
	sinphi = ypos/sqrt(xpos*xpos+ypos*ypos);

	Velx = Vely = Velz = 0.0;

	// default values:
	Density = CloudDensity;
	eint = CloudInternalEnergy;


	if (r < CloudRadius) {

	  /* Type 0: uniform cloud, 7: uniform cloud (only turb k=1-2) */

	  if (CloudType == 0 || CloudType == 7) {
	    Density = CloudDensity;
	    eint = CloudInternalEnergy;
	  }

	  /* Type 1: cloud density profile is the same as that of Nakamura & Li, ApJ, 662, 395 */

	  if (CloudType == 1) {
	    Density = CloudDensity / (1.0 + pow(3.0*r/CloudRadius,2));
	    eint = CloudInternalEnergy;
	    Velx = -CloudAngularVelocity * ypos;
	    Vely =  CloudAngularVelocity * xpos;
	    /*else {
	      Velx = -CloudAngularVelocity * (CloudRadius-r) * 1.5 / CloudRadius * ypos;
	      Vely =  CloudAngularVelocity * (CloudRadius-r) * 1.5 / CloudRadius * xpos;
	      }*/
	  }

	  /* Type 2: cloud density is uniform sphere embeded in low density medium */

	  if (CloudType == 2) {
	    Density = CloudDensity;
	    eint = CloudInternalEnergy;
	    Velx = -CloudAngularVelocity * ypos;
	    Vely =  CloudAngularVelocity * xpos;
	  }

	  /* Type 3: flattened 1/r^2 profile without ambient medium. */

	  if (CloudType == 3) {
	    Density = CloudDensity / (1.0 + pow(6.0*r/CloudRadius,2));
	    eint = CloudInternalEnergy;
	  }
 
	  /* Type 4: 1/r^2 profile with a smaller core.
	     This is a model for massive star formation with a seed
	     protostar in the center */

	  if (CloudType == 4) {
	    Density = 4.2508525*CloudDensity / (1.0 + pow(9.0*r/CloudRadius,2));
	    eint = CloudInternalEnergy;
	    Velx = -CloudAngularVelocity * ypos;
	    Vely =  CloudAngularVelocity * xpos;
	  }

	  if (CloudType == 5) {
	    float drho = rand();
	    drho = 0.0;//drho/RAND_MAX;
	    if (r < 5) {
	      Density = CloudDensity*(1.0+0.1*drho);
	    } else {
	      Density = CloudDensity/pow(r/0.05, 2)*(1.0+0.1*drho);
	    }
	    eint = CloudInternalEnergy/(1.0+0.1*drho);
	  }

	  /* Type 6: flattened 1/r^2 profile with large core and with ambient medium. */

	  if (CloudType == 6) {
	    Density = 1.0522054*CloudDensity / (1.0 + pow(4.0*r/CloudRadius,2));
	    eint = CloudInternalEnergy;
	  }


	} else {

	  if (CloudType == 0 || CloudType == 7 ) {
	    Density = CloudDensity/100.0;
	    eint = CloudInternalEnergy*100.;
	  }

	  if (CloudType == 1) {
	    Density = CloudDensity/20.0;
	    eint = CloudInternalEnergy;
	  }

	  if (CloudType == 2) {
	    Density = CloudDensity / 10.0;
	    eint = CloudInternalEnergy * 10.0;
	  }

          if (CloudType ==3) {
	    Density = max(DensityUnits, CloudDensity/(1.0 + pow(6.0*r/CloudRadius,2)));
	    eint = CloudInternalEnergy;
	  }

	  if (CloudType == 4) {
	    Density = max(DensityUnits,0.5*4.25*CloudDensity/(1.0 + pow(9.0*r/CloudRadius,2)));
	    eint = CloudInternalEnergy*200.0; //400.0;
	  }


          if (CloudType ==6) {
	    //Density = max(DensityUnits, 0.5*CloudDensity/(1.0 + pow(4.0*r/CloudRadius,2)));
	    Density = 0.1*CloudDensity/(1.0 + pow(4.0,2));
	    eint = CloudInternalEnergy*100.0; //400.0;
	  }

	}

	BaryonField[iden ][n] = Density;
	//	BaryonField[ColourNum][n] = Density*0.018477;
	BaryonField[ivx  ][n] = Velx;
	BaryonField[ivy  ][n] = Vely;
	BaryonField[ivz  ][n] = Velz;
	BaryonField[ietot][n] = eint;
	if (HydroMethod != Zeus_Hydro)
	  BaryonField[ietot][n] += 0.5*(Velx*Velx + Vely*Vely + Velz*Velz);
	if (DualEnergyFormalism) {
	  BaryonField[ieint][n] = eint;
	}
	
	if(Velx != 0.0) 
	  printf("    PROBLEM!!!! eint = %g, Velx = %g, Vely = %g, Velz = %g \n", eint, Velx, Vely, Velz);

	if (HydroMethod == MHD_RK) {
	  BaryonField[iBx  ][n]  = 0.0;//-InitialBField*sinphi;
	  BaryonField[iBy  ][n]  = 0.0;//InitialBField*sintheta;
	  BaryonField[iBz  ][n]  = InitialBField;
	  BaryonField[iPhi ][n]  = 0.0;
	  BaryonField[ietot][n] += 0.5 * pow(InitialBField,2) / Density;
	}

	if (UseDrivingField) {
	  BaryonField[idrivex][n] = 0.0;
	  BaryonField[idrivey][n] = 0.0;
	  BaryonField[idrivez][n] = 0.0;
	}
 
	if (MultiSpecies > 0) {

	  BaryonField[HIINum][n] = InitialFractionHII *
	    CoolData.HydrogenFractionByMass * BaryonField[iden][n];
	  BaryonField[HeIINum][n] = InitialFractionHeII*
	    BaryonField[iden][n] * 4.0 * (1.0-CoolData.HydrogenFractionByMass);
	  BaryonField[HeIIINum][n] = InitialFractionHeIII*
	    BaryonField[iden][n] * 4.0 * (1.0-CoolData.HydrogenFractionByMass);
	  BaryonField[HeINum][n] =
	    (1.0 - CoolData.HydrogenFractionByMass)*BaryonField[iden][n] -
	    BaryonField[HeIINum][n] - BaryonField[HeIIINum][n];

          /*if (MultiSpecies > 1) {
            BaryonField[HMNum][n] = InitialFractionHM*
              BaryonField[HIINum][n]* pow(temperature,float(0.88));
            BaryonField[H2IINum][n] = InitialFractionH2II*
              2.0*BaryonField[HIINum][n]* pow(temperature,float(1.8));
            BaryonField[H2INum][n] = InitialFractionH2I*
              BaryonField[0][n]*CoolData.HydrogenFractionByMass*pow(301.0,5.1);
	      }*/
	  
	  BaryonField[HINum][n] =
	    CoolData.HydrogenFractionByMass*BaryonField[iden][n]
	    - BaryonField[HIINum][n];
	  
          //if (MultiSpecies > 1)
	  // BaryonField[HINum][n] -= BaryonField[HMNum][n]
	  //  + BaryonField[H2IINum][n]
	  //  + BaryonField[H2INum][n];

	  /* electron "density": n_e * m_p */
	  
	  BaryonField[DeNum][n] = BaryonField[HIINum][n] +
	    0.25*BaryonField[HeIINum][n] + 0.5*BaryonField[HeIIINum][n];
	  
          //if (MultiSpecies > 1)
	  //BaryonField[DeNum][n] += 0.5*BaryonField[H2IINum][n] -
	  //  BaryonField[HMNum][n];

	}

      }
    }
  }


  /* Initialize turbulent velocity field */

  if (SetTurbulence) {

    float k1, k2, dk;
    if (CloudType == 0) {
      k1 = 2.0;
      k2 = 10.0;
      dk = 1.0;
    }
    if (CloudType == 1) {
      k1 = 2, k2 = 4, dk = 1;
      k1 = int(1.0/CloudRadius);
      dk = int(0.5/CloudRadius);
      k2 = k1 + 8.0*dk;
    }
    if (CloudType == 2) {
      k1 = 2.;
      k2 = 37;
      dk = 5;
    }
    if (CloudType == 3) {
      k1 = 2.0;
      k2 = 10.0;
      dk = 1.0;
    }
    if (CloudType == 4 || CloudType == 6) {
      k1 = 2.0;
      k2 = min(34.0, int(GridDimension[0]/10));
      printf("                GridDimension[0] = %"ISYM"\n",GridDimension[0] );
      dk = max(1.0,int((k2-k1)/10));
    }
    if (CloudType == 7) {
      k1 = 1.0;
      k2 = 2.0;
      dk = 0.5;
    }


    printf("Begin generating turbulent velocity spectrum...\n");
    Turbulence_Generator(TurbulenceVelocity, GridDimension[0], 
			 GridDimension[1],
			 GridDimension[2],
			 4.0, k1, k2, dk,
			 CellLeftEdge, CellWidth, TurbulenceSeed);    
    /*Turbulence_Generator(TurbulenceVelocity, GridDimension[0]-2*NumberOfGhostZones, 
			 GridDimension[1]-2*NumberOfGhostZones,
			 GridDimension[2]-2*NumberOfGhostZones,
			 4.0, k1, k2, dk,
			 CellLeftEdge, CellWidth, TurbulenceSeed);    */

    printf("Turbulent spectrum generated\n");

    float VelocityNormalization = 1;
// for level > 0 grids the CloudMachNumber passed in is actuall the Velocity normalization factor
  if (level > 0) VelocityNormalization = CloudMachNumber; 
  printf("Cloud Mach Number = %"GSYM" \n",CloudMachNumber);
  for (i = 0; i < 3; i++) {
    for (n = 0; n < activesize; n++) {
      TurbulenceVelocity[i][n] *= VelocityNormalization;
    }
  }


    /* Set turbulent velocity field */

    n = 0;
    for (k = 0; k < GridDimension[2]; k++) {
      for (j = 0; j < GridDimension[1]; j++) {
	for (i = 0; i < GridDimension[0]; i++, n++) {
	  /*for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {*/
	  igrid = i + GridDimension[0]*(j+k*GridDimension[1]);
	  x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	  y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	  z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];
	  
	  r = sqrt(pow(fabs(x-xc),2)+pow(fabs(y-yc),2)+pow(fabs(z-zc),2));
	  r = max(r, 0.1*CellWidth[0][0]);
	  
	  if (r < CloudRadius) {
	    BaryonField[ivx][igrid] += TurbulenceVelocity[0][n];
	    BaryonField[ivy][igrid] += TurbulenceVelocity[1][n];
	    BaryonField[ivz][igrid] += TurbulenceVelocity[2][n];
	    BaryonField[ietot][igrid] += 
	      0.5 * (pow(TurbulenceVelocity[0][n],2) + 
		     pow(TurbulenceVelocity[1][n],2) + 
		     pow(TurbulenceVelocity[2][n],2));
	  }
	} 
      }
    }    

  }

  /* Initialize driving force field = efficiency * density * velocity / t_ff*/


  if (UseDrivingField) {
    float k1, k2, dk;
    k1 = 3.0;
    k2 = 4.0;
    dk = 1.0;
    printf("Begin generating driving force field ...\n");
    Turbulence_Generator(DrivingField, GridDimension[0]-2*NumberOfGhostZones, 
			 GridDimension[1]-2*NumberOfGhostZones,
			 GridDimension[2]-2*NumberOfGhostZones,
			 4.0, k1, k2, dk,
			 CellLeftEdge, CellWidth, TurbulenceSeed);
    printf("Driving force field generated\n");


    /* Renormalize to ensure <F>=0 */

    double Fx = 0.0, Fy = 0.0, Fz = 0.0;
    for (n = 0; n < activesize; n++) {
      Fx += DrivingField[0][n];
      Fy += DrivingField[1][n];
      Fz += DrivingField[2][n];
    }

    Fx /= activesize;
    Fy /= activesize;
    Fz /= activesize;
    
    for (n = 0; n < activesize; n++) {
      DrivingField[0][n] -= Fx;
      DrivingField[1][n] -= Fy;
      DrivingField[2][n] -= Fz;
    }
    
    /* Renormalize the mass-weighted 3D rms velocity inside the cloud */
    /*
    double VelRMS = 0.0, Mass = 0.0;
    n = 0;
    for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {
	  igrid = i + GridDimension[0]*(j+k*GridDimension[1]);
	  x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	  y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	  z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];
	  
	  r = sqrt(pow(fabs(x-xc),2)+pow(fabs(y-yc),2)+pow(fabs(z-zc),2));
	  r = max(r, 0.1*CellWidth[0][0]);
	  
	  if (r <= CloudRadius) {
	    VelRMS += BaryonField[iden][igrid] * 
	      sqrt(pow(DrivingField[0][n],2) + pow(DrivingField[1][n],2) + pow(DrivingField[2][n],2));
	    Mass += BaryonField[iden][igrid];
	  }
	  
	}
      }
    }
    */
    //printf("Grid_TubInit: Mass = %"FSYM"\n",Mass);
    /* VelRMS /= Mass;
    double t_ff = sqrt(32.0/(3.0*M_PI*CloudDensity));
    double NormFactor = CloudMachNumber * CloudSoundSpeed / VelRMS / t_ff;
    for (dim = 0; i < GridRank; dim++) {
      for (n = 0; n < activesize; n++) {
	DrivingField[dim][n] *= NormFactor;
      }
      }*/

    n = 0;
    for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {
	  igrid = i + GridDimension[0]*(j+k*GridDimension[1]);
	  BaryonField[idrivex][igrid] = DrivingField[0][n];
	  BaryonField[idrivey][igrid] = DrivingField[1][n];
	  BaryonField[idrivez][igrid] = DrivingField[2][n];
	}
      }
    }

  }

  for (dim = 0; dim < GridRank; dim++) {
    delete [] TurbulenceVelocity[dim];
    delete [] DrivingField[dim];
  }    

  /* Put a sink particle if we are studying massive star formation */

  //  if (PutSink == 1 && level == MaximumRefinementLevel) {
  if (PutSink == 1 && level == 0) {  // set it up on level zero and make it mustrefine

    //    double mass_p = 20.0*1.989e33;
    double mass_p = 3.415*1.989e33;
    mass_p /= MassUnits;
    double dx = CellWidth[0][0];
    double den_p = mass_p / pow(dx,3);
    double t_dyn = sqrt(3*M_PI/(6.672e-8*den_p*DensityUnits));
    t_dyn /= TimeUnits;

    double dxm = dx / pow(RefineBy, MaximumRefinementLevel);

    printf("Adding a Sink Particle. \n");

    NumberOfParticles = 1;
    NumberOfStars = 1;    
    NumberOfParticleAttributes = 3;
    //    MaximumParticleNumber = 1;
    if (StellarWindFeedback) NumberOfParticleAttributes = 6;
    this->AllocateNewParticles(NumberOfParticles);

    ParticleMass[0] = den_p;
    ParticleNumber[0] = 0;
    ParticleType[0] = PARTICLE_TYPE_MUST_REFINE;
    ParticlePosition[0][0] = 0.35+0.5*dxm;
    ParticlePosition[1][0] = 0.45+0.5*dxm;
    ParticlePosition[2][0] = 0.70+0.5*dxm;
    ParticleVelocity[0][0] = 0.0;
    ParticleVelocity[1][0] = 0.0;
    ParticleVelocity[2][0] = 0.0;

    ParticleAttribute[0][0] = 0.0; // creation time             
    ParticleAttribute[1][0] = t_dyn;
    ParticleAttribute[2][0] = mass_p;

    if (StellarWindFeedback) {
      ParticleAttribute[3][0] = 1.0;  
      ParticleAttribute[4][0] = 0.0;
      ParticleAttribute[5][0] = 0.0;
    }

    for (i = 0; i< MAX_DIMENSION+1; i++){
      ParticleAcceleration[i] = NULL;
    }
    this->ClearParticleAccelerations();

  }


  if (PutSink == 2 && level == 0) {  // set it up on level zero and make it mustrefine

    NumberOfParticles = 64;
    NumberOfStars = 64;

    printf("Adding Sink Particles. \n");
    NumberOfParticleAttributes = 3;
    if (StellarWindFeedback) NumberOfParticleAttributes = 6;
    this->AllocateNewParticles(NumberOfParticles);

    //    double mass_p = 20.0*1.989e33;
    double mass_m = 3.415*1.989e33; //Mass of massive stars
    double mass_s = 0.01*1.989e33; //Mass of small stars
    mass_m /= MassUnits;
    mass_s /= MassUnits;
    double dx = CellWidth[0][0];
    double den_m = mass_m / pow(dx,3);
    double den_s = mass_s / pow(dx,3);
    double t_dyn_m = sqrt(3*M_PI/(6.672e-8*den_m*DensityUnits));
    double t_dyn_s = sqrt(3*M_PI/(6.672e-8*den_s*DensityUnits));
    t_dyn_m /= TimeUnits;
    t_dyn_s /= TimeUnits;
    double dxm = dx / pow(RefineBy, MaximumRefinementLevel);

    //    MaximumParticleNumber = 1;


    for (k=0; k<4; k++){
      for (j=0; j<4; j++){
	for (i=0; i<4; i++){
	  l = i+4*j+16*k;
	  printf("Creating particle %i \n",l);
	  ParticleMass[l] = den_m;
	  ParticleNumber[l] = l;
	  ParticleType[l] = PARTICLE_TYPE_MUST_REFINE;
	  ParticlePosition[0][l] = 0.125+0.25*i+0.5*dxm;
	  ParticlePosition[1][l] = 0.125+0.25*j+0.5*dxm;
	  ParticlePosition[2][l] = 0.125+0.25*k+0.5*dxm;
	  ParticleVelocity[0][l] = 0.0;
	  ParticleVelocity[1][l] = 0.0;
	  ParticleVelocity[2][l] = 0.0;
	  ParticleAcceleration[0] = NULL;
	  ParticleAcceleration[1] = NULL;
	  ParticleAcceleration[2] = NULL;

	  ParticleAttribute[0][l] = 0.001; // creation time             
	  ParticleAttribute[1][l] = t_dyn_m; // t_dyn
	  ParticleAttribute[2][l] = mass_m;

	  if (StellarWindFeedback) {
	    ParticleAttribute[3][l] = 1.0;  
	    ParticleAttribute[4][l] = 0.0;
	    ParticleAttribute[5][l] = 0.0;
	  }
	  /*for (m = 0; m< MAX_DIMENSION+1; m++){
	    ParticleAcceleration[m][l] = NULL;
	    }*/

	  this->ClearParticleAccelerations();
	  printf("Completed particle %i, position %g,%g,%g \n",l,ParticlePosition[0][l],ParticlePosition[1][l],ParticlePosition[2][l]);
	  //printf("Domain Right Edge = %g %g %g, dx = %g\n",DomainRightEdge[0],DomainRightEdge[1],DomainRightEdge[2],dx);
	}
      }
    }

  }

  if (PutSink == 3 && level == 0) {  // set it up on level zero and make it mustrefine

    NumberOfParticles = 64;
    NumberOfStars = 64;

    printf("Adding Dummy Sink Particles. \n");
    NumberOfParticleAttributes = 3;
    if (StellarWindFeedback) NumberOfParticleAttributes = 6;
    this->AllocateNewParticles(NumberOfParticles);

    //    double mass_p = 20.0*1.989e33;
    double mass_m = 1.989e13; //Mass of massive stars
    mass_m /= MassUnits;
    double dx = CellWidth[0][0];
    double den_m = mass_m / pow(dx,3);
    double t_dyn_m = sqrt(3*M_PI/(6.672e-8*den_m*DensityUnits));
    t_dyn_m /= TimeUnits;
    double dxm = dx / pow(RefineBy, MaximumRefinementLevel);

    for (k=0; k<4; k++){
      for (j=0; j<4; j++){
	for (i=0; i<4; i++){
	  l = i+4*j+16*k;
	  printf("Creating particle %i \n",l);
	  ParticleMass[l] = den_m;
	  ParticleNumber[l] = l;
	  ParticleType[l] = PARTICLE_TYPE_MUST_REFINE;
	  ParticlePosition[0][l] = 1.125+0.25*i;//+0.5*dx;
	  ParticlePosition[1][l] = 1.125+0.25*j;//+0.5*dx;
	  ParticlePosition[2][l] = 1.125+0.25*k;//+0.5*dx;
	  ParticleVelocity[0][l] = 0.0;
	  ParticleVelocity[1][l] = 0.0;
	  ParticleVelocity[2][l] = 0.0;
	  ParticleAcceleration[0] = NULL;
	  ParticleAcceleration[1] = NULL;
	  ParticleAcceleration[2] = NULL;

	  ParticleAttribute[0][l] = 0.0; // creation time             
	  ParticleAttribute[1][l] = 0.0; //t_dyn_m; // t_dyn
	  ParticleAttribute[2][l] = mass_m;

	  if (StellarWindFeedback) {
	    ParticleAttribute[3][l] = 1.0;  
	    ParticleAttribute[4][l] = 0.0;
	    ParticleAttribute[5][l] = 0.0;
	  }
	  /*for (m = 0; m< MAX_DIMENSION+1; m++){
	    ParticleAcceleration[m][l] = NULL;
	    }*/

	  this->ClearParticleAccelerations();
	  printf("Completed particle %i, position %g,%g,%g \n",l,ParticlePosition[0][l],ParticlePosition[1][l],ParticlePosition[2][l]);
	  //printf("Domain Right Edge = %g %g %g, dx = %g\n",DomainRightEdge[0],DomainRightEdge[1],DomainRightEdge[2],dx);
	}
      }
    }

  }




  /*  printf("XXX PutSinkParticle = %"ISYM"\n", PutSinkParticle);
  int PutSinkParticle = 0;
  printf("XXX PutSinkParticle = %"ISYM"\n", PutSinkParticle);
  if (PutSinkParticle == 1 && level == 0) {
    NumberOfParticleAttributes = 6;
    double mass_p = 1.1*1.989e33;
    mass_p /= MassUnits;
    double dx = CellWidth[0][0];
    double den_p = mass_p / pow(dx,3);
    double t_dyn = sqrt(3*Pi/(6.672e-8*den_p*DensityUnits));
    t_dyn /= TimeUnits;

    NumberOfParticles = 1;
    NumberOfStarParticles = 1;
    MaximumParticleNumber = 1;
    this->AllocateNewParticles(NumberOfParticles);
    ParticleMass[0] = den_p;
    ParticleNumber[0] = 0;
    ParticleType[0] = PARTICLE_TYPE_MUST_REFINE;
    ParticlePosition[0][0] = 0.50;                         
    ParticlePosition[1][0] = 0.50;
    ParticlePosition[2][0] = 0.50;

    ParticleVelocity[0][0] = 0.0;
    ParticleVelocity[1][0] = 0.0;
    ParticleVelocity[2][0] = 0.0;
    LibbyTouched = 1;
    for (i = 0; i< MAX_DIMENSION+1; i++){
      ParticleAcceleration[i] = NULL;
    }
    this->ClearParticleAccelerations();

    ParticleAttribute[0][0] = 0.0; // creation time    
    ParticleAttribute[1][0] = t_dyn; // dynamical time                                                                
    ParticleAttribute[2][0] = mass_p; //                                                                                 
    printf("XXX Sink Particle in, NumberOfParticles = %"ISYM" \n",NumberOfParticles);
    }*/


  int TestMerge = 0;
  if (TestMerge == 1 && level == 0) {

    double mass_p = 0.05*1.989e33;
    mass_p /= MassUnits;
    double dx = CellWidth[0][0];
    double den_p = mass_p / pow(dx,3);
    double t_dyn = sqrt(3*M_PI/(6.672e-8*den_p*DensityUnits));
    t_dyn /= TimeUnits;

    double dxm = dx / pow(2.0, MaximumRefinementLevel);

    NumberOfParticles = 4;
    NumberOfStars = 4;
    //    MaximumParticleNumber = 4;
    this->AllocateNewParticles(NumberOfParticles);
    for (int i = 0; i < 4; i++) {
      ParticleMass[i] = den_p;
      ParticleNumber[i] = i;
      ParticleType[i] = PARTICLE_TYPE_MUST_REFINE;
      for (int dim = 0; dim < 3; dim++)
	ParticleVelocity[dim][i] = 0.0;
      ParticleAttribute[0][i] = 0.0;
      ParticleAttribute[1][i] = 0.0;
      ParticleAttribute[2][0] = mass_p;
    }
    
    //ParticleMass[0] = 10.0*den_p;

    ParticlePosition[0][0] = 0.501;
    ParticlePosition[1][0] = 0.501;
    ParticlePosition[2][0] = 0.49;

    ParticlePosition[0][1] = 0.501;
    ParticlePosition[1][1] = 0.501;
    ParticlePosition[2][1] = 0.51;

    ParticlePosition[0][2] = 0.501;
    ParticlePosition[1][2] = 0.501;
    ParticlePosition[2][2] = 0.92;

    ParticlePosition[0][3] = 0.501;
    ParticlePosition[1][3] = 0.501;
    ParticlePosition[2][3] = 0.53;


  }

  return SUCCESS;
}

