#define DEBUG 0
/***********************************************************************
  /
  /  GRID CLASS (WALK PHOTON PACKAGES ACROSS GRID)
  /
  /  written by: Tom Abel
  /  date:       August, 2003
  /  modified1:
  /
  /  PURPOSE: This is the heart of the radiative transfer algorithm.
  /    All the work is done here. Trace particles, split them, compute
  /    photo and heating rates, communicate PhotonPackages if they are on
  /    the same processor
  /
  /  RETURNS: FAIL or SUCCESS
  /
 ************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "ExternalBoundary.h"
#include "Fluxes.h"
#include "GridList.h"
#include "Grid.h"
#include "CosmologyParameters.h"
#include "RadiativeTransferHealpixRoutines.h"

#ifdef CONFIG_BFLOAT_4
#define ROUNDOFF 1e-6
#endif
#ifdef CONFIG_BFLOAT_8
#define ROUNDOFF 1e-12
#endif
#ifdef CONFIG_BFLOAT_16
#define ROUNDOFF 1e-16
#endif

#define MAX_HEALPIX_LEVEL 13
#define MAX_COLUMN_DENSITY 1e25

int FindRootGrid(int &dummy, grid **Grids0, int nGrids0, 
        FLOAT rx, FLOAT ry, FLOAT rz, FLOAT ux, FLOAT uy, FLOAT uz);
int SplitPhotonPackage(PhotonPackageEntry *PP);
FLOAT FindCrossSection(int type, float energy);
float ComputeInterpolatedValue(float *vc, int vci, int vcj, int vck, 
        float mx, float my, float mz);

const float erg_eV = 1.602176e-12;
const float c_cgs = 2.99792e10;
const float EnergyThresholds[] = {13.6, 24.6, 54.4, 11.2};
const double PopulationFractions[] = {1., 0.25, 0.25, 1.};
enum species {
    iHI, iHeI, iHeII, iH2I, iHII
};

int grid::PhotonBoundaryConditions(grid **Grids0, int nGrids0, FLOAT *r, const FLOAT *u, PhotonPackageEntry *PP,
        grid* &MoveToGrid, int &DeltaLevel, const FLOAT *DomainWidth, int &DeleteMe, grid *ParentGrid) {
    int dummy, i;
    if (    (r[0] <= GridLeftEdge[0]) || (GridRightEdge[0] <= r[0])  ||
            (r[1] <= GridLeftEdge[1]) || (GridRightEdge[1] <= r[1])  ||
            (r[2] <= GridLeftEdge[2]) || (GridRightEdge[2] <= r[2])) {
        switch (GravityBoundaryType) {
            case TopGridPeriodic: 
                FindRootGrid(dummy, Grids0, nGrids0, r[0], r[1], r[2], u[0], u[1], u[2]);
                if (dummy >= 0) {
                    MoveToGrid = Grids0[dummy];
                    DeltaLevel = 0;
                    PP->Radius += ROUNDOFF;
                } else {
                    // Wrap the photon around the boundary
                    if (RadiativeTransferPeriodicBoundary) {
                        for (i=0; i<GridRank; i++) {
                            if (r[i] < DomainLeftEdge[i]) {
                                PP->SourcePosition[i] += DomainWidth[i];
                                r[i] += DomainWidth[i];
                            } else if (r[i] > DomainRightEdge[i]) {
                                PP->SourcePosition[i] -= DomainWidth[i];
                                r[i] -= DomainWidth[i];
                            }
                        }
                        FindRootGrid(dummy, Grids0, nGrids0, r[0], r[1], r[2], u[0], u[1], u[2]);
                        MoveToGrid = Grids0[dummy];
                        DeltaLevel = 0;
                        PP->Radius += ROUNDOFF;
                    } else {
                        // PhotonPackage left the box
                        PP->Photons=-1;
                        DeleteMe = TRUE;
                    }
                }
                return 1; // WalkPhotonPackage should return success
            case TopGridIsolated: 
                FindRootGrid(dummy, Grids0, nGrids0, r[0], r[1], r[2], u[0], u[1], u[2]);
                if (dummy >= 0) {
                    MoveToGrid = Grids0[dummy];
                    DeltaLevel = 0;
                    PP->Radius += ROUNDOFF;
                } else {
                    // PhotonPackage left the box
                    PP->Photons=-1;
                    DeleteMe = TRUE;
                }
                return 1; // WalkPhotonPackage should return success
            case SubGridIsolated:
                MoveToGrid = ParentGrid;
                DeltaLevel = -1;
                // PP->Radius += 0.01*CellWidth[0][0];
                PP->Radius += ROUNDOFF;
                if (DEBUG) 
                    fprintf(stdout, "Walk: left grid: sent photon to grid %x\n", ParentGrid);
                return 1; // WalkPhotonPackage should return success
            case GravityUndefined:
            default:
                fprintf(stdout, "grid::WalkPhotonPackage: "
                        "GravityBoundaryType = RadiationBoundary undefined %"ISYM".\n",
                        GravityBoundaryType);
                return 2; // WalkPhotonPackage should fail
        }
    }
    return 0; //WalkPhotonPackage should continue
}

void grid::PhotonDoAbsorption(int i, float* field, int kphNum, int gammaNum, int index, float conversion, FLOAT sigma, FLOAT dr, FLOAT &photons, FLOAT rate_kph, FLOAT rate_gamma) {
    FLOAT dP, n, tau;
    // Get absorber density -- no interpolation
    n = PopulationFractions[i]*field[index]*conversion;
    tau = n*sigma*dr; // optical depth of ray segment
    // at most use all photons for photo-ionizations
    if (tau > 2.e1) dP = (1.0+ROUNDOFF) * photons;
    else if (tau > 1.e-4) dP = min(photons*(1-expf(-tau)), photons);
    else dP = min(photons*tau, photons);
    // contributions to the photoionization rate is over whole timestep
    BaryonField[kphNum][index] += dP*rate_kph;
    // the heating rate is just the number of photo ionizations times
    // the excess energy units here are  eV/s/cm^3 *TimeUnits. 
    BaryonField[gammaNum][index] += dP*rate_gamma;
    photons -= dP;
}

double grid::PhotonRadiationPressureConversion(float DensityUnits, float VelocityUnits, FLOAT CellVolume) {
    /* Calculate conversion factor for radiation pressure.  In cgs, we
       actually calculate acceleration due to radiation pressure, (all
       unit conversions are in []),
       dA = dMomentum / Mass 
       = (N_absorbed * Energy / c) / (CellVolume * Density) * r_hat
       ---  N_absorbed = dP * [L^3] / [t]
       dA = (dP * [L^3] / [t] * Energy / c) / (CellVolume * Density) * r_hat
       LHS = [~] * [v] / [t]
       RHS = [L^3] / [t] * (erg_eV) / (c_cgs) / [L^3] / [rho]
       (unit conversion for acceleration will be [~])
       [~] = RHS / ([v] / [t])
       = (erg_eV / c_cgs) / ([t] * [rho]) / ([v] / [t])
       = (erg_eV / c_cgs) / [rho] / [v]
       Therefore,
       dA = [~] * dP * Energy / (CellVolume * Density) * r_hat
       Since CellVolume is constant for the grid, incorporate it into [~].
       dA = [~] * dP * Energy / Density * r_hat
     */
    return erg_eV / c_cgs / DensityUnits / VelocityUnits / CellVolume;
}

int grid::WalkPhotonPackage(PhotonPackageEntry **PP, 
        grid **MoveToGrid, grid *ParentGrid, grid *CurrentGrid, 
        grid **Grids0, int nGrids0, int DensNum, 
        int HINum, int HeINum, int HeIINum, int H2INum, 
        int kphHINum, int gammaHINum, int kphHeINum, 
        int gammaHeINum, int kphHeIINum, int gammaHeIINum, 
        int kdissH2INum, int RPresNum1, int RPresNum2, 
        int RPresNum3, int &DeleteMe, int &PauseMe, 
        int &DeltaLevel,
        float DensityUnits, float TemperatureUnits,
        float VelocityUnits, float LengthUnits,
        float TimeUnits) {

    // Check for early termination
    if ((*PP)->Photons <= 0) {
        (*PP)->Photons=-1;
        DeleteMe = TRUE;
        if (DEBUG)
            fprintf(stdout, "called WalkPhotonPackge with empty PhotonPackage "
                    "%x %x %x %x %"GSYM"\n",  (*PP), 
                    (*PP)->PreviousPackage, 
                    (*PP)->NextPackage,  PhotonPackages, 
                    -1);
        return SUCCESS;
    }
    if (((*PP) == NULL) || ((*PP)->PreviousPackage->NextPackage != (*PP))) {
        fprintf(stdout, "called WalkPhotonPackge with invalid pointer "
                "%x %x %x %x\n",  (*PP), 
                (*PP)->PreviousPackage, 
                (*PP)->PreviousPackage->NextPackage, 
                PhotonPackages);
        ENZO_FAIL("");
    }

    float ConvertToProperNumberDensity = DensityUnits/1.673e-24f;
    FLOAT DomainWidth[GridRank];
    FLOAT dx = CellWidth[0][0], dx2 = dx*dx;
    float m[GridRank], slice_factor, energy = (*PP)->Energy;
    FLOAT CellVolume=1, r[GridRank], radius, oldr, prev_radius;
    // This controls the splitting condition, where this many rays must exist in each cell
    FLOAT SplitCriteron = dx2 / RadiativeTransferRaysPerCell;
    FLOAT SplitWithinRadius = (RadiativeTransferSplitPhotonRadius > 0) ?
        RadiativeTransferSplitPhotonRadius * (3.086e21 / LengthUnits) : 2.0;
    FLOAT SplitCriteronIonized = dx2, PauseRadius = huge_number;
    // Convert escape fraction radius into code units
    float escapeRadius = RadiativeTransferPhotonEscapeRadius * (3.086e21 / LengthUnits);
    float PhotonEscapeRadius[3] = {.5 * escapeRadius, escapeRadius, 2.0 * escapeRadius};

    FLOAT photons = (*PP)->Photons;
    int i=0, index=0, type = (*PP)->Type;
    int kphNum[3] = {kphHINum, kphHeINum, kphHeIINum}, gammaNum[3] = {gammaHINum, gammaHeINum, gammaHeIINum};
    DeltaLevel = 0;
    FLOAT dir_vec[GridRank];
    if (pix2vec_nest((long) (1 << (*PP)->level), (*PP)->ipix, dir_vec)==FAIL) {
        fprintf(stdout,"grid::WalkPhotonPackage:  pix2vec_nest outor %"ISYM" %"ISYM" %"GSYM" %x\n",
                (long) (1 << (*PP)->level), (*PP)->ipix, photons, 
                (*PP) );
        (*PP)->Photons=0;
        ENZO_FAIL("");
    }

    // Cell geometry
    for (i=0; i<GridRank; i++) {
        DomainWidth[i] = DomainRightEdge[i] - DomainLeftEdge[i];
        CellVolume *= CellWidth[i][0];
    }

    /* Compute the photon distance that corresponds to a distance =
       R_merge away from the super source. 
       | (PauseRadius * dir_vec) + (vector_between_source_and_super) | = R_merge
       solve for PauseRadius, which results in a quadratic equation.

       PauseRadius = -u1 +/- sqrt(u1^2 - d^2 + R_merge^2), where
       u1 = dir_vec \dot d
       d := vector between current and super source
     */

    if (RadiativeTransferSourceClustering && (*PP)->CurrentSource != NULL) {
        FLOAT r_merge = 2*RadiativeTransferPhotonMergeRadius * (*PP)->CurrentSource->ClusteringRadius;
        FLOAT d2_ss = 0.0, u_dot_d = 0.0;
        for (i=0; i<GridRank; i++) {
            FLOAT d_ss = (*PP)->SourcePosition[i] - (*PP)->CurrentSource->Position[i];
            d2_ss += d_ss * d_ss;
            u_dot_d += dir_vec[i] * d_ss;
        }
        FLOAT sqrt_term = sqrt(u_dot_d*u_dot_d - d2_ss + r_merge*r_merge);
        PauseRadius = -u_dot_d + (sqrt_term > u_dot_d) ? sqrt_term : -sqrt_term;
    }

    // solid angle associated with package (= 4 Pi/N_package[on this level]) 
    FLOAT dP, EndTime = PhotonTime+dtPhoton, omega_package=4*M_PI/(12.0*pow(4.0,(*PP)->level));
    float dtheta = sqrt(omega_package);
    int direction, keep_walking=1, dummy, g[GridRank], count=0, u_dir[GridRank];
    FLOAT s[3] = {(*PP)->SourcePosition[0], (*PP)->SourcePosition[1], (*PP)->SourcePosition[2]};
    FLOAT u[3] = {dir_vec[0], dir_vec[1], dir_vec[2]};
    //if (u[0] == 0) u[0] = ROUNDOFF;
    if (fabs(u[1]) < ROUNDOFF) u[1] = sign(u[1])*ROUNDOFF; // zeros in y direction possible
    if (fabs(u[2]) < ROUNDOFF) u[2] = sign(u[2])*ROUNDOFF; // zeros in z direction possible

    for (i=0; i<GridRank; i++) {
        u_dir[i] = (sign(u[i])+1)/2; // Ray direction
        r[i] = s[i] + (*PP)->Radius * u[i]; // My Position in coordinates [0..1] 		     
        g[i] = (int)((r[i]-GridLeftEdge[i])/CellWidth[i][0]); // current cell
        // on cell boundaries the index will change in negative directions
        if (r[i] == GridLeftEdge[i] + (FLOAT)g[i] * CellWidth[i][0])
            g[i] += (sign(u[i])-1)/2;
    }
    // speed of light in code units. note this one is independent of a(t)
    // Modify the photon propagation speed by this parameter
    double c = c_cgs/VelocityUnits*RadiativeTransferPropagationSpeedFraction, RadiationPressureConversion;

    FLOAT sigma[4], factor2[3], emission_dt_inv = 1.0 / (*PP)->EmissionTimeInterval;
    int nSecondaryHII, nSecondaryHeIII;
    float *density = BaryonField[DensNum];
    float *fields[5] = { BaryonField[HINum], BaryonField[HeINum], BaryonField[HeIINum], 
        MultiSpecies > 1 ? BaryonField[H2INum] : 0, BaryonField[HINum+1]};

    // find relevant cross-sections premultiplied by the unit conversion factor for tau
    if (type == iHI || type == iHeI || type == iHeII || type == iH2I)
        sigma[type] = (*PP)->CrossSection * LengthUnits;
    else if (type == 4) {
        // For X-ray photons, we do heating and ionization for HI/HeI/HeII in one shot.
        for (i=0; i<3; i++) {
            sigma[i]   = FindCrossSection(i, energy) * LengthUnits;
            factor2[i]   = emission_dt_inv * (energy - EnergyThresholds[i]);
        }
        // Secondary ionizations from X-rays
        nSecondaryHII = (int) floor(energy / EnergyThresholds[iHI] - 1);
        nSecondaryHeIII = (int) floor(energy / EnergyThresholds[iHeII] - 1);
    }

    RadiationPressureConversion = PhotonRadiationPressureConversion(DensityUnits, VelocityUnits, CellVolume);

    while (keep_walking) {
        index = GRIDINDEX(g[0],g[1],g[2]);
        oldr = (*PP)->Radius;
        FLOAT min_dr = 1e20, dr;

        switch (PhotonBoundaryConditions(Grids0, nGrids0, r, u, *PP, *MoveToGrid, DeltaLevel, DomainWidth, DeleteMe, ParentGrid)) {
            case 1:
                return SUCCESS;
            case 2:
                ENZO_FAIL("");
        }

        // closest cell edge crossing radius
        for (i=0; i<GridRank; i++) {
            dr = (GridLeftEdge[i] + ((FLOAT) (g[i] + u_dir[i])) * CellWidth[i][0] - s[i]) / u[i];
            if (dr < min_dr) {
                direction = i;
                min_dr = dr;
            }
        }

        radius = min_dr + ROUNDOFF;
        dr = radius - oldr;

        if (dr < 0) {
            printf("dr < 0:   %"GSYM" %"GSYM" %"GSYM"\n", dr, min_dr, oldr);
            (*PP)->Photons = -1;
            DeleteMe = TRUE;
            return SUCCESS;
        }

        for (i=0; i<GridRank; i++)
            r[i] = s[i] + radius*u[i]; // My Position in coordinates [0..1]	     

        if (SubgridMarker[index] != CurrentGrid) {
            if (SubgridMarker[index] == NULL) {
                (*MoveToGrid) = ParentGrid;
                DeltaLevel = -1;
            } else {
                (*MoveToGrid) = SubgridMarker[index];
                DeltaLevel = 1;
                if (DEBUG) 
                    printf("different grid subgrid marker %x %x %"ISYM" %"ISYM" %"ISYM
                            " %"ISYM" %"ISYM" %"FSYM" %"FSYM" %"FSYM"\n",
                            SubgridMarker[index], CurrentGrid, index, count,
                            g[0], g[0], g[0], r[0], r[1], r[2]);
            }
            // move it at least a tiny fraction of the grid cell to not have
            // to worry about round off errors shifting photons back and
            // forth between two grids without doing anything (*PP)->Radius
            // *= (1. + 0.001 * CellWidth[0][0]);
            //      (*PP)->Radius += 0.01*CellWidth[0][0];
            (*PP)->Radius += ROUNDOFF;
            return SUCCESS;
        }

        // splitting condition split package if its associated area is
        // larger than dx^2/RadiativeTransferRaysPerCell but only as long as the radius is
        // less than one box length.  Don't split beyond SplitWithinRadius
        // if the radiation is optically thin (Xray, LW)

        // One ray per cell in ionized cells
        //    if (fields[iHII][index] > 0.5*CoolData.HydrogenFractionByMass*density[index])
        //      splitMe = omega_package*radius*radius > SplitCriteronIonized;
        //    else
        float solid_angle = radius * radius * omega_package;

        if (solid_angle > SplitCriteron && radius < SplitWithinRadius && (*PP)->level < MAX_HEALPIX_LEVEL) {
            // split the package and discontinue parent ray
            int return_value = SplitPhotonPackage(*PP);
            (*PP)->Photons = -1;
            DeleteMe = TRUE;
            NumberOfPhotonPackages += 4;
            return return_value;
        }

        if (DEBUG > 1) 
            fprintf(stdout, "%x %"ISYM" %"ISYM" %"ISYM" %"GSYM" %"GSYM"\t|\n",
                    (*PP), g[0], g[1], g[2], (*PP)->Radius, dr);

        // nor do we want transport longer than the grid timestep
        dr = min(dr, c*(EndTime-(*PP)->CurrentTime));
        FLOAT cdt = dr / c;

        // Check for ray merging, only consider a fraction of the ray to
        // make radius=PauseRadius and return.
        float fraction;
        if ((*PP)->Radius+dr > PauseRadius) {
            fraction = (PauseRadius-(*PP)->Radius) / dr;
            //fraction = min(fraction,0.1);
            //fraction = 1.0;
            dr *= fraction;
            cdt *= fraction;
            if (DEBUG > 1) 
                fprintf(stderr, "PAUSE: PP->Photons: %"GSYM"  PP->Radius: %"GSYM"\n",
                        photons, (*PP)->Radius);
            PauseMe = TRUE;
        }

        for (i=0; i<GridRank; i++)
            m[i] = fabs(s[i] + (oldr + 0.5*dr - ROUNDOFF) * u[i] - (GridLeftEdge[i] + (g[i]+0.5)*CellWidth[i][0]));
        float nearest_edge = max(max(m[0], m[1]), m[2]);
        slice_factor = min(0.5 + (0.5*dx-nearest_edge) / (dtheta*radius), 1);

        // If requested, keep track of hydrogen ionizing photons passing certain radii (for photon escape fractions)
        if (RadiativeTransferPhotonEscapeRadius > 0 && type == iHI)
            for (i = 0; i < 3; i++)
                if (radius > PhotonEscapeRadius[i] && oldr < PhotonEscapeRadius[i])
                    EscapedPhotonCount[i+1] += photons;

        switch (type) {
            case iHI: // HI/HeI/HeII ionization and heating
            case iHeI:
            case iHeII:
                PhotonDoAbsorption(type, fields[type], kphNum[type], gammaNum[type], index, ConvertToProperNumberDensity,
                        sigma[type], dr, photons, emission_dt_inv*slice_factor*slice_factor, 
                        emission_dt_inv*slice_factor*slice_factor*(energy-EnergyThresholds[type]));
                break;
            case iH2I:  // Lyman-Werner
                // We treat H2 dissociation with the shielding function from Draine & Bertoldi (1996)
                FLOAT n = PopulationFractions[type]*fields[type][index]*ConvertToProperNumberDensity;
                FLOAT oldColumnDensity = (*PP)->ColumnDensity, newColumnDensity = oldColumnDensity + n * dr * LengthUnits;
                if (oldColumnDensity < 1e14)
                    dP = (newColumnDensity < 1e14) ? 0 : photons * (1 - pow(newColumnDensity / 1e14, -0.75));
                else
                    dP = photons * (1 - pow(newColumnDensity / oldColumnDensity, -0.75));
                (*PP)->ColumnDensity = newColumnDensity;
                BaryonField[kdissH2INum][index] += photons * (sigma[type] / dx2) * emission_dt_inv * (dr * solid_angle / CellVolume);
                photons -= dP;
                break;
            case 4: // X-rays (HI/HeI/HeII all in one!)
                // Secondary ionizations (minus heating)
                double heat_factor = 1.f, kph_factor[3] = {1.f, 1.f, 1.f};
                if (RadiationXRaySecondaryIon) {
                    float xx = max(fields[iHII][index] / (fields[iHI][index] + fields[iHII][index]), 1e-4);
                    heat_factor       = 0.9971 * (1 - powf(1 - powf(xx, 0.2663f), 1.3163));
                    kph_factor[iHI]   = 0.3908 * powf(1 - powf(xx, 0.4092f), 1.7592f) * nSecondaryHII;
                    kph_factor[iHeII] = 0.0554 * powf(1 - powf(xx, 0.4614f), 1.6660f) * nSecondaryHeIII;
                }
                for (i=iHI; i<=iHeII; i++)
                    PhotonDoAbsorption(i, fields[i], kphNum[i], gammaNum[i], index, ConvertToProperNumberDensity,
                            sigma[i], dr, photons, emission_dt_inv*kph_factor[i], heat_factor*factor2[i]);
                break;
        }

        dP = (*PP)->Photons - photons;
        // acceleration due to radiation pressure: dA = [~] * dP * Energy / Density * r_hat
        if (RadiationPressure && (*PP)->Radius > (*PP)->SourcePositionDiff) {
            double kick = RadiationPressureConversion * dP * energy / density[index];
            for (i = 0; i < MAX_DIMENSION; i++)
                BaryonField[RPresNum1+i][index] += kick * dir_vec[i];
        }
        (*PP)->Photons = photons;
        // Indicator field to count photons per cell:
        BaryonField[kphNum[iHeII]][index] += 1;

        (*PP)->CurrentTime += cdt;
        (*PP)->Radius      += dr;

        // return in case we're pausing to merge or we're done
        if (PauseMe || (*PP)->CurrentTime >= EndTime)
            return SUCCESS;

        // return in case we're out of photons
        if (photons < tiny_number) {
            if (DEBUG>1)
                fprintf(stderr, "PP-Photons: %"GSYM"  PP->Radius: %"GSYM" PP->CurrentTime: \n"
                        "\tdP: %"GSYM"\tddr: %"GSYM"\t cdt: %"GSYM"\n",
                        photons, (*PP)->Radius, (*PP)->CurrentTime, dP, dr, cdt);
            (*PP)->Photons = -1;
            DeleteMe = TRUE;
            return SUCCESS;
        }

        count++;
        g[direction] += sign(u[direction]);
    } // while keep walking

    return SUCCESS;
}
