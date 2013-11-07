/***********************************************************************
/
/  GRID CLASS (SPLIT PARTICLES INTO CHILDREN PARTICLES)
/
/  written by: Ji-hoon Kim
/  date:       October, 2009
/  modified1:  John Regan October 2013
/              Converted from Fortran to C++
/
/  PURPOSE: This routine splits particles into 13 (=12+1) children particles 
/           when requested.  See Kitsionas & Whitworth (2002) for the
/           technical details of particle splitting,  which was already 
/           implemented and used in SPH/Gadget.
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "ErrorExceptions.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "CosmologyParameters.h"


int grid::CreateChildParticles(float dx, int NumberOfParticles, float *ParticleMass,
			       int *ParticleType, float *ParticlePosition[],
			       float *ParticleVelocity[], float *ParticleAttribute[],
			       int NumberAttributes, float *RefineRegionLeftEdge, 
			       float *RefineRegionRightEdge, float *CellLeftEdge[],
			       int *GridDimension, int MaximumNumberOfNewParticles,
			       int *NumberOfNewParticles)
			 
{
  int partnum = 0, i = 0, child = 0, m = 0, numpart = 0, innerchild = 0;
  int xindex=0, yindex=0, zindex=0;
  float rad = 0.0, sin60 = 0.0;
  float NewPos[3][CHILDRENPERPARENT];
  int ChildrenPerParent = CHILDRENPERPARENT, total_children = 0;
  float alpha[3];
  float l11 = 0.0, l12 = 0.0, l13 = 0.0;
  float l21 = 0.0, l22 = 0.0, l23 = 0.0;
  float l31 = 0.0, l32 = 0.0, l33 = 0.0;
  float  separation =  ParticleSplitterChildrenParticleSeparation;
  int iterations = ParticleSplitterIterations;
  
  /*
   *   The distance to the children particles rad is the same for all. It is 
   *   currently set to be dx (=CellWidth), but can be changed. 
   *   (i.e. ParticleSplitterChildrenParticleSeparation = 1 by default)
   */
  
  rad = dx * pow(separation, double(iterations));
  srand48(time(NULL));
  
  
  /* 
   * Loop over existing (parent) particles; It implicitly assumes that 
   * only DM and conventional star particles  get splitted.  Other 
   * particles - which usually become Star class particles - doesn't 
   * seem to have any reason to be splitted.  (as of Oct.2009)
   */
  
  for(partnum = 0; partnum < NumberOfParticles; partnum++)
    {
      if(ParticleMass[partnum] > 0.0 && ParticleType[partnum] <= 2)
	{
	  /* 
	   * Check that particle is within the most refined region.
	   * We only care about those particles. We're going to
	   * blatantly discriminate against non-refined types around 
	   * here.
	   */
	  if(ParticlePosition[0][partnum] < RefineRegionLeftEdge[0]  ||
	     ParticlePosition[0][partnum] > RefineRegionRightEdge[0] ||
	     ParticlePosition[1][partnum] < RefineRegionLeftEdge[1]  ||
	     ParticlePosition[1][partnum] > RefineRegionRightEdge[1] ||
	     ParticlePosition[2][partnum] < RefineRegionLeftEdge[2]  ||
	     ParticlePosition[2][partnum] > RefineRegionRightEdge[2])
	    {
	      //fprintf(stdout, "grid::PS: Particle outside RR - ignoring\n"); 
	      continue;
	    }

	  /*
	   *  Compute index of the cell that the parent particle resides in.
	   */
	  xindex = (int)((ParticlePosition[0][partnum] - CellLeftEdge[0][0]) / dx);
	  yindex = (int)((ParticlePosition[1][partnum] - CellLeftEdge[1][0]) / dx); 
	  zindex = (int)((ParticlePosition[2][partnum] - CellLeftEdge[2][0]) / dx); 

	  if (xindex < 0 || xindex >= GridDimension[0] || 
	      yindex < 0 || yindex >= GridDimension[1] || 
	      zindex < 0 || zindex >= GridDimension[2])
	    {
	      fprintf(stdout, "grid::PS: parent particle out of grid (C level): \n");
	      fprintf(stdout, "xind, yind, zind = %ld, %ld, %ld\n", xindex, yindex, zindex); 
	      continue;
	    }

	  /*
	   * CREATE CHILDREN PARTICLES 
	   */
	  /* First reduce the mass of the parent down to 1/(children+parent) = 1/13 */
	  ParticleMass[partnum] = ParticleMass[partnum] / (float)(ChildrenPerParent + 1.0);
	  /*
	   * ===================================
	   * Step I - Positioning
	   * ===================================
  
	   * Now carefully set the position.  Populate 12 children particles 
	   * around the parent at (0,0,0) assuming HCP (hexagonal closed-packed) 
	   * structure in an xyz space.  Keep in mind that Any three neighboring 
	   * children and the parent will form a tetrahedron of edge length rad.  


	   *  Step[I-1]: 6 children on the same x-y plane as the parent 
	   *  (counter-clockwise; again, here we assume the parent is at the origin)
	   */
	  sin60 = sin(M_PI/3.0);
	  NewPos[0][0] = rad;
	  NewPos[1][0] = 0.0;
	  NewPos[2][0] = 0.0;
	  
	  NewPos[0][1] = 0.5*rad;
	  NewPos[1][1] = sin60*rad;
	  NewPos[2][1] = 0.0;

	  NewPos[0][2] = -0.5*rad;
	  NewPos[1][2] = sin60*rad;
	  NewPos[2][2] = 0.0;

	  NewPos[0][3] = -1.0*rad;
	  NewPos[1][3] = 0.0;
	  NewPos[2][3] = 0.0;
	  
	  NewPos[0][4] = -0.5*rad;
	  NewPos[1][4] = -1.0*sin60*rad;
	  NewPos[2][4] = 0.0;

	  NewPos[0][5] = 0.5*rad;
	  NewPos[1][5] = -1.0*sin60*rad;
	  NewPos[2][5] = 0.0;

	  /*
	   * Step[I-2]: 3 children above the parent plane (x-y)
	   * (the height of the tetrahedron of edge r = sqrt(2/3)*r)
	   */
	  NewPos[0][6] = 0.5*rad;
	  NewPos[1][6] = 0.2887*rad;
	  NewPos[2][6] = sqrt(2.0/3.0)*rad;
	    
	  NewPos[0][7] = -0.5*rad;
	  NewPos[1][7] = 0.2887*rad;
	  NewPos[2][7] = sqrt(2.0/3.0)*rad;

	  NewPos[0][8] = 0.0;
	  NewPos[1][8] = -0.5774*rad;
	  NewPos[2][8] = sqrt(2.0/3.0)*rad;

	  /*
	   * Step[I-3]: 3 children below the parent plane (x-y)
	   */
	  NewPos[0][9] = 0.5*rad;
	  NewPos[1][9] = 0.2887*rad;
	  NewPos[2][9] = -1.0*sqrt(2.0/3.0)*rad;
	  
	  NewPos[0][10] = -0.5*rad;
	  NewPos[1][10] = 0.2887*rad;
	  NewPos[2][10] = -1.0*sqrt(2.0/3.0)*rad;

	  NewPos[0][11] = 0.0;
	  NewPos[1][11] = -0.5774*rad;
	  NewPos[2][11] = -1.0*sqrt(2.0/3.0)*rad;

	  /*
	   *  ===================================
	   *       Step II - Rotation
	   *  ===================================
	   *
	   *  Provide the particle with Euler rotation (z->x->z) to an arbitrary 
	   *  orientation.  See Eq.(11.99) of Marion & Thornton (1995)
	   *  -----------------------------------------------------------------   
	   *             lambda_11 =  cosCcosA - cosBsinAsinC
	   *             lambda_21 = -sinCcosA - cosBsinAcosC
	   *             lambda_31 =  sinBsinA
	   *             lambda_12 =  cosCsinA + cosBcosAsinC
	   *             lambda_22 = -sinCsinA + cosBcosAcosC
	   *             lambda_32 = -sinBcosA
	   *             lambda_13 =  sinCsinB
	   *             lambda_23 =  cosCsinB
	   *             lambda_33 =  cosB
           *  -----------------------------------------------------------------
	   */

	  for(m = 0; m < 3; m++)
	  {
	      alpha[m] = (float)(drand48()*M_PI*2.0);
	  }
	  
	  /* 
	   * Calculate the Euler Angles
	   */
	  l11 = cos(alpha[2])*cos(alpha[0]) - cos(alpha[1])*sin(alpha[0])*sin(alpha[2]);
	  l21 = -1.0*sin(alpha[2])*cos(alpha[0]) -1.0*cos(alpha[1])*sin(alpha[0])*cos(alpha[2]);
	  l31 = sin(alpha[1])*sin(alpha[0]);

	  l12 = cos(alpha[2])*sin(alpha[0]) + cos(alpha[1])*cos(alpha[0])*sin(alpha[2]);
	  l22 = -1.0*sin(alpha[2])*sin(alpha[0]) + cos(alpha[1])*cos(alpha[0])*cos(alpha[2]);
	  l32 = -1.0*sin(alpha[1])*cos(alpha[0]);
	  
	  l13 = sin(alpha[2])*sin(alpha[1]);
	  l23 = cos(alpha[2])*sin(alpha[1]);
	  l33 = cos(alpha[1]);

	  
	  /*          
	   * Step[II-3]: Use NewPos and Euler angles to generate ChildPos.
	   * Here we add the initial position of the parent [xyz]pos.		 
	   *
	   * Step III
	   * Copy other particle attributes from parent to child
	   */

	  for(child = total_children, innerchild = 0 ; 
	      child < total_children + CHILDRENPERPARENT; 
	      child++, innerchild++)
	    {
	      this->ParticlePosition[0][child] = ParticlePosition[0][partnum] +
		l11*NewPos[0][innerchild] + l12*NewPos[1][innerchild] + l13*NewPos[2][innerchild];
	      this->ParticlePosition[1][child] = ParticlePosition[1][partnum] +
		l21*NewPos[0][innerchild] + l22*NewPos[1][innerchild] + l23*NewPos[2][innerchild];
	      this->ParticlePosition[2][child] = ParticlePosition[2][partnum] +
		l31*NewPos[0][innerchild] + l32*NewPos[1][innerchild] + l33*NewPos[2][innerchild];

	      
	      if(this->ParticlePosition[0][child] < 0.0  ||
		 this->ParticlePosition[0][child] > 1.0  ||
		 this->ParticlePosition[1][child] < 0.0  ||
		 this->ParticlePosition[1][child] > 1.0  ||
		 this->ParticlePosition[2][child] < 0.0  ||
		 this->ParticlePosition[2][child] > 1.0)
		{
		  fprintf(stderr, "WARNING - Child kicked outside domain\n");
		  fprintf(stderr, "OldPos[%d] = (%f %f %f)\n NewPos[%d] = (%f %f %f)\n\n", 
			 partnum, ParticlePosition[0][partnum],
			 ParticlePosition[1][partnum], ParticlePosition[2][partnum], 
			 child, this->ParticlePosition[0][child],
			 this->ParticlePosition[1][child], this->ParticlePosition[2][child]);
		  fprintf(stderr, "alpha = (%f, %f, %f)\n", alpha[0], alpha[1], alpha[2]);
		  fprintf(stderr, "l11 = %f\n", l11); 
		  fprintf(stderr, "l12 = %f\n", l12); 
		  fprintf(stderr, "l13 = %f\n", l13);
		  for(i = 0; i < 3; i++)
		    {
		      fprintf(stderr, "NewPos[%d][%d] = %f\n", i, child, NewPos[0][innerchild]);
		    }
		  fprintf(stderr, "l11*NewPos[%d][%d] = %f\n", 0, child, l11*NewPos[0][innerchild]);
		  fprintf(stderr, "l12*NewPos[%d][%d] = %f\n", 1, child, l12*NewPos[1][innerchild]);
		  fprintf(stderr, "l13*NewPos[%d][%d] = %f\n", 2, child, l13*NewPos[2][innerchild]);
		  fprintf(stderr, "l21*NewPos[%d][%d] = %f\n", 0, child, l21*NewPos[0][innerchild]);
		  fprintf(stderr, "l22*NewPos[%d][%d] = %f\n", 1, child, l22*NewPos[1][innerchild]);
		  fprintf(stderr, "l23*NewPos[%d][%d] = %f\n", 2, child, l23*NewPos[2][innerchild]);
		 
		  return FAIL;
		  
		}
	 
	      this->ParticleMass[child] = ParticleMass[partnum];
	      for(i = 0; i < 3; i++)
		this->ParticleVelocity[i][child] = ParticleVelocity[i][partnum];
	      this->ParticleType[child] = ParticleType[partnum];
	      for(i = 0; i < NumberAttributes; i++)
		this->ParticleAttribute[i][child] = ParticleAttribute[i][numpart];
	    }

	  /* Loop forward by CHILDRENPERPARENT each time. */
	  total_children = total_children + CHILDRENPERPARENT;

	  if(total_children > MaximumNumberOfNewParticles)
	    {
	       fprintf(stdout, "Total number of Children (%ld) exceeded the maximum (%ld)\n", 
		       total_children, MaximumNumberOfNewParticles);
	       return FAIL;
	    }

	}
    }
#ifdef DEBUG_PS
  fprintf(stdout,"%d new child particles created in this grid from %d parents", 
	  total_children, NumberOfParticles);
#endif
 
  *NumberOfNewParticles = total_children;
  return SUCCESS;
}
