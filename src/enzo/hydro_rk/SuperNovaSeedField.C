//
//  SuperNovaSeedField.cpp
//
//  Defines the SuperNova class.
//
//  getSourceTerms method calculates the time derivative of the magnetic field
//  to be injected.
//
//  Created by Iryna Butsky on 1/5/15.
//
//
#define USE
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include "macros_and_parameters.h"
#include "typedefs.h"

#include "SuperNova.h"


SuperNova::SuperNova(){
	// Set zero values (which result in no supernova feedback) in case SuperNova 

  // object is initialized but not given proper values
	
  zHat[0] = FLOAT();
  zHat[1] = FLOAT();
  zHat[2] = FLOAT();
	
  location[0] = FLOAT();
  location[1] = FLOAT();
  location[2] = FLOAT();
	
  characteristicLength = FLOAT();
  timeStarted = FLOAT();
  characteristicTime = FLOAT();
	
  totalEnergy = float();
  sigma = float();
	
}


void SuperNova::setValues(FLOAT phi_x, FLOAT phi_y, FLOAT phi_z,  FLOAT x, FLOAT y, FLOAT z, \
			 FLOAT radius, FLOAT time_started, FLOAT duration, float energy, float sigma_sn){
	zHat[0] = phi_x;
	zHat[1] = phi_y;
	zHat[2] = phi_z;
	
	location[0] = x;
	location[1] = y;
	location[2] = z;
	
	characteristicLength = radius/3.0;
	timeStarted = time_started;
	characteristicTime = duration/5.0;
	
	totalEnergy = energy;
	sigma = sigma_sn;
	
}

FLOAT* SuperNova::getPosition(){
  return location;
}

static FLOAT Magnitude(FLOAT x[3]){
	
	return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
}

static void Rotate_Vector(FLOAT v_z[3], FLOAT v_sn[3], FLOAT v_to_rotate[3]) {
	/* Rotates vector v_to_rotate (actually changes that vector, so be careful!)
	 FROM the direction of v_z TO the direction of v_sn */
	
	float v_z_norm[3], v_sn_norm[3], u[3], cross[3];
	FLOAT copy[3];
	float R[3][3];
	
	// normalizing v_z and v_sn in case they aren't already
	// copy[3] will be a helper vector when reassigning values in v_to_rotate
	// temp stores the normalized v_sn, which needs to undergo more
	for (int i = 0; i < 3; ++i) {
		copy[i] = v_to_rotate[i];
		v_z_norm[i] = v_z[i]/Magnitude(v_z);
		v_sn_norm[i] = v_sn[i]/Magnitude(v_sn);
	}
	
	// calculating the cross product and normalizing it
	cross[0] = v_z_norm[1]*v_sn_norm[2] - v_z_norm[2]*v_sn_norm[1];
	cross[1] = -v_z_norm[0]*v_sn_norm[2] + v_z_norm[2]*v_sn_norm[0];
	cross[2] = v_z_norm[0]*v_sn_norm[1] - v_z_norm[1]*v_sn_norm[0];
	
	float mag_cross = sqrt(cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2]);
	if ( mag_cross == 0) { mag_cross += 1;}
	
	u[0] = cross[0] / mag_cross;
	u[1] = cross[1] / mag_cross;
	u[2] = cross[2] / mag_cross;
	
	// Finding angle with dot product. Mag_z and Mag_u should be = 1
	float cos_theta = (v_z_norm[0]* v_sn_norm[0] + v_z_norm[1]*v_sn_norm[1] \
				 + v_z_norm[2]*v_sn_norm[2]);
	float sin_theta = sin(acos(cos_theta));
	
	// Rotation matrix s.t. u = R*z
	// Handy algorithm on wikipedia page on rotation matrices
	R[0][0] = cos_theta + POW(u[0], 2)*(1-cos_theta);
	R[0][1] = u[0]*u[1]*(1-cos_theta) - u[2]*sin_theta;
	R[0][2] = u[0]*u[2]*(1-cos_theta) + u[1]*sin_theta;
	
	R[1][0] = u[1]*u[0]*(1-cos_theta) + u[2]*sin_theta;
	R[1][1] = cos_theta + POW(u[1], 2)*(1-cos_theta);
	R[1][2] = u[1]*u[2]*(1-cos_theta) - u[0]*sin_theta;
	
	R[2][0] = u[2]*u[0]*(1-cos_theta) - u[1]*sin_theta;
	R[2][1] = u[2]*u[1]*(1-cos_theta) + u[0]*sin_theta;
	R[2][2] = cos_theta + POW(u[2], 2)*(1-cos_theta);
	
	// Multiplying R*v_to_rotate (essentially rotating v_to_rotate)
	v_to_rotate[0] = R[0][0]*copy[0] + R[0][1]*copy[1] + R[0][2]*copy[2];
	v_to_rotate[1] = R[1][0]*copy[0] + R[1][1]*copy[1] + R[1][2]*copy[2];
	v_to_rotate[2] = R[2][0]*copy[0] + R[2][1]*copy[1] + R[2][2]*copy[2];
	
	
}

snsf_source_terms SuperNova::getSourceTerms(double dx, double dy, double dz, double enzoTime){
	// ------------------------------------------------------------------- //
	// Inputs: (1) distances to current cell, dx, dy, dz                             //
	//         (2) current time of the simulation                          //
	//                                                                     //
	// Output: (1) A struct of source terms corresponding to the given	   //
	//             supernova							               //
	//          (i)   x,y,z components of time derivative of B-field   //
	//          (ii)  time derivative of magnetic energy density	   //
	//    (iii) time derivative of total energy density		   //
	// ------------------------------------------------------------------- //
	
	FLOAT phi[3], zhat_cell[3];
	FLOAT rotated_coords[3];
	
	snsf_source_terms S;

	// the z_hat direction in the cell.
	// chosen to be the z-axis
	zhat_cell[0] = 0;
	zhat_cell[1] = 0;
	zhat_cell[2] = 1;
	
	
	// distance to supernova and time since supernova
	FLOAT r_s = sqrt(dx*dx + dy*dy + dz*dz);
	FLOAT t_s = enzoTime-timeStarted;

	
	if((t_s < 0) || (t_s > 5*characteristicTime)){
	  S.dbx = 0;
	  S.dby = 0;
	  S.dbz = 0;
	  S.dUb = 0;
	  S.dUtot = 0;
	  return S;
	}
	// these coordinates will be rotated in the next step
	// in order to calculate B s.t. z_hat is rotated to be z_hat_sn
	rotated_coords[0] = dx;
	rotated_coords[1] = dy;
	rotated_coords[2] = dz;
	
	
	// Rotated x,y,z coordinates to frame where z axis = direction of supernova
	Rotate_Vector(zhat_cell, zHat, rotated_coords);
	
	
	// how our source terms will scale with radius and time
	// r_cyl is the cylindrical radius in the rotated frame
	FLOAT r_cyl = sqrt(rotated_coords[0]*rotated_coords[0] + \
				 rotated_coords[1]*rotated_coords[1]);
	
	FLOAT r_scale = (r_cyl/characteristicLength)*\
				exp(- r_s*r_s /(characteristicLength*characteristicLength));
	FLOAT t_exp = exp(-t_s/characteristicTime);
	
	FLOAT db_t_exp = t_exp / characteristicTime;
	FLOAT b_t_exp = 1 - t_exp;
	
	// Compute phi_hat in rotated frame
        int z_direction;
	if (zHat[2] < 0) z_direction = -1;
	else z_direction = 1;

	if (r_cyl != 0){
	  phi[0] = z_direction*(-rotated_coords[1]/r_cyl);
	  phi[1] = z_direction*(rotated_coords[0]/r_cyl);
	  phi[2] = 0;

	}
        else { phi[0] = 1; 
	  phi[1] = 1;
	  phi[2] = 0;}
	// Rotate phi back to the orientation of the cell
	Rotate_Vector(zHat, zhat_cell, phi);
	
	FLOAT norm_factor = 4*totalEnergy / (POW(characteristicLength, 3) * M_PI*M_PI);
	FLOAT dB_scale = sqrt(sigma*norm_factor)*sqrt(r_scale)*db_t_exp;
	FLOAT B_scale = sqrt(sigma*norm_factor)*sqrt(r_scale)*b_t_exp;

	
	// the x,y,z components of dB in the cell's frame of reference
	S.dbx = dB_scale * phi[0];
	S.dby = dB_scale * phi[1];
	S.dbz = dB_scale * phi[2];
	
	S.bx = B_scale * phi[0];
	S.by = B_scale * phi[1];
	S.bz = B_scale * phi[2];
       
	S.dUtot = norm_factor*r_scale*db_t_exp*b_t_exp;
	S.dUb = (S.bx*S.dbx + S.by*S.dby + S.bz*S.dbz);
	

	return S;
	
}

