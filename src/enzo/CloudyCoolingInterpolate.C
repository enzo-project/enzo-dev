/***********************************************************************
/
/  CLOUDY COOLING INTERPOLATION ROUTINES
/
/  written by: Britton Smith
/  date:       November, 2005
/  modified1:  May, 2009
/              Added 4d and 5d interpolation.
/
/  PURPOSE:  Perform interpolation over Cloudy heatin, cooling, and 
/            mean molecular weight data.
/
/  Rank = 1: interpolate over temperature.
/  Rank = 2: interpolate over density and temperature.
/  Rank = 3: interpolate over density, metallicity, and temperature.
/  Rank = 4: interpolate over density, metallicity, electron fraction, 
/            and temperature.
/  Rank = 5: interpolate over density, metallicity, electron fraction, 
/            redshift (for Haardt Madau background), and temperature.
/
/  RETURNS:
/    Interpolated value.
/
************************************************************************/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "hdf5.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "CosmologyParameters.h"

// Interpolate over temperature.

float coolingGridInterpolate1D(float temperature,float *dataField)
{
  int temperatureIndex;
  int midpoint,highpoint;
  float slope,temperatureValue;

  // find temperature index

  if (temperature <= CloudyCoolingData.CloudyCoolingGridParameters[0][0]) {
    temperatureIndex = 0;
  }
  else if (temperature >= CloudyCoolingData.CloudyCoolingGridParameters[0][CloudyCoolingData.CloudyCoolingGridDimension[0] - 1]) {
    temperatureIndex = CloudyCoolingData.CloudyCoolingGridDimension[0] - 2;
  }
  else {
    temperatureIndex = 0;
    highpoint = CloudyCoolingData.CloudyCoolingGridDimension[0] - 1;
    while (highpoint - temperatureIndex > 1) {
      midpoint = (highpoint + temperatureIndex) >> 1;
      if (temperature >= CloudyCoolingData.CloudyCoolingGridParameters[0][midpoint]) temperatureIndex = midpoint;
      else highpoint = midpoint;
    }
  }

  // interpolate over temperature
    
  slope = (dataField[temperatureIndex+1] - dataField[temperatureIndex]) /
    (CloudyCoolingData.CloudyCoolingGridParameters[0][temperatureIndex+1] - 
     CloudyCoolingData.CloudyCoolingGridParameters[0][temperatureIndex]);

  temperatureValue = (temperature - CloudyCoolingData.CloudyCoolingGridParameters[0][temperatureIndex]) * slope +
    dataField[temperatureIndex];

    return temperatureValue;
}

// Interpolate over density and temperature.

float coolingGridInterpolate2D(float parameter1,float temperature,float *dataField)
{
  int parameter1Index, temperatureIndex;
  int midpoint,highpoint;
  int interpolateIndex,q;
  float slope,temperatureValue[2],p1Value;

  // find parameter1 index

  if (parameter1 <= CloudyCoolingData.CloudyCoolingGridParameters[0][0]) {
    parameter1Index = 0;
  }
  else if (parameter1 >= CloudyCoolingData.CloudyCoolingGridParameters[0][CloudyCoolingData.CloudyCoolingGridDimension[0] - 1]) {
    parameter1Index = CloudyCoolingData.CloudyCoolingGridDimension[0] - 2;
  }
  else {
    parameter1Index = 0;
    highpoint = CloudyCoolingData.CloudyCoolingGridDimension[0] - 1;
    while (highpoint - parameter1Index > 1) {
      midpoint = (highpoint + parameter1Index) >> 1;
      if (parameter1 >= CloudyCoolingData.CloudyCoolingGridParameters[0][midpoint]) parameter1Index = midpoint;
      else highpoint = midpoint;
    }
  }

  // find temperature index

  if (temperature <= CloudyCoolingData.CloudyCoolingGridParameters[1][0]) {
    temperatureIndex = 0;
  }
  else if (temperature >= CloudyCoolingData.CloudyCoolingGridParameters[1][CloudyCoolingData.CloudyCoolingGridDimension[1] - 1]) {
    temperatureIndex = CloudyCoolingData.CloudyCoolingGridDimension[1] - 2;
  }
  else {
    temperatureIndex = 0;
    highpoint = CloudyCoolingData.CloudyCoolingGridDimension[1] - 1;
    while (highpoint - temperatureIndex > 1) {
      midpoint = (highpoint + temperatureIndex) >> 1;
      if (temperature >= CloudyCoolingData.CloudyCoolingGridParameters[1][midpoint]) temperatureIndex = midpoint;
      else highpoint = midpoint;
    }
  }

  // interpolate over parameter 1 (density)

  for (q = 0;q < 2;q++) {

    // interpolate over temperature

    interpolateIndex = (q+parameter1Index) * CloudyCoolingData.CloudyCoolingGridDimension[1] + temperatureIndex;
    
    slope = (dataField[interpolateIndex+1] - dataField[interpolateIndex]) /
      (CloudyCoolingData.CloudyCoolingGridParameters[1][temperatureIndex+1] - 
       CloudyCoolingData.CloudyCoolingGridParameters[1][temperatureIndex]);
    
    temperatureValue[q] = (temperature - CloudyCoolingData.CloudyCoolingGridParameters[1][temperatureIndex]) * slope +
      dataField[interpolateIndex];
  }
  slope = (temperatureValue[1]-temperatureValue[0]) /
    (CloudyCoolingData.CloudyCoolingGridParameters[0][parameter1Index+1]-
     CloudyCoolingData.CloudyCoolingGridParameters[0][parameter1Index]);

  p1Value = (parameter1-CloudyCoolingData.CloudyCoolingGridParameters[0][parameter1Index]) * slope + temperatureValue[0];
  
  return p1Value;
}

// Interpolate over density, metallicity, and temperature.

float coolingGridInterpolate3D(float parameter1,float parameter2,float temperature,float *dataField)
{
  int parameter1Index, parameter2Index, temperatureIndex;
  int midpoint,highpoint;
  int interpolateIndex,q,w;
  float slope,temperatureValue[2],p2Value[2],p1Value;

  // find parameter1 index

  if (parameter1 <= CloudyCoolingData.CloudyCoolingGridParameters[0][0]) {
    parameter1Index = 0;
  }
  else if (parameter1 >= CloudyCoolingData.CloudyCoolingGridParameters[0][CloudyCoolingData.CloudyCoolingGridDimension[0] - 1]) {
    parameter1Index = CloudyCoolingData.CloudyCoolingGridDimension[0] - 2;
  }
  else {
    parameter1Index = 0;
    highpoint = CloudyCoolingData.CloudyCoolingGridDimension[0] - 1;
    while (highpoint - parameter1Index > 1) {
      midpoint = (highpoint + parameter1Index) >> 1;
      if (parameter1 >= CloudyCoolingData.CloudyCoolingGridParameters[0][midpoint]) parameter1Index = midpoint;
      else highpoint = midpoint;
    }
  }

  // find parameter2 index

  if (parameter2 <= CloudyCoolingData.CloudyCoolingGridParameters[1][0]) {
    parameter2Index = 0;
  }
  else if (parameter2 >= CloudyCoolingData.CloudyCoolingGridParameters[1][CloudyCoolingData.CloudyCoolingGridDimension[1] - 1]) {
    parameter2Index = CloudyCoolingData.CloudyCoolingGridDimension[1] - 2;
  }
  else {
    parameter2Index = 0;
    highpoint = CloudyCoolingData.CloudyCoolingGridDimension[1] - 1;
    while (highpoint - parameter2Index > 1) {
      midpoint = (highpoint + parameter2Index) >> 1;
      if (parameter2 >= CloudyCoolingData.CloudyCoolingGridParameters[1][midpoint]) parameter2Index = midpoint;
      else highpoint = midpoint;
    }
  }

  // find temperature index

  if (temperature <= CloudyCoolingData.CloudyCoolingGridParameters[2][0]) {
    temperatureIndex = 0;
  }
  else if (temperature >= CloudyCoolingData.CloudyCoolingGridParameters[2][CloudyCoolingData.CloudyCoolingGridDimension[2] - 1]) {
    temperatureIndex = CloudyCoolingData.CloudyCoolingGridDimension[2] - 2;
  }
  else {
    temperatureIndex = 0;
    highpoint = CloudyCoolingData.CloudyCoolingGridDimension[2] - 1;
    while (highpoint - temperatureIndex > 1) {
      midpoint = (highpoint + temperatureIndex) >> 1;
      if (temperature >= CloudyCoolingData.CloudyCoolingGridParameters[2][midpoint]) temperatureIndex = midpoint;
      else highpoint = midpoint;
    }
  }

  // interpolate over parameter 1

  for (q = 0;q < 2;q++) {

    // interpolate over parameter 2

    for (w = 0;w < 2;w++) {

      // interpolate over temperature

      interpolateIndex = ((q+parameter1Index)*CloudyCoolingData.CloudyCoolingGridDimension[1] + (w+parameter2Index))
	* CloudyCoolingData.CloudyCoolingGridDimension[2] + temperatureIndex;

      slope = (dataField[interpolateIndex+1] - dataField[interpolateIndex]) /
	(CloudyCoolingData.CloudyCoolingGridParameters[2][temperatureIndex+1] - 
	 CloudyCoolingData.CloudyCoolingGridParameters[2][temperatureIndex]);

      temperatureValue[w] = (temperature - CloudyCoolingData.CloudyCoolingGridParameters[2][temperatureIndex]) * slope +
	dataField[interpolateIndex];
    }
    slope = (temperatureValue[1]-temperatureValue[0]) /
      (CloudyCoolingData.CloudyCoolingGridParameters[1][parameter2Index+1]-CloudyCoolingData.CloudyCoolingGridParameters[1][parameter2Index]);
    p2Value[q] = (parameter2-CloudyCoolingData.CloudyCoolingGridParameters[1][parameter2Index]) * slope + temperatureValue[0];
  }
  slope = (p2Value[1]-p2Value[0]) /
    (CloudyCoolingData.CloudyCoolingGridParameters[0][parameter1Index+1]-CloudyCoolingData.CloudyCoolingGridParameters[0][parameter1Index]);
  p1Value = (parameter1-CloudyCoolingData.CloudyCoolingGridParameters[0][parameter1Index]) * slope + p2Value[0];

  return p1Value;
}

// Interpolate over density, metallicity, electron fraction, and temperature.

float coolingGridInterpolate4D(float parameter1,float parameter2,float parameter3,float temperature,float *dataField)
{
  int parameter1Index, parameter2Index, parameter3Index, temperatureIndex;
  int midpoint,highpoint;
  int interpolateIndex,q,w,e;
  float slope,temperatureValue[2],p3Value[2],p2Value[2],p1Value;

  // find parameter1 index

  if (parameter1 <= CloudyCoolingData.CloudyCoolingGridParameters[0][0]) {
    parameter1Index = 0;
  }
  else if (parameter1 >= CloudyCoolingData.CloudyCoolingGridParameters[0][CloudyCoolingData.CloudyCoolingGridDimension[0] - 1]) {
    parameter1Index = CloudyCoolingData.CloudyCoolingGridDimension[0] - 2;
  }
  else {
    parameter1Index = 0;
    highpoint = CloudyCoolingData.CloudyCoolingGridDimension[0] - 1;
    while (highpoint - parameter1Index > 1) {
      midpoint = (highpoint + parameter1Index) >> 1;
      if (parameter1 >= CloudyCoolingData.CloudyCoolingGridParameters[0][midpoint]) parameter1Index = midpoint;
      else highpoint = midpoint;
    }
  }

  // find parameter2 index

  if (parameter2 <= CloudyCoolingData.CloudyCoolingGridParameters[1][0]) {
    parameter2Index = 0;
  }
  else if (parameter2 >= CloudyCoolingData.CloudyCoolingGridParameters[1][CloudyCoolingData.CloudyCoolingGridDimension[1] - 1]) {
    parameter2Index = CloudyCoolingData.CloudyCoolingGridDimension[1] - 2;
  }
  else {
    parameter2Index = 0;
    highpoint = CloudyCoolingData.CloudyCoolingGridDimension[1] - 1;
    while (highpoint - parameter2Index > 1) {
      midpoint = (highpoint + parameter2Index) >> 1;
      if (parameter2 >= CloudyCoolingData.CloudyCoolingGridParameters[1][midpoint]) parameter2Index = midpoint;
      else highpoint = midpoint;
    }
  }

  // find parameter3 index

  if (parameter3 <= CloudyCoolingData.CloudyCoolingGridParameters[2][0]) {
    parameter3Index = 0;
  }
  else if (parameter3 >= CloudyCoolingData.CloudyCoolingGridParameters[2][CloudyCoolingData.CloudyCoolingGridDimension[2] - 1]) {
    parameter3Index = CloudyCoolingData.CloudyCoolingGridDimension[2] - 2;
  }
  else {
    parameter3Index = 0;
    highpoint = CloudyCoolingData.CloudyCoolingGridDimension[2] - 1;
    while (highpoint - parameter3Index > 1) {
      midpoint = (highpoint + parameter3Index) >> 1;
      if (parameter3 >= CloudyCoolingData.CloudyCoolingGridParameters[2][midpoint]) parameter3Index = midpoint;
      else highpoint = midpoint;
    }
  }

  // find temperature index

  if (temperature <= CloudyCoolingData.CloudyCoolingGridParameters[3][0]) {
    temperatureIndex = 0;
  }
  else if (temperature >= CloudyCoolingData.CloudyCoolingGridParameters[3][CloudyCoolingData.CloudyCoolingGridDimension[3] - 1]) {
    temperatureIndex = CloudyCoolingData.CloudyCoolingGridDimension[3] - 2;
  }
  else {
    temperatureIndex = 0;
    highpoint = CloudyCoolingData.CloudyCoolingGridDimension[3] - 1;
    while (highpoint - temperatureIndex > 1) {
      midpoint = (highpoint + temperatureIndex) >> 1;
      if (temperature >= CloudyCoolingData.CloudyCoolingGridParameters[3][midpoint]) temperatureIndex = midpoint;
      else highpoint = midpoint;
    }
  }

  // interpolate over parameter 1

  for (q = 0;q < 2;q++) {

    // interpolate over parameter 2

    for (w = 0;w < 2;w++) {

      // interpolate over parameter 3

      for (e = 0;e < 2;e++) {

      // interpolate over temperature

	interpolateIndex = (((q+parameter1Index) * CloudyCoolingData.CloudyCoolingGridDimension[1] +
			     (w+parameter2Index)) * CloudyCoolingData.CloudyCoolingGridDimension[2] +
			    (e+parameter3Index)) * CloudyCoolingData.CloudyCoolingGridDimension[3] + temperatureIndex;

	slope = (dataField[interpolateIndex+1] - dataField[interpolateIndex]) /
	  (CloudyCoolingData.CloudyCoolingGridParameters[3][temperatureIndex+1] - 
	   CloudyCoolingData.CloudyCoolingGridParameters[3][temperatureIndex]);

	temperatureValue[e] = (temperature - CloudyCoolingData.CloudyCoolingGridParameters[3][temperatureIndex]) * slope +
	  dataField[interpolateIndex];

      }

      slope = (temperatureValue[1]-temperatureValue[0]) /
	(CloudyCoolingData.CloudyCoolingGridParameters[2][parameter3Index+1] -
	 CloudyCoolingData.CloudyCoolingGridParameters[2][parameter3Index]);

      p3Value[w] = (parameter3-CloudyCoolingData.CloudyCoolingGridParameters[2][parameter3Index]) * slope + temperatureValue[0];

    }

    slope = (p3Value[1]-p3Value[0]) /
      (CloudyCoolingData.CloudyCoolingGridParameters[1][parameter2Index+1] -
       CloudyCoolingData.CloudyCoolingGridParameters[1][parameter2Index]);

    p2Value[q] = (parameter2-CloudyCoolingData.CloudyCoolingGridParameters[1][parameter2Index]) * slope + p3Value[0];

  }

  slope = (p2Value[1]-p2Value[0]) /
    (CloudyCoolingData.CloudyCoolingGridParameters[0][parameter1Index+1] -
     CloudyCoolingData.CloudyCoolingGridParameters[0][parameter1Index]);

  p1Value = (parameter1-CloudyCoolingData.CloudyCoolingGridParameters[0][parameter1Index]) * slope + p2Value[0];

  return p1Value;
}

// Interpolate over density, metallicity, electron fraction, radiation background, and temperature.

float coolingGridInterpolate5D(float parameter1,float parameter2,float parameter3,float parameter4,float temperature,float *dataField)
{
  int parameter1Index, parameter2Index, parameter3Index, parameter4Index, temperatureIndex;
  int midpoint,highpoint;
  int interpolateIndex,q,w,e,r;
  float slope,temperatureValue[2],p4Value[2],p3Value[2],p2Value[2],p1Value;

  // find parameter1 index

  if (parameter1 <= CloudyCoolingData.CloudyCoolingGridParameters[0][0]) {
    parameter1Index = 0;
  }
  else if (parameter1 >= CloudyCoolingData.CloudyCoolingGridParameters[0][CloudyCoolingData.CloudyCoolingGridDimension[0] - 1]) {
    parameter1Index = CloudyCoolingData.CloudyCoolingGridDimension[0] - 2;
  }
  else {
    parameter1Index = 0;
    highpoint = CloudyCoolingData.CloudyCoolingGridDimension[0] - 1;
    while (highpoint - parameter1Index > 1) {
      midpoint = (highpoint + parameter1Index) >> 1;
      if (parameter1 >= CloudyCoolingData.CloudyCoolingGridParameters[0][midpoint]) parameter1Index = midpoint;
      else highpoint = midpoint;
    }
  }

  // find parameter2 index

  if (parameter2 <= CloudyCoolingData.CloudyCoolingGridParameters[1][0]) {
    parameter2Index = 0;
  }
  else if (parameter2 >= CloudyCoolingData.CloudyCoolingGridParameters[1][CloudyCoolingData.CloudyCoolingGridDimension[1] - 1]) {
    parameter2Index = CloudyCoolingData.CloudyCoolingGridDimension[1] - 2;
  }
  else {
    parameter2Index = 0;
    highpoint = CloudyCoolingData.CloudyCoolingGridDimension[1] - 1;
    while (highpoint - parameter2Index > 1) {
      midpoint = (highpoint + parameter2Index) >> 1;
      if (parameter2 >= CloudyCoolingData.CloudyCoolingGridParameters[1][midpoint]) parameter2Index = midpoint;
      else highpoint = midpoint;
    }
  }

  // find parameter3 index

  if (parameter3 <= CloudyCoolingData.CloudyCoolingGridParameters[2][0]) {
    parameter3Index = 0;
  }
  else if (parameter3 >= CloudyCoolingData.CloudyCoolingGridParameters[2][CloudyCoolingData.CloudyCoolingGridDimension[2] - 1]) {
    parameter3Index = CloudyCoolingData.CloudyCoolingGridDimension[2] - 2;
  }
  else {
    parameter3Index = 0;
    highpoint = CloudyCoolingData.CloudyCoolingGridDimension[2] - 1;
    while (highpoint - parameter3Index > 1) {
      midpoint = (highpoint + parameter3Index) >> 1;
      if (parameter3 >= CloudyCoolingData.CloudyCoolingGridParameters[2][midpoint]) parameter3Index = midpoint;
      else highpoint = midpoint;
    }
  }

  // find parameter4 index

  if (parameter4 <= CloudyCoolingData.CloudyCoolingGridParameters[3][0]) {
    parameter4Index = 0;
  }
  else if (parameter4 >= CloudyCoolingData.CloudyCoolingGridParameters[3][CloudyCoolingData.CloudyCoolingGridDimension[3] - 1]) {
    parameter4Index = CloudyCoolingData.CloudyCoolingGridDimension[3] - 2;
  }
  else {
    parameter4Index = 0;
    highpoint = CloudyCoolingData.CloudyCoolingGridDimension[3] - 1;
    while (highpoint - parameter4Index > 1) {
      midpoint = (highpoint + parameter4Index) >> 1;
      if (parameter4 >= CloudyCoolingData.CloudyCoolingGridParameters[3][midpoint]) parameter4Index = midpoint;
      else highpoint = midpoint;
    }
  }

  // find temperature index

  if (temperature <= CloudyCoolingData.CloudyCoolingGridParameters[4][0]) {
    temperatureIndex = 0;
  }
  else if (temperature >= CloudyCoolingData.CloudyCoolingGridParameters[4][CloudyCoolingData.CloudyCoolingGridDimension[4] - 1]) {
    temperatureIndex = CloudyCoolingData.CloudyCoolingGridDimension[4] - 2;
  }
  else {
    temperatureIndex = 0;
    highpoint = CloudyCoolingData.CloudyCoolingGridDimension[4] - 1;
    while (highpoint - temperatureIndex > 1) {
      midpoint = (highpoint + temperatureIndex) >> 1;
      if (temperature >= CloudyCoolingData.CloudyCoolingGridParameters[4][midpoint]) temperatureIndex = midpoint;
      else highpoint = midpoint;
    }
  }

  // interpolate over parameter 1

  for (q = 0;q < 2;q++) {

    // interpolate over parameter 2

    for (w = 0;w < 2;w++) {

      // interpolate over parameter 3

      for (e = 0;e < 2;e++) {

	// interpolate over parameter 4

	for (r = 0;r < 2;r++) {

	  // interpolate over temperature

	  interpolateIndex = ((((q+parameter1Index) * CloudyCoolingData.CloudyCoolingGridDimension[1] +
				(w+parameter2Index)) * CloudyCoolingData.CloudyCoolingGridDimension[2] +
			       (e+parameter3Index)) * CloudyCoolingData.CloudyCoolingGridDimension[3] + 
			      (r+parameter4Index)) * CloudyCoolingData.CloudyCoolingGridDimension[4] + temperatureIndex;

	  slope = (dataField[interpolateIndex+1] - dataField[interpolateIndex]) /
	    (CloudyCoolingData.CloudyCoolingGridParameters[4][temperatureIndex+1] - 
	     CloudyCoolingData.CloudyCoolingGridParameters[4][temperatureIndex]);

	  temperatureValue[r] = (temperature - CloudyCoolingData.CloudyCoolingGridParameters[4][temperatureIndex]) * slope +
	    dataField[interpolateIndex];

	}

	slope = (temperatureValue[1]-temperatureValue[0]) /
	  (CloudyCoolingData.CloudyCoolingGridParameters[3][parameter4Index+1] -
	   CloudyCoolingData.CloudyCoolingGridParameters[3][parameter4Index]);

	p4Value[e] = (parameter4-CloudyCoolingData.CloudyCoolingGridParameters[3][parameter4Index]) * slope + temperatureValue[0];

      }

      slope = (p4Value[1]-p4Value[0]) /
	(CloudyCoolingData.CloudyCoolingGridParameters[2][parameter3Index+1] -
	 CloudyCoolingData.CloudyCoolingGridParameters[2][parameter3Index]);

      p3Value[w] = (parameter3-CloudyCoolingData.CloudyCoolingGridParameters[2][parameter3Index]) * slope + p4Value[0];

    }

    slope = (p3Value[1]-p3Value[0]) /
      (CloudyCoolingData.CloudyCoolingGridParameters[1][parameter2Index+1] -
       CloudyCoolingData.CloudyCoolingGridParameters[1][parameter2Index]);

    p2Value[q] = (parameter2-CloudyCoolingData.CloudyCoolingGridParameters[1][parameter2Index]) * slope + p3Value[0];

  }

  slope = (p2Value[1]-p2Value[0]) /
    (CloudyCoolingData.CloudyCoolingGridParameters[0][parameter1Index+1] -
     CloudyCoolingData.CloudyCoolingGridParameters[0][parameter1Index]);

  p1Value = (parameter1-CloudyCoolingData.CloudyCoolingGridParameters[0][parameter1Index]) * slope + p2Value[0];

  return p1Value;
}
