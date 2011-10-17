/***********************************************************************
/
/  GLOBAL DATA DECLARATIONS FOR THE 2D MHD TESTS 
/
/   Needed for the Wengen Coliding flow problem to specify inflow boundaries
/
/  written by: Tom Abel
/  date:       Oct 2010
/  modified1:
/
/  PURPOSE:
/    This is global data that pertains only to the Wengen Coliding flow test problem.
/
************************************************************************/

#ifdef DEFINE_STORAGE
# define MHD2DEXTERN
#else /* DEFINE_STORAGE */
# define MHD2DEXTERN extern
#endif /* DEFINE_STORAGE */

/* Smoothing width  */

MHD2DEXTERN float RampWidth;

MHD2DEXTERN float LowerDensity;
MHD2DEXTERN float UpperDensity;
MHD2DEXTERN float LowerVelocityX;
MHD2DEXTERN float UpperVelocityX;
MHD2DEXTERN float LowerVelocityY;
MHD2DEXTERN float UpperVelocityY;
MHD2DEXTERN float LowerPressure; 
MHD2DEXTERN float UpperPressure;
MHD2DEXTERN float LowerBx; 
MHD2DEXTERN float UpperBx;
MHD2DEXTERN float LowerBy;
MHD2DEXTERN float UpperBy;
MHD2DEXTERN int UseColour;
