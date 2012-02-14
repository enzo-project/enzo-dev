/***********************************************************************
/
/  INTERPRET COMMAND-LINE OPTIONS
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
 
/* function prototypes */
 
void PrintUsage(char *myname);
void my_exit(int status);
void auto_show_compile_options(void);
void WriteConfigure(FILE *fp);

Eint32 hide_isdigit(Eint32 c);

int InterpretCommandLine(int argc, char *argv[], char *myname,
			 int &restart, int &debug, int &extract,
			 int &InformationOutput,
			 int &OutputAsParticleData,
			 int &project, int &ProjectionDimension,
			 int &ProjectionSmooth, int &velanyl,
			 char *ParameterFile[],
			 int RegionStart[], int RegionEnd[],
			 FLOAT RegionStartCoordinate[],
			 FLOAT RegionEndCoordinate[],
			 int &RegionLevel, int &HaloFinderOnly,
			 int &WritePotentialOnly,
			 int &SmoothedDarkMatterOnly,
			 int &WriteCoolingTimeOnly,
			 int MyProcessorNumber)
{
 
  int dim;
 
  /* Initialize Region start/end. */
 
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    RegionStart[dim] = INT_UNDEFINED;
    RegionEnd[dim]   = INT_UNDEFINED;
    RegionStartCoordinate[dim] = FLOAT_UNDEFINED;
    RegionEndCoordinate[dim]   = FLOAT_UNDEFINED;
  }
  RegionLevel = INT_UNDEFINED;
 
  /* Interpret command-line arguments */
 
  char c;
 
  while (--argc > 0 && (*++argv)[0] == '-')
    while ((c = *++argv[0]))
      switch (c) {
 
	/* get beginning of region selection (float coordinate). */
 
      case 'b':
	dim = 0;
	while (dim < MAX_DIMENSION && argc > 1 && hide_isdigit(*argv[1])) {
	  argc--;
	  if (sscanf((*++argv), "%"PSYM, &RegionStartCoordinate[dim++]) != 1) {
	    if (MyProcessorNumber == ROOT_PROCESSOR)
	    ENZO_VFAIL("%s: error reading Begin coordinates\n", myname)
	  }
	  while (*(argv[0]+1))
	    ++argv[0];
	}
	break;
 
	/* Add cooling time to data */

      case 'C':
	WriteCoolingTimeOnly = TRUE;
	break;

	/* debug */
 
      case 'd':
	if (MyProcessorNumber == ROOT_PROCESSOR)
	  debug = TRUE;
	break;
 
	/* get end of region selection (integer index). */
 
      case 'e':
	dim = 0;
	while (dim < MAX_DIMENSION && argc > 1 && hide_isdigit(*argv[1])) {
	  argc--;
	  if (sscanf((*++argv), "%"ISYM, &RegionEnd[dim++]) != 1) {
	    if (MyProcessorNumber == ROOT_PROCESSOR)
	    ENZO_VFAIL("%s: error reading End indexes.\n", myname)
	  }
	  while (*(argv[0]+1))
	    ++argv[0];
	}
	break;
 
	/* get finish of region selection (float coordinate). */
 
      case 'f':
	dim = 0;
	while (dim < MAX_DIMENSION && argc > 1 && hide_isdigit(*argv[1])) {
	  argc--;
	  if (sscanf((*++argv), "%"PSYM, &RegionEndCoordinate[dim++]) != 1) {
	    if (MyProcessorNumber == ROOT_PROCESSOR)
	    ENZO_VFAIL("%s: error reading Finish coordinates\n",myname)
	  }
	  while (*(argv[0]+1))
	    ++argv[0];
	}
	break;

	/* FOF halo finder only */

      case 'F':
	HaloFinderOnly = TRUE;
	break;

      case 'g':
	WritePotentialOnly = TRUE;
	break;
 
	/* help */
 
      case 'h':
	if (MyProcessorNumber == ROOT_PROCESSOR)
	  PrintUsage(myname);
	my_exit(EXIT_SUCCESS);
	break;

    /* */

      case 'V':
	if (MyProcessorNumber == ROOT_PROCESSOR) {
	  WriteConfigure(stdout);
	  auto_show_compile_options();
	}
	my_exit(EXIT_SUCCESS);
	break;
 
	/* Information output */
 
      case 'i':
	InformationOutput = TRUE;
	break;
 
	/* level of region selection. */
 
      case 'l':
	if (--argc > 0) {
	  if (sscanf((*++argv), "%"ISYM, &RegionLevel) != 1) {
	    if (MyProcessorNumber == ROOT_PROCESSOR)
	    ENZO_VFAIL("%s: error reading level.\n", myname)
	  }
	  while (*(argv[0]+1))
	    ++argv[0];
	}
	else {
	  if (MyProcessorNumber == ROOT_PROCESSOR)
	  ENZO_VFAIL("%s: Need to specify level.\n", myname)
	}
	break;
 
	/* Smooth projection */
 
      case 'm':
	ProjectionSmooth = TRUE;
	break;

	/* Write smoothed dark matter field */

      case 'M':
	SmoothedDarkMatterOnly = TRUE;
	break;
 
	/* Output as particle data */
 
      case 'o':
	OutputAsParticleData = TRUE;
	break;
 
	/* Project to a plane. */
 
      case 'p':
	if (--argc > 0) {
	  if (sscanf((*++argv), "%"ISYM, &ProjectionDimension) != 1) {
	    if (MyProcessorNumber == ROOT_PROCESSOR)
	    ENZO_VFAIL("%s: error reading ProjectionDimension.\n",
		      myname)
	  }
	  while (*(argv[0]+1))
	    ++argv[0];
	}
	else {
	  if (MyProcessorNumber == ROOT_PROCESSOR)
	  ENZO_VFAIL("%s: Need to specify dimension.\n", myname)
	}
	project = 1;
	break;

      case 'P':
	if (--argc > 0) {
	  if (sscanf((*++argv), "%"ISYM, &ProjectionDimension) != 1) {
	    if (MyProcessorNumber == ROOT_PROCESSOR)
	    ENZO_VFAIL("%s: error reading ProjectionDimension.\n",
		      myname)
	  }
	  while (*(argv[0]+1))
	    ++argv[0];
	}
	else {
	  if (MyProcessorNumber == ROOT_PROCESSOR)
	  ENZO_VFAIL("%s: Need to specify dimension.\n", myname)
	}
	project = 2;
	break;
 
	/* restart file. */
 
      case 'r':
	restart = TRUE;
	break;
 
	/* get start of region selection. */
 
      case 's':
	dim = 0;
	while (dim < MAX_DIMENSION && argc > 1 && hide_isdigit(*argv[1])) {
	  argc--;
	  if (sscanf((*++argv), "%"ISYM, &RegionStart[dim++]) != 1) {
	    if (MyProcessorNumber == ROOT_PROCESSOR)
	    ENZO_VFAIL("%s: error reading Start indexes.\n", myname)
	  }
	  while (*(argv[0]+1))
	    ++argv[0];
	}
	break;
 
	/* Extract section. */
 
      case 'x':
	extract = TRUE;
	break;
 
      case 'v':
	velanyl = TRUE;
	break;

	/* Unknown */
 
      default:
	if (MyProcessorNumber == ROOT_PROCESSOR)
	ENZO_VFAIL("%s: unknown command-line option: -%s.\n",myname,&c)
	
      } // end of switch(c)
 
  /* Error check for number of parameters, and set parameter file. */
 
  if (argc != 1) {
    if (MyProcessorNumber == ROOT_PROCESSOR)

      PrintUsage(myname);
    my_exit(EXIT_SUCCESS);
  }
  *ParameterFile = argv[0];
 
  return SUCCESS;
}
 
 
/* Explain how to run the program. */
 
void PrintUsage(char *myname)
{
  fprintf(stderr, "usage: %s [options] param_file\n"
	          "   options are:\n"
	          "      -d(ebug)\n"
	          "      -r(estart)\n"
	          "      -x(extract)\n"
	          "         -l(evel_of_extract) level\n"
	          "      -p(roject_to_plane) dimension\n"
	          "      -P(roject_to_plane version 2) dimension\n"
                  "         -m(smooth projection)\n"
	          "      -o(utput as particle data)\n"
	          "      -g (Write Potential field only)\n"
	          "      -M (Write smoothed DM field only)\n"
	          "      -F(riends-of-friends halo finder only)\n"
	          "      -C(ooling time write only)\n"
                  "      -h(elp)\n"
	          "      -i(nformation output)\n"
	          "      -V (show compiler options and flags)\n"
	          "      -s(tart  index region) dim0 [dim1] [dim2]\n"
	          "      -e(nd    index region) dim0 [dim1] [dim2]\n"
	          "      -b(egin  coordinate region) dim0 [dim1] [dim2]\n"
	          "      -f(inish coordinate region) dim0 [dim1] [dim2]\n"
          ,myname);
}
