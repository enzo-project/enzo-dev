/***********************************************************************
/
/  INTERPRET COMMAND-LINE OPTIONS
/
/  written by: Greg Bryan
/  date:       June, 1997
/  modified1:  Robert Harkness
/  date:       November, 2003
/  modified2:  Stephen Skory
/  date        July, 2007
/
/  PURPOSE:
/
************************************************************************/
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
 
// function prototypes
 
void PrintUsage(char *myname);
 
 
 
 
int InterpretCommandLine(int argc, char *argv[], char *myname,
						int *debug, int *nostars, int *nopacked, char *inputfilename[], float *HopDensityThreshold)
{
 
  // Interpret command-line arguments
 
  char c;
 
  while (--argc > 0 && (*++argv)[0] == '-')
    while (c = *++argv[0])
      switch (c) {
 
	// debug
 
      case 'd':
	*debug = 1;
	break;
 
	// debug
 
      case 's':
	    *nostars = 0;
	    fprintf(stderr, "yes to stars set\n");
	  break;
	  
	  case 'p':
	  	*nopacked = 0;
	  	fprintf(stderr, "packed HDF set\n");
	  break;
 
	// Unknown
 
      default:
	fprintf(stderr, "%s: unknown command-line option: -%s.\n", myname, &c);
	return 1;
	
      } // end of switch(c)
 
  // Error check for number of parameters, and set parameter file
 
  if (argc != 2) {
    PrintUsage(myname);
    exit(1);
  }
 
  *HopDensityThreshold = atof(argv[0]);
  *inputfilename = argv[1];
 
  return 0;
}
 
 
// Explain how to run the program
 
void PrintUsage(char *myname)
{
  fprintf(stderr, "usage: %s [options] HopDensityThreshold enzo_restart_file\n"
              "example: %s -p 80.0 data0100\n"
	          "   options are:\n"
	          "      -d(ebug)\n"
	          "      -s (include stars, default is dm-only grouping)\n"
	          "      -p (for packed HDF enzo files)\n"
          ,myname,myname);
}
