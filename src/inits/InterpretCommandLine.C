/***********************************************************************
/
/  INTERPRET COMMAND-LINE OPTIONS
/
/  written by: Greg Bryan
/  date:       June, 1997
/  modified1:  Robert Harkness
/  date:       November, 2003
/
/  PURPOSE:
/
************************************************************************/
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
 
#include "macros_and_parameters.h"
#include "global_data.h"
 
// function prototypes
 
void PrintUsage(char *myname);
 
 
 
 
int InterpretCommandLine(int argc, char *argv[], char *myname,
			 char *ParameterFile[], char *SubGridParameterFile[])
{
 
  // Interpret command-line arguments
 
  char c;
 
  while (--argc > 0 && (*++argv)[0] == '-')
    while (c = *++argv[0])
      switch (c) {
 
	// debug
 
      case 'd':
	debug = TRUE;
	break;
 
	// debug
 
      case 's':
	if (--argc > 0) {
	  *SubGridParameterFile = new char[strlen(argv[1])+1];
	  strcpy(*SubGridParameterFile, *++argv);
	  while (*(argv[0]+1))
	    ++argv[0];
	}
	else {
	  fprintf(stderr, "%s: must specify subgrid parameter file\n", myname);
	  return FAIL;
	}
	break;
 
	// Unknown
 
      default:
	fprintf(stderr, "%s: unknown command-line option: -%s.\n", myname, &c);
	return FAIL;
	
      } // end of switch(c)
 
  // Error check for number of parameters, and set parameter file
 
  if (argc != 1) {
    PrintUsage(myname);
    exit(FAIL);
  }
 
  *ParameterFile = argv[0];
 
  return SUCCESS;
}
 
 
// Explain how to run the program
 
void PrintUsage(char *myname)
{
  fprintf(stderr, "usage: %s [options] param_file\n"
	          "   options are:\n"
	          "      -d(ebug)\n"
	          "      -s(ubgrid) param_file\n"
          ,myname);
}
