/***********************************************************************
/
/  OUTPUT CONFIGURATION
/
/  written by: James Bordner
/  date:       July 2008
/  modified1:
/
/  PURPOSE: Write configuration and version settings to data dump file
/
************************************************************************/
 
#include <stdio.h>

void auto_show_version (FILE *fp);
void auto_show_config (FILE *fp);
void auto_show_flags (FILE *fp);
 
#define DIVIDER "========================================================================\n"

void WriteConfigure(FILE *fp)
{
  // Print output of "make show-version"

  fprintf (fp, DIVIDER "make show-version\n" DIVIDER);
  auto_show_version (fp);

  // Print output of "make show-config"

  fprintf (fp, DIVIDER "gmake show-config\n" DIVIDER);
  auto_show_config (fp);

  // Print output of "make show-flags"

  fprintf (fp, DIVIDER "gmake show-flags\n" DIVIDER);
  auto_show_flags (fp);
}
