//----------------------------------------------------------------------
//
// File: jb-region.C
//
// Usage: jb-region region1 [region2 ...]
//
// Description: Filters records, leaving only those with the specified region
// Description: names.
//
//----------------------------------------------------------------------
//
// Copyright 2005 James Bordner
// Copyright 2005 Laboratory for Computational Astrophysics
// Copyright 2005 Regents of the University of California
//
//----------------------------------------------------------------------

#include "jb.h"

main(int argc, char **argv)

{

  map<string,string> Globals;
  vector<string> Attribute_names;
  vector<string> Counter_names;
  vector<string> Counter_types;
  string augregion;
  vector<string> augregion_split;
  vector<string> counters;

  // Log start of jb-merge

  jb_log("start", argv[0]);

  FILE *fp = stdin;

  // Read and write globals

  jb_read_globals(fp, Globals);
  jb_print_globals (Globals);

  // Read and write data header

  jb_read_header(fp, Attribute_names, Counter_names, Counter_types);
  jb_print_header (Attribute_names, Counter_names, Counter_types);

  // Now print only selected records based on region names:

  while (jb_read_next_record(fp,Attribute_names.size(),Counter_names.size(),
			     augregion, counters)) {

    jb_split (augregion, augregion_split);

    // Determine whether to print the record
    
    for (int i=1; i<argc; i++) {
      if (strcmp(argv[i],augregion_split[0].c_str())==0) {
	jb_print_record (augregion, counters);
	break;
      }
    }

    // Get next record
  }
  jb_log("stop", argv[0]);
}
