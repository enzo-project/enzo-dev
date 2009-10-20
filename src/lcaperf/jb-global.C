//----------------------------------------------------------------------
//
// File: jb-global.C
//
// Usage: jb-global <global1> <value1> [<global2> <value2> ... ]
//
// Description: (Re)defines the listed global variables.  Defining a global to be "*"
// Description: serves to undefine it (remember the quotes!)
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

  if (argc == 2) {

    // Print existing global value

    if (Globals.find(argv[1]) != Globals.end()) {
      printf ("%s\n",Globals[argv[1]].c_str());
    }
      
  } else {
    for (int i=1; i<argc; i+=2) {
      const char * new_global = argv[i];
      const char * new_value = argv[i+1];
      if (strcmp(new_value,"*") != 0) {
	Globals[new_global] = new_value;
      } else {
	if (Globals.find(new_global) != Globals.end()) {
	  Globals.erase(new_global);
	}
      }
    }
    jb_print_globals (Globals);

    // Read and write data header

    jb_read_header(fp, Attribute_names, Counter_names, Counter_types);
    jb_print_header (Attribute_names, Counter_names, Counter_types);

    // Print records

    while (jb_read_next_record(fp,Attribute_names.size(),Counter_names.size(),
			       augregion, counters)) {
      jb_split (augregion, augregion_split);
      jb_print_record (augregion, counters);
    }
  }
  jb_log("stop", argv[0]);
}
