//----------------------------------------------------------------------
//
// File: jb-promote.C
//
// Usage: jb-promote <global1> [ <global2> ... ] 
//
// Description: Converts the listed global variables to attributes in all
// Description: records
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
  int i;

  map<string,string> Globals;
  vector<string> Attribute_names;
  vector<string> Counter_names;
  vector<string> Counter_types;
  string augregion;
  vector<string> augregion_split;
  vector<string> counters;

  vector<string> new_attributes;
  string attribute_string;
  // Log start of jb-merge

  jb_log("start", argv[0]);

  FILE *fp = stdin;

  // Read and write globals

  jb_read_globals(fp, Globals);

  new_attributes.clear();
  attribute_string = "";
  for (i=1; i<argc; i++) {
    if (Globals.find(argv[i]) != Globals.end()) {
      new_attributes.push_back(argv[i]);
      attribute_string = attribute_string + " " + Globals[argv[i]];
      Globals.erase(argv[i]);
    }
  }

  jb_print_globals (Globals);

  // Read and write data header

  jb_read_header(fp, Attribute_names, Counter_names, Counter_types);

  int num_attributes = Attribute_names.size();
  for (i=0; i<new_attributes.size(); i++) {
    Attribute_names.push_back(new_attributes[i]);
  }

  jb_print_header (Attribute_names, Counter_names, Counter_types);

  // Now print only selected records based on region names:

  while (jb_read_next_record(fp,num_attributes,Counter_names.size(),augregion,counters)) {

    jb_split (augregion, augregion_split);

    // Determine whether to print the record
    
    augregion = augregion + attribute_string;

    jb_print_record (augregion, counters);

    // Get next record
  }
  jb_log("stop", argv[0]);
}
