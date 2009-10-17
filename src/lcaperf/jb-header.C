//----------------------------------------------------------------------
//
// File: jb-header.C
//
// Usage: jb-header
//
// Description: Displays the header information (header record)
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

  // Read and write data header

  jb_read_header(fp, Attribute_names, Counter_names, Counter_types);
  jb_print_header (Attribute_names, Counter_names, Counter_types);

  jb_log("stop", argv[0]);
}
