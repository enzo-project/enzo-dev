//----------------------------------------------------------------------
//
// Usage: jb-rename <name1-from> <name1-to> [<name*-from> <name*-to>] < <file>
//
// Description: Renames globals, attributes, and counters.
//
// Bugs: 1. Does not do error-checking--it doesn't check whether 
// Bugs:    an attribute etc. already exists with the new name.
// Bugs: 2. "region" does not get remapped: currently it's hard-coded in jb.C
// Bugs:    for file type backwards compatibility. 
//
//----------------------------------------------------------------------
//
// Copyright 2006 James Bordner
// Copyright 2006 Laboratory for Computational Astrophysics
// Copyright 2006 Regents of the University of California
//
//----------------------------------------------------------------------

#include "jb.h"

main(int argc, char **argv)
{
  int i,j;
  vector<int> filtered_counter_indices;
  map<string,string> Globals;
  map<string,string> Globals_new;
  vector<string> Attribute_names;
  vector<string> Counter_names;
  vector<string> Counter_types;
  vector<string> Counter_names_new;
  vector<string> Counter_types_new;
  string augregion;
  string new_augregion;
  vector<string> augregion_split;
  vector<string> counters;

  map<string,string> map_name;

  // Log start

  jb_log("start", argv[0]);

  // Read from stdin or named file
  
  FILE *fp = stdin;

  // Check number of arguments: exit if 0 or odd number

  if ((argc-1) == 0 || (argc-1) % 2 != 0) {
    fprintf (stderr,"%s:%d ERROR: wrong number of arguments!\n",
	     __FILE__,__LINE__);
    exit(1);
  }

  for (i=1; i<argc-1; i+=2) {
    map_name[argv[i]] = argv[i+1];
  }

  //--------------------------------------------------------------------
  // Globals
  //--------------------------------------------------------------------

  // Read Globals

  jb_read_globals(fp, Globals);

  // Remap Globals --> Globals_new

  map<string,string>::iterator global_var;

  for (global_var = Globals.begin();
       global_var != Globals.end();
       ++global_var) {
    if (map_name[global_var->first] != "") {
      Globals_new[map_name[global_var->first]] = global_var->second;
    } else {
      Globals_new[global_var->first] = global_var->second;
    }
  }

  // Write Globals_new

  jb_print_globals(Globals_new);

  //--------------------------------------------------------------------
  // Header
  //--------------------------------------------------------------------

  // Read Header

  jb_read_header(fp, Attribute_names, Counter_names, Counter_types);
  int num_attributes = Attribute_names.size();
  int num_counters   = Counter_names.size();

  // Remap Header: attributes

  for (i=0; i<Attribute_names.size(); i++) {
    if (map_name[Attribute_names[i]] != "") {
      Attribute_names[i] = map_name[Attribute_names[i]];
    }
  }

  // Remap Header: counters

  for (i=0; i<Counter_names.size(); i++) {
    if (map_name[Counter_names[i]] != "") {
      Counter_names[i] = map_name[Counter_names[i]];
    }
  }

  // Write Header

  jb_print_header (Attribute_names, Counter_names, Counter_types);

  //--------------------------------------------------------------------
  // Records
  //--------------------------------------------------------------------

  while (jb_read_next_record(fp,num_attributes,num_counters, augregion, counters)) {

    // Split augregion --> augregion_split

    jb_split (augregion, augregion_split);

    // Remap augregion_split

    for (i=0; i<augregion_split.size(); i++) {
      if (map_name[augregion_split[i]] != "") {
	augregion_split[i] = map_name[augregion_split[i]];
      }
    }

    new_augregion = augregion_split[0];

    for (i=1; i<augregion_split.size(); i++) {
      new_augregion = new_augregion + " " + augregion_split[i];
    }

    // Write record

    jb_print_record (new_augregion, counters);
  }


  jb_log("stop", argv[0]);
}

