//----------------------------------------------------------------------
//
// Usage: jb-delete <name1> [<name2> ...]
//
// Description: Deletes the attribute, global, or counter
//
// Bugs: 1. Attributes can be deleted such that augmented regions
// Bugs:    are no longer unique, resulting in an ill-formed file
// Bugs: 2. "region" cannot be deleted:  currently it's hard-coded in jb.C
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
  vector<string> Attribute_names_new;
  vector<string> Counter_names;
  vector<string> Counter_types;
  vector<string> Counter_names_new;
  vector<string> Counter_types_new;
  string augregion;
  string new_augregion;
  vector<string> augregion_split;
  vector<string> counters;
  vector<string> new_counters;

  map<string,int> delete_name;
  vector<int> attribute_indices;
  vector<int> counter_indices;

  // Log start

  jb_log("start", argv[0]);

  // Read from stdin or named file
  
  FILE *fp = stdin;

  // Check number of arguments: exit if 0

  if ((argc-1) == 0) {
    fprintf (stderr,"%s:%d ERROR: wrong number of arguments!\n",
	     __FILE__,__LINE__);
    exit(1);
  }

  for (i=1; i<argc; i++) {
    delete_name[argv[i]] = 1;
  }

  //--------------------------------------------------------------------
  // Globals
  //--------------------------------------------------------------------

  // Read Globals

  jb_read_globals(fp, Globals);


  // Delete Globals

  map<string,string>::iterator global_var;

  for (global_var = Globals.begin();
       global_var != Globals.end();
       ++global_var) {
    if (delete_name[global_var->first] != 1) {
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
    if (delete_name[Attribute_names[i]] != 1) {
      attribute_indices.push_back(i);
      Attribute_names_new.push_back(Attribute_names[i]);
    }
  }

  // Remap Header: counters

  for (i=0; i<Counter_names.size(); i++) {
    if (delete_name[Counter_names[i]] != 1) {
      counter_indices.push_back(i);
      Counter_names_new.push_back(Counter_names[i]);
      Counter_types_new.push_back(Counter_types[i]);
    }
  }

  // Write Header

  jb_print_header (Attribute_names_new, Counter_names_new, Counter_types_new);

  //--------------------------------------------------------------------
  // Records
  //--------------------------------------------------------------------

  while (jb_read_next_record(fp,num_attributes,num_counters, augregion, counters)) {

    // Split augregion --> augregion_split

    jb_split (augregion, augregion_split);

    // Remap augregion_split

    new_augregion = augregion_split[0];

    for (i=0; i<attribute_indices.size(); i++) {
      new_augregion = new_augregion + " " + augregion_split[attribute_indices[i]+1];
    }

    new_counters.clear();
    for (i=0; i<counter_indices.size(); i++) {
      new_counters.push_back(counters[counter_indices[i]]);
    }

    // Write record

    jb_print_record (new_augregion, new_counters);
  }


  jb_log("stop", argv[0]);
}

