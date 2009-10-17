//----------------------------------------------------------------------
//
// File: jb-attribute.C
//
// Usage: jb-attribute  attribute-name1 attribute-value1
// Usage:              [attribute-name2 attribute-value2 ...]
//
// Description: Filters records according to attribute values
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
  vector<string> Attribute_names;
  vector<string> Counter_names;
  vector<string> Counter_types;
  vector<string> Counter_names_new;
  vector<string> Counter_types_new;
  string augregion;
  vector<string> augregion_split;
  vector<string> counters;
  vector<string> counters_new;

  // Log start

  jb_log("start", argv[0]);

  FILE *fp = stdin;

  // Read and write globals

  jb_read_globals(fp, Globals);
  jb_print_globals (Globals);


  // Read and write data header

  jb_read_header(fp, Attribute_names, Counter_names, Counter_types);
  jb_print_header (Attribute_names, Counter_names, Counter_types);

  // First determine Attribute indices 

  vector<int>    indices;
  vector<string> values;

  for (i=1; i<argc; i+=2) {

    string argv_attribute = argv[i];
    string argv_value     = argv[i+1];

    bool any_found = false;
    for (j=0; j<Attribute_names.size(); j++) {
      if (argv_attribute == Attribute_names[j]) {
	indices.push_back(j);
	values.push_back(argv_value);
	any_found = true;
	break;
      }
    }
    if (! any_found) {
      fprintf (stderr,"Error: Attribute '%s' not in file!\n",argv_attribute.c_str());
    }
  }

  // Now print only selected records based on attribute values

  while (jb_read_next_record(fp,Attribute_names.size(),Counter_names.size(),
			     augregion, counters)) {

    jb_split (augregion, augregion_split);

    bool print_record = true;

    for (i=0; i<indices.size(); i++) {
      if (augregion_split[indices[i]+1] != values[i]) print_record = false;
    }

    if (print_record) jb_print_record (augregion, counters);
  }
  jb_log("stop", argv[0]);
}
