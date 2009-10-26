//----------------------------------------------------------------------
//
// File: jb-counter.C
//
// Usage: jb-counter counter-name1 [counter-name2 ...]
//
// Description: Filter counters from records, leaving only those named
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
  int i,j;
  vector<int> filtered_counter_indices;
  map<string,string> Globals;
  vector<string> Attribute_names;
  vector<string> Counter_names;
  vector<string> Counter_types;
  vector<string> Counter_names_new;
  vector<string> Counter_types_new;
  string augregion;
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

  filtered_counter_indices.clear();
  Counter_names_new.clear();
  Counter_types_new.clear();
  for (i=1; i<argc; i++) {
    for (j=0; j<Counter_names.size(); j++) {
      if (Counter_names[j] == argv[i]) {
	filtered_counter_indices.push_back(j);
	Counter_names_new.push_back(Counter_names[j]);
	Counter_types_new.push_back(Counter_types[j]);
      }
    }
  }
  
  jb_print_header (Attribute_names, Counter_names_new, Counter_types_new);

  // Now print records with only given counters

  while (jb_read_next_record(fp,Attribute_names.size(),Counter_names.size(),
			     augregion, counters)) {
    counters_new.clear();
    for (i=0; i<filtered_counter_indices.size(); i++) {
      counters_new.push_back(counters[filtered_counter_indices[i]]);
    }
    jb_print_record (augregion, counters_new);
  }

  jb_log("stop", argv[0]);
}
