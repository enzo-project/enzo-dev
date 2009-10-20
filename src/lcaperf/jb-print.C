//----------------------------------------------------------------------
//
// Usage: jb-print < <file>
// Usage: jb-print <file>
//
//----------------------------------------------------------------------
//
// Copyright 2004 James Bordner
// Copyright 2004 Laboratory for Computational Astrophysics
// Copyright 2004 Regents of the University of California
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

  // Log start

  jb_log("start", argv[0]);

  // Read from stdin or named file
  
  FILE *fp = (argc == 1) ? stdin : jb_file_open(argv[1]);

  // Read globals and header

  jb_read_globals(fp, Globals);
  jb_read_header(fp, Attribute_names, Counter_names, Counter_types);

  int num_attributes = Attribute_names.size();
  int num_counters   = Counter_names.size();

  while (jb_read_next_record(fp,num_attributes,num_counters, augregion, counters)) {

    printf ("\n");

    jb_split (augregion, augregion_split);

    printf ("Region %s\n", augregion_split[0].c_str());
    for (i=0; i<num_attributes; i++) {
      printf ("Attribute %s %s\n", Attribute_names[i].c_str(), augregion_split[i+1].c_str());
    }
    for (i=0; i<num_counters; i++) {
      printf ("Counter %s %s\n", Counter_names[i].c_str(), counters[i].c_str());
    }
  }


  jb_log("stop", argv[0]);
}

