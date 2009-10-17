//----------------------------------------------------------------------
//
// File: jb-columns.C
//
// Usage: jb-columns "x_1 [x_2 ... x_N]"
//
// Description: Prints the given attributes or counter values in columns
// Description: that are easy for plotting programs to read (e.g. gnuplot)
//
// Caveats: The output is unsorted, so should be run through "sort -n"
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

  vector<int> indices;
  vector<bool> isAttribute;
  vector<bool> isCounter;
  
  // Log start of jb-merge

  jb_log("start", argv[0]);

  FILE *fp = stdin;

  // Read globals and data header

  jb_read_globals(fp, Globals);
  jb_read_header(fp, Attribute_names, Counter_names, Counter_types);

  indices.clear();
  isAttribute.clear();
  isCounter.clear();

  // Find indices of columns, and whether they're Attributes or Counters

  for (i=0; i<argc-1; i++) {
    bool found = false;
    char * column = argv[i+1];
    int j;
    if (strcmp(column,"region")==0) {
       isAttribute.push_back(1);
       isCounter.push_back  (0);
       indices.push_back    (-1);
       found = true;
    }
    for (j=0; j<Attribute_names.size(); j++) {
      if (Attribute_names[j] == column) {
        isAttribute.push_back(1);
        isCounter.push_back  (0);
        indices.push_back    (j);
        found = true;
        break;
      }
    }
    for (j=0; j<Counter_names.size(); j++) {
      if (Counter_names[j] == column) {
        isAttribute.push_back(0);
        isCounter.push_back  (1);
        indices.push_back    (j);
        found = true;
        break;
      }
    }
    if (! found) {
      fprintf (stderr, "%s Error: Attribute or Counter name '%s' not found!\n",
               argv[0],column);
      exit(1);
    }
  }
        
  while (jb_read_next_record(fp,Attribute_names.size(),Counter_names.size(),
                             augregion, counters)) {

    jb_split (augregion, augregion_split);

    // Determine whether to print the record
    
    string line = "";

    for (int i=0; i<argc-1; i++) {
      string value;
      if (isAttribute[i]) {
        value = augregion_split[indices[i]+1];
      } else if (isCounter[i]) {
        value = counters[indices[i]];
      } else {
        fprintf (stderr, "%s Error: index is neither Attribute nor Counter name!\n",
                argv[0]);
        exit(1);
      }

      line = (i==0) ? value : line + " " + value;
    }
    printf ("%s\n",line.c_str());

    // Get next record
  }
  jb_log("stop", argv[0]);
}
