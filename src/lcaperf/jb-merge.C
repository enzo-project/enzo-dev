//----------------------------------------------------------------------
//
// File: jb-merge.C
//
// Usage: jb-merge <file1> [ <file2> ... ] 
//
// Description: Merges the listed jbPerf files into a single file.
// Description: Promotes globals to attributes as needed.

// Assumptions: Merged files that promote globals must have different
// Assumptions: global values for each file, e.g. processor-rank.  If no globals
// Assumptions: are promoted, then a "slow" version is used that accumulates
// Assumptions: records before writing them.
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

  // Exit if there is nothing to do

  if (argc == 1) exit(1);

  // Log start of jb-merge

  jb_log("start", argv[0]);

  // Get command-line arguments
  
  vector<string> infiles;
  infiles.resize(argc-1);
  for (i=1; i<argc; i++) infiles[i-1] = argv[i];

  // (DEBUG: print out command-line arguments)
  
  map<string,string> old_globals;
  map<string,string> new_globals;
  map<string,string> curr_globals;
  map<string,string> Globals; // @@@ Might not need--can re-use curr_globals?
  map<string,string> new_promoted_attributes;
  map<string,int> new_attribute_index;
  map<string,int> new_counter_index;
  vector<string> new_counter_names;
  vector<string> new_counter_types;

  vector<string> Attribute_names;
  vector<string> Counter_names;
  vector<string> Counter_types;
  map<int,int> attribute_remap;
  map<int,int> counter_remap;
  map<string,vecstr> new_records;
  string old_augregion;
  vector<string> old_record_counters;
  vector<string> augregion_split;
  vector<string> new_augregion_split;
  string new_augregion;
  string attribute_value;
  string counter_value;
  vector<string> new_counter_values;

  bool are_disjoint = false;

  //------------------------------------------------------------------------
  // READ GLOBALS AND HEADER RECORD
  //------------------------------------------------------------------------

  int ifile=0;
  for (ifile=0; ifile<infiles.size(); ifile++) {

    string infile = infiles[ifile];

    // First figure out which globals to promote to locals
    
    FILE *fp = jb_file_open (infile);

    // Read globals

    jb_read_globals(fp, curr_globals);

    // if first file, set new and old globals to its globals...
    
    if (infile==infiles[0]) {

      old_globals = curr_globals;
      new_globals = curr_globals;

    // ...otherwise update old and new globals accordingly

    } else {

      if (old_globals.size() != curr_globals.size()) {
	printf ("ERROR: '%s' and '%s' have different numbers of globals--exiting!",
		infiles[0].c_str(), infile.c_str());
      } else {
	map<string,string>::iterator global_var;
	for (global_var = curr_globals.begin();
	     global_var != curr_globals.end();
	     ++global_var) {
	  if (strcmp(curr_globals[global_var->first].c_str(),
		     old_globals[global_var->first].c_str()) != 0 &&
	      new_globals.find(global_var->first) != new_globals.end()) {
	    new_globals.erase(global_var->first);
	    new_promoted_attributes[global_var->first] = curr_globals[global_var->first];
	  }
	}
      }
    }

    // Which attributes...

    jb_read_header(fp, Attribute_names, Counter_names, Counter_types);

    int i;

    for (i=0; i<Attribute_names.size(); i++) {
      if (new_attribute_index.find(Attribute_names[i]) 
	  == new_attribute_index.end()) {
	new_attribute_index[Attribute_names[i]] = new_attribute_index.size();
      }
    }

    // And which counter names and types ...

    for (i=0; i<Counter_names.size(); i++) {
      if (new_counter_index.find(Counter_names[i]) 
	  == new_counter_index.end()) {
	new_counter_index[Counter_names[i]] = new_counter_index.size() - 1;
	new_counter_types.push_back(Counter_types[i]);
	new_counter_names.push_back(Counter_names[i]);
      }
    }

    jb_file_close(fp);

  }

  
  //------------------------------------------------------------------------
  // DETERMINE WHETHER RECORDS ARE DISJOINT OR NOT (SEE ASSUMPTIONS IN HEADER COMMENT!)
  //------------------------------------------------------------------------

  are_disjoint = (new_promoted_attributes.size() != 0);

  //------------------------------------------------------------------------
  // WRITE GLOBAL RECORD
  //------------------------------------------------------------------------

  jb_print_globals (new_globals);

  //------------------------------------------------------------------------
  // WRITE HEADER RECORD
  //------------------------------------------------------------------------

  Attribute_names.clear();

  //    region and existing attributes

  for (map<string,int>::iterator pattribute = new_attribute_index.begin();
       pattribute != new_attribute_index.end();
       ++pattribute) {
    Attribute_names.push_back (pattribute->first);
  }
  //   new attributes promoted from globals
  for (map<string,string>::iterator pattribute = new_promoted_attributes.begin();
       pattribute != new_promoted_attributes.end();
       ++pattribute) {
    Attribute_names.push_back (pattribute->first);
  }
  //    counters
  Counter_names.clear();
  Counter_types.clear();
  for (map<string,int>::iterator pcounter = new_counter_index.begin();
       pcounter != new_counter_index.end();
       ++pcounter) {
    Counter_types.push_back (new_counter_types[pcounter->second]);
    Counter_names.push_back (pcounter->first);
  }

  jb_print_header (Attribute_names, Counter_names, Counter_types);

  //------------------------------------------------------------------------
  // READ DATA RECORDS
  //------------------------------------------------------------------------

  new_records.clear();

  for (ifile=0; ifile<infiles.size(); ifile++) {

    string infile = infiles[ifile];

    FILE *fp = jb_file_open (infile);

    // Read globals

    jb_read_globals(fp, Globals);

    // Read header

    jb_read_header(fp, Attribute_names, Counter_names, Counter_types);

    // map attributes from new-order to file-order

    
    attribute_remap.clear();
    int num_attributes = Attribute_names.size();
    for (i=0; i<num_attributes; i++) {
      attribute_remap [new_attribute_index[Attribute_names[i]]] = i + 1;
    }

    // map counters from new-order to file-order

    counter_remap.clear();
    int num_counters = Counter_names.size();
    for (i=0; i<num_counters; i++) {
      counter_remap [new_counter_index[Counter_names[i]]] = i;
    }

    // Read data records

    while (jb_read_next_record(fp,num_attributes,num_counters,
			       old_augregion, old_record_counters)) {
    
      // Generate new augmented region for current record

      jb_split (old_augregion, augregion_split);

      new_augregion_split.clear();
      new_augregion_split.push_back(augregion_split[0]);
      
      for (map<string,int>::iterator pattribute = new_attribute_index.begin();
	   pattribute != new_attribute_index.end();
	   ++ pattribute) {
	int index = pattribute->second;
	if (attribute_remap.find(index) != attribute_remap.end()) {
	  attribute_value = augregion_split[attribute_remap[index]];
	} else {
	  attribute_value = "*";
	}
	new_augregion_split.push_back(attribute_value);
      }

      for (map<string,string>::iterator pattribute = new_promoted_attributes.begin();
	   pattribute != new_promoted_attributes.end();
	   ++pattribute) {
	new_augregion_split.push_back(Globals[pattribute->first]);
      }

      // (un-split new_augregion_split)

      new_augregion = new_augregion_split[0];

      for (i=1; i<new_augregion_split.size(); i++) {
	new_augregion = new_augregion + " " + new_augregion_split[i];
      }

      // Generate new counters for current record

      new_counter_values.clear();
      for (map<string,int>::iterator pcounter = new_counter_index.begin();
	   pcounter != new_counter_index.end();
	   ++ pcounter) {
	int index = pcounter->second;
	if (counter_remap.find(index) != counter_remap.end()) {
	  counter_value = old_record_counters[counter_remap[index]];
	} else {
	  counter_value = "*";
	}
	new_counter_values.push_back(counter_value);
      }

      // merge new augmented region and counter values those accumulated
      // so far

      if (are_disjoint) {

	jb_print_record (new_augregion, new_counter_values);

      } else {

	if (new_records.find(new_augregion) != new_records.end()) {
	  jb_counters_merge 
	    (new_records[new_augregion], new_counter_values,new_counter_values);
	}
      }
      
      // @@@ new_augregion too long
      new_records[new_augregion] = new_counter_values;

    } // while record

    jb_file_close(fp);

  } // for files

  //------------------------------------------------------------------------
  // WRITE DATA RECORDS
  //------------------------------------------------------------------------

  if (! are_disjoint) {

    for (map<string,vecstr>::iterator augregion = new_records.begin();
	 augregion != new_records.end();
	 ++augregion) {

      jb_split (augregion->first, augregion_split);
      for (i=0; i<augregion_split.size(); i++) {
	printf ("%s\n",augregion_split[i].c_str());
      }

      for (i=0; i<(new_records[augregion->first]).size(); i++) {
	printf ("%s\n",((new_records[augregion->first])[i]).c_str());
      }

      printf ("\n");
    }

  }

  jb_log("stop", argv[0]);

}
