//----------------------------------------------------------------------
//
// File: jb.h
//
// Description: Descriptiors for functions used by jb-*.C utilities.  See jb.C
// Description: for the function implementations.
//
//----------------------------------------------------------------------
//
// Copyright 2005 James Bordner
// Copyright 2005 Laboratory for Computational Astrophysics
// Copyright 2005 Regents of the University of California
//
//----------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include <map>
#include <string>
#include <vector>

using namespace std;

#define MAX_LINE_LENGTH 80
#define MAX_WORD_LENGTH 40
#define MAX_NUM_WORDS 10
#define TRACE printf ("%s:%d TRACE\n",__FILE__,__LINE__); fflush(stdout);

// float constants for "*" and "#" equivalent.  Using these should not be
// permanent.

#define STAR_FLOAT (-999999)
#define HASH_FLOAT (-999998)
#define STAR_FLOAT_STRING "-999999"
#define HASH_FLOAT_STRING "-999998"

#define _TRACE_ printf ("%s:%d TRACE\n",__FILE__,__LINE__); fflush(stdout);


//----------------------------------------------------------------------

typedef vector<string> vecstr;

struct reg_struct {
  string sorted;
  string original;
  double y;
  double x;
} ;

//----------------------------------------------------------------------

void jb_log(const char *, const char *);
FILE * jb_file_open (string);
FILE * jb_file_open_stdin ();
void jb_file_close (FILE *);
int  jb_read_line(FILE *, string & line, char words[][MAX_WORD_LENGTH]);

void jb_read_globals(FILE *, map<string,string> & globals);
void jb_print_globals (map<string,string> Globals);

void jb_read_header(FILE *, 
		    vector<string> & Attribute_name, 
		    vector<string> & Counter_name, 
		    vector<string> & Counter_type);
void jb_read_header0(FILE *, 
		    vector<string> & Attribute_name, 
		    vector<string> & Counter_name, 
		    vector<string> & Counter_type);
void jb_read_header1(FILE *, 
		    vector<string> & Attribute_name, 
		    vector<string> & Counter_name, 
		    vector<string> & Counter_type);

void jb_print_header (vector<string> & Attribute_name, 
		      vector<string> & Counter_name, 
		      vector<string> & Counter_type);
void jb_print_header0 (vector<string> & Attribute_name, 
		      vector<string> & Counter_name, 
		      vector<string> & Counter_type);
void jb_print_header1 (vector<string> & Attribute_name, 
		      vector<string> & Counter_name, 
		      vector<string> & Counter_type);

bool jb_read_next_record (FILE *, 
			 int num_attributes, 
			 int num_counters, 
			 string & augregion,
			 vector<string> & record_counters);

void jb_read_records (FILE *, 
		      int num_attributes, 
		      int num_counters, 
		      map<string,vecstr> & Regions);

void jb_print_records (map<string,vecstr> & Records);
void jb_print_record (const string & augregion, const vector<string> & counters);

void jb_split (string augregion, vector<string> & augregion_split);
void jb_counters_merge (vector<string> counters1, 
			vector<string> counters2, 
			vector<string> & new_counters);


void jb_augregion_rotate_ (string augregion, 
			   int index, 
			   string & augregion_rotated);

string jb_augregion_root_ (string augregion);



