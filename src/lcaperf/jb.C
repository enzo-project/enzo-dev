//----------------------------------------------------------------------
//
// File: jb.C
//
// Description: Functions used by jb-*.C utilities.  See jb.h for the list
// Description: of functions
//
//----------------------------------------------------------------------
//
// Copyright 2005 James Bordner
// Copyright 2005 Laboratory for Computational Astrophysics
// Copyright 2005 Regents of the University of California
//
//----------------------------------------------------------------------

#include "jb.h"

static map<FILE *,string> global_version;
#define global_version_write "1.0"

//========================================================================
void jb_log(const char * action, const char * name)
//------------------------------------------------------------------------
// Log a message (e.g. "start jb-merge") with the date in the 
// file "JBPERF.log"
//------------------------------------------------------------------------
{
#ifdef LOG_PROGRESS
  time_t timer;
  struct tm *tp;
  time (&timer);
  tp = localtime (&timer);

//   FILE *fp = fopen ("JBPERF.log","a");
//   fprintf (fp,"%-6s %-16s %04d %02d %02d %02d %02d %02d\n", action, name,
// 	   1900+tp->tm_year,
// 	   1+tp->tm_mon,
// 	   tp->tm_mday,
// 	   tp->tm_hour,
// 	   tp->tm_min,
// 	   tp->tm_sec);
//   fclose (fp);

  FILE *fp = fopen ("JBPERF.log","a");
  fprintf (fp,"%02d %02d %02d %-5s        %s\n",tp->tm_hour,tp->tm_min,tp->tm_sec,action,name);
  fclose (fp);
#endif
}
//========================================================================
FILE * jb_file_open (string filename)
{
  return fopen (filename.c_str(),"r");
}
//========================================================================
FILE * jb_file_open_stdin ()
{
  return stdin;
}
//========================================================================
void jb_file_close (FILE *fp)
{
  fclose (fp);
}
//========================================================================
int jb_read_line (FILE * fp, string & line, char words[][MAX_WORD_LENGTH])
{
  int c;
  char cline[MAX_LINE_LENGTH];
  int i=0;  // Index into line string
  int j=0;  // Index into word string
  int c0 = c = getc(fp);

  // Read line

  while ((i < MAX_LINE_LENGTH) && c != EOF && c != '\n') {
    cline[i++] = c;
    c = getc(fp);
  }
  cline[i]='\0';
  int n = i;

  if (i==MAX_LINE_LENGTH) {
    fprintf (stderr,"Line too long!");
    exit(1);
  }

  // Read words

  int iword=0;

  for (i=0; i<n; i++) {
    if (cline[i] == ' ') {
      words[iword][j] = '\0';
      ++iword;
      j=0;
    } else {
      if (j < MAX_WORD_LENGTH) {
	words[iword][j++] = cline[i];
      } else {
	fprintf (stderr,"Word too long!\n");
	words[iword][MAX_WORD_LENGTH-1] = '\0';
	fprintf (stderr,"%s\n",words[iword]);
	exit(1);
      }
    }
  }
  words[iword][j] = '\0';

  for (j=iword+1; j<MAX_NUM_WORDS; j++) words[j][0] = '\0';

  line = cline;
  return c0;
}

//========================================================================
void jb_read_globals(FILE *fp, map<string,string> & Globals)
// Reads the global definitions from the jbPerf data file fp
{
  string line;
  char words[MAX_NUM_WORDS][MAX_WORD_LENGTH];
  int c = jb_read_line (fp,line,words);

  while (c != '\n') {
    // Check that we're really in a valid Globals section
    if (strcmp(words[0],"global")!=0) {
      fprintf (stderr, "Error: Missing 'global' in global record!\n");
      exit(1);
    } else  {
      Globals[words[1]] = words[2];
      for (int i=3; strlen(words[i]) > 0; i++) {
	Globals[words[1]] = Globals[words[1]] + string (" ") + string (words[i]);
      }
    }
    c = jb_read_line (fp,line,words);
  }
  if (Globals.find("lcaperf-version") != Globals.end()) {
    global_version[fp] = Globals["lcaperf-version"];
  } else {
    global_version[fp] = "0.0";
  }
}
//========================================================================
void jb_read_header(FILE * fp, 
		    vector<string> & Attribute_name, 
		    vector<string> & Counter_name, 
		    vector<string> & Counter_type)
{
  if (global_version[fp] == "0.0") {
    jb_read_header0 (fp,Attribute_name, Counter_name, Counter_type);
  } else if (global_version[fp] == "1.0") {
    jb_read_header1 (fp,Attribute_name, Counter_name, Counter_type);
  } else {
    fprintf (stderr, "Error: Unknown version %s in jb_read_header!\n",
	     global_version[fp].c_str());
      exit(1);
  }
}

//========================================================================

void jb_read_header0(FILE * fp, 
		    vector<string> & Attribute_name, 
		    vector<string> & Counter_name, 
		    vector<string> & Counter_type)
{
  string line;
  char words[MAX_NUM_WORDS][MAX_WORD_LENGTH];

  int c = jb_read_line (fp,line,words);

  if (c != '\n') {
    if (strcmp(words[0],"region") != 0) {
      fprintf (stderr, "jb_read_header0: Error in Header section: missing 'region'!\n");
      fprintf (stderr, "%s\n",line.c_str());
      exit(1);
    }
  }
      
  c = jb_read_line (fp,line,words);

  Attribute_name.clear();
  Counter_name.clear();
  Counter_type.clear();

  while (c != '\n') {
    if (words[1][0] == '\0' || words[2][0] != '\0') {
      // not two words in line 
      fprintf (stderr, "jb_read_header0: Error in Header section: wrong line length!\n");
      fprintf (stderr, "%s\n",line.c_str());
      exit(1);
    } else {
      if (strcmp(words[0],"attribute") == 0) {
	Attribute_name.push_back(words[1]); 
      } else if (strcmp(words[0],"basic") == 0 ||
		 strcmp(words[0],"papi") == 0 ||
		 strcmp(words[0],"user") == 0 ||
		 strcmp(words[0],"derived") == 0) {
	Counter_type.push_back(words[0]);
	Counter_name.push_back(words[1]);
      } else {
	fprintf (stderr,"jb_read_header0: Error in Header section: unknown keyword \n%s\n",
		 line.c_str());
      }
    }
    c = jb_read_line (fp,line,words);
  }
}

//========================================================================

void jb_read_header1(FILE * fp, 
		    vector<string> & Attribute_name, 
		    vector<string> & Counter_name, 
		    vector<string> & Counter_type)
{
  string line;
  char words[MAX_NUM_WORDS][MAX_WORD_LENGTH];
  int c = jb_read_line (fp,line,words);


  Attribute_name.clear();
  Counter_name.clear();
  Counter_type.clear();
  while (c != '\n') {
    if (words[1][0] == '\0' || words[2][0] != '\0') {
      // not two words in line 
      fprintf (stderr, "jb_read_header1: Error in Header section: wrong line length!\n");
      fprintf (stderr, "%s\n",line.c_str());
      exit(1);
    } else {
      if (strcmp(words[0],"attribute") == 0) {
	// ???  (copying from jb.py)
	if (strcmp(words[1],"region") != 0) {
	  Attribute_name.push_back(words[1]); 
	}
      } else if (strcmp(words[0],"basic") == 0 ||
		 strcmp(words[0],"papi") == 0 ||
		 strcmp(words[0],"user") == 0 ||
		 strcmp(words[0],"derived") == 0) {
	Counter_type.push_back(words[0]);
	Counter_name.push_back(words[1]);
      } else {
	fprintf (stderr,"jb_read_header1: Error in Header section: unknown keyword \n%s\n",
		 line.c_str());
      }
    }
    c = jb_read_line (fp,line,words);
  }
}

//========================================================================

void jb_print_globals (map<string,string> Globals)
{

  // Prints the global definitions


  Globals["lcaperf-version"] = global_version_write;
  for (map<string,string>::iterator pglobal=Globals.begin();
       pglobal != Globals.end();
       ++pglobal) {
    printf ("global %s %s\n",
	    (pglobal->first).c_str(),
	    (pglobal->second).c_str());
  }
  
  printf("\n");

  return;
}

//========================================================================

void jb_print_header (vector<string> & Attribute_name, 
		      vector<string> & Counter_name, 
		      vector<string> & Counter_type)
{
  if (global_version_write == "0.0") {
    jb_print_header0 (Attribute_name, Counter_name, Counter_type);
  } else if (global_version_write == "1.0") {
    jb_print_header1 (Attribute_name, Counter_name, Counter_type);
  } else {
    fprintf (stderr, "Error: Unknown version %s in jb_print_header!\n",
	     global_version_write);
      exit(1);
  }
}

//========================================================================

void jb_print_header0 (vector<string> & Attribute_name, 
		      vector<string> & Counter_name, 
		      vector<string> & Counter_type)
{
  printf ("region\n");
  int i;
  for (i=0; i<Attribute_name.size(); i++) {
    printf ("attribute %s\n",Attribute_name[i].c_str());
  }

  for (i=0; i<Counter_name.size(); i++) {
    printf ("%s %s\n",Counter_type[i].c_str(),Counter_name[i].c_str());
  }

  printf ("\n");
  
}

//========================================================================

void jb_print_header1 (vector<string> & Attribute_name, 
		      vector<string> & Counter_name, 
		      vector<string> & Counter_type)
{
  printf ("attribute region\n");

  int i;

  for (i=0; i<Attribute_name.size(); i++) {
    printf ("attribute %s\n",Attribute_name[i].c_str());
  }

  for (i=0; i<Counter_name.size(); i++) {
    printf ("%s %s\n",Counter_type[i].c_str(),Counter_name[i].c_str());
  }

  printf ("\n");
}

//========================================================================

bool jb_read_next_record (FILE * fp, 
			   int num_attributes, 
			   int num_counters, 
			   string & augregion,
			   vector<string> & record_counters)
{
  record_counters.clear();

  string line;
  char words[MAX_NUM_WORDS][MAX_WORD_LENGTH];
  int c = jb_read_line (fp,line,words);

  // Read region name

  if (c == EOF) return false;
  if (words[0][0] == '\0' || words[1][0] != '\0') {
    fprintf (stderr,"Error in Records section: expecting region name\n");
    fprintf (stderr,"%s\n",line.c_str());
    exit(1);
  }

  augregion = words[0];

  // Read Attributes

  int i;
  for (i=0; i<num_attributes; i++) {
    c = jb_read_line (fp,line,words);
    if (words[0][0] == '\0' || words[1][0] != '\0') {
      fprintf (stderr,"Error in Records section: expecting attribute value\n");
      fprintf (stderr,"%s\n",line.c_str());
      exit(1);
    } else {
      augregion = augregion + " " + words[0];
    }
  }

  // Read Counters

  record_counters.clear();
  for (i=0; i<num_counters; i++) {
    c = jb_read_line (fp,line,words);
    if (words[0][0] == '\0' || words[1][0] != '\0') {
      fprintf (stderr,"Error in Records section: expecting counter value\n");
      fprintf (stderr,"%s\n",line.c_str());
      exit(1);
    } else {
      record_counters.push_back (words[0]);
    }
  }

  c = jb_read_line (fp,line,words);

  if (c != '\n' && c != EOF) {
      fprintf (stderr,"Error in Records section: expecting newline\n");
      fprintf (stderr,"%s\n",line.c_str());
      exit(1);
  }

  return (true);
}

//========================================================================

void jb_read_records (FILE * fp, 
		      int num_attributes, 
		      int num_counters, 
		      map<string,vecstr> & Records)
{
  string augregion;
  vector<string> counters;

  while (jb_read_next_record(fp,num_attributes,num_counters,
			     augregion, counters)) {
    Records[augregion] = counters;
  }
}

//========================================================================

void jb_print_records (map<string,vecstr> & Records)
{
  map<string,vecstr>::iterator it_record;
  for (it_record = Records.begin();
       it_record != Records.end();
       ++it_record) {
    jb_print_record (it_record->first, it_record->second);
  }
}

//========================================================================

void jb_print_record (const string & augregion, const vector<string> & counters)
{

  // @@@ buffer size enough for about 100 attributes and counters (max 20 digits)
  const int maxstr = 2000;

  char buffer[maxstr+1];
  int c=0;
  vector<string> augregion_split;

  jb_split (augregion, augregion_split);
        
  // region and attributes

  int i;
  for (i=0; i<augregion_split.size(); i++) {
    //    printf ("%s\n",augregion_split[i].c_str());
    c += snprintf (buffer+c,maxstr-c,"%s\n",augregion_split[i].c_str());
    if (c>maxstr) {
      fprintf (stderr,"jb_print_record(): out of buffer space! : \n   %d\n   %d\n",
	       c, maxstr);
      exit(1);
    }
  }

  // counters

  for (i=0; i<counters.size(); i++) {
    //    printf ("%s\n",augregion_split[i].c_str());
    std::string s = counters[i];
    if (s == STAR_FLOAT_STRING) s = "*";
    if (s == HASH_FLOAT_STRING) s = "#";
    c += snprintf (buffer+c,maxstr-c,"%s\n",s.c_str());
    if (c>maxstr) {
      fprintf (stderr,"jb_print_record(): out of buffer space! : \n   %d\n   %d\n",
	       c, maxstr);
      exit(1);
    }
  }

  printf ("%s\n",buffer);

}

//========================================================================

void jb_split (string augregion, vector<string> & augregion_split)
{
  int i;

  augregion_split.clear();
  i = augregion.find(" ");
  while (i != string::npos) {
    augregion_split.push_back(augregion.substr(0,i));
    augregion.erase (0,i+1);
    i = augregion.find(" ");
  }
  augregion_split.push_back(augregion);
}

//========================================================================

void jb_counters_merge (vector<string> counters1, 
		   vector<string> counters2, 
		   vector<string> & new_counters)
{
  if (counters1.size() != counters2.size()) {
    fprintf (stderr,"jb_counters_merge(): length mismatch! : \n   %d\n   %d\n",
	     int(counters1.size()), int(counters2.size()));
    exit(1);
  }

  new_counters.clear();

  string new_counter;
  for (int i=0; i<counters1.size(); i++) {
    // Counter is already a conflict: keep it a conflict
    if (counters1[i]=="#" || counters2[i]=="#") {
      new_counter = "#";
      // If one value is undefined, use the other value
      // (which may also be undefined)
    } else if (counters1[i]=="*") {
      new_counter = counters2[i];
    } else if (counters2[i]=="*") {
      new_counter = counters1[i];
      // If counters don't match, it's a new conflict
    } else if (counters1[i] != counters2[i]) {
      new_counter = "#";
      // If counters do match, use that value
    } else if (counters1[i] == counters2[i]) {
      new_counter = counters1[i];
    }
    new_counters.push_back(new_counter);
  }
  return;
}

//------------------------------------------------------------------------

void jb_augregion_rotate_ (string augregion, int index, string & augregion_rotated) {
  // Note: index is relative to augregion, not attribute
  
  vector<string> augregion_split;
  jb_split (augregion, augregion_split);
  int last = augregion_split.size() - 1;
  string t = augregion_split[index];
  augregion_split[index] = augregion_split[last];
  augregion_split[last] = t;
  augregion_rotated = augregion_split[0];
  for (int i=1; i<augregion_split.size(); i++) {
    augregion_rotated = augregion_rotated + " " + augregion_split[i];
  }
}
    
//------------------------------------------------------------------------

string jb_augregion_root_ (string augregion)
{
  return augregion.substr(0,augregion.rfind(" "));
}
