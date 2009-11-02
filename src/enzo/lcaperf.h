#ifndef LCAPERF_H
#define LCAPERF_H

//======================================================================
//
// File:        lcaperf.h
//
// Description: Parallel performance monitoring class
//
//----------------------------------------------------------------------
//
// Classes:     LcaPerf
//
//----------------------------------------------------------------------
//
// Copyright 2004-2005 James Bordner
// Copyright 2004-2005 Laboratory for Computational Astrophysics
// Copyright 2004-2005 Regents of the University of California
//
//======================================================================

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include <map>
#include <string>
#include <set>
#include <vector>
#include <stack>

#include "lcaperf.def"

#ifdef USE_PAPI
#  include "papi.h"
#else
#  define PAPI_NULL 0
#endif

// *********************************************************************
// <<< OLD

typedef std::vector<long long> OldCountersType;
typedef std::map<std::string,OldCountersType> RegionCountersType;
typedef std::stack<std::string> FrameType;
typedef std::map<std::string,int> StateType;

// OLD >>>
// *********************************************************************


// *********************************************************************
// <<< NEW

typedef std::map<std::string,class Counter *> CountersType;

// NEW >>>
// *********************************************************************

extern long long lcampi_global_calls_;
extern long long lcampi_global_recv_bytes_;
extern long long lcampi_global_send_bytes_;
extern long long lcampi_global_synch_;
extern long long lcampi_recv_bytes_;
extern long long lcampi_recv_calls_;
extern long long lcampi_send_bytes_;
extern long long lcampi_send_calls_;
extern long long lcampi_time_;

//======================================================================
// LcaPerf base class
//======================================================================

class LcaPerf {

  //--------------------------------------------------------------
  // Public member functions
  //--------------------------------------------------------------

 public:

  LcaPerf();
  ~LcaPerf();

  // Initialization region functions

  void new_attribute (const char * attribute_name, int attribute_type);
  void new_counter  (const char * counter_name, int counter_type);

  void delete_attribute (const char * attribute_name);
  void delete_counter (const char * counter_name, int counter_type);

  void begin (const char * suffix);
  void end   (const char * suffix);

  // Global variables

  void global    (const char * global_name,
		  long long    value);

  // Active region functions using names

  void start     (const char * region_name);

  void stop      (const char * region_name);      

  void increment (const char * counter_name,   
		  long long    value);

  void assign    (const char * counter_name,   
		  long long    value);

  void attribute (const char * attribute_name, 
		  void       * value_pointer, 
		  int          value_type);

  //  void flush ();

  //--------------------------------------------------------------
  // Public member functions -- TO ADD
  //--------------------------------------------------------------

  long long counter (const char * counter_name, int counter_type);

  //--------------------------------------------------------------
  // Public member functions -- TO DELETE
  //--------------------------------------------------------------

  //  void new_category (const char * category_name, int category_type);
  //  void category (const char * category_name, 
  //		 void       * value_pointer, 
  //		 int          value_type);

  //--------------------------------------------------------------
  // Private member functions
  //--------------------------------------------------------------

 private:

  // Counters

  void create_directories_ () ;
  inline long long * create_counters_ () const;
  inline void delete_counters_ (long long *) const;
  void start_counters_ (long long *);
  void stop_counters_ (long long *);
  void update_counters_global_ (std::string,long long *);
  void update_counters_parent_ (long long *,long long *);

  //    basic counters

  inline long long get_real_ () const; // Return real (wallclock) time
#ifdef USE_PAPI
  inline long long get_virt_ () const; // Return virt (CPU) time
#endif

  //    papi counters

  std::string papi_name_ (std::string); // Return e.g. PAPI_FP_INS given fp-ins

  void new_papi_counter_  (const char * counter_name, int counter_type);
  void delete_papi_counter_ (const char * attribute_name, int counter_type);
  inline void read_papi_ (); // Read PAPI event values into papi_counter_values_
  void start_papi_counters_ (); // Starts PAPI event set
  void stop_papi_counters_ ();  // Stops PAPI event set

  void clear_papi_eventset_();
  void create_papi_eventset_();
  void delete_papi_eventset_();

  int insert_papi_event_(const char * counter_name);
  int delete_papi_event_(const char * counter_name);

  //    user counters

  void new_user_counter_  (const char * counter_name, int counter_type);
  void create_user_counters_ ();
  void delete_user_counter_ (const char * attribute_name, int counter_type);

  // Regions

  std::string generate_augmented_region_ (std::string) const;
  std::string merge_augmented_regions_ (std::string, std::string) const;
  inline void delete_global_regions_ ();

  // Counter Indices

  void adjust_indices_ (int counter_type, int number);
  void clear_indices_ ();
  void check_indices_() const;

  // Files

  void write_ (std::string suffix = "");
  void write_globals_ (FILE *) const;                          // Write file globals
  void write_header_ (FILE *) const;                           // Write file header
  void write_field_ (FILE *, std::string, long long *) const;  // Write file field
  void pack_field_ (long long *, std::string, long long *) const;  
                                                // Pack long long array with data

  void error_ (const char * message, const char * file_name, int line_number);
  long long region_index_ (const char * region_name); 
                                                   // return index for region name
      

#ifdef OVERHEAD
  inline void overhead_init_() ;
  inline void overhead_start_(int i) ;
  inline void overhead_stop_(int i) ;
  inline void overhead_write_() const;
  inline void overhead_summary_() const;
#endif

  //--------------------------------------------------------------
  // Private member data
  //--------------------------------------------------------------

 private:

  // Counters

  long long rtime_begin_;    // Real time at start of lcaperf
#ifdef USE_PAPI
  long long vtime_begin_;    // CPU time at start of lcaperf
#endif

  int num_basic_counters_;                            // Number of defined papi counters
  std::map<std::string, int> index_;                  // index into counters arrays
  std::stack<long long *>           counters_stack_;  // stack of local counters
  std::map<std::string,long long *> counters_global_; // All global counters

  int ifirst_basic_;   // indices of first and last basic counters in Counters
  int ilast_basic_;
  int ifirst_papi_;    // indices of first and last PAPI counters in Counters
  int ilast_papi_;
  int ifirst_user_;    // indices of first and last user counters in Counters
  int ilast_user_;

  // User counters

  int num_user_counters_;                         // Number of defined user counters
  std::vector<long long> user_counter_values_;    // Array of user counters
  std::vector<int> user_counter_types_;           // Array of user counter types 
  //                                                (absolute or relative (interval))
  std::vector<std::string> user_counter_names_;    // User counter name given index
  std::map<std::string,int> user_counter_indices_; // User counter index given name

  // PAPI counters

  bool papi_is_active_;
  int papi_event_set_;                               // handle to PAPI event set
  int num_papi_counters_;                            // Number of defined papi counters
  std::vector<long long>    papi_counter_values_;    // Array of papi counters
  std::vector<std::string>  papi_counter_names_;     // Papi counter name given index
  std::map<std::string,int> papi_counter_indices_;   // Papi counter index given name
  std::vector<bool>         papi_counter_is_active_; // Whether individual PAPI counters are active

  // Attributes

  int num_attributes_;                         // Number of attributes defined
  std::map<std::string,int> attribute_indices_; // Attribute index given name
  std::vector<int>          attribute_types_;   // Attribute types given index
  std::vector<std::string>  attribute_names_;   // Attribute name given index
  std::vector<std::string>  attribute_values_;  // Attribute values given index

  // Globals

  std::map<std::string, long long> global_values_; // User-defined global variables

  // Files

  bool lDirCreated_;   // Whether the LCAPERF subdirectory has been created yet
  std::string suffix_; // Current file suffix
  std::map<std::string,bool> suffix_list_;     // List of all suffixes so far
  std::map<std::string,FILE *> fp_;            // Mapping of suffix to file pointer
  std::map<std::string,std::string> filename_; // Mapping of suffix to file name

#ifdef OVERHEAD
  long long overhead_initial_time_; // initial time
  long long overhead_total_;        // total overhead
  long long overhead_[20];          // overhead of individual functions
  int call_depth_;  // depth of LcaPerf function calls for overhead_total_
#endif

  int np_;  // Processor count
  int ip_;  // Processor rank

  bool isActive_;            // Whether we're between begin() and end()

  char lcaperf_dir_[12];       // LCAPERF.###


  // **************************************************************************

  char path_ [LCAPERF_MAX_PATH_LENGTH]; // Path to LCAPERF directory
  char path_prev_ [LCAPERF_MAX_PATH_LENGTH]; // Saved path from application
  FrameType frame_;         // Stack frame--required for exclusive counts

};

// *********************************************************************

inline long long * LcaPerf::create_counters_ () const
{
  long long * counters = new long long [MAX_NUM_COUNTERS];
  for (int i=0; i<MAX_NUM_COUNTERS; i++) counters[i] = 0;
  return counters;
}

//----------------------------------------------------------------------
inline void LcaPerf::delete_counters_ (long long * counters) const
{
  delete [] counters;
  counters = 0;
}

//----------------------------------------------------------------------

inline void LcaPerf::delete_global_regions_ ()
{
  counters_global_.clear();
}

//----------------------------------------------------------------------

inline long long LcaPerf::get_real_ () const
{
  struct timeval tv;
  struct timezone tz;
  gettimeofday(&tv,&tz);
  return (long long) (1000000) * tv.tv_sec + tv.tv_usec;
}

//----------------------------------------------------------------------

#ifdef USE_PAPI
inline long long LcaPerf::get_virt_ () const
{
  //  return PAPI_get_virt_cyc(); // More precise than PAPI_get_virt_usec()
  return PAPI_get_virt_usec(); // Doesn't roll-over like PAPI_get_virt_cyc()
}
#endif

//----------------------------------------------------------------------

inline void LcaPerf::read_papi_ ()
{
#ifdef USE_PAPI
  if (papi_event_set_ != PAPI_NULL) {
    int retval;
    for (int k=0; k<num_papi_counters_; k++) papi_counter_values_[k]=0;
    _CALL_PAPI(PAPI_read(papi_event_set_,&papi_counter_values_[0]),"read",PAPI_OK,retval);
  }
#endif
}

//========================================================================
// OVERHEAD ROUTINES
//========================================================================
//========================================================================

#ifdef OVERHEAD
inline void LcaPerf::overhead_init_()
{
  call_depth_ = 0;
  overhead_total_ = 0;
  for (int i=0; i<20; i++) {
    overhead_[i] = 0;
  }
  overhead_initial_time_ = get_real_();
}
#endif

//------------------------------------------------------------------------

#ifdef OVERHEAD
inline void LcaPerf::overhead_start_(int index)
{
  long long current_time = get_real_ ();
  if (call_depth_++ == 0) {
    overhead_total_ = current_time - overhead_total_;
  }
  overhead_[index] = current_time - overhead_[index];
}
#endif

//------------------------------------------------------------------------

#ifdef OVERHEAD
inline void LcaPerf::overhead_stop_(int index)
{
  long long current_time = get_real_ ();
  if (--call_depth_ == 0) {
    overhead_total_ = current_time - overhead_total_;
  }
  overhead_[index] = current_time - overhead_[index];
}
#endif

//------------------------------------------------------------------------

#ifdef OVERHEAD
inline void LcaPerf::overhead_write_() const
{
  const char * names[] = {"LcaPerf()",           // 0
			  "~LcaPerf()",          // 1
			  "new_attribute()",    // 2
			  "new_counter()",      // 3
			  "delete_attribute()", // 4
			  "delete_counter()",   // 5
			  "begin()",            // 6
			  "end()",              // 7
			  "start()",            // 8
			  "stop()",             // 9
			  "increment()",        // 10
			  "assign()",           // 11
			  "attribute()",        // 12
			  "counter()",          // 13
			  "global()",           // 14
			  0};
  
  if (ip_==0) {
    long long total_time = get_real_ () - overhead_initial_time_;
    for (int i=2; names[i]; i++) {
      printf ("lcaperf: %18s: %10lld %5.2f%%\n", names[i],
	      overhead_[i],100.0*(double)(overhead_[i])/total_time);
    }
    printf ("lcaperf: %18s: %10lld %5.2f%%\n","Total",
	    overhead_total_,100.0*(double)(overhead_total_)/total_time);
  }
}
#endif

//------------------------------------------------------------------------

#ifdef OVERHEAD
inline void LcaPerf::overhead_summary_() const
{
  if (ip_==0) {
    long long total_time = get_real_ () - overhead_initial_time_;
    printf ("lcaperf: %18s: %10lld %5.2f%%\n","Total",
	    overhead_total_,100.0*(double)(overhead_total_)/total_time);
  }
}
#endif

// *********************************************************************

extern LcaPerf lcaperf;

// *********************************************************************
// Fortran bindings
// *********************************************************************


//----------------------------------------------------------------------
extern "C" void LCAPERF_FORTRAN(lca_new_attribute) 
  (char *     attribute_name,
   int        *pattribute_type);

//----------------------------------------------------------------------
extern "C" void LCAPERF_FORTRAN(lca_new_counter) 
  (char *     counter_name, 
   int        *pcounter_type,
   int        *perror_code);

//----------------------------------------------------------------------
extern "C" void LCAPERF_FORTRAN(lca_begin) 
  (char *     suffix);  /* WARNING--function has default value */

//----------------------------------------------------------------------
extern "C" void LCAPERF_FORTRAN(lca_end) 
  ();

//----------------------------------------------------------------------
extern "C" void LCAPERF_FORTRAN(lca_start) 
  (char *     region_name);

//----------------------------------------------------------------------
extern "C" void LCAPERF_FORTRAN(lca_stop) 
  (char *     region_name);

//----------------------------------------------------------------------
extern "C" void LCAPERF_FORTRAN(lca_attribute) 
  (char *     attribute_name, 
   void       *pvalue_pointer,  /* WARNING--passing pointer as integer */
   int        *pvalue_type);

//----------------------------------------------------------------------
extern "C" void LCAPERF_FORTRAN(lca_increment) 
  (char *      counter_name,   
   long long  *pvalue);
//----------------------------------------------------------------------
extern "C" void LCAPERF_FORTRAN(lca_assign) 
  (char *      counter_name,   
   long long  *pvalue);


#endif
