//======================================================================
//
// File:        lcaperf.C
//
// Description: Parallel performance monitoring class
//
//----------------------------------------------------------------------
//
// Copyright 2004-2005 James Bordner
// Copyright 2004-2005 Laboratory for Computational Astrophysics
// Copyright 2004-2005 Regents of the University of California
//
//======================================================================

// #define TRACE_LCAPERF
// #define TRACE_STATE
// #define TRACE_PAPI
// #define TRACE_FRAME
// #define TRACE_COUNTERS_STACK
// #define TRACE_LCAMEM

#ifdef TRACE_LCAMEM
#  define _LCAMEM_TRACE_ \
  if (ip_==0) printf ("DEBUG %s:%d %lld %lld %lld %lld %lld %lld\n", \
		      __FILE__,__LINE__, \
		      lcamem::bytes_, \
		      lcamem::bytesHigh_, \
		      lcamem::deleteBytes_, \
		      lcamem::deleteCalls_, \
		      lcamem::newBytes_, \
		      lcamem::newCalls_);
#else
#  define _LCAMEM_TRACE_
#endif

//----------------------------------------------------------------------
// Debugging level
//----------------------------------------------------------------------

#define OUTPUT_NOCHECK 0
#define OUTPUT_ERROR   1
#define OUTPUT_WARNING 2
#define OUTPUT_DEBUG   3

#define OUTPUT_LEVEL OUTPUT_WARNING

//----------------------------------------------------------------------

#define CHECK_PAPI
// #define DEBUG

#include <assert.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <errno.h>
#include <string.h>
#include <unistd.h>

#ifdef USE_MPI
#  include <mpi.h>
#endif

#include "lcaperf.h"
#include "lcamem.h"


// *********************************************************************
// namespace lca {
// *********************************************************************

//======================================================================
// PUBLIC MEMBER FUNCTIONS
//======================================================================

//----------------------------------------------------------------------
LcaPerf::LcaPerf () // + Initialize the LcaPerf object
//----------------------------------------------------------------------
    // Counters
  : rtime_begin_(),
#ifdef USE_PAPI
    vtime_begin_(),
#endif
    num_basic_counters_(0),
    index_(),
    counters_stack_(),
    counters_global_(),
    ifirst_basic_(0),
    ilast_basic_(0),
    ifirst_papi_(0),
    ilast_papi_(0),
    ifirst_user_(0),
    ilast_user_(0),
    // User counters
    num_user_counters_(0),
    user_counter_values_(),
    user_counter_types_(),
    user_counter_names_(),
    user_counter_indices_(),
    // Papi counters
    papi_is_active_(false),
    papi_event_set_(PAPI_NULL),
    num_papi_counters_(0),
    papi_counter_values_(),
    papi_counter_names_(),
    papi_counter_indices_(),
    papi_counter_is_active_(),
    // Attributes
    num_attributes_(0),
    attribute_indices_(),
    attribute_types_(),
    attribute_names_(),
    attribute_values_(),
    // Files
    lDirCreated_(false),
    suffix_(),
    suffix_list_(),
    fp_(),
    filename_(),
#ifdef OVERHEAD
    overhead_initial_time_(0),
    overhead_total_(0),
    overhead_(),
    call_depth_(0),
#endif
    // OTHER
    isActive_(false),
    np_(1),
    ip_(0),
    path_(),
    path_prev_(),
    frame_()
#ifdef UNDER_CONSTRUCTION
  , counters_()
#endif
{
  _TRACE_LCAPERF("LcaPerf");
  _OVERHEAD_INIT;

  // Initialize counter indices

  clear_indices_();

  // Initialize PAPI

  int retval;
  _CALL_PAPI(PAPI_library_init(PAPI_VER_CURRENT),"library_init",
	     PAPI_VER_CURRENT,retval);
  _CALL_PAPI(PAPI_set_debug(PAPI_VERB_ECONT),"set_debug",PAPI_OK,retval);
  if (retval) {
    fprintf (stderr, _ERROR 
	     "PAPI_library_init() returned 'retval=%d'\n",
	     __FILE__,__LINE__,retval);
    fflush(stderr);
  }

  rtime_begin_ = get_real_();
#ifdef USE_PAPI
  vtime_begin_ = get_virt_();
#endif

  //  _CALL_PAPI(PAPI_thread_init((unsigned long (*)(void))(pthread_self), 0),
  //             "thread_init",PAPI_OK,retval),

#ifdef LCAMEM
  printf ("lcamem ARRAYS_ONLY = %s\n", 
	  lcamem_define_arrays_only ? "True" : "False");
#endif

  _TRACE_LCAPERF("LcaPerf");

}

//----------------------------------------------------------------------
LcaPerf::~LcaPerf () // + Finalize the LcaPerf object
//----------------------------------------------------------------------
{
  _TRACE_LCAPERF("~LcaPerf");

  // Close files

  std::map<std::string,FILE *>::iterator pfp;
  for (pfp  = fp_.begin(); pfp != fp_.end(); ++pfp) {
    fclose ((*pfp).second);
  }

  // Clean up PAPI

  delete_papi_eventset_();

  _OVERHEAD_WRITE;
  _TRACE_LCAPERF("~LcaPerf");
}

//----------------------------------------------------------------------
void LcaPerf::new_attribute  // + Declare a new attribute
   (const char * attribute_name, 
    int          attribute_type)
//----------------------------------------------------------------------
{
  _OVERHEAD_START(2); // new_attribute
  _TRACE_LCAPERF1("new_attribute",attribute_name);
  if (attribute_indices_.find(attribute_name) != attribute_indices_.end()) {
    fprintf (stderr, _ERROR 
	     "new_attribute() called for pre-existing attribute %s\n",
	     __FILE__,__LINE__,attribute_name);
    fflush(stderr);
  }

  int index = num_attributes_;

  ++num_attributes_;

  attribute_types_.resize(num_attributes_);
  attribute_values_.resize(num_attributes_);
  attribute_names_.resize(num_attributes_);

  attribute_types_[index] = attribute_type;
  attribute_values_[index] = LCAPERF_NULL_STRING;
  attribute_names_[index] = attribute_name;
  attribute_indices_[attribute_name] = index;
  
  _TRACE_LCAPERF1("new_attribute",attribute_name);
  _OVERHEAD_STOP(2); // new_attribute
}

//----------------------------------------------------------------------
void LcaPerf::new_counter (const char * counter_name, int counter_type)
//----------------------------------------------------------------------
{
  _OVERHEAD_START(3); // new_counter
  _TRACE_LCAPERF1("new_counter",counter_name);

  switch (counter_type) {

  case LCAPERF_COUNTER_TYPE_PAPI:

    new_papi_counter_ (counter_name, counter_type);
    break;

  case LCAPERF_COUNTER_TYPE_USER_REL:
  case LCAPERF_COUNTER_TYPE_USER_ABS:

    new_user_counter_ (counter_name, counter_type);
    break;

  default:

    fprintf (stderr, _ERROR "Bad counter_type in new_counter(%s,%d)!\n",
	     __FILE__,__LINE__,counter_name, counter_type);
    break;

  }

  _TRACE_LCAPERF1("new_counter",counter_name);
  _OVERHEAD_STOP(3); // new_counter
}

//----------------------------------------------------------------------
void LcaPerf::delete_attribute (const char * attribute_name)
//----------------------------------------------------------------------
{
  _OVERHEAD_START(4); // delete_attribute
  _TRACE_LCAPERF1("delete_attribute",attribute_name);

  // Error checking

  if (isActive_) {
    fprintf (stderr, _ERROR 
	     "delete_attribute(%s) called while in active mode!\n",
	     __FILE__,__LINE__,attribute_name);
    fflush(stderr);
    return;
  }

  if (strcmp(attribute_name,"*") == 0) {

    // Delete all attributes

    num_attributes_ = 0;
    attribute_indices_.clear();
    attribute_types_.clear();
    attribute_names_.clear();
    attribute_values_.clear();

  } else {

    if (attribute_indices_.find(attribute_name) != 
	attribute_indices_.end()) {

      // Delete specified attribute if it exists

      int attribute_index = attribute_indices_[attribute_name];
      std::map<std::string,int>::iterator ci;
      -- num_attributes_;
      attribute_indices_.erase(attribute_indices_.find(attribute_name));
      attribute_types_.erase  (attribute_types_.begin() +attribute_index);
      attribute_names_.erase  (attribute_names_.begin() +attribute_index);
      attribute_values_.erase (attribute_values_.begin()+attribute_index);

      // and update indicies to fill in the new vacancy

      for (ci = attribute_indices_.begin();
	   ci != attribute_indices_.end();
	   ++ci) {
	if (ci->second > attribute_index) {
	  --ci->second;
	}
      }

    } else {

      // Print error message if attribute does not exist

      fprintf (stderr, _ERROR 
	       "delete_attribute(%s) called with non-existing attribute!\n",
	       __FILE__,__LINE__,attribute_name);
      fflush(stderr);
    }
  }
  _TRACE_LCAPERF1("delete_attribute",attribute_name);
  _OVERHEAD_STOP(4); // delete_attribute
}

//----------------------------------------------------------------------
void LcaPerf::delete_counter (const char * counter_name, int counter_type)
//----------------------------------------------------------------------
{
  _OVERHEAD_START(5);  // delete_counter
  _TRACE_LCAPERF1("delete_counter",counter_name);

  // Error checking

  if (isActive_) {
    fprintf (stderr, _ERROR 
	     "delete_counter(%s,%d) called while in active mode!\n",
	     __FILE__,__LINE__,counter_name,counter_type);
    fflush(stderr);
    return;
  }

  bool deleteAll = (strcmp(counter_name,"*") == 0);

  switch (counter_type) {

  case LCAPERF_COUNTER_TYPE_PAPI:

    delete_papi_counter_(deleteAll ? "\0" : counter_name, counter_type);
    break;

  case LCAPERF_COUNTER_TYPE_USER_REL:
  case LCAPERF_COUNTER_TYPE_USER_ABS:

    delete_user_counter_(deleteAll ? "\0" : counter_name, counter_type);
    break;

  default:

    fprintf (stderr, _ERROR 
	     "Bad counter_type %d in delete_counter(%s,%d)!\n",
	     __FILE__,__LINE__,counter_type,counter_name,counter_type);
    fflush(stderr);
  }
    
  _TRACE_LCAPERF1("delete_counter",counter_name);
  _OVERHEAD_STOP(5);  // delete_counter
}

//----------------------------------------------------------------------
void LcaPerf::begin (const char * suffix)
//----------------------------------------------------------------------
{

  _OVERHEAD_START(6);  // begin
  _TRACE_LCAPERF1("begin", suffix);

  if (!isActive_) {

    // Change state to Active mode

    isActive_ = true;

    // Store suffix for later

    suffix_ = suffix;
    suffix_list_[suffix_] = true;

    // Get processor ID if needed
    // (not in lcaperf() because MPI will not have been initialized yet

#ifdef USE_MPI
    MPI_Comm_size (MPI_COMM_WORLD, &np_);
    MPI_Comm_rank (MPI_COMM_WORLD, &ip_);
#else
    np_ = 1;
    ip_ = 0;
#endif

    // Create directory structure if needed (root processor only!)

    if (! lDirCreated_) create_directories_ ();

    // Start counters

    start_papi_counters_();

    // Start a top-level "dummy" region: simplifies updating exclusize counts

    start (LCAPERF_NULL_STRING);
    
  } else {
    fprintf (stderr, _ERROR "begin(%s) called while in active mode!\n",
	     __FILE__,__LINE__,suffix);
    fflush(stderr);

  }

  _TRACE_LCAPERF1("begin", suffix);
  _OVERHEAD_STOP(6);  // begin
}

//----------------------------------------------------------------------
void LcaPerf::end (const char * suffix)
//----------------------------------------------------------------------
{
  _OVERHEAD_START(7); // end
  _TRACE_LCAPERF1("end",suffix);

  // Stop the top blank region

  stop (LCAPERF_NULL_STRING);
  if (!frame_.empty()) {
    fprintf (stderr,_ERROR "function stack not empty\n",
             __FILE__,__LINE__);
    fflush(stderr);
  }

  if (isActive_) {

    isActive_ = false;

    // Write events 
    // @@@ Ideally this could be turned off, if e.g. an optional parameter
    // @@@ were included (LcaPerf::end(suffix,flush=1)).  But
    // @@@ optional parameters are difficult with Fortran interfaces.

    write_ ();

    // Stop PAPI eventset for delete_counter() and new_counter() with PAPI counters

    stop_papi_counters_ ();

    // Delete all regions

    delete_global_regions_();

  } else {
    fprintf (stderr,_ERROR "end() called after event counting stopped\n",
	     __FILE__,__LINE__);
    fflush(stderr);
  }
  _TRACE_LCAPERF1("end",suffix);
  _OVERHEAD_STOP(7)  // end
  _OVERHEAD_SUMMARY;
}

//----------------------------------------------------------------------
void LcaPerf::global (const char * global_name,
		     long long global_value)
//----------------------------------------------------------------------
{
  _OVERHEAD_START(14); // global
  _TRACE_LCAPERF1("global",global_name);

  // Generate Augmented Region name given Region

  global_values_[global_name] = global_value;

  _TRACE_LCAPERF1("global",global_name);
  _OVERHEAD_STOP(14);  // global
}

//----------------------------------------------------------------------
void LcaPerf::start (const char * region_base)
//----------------------------------------------------------------------
{
  _OVERHEAD_START(8); // start
  _TRACE_LCAPERF1("start",region_base);

  // Make sure we are in active state

  if (!isActive_) {
    fprintf (stderr,_ERROR "start(%s) called before begin()\n",
             __FILE__,__LINE__,region_base);
    fflush(stderr);
  }
  assert (isActive_);

  // Generate Augmented Region name given Region

  std::string region = generate_augmented_region_(region_base);

  // Save region on stack

  _TRACE_STACK_PUSH(region);
  frame_.push(region);

  // Create local counters

  long long * counters = create_counters_();

  // Start local counters

  start_counters_ (counters);

  // Push local counters on stack

  if (counters == 0) printf (_DEBUG "Pushing 0!!!\n",__FILE__,__LINE__);
  _TRACE_COUNTERS_STACK_PUSH(region,counters);
  counters_stack_.push(counters);

  _TRACE_LCAPERF1("start",region_base);
  _OVERHEAD_STOP(8);  // start
}

//----------------------------------------------------------------------
void LcaPerf::stop (const char * region_base)
//----------------------------------------------------------------------
{
  _OVERHEAD_START(9); // stop
  _TRACE_LCAPERF1("stop",region_base);

  // Make sure we are in active state

  if (!isActive_) {
    fprintf (stderr,_ERROR "stop(%s) called before begin()\n",
             __FILE__,__LINE__,region_base);
    fflush(stderr);
  }

  // Pop local counters from stack

  long long * counters = counters_stack_.top();
  _TRACE_COUNTERS_STACK_POP;
  counters_stack_.pop();

  // Stop local counters

  stop_counters_ (counters);

  // Update parent counters

  if (counters_stack_.size() > 0) {
    update_counters_parent_ (counters_stack_.top(), counters);
  }

  // Generate Augmented Region name given Region
  
  std::string region       = generate_augmented_region_(region_base);

  // Check that attributes haven't changed, and merge if they have

  std::string region_start = frame_.top();
  frame_.pop();

  if (region_start != region) {
    region = merge_augmented_regions_(region_start,region);
  }

  update_counters_global_ (region, counters);

  delete_counters_ (counters);

  _TRACE_LCAPERF1("stop",region_base);
  _OVERHEAD_STOP(9) // stop
}

//----------------------------------------------------------------------
void LcaPerf::increment (const char *counter_name, 
			long long value)
//----------------------------------------------------------------------
{
  _OVERHEAD_START(10); // increment
  _TRACE_LCAPERF1("increment",counter_name);

  int index = user_counter_indices_[counter_name];

#if (OUTPUT_LEVEL >= OUTPUT_ERROR)

  // Check for incrementing non-existent counter

  if (user_counter_indices_.find(counter_name) 
      == user_counter_indices_.end()) {
    printf (_ERROR "incrementing nonexistent counter %s!\n",
	    __FILE__,__LINE__,counter_name);
  }

#endif

#if (OUTPUT_LEVEL >= OUTPUT_WARNING)

  // Check for calling increment() with absolute counters

  if (user_counter_types_[index] != LCAPERF_COUNTER_TYPE_USER_REL) {
    printf (_ERROR "incrementing counter of type %d != %d!\n",
	    __FILE__,__LINE__,user_counter_types_[index],LCAPERF_COUNTER_TYPE_USER_REL);
  }

#endif

  _LCAMEM_TRACE_;
  user_counter_values_[index] += value;
  _LCAMEM_TRACE_;

  _TRACE_LCAPERF1("increment",counter_name);
  _OVERHEAD_STOP(10); // increment
}

//----------------------------------------------------------------------
void LcaPerf::assign (const char *counter_name, 
		     long long value)
//----------------------------------------------------------------------
{
  _OVERHEAD_START(11); // assign
  _TRACE_LCAPERF1("assign",counter_name);

    int index = user_counter_indices_[counter_name];

#if (OUTPUT_LEVEL >= OUTPUT_ERROR)
  
  // Check for assigning to nonexistent counter

  if (user_counter_indices_.find(counter_name) 
      == user_counter_indices_.end()) {
    printf (_ERROR "trying to access nonexistent counter %s!\n",
	    __FILE__,__LINE__,counter_name);
  }

#endif

#if (OUTPUT_LEVEL >= OUTPUT_WARNING)

    // Check for calling assign() with relative counters

    if (user_counter_types_[index] != LCAPERF_COUNTER_TYPE_USER_ABS) {
      printf (_ERROR "assigning to counter %s of type %d != %d!\n",
	      __FILE__,__LINE__,
	      counter_name, user_counter_types_[index],LCAPERF_COUNTER_TYPE_USER_ABS);
    }

#endif

    user_counter_values_[index] = value;

  _TRACE_LCAPERF1("assign",counter_name);
  _OVERHEAD_STOP(11); // assign
}

//----------------------------------------------------------------------
void LcaPerf::attribute (const char * attribute_name, 
			void       * attribute_value, 
			int          attribute_type)
//----------------------------------------------------------------------
{
  _OVERHEAD_START(12); // attribute
  _TRACE_LCAPERF1("attribute",attribute_name);

  if (attribute_indices_.find(attribute_name) == attribute_indices_.end()) {
    fprintf (stderr,_ERROR "attribute(%s) called for non-existing attribute!\n",
	     __FILE__,__LINE__,attribute_name);
    fflush(stderr);
  } else {

    char attribute_string[LCAPERF_MAX_STRING];
    std::string attribute_format;
    switch (attribute_type) {
    case LCAPERF_STRING: 
      sprintf (attribute_string,"%s",(char *) attribute_value);
      break;
    case LCAPERF_INT:
      sprintf (attribute_string,"%d",*((int *)attribute_value));
      break;
    case LCAPERF_FLOAT:
      sprintf (attribute_string,"%f",*((float *)attribute_value));
      break;
    case LCAPERF_DOUBLE:
      sprintf (attribute_string,"%lf",*((double *)attribute_value));
      break;
    case LCAPERF_NULL:
      strcpy (attribute_string,LCAPERF_NULL_STRING);
      break;
    default:
      fprintf (stderr,_ERROR "Unknown attribute_type %d in attribute()\n",
	       __FILE__,__LINE__,attribute_type);
      fflush(stderr);
      break;
    }
    
    attribute_values_[attribute_indices_[attribute_name]] = attribute_string;
  }
  _TRACE_LCAPERF1("attribute",attribute_name);
  _OVERHEAD_STOP(12);  // attribute
}

//----------------------------------------------------------------------
long long LcaPerf::counter (const char * counter_name, 
			   int counter_type)
//----------------------------------------------------------------------
{
  _OVERHEAD_START(13); // counter

  int index = user_counter_indices_[counter_name];

#if (OUTPUT_LEVEL >= OUTPUT_ERROR)
  
  // Check for assigning to nonexistent counter

  if (user_counter_indices_.find(counter_name) 
      == user_counter_indices_.end()) {
    printf (_ERROR "trying to access nonexistent counter %s!\n",
	    __FILE__,__LINE__,counter_name);
  }

#endif

  fprintf (stderr,"LcaPerf::counter() is not implemented yet!\n");
  fflush(stderr);
  
  _OVERHEAD_STOP(13); // counter
  return 0;
}

//======================================================================
// PRIVATE MEMBER FUNCTIONS
//======================================================================

//----------------------------------------------------------------------
void LcaPerf::error_ (const char * message, const char * file_name, int line_number)
//----------------------------------------------------------------------
{
  fprintf (stderr,_ERROR "%s\n",file_name,line_number,message);
  fflush(stderr);
}

//----------------------------------------------------------------------
void LcaPerf::adjust_indices_ (int counter_type, int number)
//----------------------------------------------------------------------
{
  //
  // papi
  // user
  //
  if (counter_type == LCAPERF_COUNTER_TYPE_PAPI) {

    ilast_papi_   += 2*number;
    ifirst_user_  += 2*number;
    ilast_user_   += 2*number;
    num_papi_counters_ += number;

  } else if (counter_type == LCAPERF_COUNTER_TYPE_USER_REL ||
	     counter_type == LCAPERF_COUNTER_TYPE_USER_ABS) {

    ilast_user_   += 2*number;
    num_user_counters_ += number;

  } else {
      fprintf (stderr, _ERROR 
	       "adjust_indices(%d,%d) bad counter type!\n",
	       __FILE__,__LINE__,counter_type,number);
      fflush(stderr);
  }
}

//----------------------------------------------------------------------
void LcaPerf::clear_indices_ ()
//----------------------------------------------------------------------
{
  // call-count
  // time-real-start
  // time-real-stop
  // time-real-incl
  // time-real-excl
  num_basic_counters_ = 5;

#ifdef USE_PAPI
  // time-virt-start
  // time-virt-stop
  // time-virt-incl
  // time-virt-excl
  num_basic_counters_ += 4;
#endif

#ifdef USE_PMPI
  // MPI-global-calls
  // MPI-global-recv-bytes
  // MPI-global-send-bytes
  // MPI-global-synch
  // MPI-recv-bytes
  // MPI-recv-calls
  // MPI-send-bytes
  // MPI-send-calls
  // MPI-time
  num_basic_counters_ += 9;
#endif

#ifdef LCAMEM
  // lcamem-bytes-curr
  // lcamem-bytes-high
  // lcamem-delete-bytes
  // lcamem-delete-calls
  // lcamem-new-bytes
  // lcamem-new-calls
  num_basic_counters_ += 6;
#endif

  num_basic_counters_ = 24;
  num_papi_counters_ = 0;
  num_user_counters_ = 0;

  ifirst_basic_ = 0;
  ilast_basic_ =  ifirst_basic_ + num_basic_counters_ - 1;
  ifirst_papi_ =  ilast_basic_ + 1; 
  ilast_papi_ =   ifirst_papi_ + num_papi_counters_ - 1;
  ifirst_user_ =  ilast_papi_ + 1;
  ilast_user_ =   ifirst_user_ + num_user_counters_ - 1;

}

//----------------------------------------------------------------------
void LcaPerf::check_indices_() const
//----------------------------------------------------------------------
{
  // Check consistency of indicies
  if ((ilast_papi_ - ifirst_papi_ + 1) != 2*num_papi_counters_) {
    fprintf (stderr, _ERROR 
	     "ifirst_papi_(%d), ilast_papi_(%d) and "
	     "num_papi_counters_(%d) inconsistent!\n",
	     __FILE__,__LINE__,ifirst_papi_,ilast_papi_,num_papi_counters_);
  }
  // Check consistency of indicies
  if ((ilast_user_ - ifirst_user_ + 1) != 2*num_user_counters_) {
    fprintf (stderr, _ERROR 
	     "ifirst_user_(%d), ilast_user_(%d) and "
	     "num_user_counters_(%d) inconsistent!\n",
	     __FILE__,__LINE__,ifirst_user_,ilast_user_,num_user_counters_);
  }
}

//----------------------------------------------------------------------
int LcaPerf::insert_papi_event_(const char * counter_name)
//----------------------------------------------------------------------
{
  int retval;
  if (! papi_is_active_) {
    std::string papi_name = papi_name_ (counter_name);

    // Determine if counter_name corresponds to a PAPI event

#ifdef USE_PAPI
    int papi_code;
#endif
    _CALL_PAPI(PAPI_event_name_to_code
	       ((char *)(papi_name.c_str()),&papi_code),
	       "event_name_to_code",PAPI_OK,retval);
    _CALL_PAPI(PAPI_query_event (papi_code),"query_event",PAPI_OK,retval);

    if (retval == PAPI_OK) {

      // Create PAPI eventset if needed

      if (papi_event_set_ == PAPI_NULL) {
	create_papi_eventset_();
      }

      // Insert the event

      _CALL_PAPI(PAPI_add_event (PAPI3_ARG_ADJUST papi_event_set_,papi_code),
		 "add_event",PAPI_OK,retval);
    } else {

      fprintf (stderr, _ERROR 
	       "insert_papi_event_(%s) called with undefined PAPI event--ignoring\n",
	       __FILE__,__LINE__,counter_name);
      fflush(stderr);

    }
  } else {
      fprintf (stderr, _ERROR 
	       "insert_papi_event_(%s) called with active PAPI eventset!\n",
	       __FILE__,__LINE__,counter_name);
      fflush(stderr);
  }
  return retval;
}

//----------------------------------------------------------------------
int LcaPerf::delete_papi_event_(const char * counter_name)
//----------------------------------------------------------------------
{
  int retval;
  if (! papi_is_active_) {

    std::string papi_name = papi_name_ (counter_name);

    // Determine if counter_name corresponds to a PAPI event

#ifdef USE_PAPI
    int papi_code;
#endif
    _CALL_PAPI(PAPI_event_name_to_code
	       ((char *)(papi_name.c_str()),&papi_code),
	       "event_name_to_code",PAPI_OK,retval);
    _CALL_PAPI(PAPI_query_event (papi_code),"query_event",PAPI_OK,retval);

    // delete the event

    _CALL_PAPI(PAPI3_remove_event (PAPI3_ARG_ADJUST papi_event_set_,papi_code),
	       "rem_event",PAPI_OK,retval);

  } else {
    fprintf (stderr, _ERROR 
	     "delete_papi_event_(%s) called with active PAPI eventset!\n",
	     __FILE__,__LINE__,counter_name);
    fflush(stderr);
  }
  return retval;

}

//----------------------------------------------------------------------
void LcaPerf::clear_papi_eventset_()
//----------------------------------------------------------------------
{
  if (papi_event_set_ != PAPI_NULL) {
    delete_papi_eventset_();
    create_papi_eventset_();
  }
}

//----------------------------------------------------------------------
void LcaPerf::create_papi_eventset_()
//----------------------------------------------------------------------
{
  if (papi_event_set_ == PAPI_NULL) {
    int retval;
    _CALL_PAPI(PAPI_create_eventset(&papi_event_set_),"create_eventset",
	       PAPI_OK,retval);
    if (retval) {
      fprintf (stderr, _ERROR 
	       "PAPI_create_eventset() returned 'retval=%d'\n",
	       __FILE__,__LINE__,retval);
      fflush(stderr);
    }
  } else {
    fprintf (stderr,_DEBUG "Trying to re-create existing PAPI eventset!\n",
	     __FILE__,__LINE__);
  }
}
//----------------------------------------------------------------------
void LcaPerf::delete_papi_eventset_()
//----------------------------------------------------------------------
{
  if (papi_event_set_ != PAPI_NULL) {
    int retval;
    _CALL_PAPI(PAPI_cleanup_eventset( PAPI3_ARG_ADJUST papi_event_set_),
	       "cleanup_eventset",PAPI_OK,retval);
    if (retval) {
      fprintf (stderr, _ERROR 
	       "PAPI_cleanup_eventset() returned 'retval=%d'\n",
	       __FILE__,__LINE__,retval);
      fflush(stderr);
    }
    _CALL_PAPI(PAPI_destroy_eventset(&papi_event_set_),
	       "destroy_eventset", PAPI_OK,retval);
    if (retval) {
      fprintf (stderr, _ERROR 
	       "PAPI_destroy_eventset() returned 'retval=%d'\n",
	       __FILE__,__LINE__,retval);
      fflush(stderr);
    }
    papi_event_set_ = PAPI_NULL;
  }
}

//----------------------------------------------------------------------
void LcaPerf::delete_papi_counter_(const char * counter_name, int counter_type)
//----------------------------------------------------------------------
{
  if (strlen(counter_name) > 0) {
    if (papi_counter_indices_.find(counter_name) != 
	papi_counter_indices_.end()) {

      // Delete given papi_counter if it exists

      int papi_counter_index = papi_counter_indices_[counter_name];
      std::map<std::string,int>::iterator ci;

      adjust_indices_ (LCAPERF_COUNTER_TYPE_PAPI, -1);

      papi_counter_indices_.erase   (papi_counter_indices_.find(counter_name));
      papi_counter_names_.erase     (papi_counter_names_.begin() +
				     papi_counter_index);
      papi_counter_values_.erase    (papi_counter_values_.begin() +
				     papi_counter_index);
      papi_counter_is_active_.erase (papi_counter_is_active_.begin() +
				     papi_counter_index);

      // Remove it from the PAPI eventset

      delete_papi_event_(counter_name);

      // and remember to update indicies to fill in the hole!

      for (ci = papi_counter_indices_.begin();
	   ci != papi_counter_indices_.end();
	   ++ci) {
	if (ci->second > papi_counter_index) {
	  --ci->second;
	}
      }
    } else {
      // Error message if papi_counter does not exist
      fprintf (stderr, _ERROR 
	       "delete_papi_counter_(): %s not being counted!\n",
	       __FILE__,__LINE__,counter_name);
      fflush(stderr);
    }
    
  } else {

    // Delete ALL papi counters

    adjust_indices_(counter_type, -num_papi_counters_);

    papi_counter_values_.clear();
    papi_counter_names_.clear();
    papi_counter_indices_.clear();
    papi_counter_is_active_.clear();

    if (papi_event_set_ != PAPI_NULL) {
      delete_papi_eventset_();
    }
  }
}

//----------------------------------------------------------------------
void LcaPerf::start_papi_counters_ ()
//----------------------------------------------------------------------
{
  if (num_papi_counters_ > 0) {
    if (! papi_is_active_) {
      int retval;
      _CALL_PAPI(PAPI_start(papi_event_set_),"start",PAPI_OK,retval);
      if (retval) {
	fprintf (stderr, _ERROR 
		 "PAPI_start() returned 'retval=%d'\n",
		 __FILE__,__LINE__,retval);
	fflush(stderr);
      }
      papi_is_active_ = true;
    } else {
      fprintf (stderr, _ERROR 
	       "start_papi_counters_() called with active event set!\n",
	       __FILE__,__LINE__);
      fflush(stderr);
    }
  }
}

//----------------------------------------------------------------------
void LcaPerf::stop_papi_counters_ ()
//----------------------------------------------------------------------
{
  if (num_papi_counters_ > 0) {
    if (papi_is_active_) {
      int retval;
#ifdef USE_PAPI
      long long value = 0;
#endif
      papi_is_active_ = false;
      _CALL_PAPI(PAPI_stop(papi_event_set_, &value),
		 "stop",PAPI_OK,retval);
      if (retval) {
	fprintf (stderr, _ERROR 
		 "PAPI_stop() returned 'retval=%d'\n",
		 __FILE__,__LINE__,retval);
	fflush(stderr);
      }
    } else {
      fprintf (stderr, _ERROR 
	       "stop_papi_counters_() called with inactive event set!\n",
	       __FILE__,__LINE__);
      fflush(stderr);
    }
  }
}

//----------------------------------------------------------------------
void LcaPerf::delete_user_counter_(const char * counter_name, int counter_type)
//----------------------------------------------------------------------
{
  if (strlen(counter_name) > 0) {
    if (user_counter_indices_.find(counter_name) != 
	user_counter_indices_.end()) {

      // Delete given user_counter if it exists

      int user_counter_index = user_counter_indices_[counter_name];
      std::map<std::string,int>::iterator ci;

      adjust_indices_ (user_counter_types_[user_counter_index], -1);

      user_counter_indices_.erase(user_counter_indices_.find(counter_name));
      user_counter_names_.erase  (user_counter_names_.begin() + 
				  user_counter_index);
      user_counter_types_.erase  (user_counter_types_.begin() + 
				  user_counter_index);
      user_counter_values_.erase (user_counter_values_.begin() + 
				  user_counter_index);

      // and remember to update indicies to fill in the hole!

      for (ci = user_counter_indices_.begin();
	   ci != user_counter_indices_.end();
	   ++ci) {
	if (ci->second > user_counter_index) {
	  --ci->second;
	}
      }

    } else {
      fprintf (stderr, _ERROR 
	       "delete_user_counter_(%s) called with non-existing user_counter!\n",
	       __FILE__,__LINE__,counter_name);
      fflush(stderr);
    }
    
  } else {
    // Delete ALL user counters

    adjust_indices_(counter_type, -num_user_counters_);

    user_counter_indices_.clear();
    user_counter_names_.clear();
    user_counter_values_.clear();
    user_counter_types_.clear();
    
  }
}

//----------------------------------------------------------------------
std::string LcaPerf::papi_name_ (std::string counter_name)
//----------------------------------------------------------------------
{
  std::string papi_name;

  papi_name = counter_name;
  papi_name.replace(papi_name.find("-"),1,"_");
  papi_name = "PAPI_" + papi_name;

  return papi_name;
}

//----------------------------------------------------------------------
void LcaPerf::new_papi_counter_ (const char * counter_name, int counter_type)
//----------------------------------------------------------------------
{
  int papi_retval = insert_papi_event_(counter_name);

  // INSERT PAPI EVENT

  if (papi_retval == PAPI_OK) {
    int new_size = num_papi_counters_ + 1 ;
    papi_counter_names_.resize(new_size);
    papi_counter_values_.resize(new_size);
    papi_counter_is_active_.resize(new_size);

    papi_counter_indices_  [counter_name]       = num_papi_counters_;
    papi_counter_names_    [num_papi_counters_] = counter_name;
    papi_counter_values_   [num_papi_counters_] = 0;
    papi_counter_is_active_[num_papi_counters_] = true;

    adjust_indices_ (counter_type, +1);
  }
}

//----------------------------------------------------------------------
void LcaPerf::new_user_counter_ (const char * counter_name, int counter_type)
//----------------------------------------------------------------------
{

  // INSERT USER EVENT

  int new_size = num_user_counters_ + 1 ;
  user_counter_names_.resize(new_size);
  user_counter_values_.resize(new_size);
  user_counter_types_.resize(new_size);

  user_counter_indices_[counter_name] = num_user_counters_;
  user_counter_names_[num_user_counters_] = counter_name;
  user_counter_values_[num_user_counters_] = 0;
  user_counter_types_[num_user_counters_] = counter_type;

  adjust_indices_ (counter_type, +1);

}

//----------------------------------------------------------------------
void LcaPerf::write_ (std::string suffix)
//----------------------------------------------------------------------
{

  // Use current suffix_ if parameter is not included

  if (suffix == "") {
    suffix = suffix_;
  }
    
  // Save current path and change to the LcaPerf directory

  getcwd(path_prev_, LCAPERF_MAX_PATH_LENGTH);
  chdir (path_);

  // Open output file

  if (filename_[suffix] == "") {
    // Store filename for future begin(suffix) calls
    char filename[40];
    if (suffix == "") {
      sprintf (filename,"%d",ip_ % 10);
    } else {
      sprintf (filename,"%d.%s",ip_ % 10,suffix.c_str());
    }
    filename_[suffix] = filename;
    fp_[suffix] = fopen (filename, "a");

    write_globals_(fp_[suffix]);
    write_header_(fp_[suffix]);
  }

  // @@@ WARNING: counters may have changed between calls (when 
  // @@@ filename_[suffix] != "") but not checked--file will be invalid!
  // @@@ Currently user must ensure this

  FILE *fp = fp_[suffix];

  // Change directory back to original path

  chdir (path_prev_);

  std::map<std::string,long long *>::iterator pcounters;
  for (pcounters  = counters_global_.begin();
       pcounters != counters_global_.end();
       ++pcounters) {
    write_field_(fp,(*pcounters).first,(*pcounters).second);
  }

  fflush (fp);

}

//----------------------------------------------------------------------
void LcaPerf::start_counters_(long long * counters)
//----------------------------------------------------------------------
{
  int i;
  long long rtime = get_real_ ();
#ifdef USE_PAPI
  long long vtime = get_virt_ ();
#endif

  // Basic counters

  counters[IND_TIME_BEGIN]       = rtime; 
  counters[IND_TIME_ACTIVE_INCL] = rtime;
  counters[IND_TIME_ACTIVE_EXCL] = rtime;

#ifdef USE_PAPI
  counters[IND_VTIME_BEGIN]       = vtime;
  counters[IND_VTIME_ACTIVE_INCL] = vtime;
  counters[IND_VTIME_ACTIVE_EXCL] = vtime;
#endif

#ifdef USE_PMPI
  counters[IND_MPI_GLOBAL_CALLS]      = lcampi_global_calls_;
  counters[IND_MPI_GLOBAL_RECV_BYTES] = lcampi_global_recv_bytes_;
  counters[IND_MPI_GLOBAL_SEND_BYTES] = lcampi_global_send_bytes_;
  counters[IND_MPI_GLOBAL_SYNCH]      = lcampi_global_synch_;
  counters[IND_MPI_RECV_BYTES]        = lcampi_recv_bytes_;
  counters[IND_MPI_RECV_CALLS]        = lcampi_recv_calls_;
  counters[IND_MPI_SEND_BYTES]        = lcampi_send_bytes_;
  counters[IND_MPI_SEND_CALLS]        = lcampi_send_calls_;
  counters[IND_MPI_TIME]              = lcampi_time_;
#endif

#ifdef LCAMEM
  _LCAMEM_TRACE_;
  counters[IND_LCAMEM_BYTES_CURR]   = 0;
  counters[IND_LCAMEM_BYTES_HIGH]   = 0;
  counters[IND_LCAMEM_DELETE_BYTES] = lcamem::deleteBytes_;
  counters[IND_LCAMEM_DELETE_CALLS] = lcamem::deleteCalls_;
  counters[IND_LCAMEM_NEW_BYTES]    = lcamem::newBytes_;
  counters[IND_LCAMEM_NEW_CALLS]    = lcamem::newCalls_;
#endif

  check_indices_();

  // PAPI counters

  read_papi_ ();

  for (i=0; i<num_papi_counters_; i++) {
    int k = ifirst_papi_+2*i;
    counters[k]   = papi_counter_values_[i]; // inclusive counter
    counters[k+1] = papi_counter_values_[i]; // exclusive counter
  }

  // User counters

  for (i=0; i<num_user_counters_; i++) {
    int k = ifirst_user_+2*i;
    if (user_counter_types_[i] == LCAPERF_COUNTER_TYPE_USER_ABS) {
      counters[k]   = 0;                       // inclusive absolute counter
      counters[k+1] = 0;                       // exclusive absolute counter
    } else if (user_counter_types_[i] == LCAPERF_COUNTER_TYPE_USER_REL) {
      counters[k]   = user_counter_values_[i]; // "inclusive" relative counter
      counters[k+1] = user_counter_values_[i]; // "exclusive" relative counter
    } else {
      fprintf (stderr, _ERROR 
	       "start_counters_() called with unknown user counter type %d for user counter %d!\n",
	       __FILE__,__LINE__,user_counter_types_[i], i);
      fflush(stderr);
    }
  }
}

//----------------------------------------------------------------------
void LcaPerf::stop_counters_(long long * counters)
//----------------------------------------------------------------------
{
  _LCAMEM_TRACE_;
  int i;
  long long rtime = get_real_ ();
#ifdef USE_PAPI
  long long vtime = get_virt_ ();
#endif

  // Basic counters

  counters[IND_TIME_END]         = rtime;
  counters[IND_TIME_ACTIVE_INCL] = rtime - counters[IND_TIME_ACTIVE_INCL];
  counters[IND_TIME_ACTIVE_EXCL] = rtime - counters[IND_TIME_ACTIVE_EXCL];

#ifdef USE_PAPI
  counters[IND_VTIME_END]         = vtime;
  counters[IND_VTIME_ACTIVE_INCL] = vtime - counters[IND_VTIME_ACTIVE_INCL];
  counters[IND_VTIME_ACTIVE_EXCL] = vtime - counters[IND_VTIME_ACTIVE_EXCL];
#endif

#ifdef USE_PMPI
  counters[IND_MPI_GLOBAL_CALLS]      = lcampi_global_calls_ - counters[IND_MPI_GLOBAL_CALLS];
  counters[IND_MPI_GLOBAL_RECV_BYTES] = lcampi_global_recv_bytes_ - counters[IND_MPI_GLOBAL_RECV_BYTES];
  counters[IND_MPI_GLOBAL_SEND_BYTES] = lcampi_global_send_bytes_ - counters[IND_MPI_GLOBAL_SEND_BYTES];
  counters[IND_MPI_GLOBAL_SYNCH]      = lcampi_global_synch_ - counters[IND_MPI_GLOBAL_SYNCH];
  counters[IND_MPI_RECV_BYTES]        = lcampi_recv_bytes_   - counters[IND_MPI_RECV_BYTES];
  counters[IND_MPI_RECV_CALLS]        = lcampi_recv_calls_   - counters[IND_MPI_RECV_CALLS];
  counters[IND_MPI_SEND_BYTES]        = lcampi_send_bytes_   - counters[IND_MPI_SEND_BYTES];
  counters[IND_MPI_SEND_CALLS]        = lcampi_send_calls_   - counters[IND_MPI_SEND_CALLS];
  counters[IND_MPI_TIME]              = lcampi_time_         - counters[IND_MPI_TIME];
#endif

#ifdef LCAMEM
  _LCAMEM_TRACE_;
  counters[IND_LCAMEM_BYTES_CURR]   = lcamem::bytes_;
  counters[IND_LCAMEM_BYTES_HIGH]   = lcamem::bytesHigh_;
  counters[IND_LCAMEM_DELETE_BYTES] = lcamem::deleteBytes_ - counters[IND_LCAMEM_DELETE_BYTES];
  counters[IND_LCAMEM_DELETE_CALLS] = lcamem::deleteCalls_ -  counters[IND_LCAMEM_DELETE_CALLS];
  counters[IND_LCAMEM_NEW_BYTES]    = lcamem::newBytes_ -  counters[IND_LCAMEM_NEW_BYTES]   ;
  counters[IND_LCAMEM_NEW_CALLS]    = lcamem::newCalls_ -  counters[IND_LCAMEM_NEW_CALLS]   ;
#endif

  check_indices_();

  // PAPI counters

  read_papi_ ();

  for (i=0; i<num_papi_counters_; i++) {
    int k = ifirst_papi_+2*i;
    counters[k]   = papi_counter_values_[i] - counters[k];
    counters[k+1] = papi_counter_values_[i] - counters[k+1];
  }

  // User counters

  for (i=0; i<num_user_counters_; i++) {
    int k = ifirst_user_+2*i;
    if (user_counter_types_[i] == LCAPERF_COUNTER_TYPE_USER_ABS) {
      counters[k]   = user_counter_values_[i];
    } else if (user_counter_types_[i] == LCAPERF_COUNTER_TYPE_USER_REL) { 
      counters[k]   = user_counter_values_[i] - counters[k];
      counters[k+1] = user_counter_values_[i] - counters[k+1];
    } else {
      fprintf (stderr, _ERROR 
	       "stop_counters_() called with unknown user counter type %d for user counter %d!\n",
	       __FILE__,__LINE__,user_counter_types_[i], i);
      fflush(stderr);
      
    }
  }
}

//----------------------------------------------------------------------
void LcaPerf::update_counters_global_(std::string region,
				     long long * counters_local)
//----------------------------------------------------------------------
{
  int i, isNew = false;

  if (counters_global_[region]==0) {
    counters_global_[region] = create_counters_();
    isNew = true;
  }

  long long * c = counters_global_[region];

  if (isNew) {
    c[IND_TIME_BEGIN]  = counters_local[IND_TIME_BEGIN];
#ifdef USE_PAPI
    c[IND_VTIME_BEGIN] = counters_local[IND_VTIME_BEGIN];
#endif
  }

  c[IND_TIME_END]          = counters_local[IND_TIME_END];
  c[IND_TIME_ACTIVE_INCL] += counters_local[IND_TIME_ACTIVE_INCL];
  c[IND_TIME_ACTIVE_EXCL] += counters_local[IND_TIME_ACTIVE_EXCL];
  ++c[IND_CALL_COUNT];

#ifdef USE_PAPI
  c[IND_VTIME_END]          = counters_local[IND_VTIME_END];
  c[IND_VTIME_ACTIVE_INCL] += counters_local[IND_VTIME_ACTIVE_INCL];
  c[IND_VTIME_ACTIVE_EXCL] += counters_local[IND_VTIME_ACTIVE_EXCL];
#endif

#ifdef USE_PMPI
  c[IND_MPI_GLOBAL_CALLS]      += counters_local[IND_MPI_GLOBAL_CALLS];
  c[IND_MPI_GLOBAL_RECV_BYTES] += counters_local[IND_MPI_GLOBAL_RECV_BYTES];
  c[IND_MPI_GLOBAL_SEND_BYTES] += counters_local[IND_MPI_GLOBAL_SEND_BYTES];
  c[IND_MPI_GLOBAL_SYNCH]      += counters_local[IND_MPI_GLOBAL_SYNCH];
  c[IND_MPI_RECV_BYTES]        += counters_local[IND_MPI_RECV_BYTES];
  c[IND_MPI_RECV_CALLS]        += counters_local[IND_MPI_RECV_CALLS];
  c[IND_MPI_SEND_BYTES]        += counters_local[IND_MPI_SEND_BYTES];
  c[IND_MPI_SEND_CALLS]        += counters_local[IND_MPI_SEND_CALLS];
  c[IND_MPI_TIME]              += counters_local[IND_MPI_TIME];
#endif

#ifdef LCAMEM
  _LCAMEM_TRACE_;
 c[IND_LCAMEM_BYTES_CURR]   = lcamem::bytes_;
 c[IND_LCAMEM_BYTES_HIGH]   = lcamem::bytesHigh_;
 c[IND_LCAMEM_DELETE_BYTES] += counters_local[IND_LCAMEM_DELETE_BYTES];
 c[IND_LCAMEM_DELETE_CALLS] += counters_local[IND_LCAMEM_DELETE_CALLS];
 c[IND_LCAMEM_NEW_BYTES]    += counters_local[IND_LCAMEM_NEW_BYTES]   ;
 c[IND_LCAMEM_NEW_CALLS]    += counters_local[IND_LCAMEM_NEW_CALLS]   ;
#endif

  // Update PAPI counters

  for (i = ifirst_papi_; i <= ilast_papi_; i++) {
    c[i] += counters_local[i];
  }

  // Update User counters

  for (i=0; i<num_user_counters_; i++) {
    int k = ifirst_user_+2*i;
    if (user_counter_types_[i] == LCAPERF_COUNTER_TYPE_USER_ABS) {
      c[k] = counters_local[k];
    } else if (user_counter_types_[i] == LCAPERF_COUNTER_TYPE_USER_REL) { 
      c[k] += counters_local[k];
      c[k+1] += counters_local[k+1];
    } else {
      fprintf (stderr, _ERROR 
	       "update_counters_global_() called with unknown user counter type %d for user counter %d!\n",
	       __FILE__,__LINE__,user_counter_types_[i], i);
      fflush(stderr);
    }
  }

}

//----------------------------------------------------------------------
 void LcaPerf::update_counters_parent_(long long * counters_parent,
				      long long * counters_child)
//----------------------------------------------------------------------
{
  int i;

  // Update Basic counters

  counters_parent[IND_TIME_ACTIVE_EXCL] += counters_child[IND_TIME_ACTIVE_EXCL];
#ifdef USE_PAPI
  counters_parent[IND_VTIME_ACTIVE_EXCL] += counters_child[IND_VTIME_ACTIVE_EXCL];
#endif
  
  // Update PAPI counters

  for (i = ifirst_papi_+1; i <= ilast_papi_; i+=2) {
    counters_parent[i] += counters_child[i];
  }

  // Update User counters

  for (i = ifirst_user_+1; i <= ilast_user_; i+=2) {
    counters_parent[i] += counters_child[i];
  }

}

//----------------------------------------------------------------------
std::string LcaPerf::generate_augmented_region_ (std::string region) const
//----------------------------------------------------------------------
{
  std::string augmented_region;

  augmented_region = region + ":";

  // Set augmented_region = region + [":" + attribute-value]

  for (int i=0; i<num_attributes_; i++) {
    augmented_region = augmented_region + attribute_values_[i] + ":";
  }

  return augmented_region;
}

//----------------------------------------------------------------------
std::string LcaPerf::merge_augmented_regions_ (std::string region1, 
					      std::string region2) const
//----------------------------------------------------------------------
{
  std::string region;

  // Check that region portions match

  if (region1.substr(0,region1.find(":")) != region2.substr(0,region2.find(":"))) {
    fprintf (stderr, _ERROR 
	     "merge_augmented_regions_() regions %s and %s do not match!\n",
	     __FILE__,__LINE__,
	     (region1.substr(0,region1.find(":"))).c_str(),
	     (region2.substr(0,region2.find(":"))).c_str());
    fflush(stderr);
  } else {
    region = region1.substr(0,region1.find(":"));
  }

  // Change differing Attribute values to LCAPERF_NULL_STRING

  std::string tail1 
    = region1.substr(region1.find(":")+1,region1.size()-region1.find(":"));
  std::string tail2 
    = region2.substr(region2.find(":")+1,region2.size()-region2.find(":"));

  while (tail1 != "") {
    if (tail1.substr(0,tail1.find(":")) != tail2.substr(0,tail2.find(":"))) {
      region = region + ":" + LCAPERF_NULL_STRING;
    } else {
      region = region + ":" + tail1.substr(0,tail1.find(":"));
    }

    tail1 = tail1.substr(tail1.find(":")+1,tail1.size()-tail1.find(":"));
    tail2 = tail2.substr(tail2.find(":")+1,tail2.size()-tail2.find(":"));
  }

  return region1;
}

//----------------------------------------------------------------------
void LcaPerf::write_globals_ (FILE *fp) const
//----------------------------------------------------------------------
{
  fprintf (fp,"global processor-count %d\n",np_);
  fprintf (fp,"global processor-rank %d\n",ip_);
  fprintf (fp,"global lcaperf-version %d.%d\n",
	   LCAPERF_VERSION_MAJOR,LCAPERF_VERSION_MINOR);
  std::map<std::string,long long>::const_iterator g = global_values_.begin();
  while (g != global_values_.end()) {
    fprintf (fp,"global %s %lld\n",(g->first).c_str(),(g->second));
    ++g;
  }

  fprintf (fp,"\n");
}

//----------------------------------------------------------------------
void LcaPerf::write_header_ (FILE *fp) const
//----------------------------------------------------------------------
{
  
  int i;
  fprintf (fp,"attribute region\n");
  for (i=0; i<num_attributes_; i++) {
    fprintf (fp,"attribute %s\n",attribute_names_[i].c_str());
  }

  // Basic counters

  fprintf (fp,"basic call-count\n");

  fprintf (fp,"basic time-real-start\n");
  fprintf (fp,"basic time-real-stop\n");
  fprintf (fp,"basic time-real-incl\n");
  fprintf (fp,"basic time-real-excl\n");

  // Basic counters available for free with PAPI

#ifdef USE_PAPI
  fprintf (fp,"basic time-virt-start\n");
  fprintf (fp,"basic time-virt-stop\n");
  fprintf (fp,"basic time-virt-incl\n");
  fprintf (fp,"basic time-virt-excl\n");
#endif

  // Basic counters available with USE_PMPI

#ifdef USE_PMPI
  fprintf (fp,"basic MPI-global-calls\n");
  fprintf (fp,"basic MPI-global-recv-bytes\n");
  fprintf (fp,"basic MPI-global-send-bytes\n");
  fprintf (fp,"basic MPI-global-synch\n");
  fprintf (fp,"basic MPI-recv-bytes\n");
  fprintf (fp,"basic MPI-recv-calls\n");
  fprintf (fp,"basic MPI-send-bytes\n");
  fprintf (fp,"basic MPI-send-calls\n");
  fprintf (fp,"basic MPI-time\n");
#endif

  // Basic counters available with lcamem

#ifdef LCAMEM
  fprintf (fp,"basic lcamem-bytes-curr\n");
  fprintf (fp,"basic lcamem-bytes-high\n");
  fprintf (fp,"basic lcamem-delete-bytes\n");
  fprintf (fp,"basic lcamem-delete-calls\n");
  fprintf (fp,"basic lcamem-new-bytes\n");
  fprintf (fp,"basic lcamem-new-calls\n");
#endif

  // Papi counters

#ifdef USE_PAPI
  for (i=0; i<num_papi_counters_; i++) {
    fprintf (fp,"papi %s-incl\n",papi_counter_names_[i].c_str());
    fprintf (fp,"papi %s-excl\n",papi_counter_names_[i].c_str());
  }
#endif
	    
  // User counters

  for (i=0; i<num_user_counters_; i++) {
    if (user_counter_types_[i] == LCAPERF_COUNTER_TYPE_USER_ABS) {
      fprintf (fp,"user %s-abs\n",user_counter_names_[i].c_str());
    } else if (user_counter_types_[i] == LCAPERF_COUNTER_TYPE_USER_REL) {
      fprintf (fp,"user %s-incl\n",user_counter_names_[i].c_str());
      fprintf (fp,"user %s-excl\n",user_counter_names_[i].c_str());
    } else {
      fprintf (stderr, _ERROR 
	       "write_header_() called with unknown user counter type %d for user counter %d!\n",
	       __FILE__,__LINE__,user_counter_types_[i], i);
      fflush(stderr);
    }
  }

  fprintf (fp,"\n");
}

//----------------------------------------------------------------------
void LcaPerf::write_field_ (FILE *fp, 
			   std::string augmented_region, 
			   long long * counters) const
//----------------------------------------------------------------------
{
  int i;
  std::string ar = augmented_region;

  if (ar.substr(0,ar.find(":")) != LCAPERF_NULL_STRING) {

    // write region name

    fprintf (fp,"%s\n",ar.substr(0,ar.find(":")).c_str());

    // write attribute values

    ar = ar.substr(ar.find(":")+1,ar.size());
    while (ar != "") {
      fprintf (fp,"%s\n",ar.substr(0,ar.find(":")).c_str());
      ar = ar.substr(ar.find(":")+1,ar.size());
    }

    // write basic counters

    fprintf (fp,"%lld\n",counters[IND_CALL_COUNT]);

    fprintf (fp,"%lld\n",counters[IND_TIME_BEGIN] - rtime_begin_);
    fprintf (fp,"%lld\n",counters[IND_TIME_END]   - rtime_begin_);
    fprintf (fp,"%lld\n",counters[IND_TIME_ACTIVE_INCL]);
    fprintf (fp,"%lld\n",counters[IND_TIME_ACTIVE_EXCL]);

    // write basic counters if PAPI is available

#ifdef USE_PAPI
    fprintf (fp,"%lld\n",counters[IND_VTIME_BEGIN] - vtime_begin_);
    fprintf (fp,"%lld\n",counters[IND_VTIME_END]   - vtime_begin_);
    fprintf (fp,"%lld\n",counters[IND_VTIME_ACTIVE_INCL]);
    fprintf (fp,"%lld\n",counters[IND_VTIME_ACTIVE_EXCL]);
#endif

#ifdef USE_PMPI
    fprintf (fp,"%lld\n",counters[IND_MPI_GLOBAL_CALLS]);
    fprintf (fp,"%lld\n",counters[IND_MPI_GLOBAL_RECV_BYTES]);
    fprintf (fp,"%lld\n",counters[IND_MPI_GLOBAL_SEND_BYTES]);
    fprintf (fp,"%lld\n",counters[IND_MPI_GLOBAL_SYNCH]);
    fprintf (fp,"%lld\n",counters[IND_MPI_RECV_BYTES]);
    fprintf (fp,"%lld\n",counters[IND_MPI_RECV_CALLS]);
    fprintf (fp,"%lld\n",counters[IND_MPI_SEND_BYTES]);
    fprintf (fp,"%lld\n",counters[IND_MPI_SEND_CALLS]);
    fprintf (fp,"%lld\n",counters[IND_MPI_TIME]);
#endif

#ifdef LCAMEM
    fprintf (fp,"%lld\n",counters[IND_LCAMEM_BYTES_CURR]);
    fprintf (fp,"%lld\n",counters[IND_LCAMEM_BYTES_HIGH]);
    fprintf (fp,"%lld\n",counters[IND_LCAMEM_DELETE_BYTES]);
    fprintf (fp,"%lld\n",counters[IND_LCAMEM_DELETE_CALLS]);
    fprintf (fp,"%lld\n",counters[IND_LCAMEM_NEW_BYTES]);
    fprintf (fp,"%lld\n",counters[IND_LCAMEM_NEW_CALLS]);
#endif

#ifdef USE_PAPI

    // write papi counters

    for (i=0; i<num_papi_counters_; i++) {
      int k = ifirst_papi_ + 2*i;
      if (papi_counter_is_active_[i]) {
	fprintf (fp,"%lld\n",counters[k]);
	fprintf (fp,"%lld\n",counters[k+1]);
      }
    }
#endif

    // write user counters

    for (i=0; i<num_user_counters_; i++) {
      int k = ifirst_user_ + 2*i;
      fprintf (fp,"%lld\n",counters[k]);
      if (user_counter_types_[i] == LCAPERF_COUNTER_TYPE_USER_REL) {
	fprintf (fp,"%lld\n",counters[k+1]);
      }
    }

    fprintf (fp,"\n");
  }
}

//----------------------------------------------------------------------
void LcaPerf::create_directories_()
//----------------------------------------------------------------------

{
  char path_base [LCAPERF_MAX_PATH_LENGTH];  // LCAPERF directory
  char path_rest [LCAPERF_MAX_PATH_LENGTH];  // New directory to create
  int ip,ip0,ip1,ip2,ip3;
      
  // Expand np_ digits vwxyz into (np0=v,np1=w,np2=x,np3=y)

  int n,n0,n1,n2,n3; 

  n = np_;
  n0 = 0;
  n1 = 0;
  n2 = 0;
  n3 = 0;
  if (n >= 10) {
    n3 = n / 10;
  }	
  if (n >= 100) {
    n3 = 9;
    n2 = n / 100;
  }
  if (n >= 1000) {
    n2 = 9;
    n1 = n / 1000;
  }
  if (n >= 10000) {
    n1 = 9;
    n0 = n / 10000;
  }

  // Create subdirectory structure

  if (ip_ == 0) {

    // create new LCAPERF.# directory

    int dir_count = 1;
    sprintf (lcaperf_dir_,"%s.%03d",LCAPERF_DIR,dir_count);
    int mkdir_error = mkdir (lcaperf_dir_,0755);
    while (mkdir_error && errno == EEXIST) {
      ++ dir_count;
      sprintf (lcaperf_dir_,"%s.%03d",LCAPERF_DIR,dir_count);
      mkdir_error = mkdir (lcaperf_dir_,0755);      
    }

    if (mkdir_error) {
      fprintf (stderr,"Error %d creating %s directory: %s!\n",
	       errno, LCAPERF_DIR,strerror(errno));
      fflush(stderr);
    }
  }

#ifdef USE_MPI
  // Make sure all processors know which directory to use
  MPI_Bcast (lcaperf_dir_, strlen("LCAPERF.###")+1, MPI_CHAR, 0, MPI_COMM_WORLD);
#endif
  
  // Save path base "/.../LCAPERF"

  getcwd(path_base, LCAPERF_MAX_PATH_LENGTH);
  strcat (path_base,"/");
  strcat (path_base,lcaperf_dir_);

  if (ip_ == 0) {

    std::string path;
    int i0,i1,i2,i3;

    // Create subdirectory structure

    bool error_mkdir = false;
    for (i0 = 0; i0 <= n0; i0++) {
      sprintf (path_rest,"/%1d",i0);
      path = (std::string(path_base) + path_rest);
      error_mkdir |= mkdir (path.c_str(),0755);
      for (i1 = 0; i1 <= n1; i1++) {
	sprintf (path_rest,"/%1d/%1d",i0,i1);
	path = (std::string(path_base) + path_rest);
	error_mkdir |= mkdir (path.c_str(),0755);
	for (i2 = 0; i2 <= n2; i2++) {
	  sprintf (path_rest,"/%1d/%1d/%1d",i0,i1,i2);
	  path = (std::string(path_base) + path_rest);
	  error_mkdir |= mkdir (path.c_str(),0755);
	  for (i3 = 0; i3 <= n3; i3++) {
	    sprintf (path_rest,"/%1d/%1d/%1d/%1d",i0,i1,i2,i3);
	    path = (std::string(path_base) + path_rest);
	    error_mkdir |= mkdir (path.c_str(),0755);
	  }
	}
      }
    }

    if (error_mkdir) {
      fprintf (stderr,_ERROR "mkdir error!\n",__FILE__,__LINE__);
    }
  }

  // Barrier to ensure non-root processors don't access
  // nonexistent subdirectories!

#ifdef USE_MPI
  if (!lDirCreated_) {
    MPI_Barrier(MPI_COMM_WORLD);
  }
#endif

  // Save the directory path for this processor

  ip = ip_;
  ip0 = ip /10000;
  ip -= ip0*10000;
  ip1 = ip /1000;
  ip -= ip1*1000;
  ip2 = ip /100;
  ip -= ip2*100;
  ip3 = ip /10;
  ip -= ip3*10;

  sprintf (path_rest,"/%1d/%1d/%1d/%1d",ip0,ip1,ip2,ip3);
  strcpy (path_,path_base);
  strcat (path_,path_rest);

  // Mark directory as created

  lDirCreated_ = true;

}

//======================================================================
// FORTRAN BINDINGS
//======================================================================

static char cstr[LCAPERF_MAX_STRING];
void copy_string (char *str1, char *str2, int len)
{
  assert (abs(len+1) < LCAPERF_MAX_STRING);
  for (int i=0; i<len && !isspace(str2[i]); i++) str1[i] = str2[i];
  str1[len] = 0;
}

//----------------------------------------------------------------------

extern "C" void 
LCAPERF_FORTRAN(lcac_new_attribute) 
  (char *     attribute_name,
   int        *len,
   int        *pattribute_type)
{
  copy_string(cstr,attribute_name,*len);
  lcaperf.new_attribute(cstr,*pattribute_type);
}

//----------------------------------------------------------------------

extern "C" void 
LCAPERF_FORTRAN(lcac_new_counter) 
  (char *     counter_name, 
   int        *len,
   int        *pcounter_type)
{
  copy_string(cstr,counter_name,*len);
  lcaperf.new_counter(cstr,*pcounter_type);
}

//----------------------------------------------------------------------

extern "C" void 
LCAPERF_FORTRAN(lcac_delete_attribute) 
  (char *     attribute_name,    
   int        *len) 
{
  copy_string(cstr,attribute_name,*len);
  lcaperf.delete_attribute(cstr);
}

//----------------------------------------------------------------------

extern "C" void 
LCAPERF_FORTRAN(lcac_delete_counter) 
  (char *     counter_name,    
   int        *len,
   int        *pcounter_type) 
{
  copy_string(cstr,counter_name,*len);
  lcaperf.delete_counter(cstr,*pcounter_type);
}

//----------------------------------------------------------------------

extern "C" void 
LCAPERF_FORTRAN(lcac_begin) 
  (char *     segment,
   int        *len) 
{
  copy_string(cstr,segment,*len);
  lcaperf.begin(cstr);
}

//----------------------------------------------------------------------

extern "C" void 
LCAPERF_FORTRAN(lcac_end)
  (char *     segment,
   int        *len) 
{
  copy_string(cstr,segment,*len);
  lcaperf.end(cstr);
}

//----------------------------------------------------------------------

extern "C" void 
LCAPERF_FORTRAN(lcac_global) 
  (char *      global_name,   
   int         *len,   
   long long   *pvalue)
{
  copy_string(cstr,global_name,*len);
  lcaperf.global(cstr,*pvalue);
}

//----------------------------------------------------------------------

extern "C" void 
LCAPERF_FORTRAN(lcac_start) 
  (char *     region_name,
   int        *len)
{
  copy_string(cstr,region_name,*len);
  lcaperf.start(cstr);
}

//----------------------------------------------------------------------

extern "C" void 
LCAPERF_FORTRAN(lcac_stop) 
  (char *     region_name, 
   int        *len)
{
  copy_string(cstr,region_name,*len);
  lcaperf.stop(cstr);
}

//----------------------------------------------------------------------

extern "C" void 
LCAPERF_FORTRAN(lcac_attribute) 
  (char *     attribute_name, 
   int        *len, 
   void       *pvalue_pointer,  /* WARNING--passing pointer as integer */
   int        *pvalue_type)
{
  copy_string(cstr,attribute_name,*len);
  lcaperf.attribute(cstr,pvalue_pointer, *pvalue_type);
}

//----------------------------------------------------------------------

extern "C" void 
LCAPERF_FORTRAN(lcac_increment) 
  (char *      counter_name,   
   int         *len,   
   long long   *pvalue)
{
  copy_string(cstr,counter_name,*len);
  lcaperf.increment(cstr,*pvalue);
}

//----------------------------------------------------------------------

extern "C" void 
LCAPERF_FORTRAN(lcac_assign) 
  (char *      counter_name,   
   int         *len,   
   long long   *pvalue)
{
  copy_string(cstr,counter_name,*len);
  lcaperf.assign(cstr,*pvalue);
}

//----------------------------------------------------------------------

LcaPerf lcaperf;

