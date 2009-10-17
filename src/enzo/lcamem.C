/*****************************************************************************
 *                                                                           *
 * Copyright 2004 James Bordner                                              *
 * Copyright 2004 Laboratory for Computational Astrophysics                  *
 * Copyright 2004 Board of Trustees of the University of Illinois            *
 * Copyright 2004 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <new>

#include "lcamem.h"

#undef TRACE

/* don't trace new/delete() just new[]()/delete[]() */

#undef ARRAYS_ONLY

#ifdef ARRAYS_ONLY
  const bool lcamem_define_arrays_only = true;
#else
  const bool lcamem_define_arrays_only = false;
#endif

namespace lcamem {
  void *new_(size_t bytes) throw (std::bad_alloc);
  void delete_(void *p) throw ();
  long long bytes_ = 0;
  long long bytesHigh_ = 0;
  long long newCalls_ = 0;
  long long newBytes_ = 0;
  long long deleteCalls_ = 0;
  long long deleteBytes_ = 0;
}

//----------------------------------------------------------------------

void lcamem::clear_ ()
{
  lcamem::bytes_ = 0;
  lcamem::bytesHigh_ = 0;
  lcamem::newCalls_ = 0;
  lcamem::newBytes_ = 0;
  lcamem::deleteCalls_ = 0;
  lcamem::deleteBytes_ = 0;
  return;
}

//----------------------------------------------------------------------

void *lcamem::new_(size_t bytes) throw (std::bad_alloc)
{
  if (bytes==0) return NULL;

  long long *pi;

  pi = (long long *)(malloc(bytes+sizeof(long long)));
  if (pi==0) throw std::bad_alloc();

  pi[0] = bytes;

  ++ lcamem::newCalls_ ;
  lcamem::newBytes_ += bytes;
  lcamem::bytes_ += bytes;
  if (lcamem::bytes_ > lcamem::bytesHigh_) lcamem::bytesHigh_ = lcamem::bytes_;

  return (void *)(pi+1);
}

//----------------------------------------------------------------------

void lcamem::delete_ (void *p) throw ()

{
  // Assumes p != 0

  long long *pi = (long long *)(p) - 1;

  ++ lcamem::deleteCalls_ ;
  lcamem::deleteBytes_ += pi[0];
  lcamem::bytes_ -= pi[0];

  free(pi);
}

//----------------------------------------------------------------------

void *operator new (size_t bytes) throw (std::bad_alloc)
{

#ifdef TRACE
  printf ("Entering lcamem new (%d) (%lld,%lld)\n",(int) bytes,
	  lcamem::bytes_, lcamem::bytesHigh_);
#endif

  long long * p = (long long *) lcamem::new_(bytes);

#ifdef TRACE
  printf ("Exiting lcamem new (%p) (%lld,%lld)\n",p,
	  lcamem::bytes_, lcamem::bytesHigh_);
#endif

  // Return pointer to new storage

  return (void *) p;
}

//----------------------------------------------------------------------
#ifndef ARRAYS_ONLY
//----------------------------------------------------------------------

void operator delete (void *p) throw ()
{
  if (p==0) return;

#ifdef TRACE
  printf ("Entering lcamem delete (%lld,%p) (%lld,%lld)\n",*((long long *)(p)-1),p,
	  lcamem::bytes_, lcamem::bytesHigh_);
#endif

  lcamem::delete_(p);

#ifdef TRACE
  printf ("Exiting lcamem delete (%lld) (%lld,%lld)\n",*((long long *)(p)-1),
	  lcamem::bytes_, lcamem::bytesHigh_);
#endif
}

//----------------------------------------------------------------------
#endif
//----------------------------------------------------------------------

void *operator new [] (size_t bytes) throw (std::bad_alloc)
{
#ifdef TRACE
  printf ("Entering lcamem new [] (%d) (%lld,%lld)\n",(int)bytes,
	  lcamem::bytes_, lcamem::bytesHigh_);
#endif

  long long * p = (long long *) lcamem::new_(bytes);

#ifdef TRACE
  printf ("Exiting lcamem new [] (%p) (%lld,%lld)\n",p,
	  lcamem::bytes_, lcamem::bytesHigh_);
#endif

  // Return pointer to new storage

  return (void *)(p);

}

//----------------------------------------------------------------------

void operator delete [] (void *p) throw ()
{
  if (p==0) return;

#ifdef TRACE
  printf ("Entering lcamem delete [] (%lld,%p) (%lld,%lld)\n",*((long long *)(p)-1),p,
	  lcamem::bytes_, lcamem::bytesHigh_);
#endif

  lcamem::delete_(p);

#ifdef TRACE
  printf ("Exiting lcamem delete [] (%lld) (%lld,%lld)\n",*((long long *)(p)-1),
	  lcamem::bytes_, lcamem::bytesHigh_);
#endif
}



