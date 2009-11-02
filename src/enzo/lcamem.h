#ifndef LCAMEM_H
#define LCAMEM_H
//======================================================================
//
// File:        lcamem.h
//
// Description: Overrides new() and delete [] () to track memory allocation
//
//----------------------------------------------------------------------
//
// Classes:     lcamem
//
//----------------------------------------------------------------------
//
// Copyright 2004 James Bordner
// Copyright 2004 Laboratory for Computational Astrophysics
// Copyright 2004 Regents of the University of California
//
//======================================================================

namespace lcamem {
  void *new_(size_t bytes) throw (std::bad_alloc);
  void delete_(void *p) throw ();
  void clear_();
  extern long long bytes_;
  extern long long bytesHigh_;
  extern long long newCalls_;
  extern long long newBytes_;
  extern long long deleteCalls_;
  extern long long deleteBytes_;
}  

extern const bool lcamem_define_arrays_only;

#endif
