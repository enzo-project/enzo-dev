/***********************************************************************
/
/  MEMORY POOL - MEMORY CHUNK STRUCTURE
/
/  written by: John Wise
/  date:       February, 2010
/
/  PURPOSE:    Custom memory manager to avoid heap fragmentation, 
/              especially with photon packages.
/
************************************************************************/
#ifndef __MEMORYCHUNK_H
#define __MEMORYCHUNK_H

#include "MemoryBlock.h"

namespace MPool
{

  /* 
     This struct will hold and manage the actual allocated
     memory. Every MemoryChunk will point to a MemoryBlock and the
     next MemoryChunk in the linked list. 
  */

  typedef struct MemoryChunk
  {
    TByte* Data;  // data
    size_t DataSize;  // size of data
    size_t UsedSize;  // actual used size

    // true when this MemoryChunk's Data can be deallocated with free()
    bool IsAllocationChunk;

    // Next block in the linked list
    MemoryChunk* Next;
  } MemoryChunk;

}

#endif
