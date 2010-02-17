/***********************************************************************
/
/  MEMORY POOL ABSTRACT CLASS
/
/  written by: John Wise
/  date:       February, 2010
/
/  PURPOSE:    Custom memory manager to avoid heap fragmentation, 
/              especially with photon packages.
/
************************************************************************/
#ifndef __MEMORYBLOCK_H
#define __MEMORYBLOCK_H

namespace MPool
{

  // short typedef for a byte
  typedef unsigned char TByte;

  /* Abstract base class (interface) for the memory pool. */

  class MemoryBlock
  {
  public:
    virtual ~MemoryBlock(void) {};
    virtual void* GetMemory(const size_t &MemorySize)=0;
    virtual void FreeMemory(void* sMemoryBlock)=0;
  };

}

#endif
