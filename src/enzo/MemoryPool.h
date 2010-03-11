/***********************************************************************
/
/  MEMORY POOL CLASS
/
/  written by: John Wise
/  date:       February, 2010
/
/  PURPOSE:    Custom memory manager to avoid heap fragmentation, 
/              especially with photon packages.
/
************************************************************************/
#ifndef __MEMORYPOOL_H
#define __MEMORYPOOL_H

#include "MemoryBlock.h"
#include "MemoryChunk.h"

namespace MPool
{

  static const size_t DEFAULT_MEMORY_POOL_SIZE = 1024;  // Bytes
  static const size_t DEFAULT_MEMORY_CHUNK_SIZE = 16;   // Bytes
  static const size_t DEFAULT_MEMORY_ALLOCATION = 16*DEFAULT_MEMORY_CHUNK_SIZE;

  class MemoryPool : public virtual MemoryBlock
  {

  private:
    
    // Allocates MemorySize bytes from the OS
    bool AllocateMemory(const size_t &MemorySize);

    // Frees all allocated memory back to the OS
    void FreeAllAllocatedMemory(void);

    // Returns number of memory chunks to hold MemorySize bytes
    size_t CalculateNeededChunks(const size_t &MemorySize);

    // Returns the amount of memory that's best managed by the memory chunks
    size_t CalculateBestMemoryBlockSize
    (const size_t &RequestedMemoryBlockSize);
    
    // Returns a chunk that can hold MemorySize bytes
    MemoryChunk* FindChunkSuitableToHoldMemory(const size_t &MemorySize);
    
    // Finds a chunk whose data points to this pointer.  If not found,
    // return NULL.
    MemoryChunk* FindChunkHoldingPointer(void *sMemoryBlock);

    // Skip the given amount of chunks, starting at the given starting
    // pointer
    MemoryChunk* SkipChunks(MemoryChunk* StartChunk, size_t ChunksToSkip);
    
    // Set default values to the Chunk
    MemoryChunk* SetChunkDefaults(MemoryChunk* Chunk);

    // Frees memory linked to the given chunk back to the MemoryPool again.
    void FreeChunks(MemoryChunk* Chunk);

    // Deallocates all memory chunks back to the OS
    void DeallocateAllChunks(void);

    // Link the given block to the linked list of memory chunks
    bool LinkChunksToData(MemoryChunk* NewChunk, size_t ChunkCount,
			  TByte* NewMemBlock);

    // Sets the UsedSize of this chunk pointer to MemBlockSize
    void SetMemoryChunkValues(MemoryChunk* Chunk,
			      const size_t &MemBlockSize);

    // Recalculates the DataSize of all chunks when the memory pool
    // grows by AllocateMemory()
    bool RecalcChunkMemorySize(MemoryChunk* Chunks, size_t ChunkCount);

    // Pointer to the first and last chunk
    MemoryChunk* FirstChunk;
    MemoryChunk* LastChunk;

    // Cursor chunk that speeds up navigation in the linked list
    MemoryChunk* CursorChunk;

    // Total size of the memory pool in bytes
    size_t TotalMemoryPoolSize;
    
    // Used and Free memory in bytes
    size_t UsedMemoryPoolSize;
    size_t FreeMemoryPoolSize;
    
    // Amount of memory in a single chunk
    size_t MemoryChunkSize;

    // Total amount of memory chunks
    size_t MemoryChunkCount;
    
    // Counter for GetMemory() and FreeMemory() calls.  Indirectly
    // counts the number of objects inside the memory pool.
    size_t ObjectCount;

    // Set to true, if all (de)allocated memory is set to a predefined
    // value through memset().  Useful for debugging.
    bool SetMemoryData;

    // Minimal amount of memory that can allocated by AllocateMemory()
    size_t MinimalMemorySizeToAllocate;

  public:

    /* CONSTRUCTOR:
       InitializeMemoryPoolSize: initial size (in bytes) of the memory pool

       MemoryChunkSize: size of each memory chunk.  Low values
                        increases runtime but decreases memory
                        overhead and fragmentation.

       MinimalMemorySizeToAllocate: minimum size to allocate when requesting
                                    more memory from the OS

       SetMemoryData: when true, you set the allocated/freed memory to some
                      value, which is useful for debugging.  Bad for runtimes.
     */

    MemoryPool(const size_t &sInitialMemoryPoolSize = 
	       DEFAULT_MEMORY_POOL_SIZE,
	       const size_t &sMemoryChunkSize = DEFAULT_MEMORY_CHUNK_SIZE,
	       const size_t &sMinimalMemorySizeToAllocate =
	       DEFAULT_MEMORY_ALLOCATION,
	       bool bSetMemoryData = false);

    /* DESTRUCTOR */
    ~MemoryPool(void);

    // Gets MemorySize bytes from memory pool
    void* GetMemory(const size_t &MemorySize);

    // Free the allocated memory
    void FreeMemory(void *sMemoryBlock);

    // Checks if the pointer is in the memory pool
    bool IsValidPointer(void* Pointer);

  };
}

#endif
