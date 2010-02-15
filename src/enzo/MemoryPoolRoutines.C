/***********************************************************************
/
/  MEMORY POOL CLASS ROUTINES
/
/  written by: John Wise
/  date:       February, 2010
/
/  PURPOSE:    Custom memory manager to avoid heap fragmentation, 
/              especially with photon packages.
/
************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "MemoryPool.h"

namespace MPool
{

  // For debugging (SetMemoryData)
  static const int FREED_MEMORY = 0xAA;
  static const int NEW_ALLOCATED_MEMORY = 0xFF;

  /******************************/
  /* CONSTRUCTOR AND DESTRUCTOR */
  /******************************/

  MemoryPool::MemoryPool(const size_t &sInitialMemoryPoolSize,
			 const size_t &sMemoryChunkSize,
			 const size_t &sMinimalMemorySizeToAllocate,
			 bool bSetMemoryData)
  {
    FirstChunk = NULL;
    LastChunk = NULL;
    CursorChunk = NULL;
    TotalMemoryPoolSize = 0;
    UsedMemoryPoolSize = 0;
    FreeMemoryPoolSize = 0;
    MemoryChunkSize = sMemoryChunkSize;
    MemoryChunkCount = 0;
    ObjectCount = 0;
    SetMemoryData = bSetMemoryData;
    MinimalMemorySizeToAllocate = sMinimalMemorySizeToAllocate;
    
    this->AllocateMemory(sInitialMemoryPoolSize);
    
  }
  
  MemoryPool::~MemoryPool(void)
  {
    FreeAllAllocatedMemory();
    DeallocateAllChunks();
    //assert(ObjectCount == 0);  // Check if everything was deallocated
  }

  /**************************************************/
  /*                ALL OTHER ROUTINES              */
  /**************************************************/

  void* MemoryPool::GetMemory(const size_t &MemorySize)
  {
    size_t BestMemBlockSize = CalculateBestMemoryBlockSize(MemorySize);
    size_t AllocationBlockSize;
    MemoryChunk* Chunk = NULL;

    while (Chunk == NULL) {

      // Are there chunks available to hold the memory?
      Chunk = FindChunkSuitableToHoldMemory(BestMemBlockSize);
      if (Chunk == NULL) {
	// No chunk was found -> Memory pool is too small.  Request
	// more memory from the OS.
	AllocationBlockSize = 
	  max(BestMemBlockSize, 
	      CalculateBestMemoryBlockSize(MinimalMemorySizeToAllocate));
	AllocateMemory(AllocationBlockSize);
      } // ENDIF

    } // ENDWHILE

    /* After we found a suitable chunk, adjust the internal variables */
    
    UsedMemoryPoolSize += BestMemBlockSize;
    FreeMemoryPoolSize -= BestMemBlockSize;
    ObjectCount++;
    SetMemoryChunkValues(Chunk, BestMemBlockSize);

    // Return the chunk pointer
    return ((void*) Chunk->Data);

  }

  /**********************************************************************/
  
  void MemoryPool::FreeMemory(void* sMemoryBlock)
  {
    MemoryChunk* Chunk = FindChunkHoldingPointer(sMemoryBlock);
    if (Chunk != NULL)
      FreeChunks(Chunk);
    else
      assert(false && "Error: requested pointer not in memory pool");
    assert(ObjectCount > 0 && 
	   "Error: requested to delete more memory than allocated.");
    ObjectCount--;
  }

  /**********************************************************************/

  bool MemoryPool::AllocateMemory(const size_t &MemorySize)
  {
    unsigned int NeededChunks = CalculateNeededChunks(MemorySize);
    size_t BestMemBlockSize = CalculateBestMemoryBlockSize(MemorySize);
    
    // Allocate memory from the OS
    TByte* NewMemBlock = (TByte*) malloc(BestMemBlockSize);

    // Allocate chunk array to manage the memory
    MemoryChunk* NewChunks = 
      (MemoryChunk*) malloc(NeededChunks * sizeof(MemoryChunk));
    assert(((NewMemBlock) && (NewChunks)) && "Error: out of memory?");

    printf("AllocateMemory:\n"
	   "\tNeededChunks = %d\n"
	   "\tBestMemBlockSize = %d\n"
	   "\tNewChunks = %x\n",
	   NeededChunks, BestMemBlockSize, NewChunks);

    // Adjust the internal values
    TotalMemoryPoolSize += BestMemBlockSize;
    FreeMemoryPoolSize += BestMemBlockSize;
    MemoryChunkCount += NeededChunks;

    if (SetMemoryData)
      memset( ((void*) NewMemBlock), NEW_ALLOCATED_MEMORY, BestMemBlockSize);

    return LinkChunksToData(NewChunks, NeededChunks, NewMemBlock);
    
  }

  /**********************************************************************/

  unsigned int MemoryPool::CalculateNeededChunks(const size_t &MemorySize)
  {
    float f = ((float) MemorySize) / ((float) MemoryChunkSize);
    return ((unsigned int) ceil(f));
  }
  
  /**********************************************************************/
  
  unsigned int MemoryPool::CalculateBestMemoryBlockSize
  (const size_t &RequestedMemoryBlockSize)
  {
    unsigned int NeededChunks =
      CalculateNeededChunks(RequestedMemoryBlockSize);
    return ((unsigned int)(NeededChunks * MemoryChunkSize));
  }

  /**********************************************************************/

  void MemoryPool::FreeChunks(MemoryChunk* Chunk)
  {
    MemoryChunk* CurrentChunk = Chunk;
    unsigned int ChunkCount = CalculateNeededChunks(CurrentChunk->UsedSize);
    unsigned int i;
    for (i = 0; i < ChunkCount; i++)
      if (CurrentChunk) {
	// Optional: set allocated memory to FREED_MEMORY
	if (SetMemoryData)
	  memset( ((void*) CurrentChunk->Data), FREED_MEMORY, MemoryChunkSize);
	
	// Set used size to zero and adjust memory pool values.  Then
	// go to the next chunk.
	CurrentChunk->UsedSize = 0;
	UsedMemoryPoolSize -= MemoryChunkSize;
	FreeMemoryPoolSize += MemoryChunkSize;
	CurrentChunk = CurrentChunk->Next;
      } // ENDIF CurrentChunk
  }

  /**********************************************************************/

  MemoryChunk* MemoryPool::FindChunkSuitableToHoldMemory
  (const size_t &MemorySize)
  {

    unsigned int i, j, NeededChunks, ChunksToSkip = 0;
    bool EverythingFree = true;
    MemoryChunk* Chunk = CursorChunk;
    MemoryChunk* ThisChunk;
    NeededChunks = CalculateNeededChunks(MemorySize);
    for (i = 0; i < MemoryChunkCount; i++) {
      if (Chunk->DataSize >= MemorySize) {
	// Check if memory is free in the following "NeededChunks" chunks
	EverythingFree = true;
	for (j = 0, ThisChunk = Chunk; j < NeededChunks; j++) {
	  if (ThisChunk->UsedSize == 0) {
	    ThisChunk = ThisChunk->Next;
	  } else {
	    EverythingFree = false;
	    break;
	  }
	} // ENDFOR j
	if (EverythingFree) {
	  CursorChunk = Chunk;
	  return Chunk;
	}
      } // ENDIF DataSize >= MemorySize
      ChunksToSkip = CalculateNeededChunks(Chunk->UsedSize);
      if (ChunksToSkip == 0) ChunksToSkip = 1;
      Chunk = SkipChunks(Chunk, ChunksToSkip);
//      if (Chunk == NULL) // End of list -> start over at the beginning
//	Chunk = FirstChunk;
    } // ENDFOR

    return NULL;

  }

  /**********************************************************************/

  MemoryChunk* MemoryPool::SkipChunks(MemoryChunk* StartChunk,
				      unsigned int ChunksToSkip)
  {

    unsigned int i;
    MemoryChunk* CurrentChunk = StartChunk;
    for (i = 0; i < ChunksToSkip; i++)
      if (CurrentChunk) {
	CurrentChunk = CurrentChunk->Next;
	if (CurrentChunk == NULL) // End of the list
	  CurrentChunk = FirstChunk;
      } 
      else {
	// Error if tries to skip more chunks than available from the
	// start chunk
	assert(false && "SkipChunks: ERROR.  NULL pointed not expected.");
	break;
      }
    return CurrentChunk;
  }

  /**********************************************************************/

  void MemoryPool::SetMemoryChunkValues(MemoryChunk* Chunk,
					const size_t &MemBlockSize)
  {
    if (Chunk)
      Chunk->UsedSize = MemBlockSize;
    else
      assert(false && "SetMemoryChunkValues: ERROR.  NULL pointer?");
  }

  /**********************************************************************/

  bool MemoryPool::LinkChunksToData(MemoryChunk* NewChunks,
				    unsigned int ChunkCount,
				    TByte* NewMemBlock)
  {
    MemoryChunk* NewChunk = NULL;
    unsigned int i, MemOffset = 0;
    bool AllocationChunkAssigned = false;
    for (i = 0; i < ChunkCount; i++) {
      if (!FirstChunk) {
	FirstChunk = SetChunkDefaults(&(NewChunks[0]));
	LastChunk = FirstChunk;
	CursorChunk = FirstChunk;
      } else {
	NewChunk = SetChunkDefaults(&(NewChunks[i]));
	LastChunk->Next = NewChunk;
	LastChunk = NewChunk;
      } // ENDELSE !FirstChunk

      MemOffset = i * ((unsigned int) MemoryChunkSize);
      LastChunk->Data = &(NewMemBlock[MemOffset]);

      /* The first chunk assigned to the new memory block will be an
	 "allocation chunk".  Meaning that this chunk stores the
	 original pointer to the MemBlock and is responsible for
	 freeing the memory later. */
      if (!AllocationChunkAssigned) {
	LastChunk->IsAllocationChunk = true;
	AllocationChunkAssigned = true;
      }
    } // ENDFOR i
    //return RecalcChunkMemorySize(FirstChunk, MemoryChunkCount);
    return RecalcChunkMemorySize(NewChunks, ChunkCount);
  }

  /**********************************************************************/

  bool MemoryPool::RecalcChunkMemorySize(MemoryChunk* Chunk,
					 unsigned int ChunkCount)
  {
    unsigned int i, MemOffset = 0;
    unsigned int TotalChunkMemory = ChunkCount * ((unsigned int) MemoryChunkSize);
    for (i = 0; i < ChunkCount; i++)
      if (Chunk) {
	MemOffset = i * ((unsigned int) MemoryChunkSize);
	Chunk->DataSize = TotalChunkMemory - MemOffset;
	Chunk = Chunk->Next;
      } else {
	assert(false && "RecalcChunkMemorySize: ERROR.  NULL pointer?");
	return false;
      }
    return true;
  }

  /**********************************************************************/

  MemoryChunk* MemoryPool::SetChunkDefaults(MemoryChunk* Chunk)
  {
    if (Chunk) {
      Chunk->Data = NULL;
      Chunk->DataSize = 0;
      Chunk->UsedSize = 0;
      Chunk->IsAllocationChunk = false;
      Chunk->Next = NULL;
    }
    return Chunk;
  }

  /**********************************************************************/

  MemoryChunk* MemoryPool::FindChunkHoldingPointer(void* sMemoryBlock)
  {
    MemoryChunk* TempChunk = FirstChunk;
    while (TempChunk) {
      if (TempChunk->Data == ((TByte*) sMemoryBlock))
	break;
      TempChunk = TempChunk->Next;
    }
    return TempChunk;
  }

  /**********************************************************************/

  void MemoryPool::FreeAllAllocatedMemory(void)
  {
    MemoryChunk* Chunk = FirstChunk;
    while (Chunk) {
      if (Chunk->IsAllocationChunk)
	free((void*) Chunk->Data);
      Chunk = Chunk->Next;
    } // ENDWHILE
  }

  /**********************************************************************/

  void MemoryPool::DeallocateAllChunks(void)
  {
    MemoryChunk* Chunk = FirstChunk;
    MemoryChunk* ChunkToDelete = NULL;
    while (Chunk) {
      if (Chunk->IsAllocationChunk) {
	if (ChunkToDelete)
	  free((void*) ChunkToDelete);
	ChunkToDelete = Chunk;
      }
      Chunk = Chunk->Next;
    } // ENDWHILE
    if (ChunkToDelete)
      free((void*) ChunkToDelete);
  }

  /**********************************************************************/

  bool MemoryPool::IsValidPointer(void* Pointer)
  {
    MemoryChunk* Chunk = FirstChunk;
    while (Chunk) {
      if (Chunk->Data == ((TByte*) Pointer))
	return true;
      Chunk = Chunk->Next;
    }
    return false;
  }

  /**********************************************************************/

}
