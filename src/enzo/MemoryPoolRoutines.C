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
#include "global_data.h"
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

    /* Chunks are allocated with the "Data" memory, so it's
       deallocated with it, too. */
    //DeallocateAllChunks();
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
    size_t NeededChunks = CalculateNeededChunks(MemorySize);
    size_t BestMemBlockSize = CalculateBestMemoryBlockSize(MemorySize);
    
    // Allocate memory from the OS for data and chunks
    void* MemBlock = malloc(BestMemBlockSize +
			    NeededChunks * sizeof(MemoryChunk));
    //TByte* NewMemBlock = (TByte*) malloc(BestMemBlockSize);
    TByte* NewMemBlock = (TByte*) MemBlock;

    // Allocate chunk array to manage the memory
    //MemoryChunk* NewChunks = 
    //  (MemoryChunk*) malloc(NeededChunks * sizeof(MemoryChunk));
    MemoryChunk* NewChunks = (MemoryChunk*) ((TByte*)(MemBlock)+BestMemBlockSize);
    assert(((NewMemBlock) && (NewChunks)) && "Error: out of memory?");

    // Adjust the internal values
    TotalMemoryPoolSize += BestMemBlockSize;
    FreeMemoryPoolSize += BestMemBlockSize;
    MemoryChunkCount += NeededChunks;

#ifdef MEM_TRACE
    printf("P%d: AllocateMemory: %0.3f MB, (+%0.3f MB)\n",
	   MyProcessorNumber, TotalMemoryPoolSize/1048576.0,
	   BestMemBlockSize/1048576.0);
#endif

    if (SetMemoryData)
      memset( ((void*) NewMemBlock), NEW_ALLOCATED_MEMORY, BestMemBlockSize);

    return LinkChunksToData(NewChunks, NeededChunks, NewMemBlock);
    
  }

  /**********************************************************************/

  size_t MemoryPool::CalculateNeededChunks(const size_t &MemorySize)
  {
    float f = ((float) MemorySize) / ((float) MemoryChunkSize);
    return ((size_t) ceil(f));
  }
  
  /**********************************************************************/
  
  size_t MemoryPool::CalculateBestMemoryBlockSize
  (const size_t &RequestedMemoryBlockSize)
  {
    size_t NeededChunks =
      CalculateNeededChunks(RequestedMemoryBlockSize);
    return ((size_t)(NeededChunks * MemoryChunkSize));
  }

  /**********************************************************************/

  void MemoryPool::FreeChunks(MemoryChunk* Chunk)
  {
    MemoryChunk* CurrentChunk = Chunk;
    size_t ChunkCount = CalculateNeededChunks(CurrentChunk->UsedSize);
    size_t i;
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

    size_t i, j, NeededChunks, ChunksToSkip = 0;
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
				      size_t ChunksToSkip)
  {

    size_t i;
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
				    size_t ChunkCount,
				    TByte* NewMemBlock)
  {
    MemoryChunk* NewChunk = NULL;
    size_t i, MemOffset = 0;
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

      MemOffset = i * ((size_t) MemoryChunkSize);
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
					 size_t ChunkCount)
  {
    size_t i, MemOffset = 0;
    size_t TotalChunkMemory = ChunkCount * ((size_t) MemoryChunkSize);
    for (i = 0; i < ChunkCount; i++)
      if (Chunk) {
	MemOffset = i * ((size_t) MemoryChunkSize);
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

  inline MemoryChunk* MemoryPool::FindChunkHoldingPointer(void* sMemoryBlock)
  {
    /* Memory blocks and chunks were allocated together, so we can do
       some pointer arithemtic to find the chunk. */
    
    size_t ThisDataSize, nChunks, ChunkNumber;
    TByte *mFirstBlock, *mLastBlock, *mem;
    MemoryChunk *FirstChunkInAlloc, *LastChunkInAlloc;
    mem = (TByte*) sMemoryBlock;

    // Find the memory allocation "big" block
    mFirstBlock = FirstChunk->Data;
    //ThisDataSize = FirstChunk->DataSize;
    nChunks = FirstChunk->DataSize / MemoryChunkSize;
    mLastBlock = mFirstBlock + FirstChunk->DataSize;

    /* Go to the next allocation block.  First go directly to the last
       MemoryChunk in the current allocation block.  Then use the Next
       pointer.  Remember that the allocation blocks aren't contigious. */

    while (mem < mFirstBlock || mem >= mLastBlock) {

      LastChunkInAlloc = (MemoryChunk*) (mLastBlock + 
					 (nChunks-1) * sizeof(MemoryChunk));
      FirstChunkInAlloc = LastChunkInAlloc->Next;
      assert((FirstChunkInAlloc) && "Cannot find chunk.");

      mFirstBlock = FirstChunkInAlloc->Data;
      //ThisDataSize = FirstChunkInAlloc->DataSize;
      nChunks = FirstChunkInAlloc->DataSize / MemoryChunkSize;
      mLastBlock = mFirstBlock + FirstChunkInAlloc->DataSize;
      
    } // ENDWHILE

    /* Step 1: Find the MemoryChunk number of sMemoryBlock.
       Step 2: Go to the address of that MemoryChunk */

    ChunkNumber = (size_t) (mem - mFirstBlock) / MemoryChunkSize;
    MemoryChunk* TempChunk = (MemoryChunk*)
      (mLastBlock + ChunkNumber * sizeof(MemoryChunk));

#ifdef OLD_WAY
    MemoryChunk* TempChunk = FirstChunk;
    while (TempChunk) {
      if (TempChunk->Data == ((TByte*) sMemoryBlock))
	break;
      TempChunk = TempChunk->Next;
    }
#endif /* OLD_WAY */

    return TempChunk;
  }

  /**********************************************************************/

  void MemoryPool::FreeAllAllocatedMemory(void)
  {
    MemoryChunk* Chunk = FirstChunk;
    void* DataToDelete = NULL;
    while (Chunk) {
      if (Chunk->IsAllocationChunk) {
	if (DataToDelete != NULL)
	  free(DataToDelete);
	DataToDelete = (void*) Chunk->Data;
      }
      Chunk = Chunk->Next;
    } // ENDWHILE
    if (DataToDelete != NULL)
      free(DataToDelete);
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
