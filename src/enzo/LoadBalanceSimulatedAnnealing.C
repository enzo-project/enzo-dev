/***********************************************************************
/
/  COMMUNICATION ROUTINE: LOAD BALANCE BY SIMULATED ANNEALING
/
/  written by: John Wise
/  date:       July, 2009
/  modified1:
/
/  NOTES: Based on Chapter 11 of Parallel Computing Works! (Fox et al.)
/         http://www.netlib.org/utk/lsi/pcwLSI/text/
/
/         In each iteration, a random grid is moved to a random processor.  
/         If this improves the balance (quantified by the cost function, H), 
/         the change is kept.  If not, there is a probability exp(-dH/T) 
/         the change is kept, where dH is the change in the cost function 
/         and T is the current "temperature".  The temperature is steadily
/         decreased in the process.
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "TopGridData.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "communication.h"
#include "CommunicationUtilities.h"

void fpcol(Eflt64 *x, int n, int m, FILE *log_fptr);

int anneal(int ngrids, int nnodes, int* &proc, double* cells, double* scells, 
	   float temperature, float node_weight, float same_node_prob);
double TotalCostFunction(int ngrids, int nnodes, int* proc, double* cells, 
			 double* scells, float node_weight);
double ProcessorCostFunction(int ngrids, int* proc, double* cells);
double NodeCostFunction(int ngrids, int nnodes, int* proc, double* scells, 
			float weight);
int DetermineNodeNumber(int pp, int nnodes);

#define NO_LB_DEBUG
//#define LB_DEBUG

int LoadBalanceSimulatedAnnealing(int NumberOfGrids, int NumberOfNodes, 
				  int* &ProcessorNumbers, double* NumberOfCells, 
				  double* NumberOfSubcells)
{

  /* Tunable parameters */

  // The higher the initial temperature, the more the initial
  // configuration is forgotten.  1e-3 is a good value because we want
  // to keep the grid transfers to a minimum.  This parameter is
  // expressed in the fraction of the initial cost function.
  const float initial_temperature = 1e-3;

  // Preference to which the nodes are balanced with the number of
  // subgrid cells over the balance of the cells on processors.  4 is
  // a good value.
  const float node_balance_weight = 1.0;

  // In the random choice of a new processor, this is the probability
  // that this processor will be on the same node as before.
  // Currently I haven't tested how good this parameter works, so I'm
  // setting it to zero.
  const float same_node_prob = 0.0;

  // Number of iterations (i.e. temperature reductions)
  const int NumberOfIterations = 100;

  // Number of attempted moves at a single temperature
  const int NumberOfTries = NumberOfGrids;

  // Maximum number of moves at a temperature
  int MoveLimit = NumberOfTries / 5;

  // Factor to reduce the temperature each iteration
  const float TemperatureReduce = 0.95;

  /************************************************************************/

  if (LoadBalancing < 2 || LoadBalancing > 3)
    ENZO_FAIL("Calling LBSimulatedAnnealing with a wrong LoadBalancing setting.");

  if (rand_init == 0) {
    srand( time(NULL) );
    rand_init = 1;
  }

  int i, j, node, nsucc;
  float temperature, H;
  double *ProcWork, *NodeWork;

  ProcWork = new double[NumberOfProcessors];
  NodeWork = new double[NumberOfNodes];

  for (i = 0; i < NumberOfProcessors; i++)
    ProcWork[i] = 0;
  for (i = 0; i < NumberOfNodes; i++)
    NodeWork[i] = 0;

  for (i = 0; i < NumberOfGrids; i++) {
    node = DetermineNodeNumber(ProcessorNumbers[i], NumberOfNodes);
    ProcWork[ProcessorNumbers[i]] += NumberOfCells[i];
    NodeWork[node] += NumberOfSubcells[i];
  }

  if (debug) {
#ifdef LB_DEBUG
    printf("BEFORE LB proc:");
    fpcol(ProcWork, NumberOfProcessors, 16, stdout);
#endif /* LB_DEBUG */
    printf("BEFORE LB nodes:");
    fpcol(NodeWork, NumberOfNodes, 16, stdout);
  }

  temperature = initial_temperature *
    TotalCostFunction(NumberOfGrids, NumberOfNodes, ProcessorNumbers, NumberOfCells,
		      NumberOfSubcells, node_balance_weight);
     
  for (i = 0; i < NumberOfIterations; i++) {
    nsucc = 0;
    for (j = 0; j < NumberOfTries; j++) {

      if (anneal(NumberOfGrids, NumberOfNodes, ProcessorNumbers, NumberOfCells, 
		 NumberOfSubcells, temperature, node_balance_weight, 
		 same_node_prob) == TRUE)
	nsucc++;
      if (nsucc > MoveLimit)
	break;

    } // ENDFOR j

#ifdef LB_DEBUG
    H = TotalCostFunction(NumberOfGrids, NumberOfNodes, ProcessorNumbers, 
			  NumberOfCells, NumberOfSubcells, node_balance_weight); 

    if (debug)
      printf("iter %3d / %3d: T = %8.3g, H = %8.4g, nsucc = %d\n",
	     i, j, temperature, H, nsucc);
#endif /* LB_DEBUG */

    temperature *= TemperatureReduce;

  } // ENDFOR i
  
  for (i = 0; i < NumberOfProcessors; i++)
    ProcWork[i] = 0;
  for (i = 0; i < NumberOfNodes; i++)
    NodeWork[i] = 0;

  for (i = 0; i < NumberOfGrids; i++) {
    node = DetermineNodeNumber(ProcessorNumbers[i], NumberOfNodes);
    ProcWork[ProcessorNumbers[i]] += NumberOfCells[i];
    NodeWork[node] += NumberOfSubcells[i];
  }

  if (debug) {
#ifdef LB_DEBUG
    printf("AFTER LB proc:");
    fpcol(ProcWork, NumberOfProcessors, 16, stdout);
#endif /* LB_DEBUG */
    printf("AFTER LB nodes:");
    fpcol(NodeWork, NumberOfNodes, 16, stdout);
  }

  delete [] ProcWork;
  delete [] NodeWork;

  return SUCCESS;
}

/************************************************************************/

int anneal(int ngrids, int nnodes, int* &proc, double* cells, double* scells, 
	   float temperature, float node_weight, float same_node_prob)
{

  int change_grid, old_proc, new_proc, old_node, new_node;
  int random_proc_on_node;
  float prob;
  double H0, H1, dH;

  // Choose a random grid
  change_grid = rand() % ngrids;

  old_proc = proc[change_grid];
  old_node = DetermineNodeNumber(old_proc, nnodes);

  // Choose a random processor
  new_proc = old_proc;
  while (new_proc == old_proc) {
    new_proc = rand() % NumberOfProcessors;
    if (same_node_prob > 0) {
      if (LoadBalancing == 2)
	random_proc_on_node = (rand() % CoresPerNode) + old_node * CoresPerNode;
      else if (LoadBalancing == 3)
	random_proc_on_node = (rand() % CoresPerNode) * nnodes + old_node;
      if (float(rand()) / RAND_MAX > same_node_prob)
	new_proc = random_proc_on_node;
    } // ENDIF same_node_prob > 0
  } // ENDWHILE

  // Original cost function
  H0 = TotalCostFunction(ngrids, nnodes, proc, cells, scells, node_weight);

  // Move grid to new processor and evaulate cost function again
  proc[change_grid] = new_proc;
  H1 = TotalCostFunction(ngrids, nnodes, proc, cells, scells, node_weight);

  // If the change reduces the cost function, keep the change.  If
  // not, there is a probability of exp(-dH/T) that the change is
  // kept.
  dH = H1 - H0;
  if (dH < 0)
    prob = 1;
  else
    prob = exp(-dH / temperature);
  
  if (float(rand()) / RAND_MAX < prob) {
    // Accept and return true
    return TRUE;
  } else {
    // Reject, return to original state, and return false
    proc[change_grid] = old_proc;
    return FALSE;
  }

}

/************************************************************************/
/*                            COST FUNCTIONS                            */
/************************************************************************/

double TotalCostFunction(int ngrids, int nnodes, int* proc, double* cells, 
			 double* scells, float node_weight)
{
  return ProcessorCostFunction(ngrids, proc, cells) +
    NodeCostFunction(ngrids, nnodes, proc, scells, node_weight);
}

double ProcessorCostFunction(int ngrids, int* proc, double* cells)
{

  // Least squares of the work on each processor.

  int i;
  double factor, result = 0, TotalWork = 0;
  double *WorkPerProc = new double[NumberOfProcessors];

  for (i = 0; i < NumberOfProcessors; i++)
    WorkPerProc[i] = 0;

  for (i = 0; i < ngrids; i++) {
    TotalWork += cells[i];
    WorkPerProc[proc[i]] += cells[i];
  }
  
  for (i = 0; i < NumberOfProcessors; i++)
    result += WorkPerProc[i] * WorkPerProc[i];

  factor = NumberOfProcessors / TotalWork;
  result *= factor*factor;
  
  delete [] WorkPerProc;

  return result;

}

double NodeCostFunction(int ngrids, int nnodes, int* proc, double* scells, 
			float weight)
{

  int i, node;
  double factor, result = 0, TotalWork = 0;
  double *WorkPerNode = new double[nnodes];

  for (i = 0; i < nnodes; i++)
    WorkPerNode[i] = 0;

  for (i = 0; i < ngrids; i++) {
    TotalWork += scells[i];
    node = DetermineNodeNumber(proc[i], nnodes);
    WorkPerNode[node] += scells[i];
  }

  for (i = 0; i < nnodes; i++)
    result += WorkPerNode[i] * WorkPerNode[i];

  factor = nnodes / TotalWork;
  result *= weight * factor * factor;

  delete [] WorkPerNode;

  return result;

}

int DetermineNodeNumber(int pp, int nnodes)
{
  int result = 0;
  if (LoadBalancing == 2)
    result = pp / CoresPerNode;
  else if (LoadBalancing == 3)
    result = pp % nnodes;
  else
    ENZO_FAIL("DetermineNodeNumber isn't valid for this type of load balancing.");
  return result;
}
