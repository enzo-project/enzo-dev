#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "phys_constants.h"
#include "Star.h"

/* Handle *when* an SN occurs in CheckTimeAction:

   3 actions:
     1) SNIa at constant rate (easiest to do)

     2) Isolated, single SNII at variable rate:
        - Fix initial rate given cold gas mass to get SNR and SNR_individual
        - After every SN explosion, recompute cold gas mass
          and determine next explosion time (store SNR as global)
             T_next = T_now + 1 / SNR_individual

          where SNR_individual = SNR * (1.0 - f_cluster)

     3) Some fraction of SNII will be clustered (most?):
        - Do constant clustering size, say n_cluster = 8 (pick reasonable number for avg cluster size)
              - argue this is to remove stochasicity and improve reliablity
        - Fix initial rate given cold gas mass (f_clust / n_cluster * SNR)
        - After every clustered SN explosion, recompute when the next one
          goes off: T_next = T_now + 1 / SNR_cluster with
          SNR_cluster = (f_clust / n_cluster) * SNR

     4) Pick random positions for all SN as some position in initial galaxy disk size. Bias
        SNII towards center, keep SNIa at initial disk without bias.
*/


/* function prototypes */
int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);

