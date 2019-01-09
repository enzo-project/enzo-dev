#ifdef TASKMAP

#ifdef USE_LL

#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include "llapi.h"

#include "macros_and_parameters.h"

 
Eint64 FreeRealMem( char *node )
{

  LL_element *queryObject;
  LL_element *job, *step;

  Eint32 rc, obj_count, err_code, ok;

  Eint64 free;
  Eint64 freemem;
  Eint64 def_mem;

  def_mem = (Eint64)(12.0); /* minimum memory on DataStar */

  Eint32 node_number;

#ifdef SEABORG
    sscanf(node, "s%5d", &node_number);
    fprintf(stderr, "Seaborg Node %s number %d\n", node, node_number);

    if ( 701 <= node_number && node_number <= 713 ) {
      def_mem = (Eint64)(56.0);
    } else {

      if ( 3701 <= node_number && node_number <= 5213 ) {
        def_mem = (Eint64)(28.0);
      } else {
        def_mem = (Eint64)(14.0);
      }

    }
#endif

#ifdef DATASTAR
    sscanf(node, "ds%3d", &node_number);
    fprintf(stderr, "DataStar Node %s number %d\n", node, node_number);

    if ( 4 <= node_number && node_number <= 11 ) {
      def_mem = (Eint64)(120.0);
    } else {

      if ( 300 <= node_number && node_number <= 395 ) {
        def_mem = (Eint64)(27.0);
      } else {
        def_mem = (Eint64)(12.0);
      }

    }
#endif

  def_mem = def_mem * 1024.0;
  fprintf(stderr,"Default node memory on node %s [%d] is %"ISYM" MBytes\n", node, node_number, def_mem);

  ok = 0;

  /* Initialize the query */

  queryObject = ll_query(MACHINES); 

  if (!queryObject) {
    fprintf(stderr, "Query MACHINES: ll_query() returns NULL.\n");
    freemem = def_mem;
    return(freemem);
  }

  char **nodelist;
  nodelist = (char **) malloc(2*sizeof(char *));
  char nodename[6];
  strcpy(nodename, node);
 
  nodelist[0] = nodename;
  nodelist[1] = NULL; 

  rc = ll_set_request(queryObject, QUERY_HOST, nodelist, ALL_DATA);

  if (rc) {
    fprintf(stderr, "Query MACHINES: ll_set_request() return code is non-zero.\n");
    freemem = def_mem;
    return(freemem);
  }

  job = ll_get_objs(queryObject, LL_CM, NULL, &obj_count, &err_code);
 
  /* Process the job objects */

  rc = ll_get_data(job, LL_MachineFreeRealMemory64, &free);

  if (rc) {
    fprintf(stderr, "Get data: ll_get_data() return code is non-zero: %d\n", rc);
    freemem = def_mem;
    return(freemem);
  }

/*
  fprintf(stderr, "Free: %lld MBytes\n", free);
*/

  ll_free_objs(queryObject);
  ll_deallocate(queryObject);

  freemem = (Eint64)(free);

  return(freemem);
}

#else

#include <stdlib.h>
#include <stdio.h>

#include "macros_and_parameters.h"

Eint64 FreeRealMem( char *node )
{

  Eint64 freemem;
  Eint64 def_mem;

  def_mem = (Eint64)(12.5); /* minimum memory of DataStar */

  Eint32 node_number;


#ifdef SEABORG
    sscanf(node, "s%5d", &node_number);
    fprintf(stderr, "Seaborg Node %s number %d\n", node, node_number);

    if ( 701 <= node_number && node_number <= 713 ) {
      def_mem = (Eint64)(56.0);
    } else {

      if ( 3701 <= node_number && node_number <= 5213 ) {
        def_mem = (Eint64)(28.0);
      } else {
        def_mem = (Eint64)(14.0);
      }

    }
#endif

#ifdef DATASTAR
    sscanf(node, "ds%3d", &node_number);
    fprintf(stderr, "DataStarNode %s number %d\n", node, node_number);

    if ( 4 <= node_number && node_number <= 11 ) {
      def_mem = (Eint64)(120.0);
    } else {

      if ( 300 <= node_number && node_number <= 395 ) {
        def_mem = (Eint64)(27.0);
      } else {
        def_mem = (Eint64)(12.0);
      }

    }
#endif

  def_mem = def_mem * 1024.0;
  fprintf(stderr,"Default node memory on node %s [%d] is %"ISYM" MBytes\n", node, node_number, def_mem);

  freemem = def_mem;

  return(freemem);
}

#endif

#else

#include <stdlib.h>
#include <stdio.h>

#include "macros_and_parameters.h"

Eint64 FreeRealMem( char *node )
{

  Eint64 freemem;

  fprintf(stderr, "Error in FreeRealMem - the TASKMAP macro is not defined! Default to 12.0 GBytes.\n");

  freemem = (Eint64)(12.0*1024.0); /* minimum memory on DataStar */

  return(freemem);
}

#endif
