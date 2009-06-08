/***********************************************************************
/
/  GLOBAL DATA DECLARATIONS
/
/  written by: Greg Bryan
/  date:       June, 1997
/  modified1:
/
/  PURPOSE:
/    This is the global data, which should be held to a minimum.
/
************************************************************************/

#ifdef DEFINE_STORAGE
# define EXTERN
#else /* DEFINE_STORAGE */
# define EXTERN extern
#endif

/* debugging flag */

EXTERN int debug;
