#include <stdio.h>
#include <time.h>
#include <sys/times.h>
#include <sys/time.h>

#if defined(IRIS4) || defined(CONVEX) || defined(COMPAQ) || defined(SUN) || defined(LINUX) || defined(IA64) || defined(CRAYX1) || defined(XT3)
double wall_clock_(void);

double wall_clock_(void)
#endif
#if defined(SP2) || defined(HP)
double wall_clock(void);

double wall_clock(void)
#endif
{
  struct timeval timestr;
  void *Tzp=0;
  gettimeofday(&timestr, Tzp);
  return (double)timestr.tv_sec+1.0E-06*(double)timestr.tv_usec;
}
