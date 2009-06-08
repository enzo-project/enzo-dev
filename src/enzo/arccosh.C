#include <stdio.h>
#include <math.h>
 
#ifdef CRAYX1
 
double arccosh(double x)
{
  if( x == 1.0 ) {
    return (0.0);
  }
  if( x < 1.0 ) {
    printf("ArcCosh argument < 1; x = %16.9e\n", x);
    return (0.0);
  }
 
  return (log(x+sqrt(x*x-1.0)));
}
 
#else
 
double arccosh(double x)
{
  return (acosh(x));
}
 
#endif
