#include<stdio.h>
#include<ctype.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"

Eint32 hide_isdigit(Eint32 c)
{
  return (c >= '0' && c <= '9');
}
