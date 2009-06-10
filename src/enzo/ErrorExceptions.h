/***********************************************************************
/
/  EXCEPTION CLASS
/
/  written by: Matthew Turk
/  date:       June, 2009
/  modified:   
/
/  PURPOSE:
/
************************************************************************/
#ifndef __ENZO_EXCEPTIONS__
#define __ENZO_EXCEPTIONS__

#include <execinfo.h>
#include <signal.h>
#include <stdio.h>

#include <exception>

// Example from 
//
// http://www.ibm.com/developerworks/linux/library/l-cppexcep.html
//

 class EnzoFatalException
 {
 public:
     EnzoFatalException(const char * error_msg)
     {
         void * array[25];
         int nSize = backtrace(array, 25);
         char ** symbols = backtrace_symbols(array, nSize);
         fprintf(stderr, "Caught fatal exception:\n\n");
         fprintf(stderr, "   '%s'\n\n", error_msg);
         fprintf(stderr, "Backtrace:\n\n");

         for (int i = 0; i < nSize; i++)
         {
             fprintf(stderr, "BT symbol: %s\n", symbols[i]);
         }

         delete [] symbols;
     }
    void WriteDebuggingOutput()
    {

    }
 };

#endif
