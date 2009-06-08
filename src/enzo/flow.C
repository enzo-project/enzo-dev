#include <stdio.h>
 
#include "flowdefs.h"
 
 
void flow_trace1( const char *name )
{
  void print_flow_trace( const char *io, const char *name );
 
  flow_trace_level = flow_trace_level + 1;
  print_flow_trace("> ", name);
}
 
void flow_trace2( const char *name )
{
  void print_flow_trace( const char *io, const char *name );
 
  print_flow_trace("< ", name);
  flow_trace_level = flow_trace_level - 1;
}
 
 
void print_flow_trace( const char *io, const char *name )
{
 
/*
char line[132];
const char pad[2] = ".";
int ll;
 
  if (strlen(name)+strlen(io)+flow_trace_level < 132)
  {
 
    line[0] = NULL;
 
    for (ll = 0; ll < flow_trace_level; ll++)
    {
      strcat(line, pad);
    }
 
    strcat(line, io);
    strcat(line, name);
 
    fprintf(flow_trace_fptr, "%8d : %s\n", flow_trace_level, line);
  }
*/
 
fprintf(flow_trace_fptr, "%8d %s\n", flow_trace_level, name);
 
}
