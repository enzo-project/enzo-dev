#ifdef DEFINE_STORAGE
# define EXTERN
#else /* DEFINE_STORAGE */
# define EXTERN extern
#endif

EXTERN int flow_trace_on;
EXTERN int flow_trace_level;
EXTERN FILE *flow_trace_fptr;


void flow_trace1( const char *name );
void flow_trace2( const char *name );
void print_flow_trace( const char *io, const char *name );

