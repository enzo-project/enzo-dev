#ifndef __FOF_PROTO_H
#define __FOF_PROTO_H

void   allocate_memory(FOFData &AllVars);
void   check_cell(FOFData &AllVars, int p, int i, int j, int k);
void   compile_group_catalogue(FOFData &AllVars);

Eint32 comp_func(void const *a, void const *b);
Eint32 comp_func2(void const *a, void const *b);
Eint32 comp_func_partcoord(void const *a, void const *b);
Eint32 comp_func_gr(void const *a, void const *b);
Eint32 comp_func_partminid(void const *a, void const *b);
Eint32 compare_slab(const void *a, const void *b);

int    coarse_binning(FOFData &AllVars);
void   deallocate_all_memory(FOFData &D);
void   density(FOFData &A);
int    do_subfind_in_group(FOFData &D, FOF_particle_data *pbuf, int grlen, 
			   int *sublen, int *suboffset);
void   exchange_shadow(FOFData &AllVars, int TopGridResolution, bool SmoothData);
void   find_groups(FOFData &AllVars);
void   find_minids(FOFData &AllVars);
void   find_subgroups(FOFData &D);
int    get_particles(int dest, int minid, int len, FOF_particle_data *buf, 
		     FOFData &AllVars);
void   get_properties(FOFData D, FOF_particle_data *p, int len, bool subgroup, float *pcm, 
		      float *pcmv, float *pmtot, float *pmstars, float *pmvir,
		      float *prvir, float *pL, float *pvrms, float *pspin);
void   iindexx(int n, int arr[], int indx[]);
void   indexx(int n, float arr[], int indx[]);
void   init_coarse_grid(FOFData &AllVars);
int    link_across(FOFData &AllVars);
void   linkit(int p, int s, FOFData &AllVars);
void   link_local_slab(FOFData &AllVars);
void   marking(FOFData &AllVars);
int    number_of_unbound(FOFData &D, int head, int len);
void   order_subgroups_by_potential(FOFData &D);
double FOF_periodic_wrap(double x, double boxsize);
double FOF_periodic(double x, double boxsize);
void   save_groups(FOFData &AllVars, int CycleNumber, FLOAT EnzoTime);
float  selectb(unsigned long k, unsigned long n, float arr[], int ind[]);
void   sort_int(unsigned long n, int arr[]);
void   sort2_flt_int(unsigned long n, float arr[], int brr[]);
void   sort2_int(unsigned long n, int arr[], int brr[]);
void   set_sph_kernel(FOFData &A);
void   set_units(FOFData &AllVars);
void   stitch_together(FOFData &AllVars);
void   subfind(FOFData &D, int CycleNumber, FLOAT EnzoTime);
int    unbind(FOFData &D, int head, int len);
void   unbind_node(FOFData &D, int k);
void   walk_tree_and_unbind(FOFData &D);

#endif
