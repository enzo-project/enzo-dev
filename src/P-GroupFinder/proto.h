#ifdef USE_HDF4
#include "hdf.h"
#endif
#ifdef USE_HDF5
#include "hdf5.h"
#endif

void   allocate_memory(void);
void   check_cell(int p, int i, int j, int k);
void   compile_group_catalogue(void);
int    comp_func(void const *a, void const *b);
int    comp_func2(void const *a, void const *b);
int    comp_func_partcoord(void const *a, void const *b);
int    comp_func_gr(void const *a, void const *b);
int    comp_func_partminid(void const *a, void const *b);
void   count_local_particles(char *fname, int files);
int    course_binning(void);
void   density(void);
int    do_subfind_in_group(struct particle_data *pbuf, int grlen, int *sublen, int *suboffset);
void   exchange_shadow(void);
int    find_files(char *fname);
void   find_groups(void);
void   find_minids(void);
void   find_subgroups(void);
int    get_particles(int dest, int minid, int len, struct particle_data *buf);
void get_properties(struct particle_data *p, int len, float *pcm, float *pmtot, 
		    float *pmgas, float *pmstars, float *pmsfr, float *pmcold,
		    int subgroup, float *pcmv, float *pmvir, float *prvir,
		    float *pL, float *pvrms, float *pspin);
int    get_slab(int index);
void   iindexx(unsigned int n, int arr[], unsigned int indx[]);
void   indexx(unsigned long n, float arr[], int indx[]);
void   init_coarse_grid(void);
int    link_accross(void);
void   linkit(int p, int s);
void   link_local_slab(void);
void   loadpositions(char *fname, int files);
void   marking(void);
void   *mymalloc(size_t size);
int    number_of_unbound(int head, int len);
void   order_subgroups_by_potential(void);
double periodic_wrap(double x);
double periodic(double x);
void save_groups(char *particles_fname, char *particles_fname5, 
		 char *catalogue_fname, char *parttypes_fname, 
		 char *partids_fname, char *cataloguetxt);
float  selectb(unsigned long k, unsigned long n, float arr[], int ind[]);
void   sort2_flt_int(unsigned long n, float arr[], int brr[]);
void   sort2_int(unsigned long n, int arr[], int brr[]);
void   set_sph_kernel(void);
void   stitch_together(void);
void subfind(char *particles_fname, char *catalogue_fname, 
	     char *subhalo_fname, char *parttypes_fname, char *partids_fname, 
	     char *subprop_fname, char *prop_fname, char *sparticles_fname5,
	     char *scataloguetxt);
int    unbind(int head, int len);
void   unbind_node(int k);
void   walk_tree_and_unbind(void);

void enzoLoadPositions (char *fname, int files);
int enzoFindFiles (char *fname);
void enzoCountLocalParticles (char *fname, int files);
void MarkInteriorParticles (char *fname, int files);

#ifdef USE_HDF4
void ReadParticleField_FLOAT64 (int32 sd_id, char *label, int nPart, 
				double **data);
void ReadParticleField_FLOAT32 (int32 sd_id, char *label, int nPart, 
				float **data);
void ReadParticleField_INT (int32 sd_id, char *label, int nPart, int **data);
#endif

#ifdef USE_HDF5
void ReadParticleFieldHDF5_INT(hid_t group_id, char *label, int nPart, PINT **data);
void ReadParticleFieldHDF5_FLOAT(hid_t group_id, char *label, int nPart, float **data);
void ReadParticleFieldHDF5_DOUBLE(hid_t group_id, char *label, int nPart, double **data);
#endif

int compare_slab(const void *a, const void *b);
void write_ascii_catalog(char *catalogue_fname, char *fofprop_fname, 
			 char *cataloguetxt);
