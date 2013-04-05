void WriteSingleCube(float * array, int Dims[], char* string, int dNum, int gNum, char * label);

int TracerParticlesAddToRestart_DoIt(char * filename, HierarchyEntry *TopGrid,
				    TopGridData *MetaData);

#ifdef DEFINE_STORAGE
# define EXTERN
#else /* DEFINE_STORAGE */
# define EXTERN extern
#endif


