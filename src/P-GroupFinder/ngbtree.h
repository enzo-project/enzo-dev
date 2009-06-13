
void ngb_treebuild(int Npart);

void ngb_treefree(void);
void ngb_treeallocate(int npart,int maxnodes);  /* usually maxnodes=2*npart is suffiecient */

float ngb_treefind(double xyz[3], int desngb, float hguess,int **ngblistback, float **r2listback);


float ngb_treetest(float xyz[3], int desngb, float hguess,int **ngblistback, float **r2listback);
