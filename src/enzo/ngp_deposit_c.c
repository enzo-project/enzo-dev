#ifdef CONFIG_PFLOAT_16

#define ENZO_PYTHON_IMPORTED
#include "macros_and_parameters.h"
#undef ENZO_PYTHON_IMPORTED

// 10/16/2013 Sam Skillman
// Converted from ngp_deposit.F by hand.

int ngp_deposit_c(FLOAT *posx, FLOAT *posy, FLOAT *posz,
    int *ndim, int *npositions, float *mass, float *field,
    FLOAT *leftedge, int *dim1, int *dim2, int *dim3,
    float *cellsize)
{
    /* Local variables */
    int index;

    static int n, i1, j1, k1;
    const FLOAT fact = ((FLOAT) 1.0) / *cellsize;
    const FLOAT half = 0.5001;
    static FLOAT xpos, ypos, zpos;
    static FLOAT edge1, edge2, edge3;

    edge1 = (FLOAT) (*dim1) - half;
    edge2 = (FLOAT) (*dim2) - half;
    edge3 = (FLOAT) (*dim3) - half;

    int field_dim1, field_dim2, field_offset;
    /* Parameter adjustments */
    --mass;
    --posz;
    --posy;
    --posx;
    --leftedge;
    field_dim1 = *dim1;
    field_dim2 = *dim2;
    field_offset = 1 + field_dim1 * (1 + field_dim2);
    field -= field_offset;

    int np = *npositions;
    if (*ndim == 1){
        for (n=0; n<np; n++){
            xpos =  min(max((posx[n] - leftedge[0])*fact, half), edge1);

            i1 = (int) (xpos + ((FLOAT) 0.5));

            index = i1 + j1*(*dim1);
            field[index] += mass[n];
        }
    }

    if (*ndim == 2){
        for (n=0; n<np; n++){
            xpos =  min(max((posx[n] - leftedge[0])*fact, half), edge1);
            ypos =  min(max((posy[n] - leftedge[1])*fact, half), edge2);

            i1 = (int) (xpos + 0.5);
            j1 = (int) (ypos + 0.5);

            index = i1 + (j1 + k1*field_dim2)*field_dim1;
            field[index] += mass[n];
        }
    }

    if (*ndim == 3){
        for (n=0; n<np; n++){
            xpos =  min(max((posx[n] - leftedge[0])*fact, half), edge1);
            ypos =  min(max((posy[n] - leftedge[1])*fact, half), edge2);
            zpos =  min(max((posz[n] - leftedge[2])*fact, half), edge3);

            i1 = (int) (xpos + 0.5);
            j1 = (int) (ypos + 0.5);
            k1 = (int) (zpos + 0.5);

            index = i1 + (j1 + k1*field_dim2)*field_dim1;
            field[index] += mass[n];
        }
    }

    return 0;
}


#endif
