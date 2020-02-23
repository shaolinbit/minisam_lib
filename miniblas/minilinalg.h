#ifndef __MINILINALG_H__
#define __MINILINALG_H__

#include <stdlib.h>
#include "minivector_double.h"
#include "minimatrix_double.h"
#include "miniblas.h"
#include <stdexcept>



// LU Decomposition

void minilinalg_LU_decompmat(minimatrix* A,const minimatrix& source);
double minilinalg_LU_det (minimatrix * LU, int signum);
void minilinalg_luMatInverse(minimatrix *ainverse,minimatrix *a);
void nr_ludcmp(minimatrix* a,int* indx, int *d);
void nr_lubksb(minimatrix* a,int* indx, double* b);

/* QR decomposition */

int minilinalg_Golub_QR_decomp (minimatrix * A,
                          minivector * tau);

int minilinalg_R_solve (const minimatrix& R,
                        const minivector& b,
                        minivector * x);

/* Cholesky Decomposition */
int minilinalg_nr_cholesky_decomp(minimatrix * A,bool lowerorupper=true);





#endif /* __MINILINALG_H__ */
