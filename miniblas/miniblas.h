#ifndef MINIBLAS_H
#define MINIBLAS_H

#include "minivector_double.h"
#include "minimatrix_double.h"
#include <math.h>

enum MINIBLAS_TRANS {blasNoTrans=11, blasTrans=12};
enum MINIBLAS_UPORLOWER {blasUpper=21, blasLower=22};
enum MINIBLAS_UNIT {blasNonUnit=31, blasUnit=32};

/*
 * Level 1
 */

int miniblas_ddot (const minivector& X,
                    const minivector& Y,
                    double * result);

double miniblas_vector_ddot(const minivector& X,
                       const minivector& Y);


double miniblas_dnrm2  (const minivector& X);



int  miniblas_daxpy (double alpha,
                     const minivector& X,
                     minivector * Y);

void miniblas_dscal  (double alpha, minivector * X);

/*
 * Level 2
 */

int  miniblas_dgemv (MINIBLAS_TRANS TransA,
                      double alpha,
                      const minimatrix& A,
                      const minivector& X,
                      double beta,
                      minivector* Y);


int  miniblas_dtrsv (MINIBLAS_UPORLOWER Uplo,
                     MINIBLAS_TRANS TransA, MINIBLAS_UNIT Diag,
                     const minimatrix& A,
                     minivector * X);


int  miniblas_dger (double alpha,
                    const minivector& X,
                    const minivector& Y,
                    minimatrix * A);

int  miniblas_dsyr (MINIBLAS_UPORLOWER Uplo,
                    double alpha,
                    const minivector& X,
                    minimatrix * A);


/*
 * level 3
 */


int  miniblas_dgemm (MINIBLAS_TRANS TransA,
                      MINIBLAS_TRANS TransB,
                      double alpha,
                      const minimatrix& A,
                      const minimatrix& B,
                      double beta,
                      minimatrix * C);

int  miniblas_dsyrk (MINIBLAS_UPORLOWER Uplo,
                     MINIBLAS_TRANS Trans,
                     double alpha,
                     const minimatrix& A,
                     double beta,
                     minimatrix * C);


#endif
