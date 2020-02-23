#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

#include <vector>
#include <map>
#include <list>
#include "../miniblas/minimatrix_double.h"
#include "../miniblas/minivector_double.h"
#include "../miniblas/minilinalg.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2 M_PI * 0.5
#endif

namespace minisam
{
///////////////////////
void zeroBelowDiagonal(minimatrix *A, size_t cols = 0);
double norm2d(const minivector *v1);
double norm2d(const minivector& v1);
double norm2d(const minivector* p, minimatrix* H);
double norm2d(const minivector& p, minimatrix* H);

bool assert_equal(const minivector *expected, const minivector *actual, double tol);
bool assert_equal(const minivector& expected, const minivector& actual, double tol);

bool equal_with_abs_tol(const minivector *vec1, const minivector *vec2, double tol);
bool equal_with_abs_tol(const minivector& vec1, const minivector& vec2, double tol);

void vector_scale_inplace(const minivector& v, minimatrix& A, bool inf_mask = false);

minimatrix vector_scale(const minivector& v, const minimatrix&A, bool inf_mask = false);       // row
minimatrix *vector_scale(const minimatrix *A, const minivector *v, bool inf_mask = false); // column

double distance2(const minivector& v1,const minivector& v2);

inline minimatrix skewSymmetric_object(double wx, double wy, double wz)
{

    minimatrix m3d(3,3);
    m3d.data[0]=0.0;
    m3d.data[1]=-wz;
    m3d.data[2]=wy;
    m3d.data[3]=wz;
    m3d.data[4]=0.0;
    m3d.data[5]=-wx;
    m3d.data[6]=-wy;
    m3d.data[7]=wx;
    m3d.data[8]=0.0;


    return m3d;
}
inline minimatrix* skewSymmetric(double wx, double wy, double wz)
{
    minimatrix* m3d=new minimatrix(3,3);
    m3d->data[0]=0.0;
    m3d->data[1]=-wz;
    m3d->data[2]=wy;
    m3d->data[3]=wz;
    m3d->data[4]=0.0;
    m3d->data[5]=-wx;
    m3d->data[6]=-wy;
    m3d->data[7]=wx;
    m3d->data[8]=0.0;


    return m3d;
}


inline void skewSymmetric(minimatrix* m3d,double wx, double wy, double wz)
{

    if(m3d->size1!=3||m3d->size2!=3)
    {
        minimatrix_resize(m3d,3,3);
    }

    m3d->data[0]=0.0;
    m3d->data[1]=-wz;
    m3d->data[2]=wy;
    m3d->data[3]=wz;
    m3d->data[4]=0.0;
    m3d->data[5]=-wx;
    m3d->data[6]=-wy;
    m3d->data[7]=wx;
    m3d->data[8]=0.0;
}

inline minimatrix*  skewSymmetric(const minivector* w)
{
    return skewSymmetric(w->data[0], w->data[1], w->data[2]);
}
inline void skewSymmetric(minimatrix* m3d,const minivector& w)
{
    return skewSymmetric(m3d,w.data[0], w.data[1], w.data[2]);
}


bool equal_with_abs_tol(const minimatrix *A, const minimatrix *B, double tol = 1e-9);
bool equal_with_abs_tol(const minimatrix& A, const minimatrix& B, double tol = 1e-9);

void backSubstituteUpper(const minimatrix& U, const minivector& b, minivector *x);

void inplace_QR(minimatrix* A);


std::pair<int, bool> choleskyCareful(minimatrix*ATA, int order = -1);
bool choleskyPartial(minimatrix *ABC, int nFrontal, int topleft = 0);

};     // namespace minisam
#endif // MATRIX_H_INCLUDED
