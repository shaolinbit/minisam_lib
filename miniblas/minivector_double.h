#ifndef __MINIVECTOR_DOUBLE_H__
#define __MINIVECTOR_DOUBLE_H__

#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iosfwd>
#include "minimatrix_double.h"
#include "memorypool.h"

#define miniblas_max(a,b) ((a) > (b) ? (a) : (b))
#define miniblas_min(a,b) ((a) < (b) ? (a) : (b))

struct minivector:public minimatrix
{
    minivector()
    {

    }
    minivector(int n):minimatrix(n,1)
    {

    }
    minivector(const minivector& mcopy):minimatrix(mcopy.size1,1)
    {
        dimension=mcopy.dimension;
        const size_t src_prd = mcopy.prd ;

        for (size_t j = 0; j < this->size1; j++)
        {
            this->data[ j]
                = mcopy.data[ src_prd * j];
        }

    }
    minivector(const minimatrix& mcopy):minimatrix(mcopy.size1,1)
    {
        dimension=mcopy.dimension;
        if(mcopy.size2!=1)
        {
            throw std::invalid_argument("vectors must have same length");
        }
        const size_t src_prd = mcopy.prd ;
        for (size_t j = 0; j < this->size1; j++)
        {

            this->data[ j]
                = mcopy.data[ src_prd * j];
        }

    }
    minivector(const minimatrix* m_memory)
    {
        if(m_memory->size2!=1)
        {
           throw std::invalid_argument("minimatrix is not a vector");
        }
        size1=m_memory->size1;
        size2=1;
        prd=m_memory->prd;
        data=m_memory->data;
        owner=0;
        dimension=m_memory->dimension;

    }
    minivector(const minivector* m_memory)
    {
        size1=m_memory->size1;
        size2=1;
        prd=m_memory->prd;
        data=m_memory->data;
        owner=0;
        dimension=m_memory->dimension;
    }

    minivector(int n,double value):minimatrix(n,1)
    {
        for (size_t j = 0; j < this->size1; j++)
        {

            this->data[j]
                = value;
        }
    }
    minivector(double x,double y,double z):minimatrix(3,1)
    {
        data[0]=x;
        data[1]=y;
        data[2]=z;
    }
    virtual minimatrix LocalCoordinates(const minimatrix* mpose,minimatrix* H1=NULL,minimatrix* H2=NULL) const
    {
        const size_t N = size1;
        minivector result(size1);

        if((mpose->size1 != N)||(mpose->size2!=1))
        {
            throw std::invalid_argument("vectors must have same length");
        }
        else
        {
            const size_t prd_b = mpose->prd;

            size_t i;

            for (i = 0; i < N; i++)
            {
                result.data[i * prd]=-data[i*prd] +mpose->data[i * prd_b];
            }

            return result;
        }
    }
    virtual minimatrix* Retract(const minimatrix* mpose)
    {
        minivector* result=new minivector(size1);
        const size_t M = size1;

        if(mpose->size1!=M||mpose->size2!=1)
        {
            throw std::invalid_argument("matrices must have same dimensions");
        }
        else
        {
            const size_t prd_b = mpose->prd;

            size_t i;

            for (i = 0; i < M; i++)
            {
                result->data[i * prd]=data[i*prd]+mpose->data[i * prd_b];
            }

            return result;
        }
    }
    virtual minimatrix between(const minimatrix* mpose) const
    {
        size_t i, j;
        minivector result(size1);
        for (i = 0; i < size1; i++)
        {
            result.data[i * prd ] = mpose->data[i * prd]-data[i * prd];
        }
        return result;
    }
    virtual minimatrix between(const minimatrix* mpose,minimatrix& H1,minimatrix& H2) const
    {

        const size_t p = dimension ;
        const size_t q = dimension ;
        size_t i, j;
        minimatrix_resize(&H1,p,q);
        minimatrix_set_neg_identity(&H1);
        minimatrix_resize(&H2,p,q);
        minimatrix_set_identity(&H2);

        minivector result(size1);
        for (i = 0; i < size1; i++)
        {
            result.data[i * prd ] = mpose->data[i * prd]-data[i * prd];
        }
        return result;
    }




};



minivector minivector_2dim(double x,double y);

void minivector_resize(minivector* m,const size_t n1);

minivector
minivector_subvector_var(minivector *v,
                     size_t i,
                     size_t n);
minivector
minivector_subvector(const minivector& v,
                     size_t i,
                     size_t n);

void minivector_set_zero (minivector * v);
void minivector_set_all (minivector * v, double x);
int minivector_set_basis (minivector * v, size_t i);
int minivector_memcpy (minivector * dest, const minivector * src);
int minivector_memcpy (minivector * dest, const minivector& src);
int minivector_memcpy_normalized(minivector * dest, const minivector& src);
int minivector_reverse (minivector * v);

int minivector_swap (minivector * v, minivector * w);
int minivector_swap_elements (minivector * v, const size_t i, const size_t j);
double minivector_lpNormInfinity(const minivector& v);
int minivector_vectorbiggerthanthreshold(const minivector& v1,const minivector& v2);
int minivector_add (minivector * a, const minivector * b);
int minivector_add (minivector * a, const minivector& b);
int minivector_add (minivector * a, const minivector& b, const minivector& c);
int minivector_add_by(minivector * a, double a1, const minivector& b);
int minivector_add_xby(minivector * a, const minivector& b,double a1, const minivector& c);
int minivector_axby(minivector * a, double a1,const minivector& b,double a2, const minivector& c);
int minivector_sub (minivector * a, const minivector * b);
int minivector_sub (minivector * a, const minivector& b);
int minivector_sub (minivector * a, const minivector& b, const minivector& c);int minivector_sub (minivector * a, const minivector* b, const minivector* c);
int minivector_mul (minivector * a, const minivector * b);
int minivector_div (minivector * a, const minivector * b);
int minivector_scale (minivector * a, const double x);
int minivector_scale_vec(minivector * a, const double x,const minivector& b);
minivector minivector_square(const minivector& b);
void minivector_cwisesqrt(minivector * a);
void minivector_cwisemax(minivector * a,double maxnumber);
void minivector_cwisemin(minivector * a,double minnumber);
minivector minivector_cwisesqrt(const minivector& a);
minivector minivector_cwiseProduct(const minivector& a,const minivector& b);
void minivector_cwiseProduct(minivector* result,const minivector& a,const minivector& b);
minivector minivector_reciprocal(const minivector& a);
void minivector_print(const minivector& m);
void minivector_ostream(std::ostream& os,const minivector& m);
int minivector_add_constant (minivector * a, const double x);

int minivector_equal (const minivector& u,
                      const minivector& v);
int minivector_fabsvectorbiggerthanthreshold(const minivector& v1,const minivector& v2);

int minivector_isnull (const minivector * v);
int minivector_isnull (const minivector& v);
int minivector_hasnan(const minivector& v);

int minivector_ispos (const minivector * v);
int minivector_isneg (const minivector * v);
int minivector_isnonneg (const minivector * v);

double minivector_get (const minivector * v, const size_t i);
void minivector_set (minivector * v, const size_t i, double x);

minivector
minimatrix_row (minimatrix * m, const size_t i);

minivector minimatrix_row (const minimatrix& m, const size_t i);

minivector
minimatrix_column (minimatrix * m, const size_t j);

minivector
minimatrix_column (const minimatrix& m, const size_t j);

minivector
minimatrix_diagonal (const minimatrix& m);
int minimatrix_set_diagonal(minimatrix* a,const minivector& b);


minivector
minimatrix_subrow (minimatrix * m, const size_t i,
                   const size_t offset, const size_t n);

minivector
minimatrix_subcolumn (minimatrix * m, const size_t j,
                      const size_t offset, const size_t n);
minivector
minimatrix_subcolumn (const minimatrix& m, const size_t j,
                      const size_t offset, const size_t n);

minimatrix minimatrix_vector_to_matrix(const minivector& a);



int minimatrix_get_row(minivector * v, const minimatrix * m, const size_t i);
int minimatrix_get_col(minivector * v, const minimatrix * m, const size_t j);
int minimatrix_set_row(minimatrix * m, const size_t i, const minivector * v);
int minimatrix_set_col(minimatrix * m, const size_t j, const minivector * v);
void QuaternionToMatrix(const minivector& qnb,minimatrix* mat);
void QuaternionToMatrixTrans(const minivector& qnb,minimatrix* mat);
void Quaternion_Multiply(minivector* qqQ,const minivector &gQ, const minivector & hQ);

void Quaternion_inv_Multiply(minivector* qqQ,const minivector &gQ, const minivector & hQ);

void AngleAxisd_to_Quaternion(minivector* Q,const minivector& angleaxidandangle);

void AngleAxisangle_to_Quaternion(minivector* Q,const minivector& angleaxis,double angle);

minimatrix minimatrix_vector_asDiagonal(const minivector& v);
#endif /* __MINIVECTOR_DOUBLE_H__ */


