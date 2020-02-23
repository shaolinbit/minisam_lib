#ifndef MINIMATRIX_DOUBLE_H
#define MINIMATRIX_DOUBLE_H

#include <stdlib.h>
#include <iostream>
#include "memorypool.h"


struct minimatrix
{
    size_t size1;
    size_t size2;
    size_t prd;//physical row dimension
    size_t dimension;
    double * data;
    int owner;
    minimatrix():size1(0),size2(0),prd(0),data(NULL),owner(0),dimension(0)
    {

    }
    minimatrix(int m,int n):size1(m),size2(n),prd(n),owner(1),dimension(m*n)
    {
        data=(double *) malloc (m*n * sizeof (double));
    }
    minimatrix(const minimatrix& mcopy):size1(mcopy.size1),size2(mcopy.size2),
        prd(mcopy.size2),owner(1),dimension(mcopy.dimension)
    {
        const size_t src_size1 = mcopy.size1;
        const size_t src_size2 = mcopy.size2;
        data=(double *) malloc (src_size1*src_size2 * sizeof (double));

        const size_t src_prd = mcopy.prd ;
        size_t i, j;

        for (i = 0; i < src_size1 ; i++)
        {
            for (j = 0; j < src_size2; j++)
            {
                data[ prd * i + j]
                    = mcopy.data[src_prd * i + j];
            }
        }

    }
    minimatrix(minimatrix* m_memory)
    {
        size1=m_memory->size1;
        size2=m_memory->size2;
        prd=m_memory->prd;
        data=m_memory->data;
        owner=0;
        dimension=m_memory->dimension;
    }
    minimatrix(const minimatrix* m_memory)
    {
        size1=m_memory->size1;
        size2=m_memory->size2;
        prd=m_memory->prd;
        data=m_memory->data;
        owner=0;
        dimension=m_memory->dimension;
    }
    virtual double x() const
    {
        if((data!=NULL)&&(size1*size2>=1))
        {
            return data[0];
        }
    }
    virtual double y() const
    {
        if((data!=NULL)&&(size1*size2>1))
        {
            return data[1];
        }
    }
    virtual minimatrix LocalCoordinates(const minimatrix* mpose) const
    {
        minimatrix result(size1,size2);
        const size_t M = size1;
        const size_t N = size2;

        if(mpose->size1!=M||mpose->size2!=N)
        {
            throw std::invalid_argument("matrices must have same dimensions");
        }
        else
        {
            const size_t prd_a = result.prd;
            const size_t prd_b = this->prd;
            const size_t prd_c = mpose->prd;

            size_t i, j;

            for (i = 0; i < M; i++)
            {
                for (j = 0; j < N; j++)
                {
                    result.data[i * prd_a + j] = data[i * prd_b + j]-mpose->data[i * prd_c + j];
                }
            }
        }
        return result;
    }
    virtual minimatrix* Retract(const minimatrix* mpose)
    {
        minimatrix* result=new minimatrix(size1,size2);
        const size_t M = size1;
        const size_t N = size2;

        if(mpose->size1!=M||mpose->size2!=N)
        {
            throw std::invalid_argument("matrices must have same dimensions");
        }
        else
        {
            const size_t prd_a = result->prd;
            const size_t prd_b = this->prd;
            const size_t prd_c = mpose->prd;

            size_t i, j;

            for (i = 0; i < M; i++)
            {
                for (j = 0; j < N; j++)
                {
                    result->data[i * prd_a + j] = data[i * prd_b + j]+mpose->data[i * prd_c + j];
                }
            }
        }
        return result;
    }
    virtual minimatrix between(const minimatrix* mpose) const
    {
        size_t i, j;
        minimatrix result(size1,size2);
        for (i = 0; i < size1; i++)
        {
            for (j = 0; j < size2; j++)
            {
                result.data[i * prd + j] = mpose->data[i * prd + j]-data[i * prd + j];
            }
        }

        return result;
    }
    virtual minimatrix between(const minimatrix* mpose,minimatrix& H1,minimatrix& H2) const
    {

        const size_t p = dimension ;
        const size_t q = dimension ;
        size_t i, j;

        const double zero = 0.0;
        const double n_one = -1.0;
        const double one = 1.0;

        if(H1.size1!=dimension||H1.size2!=dimension)
        {

            if(H1.data!=0)
            {
                free(H1.data);
            }

            H1.data=(double *) malloc (dimension * dimension * sizeof (double));
            if (H1.data == 0)
            {
                throw std::invalid_argument("failed to allocate space for block");
            }
            H1.size1 = dimension;
            H1.size2 = dimension;
            H1.prd = dimension;
            H1.owner = 1;
        }
        for (i = 0; i < p; i++)
        {
            for (j = 0; j < q; j++)
            {
                *(double *) (H1.data +  (i * prd + j)) = ((i == j) ? n_one : zero);
            }
        }
        if(H2.size1!=dimension||H2.size2!=dimension)
        {

            if(H2.data!=0)
            {
                free(H2.data);
            }

            H2.data=(double *) malloc (dimension * dimension * sizeof (double));
            if (H2.data == 0)
            {
                throw std::invalid_argument("failed to allocate space for block");
            }
            H2.size1 = dimension;
            H2.size2 = dimension;
            H2.prd = dimension;
            H2.owner = 1;
        }
        for (i = 0; i < p; i++)
        {
            for (j = 0; j < q; j++)
            {
                *(double *) (H2.data +  (i * prd + j)) = ((i == j) ? one : zero);
            }
        }
        minimatrix result(size1,size2);
        for (i = 0; i < size1; i++)
        {
            for (j = 0; j < size2; j++)
            {
                result.data[i * prd + j] = mpose->data[i * prd + j]-data[i * prd + j];
            }
        }

        return result;
    }


    ~minimatrix()
    {
        if((owner)&&(data!=NULL))
        {
            free(data);
            owner=0;
        }
    }
} ;


minimatrix minimatrix_mat3(const double x00, const double x01,const double x02,
                           const double x10, const double x11,const double x12,
                           const double x20, const double x21,const double x22 );
void minimatrix_resize(minimatrix* m, const size_t n1, const size_t n2);
void minimatrix_print(const minimatrix * m);
void minimatrix_print(const minimatrix& m);

minimatrix
minimatrix_blockmatrix_var (minimatrix * m,
                        const size_t i, const size_t j,
                        const size_t n1, const size_t n2);

minimatrix
minimatrix_varprd_blockmatrix (minimatrix * m,
                               const size_t i, const size_t j,
                               const size_t n1, const size_t n2);

minimatrix
minimatrix_varprd_blockmatrix (const minimatrix& m,
                               const size_t i, const size_t j,
                               const size_t n1, const size_t n2);
minimatrix
minimatrix_blockmatrix (const minimatrix& m,
                        const size_t i, const size_t j,
                        const size_t n1, const size_t n2);
void minimatrix_set_zero (minimatrix * m);
void minimatrix_set_identity (minimatrix * m);
minimatrix minimatrix_identity_mat(int n);
void minimatrix_set_neg_identity (minimatrix * m);
void minimatrix_set_all (minimatrix * m, double x);
int minimatrix_memcpy(minimatrix * dest, const minimatrix* src);
int minimatrix_memcpy(minimatrix * dest, const minimatrix& src);
int minimatrix_copy_diffmatrix(minimatrix * dest, const minimatrix& src);
int minimatrix_triangularViewUpper(minimatrix* dest,const minimatrix& src);
int minimatrix_triangularViewLower(minimatrix* dest,const minimatrix& src);
int minimatrix_transpose (minimatrix * m);
minimatrix minimatrix_transmat (const minimatrix& m);
int minimatrix_transpose_memcpy (minimatrix * dest, const minimatrix * src);
int minimatrix_transpose_memcpy (minimatrix * dest, const minimatrix& src);
int minimatrix_equal (const minimatrix& a, const minimatrix&  b);
int minimatrix_add (minimatrix * a, const minimatrix& b);
int minimatrix_mat_add (minimatrix * a, const minimatrix& b, const minimatrix& c);
int minimatrix_mat_add (minimatrix * a, const minimatrix& b,double alpha, const minimatrix& c);
int minimatrix_mat_add (minimatrix * a,double alpha, const minimatrix& b,double beta, const minimatrix& c);
int minimatrix_add_by(minimatrix * a, double alpha, const minimatrix& c);
int minimatrix_add_transpose (minimatrix * a, const minimatrix& b);
int minimatrix_sub (minimatrix * a, const minimatrix * b);
int minimatrix_mat_sub (minimatrix * a, const minimatrix& b,const minimatrix& c);
int minimatrix_mat_sub (minimatrix * a, const minimatrix* b,const minimatrix* c);
int minimatrix_scale (minimatrix * a, const double x);
int minimatrix_mat_scale (minimatrix * a, const double x,const minimatrix& b);
void minimatrix_triangularView_setLowZero(minimatrix * a,unsigned int n);
minimatrix minimatrix_toprows(const minimatrix& a,int n);
minimatrix minimatrix_selfadjointview(const minimatrix& a,bool upperorlower=true);
double   minimatrix_get(const minimatrix * m, const size_t i, const size_t j);
double   minimatrix_get(const minimatrix& m, const size_t i, const size_t j);
void    minimatrix_set(minimatrix * m, const size_t i, const size_t j, const double x);

#endif
