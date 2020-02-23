#ifndef MEMORYPOOL_H
#define MEMORYPOOL_H

#include <stdlib.h>
#include <stdio.h>
#include <stdexcept>


enum
{
    MINI_SUCCESS  = 0,
    MINI_FAILURE  = -1
} ;

struct memory_pool
{
    size_t size;
    int index;
    double *data;

    memory_pool():size(0),index(0),data(NULL)
    {

    }
    memory_pool(size_t n):size(n),index(0)
    {
        data=(double *) malloc (n * sizeof (double));
        if(data==NULL)
        {
             throw std::invalid_argument("Memory pool allocated failed.");
        }
    }
    ~memory_pool()
    {
        if(data!=NULL)
        {
            free(data);
        }
    }
};
struct mini_int_vector
{
    size_t size;
    size_t *data;
    mini_int_vector():size(0),data(NULL)
    {

    }
    mini_int_vector(size_t m):size(m)
    {
        data=(size_t *) malloc (m* sizeof (size_t));
    }
    ~mini_int_vector()
    {
        if(data!=NULL)
        {
            free(data);
        }
    }

};



#endif
