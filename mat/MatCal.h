#ifndef MATCAL_H_INCLUDED
#define MATCAL_H_INCLUDED

#include "../gmfconfig.h"
#include <vector>
#include <map>
#include <list>
#include  <algorithm>
#include "Matrix.h"
#include "../miniblas/minimatrix_double.h"
#include "../miniblas/minivector_double.h"
/*
This file was added by Shaolin for matrix function in isam.
*/



namespace minisam
{



double dot3d(const minivector* p, const minivector* q,
             minimatrix* H1,
             minimatrix* H2);
double dot3d(const minivector& p, const minivector& q,
             minimatrix* H1,
             minimatrix* H2);


double dot3d(const minivector* p, const minivector* q);

double norm3d(const minivector* p);

double norm3d(const minivector& p);


double norm3d(const minivector& p, minimatrix* H);


minivector* cross3d(const minivector* p,
                    const minivector* q,
                    minimatrix* H1,
                    minimatrix* H2);
minivector cross3d(const minivector& p,
                   const minivector& q,
                   minimatrix* H1,
                   minimatrix* H2);

minivector* cross3d(const minivector* p,
                    const minivector* q,
                    minimatrix* H1);
minivector cross3d(const minivector& p,
                   const minivector& q,
                   minimatrix* H1);

minivector cross3d(const minivector& p,const minivector& q);


minivector* cross3d(const minivector* p,const minivector* q);


minivector* normalize3d(const minivector* p,
                        minimatrix* H);
minivector normalize3d(const minivector& p,
                       minimatrix* H);



void normalize3d(const minivector& p,
                 minivector* H);
void normalize3d(minivector* p);


minivector* normalize3d(const minivector* p);
minivector normalize3d(const minivector& p);

minivector VectorValuesgetvector(const std::map<int,minivector>& VectorValues,const std::vector<int>& keys);


std::map<int,minivector> VectorValuesClone(const std::map<int,minivector>& vectorvalues);
void  VectorValuesClone(const std::map<int,minivector>& vectorvalues,std::map<int,minivector>* cvv);

double VectorValuesDot(const std::map<int,minivector>& vectorvalues);
double VectorValuesDot(const std::map<int,minivector>& vectorvalues1,const std::map<int,minivector>& vectorvalues2);


void VectorValuesScaleMultiply(std::map<int,minivector> *vectorvalues,double a);


void VectorValuesaxpy(std::map<int,minivector>* y,
                      const std::map<int,minivector>& x,double a);




std::vector<int> getvectorfrom3keys(int key1,int key2,int key3);
std::vector<int> getvectorfrom5keys(int key1,int key2,int key3,int key4,int key5);
std::vector<int> getvectorfrom6keys(int key1,int key2,int key3,int key4,int key5,int key6);

std::map<int,minivector> VectorValuesZero(const std::map<int,minivector>& vvz);

void VectorValuesSetZero(std::map<int,minivector>* vvz);
std::map<int,minivector> VectorValuesRetract(const std::map<int,minivector>& vectorvalues1,
        const std::map<int,minivector>& vectorvalues2);
std::map<int,minivector> DVectorValuesRetract(const std::map<int,minivector>& vectorvalues1,
        const std::map<int,minivector>& vectorvalues2);
std::map<int,minimatrix*> ValuesRetract(const std::map<int,minimatrix*>& vectorvalues1,
                                        const std::map<int,minivector>& vectorvalues2);
void ValuesRetract_atmap(const std::map<int,minimatrix*>& vectorvalues1,
                         const std::map<int,minivector>& vectorvalues2,std::map<int,minimatrix*>& vectorvalues3);

void DVectorValuesRetract_atmap(const std::map<int,minivector>& vectorvalues1,
                                const std::map<int,minivector>& vectorvalues2,std::map<int,minivector>& vectorvalues3);

std::map<int,minivector> VectorValuesSub(const std::map<int,minivector>& vectorvalues1,
        const std::map<int,minivector>& vectorvalues2);

std::map<int,minivector> VectorValuesAdd(const std::map<int,minivector>& vectorvalues1,
        const std::map<int,minivector>& vectorvalues2);



void VectorValuesUpdate(std::map<int,minivector>* vectorvalues1,const std::map<int,minivector>& vectorvalues2);

double VectorValuesNorm(const std::map<int,minivector>& vvz);


double VectorValuesSquareNorm(const std::map<int,minivector>& vvz);


};
#endif // MATRIX_H_INCLUDED
