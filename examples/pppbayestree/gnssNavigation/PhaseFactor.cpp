/**
 *  @file   PhaseFactor.cpp
 *  @author Ryan Watson & Jason Gross
 *  @brief  Implementation file for carrier-phase factor
 **/

#include "PhaseFactor.h"

using namespace std;

namespace minisam
{

minivector PhaseFactor::evaluateError(const minimatrix* q, const minimatrix* g) const
{
    minivector h=obsMap(satXYZ_, nomXYZ_, 1);
    //double est = (h.transpose() * q) + g[0];
    double est;
    miniblas_ddot(h,minivector(q),&est);
    est+=g->data[0];

    minivector result(1,est-measured_);
    //result.data[0]=est-measured_;
    return result;

}

minivector PhaseFactor::evaluateError(const minimatrix* q, const minimatrix* g,
                                      minimatrix& H1,minimatrix& H2) const
{
    minivector h=obsMap(satXYZ_, nomXYZ_, 1);
    //double est = (h.transpose() * q) + g[0];
    double est;
    miniblas_ddot(h,minivector(q),&est);
    est+=g->data[0];

    // H1.resize(1,5);
    // H1<<h(0),h(1),h(2),h(3),h(4);
    minimatrix_resize(&H1,1,5);
    H1.data[0]=h.data[0];
    H1.data[1]=h.data[1];
    H1.data[2]=h.data[2];
    H1.data[3]=h.data[3];
    H1.data[4]=h.data[4];

    //H2.resize(1,1);
    //H2<<1.0;
    minimatrix_resize(&H2,1,1);
    H2.data[0]=1.0;

    minivector result(1,est-measured_);
    //result.data[0]=est-measured_;


    return result;


}
};
