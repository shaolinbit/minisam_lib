/**
 *  @file   GnssBetweenFactor.cpp
 *  @author Ryan Watson & Jason Gross
 *  @brief  Implementation file for GnssBetweenFactor
 **/

#include "GnssBetweenFactor.h"

using namespace std;

namespace minisam
{
minivector
GnssBetweenFactor::evaluateError(const minimatrix* q, const minimatrix* p) const
{

    //double est = GNSS_norm5(nonBiasStates(q-p));
    double est = GNSS_norm5(nonBiasStates(q->data[0]-p->data[0],q->data[1]-p->data[1],
    q->data[2]-p->data[2],q->data[3]-p->data[3],q->data[4]-p->data[4]));
    minivector result(1,est);
    //result.data[0]=est;
    return result;
    //return (minivector(1) << est ).finished();

}
minivector
GnssBetweenFactor::evaluateError(const minimatrix* q, const minimatrix* p, minimatrix &H1, minimatrix &H2) const
{

    // h.data[0]=1.0;h.data[1]=1.0;h.data[2]=1.0;h.data[3]=1.0;h.data[4]=1.0;
    // h <<1,1,1,1,1;
    // H1.resize(1,5);
    // H1<<h;
    minimatrix_resize(&H1,1,5);
    minimatrix_set_all(&H1,1.0);
    /*
    H1.data[0]=1.0;
    H1.data[1]=1.0;
    H1.data[2]=1.0;
    H1.data[3]=1.0;
    H1.data[4]=1.0;
    */
    // H2.resize(1,5);
    // H2<<h;
    //double est = GNSS_norm5(nonBiasStates(q-p));
    //return (minivector(1) << est ).finished();
    double est = GNSS_norm5(nonBiasStates(q->data[0]-p->data[0],
    q->data[1]-p->data[1],q->data[2]-p->data[2],
        q->data[3]-p->data[3],q->data[4]-p->data[4]));
    minivector result(1,est);
   // result.data[0]=est;
    return result;

}

};
