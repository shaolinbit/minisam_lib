/**
 *  @file   PseudorangeFactor.cpp
 *  @author Ryan Watson & Jason Gross
 *  @brief  Implementation file for pseudorange factor
 **/

#include "PseudorangeFactor.h"

using namespace std;

namespace minisam
{


minivector  PseudorangeFactor::evaluateError(const minimatrix* q) const
{

    minivector  h = obsMap(*satXYZ_, *nomXYZ_, 1);
    // double est = h.dot(q);
    double est;
    miniblas_ddot(h,q,&est);
    minivector result(1,est-measured_);
    //result<<est-measured_;
   // result.data[0]=est-measured_;

    return result;
}


minivector  PseudorangeFactor::evaluateError(const minimatrix* q,
        minimatrix& H) const
{
    minivector  h = obsMap(*satXYZ_, *nomXYZ_, 1);
    // H.resize(1,5);
    // H<<h(0),h(1),h(2),h(3),h(4);
    minimatrix_resize(&H,1,5);
    /*
    H.data[0]=h.data[0];
    H.data[1]=h.data[1];
    H.data[2]=h.data[2];
    H.data[3]=h.data[3];
    H.data[4]=h.data[4];*/
    minimatrix_set(&H,0,0,h.data[0]);
    minimatrix_set(&H,0,1,h.data[1]);
    minimatrix_set(&H,0,2,h.data[2]);
    minimatrix_set(&H,0,3,h.data[3]);
    minimatrix_set(&H,0,4,h.data[4]);


    double est;
    miniblas_ddot(h,minivector(q),&est);
    minivector result(1,est-measured_);
    //result<<est-measured_;
    //result.data[0]=est-measured_;

    if(DEBUGSTATE)
    {
    minimatrix_print(H);
    cout<<endl;
    cout<<"result"<<endl;
    minimatrix_print(result);
    }



    return result;

}
};
