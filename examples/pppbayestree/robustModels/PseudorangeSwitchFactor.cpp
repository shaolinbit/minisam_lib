/**
 *  @file   PseudorangeSwitchFactor.cpp
 *  @author Ryan Waton and Jason Gross
 *  @brief  Implementation file for pseudorange switchable factor
 **/

#include "PseudorangeSwitchFactor.h"

using namespace std;

namespace minisam
{


minivector PseudorangeSwitchFactor::evaluateError(const minimatrix* q,
        const minimatrix* s) const//SwitchVariableLinear is a double type.
{
    //double error = (h_.transpose()*q)-measured_;
    double error;
    miniblas_ddot(*h_,minivector(q),&error);
    error-=measured_;
    error *= minimatrix_get(s,0,0);//error *= s.data[0];
    minivector result(1);
    result.data[0]=error;
    // result<<error;
    return result;
}
minivector PseudorangeSwitchFactor::evaluateError(const minimatrix* q,
        const minimatrix* s,minimatrix& H1,minimatrix& H2) const//SwitchVariableLinear is a double type.
{

    // double error = (h_.transpose()*q)-measured_;
    // error *= s(0);
    double error;
    miniblas_ddot(*h_,minivector(q),&error);
    error-=measured_;
    error *= minimatrix_get(s,0,0);//error *= s.data[0];
    minivector result(1);
    result.data[0]=error;

    //H1.resize(1,5);
    //H1<<h_.transpose()*s(0);
    minimatrix_resize(&H1,1,5);
    /*
    H1.data[0]=h_.gnssstate_vector_.data[0]*s.data[0];
    H1.data[1]=h_.gnssstate_vector_.data[1]*s.data[0];
    H1.data[2]=h_.gnssstate_vector_.data[2]*s.data[0];
    H1.data[3]=h_.gnssstate_vector_.data[3]*s.data[0];
    H1.data[4]=h_.gnssstate_vector_.data[4]*s.data[0];
    */
    minimatrix_set(&H1,0,0,minivector_get(h_,0)*minimatrix_get(s,0,0));
    minimatrix_set(&H1,0,1,minivector_get(h_,1)*minimatrix_get(s,0,0));
    minimatrix_set(&H1,0,2,minivector_get(h_,2)*minimatrix_get(s,0,0));
    minimatrix_set(&H1,0,3,minivector_get(h_,3)*minimatrix_get(s,0,0));
    minimatrix_set(&H1,0,4,minivector_get(h_,4)*minimatrix_get(s,0,0));


    //H2.resize(1,1);
    // H2<<error;
    minimatrix_resize(&H2,1,1);
    H2.data[0]=error;

    return result;
}
};
