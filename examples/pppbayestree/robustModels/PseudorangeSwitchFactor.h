/**
 *  @file   PseudorangeSwitchFactor.h
 *  @author Ryan Watson & Jason Gross
 *  @brief  Header file for Pseudorange Switchable factor
 **/

#pragma once


#include "minisam/nonlinear/NonlinearFactor.h"
#include "../gnssNavigation/gnssStateVec.h"
//#include "../robustModels/switchVariableLinear.h"


using namespace std;

namespace minisam
{

class  PseudorangeSwitchFactor: public NoiseModelFactor2
{
private:
    typedef NoiseModelFactor2  Base;
    double measured_;
    minivector* h_;

public:

    typedef PseudorangeSwitchFactor This;

    PseudorangeSwitchFactor(): measured_(0),h_(new gnssStateVec(1.0,1.0,1.0,1.0,1.0))
    {


    }

    virtual ~PseudorangeSwitchFactor() {
    delete h_;
    }

    double measured()
    {
        return measured_;
    }

    minivector* h()
    {
        return h_;
    }


    PseudorangeSwitchFactor(int j, int k, const double deltaObs, minivector* obsMap, GaussianNoiseModel* model):
        Base(model, j,k), measured_(deltaObs),h_(obsMap)
    {
    
    }

    virtual NoiseModelFactor* clone()
    {
        PseudorangeSwitchFactor* npsf=new PseudorangeSwitchFactor(key1(),key2(),measured(),h(),noiseModel());
        return npsf;
    }


    /// vector of errors

    virtual minivector evaluateError(const minimatrix* q,
                             const minimatrix* s) const;//SwitchVariableLinear is a double type.
    virtual minivector evaluateError(const minimatrix* q,
                             const minimatrix* s,minimatrix& H1,minimatrix& H2) const;//SwitchVariableLinear is a double type.



}; // PseudorangeSwitchFactor Factor
}; // namespace
