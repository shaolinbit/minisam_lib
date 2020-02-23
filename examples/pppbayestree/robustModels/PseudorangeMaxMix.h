/**
 *  @file   PseudorangeMaxMix.h
 *  @author Ryan
 *  @brief  Header file for Pseudorange Max-Mix factor
 **/

#pragma once

#include "minisam/linear/NoiseModel.h"
#include "minisam/nonlinear/NonlinearFactor.h"
#include "../gnssNavigation/FolderUtils.h"
#include "../gnssNavigation/gnssStateVec.h"


using namespace std;

namespace minisam
{

class  PseudorangeMaxMix: public NoiseModelFactor1
{
private:
    typedef NoiseModelFactor1  Base;
    double measured_, w_, hyp_;
    minivector* h_;//gnssStateVec h_;
    GaussianNoiseModel* nullHypothesisModel_;


public:


    typedef PseudorangeMaxMix This;

    PseudorangeMaxMix(): measured_(0),h_(new gnssStateVec(1.0,1.0,1.0,1.0,1.0))
    {

    }

    virtual ~PseudorangeMaxMix() {}

    double measured()
    {
        return measured_;
    }

    minivector* h()
    {
        return h_;
    }

    double w()
    {
        return w_;
    }

    double hyp()
    {
        return hyp_;
    }

    GaussianNoiseModel* noiseModel2()
    {
        return nullHypothesisModel_;
    }


    PseudorangeMaxMix(int key, const double deltaObs, minivector* obsMap,
                      GaussianNoiseModel* model1, GaussianNoiseModel* model2,
                      const double& hypNoise, double weight): NoiseModelFactor1(model1, key),
        measured_(deltaObs),  w_(weight), hyp_(hypNoise),h_(obsMap),
        nullHypothesisModel_(model2)
    {
       

    };

    virtual NoiseModelFactor* clone() const
    {
        PseudorangeMaxMix* npmm=new PseudorangeMaxMix(key(),measured_,h_,noiseModel_,nullHypothesisModel_,hyp_,w_);
        return npmm;
    }


/// vector of errors
    virtual minivector evaluateError(const minimatrix* q) const;
    virtual minivector evaluateError(const minimatrix* q,minimatrix& H1 ) const;


}; // PseudorangeMaxMix Factor
//} // namespace
};
