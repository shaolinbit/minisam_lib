/**
 *  @file   PseudorangeMaxMix.h
 *  @author Ryan
 *  @brief  Header file for Pseudorange Max-Mix factor
 **/

#pragma once

#include <Eigen/Eigen>
#include "minisam/linear/NoiseModel.h"
#include "minisam/nonlinear/NonlinearFactor.h"
#include "../gnssNavigation/FolderUtils.h"
#include "../gnssNavigation/gnssStateVec.h"


using namespace std;

namespace minisam {

class  PseudorangeMaxMix: public NoiseModelFactor1//<gnssStateVec> {
{
private:
  typedef NoiseModelFactor1  Base;
  double measured_, w_, hyp_;
  gnssStateVec h_;
  SharedNoiseModel* nullHypothesisModel_;


public:


  typedef PseudorangeMaxMix This;

  PseudorangeMaxMix(): measured_(0) {
  h_.resize(5);
  h_<<1,1,1,1,1;
   }

  virtual ~PseudorangeMaxMix() {}

   double measured()
{
  return measured_;
}

Eigen::VectorXd h()
{
  return h_;
}

double w()
{return w_;}

double hyp()
{return hyp_;}

SharedNoiseModel* noiseModel2()
{
  return nullHypothesisModel_;
}


   PseudorangeMaxMix(int key, const double deltaObs, const Eigen::VectorXd& obsMap,
    SharedNoiseModel* model1, SharedNoiseModel* model2,
    const double& hypNoise, double weight): NoiseModelFactor1(model1, key),
    measured_(deltaObs), h_(obsMap), w_(weight), hyp_(hypNoise),
    nullHypothesisModel_(model2) {  };

 virtual NonlinearFactor* clone() const {
    PseudorangeMaxMix* npmm=new PseudorangeMaxMix(key(),measured_,h_,noiseModel_,nullHypothesisModel_,hyp_,w_);
    return npmm;
}


/// vector of errors
  Eigen::VectorXd evaluateError(const Eigen::VectorXd& q) const;
  Eigen::VectorXd evaluateError(const Eigen::VectorXd& q,Eigen::MatrixXd& H1 ) const;


}; // PseudorangeMaxMix Factor
//} // namespace
};
