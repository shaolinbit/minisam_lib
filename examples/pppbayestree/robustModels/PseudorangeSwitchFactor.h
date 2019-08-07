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

namespace minisam {

class  PseudorangeSwitchFactor: public NoiseModelFactor2
{
private:
  typedef NoiseModelFactor2  Base;
  double measured_;
  gnssStateVec h_;

public:

  typedef PseudorangeSwitchFactor This;

  PseudorangeSwitchFactor(): measured_(0) {
  h_.resize(5);
  h_<<1,1,1,1,1;

   }

  virtual ~PseudorangeSwitchFactor() {}

  double measured()
{
  return measured_;
}

Eigen::VectorXd h()
{
  return h_;
}


  PseudorangeSwitchFactor(int j, int k, const double deltaObs, const Eigen::VectorXd& obsMap, GaussianNoiseModel* model):
    Base(model, j,k), measured_(deltaObs) {h_=gnssStateVec(obsMap);}

  virtual NonlinearFactor* clone()  {
    PseudorangeSwitchFactor* npsf=new PseudorangeSwitchFactor(key1(),key2(),measured(),h(),noiseModel());
    return npsf;
}


  /// vector of errors

        Eigen::VectorXd evaluateError(const Eigen::VectorXd & q,
      const Eigen::VectorXd & s) const;//SwitchVariableLinear is a double type.
        Eigen::VectorXd evaluateError(const Eigen::VectorXd & q,
      const Eigen::VectorXd & s,Eigen::MatrixXd& H1,Eigen::MatrixXd& H2) const;//SwitchVariableLinear is a double type.



}; // PseudorangeSwitchFactor Factor
}; // namespace
