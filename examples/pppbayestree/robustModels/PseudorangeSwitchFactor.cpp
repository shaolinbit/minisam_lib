/**
 *  @file   PseudorangeSwitchFactor.cpp
 *  @author Ryan Waton and Jason Gross
 *  @brief  Implementation file for pseudorange switchable factor
 **/

#include "PseudorangeSwitchFactor.h"

using namespace std;
//using namespace vertigo;

namespace minisam {


  Eigen::VectorXd PseudorangeSwitchFactor::evaluateError(const Eigen::VectorXd & q,
      const Eigen::VectorXd & s) const//SwitchVariableLinear is a double type.
      {
     double error = (h_.transpose()*q)-measured_;
    error *= s(0);
    Eigen::VectorXd result(1);
    result<<error;
    return result;
      }
        Eigen::VectorXd PseudorangeSwitchFactor::evaluateError(const Eigen::VectorXd & q,
      const Eigen::VectorXd & s,Eigen::MatrixXd& H1,Eigen::MatrixXd& H2) const//SwitchVariableLinear is a double type.
      {

   double error = (h_.transpose()*q)-measured_;
    error *= s(0);
    Eigen::VectorXd result(1);
    H1.resize(1,5);
    H1<<h_.transpose()*s(0);
    H2.resize(1,1);
    H2<<error;
    result<<error;
    return result;
      }
};
