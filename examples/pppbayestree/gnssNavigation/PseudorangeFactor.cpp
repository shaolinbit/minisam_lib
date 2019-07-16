/**
 *  @file   PseudorangeFactor.cpp
 *  @author Ryan Watson & Jason Gross
 *  @brief  Implementation file for pseudorange factor
 **/

#include "PseudorangeFactor.h"

using namespace std;

namespace minisam {


Eigen::VectorXd  PseudorangeFactor::evaluateError(const Eigen::VectorXd& q) const
{

        Eigen::VectorXd  h = obsMap(satXYZ_, nomXYZ_, 1);
        Eigen::VectorXd result(1);
        double est = h.dot(q);
        result<<est-measured_;
        return result;
}


Eigen::VectorXd  PseudorangeFactor::evaluateError(const Eigen::VectorXd& q,
                     Eigen::MatrixXd& H) const
{
  Eigen::VectorXd  h = obsMap(satXYZ_, nomXYZ_, 1);
   H.resize(1,5);
   H<<h(0),h(1),h(2),h(3),h(4);
    double est = h.dot(q);
        Eigen::VectorXd result(1);
        result<<est-measured_;
        return result;

}
};
