/**
 *  @file   PhaseFactor.cpp
 *  @author Ryan Watson & Jason Gross
 *  @brief  Implementation file for carrier-phase factor
 **/

#include "PhaseFactor.h"

using namespace std;

namespace minisam {

Eigen::VectorXd PhaseFactor::evaluateError(const Eigen::VectorXd& q, const Eigen::VectorXd& g) const
{
    Eigen::VectorXd h=obsMap(satXYZ_, nomXYZ_, 1);
     double est = (h.transpose() * q) + g[0];
     Eigen::VectorXd result(1);
     result<<est-measured_;
     return result;

}

Eigen::VectorXd PhaseFactor::evaluateError(const Eigen::VectorXd& q, const Eigen::VectorXd& g,
                              Eigen::MatrixXd& H1,Eigen::MatrixXd& H2) const
{
    Eigen::VectorXd h=obsMap(satXYZ_, nomXYZ_, 1);
     double est = (h.transpose() * q) + g[0];
     H1.resize(1,5);
     H1<<h(0),h(1),h(2),h(3),h(4);
     H2.resize(1,1);
     H2<<1.0;
     Eigen::VectorXd result(1);
     result<<est-measured_;
     return result;


}
};
