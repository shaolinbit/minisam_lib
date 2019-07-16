/**
 *  @file   GnssBetweenFactor.cpp
 *  @author Ryan Watson & Jason Gross
 *  @brief  Implementation file for GnssBetweenFactor
 **/

#include "GnssBetweenFactor.h"

using namespace std;

namespace minisam
{
 Eigen::VectorXd
    GnssBetweenFactor::evaluateError(const Eigen::VectorXd& q, const Eigen::VectorXd& p) const
    {

        double est = GNSS_norm5(nonBiasStates(q-p));
        return (Eigen::VectorXd(1) << est ).finished();


    }
    Eigen::VectorXd
    GnssBetweenFactor::evaluateError(const Eigen::VectorXd& q, const Eigen::VectorXd& p, Eigen::MatrixXd &H1, Eigen::MatrixXd &H2) const
    {

      Eigen::VectorXd h(5);
      h <<1,1,1,1,1;
      H1.resize(1,5);
      H1<<h;
      H2.resize(1,5);
      H2<<h;
        double est = GNSS_norm5(nonBiasStates(q-p));
        return (Eigen::VectorXd(1) << est ).finished();

    }

};
