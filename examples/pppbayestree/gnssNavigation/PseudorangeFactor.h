/**
 *  @file   PseudorangeFactor.h
 *  @author Ryan Watson & Jason Gross
 *  @brief  Header file for Pseudorange factor
 **/

#pragma once

#include "../gnssNavigation/GnssTools.h"
#include "minisam/nonlinear/NonlinearFactor.h"
#include "../gnssNavigation/nonBiasStates.h"

namespace minisam {

class  PseudorangeFactor : public NoiseModelFactor1
{
private:
typedef NoiseModelFactor1 Base;
Eigen::Vector3d nomXYZ_;
Eigen::Vector3d satXYZ_;
nonBiasStates h_;
double measured_;

public:

typedef PseudorangeFactor This;

PseudorangeFactor() : measured_(0) {
         Eigen::VectorXd bh(5);
        bh<<1,1,1,1,1;

        h_=nonBiasStates(bh);
}

virtual ~PseudorangeFactor() {
}

Eigen::Vector3d satXYZ()
{
  return satXYZ_;
}

Eigen::Vector3d nomXYZ()
{
  return nomXYZ_;
}

double measured()
{
  return measured_;
}


PseudorangeFactor(int key, const double deltaObs, const Eigen::Vector3d& satXYZ,
const Eigen::Vector3d& nomXYZ, SharedNoiseModel* model) :
        Base(model, key), measured_(deltaObs), satXYZ_(satXYZ) {
        nomXYZ_=nomXYZ;
}


virtual NonlinearFactor* clone() const
{
   PseudorangeFactor* nprf=new PseudorangeFactor(key(),measured_,satXYZ_,nomXYZ_,noiseModel());
   return nprf;

}
/// vector of errors

 virtual Eigen::VectorXd unwhitenedError(const std::map<int, Eigen::VectorXd>& x) const
    {
        std::map<int, Eigen::VectorXd>::const_iterator xbegin = x.find(keys_[0]);
        Eigen::VectorXd x1 = xbegin->second;
        return evaluateError(x1);
    }

Eigen::VectorXd  evaluateError(const Eigen::VectorXd& q) const;


Eigen::VectorXd  evaluateError(const Eigen::VectorXd& q,
                     Eigen::MatrixXd& H) const;


}; // PseudorangeFactor Factor
}; // namespace
