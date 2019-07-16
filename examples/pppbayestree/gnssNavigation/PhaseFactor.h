/**
 *  @file   PhaseFactor.h
 *  @author Ryan Watson & Jason Gross
 *  @brief  Header file for Carrier-Phase Factor
 **/

#pragma once

#include "../gnssNavigation/GnssTools.h"
#include "minisam/nonlinear/NonlinearFactor.h"
#include "../gnssNavigation/nonBiasStates.h"

namespace minisam {

class  PhaseFactor : public NoiseModelFactor2
{
private:
typedef NoiseModelFactor2 Base;
Eigen::Vector3d satXYZ_;
Eigen::Vector3d nomXYZ_;
double measured_;
nonBiasStates h_;

public:

typedef PhaseFactor This;

PhaseFactor() : measured_(0) {

        Eigen::VectorXd bh(5) ;
        bh<<1,1,1,1,1;

        h_=nonBiasStates(bh);
}

virtual ~PhaseFactor() {
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

PhaseFactor(int deltaStates, int bias, const double measurement,
            const Eigen::Vector3d& satXYZ, const Eigen::Vector3d& nomXYZ,SharedNoiseModel* model) :
        Base(model, deltaStates, bias), measured_(measurement)
{
        satXYZ_=satXYZ;
        nomXYZ_=nomXYZ;
}


virtual NonlinearFactor* clone() const{
       PhaseFactor* npf=new PhaseFactor(key1(),key2(),measured_,satXYZ_,nomXYZ_,noiseModel());
       return npf;
}

Eigen::VectorXd evaluateError(const Eigen::VectorXd& q, const Eigen::VectorXd& g) const;

Eigen::VectorXd evaluateError(const Eigen::VectorXd& q, const Eigen::VectorXd& g,
                              Eigen::MatrixXd& H1,Eigen::MatrixXd& H2) const;


}; // PhaseFactor Factor
}; // namespace
