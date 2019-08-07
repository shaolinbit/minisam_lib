/**
 *  @file   GnssBetweenFactor.h
 *  @author Ryan Watson & Jason Gross
 *  @brief  Header file for GnssBetweenFactor
 **/

#pragma once


#include "../gnssNavigation/nonBiasStates.h"
#include "minisam/nonlinear/NonlinearFactor.h"

namespace minisam {

class  GnssBetweenFactor: public NoiseModelFactor2
 {

public:


  GnssBetweenFactor(int state1, int state2, GaussianNoiseModel* model):
    NoiseModelFactor2(model, state1, state2) { }

 virtual NonlinearFactor* clone() const
    {
        GnssBetweenFactor *newfactor = new GnssBetweenFactor(key1(),key2(),noiseModel());
        return newfactor;
    }


    virtual Eigen::VectorXd
    evaluateError(const Eigen::VectorXd& X1, const Eigen::VectorXd& X2) const;

    virtual Eigen::VectorXd
    evaluateError(const Eigen::VectorXd& X1, const Eigen::VectorXd& X2, Eigen::MatrixXd &H1, Eigen::MatrixXd &H2) const;


}; // PseudorangeFactor Factor
}; // namespace
