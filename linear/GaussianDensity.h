#ifndef GAUSSIANDENSITY_H
#define GAUSSIANDENSITY_H

/**
 * @file    GaussianDensity.h
 * @brief   A Gaussian Density
 * @author  Frank Dellaert
 * @date    Jan 21, 2012
 */

// \callgraph
//#pragma once

#include "../linear/GaussianConditional.h"


/**
* A Gaussian density.
*
* It is implemented as a GaussianConditional without parents.
* The negative log-probability is given by \f$ |Rx - d|^2 \f$
* with \f$ \Lambda = \Sigma^{-1} = R^T R \f$ and \f$ \mu = R^{-1} d \f$
*/
class  GaussianDensity : public GaussianConditional
{

public:

    //typedef boost::shared_ptr<GaussianDensity> shared_ptr;

    /// default constructor needed for serialization
    GaussianDensity() :
        GaussianConditional()
    {
    }

    /// Copy constructor from GaussianConditional
    GaussianDensity(const GaussianConditional& conditional) :
        GaussianConditional(conditional)
    {

        //if((conditional.nrParents()).size()!= 0)
        //  throw std::invalid_argument("GaussianDensity can only be created from a conditional with no parents");
    }

    /// constructor using d, R
    GaussianDensity(int key, const Eigen::VectorXd& d, const Eigen::MatrixXd& R, DiagonalNoiseModel* noiseModel =new DiagonalNoiseModel()) :
        GaussianConditional(key, d, R, noiseModel) {}

    virtual ~GaussianDensity();

    /// Construct using a mean and variance
    static GaussianDensity FromMeanAndStddev(int key, const Eigen::VectorXd& mean, const double& sigma);

    /// print
    // void print(const std::string& = "GaussianDensity",
    //  const KeyFormatter& formatter = DefaultKeyFormatter) const;

    /// Mean \f$ \mu = R^{-1} d \f$
    Eigen::VectorXd mean() const;

    /// Covariance matrix \f$ \Sigma = (R^T R)^{-1} \f$
    Eigen::MatrixXd covariance() const;

    int GDgetDim(std::vector<int>::const_iterator variable) const;

};

#endif // GAUSSIANDENSITY_H
