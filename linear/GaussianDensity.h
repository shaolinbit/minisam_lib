#ifndef GAUSSIANDENSITY_H
#define GAUSSIANDENSITY_H

/**
 * @file    GaussianDensity.h
 * @brief   A Gaussian Density
 */


#include "../linear/GaussianConditional.h"

namespace minisam {

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

    /// default constructor needed for serialization
    GaussianDensity() :
        GaussianConditional()
    {
    }

    /// Copy constructor from GaussianConditional
    GaussianDensity(const GaussianConditional& conditional) :
        GaussianConditional(conditional)
    {
        if(conditional.Parentsize_!= 0)
          throw std::invalid_argument("GaussianDensity can only be created from a conditional with no parents");
    }

    /// constructor using d, R
    GaussianDensity(int key,  const minimatrix& R, const minivector& d,GaussianNoiseModel* noiseModel = NULL) :
        GaussianConditional(key, R,d, noiseModel) {}

    virtual ~GaussianDensity();

    /// Construct using a mean and variance
     GaussianDensity* FromMeanAndStddev(int key, const minivector& mean, const double& sigma);


    /// Mean \f$ \mu = R^{-1} d \f$
    minivector mean() const;

    /// Covariance matrix \f$ \Sigma = (R^T R)^{-1} \f$
    minimatrix covariance() const;


};
}

#endif // GAUSSIANDENSITY_H
