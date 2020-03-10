#ifndef LEVENBERGMARQUARDTSTATE_H
#define LEVENBERGMARQUARDTSTATE_H


/**
 * @file    LevenbergMarquardtState.h
 * @brief   A LevenbergMarquardtState class containing most of the logic for Levenberg-Marquardt
 */

#include "../nonlinear/NonlinearOptimizerState.h"
#include "../linear/GaussianFactorGraph.h"
#include "../linear/JacobianFactor.h"
#include "../linear/NoiseModel.h"
#include "../nonlinear/LevenbergMarquardtParams.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>


/** Small struct to cache objects needed for damping.
  * This is used in buildDampedSystem  */
namespace minisam
{
struct CachedModel
{
    CachedModel()
    {
        model=NULL;   // default int makes zero-size matrices
        A.data=NULL;
        A.owner=0;
        A.prd=0;
        A.size1=0;
        A.size2=0;
        b.data=NULL;
        b.owner=0;
        b.prd=0;
        b.size1=0;
        b.size2=0;
    }
    CachedModel(int dim, double sigma)
        : A(minimatrix_identity_mat(dim)),
          b(minivector(dim,0.0)),//b(Eigen::VectorXd::Zero(dim)),
          model(new IsotropicNoiseModel(dim, sigma)) {}
    CachedModel(int dim, double sigma, const minivector& diagonal)
        :A(minimatrix_vector_asDiagonal(diagonal)),// A(Eigen::DiagonalMatrix<double, Eigen::Dynamic>(diagonal)),
         b(minivector(dim,0.0)),// b(Eigen::VectorXd::Zero(dim)),
         model(new IsotropicNoiseModel(dim, sigma)) {}
    minimatrix A;
    minivector b;
    GaussianNoiseModel* model;
    ~CachedModel()
    {
        if(model!=NULL)
        {
            delete model;
            model=NULL;
        }
    }
};

class LevenbergMarquardtState : public NonlinearOptimizerState
{

public:

  //  double lambda;
  //  double currentFactor;
 //   int totalNumberInnerIterations;  ///< The total number of inner iterations in the
    // optimization (for each iteration, LM tries multiple
    // inner iterations with different lambdas)
    LevenbergMarquardtState(const std::map<int,minimatrix*>& initialValues,
                            double error,double lambda, double currentFactor,
                            double delta=0.0,
                            unsigned int iterations = 0, unsigned int totalNumberInnerIterations = 0)
        : NonlinearOptimizerState(initialValues, error, iterations,delta,lambda,currentFactor,totalNumberInnerIterations)
        {


        } // lambda(lambda),
        //  currentFactor(currentFactor),
        //  totalNumberInnerIterations(totalNumberInnerIterations) {}

    // Applies policy to *increase* lambda: should be used if the current update was NOT successful
    void increaseLambda(const LevenbergMarquardtParams& params);

    // Apply policy to decrease lambda if the current update was successful
    // stepQuality not used in the naive policy)
    // Take ownsership of newValues, must be passed an rvalue
    LevenbergMarquardtState* decreaseLambda(const LevenbergMarquardtParams& params, double stepQuality,
                                            std::map<int,minimatrix*>& newValues,
                                            double newError) const;
    void decreaseLambda(const LevenbergMarquardtParams& params, double stepQuality,
                                            std::map<int,minimatrix*>& newValues,
                                            double newError,NonlinearOptimizerState* setstate) const;

    // Small cache of A|b|model indexed by dimension. Can save many mallocs.
    mutable std::vector<CachedModel*> noiseModelCache;

    CachedModel* getCachedModel(int dim) const ;

    virtual ~LevenbergMarquardtState();

    /// Build a damped system for a specific lambda, vanilla version
    GaussianFactorGraph buildDampedSystem(GaussianFactorGraph damped /* gets copied */) const ;

    /// Build a damped system, use hessianDiagonal per variable (more expensive)
    GaussianFactorGraph buildDampedSystem(GaussianFactorGraph damped,  // gets copied
                                          const std::map<int,minivector>& sqrtHessianDiagonal) const;
};
};

#endif // LEVENBERGMARQUARDTSTATE_H
