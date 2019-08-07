#ifndef LEVENBERGMARQUARDTSTATE_H
#define LEVENBERGMARQUARDTSTATE_H


/**
 * @file    LevenbergMarquardtState.h
 * @brief   A LevenbergMarquardtState class containing most of the logic for Levenberg-Marquardt
 * @author  Frank Dellaert
 * @date    April 2016
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
    }
    CachedModel(int dim, double sigma)
        : A(Eigen::MatrixXd::Identity(dim, dim)),
          b(Eigen::VectorXd::Zero(dim)),
          model(new IsotropicNoiseModel(dim, sigma)) {}
    CachedModel(int dim, double sigma, const Eigen::VectorXd& diagonal)
        : A(Eigen::DiagonalMatrix<double, Eigen::Dynamic>(diagonal)),
          b(Eigen::VectorXd::Zero(dim)),
          model(new IsotropicNoiseModel(dim, sigma)) {}
    Eigen::MatrixXd A;
    Eigen::VectorXd b;
    DiagonalNoiseModel* model;
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

    double lambda;
    double currentFactor;
    int totalNumberInnerIterations;  ///< The total number of inner iterations in the
    // optimization (for each iteration, LM tries multiple
    // inner iterations with different lambdas)
#ifdef GMF_Using_Pose3
    // Constructor version that takes ownership of values
    LevenbergMarquardtState(const std::map<int,Eigen::VectorXd>& initialValues,const std::map<int,Pose3>& initialposes, double error, double lambda, double currentFactor,
                            unsigned int iterations = 0, unsigned int totalNumberInnerIterations = 0)
        : NonlinearOptimizerState(initialValues,initialposes, error, iterations),
          lambda(lambda),
          currentFactor(currentFactor),
          totalNumberInnerIterations(totalNumberInnerIterations) {}
#else

    // Constructor version that takes ownership of values
    LevenbergMarquardtState(const std::map<int,Eigen::VectorXd>& initialValues,const std::map<int,Pose2>& initialposes, double error, double lambda, double currentFactor,
                            unsigned int iterations = 0, unsigned int totalNumberInnerIterations = 0)
        : NonlinearOptimizerState(initialValues,initialposes, error, iterations),
          lambda(lambda),
          currentFactor(currentFactor),
          totalNumberInnerIterations(totalNumberInnerIterations) {}
#endif // GMF_Using_Pose3



    // Applies policy to *increase* lambda: should be used if the current update was NOT successful
    void increaseLambda(const LevenbergMarquardtParams& params);

    // Apply policy to decrease lambda if the current update was successful
    // stepQuality not used in the naive policy)
    // Take ownsership of newValues, must be passed an rvalue
#ifdef GMF_Using_Pose3
    LevenbergMarquardtState* decreaseLambda(const LevenbergMarquardtParams& params, double stepQuality,
                                            std::map<int,Eigen::VectorXd>& newValues,
                                            std::map<int,Pose3>& newPoses,
                                            double newError) const;
#else
    LevenbergMarquardtState* decreaseLambda(const LevenbergMarquardtParams& params, double stepQuality,
                                            std::map<int,Eigen::VectorXd>& newValues,
                                            std::map<int,Pose2>& newPoses,
                                            double newError) const;
#endif // GMF_Using_Pose3



    // Small cache of A|b|model indexed by dimension. Can save many mallocs.
    mutable std::vector<CachedModel*> noiseModelCache;

    CachedModel* getCachedModel(int dim) const ;

    virtual ~LevenbergMarquardtState();

    /// Build a damped system for a specific lambda, vanilla version
    GaussianFactorGraph buildDampedSystem(GaussianFactorGraph damped /* gets copied */) const ;

    /// Build a damped system, use hessianDiagonal per variable (more expensive)
    GaussianFactorGraph buildDampedSystem(GaussianFactorGraph damped,  // gets copied
                                          const std::map<int,Eigen::VectorXd>& sqrtHessianDiagonal) const;
};
};

#endif // LEVENBERGMARQUARDTSTATE_H
