#ifndef LEVENBERGMARQUARDTOPTIMIZER_H
#define LEVENBERGMARQUARDTOPTIMIZER_H


/**
 * @file    LevenbergMarquardtOptimizer.h
 * @brief   A nonlinear optimizer that uses the Levenberg-Marquardt trust-region scheme
 * @author  Richard Roberts
 * @author  Frank Dellaert
 * @author  Luca Carlone
 * @date    Feb 26, 2012
 */

//#pragma once

#include "../nonlinear/NonlinearOptimizer.h"
#include "../nonlinear/LevenbergMarquardtParams.h"
#include <sys/time.h>

namespace minisam
{
/**
 * This class performs Levenberg-Marquardt nonlinear optimization
 */
class  LevenbergMarquardtOptimizer: public NonlinearOptimizer
{

public:
    LevenbergMarquardtParams params_; ///< LM parameters
    struct timeval startTime_;//write it later;
    void initTime();

public:
    /// @name Constructors/Destructor
    /// @{

    /** Standard constructor, requires a nonlinear factor graph, initial
     * variable assignments, and optimization parameters.  For convenience this
     * version takes plain objects instead of shared pointers, but internally
     * copies the objects.
     * @param graph The nonlinear factor graph to optimize
     * @param initialValues The initial variable assignments
     * @param params The optimization parameters
     */

#ifdef GMF_Using_Pose3
    LevenbergMarquardtOptimizer(const NonlinearFactorGraph& graph, const std::map<int,Eigen::VectorXd>& initialValues,
                                const std::map<int,Pose3>& initialPoses,const LevenbergMarquardtParams& params = LevenbergMarquardtParams());
#else
    LevenbergMarquardtOptimizer(const NonlinearFactorGraph& graph, const std::map<int,Eigen::VectorXd>& initialValues,
                                const std::map<int,Pose2>& initialPoses, const LevenbergMarquardtParams& params = LevenbergMarquardtParams());
#endif // GMF_Using_Pose3


    /** Standard constructor, requires a nonlinear factor graph, initial
     * variable assignments, and optimization parameters.  For convenience this
     * version takes plain objects instead of shared pointers, but internally
     * copies the objects.
     * @param graph The nonlinear factor graph to optimize
     * @param initialValues The initial variable assignments
     */
#ifdef GMF_Using_Pose3
    LevenbergMarquardtOptimizer(const NonlinearFactorGraph& graph, const std::map<int,Eigen::VectorXd>& initialValues,
                                const std::map<int,Pose3>& initialPoses,
                                const std::vector<int>& ordering,
                                const LevenbergMarquardtParams& params = LevenbergMarquardtParams());
#else
    LevenbergMarquardtOptimizer(const NonlinearFactorGraph& graph, const std::map<int,Eigen::VectorXd>& initialValues,
                                const std::map<int,Pose2>& initialPoses,
                                const std::vector<int>& ordering,
                                const LevenbergMarquardtParams& params = LevenbergMarquardtParams());
#endif // GMF_Using_Pose3
    /** Virtual destructor */
    virtual ~LevenbergMarquardtOptimizer()
    {
    }

    /// @}

    /// @name Standard interface
    /// @{

    /// Access the current damping value
    double lambda() const;

    /// Access the current number of inner iterations
    int getInnerIterations() const;

    /// @}

    /// @name Advanced interface
    /// @{

    /** Perform a single iteration, returning a new NonlinearOptimizer class
     * containing the updated variable assignments, which may be retrieved with
     * values().
     */
    GaussianFactorGraph iterate() override;

    /** Read-only access the parameters */
    const LevenbergMarquardtParams params() const
    {
        return params_;
    }

    /** linearize, can be overwritten */
    virtual GaussianFactorGraph linearize() const;

    /** Build a damped system for a specific lambda -- for testing only */
    GaussianFactorGraph buildDampedSystem(const GaussianFactorGraph& linear,
                                          const std::map<int,Eigen::VectorXd>& sqrtHessianDiagonal) const;

    /** Inner loop, changes state, returns true if successful or giving up */
    bool tryLambda(const GaussianFactorGraph& linear, const std::map<int,Eigen::VectorXd>& sqrtHessianDiagonal);

    /// @}

protected:

    /** Access the parameters (base class version) */
    const NonlinearOptimizerParams _params() const override
    {
        return params_;
    }
};

};

#endif // LEVENBERGMARQUARDTOPTIMIZER_H
