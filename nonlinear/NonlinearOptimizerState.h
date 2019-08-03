#ifndef NOLINEAROPTIMIZERSTATE_H
#define NOLINEAROPTIMIZERSTATE_H


/**
 * @file NonlinearOptimizerState.h
 * @brief Private class for NonlinearOptimizer state
 * @author Richard Roberts
 * @author Frank Dellaert
 * @date Sep 7, 2009
 */


#pragma once
#include <map>
#include "../base/Matrix.h"
#include "../geometry/Pose2.h"
#include "../geometry/Pose3.h"

#include "../gmfconfig.h"
using namespace std;
namespace minisam
{
/**
 * Base class for a nonlinear optimization state, including the current estimate
 * of the variable values, error, and number of iterations.  Optimizers derived
 * from NonlinearOptimizer usually also define a derived state class containing
 * additional state specific to the algorithm (for example, Dogleg state
 * contains the current trust region radius).
 */
class NonlinearOptimizerState
{
public:
    /** The current estimate of the variable values. */

    std::map<int,Eigen::VectorXd> values;
#ifdef GMF_Using_Pose3
    std::map<int,Pose3> poses;
#else
    std::map<int,Pose2> poses;
#endif // GMF_Using_Pose3

    /** The factor graph error on the current values. */
    double error;

    /** The number of optimization iterations performed. */
    int iterations;
    double delta;

    NonlinearOptimizerState() {}
    virtual ~NonlinearOptimizerState() {}

#ifdef GMF_Using_Pose3
    NonlinearOptimizerState(const std::map<int,Eigen::VectorXd>& invalues,const std::map<int,Pose3>& inposes,
                            double error,unsigned int iterations = 0,double delta=0.0)
        : values(invalues),poses(inposes), error(error), iterations(iterations),delta(delta) {}

    // Constructor version that takes ownership of values
    NonlinearOptimizerState(std::map<int,Eigen::VectorXd>& invalues,std::map<int,Pose3>& inposes,double error,
                            unsigned int iterations = 0,double delta=0.0)
        : values(std::move(invalues)),poses(std::move(inposes)), error(error), iterations(iterations),delta(delta) {}
#else
    NonlinearOptimizerState(const std::map<int,Eigen::VectorXd>& invalues,const std::map<int,Pose2>& inposes,
                            double error, unsigned int iterations = 0,double delta=0.0)
        : values(invalues),poses(inposes), error(error), iterations(iterations),delta(delta) {}

    // Constructor version that takes ownership of values
    NonlinearOptimizerState(std::map<int,Eigen::VectorXd>& invalues,std::map<int,Pose2>& inposes,
                            double error,unsigned int iterations = 0,double delta=0.0)
        : values(std::move(invalues)),poses(std::move(inposes)), error(error), iterations(iterations),delta(delta) {}
#endif // GMF_Using_Pose3
};
#ifdef GMF_Using_Pose3
NonlinearOptimizerState* GetNonlinearOptimizerState(const std::map<int,Eigen::VectorXd>& invalues,
        const std::map<int,Pose3>& inposes,
        double error, unsigned int iterations = 0,double delta=0.0);
#else
NonlinearOptimizerState* GetNonlinearOptimizerState(const std::map<int,Eigen::VectorXd>& invalues,
        const std::map<int,Pose2>& inposes,
        double error, unsigned int iterations = 0,double delta=0.0);
#endif
};
#endif // NOLINEAROPTIMIZERSTATE_H

