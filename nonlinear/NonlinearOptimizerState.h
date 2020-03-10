#ifndef NOLINEAROPTIMIZERSTATE_H
#define NOLINEAROPTIMIZERSTATE_H


/**
 * @file NonlinearOptimizerState.h
 * @brief Private class for NonlinearOptimizer state
 */


#pragma once
#include <map>
#include "../mat/Matrix.h"
#include "../geometry/Pose2.h"
#include "../geometry/Pose3.h"

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

    std::map<int,minimatrix*> values;


    /** The factor graph error on the current values. */
    double error;

    /** The number of optimization iterations performed. */
    int iterations;
    double delta;
    double lambda;
    double currentFactor;
    int totalNumberInnerIterations;  ///< The total number of inner iterations in the


    NonlinearOptimizerState() {}
    virtual ~NonlinearOptimizerState()
    {
    /*
     for(auto& kst:values)
     {
       if(kst.second!=NULL)
       {
         delete kst.second;
         kst.second=NULL;
       }
     }
     values.clear();
     */
    }
    NonlinearOptimizerState(const std::map<int,minimatrix*>& invalues,
                            double error,unsigned int iterations = 0,double delta=0.0,
                           double inlambda=0.0,
    double incurrentFactor=0.0,
    int intotalNumberInnerIterations=0)
        : values(invalues), error(error), iterations(iterations),delta(delta),
        lambda(inlambda),currentFactor(incurrentFactor),totalNumberInnerIterations(intotalNumberInnerIterations) {}

    // Constructor version that takes ownership of values
    NonlinearOptimizerState(std::map<int,minimatrix*>& invalues,double error,
                            unsigned int iterations = 0,double delta=0.0,double inlambda=0.0,
    double incurrentFactor=0.0,
    int intotalNumberInnerIterations=0)
        : values(std::move(invalues)), error(error), iterations(iterations),delta(delta),
        lambda(inlambda),currentFactor(incurrentFactor),totalNumberInnerIterations(intotalNumberInnerIterations) {}

};
NonlinearOptimizerState* GetNonlinearOptimizerState(const std::map<int,minimatrix*>& invalues,
        double error, unsigned int iterations = 0,double delta=0.0);

};
#endif // NOLINEAROPTIMIZERSTATE_H

