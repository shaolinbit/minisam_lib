#ifndef GAUSSNEWTONOPTIMIZER_H
#define GAUSSNEWTONOPTIMIZER_H


/**
 * @file    GaussNewtonOptimizer.h
 * @brief
 * @author
 * @date
 */

#pragma once

#include "../nonlinear/NonlinearOptimizer.h"
#include "../nonlinear/NonlinearFactorGraph.h"

namespace minisam
{

class GaussNewtonOptimizer;

/** Parameters for Gauss-Newton optimization, inherits from
 * NonlinearOptimizationParams.
 */
class  GaussNewtonParams : public NonlinearOptimizerParams
{
};

/**
 * This class performs Gauss-Newton nonlinear optimization
 */
class GaussNewtonOptimizer : public NonlinearOptimizer
{

protected:
    GaussNewtonParams params_;

public:
    /// @name Standard interface
    /// @{

    /** Standard constructor, requires a nonlinear factor graph, initial
     * variable assignments, and optimization parameters.  For convenience this
     * version takes plain objects instead of shared pointers, but internally
     * copies the objects.
     * @param graph The nonlinear factor graph to optimize
     * @param initialValues The initial variable assignments
     * @param params The optimization parameters
     */
    GaussNewtonOptimizer(const NonlinearFactorGraph& graph,
                         const std::map<int,minimatrix*>& initialValues,
                         const GaussNewtonParams& params = GaussNewtonParams());

    /** Standard constructor, requires a nonlinear factor graph, initial
     * variable assignments, and optimization parameters.  For convenience this
     * version takes plain objects instead of shared pointers, but internally
     * copies the objects.
     * @param graph The nonlinear factor graph to optimize
     * @param initialValues The initial variable assignments
     */
    GaussNewtonOptimizer(const NonlinearFactorGraph& graph,
                         const std::map<int,minimatrix*>& initialValues,
                         const std::vector<int>& ordering);
    /// @}

    /// @name Advanced interface
    /// @{

    /** Virtual destructor */
    virtual ~GaussNewtonOptimizer() {}

    /** Perform a single iteration, returning a new NonlinearOptimizer class
     * containing the updated variable assignments, which may be retrieved with
     * values().
     */
    GaussianFactorGraph iterate() override;

    /** Read-only access the parameters */
    const GaussNewtonParams params() const
    {
        return params_;
    }

    /// @}

protected:
    /** Access the parameters (base class version) */
    const NonlinearOptimizerParams _params() const override
    {
        return params_;
    }

    /** Internal function for computing a COLAMD ordering if no ordering is specified */
    GaussNewtonParams ensureHasOrdering(GaussNewtonParams params,
                                        const NonlinearFactorGraph& graph) const;

};

}
#endif // GAUSSNEWTONOPTIMIZER_H

