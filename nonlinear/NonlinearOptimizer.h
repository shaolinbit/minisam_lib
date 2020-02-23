#ifndef NONLINEAROPTIMIZER_H
#define NONLINEAROPTIMIZER_H


/**
 * @file NonlinearOptimizer.h
 * @brief Base class and parameters for nonlinear optimization algorithms
 */

#pragma once

#include "../nonlinear/NonlinearFactorGraph.h"
#include "../nonlinear/NonlinearOptimizerParams.h"
#include "../nonlinear/NonlinearOptimizerState.h"


namespace minisam
{
struct NonlinearOptimizerState;
/**
 * This is the abstract interface for classes that can optimize for the
 * maximum-likelihood estimate of a NonlinearFactorGraph.
 *
 * To use a class derived from this interface, construct the class with a
 * NonlinearFactorGraph and an initial Values variable assignment.  Next, call the
 * optimize() method which returns the optimized variable assignment.
 *
 * Simple and compact example:
 * \code
// One-liner to do full optimization and use the result.
Values result = DoglegOptimizer(graph, initialValues).optimize();
\endcode
 *
 * Example exposing more functionality and details:
 * \code
// Create initial optimizer
DoglegOptimizer optimizer(graph, initialValues);

// Run full optimization until convergence.
Values result = optimizer->optimize();

// The new optimizer has results and statistics
cout << "Converged in " << optimizer.iterations() << " iterations "
        "with final error " << optimizer.error() << endl;
\endcode
 *
 * Example of setting parameters before optimization:
 * \code
// Each derived optimizer type has its own parameters class, which inherits from NonlinearOptimizerParams
DoglegParams params;
params.factorization = DoglegParams::QR;
params.relativeErrorTol = 1e-3;
params.absoluteErrorTol = 1e-3;

// Optimize
Values result = DoglegOptimizer(graph, initialValues, params).optimize();
\endcode
 *
 * This interface also exposes an iterate() method, which performs one
 * iteration.  The optimize() method simply calls iterate() multiple times,
 * until the error changes less than a threshold.  We expose iterate() so that
 * you can easily control what happens between iterations, such as drawing or
 * printing, moving points from behind the camera to in front, etc.
 *
 * For more flexibility you may override virtual methods in your own derived class.
 */

class NonlinearOptimizer
{

protected:
    NonlinearFactorGraph graph_; ///< The graph with nonlinear factors

    NonlinearOptimizerState* state_; ///< PIMPL'd state


public:

    NonlinearOptimizer();//{}

    /// @name Standard interface
    /// @{

    /** Optimize for the maximum-likelihood estimate, returning a new
     * NonlinearOptimizer class containing the optimized variable assignments,
     * which may be retrieved with values().
     *
     * This function simply calls iterate() in a loop, checking for convergence
     * with check_convergence().  For fine-grain control over the optimization
     * process, you may call iterate() and check_convergence() yourself, and if
     * needed modify the optimization state between iterations.
     */
    virtual const std::map<int,minimatrix*> optimize();

    /**
     * Optimize, but return empty result if any uncaught exception is thrown
     * Intended for MATLAB. In C++, use above and catch exceptions.
     * No message is printed: it is up to the caller to check the result
     * @param optimizer a non-linear optimizer
     */
    const std::map<int,minimatrix*> optimizeSafely();

    /// return error
    double error() const;

    /// return number of iterations
    int iterations() const;

    /// return values
    const std::map<int,minimatrix*> values() const;

    /// @}

    /// @name Advanced interface
    /// @{

    /** Virtual destructor */
    virtual ~NonlinearOptimizer();

    /** Default function to do linear solve, i.e. optimize a GaussianFactorGraph */
    virtual std::map<int,minivector> solve(const GaussianFactorGraph &gfg,
                                           const NonlinearOptimizerParams& params) const;

    /** Perform a single iteration, returning a new NonlinearOptimizer class
     * containing the updated variable assignments, which may be retrieved with
     * values().
     */
    virtual GaussianFactorGraph iterate() = 0;

    /// @}

protected:
    /** A default implementation of the optimization loop, which calls iterate()
     * until checkConvergence returns true.
     */
    void defaultOptimize();

    virtual const NonlinearOptimizerParams _params() const = 0;

    /** Constructor for initial construction of base classes. Takes ownership of state. */
    NonlinearOptimizer(const NonlinearFactorGraph& graph,
                       NonlinearOptimizerState* state);
};

/** Check whether the relative error decrease is less than relativeErrorTreshold,
 * the absolute error decrease is less than absoluteErrorTreshold, <em>or</em>
 * the error itself is less than errorThreshold.
 */
bool checkConvergence(double relativeErrorTreshold,
                      double absoluteErrorTreshold, double errorThreshold,
                      double currentError, double newError, NonlinearOptimizerParams::Verbosity verbosity = NonlinearOptimizerParams::SILENT);

bool checkConvergence(const NonlinearOptimizerParams& params, double currentError,
                      double newError);

};
#endif // NONLINEAROPTIMIZER_H
