#ifndef NONLINEAROPTIMIZERPARAMS_H
#define NONLINEAROPTIMIZERPARAMS_H


/**
 * @file   NonlinearOptimizerParams.h
 * @brief  Parameters for nonlinear optimization
 */


#include <string>
#include "../inference/Ordering.h"

namespace minisam
{
/** The common parameters for Nonlinear optimizers.  Most optimizers
 * deriving from NonlinearOptimizer also subclass the parameters.
 */
class  NonlinearOptimizerParams
{
public:
    /** See NonlinearOptimizerParams::verbosity */
    enum Verbosity
    {
        SILENT, TERMINATION, ERROR, VALUES, DELTA, LINEAR
    };

    int maxIterations; ///< The maximum iterations to stop iterating (default 100)
    double relativeErrorTol; ///< The maximum relative error decrease to stop iterating (default 1e-5)
    double absoluteErrorTol; ///< The maximum absolute error decrease to stop iterating (default 1e-5)
    double errorTol; ///< The maximum total error to stop iterating (default 0.0)
    Verbosity verbosity; ///< The printing verbosity during optimization (default SILENT)
    Ordering_OrderingType orderingType;///< The method of ordering use during variable elimination (default COLAMD)

    NonlinearOptimizerParams();
    virtual ~NonlinearOptimizerParams();


    int getMaxIterations() const;
    double getRelativeErrorTol() const;
    double getAbsoluteErrorTol() const;
    double getErrorTol() const;
    std::string getVerbosity() const;

    void setMaxIterations(int value);
    void setRelativeErrorTol(double value);
    void setAbsoluteErrorTol(double value);
    void setErrorTol(double value);
    void setVerbosity(const std::string &src);

    static Verbosity verbosityTranslator(const std::string &s) ;
    static std::string verbosityTranslator(Verbosity value) ;


public:

    /** See NonlinearOptimizerParams::linearSolverType */

    enum LinearSolverType
    {
        MULTIFRONTAL_CHOLESKY,
        MULTIFRONTAL_QR,
        SEQUENTIAL_CHOLESKY,
        SEQUENTIAL_QR,
        CHOLMOD, /* Experimental Flag */
    };

    LinearSolverType linearSolverType; ///< The type of linear solver to use in the nonlinear optimizer
    std::vector<int> ordering; ///< The variable elimination ordering, or empty to use COLAMD (default: empty)
    // IterativeOptimizationParameters* iterativeParams; ///< The container for iterativeOptimization parameters. used in CG Solvers.

    NonlinearOptimizerParams& operator=(const NonlinearOptimizerParams& rObj);

    inline bool isMultifrontal() const
    {
        return (linearSolverType == MULTIFRONTAL_CHOLESKY)
               || (linearSolverType == MULTIFRONTAL_QR);
    }

    inline bool isSequential() const
    {
        return (linearSolverType == SEQUENTIAL_CHOLESKY)
               || (linearSolverType == SEQUENTIAL_QR);
    }

    inline bool isCholmod() const
    {
        return (linearSolverType == CHOLMOD);
    }


    int getEliminationFunction() const;
    std::string getLinearSolverType() const;

    void setLinearSolverType(const std::string& solver);

    void setOrdering(const std::vector<int>& ordering);
    std::string getOrderingType();

    // Note that if you want to use a custom ordering, you must set the ordering directly, this will switch to custom type
    void setOrderingType(const std::string& ordering);

private:
    std::string linearSolverTranslator(LinearSolverType linearSolverType) const;

    LinearSolverType linearSolverTranslator(const std::string& linearSolverType) const;


    std::string orderingTypeTranslator(Ordering_OrderingType type);

    Ordering_OrderingType orderingTypeTranslator(const std::string& type) const;
};


inline NonlinearOptimizerParams& NonlinearOptimizerParams::operator=(const NonlinearOptimizerParams& rObj)
{
    this->maxIterations = rObj.maxIterations;
    this->relativeErrorTol = rObj.relativeErrorTol;
    this->absoluteErrorTol = rObj.absoluteErrorTol;
    this->errorTol = rObj.errorTol;
    this->verbosity = rObj.verbosity;
    this->orderingType = rObj.orderingType;
    this->linearSolverType = rObj.linearSolverType;
    this->ordering = rObj.ordering;
    return *this;
}
};

#endif // NONLINEAROPTIMIZERPARAMS_H

