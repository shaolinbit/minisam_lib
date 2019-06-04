#ifndef NONLINEAROPTIMIZERPARAMS_H
#define NONLINEAROPTIMIZERPARAMS_H


/**
 * @file   NonlinearOptimizerParams.h
 * @brief  Parameters for nonlinear optimization
 * @author
 * @date
 */


//#include "../linear/GaussianFactorGraph.h"
//#include "../linear/SubgraphSolver.h"
//#include <boost/optional.hpp>
#include <string>
#include "../inference/Ordering.h"


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
    Ordering::OrderingType orderingType; ///< The method of ordering use during variable elimination (default COLAMD)

    NonlinearOptimizerParams() :
        maxIterations(100), relativeErrorTol(1e-5), absoluteErrorTol(1e-5), errorTol(
            0.0), verbosity(SILENT), orderingType(Ordering::COLAMD),
        linearSolverType(MULTIFRONTAL_CHOLESKY) {}

    virtual ~NonlinearOptimizerParams()
    {
    }
// virtual void print(const std::string& str = "") const;

    int getMaxIterations() const
    {
        return maxIterations;
    }
    double getRelativeErrorTol() const
    {
        return relativeErrorTol;
    }
    double getAbsoluteErrorTol() const
    {
        return absoluteErrorTol;
    }
    double getErrorTol() const
    {
        return errorTol;
    }
    std::string getVerbosity() const
    {
        return verbosityTranslator(verbosity);
    }

    void setMaxIterations(int value)
    {
        maxIterations = value;
    }
    void setRelativeErrorTol(double value)
    {
        relativeErrorTol = value;
    }
    void setAbsoluteErrorTol(double value)
    {
        absoluteErrorTol = value;
    }
    void setErrorTol(double value)
    {
        errorTol = value;
    }
    void setVerbosity(const std::string &src)
    {
        verbosity = verbosityTranslator(src);
    }

    static Verbosity verbosityTranslator(const std::string &s) ;
    static std::string verbosityTranslator(Verbosity value) ;

    // Successive Linearization Parameters

public:

    /** See NonlinearOptimizerParams::linearSolverType */

    enum LinearSolverType
    {
        MULTIFRONTAL_CHOLESKY,
        MULTIFRONTAL_QR,
        SEQUENTIAL_CHOLESKY,
        SEQUENTIAL_QR,
       // Iterative, /* Experimental Flag */
        CHOLMOD, /* Experimental Flag */
    };

    LinearSolverType linearSolverType; ///< The type of linear solver to use in the nonlinear optimizer
    Ordering ordering; ///< The variable elimination ordering, or empty to use COLAMD (default: empty)

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

    /*
    inline bool isIterative() const
    {
        return (linearSolverType == Iterative);
    }*/

    int getEliminationFunction() const
    {
        switch (linearSolverType)
        {
        case MULTIFRONTAL_CHOLESKY:
        case SEQUENTIAL_CHOLESKY:
            return 1;

        case MULTIFRONTAL_QR:
        case SEQUENTIAL_QR:
            return 0;

        default:
            throw std::runtime_error(
                "Nonlinear optimization parameter \"factorization\" is invalid");
        }
    }

    std::string getLinearSolverType() const
    {
        return linearSolverTranslator(linearSolverType);
    }

    void setLinearSolverType(const std::string& solver)
    {
        linearSolverType = linearSolverTranslator(solver);
    }

   // void setIterativeParams(const IterativeOptimizationParameters& params);

    void setOrdering(const Ordering& ordering)
    {
        this->ordering = ordering;
        this->orderingType = Ordering::CUSTOM;
    }

    std::string getOrderingType()
    {
        return orderingTypeTranslator(orderingType);
    }

    // Note that if you want to use a custom ordering, you must set the ordering directly, this will switch to custom type
    void setOrderingType(const std::string& ordering)
    {
        orderingType = orderingTypeTranslator(ordering);
    }

private:
    std::string linearSolverTranslator(LinearSolverType linearSolverType) const;

    LinearSolverType linearSolverTranslator(const std::string& linearSolverType) const;


    std::string orderingTypeTranslator(Ordering::OrderingType type);

    Ordering::OrderingType orderingTypeTranslator(const std::string& type) const;

};

// For backward compatibility:
//typedef NonlinearOptimizerParams SuccessiveLinearizationParams;

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
   // this->iterativeParams = rObj.iterativeParams;
    return *this;
}


#endif // NONLINEAROPTIMIZERPARAMS_H

