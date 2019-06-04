#ifndef DOGLEGOPITIMIZER_H
#define DOGLEGOPITIMIZER_H


/**
 * @file    DoglegOptimizer.h
 * @brief
 * @author
 * @date
 */

//#pragma once

#include "../nonlinear/NonlinearOptimizer.h"


class DoglegOptimizer;

/** Parameters for Levenberg-Marquardt optimization.  Note that this parameters
 * class inherits from NonlinearOptimizerParams, which specifies the parameters
 * common to all nonlinear optimization algorithms.  This class also contains
 * all of those parameters.
 */
class  DoglegParams : public NonlinearOptimizerParams
{
public:
    /** See DoglegParams::dlVerbosity */
    enum VerbosityDL
    {
        SILENT,
        VERBOSE
    };

    double deltaInitial; ///< The initial trust region radius (default: 1.0)
    VerbosityDL verbosityDL; ///< The verbosity level for Dogleg (default: SILENT), see also NonlinearOptimizerParams::verbosity

    DoglegParams() :
        deltaInitial(1.0), verbosityDL(SILENT) {}

    virtual ~DoglegParams() {}
    /*
    void print(const std::string& str = "") const override {
      NonlinearOptimizerParams::print(str);
      std::cout << "               deltaInitial: " << deltaInitial << "\n";
      std::cout.flush();
    }*/

    double getDeltaInitial() const
    {
        return deltaInitial;
    }
    std::string getVerbosityDL() const
    {
        return verbosityDLTranslator(verbosityDL);
    }

    void setDeltaInitial(double deltaInitial)
    {
        this->deltaInitial = deltaInitial;
    }
    void setVerbosityDL(const std::string& verbosityDL)
    {
        this->verbosityDL = verbosityDLTranslator(verbosityDL);
    }

private:
    VerbosityDL verbosityDLTranslator(const std::string& verbosityDL) const;
    std::string verbosityDLTranslator(VerbosityDL verbosityDL) const;
};

/**
 * This class performs Dogleg nonlinear optimization
 */
class  DoglegOptimizer : public NonlinearOptimizer
{

protected:
    DoglegParams params_;

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
#ifdef GMF_Using_Pose3
    DoglegOptimizer(const NonlinearFactorGraph& graph, const std::map<int,Eigen::VectorXd>& initialValues,
                    const std::map<int,Pose3>& initialPoses,
                    const DoglegParams& params = DoglegParams());
#else
    DoglegOptimizer(const NonlinearFactorGraph& graph, const std::map<int,Eigen::VectorXd>& initialValues,
                    const std::map<int,Pose2>& initialPoses,
                    const DoglegParams& params = DoglegParams());
#endif
    /** Standard constructor, requires a nonlinear factor graph, initial
     * variable assignments, and optimization parameters.  For convenience this
     * version takes plain objects instead of shared pointers, but internally
     * copies the objects.
     * @param graph The nonlinear factor graph to optimize
     * @param initialValues The initial variable assignments
     */
#ifdef GMF_Using_Pose3
    DoglegOptimizer(const NonlinearFactorGraph& graph, const std::map<int,Eigen::VectorXd>& initialValues,
                    const std::map<int,Pose3>& initialPoses,
                    const Ordering& ordering);
#else
    DoglegOptimizer(const NonlinearFactorGraph& graph, const std::map<int,Eigen::VectorXd>& initialValues,
                    const std::map<int,Pose2>& initialPoses,
                    const Ordering& ordering);
#endif

    /// @}

    /// @name Advanced interface
    /// @{

    /** Virtual destructor */
    virtual ~DoglegOptimizer() {}

    /** Perform a single iteration, returning a new NonlinearOptimizer class
     * containing the updated variable assignments, which may be retrieved with
     * values().
     */
    GaussianFactorGraph iterate() override;

    /** Read-only access the parameters */
    const DoglegParams params() const
    {
        return params_;
    }

    /** Access the current trust region radius delta */
    double getDelta() const;

    /// @}

protected:
    /** Access the parameters (base class version) */
    virtual const NonlinearOptimizerParams _params() const override
    {
        return params_;
    }

    /** Internal function for computing a COLAMD ordering if no ordering is specified */
    DoglegParams ensureHasOrdering(DoglegParams params, const NonlinearFactorGraph& graph) const;
};

//}
#endif // DOGLEGOPITIMIZER_H
