#ifndef ISAM2PARAMS_H_INCLUDED
#define ISAM2PARAMS_H_INCLUDED

/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 * @file    ISAM2Params.h
 * @brief   Parameters for iSAM 2.
 * @author  Michael Kaess, Richard Roberts, Frank Dellaert
 */

#include "../nonlinear/DoglegOptimizerImpl.h"
#include "../nonlinear/NonlinearFactorGraph.h"
namespace minisam
{
/**
 * @addtogroup ISAM2
 * Parameters for ISAM2 using Gauss-Newton optimization.  Either this class or
 * ISAM2DoglegParams should be specified as the optimizationParams in
 * ISAM2Params, which should in turn be passed to ISAM2(const ISAM2Params&).
 */
struct  ISAM2GaussNewtonParams
{
    double wildfireThreshold; ///< Continue updating the linear delta only when changes are above this threshold (default: 0.001)

    /** Specify parameters as constructor arguments */
    ISAM2GaussNewtonParams(
        double _wildfireThreshold = 0.001 ///< see ISAM2GaussNewtonParams public variables, ISAM2GaussNewtonParams::wildfireThreshold
    ) : wildfireThreshold(_wildfireThreshold) {}


    double getWildfireThreshold() const
    {
        return wildfireThreshold;
    }
    void setWildfireThreshold(double wildfireThreshold)
    {
        this->wildfireThreshold = wildfireThreshold;
    }
};
struct  ISAM2DoglegParams
{
    double initialDelta; ///< The initial trust region radius for Dogleg
    double wildfireThreshold; ///< Continue updating the linear delta only when changes are above this threshold (default: 1e-5)
    DoglegOptimizerImpl::TrustRegionAdaptationMode adaptationMode; ///< See description in DoglegOptimizerImpl::TrustRegionAdaptationMode
    bool verbose; ///< Whether Dogleg prints iteration and convergence information

    /** Specify parameters as constructor arguments */
    ISAM2DoglegParams(
        double _initialDelta = 1.0, ///< see ISAM2DoglegParams::initialDelta
        double _wildfireThreshold = 1e-5, ///< see ISAM2DoglegParams::wildfireThreshold
        DoglegOptimizerImpl::TrustRegionAdaptationMode _adaptationMode = DoglegOptimizerImpl::SEARCH_EACH_ITERATION, ///< see ISAM2DoglegParams::adaptationMode
        bool _verbose = false ///< see ISAM2DoglegParams::verbose
    ) : initialDelta(_initialDelta), wildfireThreshold(_wildfireThreshold),
        adaptationMode(_adaptationMode), verbose(_verbose) {}


    double getInitialDelta() const
    {
        return initialDelta;
    }
    double getWildfireThreshold() const
    {
        return wildfireThreshold;
    }
    std::string getAdaptationMode() const
    {
        return adaptationModeTranslator(adaptationMode);
    };
    bool isVerbose() const
    {
        return verbose;
    };

    void setInitialDelta(double initialDelta)
    {
        this->initialDelta = initialDelta;
    }
    void setWildfireThreshold(double wildfireThreshold)
    {
        this->wildfireThreshold = wildfireThreshold;
    }
    void setAdaptationMode(const std::string& adaptationMode)
    {
        this->adaptationMode = adaptationModeTranslator(adaptationMode);
    }
    void setVerbose(bool verbose)
    {
        this->verbose = verbose;
    };

    std::string adaptationModeTranslator(const DoglegOptimizerImpl::TrustRegionAdaptationMode& adaptationMode) const;
    DoglegOptimizerImpl::TrustRegionAdaptationMode adaptationModeTranslator(const std::string& adaptationMode) const;
};

/**
 * @addtogroup ISAM2
 * Parameters for the ISAM2 algorithm.  Default parameter values are listed below.
 */
struct  ISAM2Params
{

    /** Optimization parameters, this both selects the nonlinear optimization
     * method and specifies its parameters, either ISAM2GaussNewtonParams or
     * ISAM2DoglegParams.  In the former, Gauss-Newton optimization will be used
     * with the specified parameters, and in the latter Powell's dog-leg
     * algorithm will be used with the specified parameters.
     */
    ISAM2GaussNewtonParams *optimizationParamsGaussNewton;
    ISAM2DoglegParams *optimizationParamsDogleg;

    /** Only relinearize variables whose linear delta magnitude is greater than
     * this threshold (default: 0.1).  If this is a FastMap<char,Vector> instead
     * of a double, then the threshold is specified for each dimension of each
     * variable type.  This parameter then maps from a character indicating the
     * variable type to a Vector of thresholds for each dimension of that
     * variable.  For example, if Pose keys are of type TypedSymbol<'x',Pose3>,
     * and landmark keys are of type TypedSymbol<'l',Point3>, then appropriate
     * entries would be added with:
     * \code
       FastMap<char,Vector> thresholds;
       thresholds['x'] = (Vector(6) << 0.1, 0.1, 0.1, 0.5, 0.5, 0.5).finished(); // 0.1 rad rotation threshold, 0.5 m translation threshold
       thresholds['l'] = Vector3(1.0, 1.0, 1.0);                // 1.0 m landmark position threshold
       params.relinearizeThreshold = thresholds;
       \endcode
     */

    double relinearizeThresholdDouble;

    std::map<char,Eigen::VectorXd> *relinearizeThresholdMap;


    int relinearizeSkip; ///< Only relinearize any variables every relinearizeSkip calls to ISAM2::update (default: 10)

    bool enableRelinearization; ///< Controls whether ISAM2 will ever relinearize any variables (default: true)

    bool evaluateNonlinearError; ///< Whether to evaluate the nonlinear error before and after the update, to return in ISAM2Result from update()

    enum Factorization { CHOLESKY, QR };

    Factorization factorization;

    /** Whether to cache linear factors (default: true).
     * This can improve performance if linearization is expensive, but can hurt
     * performance if linearization is very cleap due to computation to look up
     * additional keys.
     */
    bool cacheLinearizedFactors;

    bool enableDetailedResults; ///< Whether to compute and return ISAM2Result::detailedResults, this can increase running time (default: false)

    /** Check variables for relinearization in tree-order, stopping the check once a variable does not need to be relinearized (default: false).
     * This can improve speed by only checking a small part of the top of the tree. However, variables below the check cut-off can accumulate
     * significant deltas without triggering relinearization. This is particularly useful in exploration scenarios where real-time performance
     * is desired over correctness. Use with caution.
     */
    bool enablePartialRelinearizationCheck;

    /// When you will be removing many factors, e.g. when using ISAM2 as a fixed-lag smoother, enable this option to
    /// add factors in the first available factor slots, to avoid accumulating NULL factor slots, at the cost of
    /// having to search for slots every time a factor is added.
    bool findUnusedFactorSlots;

    /** Specify parameters as constructor arguments */
    ISAM2Params(
        ISAM2DoglegParams *_optimizationParamsDog=NULL, ///< see ISAM2Params::optimizationParams
        ISAM2GaussNewtonParams *_optimizationParamsGaussNewton=NULL,
        double _relinearizeThresholdDouble = 0.1, ///< see ISAM2Params::relinearizeThreshold
        std::map<char,Eigen::VectorXd> *_relinearizeThresholdMap=NULL,
        int _relinearizeSkip = 10, ///< see ISAM2Params::relinearizeSkip
        bool _enableRelinearization = true, ///< see ISAM2Params::enableRelinearization
        bool _evaluateNonlinearError = false, ///< see ISAM2Params::evaluateNonlinearError
        Factorization _factorization = ISAM2Params::CHOLESKY, ///< see ISAM2Params::factorization
        bool _cacheLinearizedFactors = true ///< see ISAM2Params::cacheLinearizedFactors
                                       //  const KeyFormatter& _keyFormatter = DefaultKeyFormatter ///< see ISAM2::Params::keyFormatter
    ) :optimizationParamsDogleg(_optimizationParamsDog),
        optimizationParamsGaussNewton(_optimizationParamsGaussNewton),
        relinearizeThresholdDouble(_relinearizeThresholdDouble),
        relinearizeThresholdMap(_relinearizeThresholdMap),
        relinearizeSkip(_relinearizeSkip), enableRelinearization(_enableRelinearization),
        evaluateNonlinearError(_evaluateNonlinearError), factorization(_factorization),
        cacheLinearizedFactors(_cacheLinearizedFactors),
        enableDetailedResults(false), enablePartialRelinearizationCheck(false),
        findUnusedFactorSlots(false) {}
    /// @name Getters and Setters for all properties
    /// @{

    ISAM2DoglegParams getOptimizationParamsDogleg() const
    {
        return *optimizationParamsDogleg;
    }

    ISAM2GaussNewtonParams getOptimizationParamsGaussNewton() const
    {
        return *optimizationParamsGaussNewton;
    }

    double  getrelinearizeThresholdDouble() const
    {
        return relinearizeThresholdDouble;
    }

    std::map<char,Eigen::VectorXd> getrelinearizeThresholdMap() const
    {
        return *relinearizeThresholdMap;
    }



    int getRelinearizeSkip() const
    {
        return relinearizeSkip;
    }
    bool isEnableRelinearization() const
    {
        return enableRelinearization;
    }
    bool isEvaluateNonlinearError() const
    {
        return evaluateNonlinearError;
    }
    std::string getFactorization() const
    {
        return factorizationTranslator(factorization);
    }
    bool isCacheLinearizedFactors() const
    {
        return cacheLinearizedFactors;
    }
    bool isEnableDetailedResults() const
    {
        return enableDetailedResults;
    }
    bool isEnablePartialRelinearizationCheck() const
    {
        return enablePartialRelinearizationCheck;
    }
    void setOptimizationParamsDogleg(ISAM2DoglegParams *optimizationParamsDL)
    {
        this->optimizationParamsDogleg = optimizationParamsDL;
    }
    void setOptimizationParamsGaussNewton(ISAM2GaussNewtonParams *optimizationParamsGN)
    {
        this->optimizationParamsGaussNewton = optimizationParamsGN;
    }

    void setRelinearizeThresholdDouble(double relinearizeThresholdD)
    {
        this->relinearizeThresholdDouble = relinearizeThresholdD;
    }

    void setRelinearizeThresholdMap(std::map<char,Eigen::VectorXd> *relinearizeThresholdM)
    {
        this->relinearizeThresholdMap = relinearizeThresholdM;
    }

    void setRelinearizeSkip(int relinearizeSkip)
    {
        this->relinearizeSkip = relinearizeSkip;
    }
    void setEnableRelinearization(bool enableRelinearization)
    {
        this->enableRelinearization = enableRelinearization;
    }
    void setEvaluateNonlinearError(bool evaluateNonlinearError)
    {
        this->evaluateNonlinearError = evaluateNonlinearError;
    }
    void setFactorization(const std::string& factorization)
    {
        this->factorization = factorizationTranslator(factorization);
    }
    void setCacheLinearizedFactors(bool cacheLinearizedFactors)
    {
        this->cacheLinearizedFactors = cacheLinearizedFactors;
    }

    void setEnableDetailedResults(bool enableDetailedResults)
    {
        this->enableDetailedResults = enableDetailedResults;
    }
    void setEnablePartialRelinearizationCheck(bool enablePartialRelinearizationCheck)
    {
        this->enablePartialRelinearizationCheck = enablePartialRelinearizationCheck;
    }


    int  getEliminationFunction() const
    {
        return factorization == CHOLESKY ? 1 : 0;
    }

    /// @}

    /// @name Some utilities
    /// @{

    static Factorization factorizationTranslator(const std::string& str);
    static std::string factorizationTranslator(const Factorization& value);

    /// @}
};

/**
 * @addtogroup ISAM2
 * This struct is returned from ISAM2::update() and contains information about
 * the update that is useful for determining whether the solution is
 * converging, and about how much work was required for the update.  See member
 * variables for details and information about each entry.
 */
struct ISAM2Result
{
    /** The nonlinear error of all of the factors, \a including new factors and
     * variables added during the current call to ISAM2::update().  This error is
     * calculated using the following variable values:
     * \li Pre-existing variables will be evaluated by combining their
     * linearization point before this call to update, with their partial linear
     * delta, as computed by ISAM2::calculateEstimate().
     * \li New variables will be evaluated at their initialization points passed
     * into the current call to update.
     * \par Note: This will only be computed if ISAM2Params::evaluateNonlinearError
     * is set to \c true, because there is some cost to this computation.
     */
    double errorBefore;

    /** The nonlinear error of all of the factors computed after the current
     * update, meaning that variables above the relinearization threshold
     * (ISAM2Params::relinearizeThreshold) have been relinearized and new
     * variables have undergone one linear update.  Variable values are
     * again computed by combining their linearization points with their
     * partial linear deltas, by ISAM2::calculateEstimate().
     * \par Note: This will only be computed if ISAM2Params::evaluateNonlinearError
     * is set to \c true, because there is some cost to this computation.
     */
    double errorAfter;

    /** The number of variables that were relinearized because their linear
     * deltas exceeded the reslinearization threshold
     * (ISAM2Params::relinearizeThreshold), combined with any additional
     * variables that had to be relinearized because they were involved in
     * the same factor as a variable above the relinearization threshold.
     * On steps where no relinearization is considered
     * (see ISAM2Params::relinearizeSkip), this count will be zero.
     */
    int variablesRelinearized;

    /** The number of variables that were reeliminated as parts of the Bayes'
     * Tree were recalculated, due to new factors.  When loop closures occur,
     * this count will be large as the new loop-closing factors will tend to
     * involve variables far away from the root, and everything up to the root
     * will be reeliminated.
     */
    int variablesReeliminated;

    /** The number of factors that were included in reelimination of the Bayes' tree. */
    int factorsRecalculated;

    /** The number of cliques in the Bayes' Tree */
    int cliques;

    /** The indices of the newly-added factors, in 1-to-1 correspondence with the
     * factors passed as \c newFactors to ISAM2::update().  These indices may be
     * used later to refer to the factors in order to remove them.
     */
    std::vector<int> newFactorsIndices;

    /** A struct holding detailed results, which must be enabled with
     * ISAM2Params::enableDetailedResults.
     */
    struct DetailedResults
    {
        /** The status of a single variable, this struct is stored in
         * DetailedResults::variableStatus */
        struct VariableStatus
        {
            /** Whether the variable was just reeliminated, due to being relinearized,
             * observed, new, or on the path up to the root clique from another
             * reeliminated variable. */
            bool isReeliminated;
            bool isAboveRelinThreshold; ///< Whether the variable was just relinearized due to being above the relinearization threshold
            bool isRelinearizeInvolved; ///< Whether the variable was below the relinearization threshold but was relinearized by being involved in a factor with a variable above the relinearization threshold
            bool isRelinearized; /// Whether the variable was relinearized, either by being above the relinearization threshold or by involvement.
            bool isObserved; ///< Whether the variable was just involved in new factors
            bool isNew; ///< Whether the variable itself was just added
            bool inRootClique; ///< Whether the variable is in the root clique
            VariableStatus(): isReeliminated(false), isAboveRelinThreshold(false), isRelinearizeInvolved(false),
                isRelinearized(false), isObserved(false), isNew(false), inRootClique(false) {}
        };

        /** The status of each variable during this update, see VariableStatus.
         */
        std::map<int, VariableStatus> variableStatus;
    };

    /** Detailed results, if enabled by ISAM2Params::enableDetailedResults.  See
     * Detail for information about the results data stored here. */
    DetailedResults *detail;

    /** Getters and Setters */
    int getVariablesRelinearized() const
    {
        return variablesRelinearized;
    };
    int getVariablesReeliminated() const
    {
        return variablesReeliminated;
    };
    int getCliques() const
    {
        return cliques;
    };
};
};

#endif // ISAM2PARAMS_H_INCLUDED
