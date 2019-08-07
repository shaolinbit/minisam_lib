#ifndef DOGLEGOPITIMIZERIMPL_H
#define DOGLEGOPITIMIZERIMPL_H

/**
 * @file    DoglegOptimizerImpl.h
 * @brief   Nonlinear factor graph optimizer using Powell's Dogleg algorithm (detail implementation)
 * @author  Richard Roberts
 */

#include <iomanip>
#include "../nonlinear/NonlinearFactorGraph.h"
#include "../inference/BayesTree.h"
#include "../linear/GaussianBayesNet.h"
namespace minisam
{
class ISAM2;

/** This class contains the implementation of the Dogleg algorithm.  It is used
 * by DoglegOptimizer and can be used to easily put together custom versions of
 * Dogleg.  Each function is well-documented and unit-tested.  The notation
 * here matches that in "trustregion.pdf" in doc, see this file for further
 * explanation of the computations performed by this class.
 */
struct  DoglegOptimizerImpl
{

    struct  IterationResult
    {
        double delta;
        std::map<int,Eigen::VectorXd> dx_d;
        double f_error;
    };

    /** Specifies how the trust region is adapted at each Dogleg iteration.  If
     * this is SEARCH_EACH_ITERATION, then the trust region radius will be
     * increased potentially multiple times during one iteration until increasing
     * it further no longer decreases the error.  If this is
     * ONE_STEP_PER_ITERATION, then the step in one iteration will not exceed the
     * current trust region radius, but the radius will be increased for the next
     * iteration if the error decrease is good.  The former will generally result
     * in slower iterations, but sometimes larger steps in early iterations.  The
     * latter generally results in faster iterations but it may take several
     * iterations before the trust region radius is increased to the optimal
     * value.  Generally ONE_STEP_PER_ITERATION should be used, corresponding to
     * most published descriptions of the algorithm.
     */
    enum TrustRegionAdaptationMode
    {
        SEARCH_EACH_ITERATION,
        SEARCH_REDUCE_ONLY,
        ONE_STEP_PER_ITERATION
    };

    /**
     * Compute the update point for one iteration of the Dogleg algorithm, given
     * an initial trust region radius \f$ \delta \f$.  The trust region radius is
     * adapted based on the error of a NonlinearFactorGraph \f$ f(x) \f$, and
     * depending on the update mode \c mode.
     *
     * The update is computed using a quadratic approximation \f$ M(\delta x) \f$
     * of an original nonlinear error function (a NonlinearFactorGraph) \f$ f(x) \f$.
     * The quadratic approximation is represented as a GaussianBayesNet \f$ \bayesNet \f$, which is
     * obtained by eliminating a GaussianFactorGraph resulting from linearizing
     * the nonlinear factor graph \f$ f(x) \f$.  Thus, \f$ M(\delta x) \f$ is
     * \f[
     * M(\delta x) = (R \delta x - d)^T (R \delta x - d)
     * \f]
     * where \f$ R \f$ and \f$ d \f$ together are a Bayes' net or Bayes' tree.
     * \f$ R \f$ is upper-triangular and \f$ d \f$ is a vector, represented
     * as a BayesNet<GaussianConditional> (GaussianBayesNet) or
     * BayesTree<GaussianConditional>, containing GaussianConditional s.
     *
     * @tparam M The type of the Bayes' net or tree, currently
     * either BayesNet<GaussianConditional> (or GaussianBayesNet) or BayesTree<GaussianConditional>.
     * @tparam F For normal usage this will be NonlinearFactorGraph<VALUES>.
     * @tparam VALUES The Values or TupleValues to pass to F::error() to evaluate
     * the error function.
     * @param delta The initial trust region radius.
     * @param mode See DoglegOptimizerImpl::TrustRegionAdaptationMode
     * @param Rd The Bayes' net or tree as described above.
     * @param f The original nonlinear factor graph with which to evaluate the
     * accuracy of \f$ M(\delta x) \f$ to adjust \f$ \delta \f$.
     * @param x0 The linearization point about which \f$ \bayesNet \f$ was created
     * @param ordering The variable ordering used to create\f$ \bayesNet \f$
     * @param f_error The result of <tt>f.error(x0)</tt>.
     * @return A DoglegIterationResult containing the new \c delta, the linear
     * update \c dx_d, and the resulting nonlinear error \c f_error.
     */
#ifdef GMF_Using_Pose3

    static IterationResult Iterate(
        double delta, TrustRegionAdaptationMode mode, const std::map<int,Eigen::VectorXd>& dx_u,
        const std::map<int,Eigen::VectorXd>& dx_n,
        const BayesTree& Rd, const NonlinearFactorGraph& f, const std::map<int,Eigen::VectorXd>& x0,
        const std::map<int,Pose3>& xpose0,const double f_error, const bool verbose=false);
    static IterationResult Iterate(
        double delta, TrustRegionAdaptationMode mode, const std::map<int,Eigen::VectorXd>& dx_u,
        const std::map<int,Eigen::VectorXd>& dx_n,
        const GaussianBayesNet& Rd, const NonlinearFactorGraph& f, const std::map<int,Eigen::VectorXd>& x0,
        const std::map<int,Pose3>& xpose0, const double f_error, const bool verbose=false);
    static IterationResult Iterate(
        double delta, TrustRegionAdaptationMode mode, const std::map<int,Eigen::VectorXd>& dx_u,
        const std::map<int,Eigen::VectorXd>& dx_n,
        const ISAM2& Rd, const NonlinearFactorGraph& f, const std::map<int,Eigen::VectorXd>& x0,
        const std::map<int,Pose3>& xpose0, const double f_error, const bool verbose=false);
#else

    static IterationResult Iterate(
        double delta, TrustRegionAdaptationMode mode, const std::map<int,Eigen::VectorXd>& dx_u,
        const std::map<int,Eigen::VectorXd>& dx_n,
        const BayesTree& Rd, const NonlinearFactorGraph& f, const std::map<int,Eigen::VectorXd>& x0,
        const std::map<int,Pose2>& xpose0,const double f_error, const bool verbose=false);


    static IterationResult Iterate(
        double delta, TrustRegionAdaptationMode mode, const std::map<int,Eigen::VectorXd>& dx_u,
        const std::map<int,Eigen::VectorXd>& dx_n,
        const GaussianBayesNet& Rd, const NonlinearFactorGraph& f, const std::map<int,Eigen::VectorXd>& x0,
        const std::map<int,Pose2>& xpose0, const double f_error, const bool verbose=false);
    static IterationResult Iterate(
        double delta, TrustRegionAdaptationMode mode, const std::map<int,Eigen::VectorXd>& dx_u,
        const std::map<int,Eigen::VectorXd>& dx_n,
        const ISAM2& Rd, const NonlinearFactorGraph& f, const std::map<int,Eigen::VectorXd>& x0,
        const std::map<int,Pose2>& xpose0, const double f_error, const bool verbose=false);
#endif
    /**
     * Compute the dogleg point given a trust region radius \f$ \delta \f$.  The
     * dogleg point is the intersection between the dogleg path and the trust
     * region boundary, see doc/trustregion.pdf for more details.
     *
     * The update is computed using a quadratic approximation \f$ M(\delta x) \f$
     * of an original nonlinear error function (a NonlinearFactorGraph) \f$ f(x) \f$.
     * The quadratic approximation is represented as a GaussianBayesNet \f$ \bayesNet \f$, which is
     * obtained by eliminating a GaussianFactorGraph resulting from linearizing
     * the nonlinear factor graph \f$ f(x) \f$.  Thus, \f$ M(\delta x) \f$ is
     * \f[
     * M(\delta x) = (R \delta x - d)^T (R \delta x - d)
     * \f]
     * where \f$ R \f$ and \f$ d \f$ together are a Bayes' net.  \f$ R \f$ is
     * upper-triangular and \f$ d \f$ is a vector, represented as a
     * GaussianBayesNet, containing GaussianConditional s.
     *
     * @param delta The trust region radius
     * @param dx_u The steepest descent point, i.e. the Cauchy point
     * @param dx_n The Gauss-Newton point
     * @return The dogleg point \f$ \delta x_d \f$
     */
    static std::map<int,Eigen::VectorXd> ComputeDoglegPoint(double delta,
            const std::map<int,Eigen::VectorXd>& dx_u, const std::map<int,Eigen::VectorXd>& dx_n, const bool verbose=false);

    /** Compute the point on the line between the steepest descent point and the
     * Newton's method point intersecting the trust region boundary.
     * Mathematically, computes \f$ \tau \f$ such that \f$ 0<\tau<1 \f$ and
     * \f$ \| (1-\tau)\delta x_u + \tau\delta x_n \| = \delta \f$, where \f$ \delta \f$
     * is the trust region radius.
     * @param delta Trust region radius
     * @param x_u Steepest descent minimizer
     * @param x_n Newton's method minimizer
     */
    static std::map<int,Eigen::VectorXd> ComputeBlend(double delta,
            const std::map<int,Eigen::VectorXd>& x_u, const std::map<int,Eigen::VectorXd>& x_n, const bool verbose=false);
};




};
#endif // DOGLEGOPITIMIZERIMPL_H
