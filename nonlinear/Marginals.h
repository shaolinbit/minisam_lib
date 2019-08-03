#ifndef MARGINALS_H
#define MARGINALS_H


/**
 * @file Marginals.h
 * @brief A class for computing marginals in a NonlinearFactorGraph
 * @author Richard Roberts
 * @date May 14, 2012
 */


#include "../inference/BayesTree.h"
#include "../nonlinear/NonlinearFactorGraph.h"


namespace minisam
{


/**
 * A class for computing Gaussian marginals of variables in a NonlinearFactorGraph
 */
class  Marginals
{

public:

    /** The linear factorization mode - either CHOLESKY (faster and suitable for most problems) or QR (slower but more numerically stable for poorly-conditioned problems). */
    enum Factorization
    {
        CHOLESKY,
        QR
    };

protected:

    GaussianFactorGraph graph_;
    std::map<int,Eigen::VectorXd> values_;
#ifdef GMF_Using_Pose3
    std::map<int,Pose3> poses_;
#else
    std::map<int,Pose2> poses_;
#endif // GMF_Using_Pose3
    Factorization factorization_;
public:
    BayesTree* bayesTree_;



    /// Default constructor
    Marginals();
    ~Marginals();

    void clearlineargraph();



    /** Construct a marginals class.
     * @param graph The factor graph defining the full joint density on all variables.
     * @param solution The linearization point about which to compute Gaussian marginals (usually the MLE as obtained from a NonlinearOptimizer).
     * @param factorization The linear decomposition mode - either Marginals::CHOLESKY (faster and suitable for most problems) or Marginals::QR (slower but more numerically stable for poorly-conditioned problems).
     * @param ordering An optional variable ordering for elimination.
     */
#ifdef GMF_Using_Pose3
    Marginals(const NonlinearFactorGraph& graph, const std::map<int,Eigen::VectorXd>& solution,
              const std::map<int,Pose3>& sposes_,
              Ordering ordering, Factorization factorization = CHOLESKY);
    Marginals(const NonlinearFactorGraph& graph, const std::map<int,Eigen::VectorXd>& solution,
              const std::map<int,Pose3>& sposes_,
              Factorization factorization = CHOLESKY);
#else
    Marginals(const NonlinearFactorGraph& graph, const std::map<int,Eigen::VectorXd>& solution,
              const std::map<int,Pose2>& sposes_,
              Ordering ordering, Factorization factorization = CHOLESKY);
    Marginals(const NonlinearFactorGraph& graph, const std::map<int,Eigen::VectorXd>& solution,
              const std::map<int,Pose2>& sposes_,
              Factorization factorization = CHOLESKY);
#endif // GMF_Using_Pose3

    /** Compute the marginal factor of a single variable */
    std::pair<RealGaussianFactor*,GaussianBayesNet*>  marginalFactor(int variable);

    /** Compute the marginal information matrix of a single variable.
      * Use LLt(const Matrix&) or RtR(const Matrix&) to obtain the square-root information matrix. */
    Eigen::MatrixXd marginalInformation(int variable);

    /** Compute the marginal covariance of a single variable */
    Eigen::MatrixXd marginalCovariance(int variable);

    /** Optimize the bayes tree */
    std::map<int,Eigen::VectorXd> optimize() const;
};

};
#endif // MARGINALS_H
