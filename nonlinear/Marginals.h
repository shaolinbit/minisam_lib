#ifndef MARGINALS_H
#define MARGINALS_H


/**
 * @file Marginals.h
 * @brief A class for computing marginals in a NonlinearFactorGraph
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
    std::map<int,minimatrix*> values_;
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
    Marginals(const NonlinearFactorGraph& graph, const std::map<int,minimatrix*>& solution,
              const std::vector<int>& ordering,const Factorization& factorization= CHOLESKY);
    Marginals(const NonlinearFactorGraph& graph, const std::map<int,minimatrix*>& solution,
              const Factorization& factorization= CHOLESKY);
    /** Compute the marginal factor of a single variable */
    std::pair<RealGaussianFactor*,GaussianBayesNet*>  marginalFactor(int variable);

    /** Compute the marginal information matrix of a single variable.
      * Use LLt(const Matrix&) or RtR(const Matrix&) to obtain the square-root information matrix. */
    minimatrix marginalInformation(int variable);

    /** Compute the marginal covariance of a single variable */
    minimatrix marginalCovariance(int variable);

    /** Optimize the bayes tree */
    std::map<int,minivector> optimize() const;
};

};
#endif // MARGINALS_H
