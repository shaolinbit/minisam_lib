#ifndef MARGINALS_H
#define MARGINALS_H


/**
 * @file Marginals.h
 * @brief A class for computing marginals in a NonlinearFactorGraph
 * @author
 * @date
 */

//#pragma once


#include "../inference/BayesTreePointer.h"
#include "../nonlinear/NonlinearFactorGraph.h"

class JointMarginal;

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
    BayesTreePointer* bayesTree_;

public:

    /// Default constructor only for Cython wrapper
    Marginals() {}

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
    /** print */
    //void print(const std::string& str = "Marginals: ", const KeyFormatter& keyFormatter = DefaultKeyFormatter) const;

    /** Compute the marginal factor of a single variable */
    RealGaussianFactor* marginalFactor(int variable);

    /** Compute the marginal information matrix of a single variable.
      * Use LLt(const Matrix&) or RtR(const Matrix&) to obtain the square-root information matrix. */
    Eigen::MatrixXd marginalInformation(int variable);

    /** Compute the marginal covariance of a single variable */
    Eigen::MatrixXd marginalCovariance(int variable);

    /** Compute the joint marginal covariance of several variables */
    JointMarginal jointMarginalCovariance(const std::vector<int>& variables);

    /** Compute the joint marginal information of several variables */
    JointMarginal jointMarginalInformation(const std::vector<int>& variables);

    /** Optimize the bayes tree */
    std::map<int,Eigen::VectorXd> optimize() const;
};

/**
 * A class to store and access a joint marginal, returned from Marginals::jointMarginalCovariance and Marginals::jointMarginalInformation
 */
class  JointMarginal
{

protected:
    SVBlockMatrix blockMatrix_;
    std::vector<int> keys_;
    std::map<int, int> indices_;

public:
    /// Default constructor only for Cython wrapper
    JointMarginal() {}

    /** Access a block, corresponding to a pair of variables, of the joint
     * marginal.  Each block is accessed by its "vertical position",
     * corresponding to the variable with nonlinear Key \c iVariable and
     * "horizontal position", corresponding to the variable with nonlinear Key
     * \c jVariable.
     *
     * For example, if we have the joint marginal on a 2D pose "x3" and a 2D
     * landmark "l2", then jointMarginal(Symbol('x',3), Symbol('l',2)) will
     * return the 3x2 block of the joint covariance matrix corresponding to x3
     * and l2.
     * @param iVariable The nonlinear Key specifying the "vertical position" of the requested block
     * @param jVariable The nonlinear Key specifying the "horizontal position" of the requested block
     */
    Eigen::MatrixXd operator()(int iVariable, int jVariable) const
    {
        const auto indexI = indices_.at(iVariable);
        const auto indexJ = indices_.at(jVariable);
        return blockMatrix_.Sblock(indexI, indexJ);
    }

    /** Synonym for operator() */
    Eigen::MatrixXd at(int iVariable, int jVariable) const
    {
        return (*this)(iVariable, jVariable);
    }

    /** The full, dense covariance/information matrix of the joint marginal. */
    Eigen::MatrixXd fullMatrix() const
    {
        return blockMatrix_.selfadjointView();
    }

    /** Print */
    //void print(const std::string& s = "", const KeyFormatter& formatter = DefaultKeyFormatter) const;

protected:
    JointMarginal(const Eigen::MatrixXd& fullMatrix, const std::vector<int>& dims, const std::vector<int>& keys) :
        blockMatrix_(dims, fullMatrix),keys_(keys),indices_(Ordering(keys).invert())
    {
        //;
        //  Ordering nor(keys);
        //  indices_(nor.invert());
    }

    friend class Marginals;

};
#endif // MARGINALS_H
