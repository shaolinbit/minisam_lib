#ifndef GAUSSIANBAYESNETPOINTER_H_INCLUDED
#define GAUSSIANBAYESNETPOINTER_H_INCLUDED

/**
 * @file    GaussianBayesNet.h
 * @brief   Chordal Bayes Net, the result of eliminating a factor graph
 * @brief   GaussianBayesNet
 * @author
 */

// \callgraph

//#pragma once

#include "../linear/GaussianConditional.h"
#include "../inference/FactorPointerGraph.h"
#include "../inference/Ordering.h"
#include <Eigen/Core>

//class Ordering;
//class GaussianBayesNet;
class GaussianConditional;
/** A Bayes net made from linear-Gaussian densities */
class GaussianBayesNetPointer: public FactorPointerGraph<GaussianConditional>
{
public:
    /// @name Standard Constructors
    /// @{

    /** Construct empty factor graph */
    GaussianBayesNetPointer() {}

    /** Construct from iterator over conditionals */
    //template<typename ITERATOR>
    //GaussianBayesNetPointer(ITERATOR firstConditional, ITERATOR lastConditional) : FactorGraph<GaussianConditional>(firstConditional, lastConditional) {}

    /** Construct from container of factors (shared_ptr or plain objects) */
    //template<class CONTAINER>
    explicit GaussianBayesNetPointer(const std::vector<GaussianConditional*>& conditionals) : FactorPointerGraph<GaussianConditional>(conditionals) {}

    /** Implicit copy/downcast constructor to override explicit template container constructor */
    GaussianBayesNetPointer(FactorPointerGraph<GaussianConditional>& graph) : FactorPointerGraph<GaussianConditional>(graph) {}
    /// @}


    /// @name Standard Interface
    /// @{

    /// Solve the GaussianBayesNet, i.e. return \f$ x = R^{-1}*d \f$, by back-substitution
    std::map<int,Eigen::VectorXd> optimize() const;

    /// Version of optimize for incomplete BayesNet, needs solution for missing variables
    std::map<int,Eigen::VectorXd> optimize(const std::map<int,Eigen::VectorXd>& solutionForMissing) const;

    ///@}

    ///@name Linear Algebra
    ///@{

    /**
     * Return (dense) upper-triangular matrix representation
     */
    std::pair<Eigen::MatrixXd, Eigen::VectorXd> matrix() const;

    /**
     * Optimize along the gradient direction, with a closed-form computation to perform the line
     * search.  The gradient is computed about \f$ \delta x=0 \f$.
     *
     * This function returns \f$ \delta x \f$ that minimizes a reparametrized problem.  The error
     * function of a GaussianBayesNet is
     *
     * \f[ f(\delta x) = \frac{1}{2} |R \delta x - d|^2 = \frac{1}{2}d^T d - d^T R \delta x +
     * \frac{1}{2} \delta x^T R^T R \delta x \f]
     *
     * with gradient and Hessian
     *
     * \f[ g(\delta x) = R^T(R\delta x - d), \qquad G(\delta x) = R^T R. \f]
     *
     * This function performs the line search in the direction of the gradient evaluated at \f$ g =
     * g(\delta x = 0) \f$ with step size \f$ \alpha \f$ that minimizes \f$ f(\delta x = \alpha g)
     * \f$:
     *
     * \f[ f(\alpha) = \frac{1}{2} d^T d + g^T \delta x + \frac{1}{2} \alpha^2 g^T G g \f]
     *
     * Optimizing by setting the derivative to zero yields \f$ \hat \alpha = (-g^T g) / (g^T G g)
     * \f$.  For efficiency, this function evaluates the denominator without computing the Hessian
     * \f$ G \f$, returning
     *
     * \f[ \delta x = \hat\alpha g = \frac{-g^T g}{(R g)^T(R g)} \f] */
    std::map<int,Eigen::VectorXd> optimizeGradientSearch() const;

    /** Compute the gradient of the energy function, \f$ \nabla_{x=x_0} \left\Vert \Sigma^{-1} R x -
     * d \right\Vert^2 \f$, centered around \f$ x = x_0 \f$. The gradient is \f$ R^T(Rx-d) \f$.
     *
     * @param x0 The center about which to compute the gradient
     * @return The gradient as a VectorValues */
    std::map<int,Eigen::VectorXd> gradient(const std::map<int,Eigen::VectorXd>& x0) const;

    /** Compute the gradient of the energy function, \f$ \nabla_{x=0} \left\Vert \Sigma^{-1} R x - d
     * \right\Vert^2 \f$, centered around zero. The gradient about zero is \f$ -R^T d \f$.  See also
     * gradient(const GaussianBayesNet&, const VectorValues&).
     *
     * @param [output] g A VectorValues to store the gradient, which must be preallocated, see
     *        allocateVectorValues */
    std::map<int,Eigen::VectorXd> gradientAtZero() const;

    /** Mahalanobis norm error. */
    double error(const std::map<int,Eigen::VectorXd>& x) const;

    /**
     * Computes the determinant of a GassianBayesNet. A GaussianBayesNet is an upper triangular
     * matrix and for an upper triangular matrix determinant is the product of the diagonal
     * elements. Instead of actually multiplying we add the logarithms of the diagonal elements and
     * take the exponent at the end because this is more numerically stable.
     * @param bayesNet The input GaussianBayesNet
     * @return The determinant */
    double determinant() const;

    /**
     * Computes the log of the determinant of a GassianBayesNet. A GaussianBayesNet is an upper
     * triangular matrix and for an upper triangular matrix determinant is the product of the
     * diagonal elements.
     * @param bayesNet The input GaussianBayesNet
     * @return The determinant */
    double logDeterminant() const;

    /**
     * Backsubstitute with a different RHS vector than the one stored in this BayesNet.
     * gy=inv(R*inv(Sigma))*gx
     */
    std::map<int,Eigen::VectorXd> backSubstitute(const std::map<int,Eigen::VectorXd>& gx) const;

    /**
     * Transpose backsubstitute with a different RHS vector than the one stored in this BayesNet.
     * gy=inv(L)*gx by solving L*gy=gx.
     * gy=inv(R'*inv(Sigma))*gx
     * gz'*R'=gx', gy = gz.*sigmas
     */
    std::map<int,Eigen::VectorXd> backSubstituteTranspose(const std::map<int,Eigen::VectorXd>& gx) const;


};


GaussianFactorGraph* getgfrombnPointer(GaussianBayesNetPointer* BN);


#endif // GAUSSIANBAYESNETPOINTER_H_INCLUDED
