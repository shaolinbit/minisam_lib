#ifndef HESSIANFACTOR_H
#define HESSIANFACTOR_H
/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 * @file    HessianFactor.h
 * @brief   Contains the HessianFactor class, a general quadratic factor
 * @author  Richard Roberts
 * @date    Dec 8, 2010
 */
#include "../linear/RealGaussianFactor.h"
#include "../linear/Scatter.h"
#include "../base/SVBlockMatrix.h"
namespace minisam
{

// Forward declarations
class Ordering;
class GaussianFactorGraph;
class Scatter;

/**
 * @brief A Gaussian factor using the canonical parameters (information form)
 *
 * HessianFactor implements a general quadratic factor of the form
 * \f[ E(x) = 0.5 x^T G x - x^T g + 0.5 f \f]
 * that stores the matrix \f$ G \f$, the vector \f$ g \f$, and the constant term \f$ f \f$.
 *
 * When \f$ G \f$ is positive semidefinite, this factor represents a Gaussian,
 * in which case \f$ G \f$ is the information matrix \f$ \Lambda \f$,
 * \f$ g \f$ is the information vector \f$ \eta \f$, and \f$ f \f$ is the residual
 * sum-square-error at the mean, when \f$ x = \mu \f$.
 *
 * Indeed, the negative log-likelihood of a Gaussian is (up to a constant)
 * @f$ E(x) = 0.5(x-\mu)^T P^{-1} (x-\mu) @f$
 * with @f$ \mu @f$ the mean and  @f$ P @f$ the covariance matrix. Expanding the product we get
 * @f[
 * E(x) = 0.5 x^T P^{-1} x - x^T P^{-1} \mu + 0.5 \mu^T P^{-1} \mu
 * @f]
 * We define the Information matrix (or Hessian) @f$ \Lambda = P^{-1} @f$
 * and the information vector @f$ \eta = P^{-1} \mu = \Lambda \mu @f$
 * to arrive at the canonical form of the Gaussian:
 * @f[
 * E(x) = 0.5 x^T \Lambda x - x^T \eta + 0.5 \mu^T \Lambda \mu
 * @f]
 *
 * This factor is one of the factors that can be in a GaussianFactorGraph.
 * It may be returned from NonlinearFactor::linearize(), but is also
 * used internally to store the Hessian during Cholesky elimination.
 *
 * This can represent a quadratic factor with characteristics that cannot be
 * represented using a JacobianFactor (which has the form
 * \f$ E(x) = \Vert Ax - b \Vert^2 \f$ and stores the Jacobian \f$ A \f$
 * and error vector \f$ b \f$, i.e. is a sum-of-squares factor).  For example,
 * a HessianFactor need not be positive semidefinite, it can be indefinite or
 * even negative semidefinite.
 *
 * If a HessianFactor is indefinite or negative semi-definite, then in order
 * for solving the linear system to be possible,
 * the Hessian of the full system must be positive definite (i.e. when all
 * small Hessians are combined, the result must be positive definite).  If
 * this is not the case, an error will occur during elimination.
 *
 * This class stores G, g, and f as an augmented matrix HessianFactor::matrix_.
 * The upper-left n x n blocks of HessianFactor::matrix_ store the upper-right
 * triangle of G, the upper-right-most column of length n of HessianFactor::matrix_
 * stores g, and the lower-right entry of HessianFactor::matrix_ stores f, i.e.
 * \code
   HessianFactor::matrix_ = [ G11 G12 G13 ... g1
                                0 G22 G23 ... g2
                                0   0 G33 ... g3
                                :   :   :      :
                                0   0   0 ...  f ]
   \endcode
   Blocks can be accessed as follows:
   \code
   G11 = info(begin(), begin());
   G12 = info(begin(), begin()+1);
   G23 = info(begin()+1, begin()+2);
   g2 = linearTerm(begin()+1);
   f = constantTerm();
   .......
   \endcode
 */
class  HessianFactor : public RealGaussianFactor
{
public:

    /** default constructor for I/O */
    HessianFactor();

    HessianFactor(const RealGaussianFactor& factor);
    /** Construct a unary factor.  G is the quadratic term (Hessian matrix), g
     * the linear term (a vector), and f the constant term.  The quadratic
     * error is:
     * 0.5*(f - 2*x'*g + x'*G*x)
     */
    HessianFactor(int j, const Eigen::MatrixXd& G, const Eigen::VectorXd& g, double f);

    /** Construct a unary factor, given a mean and covariance matrix.
     * error is 0.5*(x-mu)'*inv(Sigma)*(x-mu)
    */
    HessianFactor(int j, const Eigen::VectorXd& mu, const Eigen::MatrixXd& Sigma);

    HessianFactor(const std::vector<std::pair<int, Eigen::MatrixXd>> &terms,
                  const Eigen::VectorXd &b);

    /** Combine a set of factors into a single dense HessianFactor */
    explicit HessianFactor(const GaussianFactorGraph& factors,Scatter& scatter);
    explicit HessianFactor(const std::vector<const RealGaussianFactor*>& factors,Scatter& scatter);

    /** Destructor */
    virtual ~HessianFactor() {}

    /** Clone this HessianFactor */
    virtual RealGaussianFactor* clone() const;

    /** Return the number of columns and rows of the Hessian matrix, including the information vector. */
    int rows() const;

    /** Return the constant term \f$ f \f$ as described above
     * @return The constant term \f$ f \f$
     */
    double constantTerm() const;
    /** Return the constant term \f$ f \f$ as described above
     * @return The constant term \f$ f \f$
     */
    double& constantTerm();

    /** Return the part of linear term \f$ g \f$ as described above corresponding to the requested variable.
     * @param j Which block row to get, as an iterator pointing to the slot in this factor.  You can
     * use, for example, begin() + 2 to get the 3rd variable in this factor.
     * @return The linear term \f$ g \f$ */
    Eigen::Block<const Eigen::MatrixXd>  linearTerm(std::vector<int>::const_iterator j) const;


    /** Return the complete linear term \f$ g \f$ as described above.
     * @return The linear term \f$ g \f$ */
    Eigen::Block<const Eigen::MatrixXd>  linearTerm() const;
    /** Return the complete linear term \f$ g \f$ as described above.
     * @return The linear term \f$ g \f$ */
    Eigen::Block<Eigen::MatrixXd> linearTerm();

    //S
    /// Return underlying information matrix.
    const SVBlockMatrix& info() const;

    //S
    /// Return non-const information matrix.
    SVBlockMatrix& info();
    /// Return self-adjoint view onto the information matrix (NOT augmented).
    Eigen::SelfAdjointView<Eigen::Block<const Eigen::MatrixXd>, Eigen::Upper> informationView() const;

    /// Solve the system A'*A delta = A'*b in-place, return delta as std::map<int,Eigen::VectorXd>
    std::map<int,Eigen::VectorXd> solve();

public:

    /// Allocate for given scatter pattern
    void Allocate(const Scatter& scatter);

    /// Constructor with given scatter pattern, allocating but not initializing storage.
    HessianFactor(const Scatter& scatter);

};

void _FromJacobianHelper(const RealGaussianFactor& jf,SVBlockMatrix& info);

};

#endif // HESSIANFACTOR_H
