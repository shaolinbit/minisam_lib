#ifndef JACOBIANFACTOR_H
#define JACOBIANFACTOR_H


/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 * @file    JacobianFactor.h
 * @author  Richard Roberts
 * @author  Christian Potthast
 * @author  Frank Dellaert
 * @date    Dec 8, 2010
 */


#include "../linear/RealGaussianFactor.h"
#include "../linear/GaussianFactorGraph.h"
#include "../linear/NoiseModel.h"
#include "../inference/Ordering.h"
#include "../base/SVBlockMatrix.h"
#include "../inference/VariableSlots.h"

namespace minisam
{

// Forward declarations
class GaussianFactorGraph;
class Ordering;

/**
   * A Gaussian factor in the squared-error form.
   *
   * JacobianFactor implements a
   * Gaussian, which has quadratic negative log-likelihood
   * \f[ E(x) = \frac{1}{2} (Ax-b)^T \Sigma^{-1} (Ax-b) \f]
   * where \f$ \Sigma \f$ is a \a diagonal covariance matrix.  The
   * matrix \f$ A \f$, r.h.s. vector \f$ b \f$, and diagonal noise model
   * \f$ \Sigma \f$ are stored in this class.
   *
   * This factor represents the sum-of-squares error of a \a linear
   * measurement function, and is created upon linearization of a NoiseModelFactor,
   * which in turn is a sum-of-squares factor with a nonlinear measurement function.
   *
   * Here is an example of how this factor represents a sum-of-squares error:
   *
   * Letting \f$ h(x) \f$ be a \a linear measurement prediction function, \f$ z \f$ be
   * the actual observed measurement, the residual is
   * \f[ f(x) = h(x) - z . \f]
   * If we expect noise with diagonal covariance matrix \f$ \Sigma \f$ on this
   * measurement, then the negative log-likelihood of the Gaussian induced by this
   * measurement model is
   * \f[ E(x) = \frac{1}{2} (h(x) - z)^T \Sigma^{-1} (h(x) - z) . \f]
   * Because \f$ h(x) \f$ is linear, we can write it as
   * \f[ h(x) = Ax + e \f]
   * and thus we have
   * \f[ E(x) = \frac{1}{2} (Ax-b)^T \Sigma^{-1} (Ax-b) \f]
   * where \f$ b = z - e \f$.
   *
   * This factor can involve an arbitrary number of variables, and in the
   * above example \f$ x \f$ would almost always be only be a subset of the variables
   * in the entire factor graph.  There are special constructors for 1-, 2-, and 3-
   * way JacobianFactors, and additional constructors for creating n-way JacobianFactors.
   * The Jacobian matrix \f$ A \f$ is passed to these constructors in blocks,
   * for example, for a 2-way factor, the constructor would accept \f$ A1 \f$ and \f$ A2 \f$,
   * as well as the variable indices \f$ j1 \f$ and \f$ j2 \f$
   * and the negative log-likelihood represented by this factor would be
   * \f[ E(x) = \frac{1}{2} (A_1 x_{j1} + A_2 x_{j2} - b)^T \Sigma^{-1} (A_1 x_{j1} + A_2 x_{j2} - b) . \f]
   */
class JacobianFactor : public RealGaussianFactor
{
public:
    /** default constructor for I/O */
    JacobianFactor();

    /** Convert from other GaussianFactor */ //
    explicit JacobianFactor(const RealGaussianFactor &gf);

    /** Copy constructor */
    explicit JacobianFactor(const JacobianFactor &jf) : RealGaussianFactor(jf.keys_,jf.Ab_, jf.model_, 0)
    {
    }

    /** Construct an n-ary factor
     * @tparam TERMS A container whose value type is std::pair<Key, Matrix>, specifying the
     *         collection of keys and matrices making up the factor. */
    JacobianFactor(const std::map<int, Eigen::MatrixXd> &terms, const Eigen::VectorXd &b,
                   DiagonalNoiseModel *model =new DiagonalNoiseModel());

    JacobianFactor(const std::vector<std::pair<int, Eigen::MatrixXd>> &terms,
                   const Eigen::VectorXd &b, DiagonalNoiseModel *model =new DiagonalNoiseModel());

    /** Constructor with arbitrary number keys, and where the augmented matrix is given all together
     *  instead of in block terms.  Note that only the active view of the provided augmented matrix
     *  is used, and that the matrix data is copied into a newly-allocated matrix in the constructed
     *  factor. */
    JacobianFactor(
        const std::vector<int> &keys, const SVBlockMatrix& augmentedMatrix,  DiagonalNoiseModel *sigmas =new DiagonalNoiseModel());

    /**
     * Build a dense joint factor from all the factors in a factor graph.  If a VariableSlots
     * structure computed for \c graph is already available, providing it will reduce the amount of
     * computation performed.

    explicit JacobianFactor(
        const GaussianFactorGraph &graph);*/
    explicit JacobianFactor(const GaussianFactorGraph& graph);

    explicit JacobianFactor(const GaussianFactorGraph& graph,
                            const Ordering &ordering);
    explicit JacobianFactor(const std::vector<const RealGaussianFactor*>& graph,
                            const Ordering &ordering);
    /** destructor */
    virtual ~JacobianFactor();

    /*Clone this JacobianFactor*/
    virtual RealGaussianFactor *clone() const;

    //V
    /** Return the full augmented Jacobian matrix of this factor as a VerticalBlockMatrix object. */
    const SVBlockMatrix& matrixObject() const;

    //V
    /** Mutable access to the full augmented Jacobian matrix of this factor as a VerticalBlockMatrix object. */
    SVBlockMatrix& matrixObject();

    /** is noise model constrained ? */
    bool isConstrained() const;

    /**
     * return the number of rows in the corresponding linear system
     */
    int rows() const;

    /**
     * return the number of columns in the corresponding linear system
     */
    int cols() const;

    JacobianFactor(JacobianFactor &jf);

    /** get a copy of model */
    const DiagonalNoiseModel* get_model() const;

    /** get a copy of model (non-const version) */
    DiagonalNoiseModel* get_model();
    /** Get a view of the r.h.s. vector b, not weighted by noise */

    Eigen::VectorXd getb();

    const Eigen::VectorXd getb() const;
    /** Get a view of the A matrix for the variable pointed to by the given key iterator */
    Eigen::Block<const Eigen::MatrixXd> getA(std::vector<int>::const_iterator variable) const;

    /** Get a view of the A matrix, not weighted by noise */
    Eigen::Block<const Eigen::MatrixXd> getA() const;

    /** Get a view of the A matrix for the variable pointed to by the given key iterator (non-const version) */
    Eigen::Block<Eigen::MatrixXd> getA(std::vector<int>::iterator variable);

    /** Get a view of the A matrix for the variable pointed to by the given key iterator (non-const version) */
    Eigen::Block<const Eigen::MatrixXd> getA(int vpos) const;

    /** Get a view of the A matrix */
    Eigen::Block<Eigen::MatrixXd> getA();

    /** Return A*x */
    Eigen::VectorXd operator*(const std::map<int, Eigen::VectorXd> &x) const;


    /** Return a whitened version of the factor, i.e. with unit diagonal noise model. */
    JacobianFactor& whiten() const;

    /** set noiseModel correctly */
    void setModel(bool anyConstrained, const Eigen::VectorXd &sigmas);

    JacobianFactor &operator=(const JacobianFactor &rObj);

protected:
    /// Internal function to fill blocks and set dimensions
    void fillTerms(const std::map<int, Eigen::MatrixXd> &terms, const Eigen::VectorXd &b,  DiagonalNoiseModel* noiseModel);
    void fillVecTerms(const std::vector<std::pair<int, Eigen::MatrixXd>> &terms, const Eigen::VectorXd &b,
                      DiagonalNoiseModel *noiseModel);

public:
    JacobianFactor(const std::vector<int>& keys,const std::vector<int>& dims, int m,
                   DiagonalNoiseModel* model =new DiagonalNoiseModel()) : RealGaussianFactor(keys,SVBlockMatrix(dims.begin(), dims.end(), m, true, false), model, 0)
    {
    }

    friend class NonlinearFactor;
}; // JacobianFactor
std::vector<int> JFgetvectorint_fromiterator(const Factor &nrParents_);

std::vector<JacobianFactor*> _convertOrCastToJacobians(
    const GaussianFactorGraph &factors);
std::vector<const JacobianFactor*> _convertOrCastToJacobiansConst(
    const std::vector<const RealGaussianFactor*> &factors);

};
#endif // JACOBIANFACTOR_H
