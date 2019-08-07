
/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 * @file    GaussianConditional.h
 * @brief   Conditional Gaussian Base class
 * @author  Christian Potthast
 */

#ifndef GAUSSIANCONDITIONAL_H
#define GAUSSIANCONDITIONAL_H


#include "../linear/JacobianFactor.h"
#include "../linear/NoiseModel.h"
#include "../base/SVBlockMatrix.h"

namespace minisam
{

/**
  * A conditional Gaussian functions as the node in a Bayes network
  * It has a set of parents y,z, etc. and implements a probability density on x.
  * The negative log-probability is given by \f$ \frac{1}{2} |Rx - (d - Sy - Tz - ...)|^2 \f$
  */
class GaussianConditional : public JacobianFactor
{
public:
    Factor* nrFrontals_;
    Factor* nrParents_;
public:
    /** default constructor needed for serialization */
    GaussianConditional(); // {}
    virtual ~GaussianConditional();
    GaussianConditional(const GaussianConditional &rObj);

    /** constructor with no parents |Rx-d| */
    GaussianConditional(int key, const Eigen::VectorXd &d, const Eigen::MatrixXd &R,
                        DiagonalNoiseModel *sigmas =new DiagonalNoiseModel());

    /** Constructor with arbitrary number of frontals and parents.
    *   @tparam TERMS A container whose value type is std::pair<Key, Matrix>, specifying the
    *           collection of keys and matrices making up the conditional. */
    GaussianConditional(const std::pair<int, Eigen::MatrixXd> &terms,
                        const JacobianFactor &nrFrontals, const JacobianFactor &nrParents, const Eigen::VectorXd &d,
                        DiagonalNoiseModel *sigmas =new DiagonalNoiseModel());

    /** Constructor with arbitrary number keys, and where the augmented matrix is given all together
     *  instead of in block terms.  Note that only the active view of the provided augmented matrix
     *  is used, and that the matrix data is copied into a newly-allocated matrix in the constructed
     *  factor. */
    GaussianConditional(
        const std::vector<int> &keys, int nrFrontalssize,const SVBlockMatrix& augmentedMatrix,
        DiagonalNoiseModel *sigmas =new DiagonalNoiseModel());

    /** Return a view of the upper-triangular R block of the conditional */
    Eigen::Block<const Eigen::MatrixXd> get_R() const;

    /** Get a view of the parent blocks. */
    Eigen::Block<const Eigen::MatrixXd> get_S() const;

    /** Get a view of the S matrix for the variable pointed to by the given key iterator */
    Eigen::Block<const Eigen::MatrixXd> get_S(std::vector<int>::const_iterator variable) const;
    /** Get a view of the r.h.s. vector d */
    const Eigen::VectorXd get_d() const;

    /**
    * Solves a conditional Gaussian and writes the solution into the entries of
    * \c x for each frontal variable of the conditional.  The parents are
    * assumed to have already been solved in and their values are read from \c x.
    * This function works for multiple frontal variables.
    *
    * Given the Gaussian conditional with log likelihood \f$ |R x_f - (d - S x_s)|^2 \f$,
    * where \f$ f \f$ are the frontal variables and \f$ s \f$ are the separator
    * variables of this conditional, this solve function computes
    * \f$ x_f = R^{-1} (d - S x_s) \f$ using back-substitution.
    *
    * @param parents VectorValues containing solved parents \f$ x_s \f$.
    */
    std::map<int, Eigen::VectorXd> solve(const std::map<int, Eigen::VectorXd> &parents) const;


    // operators

    GaussianConditional &operator=(const GaussianConditional &rObj);
    bool operator==(const GaussianConditional& other) const;

    int getsizenrFrontals() const;
    int getsizenrParents() const;

    std::vector<int>::const_iterator cbeginFrontals()
    const;
    std::vector<int>::const_iterator cendFrontals()
    const;

    std::vector<int>::const_iterator cbeginParents()
    const;
    std::vector<int>::const_iterator cendParents()
    const;

    /// @}
    /// @name Advanced Interface
    /// @{

    Factor* nrFrontals() const;
    Factor* nrParents() const;



}; // GaussianConditional

// operators inline

inline GaussianConditional &GaussianConditional::operator=(const GaussianConditional &rObj)
{
    JacobianFactor::operator=(rObj);
    nrFrontals_ =new Factor(*rObj.nrFrontals_);
    nrParents_ =new Factor(*rObj.nrParents_);
    return *this;
}

/**
   * Multiply all factors and eliminate the given keys from the resulting factor using a QR
   * variant that handles constraints (zero sigmas). Computation happens in noiseModel::Gaussian::QR
   * Returns a conditional on those keys, and a new factor on the separator.
   */
std::pair<GaussianConditional*, JacobianFactor*> EliminateQR(const std::vector<const RealGaussianFactor*> &factors,
        const std::vector<int>& keys);
/**
*   Densely partially eliminate with Cholesky factorization.  JacobianFactors are
*   left-multiplied with their transpose to form the Hessian using the conversion constructor
*   HessianFactor(const JacobianFactor&).
*
*   If any factors contain constrained noise models, this function will fail because our current
*   implementation cannot handle constrained noise models in Cholesky factorization.  The
*   function EliminatePreferCholesky() automatically does QR instead when this is the case.
*
*   Variables are eliminated in the order specified in \c keys.
*
*   @param factors Factors to combine and eliminate
*   @param keys The variables to eliminate and their elimination ordering
*   @return The conditional and remaining factor
*
*   \addtogroup LinearSolving */

std::pair<GaussianConditional*, HessianFactor*> EliminateCholesky(const std::vector<const RealGaussianFactor*>&factors,
        const std::vector<int>& keys);
/**
     *  In-place elimination that returns a conditional on (ordered) keys specified, and leaves
     *  this factor to be on the remaining keys (separator) only. Does dense partial Cholesky.
     */
GaussianConditional* HFeliminateCholesky(const std::vector<int>& keys, HessianFactor *HessianFactor);
GaussianConditional* splitConditional(JacobianFactor* jb, int nrFrontals);
/**
*   Densely partially eliminate with Cholesky factorization.  JacobianFactors are
*   left-multiplied with their transpose to form the Hessian using the conversion constructor
*   HessianFactor(const JacobianFactor&).
*
*   This function will fall back on QR factorization for any cliques containing JacobianFactor's
*   with constrained noise models.
*
*   Variables are eliminated in the order specified in \c keys.
*
*   @param factors Factors to combine and eliminate
*   @param keys The variables to eliminate and their elimination ordering
*   @return The conditional and remaining factor
*
*   \addtogroup LinearSolving */

std::pair<GaussianConditional*, RealGaussianFactor*>
EliminatePreferCholesky(const std::vector<const RealGaussianFactor*> &factors, const std::vector<int> &keys);
};

#endif // GAUSSIANCONDITIONAL_H
