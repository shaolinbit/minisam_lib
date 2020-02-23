#ifndef REALGAUSSIANFACTOR_H_INCLUDED
#define REALGAUSSIANFACTOR_H_INCLUDED

///Modified from GaussianFactor.h. However, this class is not a virtual class any more.

/**
 * @file    GaussianFactor.h
 * @brief   A factor with a quadratic error function - GaussianFactor
 */





#include "../inference/Factor.h"
#include "../mat/Matrix.h"
#include "../mat/MatCal.h"
#include "../mat/GaussianBlockMatrix.h"
#include "../linear/NoiseModel.h"
#include <map>

namespace minisam
{

/**
   * A base class for JacobianFactor and HessianFactor. A GaussianFactor has a
   * quadratic error function. The factor value
   * is exp(-0.5*||Ax-b||^2) */
    class RealGaussianFactor
{

public:
    GaussianBlockMatrix Ab_;         // the block view of the full matrix
    GaussianNoiseModel* model_; // Gaussian noise model with diagonal covariance matrix
    int TypeGaussianFactor;    //{Jacobian,Hessian};
    bool iswrapper_;
    std::vector<int> keys_;

public:
    /** Default constructor creates empty factor */
    RealGaussianFactor();

    /** Construct from container of keys.  This constructor is used internally from derived factor
     *  constructors, either from a container of keys or from a boost::assign::list_of. */
    RealGaussianFactor(const std::vector<int> &keys,
                       const GaussianBlockMatrix& sb,GaussianNoiseModel *model, int GaussianFactorT);
    /** Destructor */
    virtual ~RealGaussianFactor();

    RealGaussianFactor(const RealGaussianFactor &rf);

    /** Equals for testable */
    const GaussianBlockMatrix& info() const;
    int rows() const;

    minivector JFunweighted_error(const std::map<int, minivector> &c) const;

    minivector JFerror_vector(const std::map<int, minivector> &c) const;
    /** Print for testable */
    double error(const std::map<int, minivector> &c) const; /**  0.5*(A*x-b)'*D*(A*x-b) */

    /** Return the dimension of the variable pointed to by the given key iterator */
    int getDim(std::vector<int>::const_iterator variable) const;

    /**
     * Return a dense \f$ [ \;A\;b\; ] \in \mathbb{R}^{m \times n+1} \f$
     * Jacobian matrix, augmented with b with the noise models baked
     * into A and b.  The negative log-likelihood is
     * \f$ \frac{1}{2} \Vert Ax-b \Vert^2 \f$.  See also
     * GaussianFactorGraph::jacobian and GaussianFactorGraph::sparseJacobian.
     */
    minimatrix augmentedJacobian();
    minimatrix information() const;

    /// Access the factor's involved variable keys
    const std::vector<int>& keys() const;

    /** @return keys involved in this factor */
    std::vector<int>& keys();

    /** Clone a factor (make a deep copy) */
    RealGaussianFactor *clone() const;

    RealGaussianFactor &operator=(const RealGaussianFactor &rObj);


    /** Update an information matrix by adding the information corresponding to this factor
     * (used internally during elimination).
     * @param scatter A mapping from variable index to slot index in this HessianFactor
     * @param info The information matrix to be updated
     */

    ///S
    void updateHessian(const std::vector<int> &keys,
                       GaussianBlockMatrix* info) const;

    std::map<int, minivector> hessianDiagonal() const;


    void transposeMultiplyAdd(double alpha, const minivector &e,
                              std::map<int, minivector> &x) const;

    std::map<int, minivector> gradientAtZero() const;


    /** Test whether the factor is empty */
    bool empty() const;

    int size() const;

    // Determine position of a given key
    static int Slot(const std::vector<int> &keys, int key);

}; // RealGaussianFactor
void setRealGaussianFactor(RealGaussianFactor* r1,const RealGaussianFactor* r2);
};

#endif // REALGAUSSIANFACTOR_H_INCLUDED
