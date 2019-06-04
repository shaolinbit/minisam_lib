#ifndef REALGAUSSIANFACTOR_H_INCLUDED
#define REALGAUSSIANFACTOR_H_INCLUDED

/**
 * @file    RealGaussianFactor.h
 * @brief   A factor with a quadratic error function - a Gaussian
 * @brief   GaussianFactor
 * @author
 */

//#pragma once

#include "../inference/Factor.h"
#include "../base/Matrix.h"
#include "../base/SVBlockMatrix.h"
#include "../linear/NoiseModel.h"
#include <map>

//class SVBlockMatrix;

/**
   * An abstract virtual base class for JacobianFactor and HessianFactor. A GaussianFactor has a
   * quadratic error function. GaussianFactor is non-mutable (all methods const!). The factor value
   * is exp(-0.5*||Ax-b||^2) */
class RealGaussianFactor// : public Factor
{

public:
    SVBlockMatrix Ab_;         // the block view of the full matrix
    DiagonalNoiseModel* model_; // Gaussian noise model with diagonal covariance matrix
    int TypeGaussianFactor;    //{Jacobian,Hessian};
   // int cliqueindex_;
    bool iswrapper_;
    std::vector<int> keys_;

public:
    /** Default constructor creates empty factor */
    RealGaussianFactor() {model_=NULL;iswrapper_=false;}

    /** Construct from container of keys.  This constructor is used internally from derived factor
     *  constructors, either from a container of keys or from a boost::assign::list_of. */
    //template<typename CONTAINER>
    RealGaussianFactor(const std::vector<int> &keys, const SVBlockMatrix& sb,DiagonalNoiseModel *model, int GaussianFactorT) :
         keys_(keys),
        Ab_(sb), model_(model), TypeGaussianFactor(GaussianFactorT),
      iswrapper_(false) {}

    /** Destructor */
    virtual ~RealGaussianFactor();

    RealGaussianFactor(const RealGaussianFactor &rf) : keys_(rf.keys()),
        Ab_(rf.Ab_), model_(rf.model_), TypeGaussianFactor(rf.TypeGaussianFactor),
        iswrapper_(rf.iswrapper_) {}

    /** Equals for testable */
    const SVBlockMatrix& info() const
    {
        return Ab_;
    }
    int rows() const
    {
        return Ab_.rows();
    }

    Eigen::VectorXd JFunweighted_error(const std::map<int, Eigen::VectorXd> &c) const;

    Eigen::VectorXd JFerror_vector(const std::map<int, Eigen::VectorXd> &c) const;
    /** Print for testable */
    double error(const std::map<int, Eigen::VectorXd> &c) const; /**  0.5*(A*x-b)'*D*(A*x-b) */

    /** Return the dimension of the variable pointed to by the given key iterator */
    int getDim(std::vector<int>::const_iterator variable) const;

    /**
     * Return a dense \f$ [ \;A\;b\; ] \in \mathbb{R}^{m \times n+1} \f$
     * Jacobian matrix, augmented with b with the noise models baked
     * into A and b.  The negative log-likelihood is
     * \f$ \frac{1}{2} \Vert Ax-b \Vert^2 \f$.  See also
     * GaussianFactorGraph::jacobian and GaussianFactorGraph::sparseJacobian.
     */
    Eigen::MatrixXd augmentedJacobian();

      /// Access the factor's involved variable keys
    const std::vector<int>& keys() const
    {
        return keys_;
    }

    /** @return keys involved in this factor */
    std::vector<int>& keys()
    {
        return keys_;
    }

    /**
     * Return the dense Jacobian \f$ A \f$ and right-hand-side \f$ b \f$,
     * with the noise models baked into A and b. The negative log-likelihood
     * is \f$ \frac{1}{2} \Vert Ax-b \Vert^2 \f$.  See also
     * GaussianFactorGraph::augmentedJacobian and
     * GaussianFactorGraph::sparseJacobian.
     */
    std::pair<Eigen::MatrixXd, Eigen::VectorXd> jacobian() const;

    /** Return the augmented information matrix represented by this GaussianFactor.
     * The augmented information matrix contains the information matrix with an
     * additional column holding the information vector, and an additional row
     * holding the transpose of the information vector.  The lower-right entry
     * contains the constant error term (when \f$ \delta x = 0 \f$).  The
     * augmented information matrix is described in more detail in HessianFactor,
     * which in fact stores an augmented information matrix.
     */
    Eigen::MatrixXd augmentedInformation() const;

    /** Return the non-augmented information matrix represented by this
     * GaussianFactor.
     */
    Eigen::MatrixXd information() const;

    /// Return the diagonal of the Hessian for this factor
    std::map<int, Eigen::VectorXd> hessianDiagonal() const;

    /// Return the block diagonal of the Hessian for this factor
    std::map<int, Eigen::MatrixXd> hessianBlockDiagonal() const;

    /** Clone a factor (make a deep copy) */
    RealGaussianFactor *clone() const;

    RealGaussianFactor &operator=(const RealGaussianFactor &rObj)
    {
        //Factor::operator=(rObj);
        keys_=rObj.keys_;
        this->Ab_ = rObj.Ab_;
        // cout<<this->Ab_.matrix()<<endl;
        this->model_ = rObj.model_;
        this->TypeGaussianFactor = rObj.TypeGaussianFactor;
        //this->cliqueindex_ = rObj.cliqueindex_;
        this->iswrapper_=rObj.iswrapper_;
        return *this;
    }

    /**
     * Construct the corresponding anti-factor to negate information
     * stored stored in this factor.
     * @return a HessianFactor with negated Hessian matrices
     */
    RealGaussianFactor &negate() const;

    RealGaussianFactor *negate();

    /** Update an information matrix by adding the information corresponding to this factor
     * (used internally during elimination).
     * @param scatter A mapping from variable index to slot index in this HessianFactor
     * @param info The information matrix to be updated
     */

    ///S
    void updateHessian(const std::vector<int> &keys,
                        SVBlockMatrix* info) const;



    //from J
    void transposeMultiplyAdd(double alpha, const Eigen::VectorXd &e,
                              std::map<int, Eigen::VectorXd> &x) const;
    //from J
    Eigen::VectorXd MultiplyVectorValuesX(const std::map<int, Eigen::VectorXd> &x) const;

    /// y += alpha * A'*A*x
    void multiplyHessianAdd(double alpha, const std::map<int, Eigen::VectorXd> &x, std::map<int, Eigen::VectorXd> &y) const;

    /// A'*b for Jacobian, eta for Hessian
    std::map<int, Eigen::VectorXd> gradientAtZero() const;

    /// Raw memory access version of gradientAtZero
    // void gradientAtZero(double* d) const;

    /// Gradient wrt a key at any values
    Eigen::VectorXd gradient(int key, const std::map<int, Eigen::VectorXd> &x) const;

    /** Test whether the factor is empty */
    bool empty() const
    {
        return keys_.size() == 0;
    }

    int size() const
    {
        return keys_.size();
    }


    // Determine position of a given key
    // template <typename CONTAINER>
    static int Slot(const std::vector<int> &keys, int key)
    {
        return std::find(keys.begin(), keys.end(), key) - keys.begin();
    }

}; // RealGaussianFactor
void setRealGaussianFactor(RealGaussianFactor* r1,const RealGaussianFactor* r2);
#endif // REALGAUSSIANFACTOR_H_INCLUDED
