#ifndef NOISEMODEL_H_INCLUDED
#define NOISEMODEL_H_INCLUDED

#include <Eigen/Core>
#include <Eigen/Cholesky>
#include "../base/Matrix.h"

#include <iostream>
#include "../gmfconfig.h"

using namespace std;

/// All noise models live in the noiseModel namespace
// namespace noiseModel {

// Forward declaration
class GaussianNoiseModel;
class DiagonalNoiseModel;
class ConstrainedNoiseModel;
class IsotropicNoiseModel;
class UnitNoiseModel;

//---------------------------------------------------------------------------------------

/**
 * noiseModel::Base is the abstract base class for all noise models.
 *
 * Noise models must implement a 'whiten' function to normalize an error vector,
 * and an 'unwhiten' function to unnormalize an error vector.
 */
class SharedNoiseModel
{
    // public:
    //  typedef boost::shared_ptr<Base> shared_ptr;

protected:
    int dim_;
    bool isConstrained_; // default false
    bool isUnit_;        // default false

public:
    /// primary constructor @param dim is the dimension of the model
    SharedNoiseModel(int dim = 1) : dim_(dim), isConstrained_(false), isUnit_(false) {}
    virtual ~SharedNoiseModel();

    SharedNoiseModel &operator=(const SharedNoiseModel &rObj)
    {
        this->dim_ = rObj.dim_;
        this->isConstrained_ = rObj.isConstrained_;
        this->isUnit_ = rObj.isUnit_;
        return *this;
    }

    /// Dimensionality
    int dim() const
    {
        return dim_;
    }
    void setdim(int sdim)
    {
        dim_ = sdim;
    }
    /// true if a constrained noise model, saves slow/clumsy dynamic casting
    bool isConstrained() const
    {
        return isConstrained_;    // default false
    }
    /// true if a unit noise model, saves slow/clumsy dynamic casting
    bool isUnit() const
    {
        return isUnit_;    // default false
    }

    // virtual void print(const std::string& name = "") const = 0;
    // virtual bool equals(const SharedNoiseModel& expected, double tol=1e-9)
    // const = 0;

    /// Calculate standard deviations
    virtual Eigen::VectorXd sigmas() const = 0;

    /// Whiten an error vector.
    virtual Eigen::VectorXd whiten(const Eigen::VectorXd &v) const = 0;
    /// Whiten a matrix.
    virtual Eigen::MatrixXd Whiten(const Eigen::MatrixXd &H) const = 0;
    /// Unwhiten an error vector.
    virtual Eigen::VectorXd unwhiten(const Eigen::VectorXd &v) const = 0;

    virtual double distance(const Eigen::VectorXd &v) const = 0;

    virtual void WhitenSystem(std::vector<Eigen::MatrixXd> &A,
                              Eigen::VectorXd &b) const = 0;
    virtual void WhitenSystem(Eigen::MatrixXd &A, Eigen::VectorXd &b) const = 0;
    virtual void WhitenSystem(Eigen::MatrixXd &A1, Eigen::MatrixXd &A2,
                              Eigen::VectorXd &b) const = 0;
    virtual void WhitenSystem(Eigen::MatrixXd &A1, Eigen::MatrixXd &A2,
                              Eigen::MatrixXd &A3, Eigen::VectorXd &b) const = 0;

    /** in-place whiten, override if can be done more efficiently */
    void whitenInPlace(Eigen::VectorXd &v) const
    {
        v = whiten(v);
    }
    /** in-place unwhiten, override if can be done more efficiently */
    void unwhitenInPlace(Eigen::VectorXd &v) const
    {
        v = unwhiten(v);
    }

    /** in-place whiten, override if can be done more efficiently */
    void whitenInPlace(Eigen::Block<Eigen::VectorXd> &v) const
    {
        v = whiten(v);
    }
    /** in-place unwhiten, override if can be done more efficiently */
    void unwhitenInPlace(Eigen::Block<Eigen::VectorXd> &v) const
    {
        v = unwhiten(v);
    }
};

//---------------------------------------------------------------------------------------

/**
 * Gaussian implements the mathematical model
 *  |R*x|^2 = |y|^2 with R'*R=inv(Sigma)
 * where
 *   y = whiten(x) = R*x
 *   x = unwhiten(x) = inv(R)*y
 * as indeed
 *   |y|^2 = y'*y = x'*R'*R*x
 * Various derived classes are available that are more efficient.
 * The named constructors return a shared_ptr because, when the smart flag is
 * true, the underlying object might be a derived class such as Diagonal.
 */
class GaussianNoiseModel : public SharedNoiseModel
{
protected:
    /** Matrix square root of information matrix (R) */
    Eigen::MatrixXd sqrt_information_;

public:
    /**
    * Return R itself, but note that Whiten(H) is cheaper than R*H
    */
    const Eigen::MatrixXd &thisR() const
    {
        // should never happen
        // if (!sqrt_information_) throw std::runtime_error("Gaussian: has no R
        // matrix"); cout<<"sqrt_information_:"<<sqrt_information_<<endl;
        return sqrt_information_;
    }

public:
    GaussianNoiseModel(int dim = 0) : SharedNoiseModel(dim) {}

    GaussianNoiseModel(const Eigen::MatrixXd &R,bool smart = true);
    // : SharedNoiseModel(R.cols())
    // {

    // }

    /** protected constructor takes square root information matrix */
    //   GaussianNoiseModel(int dim = 1, const Eigen::MatrixXd sqrt_information2)
    //   :SharedNoiseModel(dim), sqrt_information_(sqrt_information2) {
    //   }
    /** protected constructor takes square root information matrix */

public:
    // typedef boost::shared_ptr<Gaussian> shared_ptr;

    virtual ~GaussianNoiseModel();

    GaussianNoiseModel &operator=(const GaussianNoiseModel &rObj)
    {
        SharedNoiseModel::operator=(rObj);
        this->sqrt_information_ = rObj.sqrt_information_;
        return *this;
    }

    /**
    * A Gaussian noise model created by specifying a square root information
    * matrix.
    * @param R The (upper-triangular) square root information matrix
    * @param smart check if can be simplified to derived class
    */
    static GaussianNoiseModel SqrtInformation(const Eigen::MatrixXd &R,
            bool smart = true);

    void setR(const Eigen::MatrixXd R)
    {
        sqrt_information_ = R;
    }

    /**
    * A Gaussian noise model created by specifying an information matrix.
    * @param M The information matrix
    * @param smart check if can be simplified to derived class
    */
    static GaussianNoiseModel Information(const Eigen::MatrixXd &M,
                                          bool smart = true);

    /**
    * A Gaussian noise model created by specifying a covariance matrix.
    * @param covariance The square covariance Matrix
    * @param smart check if can be simplified to derived class
    */
    static GaussianNoiseModel Covariance(const Eigen::MatrixXd &covariance,
                                         bool smart = true);
    GaussianNoiseModel* CovariancePointer(
        const Eigen::MatrixXd &covariance, bool smart = true);

    // virtual void print(const std::string& name) const;
    // virtual bool equals(const SharedNoiseModel& expected, double tol=1e-9)
    // const;
    virtual Eigen::VectorXd sigmas() const;
    virtual Eigen::VectorXd whiten(const Eigen::VectorXd &v) const;
    virtual Eigen::VectorXd unwhiten(const Eigen::VectorXd &v) const;
    /**
    * Multiply a derivative with R (derivative of whiten)
    * Equivalent to whitening each column of the input matrix.
    */
    virtual Eigen::MatrixXd Whiten(const Eigen::MatrixXd &H) const;

    /**
    * Mahalanobis distance v'*R'*R*v = <R*v,R*v>
    */
    virtual double Mahalanobis(const Eigen::VectorXd &v) const;

    virtual double distance(const Eigen::VectorXd &v) const;

    /**
    * In-place version
    */
    virtual void WhitenInPlace(Eigen::MatrixXd &H) const;
    /**
    * In-place version
    */
    virtual void WhitenInPlace(Eigen::Block<Eigen::MatrixXd> H) const;

    /**
    * Whiten a system, in place as well
    */
    virtual void WhitenSystem(std::vector<Eigen::MatrixXd> &A,
                              Eigen::VectorXd &b) const;
    virtual void WhitenSystem(Eigen::MatrixXd &A, Eigen::VectorXd &b) const;
    virtual void WhitenSystem(Eigen::MatrixXd &A1, Eigen::MatrixXd &A2,
                              Eigen::VectorXd &b) const;
    virtual void WhitenSystem(Eigen::MatrixXd &A1, Eigen::MatrixXd &A2,
                              Eigen::MatrixXd &A3, Eigen::VectorXd &b) const;

    /**
    * Apply appropriately weighted QR factorization to the system [A b]
    *               Q'  *   [A b]  =  [R d]
    * Dimensions: (r*m) * m*(n+1) = r*(n+1), where r = min(m,n).
    * This routine performs an in-place factorization on Ab.
    * Below-diagonal elements are set to zero by this routine.
    * @param Ab is the m*(n+1) augmented system matrix [A b]
    * @return Empty SharedDiagonal() noise model: R,d are whitened
    */
    // virtual DiagonalNoiseModel QR(Eigen::MatrixXd& Ab) const;

    /// Return R itself, but note that Whiten(H) is cheaper than R*H
    virtual Eigen::MatrixXd R() const;

    /// Compute information matrix
    Eigen::MatrixXd information() const
    {
        return R().transpose() * R();
    }

    /// Compute covariance matrix
    Eigen::MatrixXd covariance() const
    {
        return information().inverse();
    }

    /*
    private:
    friend class boost::serialization::access;
    template<class ARCHIVE>
    void serialize(ARCHIVE & ar, const unsigned int) {
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Base);
    ar & BOOST_SERIALIZATION_NVP(sqrt_information_);
    }
    */
    //   friend class NonlinearFactor;
    // friend class JacobianFactor;

}; // Gaussian

//---------------------------------------------------------------------------------------

/**
 * A diagonal noise model implements a diagonal covariance matrix, with the
 * elements of the diagonal specified in a Vector.  This class has no public
 * constructors, instead, use the static constructor functions Sigmas etc...
 */
class DiagonalNoiseModel : public GaussianNoiseModel
{
protected:
    /**
    * Standard deviations (sigmas), their inverse and inverse square
    * (weights/precisions) These are all computed at construction: the idea is to
    * use one shared model where computation is done only once, the common use
    * case in many problems.
    */
    Eigen::VectorXd sigmas_, invsigmas_, precisions_;

public:
    /** protected constructor - no initializations */

    int nullmodel; // 0 nonnull,1 null;

    DiagonalNoiseModel();

    DiagonalNoiseModel(const DiagonalNoiseModel &dg)
        : GaussianNoiseModel(dg.sigmas_.size()),
          sigmas_(dg.sigmas_),
          invsigmas_(dg.invsigmas_),
          precisions_(dg.precisions_),
          nullmodel(dg.nullmodel)
    {
        int n = sigmas_.size();
        isUnit_ = true;
        for (int j = 1; j < n; j++)
        {
            if (sigmas_(j) != sigmas_(0))
            {
                isUnit_ = false;
                break;
            }
        }
    }

    /** constructor to allow for disabling initialization of invsigmas */
    DiagonalNoiseModel(const Eigen::VectorXd &sigmas);

    DiagonalNoiseModel &operator=(const DiagonalNoiseModel &rObj)
    {
        GaussianNoiseModel::operator=(rObj);
        this->sigmas_ = rObj.sigmas_;
        this->invsigmas_ = rObj.sigmas_;
        this->precisions_ = rObj.precisions_;
        this->nullmodel = rObj.nullmodel;
        return *this;
    }

public:
    // typedef boost::shared_ptr<Diagonal> shared_ptr;

    virtual ~DiagonalNoiseModel();

    /**
    * A diagonal noise model created by specifying a Vector of sigmas, i.e.
    * standard deviations, the diagonal of the square root covariance matrix.
    */
    static DiagonalNoiseModel Sigmas(const Eigen::VectorXd &sigmas,
                                     bool smart = true);

    /**
    * A diagonal noise model created by specifying a Vector of variances, i.e.
    * i.e. the diagonal of the covariance matrix.
    * @param variances A vector containing the variances of this noise model
    * @param smart check if can be simplified to derived class
    */
    static DiagonalNoiseModel Variances(const Eigen::VectorXd &variances,
                                        bool smart = true);

    static DiagonalNoiseModel* VariancesPointer(const Eigen::VectorXd &variances,
                                        bool smart = true);
    /**
    * A diagonal noise model created by specifying a Vector of precisions, i.e.
    * i.e. the diagonal of the information matrix, i.e., weights
    */
    static DiagonalNoiseModel Precisions(const Eigen::VectorXd &precisions,
                                         bool smart = true)
    {
        return Variances(precisions.array().inverse(), smart);
    }
    static DiagonalNoiseModel* PrecisionsPointer(const Eigen::VectorXd &precisions,
                                         bool smart = true)
    {
        return VariancesPointer(precisions.array().inverse(), smart);
    }

    // virtual void print(const std::string& name) const;
    virtual Eigen::VectorXd sigmas() const;
    virtual Eigen::VectorXd whiten(const Eigen::VectorXd &v) const;
    virtual Eigen::VectorXd unwhiten(const Eigen::VectorXd &v) const;
    virtual Eigen::MatrixXd Whiten(const Eigen::MatrixXd &H) const;
    virtual void WhitenInPlace(Eigen::MatrixXd &H) const;
    virtual void WhitenInPlace(Eigen::Block<Eigen::MatrixXd> H) const;

    /**
    * Return standard deviations (sqrt of diagonal)
    */
    virtual double sigma(int i) const;

    /**
    * Return sqrt precisions
    */
    const Eigen::VectorXd invsigmas() const
    {
        return invsigmas_;
    }
    double invsigma(int i) const
    {
        return invsigmas_(i);
    }

    /**
    * Return precisions
    */
    const Eigen::VectorXd precisions() const
    {
        return precisions_;
    }
    double precision(int i) const
    {
        return precisions_(i);
    }

    /**
    * Return R itself, but note that Whiten(H) is cheaper than R*H
    */
    virtual Eigen::MatrixXd R() const;

    /**
    * Apply appropriately weighted QR factorization to the system [A b]
    *               Q'  *   [A b]  =  [R d]
    * Dimensions: (r*m) * m*(n+1) = r*(n+1), where r = min(m,n).
    * This routine performs an in-place factorization on Ab.
    * Below-diagonal elements are set to zero by this routine.
    * @param Ab is the m*(n+1) augmented system matrix [A b]
    * @return Empty SharedDiagonal() noise model: R,d are whitened
    */
    virtual DiagonalNoiseModel* QR(Eigen::MatrixXd &Ab) const;

    //      friend class NonlinearFactor;
    //  friend class JacobianFactor;
    /*
    private:

    friend class boost::serialization::access;
    template<class ARCHIVE>
    void serialize(ARCHIVE & ar, const unsigned int ) {
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Gaussian);
    ar & BOOST_SERIALIZATION_NVP(sigmas_);
    ar & BOOST_SERIALIZATION_NVP(invsigmas_);
    }*/
}; // Diagonal

//--------------------------------------------------------------------------------------

/**
 * A Constrained constrained model is a specialization of Diagonal which allows
 * some or all of the sigmas to be zero, forcing the error to be zero there.
 * All other Gaussian models are guaranteed to have a non-singular square-root
 * information matrix, but this class is specifically equipped to deal with
 * singular noise models, specifically: whiten will return zero on those
 * components that have zero sigma *and* zero error, unchanged otherwise.
 *
 * While a hard constraint may seem to be a case in which there is infinite
 * error, we do not ever produce an error value of infinity to allow for
 * constraints to actually be optimized rather than self-destructing if not
 * initialized correctly.
 */
class ConstrainedNoiseModel : public DiagonalNoiseModel
{
protected:
    // Sigmas are contained in the base class
    Eigen::VectorXd mu_; ///< Penalty function weight - needs to be large enough
    ///< to dominate soft constraints

    /**
    * protected constructor takes sigmas.
    * prevents any inf values
    * from appearing in invsigmas or precisions.
    * mu set to large default value (1000.0)
    */
    ConstrainedNoiseModel(
        const Eigen::VectorXd &sigmas = Eigen::Matrix<double, 1, 1>::Zero());

    /**
    * Constructor that prevents any inf values
    * from appearing in invsigmas or precisions.
    * Allows for specifying mu.
    */
    ConstrainedNoiseModel(const Eigen::VectorXd &mu,
                          const Eigen::VectorXd &sigmas);

public:
    // typedef boost::shared_ptr<Constrained> shared_ptr;

    virtual ~ConstrainedNoiseModel();

    ConstrainedNoiseModel &operator=(const ConstrainedNoiseModel &rObj)
    {
        DiagonalNoiseModel::operator=(rObj);
        this->mu_ = rObj.mu_;
        return *this;
    }

    /// Return true if a particular dimension is free or constrained
    bool constrained(int i) const;

    /// Access mu as a vector
    const Eigen::VectorXd mu() const
    {
        return mu_;
    }

    /**
    * A diagonal noise model created by specifying a Vector of
    * standard devations, some of which might be zero
    */
    static ConstrainedNoiseModel MixedSigmas(const Eigen::VectorXd &mu,
            const Eigen::VectorXd &sigmas);
    static ConstrainedNoiseModel* MixedSigmasPointer(const Eigen::VectorXd &mu,
            const Eigen::VectorXd &sigmas);
    /**
    * A diagonal noise model created by specifying a Vector of
    * standard devations, some of which might be zero
    */
    static ConstrainedNoiseModel MixedSigmas(const Eigen::VectorXd &sigmas)
    {
        return MixedSigmas(Eigen::VectorXd::Constant(sigmas.size(), 1000.0),
                           sigmas);
    }
     static ConstrainedNoiseModel* MixedSigmasPointer(const Eigen::VectorXd &sigmas)
    {
        return  MixedSigmasPointer(Eigen::VectorXd::Constant(sigmas.size(), 1000.0),
                           sigmas);
    }

    /**
    * A diagonal noise model created by specifying a Vector of
    * standard devations, some of which might be zero
    */
    static ConstrainedNoiseModel MixedSigmas(double m,
            const Eigen::VectorXd &sigmas)
    {
        return MixedSigmas(Eigen::VectorXd::Constant(sigmas.size(), m), sigmas);
    }

    /**
    * A diagonal noise model created by specifying a Vector of
    * standard devations, some of which might be zero
    */
    static ConstrainedNoiseModel MixedVariances(
        const Eigen::VectorXd &mu, const Eigen::VectorXd &variances)
    {
        ConstrainedNoiseModel anewone(mu, variances.cwiseSqrt());
        return anewone;
        // return &(new ConstrainedNoiseModel(mu, variances.cwiseSqrt()));
    }
     static ConstrainedNoiseModel* MixedVariancesPointer(
        const Eigen::VectorXd &mu, const Eigen::VectorXd &variances)
    {
        ConstrainedNoiseModel* anewone=new ConstrainedNoiseModel(mu, variances.cwiseSqrt());
        return anewone;
        // return &(new ConstrainedNoiseModel(mu, variances.cwiseSqrt()));
    }
    static ConstrainedNoiseModel MixedVariances(
        const Eigen::VectorXd &variances)
    {
        ConstrainedNoiseModel anewone(variances.cwiseSqrt());
        return anewone;
        // return &(new ConstrainedNoiseModel(variances.cwiseSqrt()));
    }
      static ConstrainedNoiseModel* MixedVariancesPointer(
        const Eigen::VectorXd &variances)
    {
        ConstrainedNoiseModel* anewone=new ConstrainedNoiseModel(variances.cwiseSqrt());
        return anewone;
        // return &(new ConstrainedNoiseModel(variances.cwiseSqrt()));
    }

    /**
    * A diagonal noise model created by specifying a Vector of
    * precisions, some of which might be inf
    */
    static ConstrainedNoiseModel MixedPrecisions(
        const Eigen::VectorXd &mu, const Eigen::VectorXd &precisions)
    {
        return MixedVariances(mu, precisions.array().inverse());
    }
    static ConstrainedNoiseModel MixedPrecisions(
        const Eigen::VectorXd &precisions)
    {
        return MixedVariances(precisions.array().inverse());
    }
    static ConstrainedNoiseModel* MixedPrecisionsPointer(
        const Eigen::VectorXd &mu, const Eigen::VectorXd &precisions)
    {
        return MixedVariancesPointer(mu, precisions.array().inverse());
    }
    static ConstrainedNoiseModel* MixedPrecisionsPointer(
        const Eigen::VectorXd &precisions)
    {
        return MixedVariancesPointer(precisions.array().inverse());
    }

    /**
    * The distance function for a constrained noisemodel,
    * for non-constrained versions, uses sigmas, otherwise
    * uses the penalty function with mu
    */
    virtual double distance(const Eigen::VectorXd &v) const;

    /** Fully constrained variations
    static ConstrainedNoiseModel *All(int dim)
    {
        ConstrainedNoiseModel anewone(Eigen::VectorXd::Constant(dim, 1000.0),
                                      Eigen::VectorXd::Constant(dim, 0));
        return anewone;
        // return &(new ConstrainedNoiseModel(Eigen::VectorXd::Constant(dim,
        // 1000.0), Eigen::VectorXd::Constant(dim,0)));
    }*/

    /** Fully constrained variations */
    static ConstrainedNoiseModel All(int dim, const Eigen::VectorXd &mu)
    {
        ConstrainedNoiseModel anewone(mu, Eigen::VectorXd::Constant(dim, 0));
        return anewone;
        // return ConstrainedNoiseModel*(new ConstrainedNoiseModel(mu,
        // Eigen::VectorXd::Constant(dim,0)));
    }

    /** Fully constrained variations with a mu parameter */
    static ConstrainedNoiseModel All(int dim, double mu)
    {
        ConstrainedNoiseModel anewone(Eigen::VectorXd::Constant(dim, mu),
                                      Eigen::VectorXd::Constant(dim, 0));
        return anewone;
        // return ConstrainedNoiseModel*(new
        // ConstrainedNoiseModel(Eigen::VectorXd::Constant(dim, mu),
        // Eigen::VectorXd::Constant(dim,0)));
    }

    // virtual void print(const std::string& name) const;

    /// Calculates error vector with weights applied
    virtual Eigen::VectorXd whiten(const Eigen::VectorXd &v) const;

    /// Whitening functions will perform partial whitening on rows
    /// with a non-zero sigma.  Other rows remain untouched.
    virtual Eigen::MatrixXd Whiten(const Eigen::MatrixXd &H) const;
    virtual void WhitenInPlace(Eigen::MatrixXd &H) const;
    virtual void WhitenInPlace(Eigen::Block<Eigen::MatrixXd> H) const;

    /**
    * Apply QR factorization to the system [A b], taking into account constraints
    *               Q'  *   [A b]  =  [R d]
    * Dimensions: (r*m) * m*(n+1) = r*(n+1), where r = min(m,n).
    * This routine performs an in-place factorization on Ab.
    * Below-diagonal elements are set to zero by this routine.
    * @param Ab is the m*(n+1) augmented system matrix [A b]
    * @return diagonal noise model can be all zeros, mixed, or not-constrained
    */
    virtual DiagonalNoiseModel* QR(Eigen::MatrixXd &Ab) const;

    /**
    * Returns a Unit version of a constrained noisemodel in which
    * constrained sigmas remain constrained and the rest are unit scaled
    */
    DiagonalNoiseModel unit() const;

    // private:
    /** Serialization function */
    //   friend class boost::serialization::access;
    //   template<class ARCHIVE>
    //   void serialize(ARCHIVE & ar, const unsigned int /*version*/) {
    //      ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Diagonal);
    //      ar & BOOST_SERIALIZATION_NVP(mu_);
    //   }

}; // Constrained

//---------------------------------------------------------------------------------------

/**
 * An isotropic noise model corresponds to a scaled diagonal covariance
 * To construct, use one of the static methods
 */

class IsotropicNoiseModel : public DiagonalNoiseModel
{
protected:
    double dsigma_, dinvsigma_;

public:
    IsotropicNoiseModel(int dim, double sigma)
        : DiagonalNoiseModel(Eigen::VectorXd::Constant(dim, sigma)),
          dsigma_(sigma),
          dinvsigma_(1.0 / sigma)
    {
        setdim(dim);
    }

    IsotropicNoiseModel()
        : DiagonalNoiseModel(Eigen::Matrix<double, 1, 1>::Constant(1.0)),
          dsigma_(1.0),
          dinvsigma_(1.0) {}

    virtual ~IsotropicNoiseModel();

    IsotropicNoiseModel &operator=(const IsotropicNoiseModel &rObj)
    {
        DiagonalNoiseModel::operator=(rObj);
        this->dsigma_ = rObj.dsigma_;
        this->dinvsigma_ = rObj.dinvsigma_;
        return *this;
    }

    // typedef boost::shared_ptr<Isotropic> shared_ptr;

    static IsotropicNoiseModel Sigma(int dim, double sigma, bool smart = true);

     IsotropicNoiseModel* SigmaPointer(int dim, double sigma, bool smart = true);

    static IsotropicNoiseModel Variance(int dim, double variance,
                                        bool smart = true);
    static IsotropicNoiseModel* VariancePointer(int dim, double variance,
                                        bool smart = true);

    static IsotropicNoiseModel Precision(int dim, double precision,
                                         bool smart = true)
    {
        return Variance(dim, 1.0 / precision, smart);
    }

    // virtual void print(const std::string& name) const;
    virtual double Mahalanobis(const Eigen::VectorXd &v) const;
    virtual Eigen::VectorXd whiten(const Eigen::VectorXd &v) const;
    virtual Eigen::VectorXd unwhiten(const Eigen::VectorXd &v) const;
    virtual Eigen::MatrixXd Whiten(const Eigen::MatrixXd &H) const;
    virtual void WhitenInPlace(Eigen::MatrixXd &H) const;
    virtual void WhitenInPlace(Eigen::Block<Eigen::MatrixXd> H) const;

    double sigma() const
    {
        return dsigma_;
    }

    // private:

    // friend class boost::serialization::access;
    // template<class ARCHIVE>
    // void serialize(ARCHIVE & ar, const unsigned int ) {
    // ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Diagonal);
    // ar & BOOST_SERIALIZATION_NVP(sigma_);
    // ar & BOOST_SERIALIZATION_NVP(invsigma_);
    // }
};

//---------------------------------------------------------------------------------------

/**
 * Unit: i.i.d. unit-variance noise on all m dimensions.
 */
/*
    class  UnitNoiseModel : public IsotropicNoiseModel {
    protected:

      UnitNoiseModel(int dim=1): IsotropicNoiseModel(dim,1.0) {}

    public:

    //  typedef boost::shared_ptr<Unit> shared_ptr;

      virtual ~UnitNoiseModel() {}


       * Create a unit covariance noise model

      static UnitNoiseModel Create(int dim) {
          UnitNoiseModel anewone(dim);
          cout<<"anewone.thisR():"<<endl<<anewone.thisR()<<endl;

          return anewone;
       // return shared_ptr(new Unit(dim));
      }

      //virtual void print(const std::string& name) const;
      virtual double Mahalanobis(const Eigen::VectorXd& v) const {return
   v.dot(v); } virtual Eigen::VectorXd whiten(const Eigen::VectorXd& v) const {
   return v; } virtual Eigen::VectorXd unwhiten(const Eigen::VectorXd& v) const
   { return v; } virtual Eigen::MatrixXd Whiten(const Eigen::MatrixXd& H) const
   { return H; } virtual void WhitenInPlace(Eigen::MatrixXd& ) const {} virtual
   void WhitenInPlace(Eigen::Block<Eigen::MatrixXd> ) const {} virtual void
   whitenInPlace(Eigen::VectorXd& ) const {} virtual void
   unwhitenInPlace(Eigen::VectorXd& ) const {} virtual void
   whitenInPlace(Eigen::Block<Eigen::VectorXd>& ) const {} virtual void
   unwhitenInPlace(Eigen::Block<Eigen::VectorXd>& ) const {}

   // private:
    //  friend class boost::serialization::access;
    //  template<class ARCHIVE>
    //  void serialize(ARCHIVE & ar, const unsigned int ) {
   //     ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Isotropic);
   //   }
    };
*/
GaussianNoiseModel* newGaussianNoiseModelCovariance(
    const Eigen::MatrixXd &covariance);
GaussianNoiseModel* NGaussianNoiseModel(const Eigen::MatrixXd &cov);

DiagonalNoiseModel* NDiagonalNoiseModelPrecision(const Eigen::VectorXd& precisions,bool smart=true);
 DiagonalNoiseModel* DiagonalNoiseModelSigmasPointer(const Eigen::VectorXd &sigmas,
                                     bool smart = true);
#endif // NOISEMODEL_H_INCLUDED

