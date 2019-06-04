#ifndef NONLINEARFACTOR_H
#define NONLINEARFACTOR_H

/**
 * @file    NonlinearFactor.h
 * @brief   Non-linear factor base classes
 * @author
 */

#include "../inference/Factor.h"
#include "../linear/RealGaussianFactor.h"
#include "../linear/NoiseModel.h"
#include "../geometry/Pose3.h"
#include "../geometry/Pose2.h"
#include "../gmfconfig.h"
//#include "../linear/HessianFactor.h"
//#include "../linear/JacobianFactor.h"
/**
 * Macro to add a standard clone function to a derived factor
 * @deprecated: will go away shortly - just add the clone function directly
 */

/* ************************************************************************* */

/**
 * Nonlinear factor base class
 *
 * \nosubgrouping
 */

//class NoiseModelFactor;
//class HessianFactor;
//class JacobianFactor;
class NonlinearFactor : public Factor
{
public:
    int factortype;

    /// @name Standard Constructors
    /// @{

    /** Default constructor for I/O only */
    NonlinearFactor(int factort = 0) : factortype(factort) {}

    /**
    * Constructor from a collection of the keys involved in this factor
    */
    //template<typename CONTAINER>
    NonlinearFactor(const std::vector<int>& keys, int factort = 0) : Factor(keys), factortype(factort) {}

    /** Destructor */
    ~NonlinearFactor() {}

    /**
    * Calculate the error of the factor
    * This is typically equal to log-likelihood, e.g. \f$ 0.5(h(x)-z)^2/sigma^2 \f$ in case of Gaussian.
    * You can override this for systems with unusual noise models.
    */
    virtual double error(const std::map<int, Eigen::VectorXd>& c) const
    {
        ///It is wrong. Fix it later.
        return 0;
    }

//#ifdef GMF_Using_Pose3
    virtual double errorPose3(const std::map<int, Pose3>& c) const
    {
        ///It is wrong. Fix it later.
        return 0;
    }
    virtual double LCFPerrorP3v(const std::map<int, Pose3>& p3, const std::map<int, Eigen::VectorXd>& v) const
    {
        ///It is wrong. Fix it later.
        return 0;
    }

    virtual double errorP3v(const std::map<int, Pose3>& p3, const std::map<int, Eigen::VectorXd>& v) const
    {
        ///It is wrong. Fix it later.
        return 0;
    }
//#else
    virtual double errorPose2(const std::map<int, Pose2>& c) const
    {
        ///It is wrong. Fix it later.
        return 0;
    }

    virtual double LCFPerrorP2v(const std::map<int, Pose2>& p3, const std::map<int, Eigen::VectorXd>& v) const
    {
        ///It is wrong. Fix it later.
        return 0;
    }

    virtual double errorP2v(const std::map<int, Pose2>& p3, const std::map<int, Eigen::VectorXd>& v) const
    {
        ///It is wrong. Fix it later.
        return 0;
    }
    virtual double errorP2v(const std::map<int, Pose2>& p3, const std::map<int, Eigen::Vector2d>& v) const
    {
        ///It is wrong. Fix it later.
        return 0;
    }
//#endif

    /** get the dimension of the factor (number of rows on linearization) */
    int dim()
    {
        return keys_.size();
    }

    /**
    * Checks whether a factor should be used based on a set of values.
    * This is primarily used to implment inequality constraints that
    * require a variable active set. For all others, the default implementation
    * returning true solves this problem.
    *
    * In an inequality/bounding constraint, this active() returns true
    * when the constraint is *NOT* fulfilled.
    * @return true if the constraint is active
    */
    virtual bool active(const std::map<int, Eigen::VectorXd>& c /*c*/) const
    {
        return true;
    }
    //
    /** linearize to a GaussianFactor*/

    //implement after jacobian factor implemented!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!!!!!!!!!!!!!!!!!
    virtual RealGaussianFactor
    linearizeVector(const std::map<int, Eigen::VectorXd>& c) const
    {
        RealGaussianFactor ng;
        return ng;
    }
    virtual RealGaussianFactor*
    linearizeVectorPointer(const std::map<int, Eigen::VectorXd>& c,int factorization) const
    {
        //RealGaussianFactor* ng=new RealGaussianFactor();
        RealGaussianFactor* ng=new RealGaussianFactor();

        return ng;
    }
//#ifdef GMF_Using_Pose3
    virtual RealGaussianFactor linearizePose(const std::map<int, Pose3>& x) const
    {
        RealGaussianFactor ng;
        return ng;
    }
    virtual RealGaussianFactor* linearizePosePointer(const std::map<int, Pose3>& x,int factorization) const
    {
        RealGaussianFactor* ng=new RealGaussianFactor();
        return ng;
    }

    virtual RealGaussianFactor linearizePV(const std::map<int, Pose3>& x1,
                                           const std::map<int, Eigen::VectorXd>& x2) const
    {
        RealGaussianFactor ng;
        return ng;
    }
    virtual RealGaussianFactor* linearizePVPointer(const std::map<int, Pose3>& x1,
            const std::map<int, Eigen::VectorXd>& x2,int factorization) const
    {
        RealGaussianFactor* ng=new RealGaussianFactor();
        return ng;
    }
//#else
    virtual RealGaussianFactor linearizePose(const std::map<int, Pose2>& x) const
    {
        RealGaussianFactor ng;
        return ng;
    }

    virtual RealGaussianFactor* linearizePosePointer(const std::map<int, Pose2>& x,int factorization) const
    {
        RealGaussianFactor* ng=new RealGaussianFactor();
        return ng;
    }


    virtual RealGaussianFactor linearizePV(const std::map<int, Pose2>& x1,
                                           const std::map<int, Eigen::VectorXd>& x2) const
    {
        RealGaussianFactor ng;
        return ng;
    }
    virtual RealGaussianFactor* linearizePVPointer(const std::map<int, Pose2>& x1,
            const std::map<int, Eigen::VectorXd>& x2,int factorization) const
    {
         RealGaussianFactor* ng=new RealGaussianFactor();
        return ng;
    }
    virtual RealGaussianFactor linearizePV(const std::map<int, Pose2>& x1,
                                           const std::map<int, Eigen::Vector2d>& x2) const
    {
        RealGaussianFactor ng;
        return ng;
    }

    virtual RealGaussianFactor* linearizePVPointer(const std::map<int, Pose2>& x1,
            const std::map<int, Eigen::Vector2d>& x2,int factorization) const
    {
        RealGaussianFactor* ng=new RealGaussianFactor();
        return ng;
    }

    virtual NonlinearFactor *clone() const
    {
        NonlinearFactor *newfactor = new NonlinearFactor(keys_);
        return newfactor;
    }
    /**
    * Creates a shared_ptr clone of the
    * factor with different keys using
    * a map from old->new keys
    */
    //NonlinearFactor rekey(const std::map<int,int> rekey_mapping) const;

    /**
    * Clones a factor and fully replaces its keys
    * @param new_keys is the full replacement set of keys
    */
    //NonlinearFactor rekey(const std::vector<int> new_keys) const;
}; // \class NonlinearFactor

/* ************************************************************************* */
/**
 * A nonlinear sum-of-squares factor with a zero-mean noise model
 * implementing the density \f$ P(z|x) \propto exp -0.5*|z-h(x)|^2_C \f$
 * Templated on the parameter type X and the values structure Values
 * There is no return type specified for h(x). Instead, we require
 * the derived class implements \f$ \mathtt{error\_vector}(x) = h(x)-z \approx A \delta x - b \f$
 * This allows a graph to have factors with measurements of mixed type.

 * The noise model is typically Gaussian, but robust and constrained error models are also supported.
 */
class NoiseModelFactor : public NonlinearFactor
{

protected:
    SharedNoiseModel *noiseModel_; /** Noise model */
    //int factortype;//0:vector,1:pose,2:pose vector

public:
    /** Default constructor for I/O only */
    NoiseModelFactor(int factort = 0) : NonlinearFactor(factort) {}

    /** Destructor */
    virtual ~NoiseModelFactor()
    {
    if(noiseModel_!=NULL)
    {
    delete noiseModel_;
    noiseModel_=NULL;
    }
    }

    /**
    * Constructor
    */
    //template<typename CONTAINER>
    NoiseModelFactor(SharedNoiseModel *noiseModel, const std::vector<int>& keys, int factort = 0) : noiseModel_(noiseModel), NonlinearFactor(keys, factort) {}

protected:
    /**
    * Constructor - only for subclasses, as this does not set keys.
    */
    NoiseModelFactor(SharedNoiseModel *noiseModel, int factort = 0) : noiseModel_(noiseModel), NonlinearFactor(factort) {}

public:
    /** get the dimension of the factor (number of rows on linearization) */
    virtual int dim() const
    {
        return noiseModel_->dim();
    }

    /// access to the noise model
    SharedNoiseModel *noiseModel() const
    {
        return noiseModel_;
    }

    /**
    * Error function *without* the NoiseModel, \f$ z-h(x) \f$.
    * Override this method to finish implementing an N-way factor.
    * If the optional arguments is specified, it should compute
    * both the function evaluation and its derivative(s) in H.
    */
    virtual Eigen::VectorXd unwhitenedError(const std::map<int, Eigen::VectorXd>& x) const
    {
        Eigen::VectorXd f;
        return f;
    }
    virtual Eigen::VectorXd unwhitenedError(const std::map<int, Eigen::VectorXd>& x,
                                            std::vector<Eigen::MatrixXd> &H) const
    {
        Eigen::VectorXd f;
        return f;
    }

///#ifdef GMF_Using_Pose3
    virtual Eigen::VectorXd unwhitenedError(const std::map<int, Pose3>& x) const
    {
        Eigen::VectorXd f;
        return f;
    }
    virtual Eigen::VectorXd unwhitenedError(const std::map<int, Pose3>& x, std::vector<Eigen::MatrixXd> &H) const
    {
        Eigen::VectorXd f;
        return f;
    }
    virtual Eigen::VectorXd unwhitenedError(const std::map<int, Pose3>& x1,
                                            const std::map<int, Eigen::VectorXd>& x2,
                                            std::vector<Eigen::MatrixXd> &H) const
    {
        Eigen::VectorXd f;
        return f;
    }

    virtual Eigen::VectorXd unwhitenedError(const std::map<int, Pose3>& x1,
                                            const std::map<int, Eigen::VectorXd>& x2) const
    {
        Eigen::VectorXd f;
        return f;
    }
//#else
    virtual Eigen::VectorXd unwhitenedError(const std::map<int, Pose2>& x) const
    {
        Eigen::VectorXd f;
        return f;
    }
    virtual Eigen::VectorXd unwhitenedError(const std::map<int, Pose2>& x, std::vector<Eigen::MatrixXd>& H) const
    {
        Eigen::VectorXd f;
        return f;
    }
    virtual Eigen::VectorXd unwhitenedError(const std::map<int, Pose2>& x1,
                                            const std::map<int, Eigen::VectorXd>& x2,
                                            std::vector<Eigen::MatrixXd> &H) const
    {
        Eigen::VectorXd f;
        return f;
    }

    virtual Eigen::VectorXd unwhitenedError(const std::map<int, Pose2>& x1,
                                            const std::map<int, Eigen::VectorXd>& x2) const
    {
        Eigen::VectorXd f;
        return f;
    }

    virtual Eigen::VectorXd unwhitenedError(const std::map<int, Pose2>& x1,const std::map<int,Eigen::Vector2d>& x2) const
    {
        Eigen::VectorXd f;
        return f;
    }

//#endif // GMF_Using_Pose3

    /**
    * Vector of errors, whitened
    * This is the raw error, i.e., i.e. \f$ (h(x)-z)/\sigma \f$ in case of a Gaussian
    */
    Eigen::VectorXd whitenedError(const std::map<int, Eigen::VectorXd>& c) const;

    /**
    * Calculate the error of the factor.
    * This is the log-likelihood, e.g. \f$ 0.5(h(x)-z)^2/\sigma^2 \f$ in case of Gaussian.
    * In this class, we take the raw prediction error \f$ h(x)-z \f$, ask the noise model
    * to transform it to \f$ (h(x)-z)^2/\sigma^2 \f$, and then multiply by 0.5.
    */
    virtual double error(const std::map<int, Eigen::VectorXd>& c) const;

#ifdef GMF_Using_Pose3
    virtual double errorPose3(const std::map<int, Pose3>& c) const;
    virtual double errorP3v(const std::map<int, Pose3>& p3, const std::map<int, Eigen::VectorXd>& v) const;
#else
    virtual double errorPose2(const std::map<int, Pose2>& c) const;
    virtual double errorP2v(const std::map<int, Pose2>& p2, const std::map<int, Eigen::VectorXd>& v) const;
    virtual double errorP2v(const std::map<int, Pose2>& p2, const std::map<int, Eigen::Vector2d>& v) const;
#endif // GMF_Using_Pose3

    /**
    * Linearize a non-linearFactorN to get a GaussianFactor,
    * \f$ Ax-b \approx h(x+\delta x)-z = h(x) + A \delta x - z \f$
    * Hence \f$ b = z - h(x) = - \mathtt{error\_vector}(x) \f$
    */

    virtual RealGaussianFactor linearizeVector(const std::map<int, Eigen::VectorXd>& x) const override;
    virtual RealGaussianFactor linearizePose(const std::map<int, Pose3>& x) const override;
    virtual RealGaussianFactor linearizePV(const std::map<int, Pose3>& x1,
                                           const std::map<int, Eigen::VectorXd>& x2) const override;
    virtual RealGaussianFactor linearizePose(const std::map<int, Pose2>& x) const override;
    virtual RealGaussianFactor linearizePV(const std::map<int, Pose2>& x1,
                                           const std::map<int, Eigen::VectorXd>& x2) const override;
    virtual RealGaussianFactor* linearizeVectorPointer(const std::map<int, Eigen::VectorXd>& x,int factorizaton) const override;
    virtual RealGaussianFactor* linearizePosePointer(const std::map<int, Pose3>& x,int factorization) const override;
    virtual RealGaussianFactor* linearizePVPointer(const std::map<int, Pose3>& x1,
            const std::map<int, Eigen::VectorXd>& x2,int factorization) const override;
    virtual RealGaussianFactor* linearizePosePointer(const std::map<int, Pose2>& x,int factorization) const override;
    virtual RealGaussianFactor* linearizePVPointer(const std::map<int, Pose2>& x1,
            const std::map<int, Eigen::VectorXd>& x2,int factorization) const override;
};     // \class NoiseModelFactor

/* ************************************************************************* */

/**
 * A convenient base class for creating your own NoiseModelFactor with 1
 * variable.  To derive from this class, implement evaluateError().
 *
 * Templated on a values structure type. The values structures are typically
 * more general than just vectors, e.g., Rot3 or Pose3,
 * which are objects in non-linear manifolds (Lie groups).
 */
//template<class VALUE>
class NoiseModelFactor1 : public NoiseModelFactor
{

    //int factortype;//0:vector,1:pose

public:
    /** Default constructor for I/O only */
    NoiseModelFactor1(int factort = 0) : NoiseModelFactor(factort) {}

    virtual ~NoiseModelFactor1() {}

    inline int key() const
    {
        return keys_.front();
    }

    /**
    *  Constructor
    *  @param noiseModel shared pointer to noise model
    *  @param key1 by which to look up X value in Values
    */
    NoiseModelFactor1(SharedNoiseModel *noiseModel, int key1, int factort = 0) : NoiseModelFactor(noiseModel, factort)
    {
        std::vector<int> key;
        keys_.push_back(key1);
    }

    /** Calls the 1-key specific version of evaluateError, which is pure virtual
    *  so must be implemented in the derived class.
    */
    virtual Eigen::VectorXd unwhitenedError(const std::map<int, Pose3>& x) const
    {
        std::map<int, Pose3>::const_iterator xbegin = x.find(keys_[0]);
        const Pose3& x1 = xbegin->second;
        return evaluateError(x1);
    }
    virtual Eigen::VectorXd unwhitenedError(const std::map<int, Pose3>& x, std::vector<Eigen::MatrixXd> &H) const
    {
        std::map<int, Pose3>::const_iterator xbegin = x.find(keys_[0]);

        const Pose3& x1 = xbegin->second;
        return evaluateError(x1, H.front());
    }
    virtual Eigen::VectorXd unwhitenedError(const std::map<int, Pose2>& x) const
    {
        cout<<keys()[0]<<endl;
        std::map<int, Pose2>::const_iterator xbegin = x.find(keys_[0]);
        const Pose2& x1 = xbegin->second;
        return evaluateError(x1);
    }
    virtual Eigen::VectorXd unwhitenedError(const std::map<int, Pose2>& x, std::vector<Eigen::MatrixXd> &H) const
    {
        std::map<int, Pose2>::const_iterator xbegin = x.find(keys_[0]);
        const Pose2& x1 = xbegin->second;
        return evaluateError(x1, H.front());
    }
//#endif // GMF_Using_Pose3

    virtual Eigen::VectorXd unwhitenedError(const std::map<int, Eigen::VectorXd>& x) const
    {
        std::map<int, Eigen::VectorXd>::const_iterator xbegin = x.find(keys_[0]);
        Eigen::VectorXd x1 = xbegin->second;
        return evaluateError(x1);
    }

    virtual Eigen::VectorXd unwhitenedError(const std::map<int, Eigen::VectorXd>& x, std::vector<Eigen::MatrixXd> &H) const
    {
        std::map<int, Eigen::VectorXd>::const_iterator xbegin = x.find(keys_[0]);
        Eigen::VectorXd x1 = xbegin->second;
        return evaluateError(x1, H.front());
    }
    /**
    *  Override this method to finish implementing a unary factor.
    *  If the optional Matrix reference argument is specified, it should compute
    *  both the function evaluation and its derivative in X.
    */

    virtual Eigen::VectorXd evaluateError(const Eigen::VectorXd& x) const
    {
        Eigen::VectorXd xb;
        return xb;
    }
    virtual Eigen::VectorXd evaluateError(const Eigen::VectorXd& x, Eigen::MatrixXd &H) const
    {
        Eigen::VectorXd xb;
        return xb;
    }

//#ifdef GMF_Using_Pose3
    virtual Eigen::VectorXd evaluateError(const Pose3& x) const
    {
        Eigen::VectorXd xb;
        return xb;

    }
    virtual Eigen::VectorXd evaluateError(const Pose3& x, Eigen::MatrixXd &H) const
    {
        Eigen::VectorXd xb;
        return xb;
    }
//#else
    virtual Eigen::VectorXd evaluateError(const Pose2& x) const
    {
        Eigen::VectorXd xb;
        return xb;
    }
    virtual Eigen::VectorXd evaluateError(const Pose2& x, Eigen::MatrixXd &H) const
    {
        Eigen::VectorXd xb;
        return xb;
    }
}; // \class NoiseModelFactor1

/* ************************************************************************* */
/** A convenient base class for creating your own NoiseModelFactor with 2
 * variables.  To derive from this class, implement evaluateError(). */
//template<class VALUE1, class VALUE2>
class NoiseModelFactor2 : public NoiseModelFactor
{
public:
    /**
    * Default Constructor for I/O
    */
    NoiseModelFactor2(int factort = 2) : NoiseModelFactor(factort) {}

    /**
    * Constructor
    * @param noiseModel shared pointer to noise model
    * @param j1 key of the first variable
    * @param j2 key of the second variable
    */
    NoiseModelFactor2(SharedNoiseModel *noiseModel, int j1, int j2, int factort = 2) : NoiseModelFactor(noiseModel, factort)
    {
        //std::vector<int> key;
        keys_.reserve(2);
        keys_.push_back(j1);
        keys_.push_back(j2);
        //NoiseModelFactor(noiseModel,key);
    }
    //:Base(noiseModel, cref_list_of<2>(j1)(j2)) {}

    virtual ~NoiseModelFactor2() {}

    /** methods to retrieve both keys */
    inline int key1() const
    {
        return keys_[0];
    }
    inline int key2() const
    {
        return keys_[1];
    }

    virtual Eigen::VectorXd unwhitenedError(const std::map<int, Eigen::VectorXd>& x) const
    {
        std::map<int, Eigen::VectorXd>::const_iterator xbegin = x.find(key1());
        Eigen::VectorXd x1 = xbegin->second;
        std::map<int, Eigen::VectorXd>::const_iterator xsecond = x.find(key2());
        Eigen::VectorXd x2 = xsecond->second;
        return evaluateError(x1, x2);
    }
    virtual Eigen::VectorXd unwhitenedError(const std::map<int, Eigen::VectorXd>& x, std::vector<Eigen::MatrixXd> &H) const
    {
        std::map<int, Eigen::VectorXd>::const_iterator xbegin = x.find(key1());
        Eigen::VectorXd x1 = xbegin->second;
        std::map<int, Eigen::VectorXd>::const_iterator xsecond = x.find(key2());
        Eigen::VectorXd x2 = xsecond->second;

        // return evaluateError(x1, x2);

        return evaluateError(x1, x2, *(H.begin()), *(H.begin() + 1));
    }
    /**
    *  Override this method to finish implementing a binary factor.
    *  If any of the optional Matrix reference arguments are specified, it should compute
    *  both the function evaluation and its derivative(s) in X1 (and/or X2).
    */
    // virtual Eigen::VectorXd
    // evaluateError(const Eigen::MatrixXd&, const X2&, std::optional<Eigen::MatrixXd&> H1 =
    //     std::nullopt, std::optional<Eigen::MatrixXd&> H2 = std::nullopt) const = 0;

    virtual Eigen::VectorXd
    evaluateError(const Eigen::VectorXd& X1, const Eigen::VectorXd& X2) const
    {
        Eigen::VectorXd VV;
        return VV;
    }

    virtual Eigen::VectorXd
    evaluateError(const Eigen::VectorXd& X1, const Eigen::VectorXd& X2, Eigen::MatrixXd &H1, Eigen::MatrixXd &H2) const
    {
        Eigen::VectorXd VV;
        return VV;
    }

}; // \class NoiseModelFactor2

/* ************************************************************************* */
static void check(SharedNoiseModel *noiseModel, int m);
#endif // NONLINEARFACTOR_H
