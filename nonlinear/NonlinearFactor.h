#ifndef NONLINEARFACTOR_H
#define NONLINEARFACTOR_H


/**
 * @file    NonlinearFactor.h
 * @brief   Non-linear factor base classes
 */

#include "../inference/Factor.h"
#include "../linear/RealGaussianFactor.h"
#include "../linear/NoiseModel.h"
#include "../geometry/Pose3.h"
#include "../geometry/Pose2.h"
#include "../gmfconfig.h"


namespace minisam
{


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

class NoiseModelFactor
{
public:
    /// Iterator over keys
    typedef std::vector<int>::iterator iterator;
    /// Const iterator over keys
    typedef std::vector<int>::const_iterator const_iterator;

    /// The keys involved in this factor
    std::vector<int> keys_;

    GaussianNoiseModel *noiseModel_; /** Noise model */

public:
    /** Default constructor for I/O only */
    NoiseModelFactor();

    NoiseModelFactor(const std::vector<int>& keys);

    /** Destructor */
    virtual ~NoiseModelFactor();

    /**
    * Constructor
    */
    NoiseModelFactor(GaussianNoiseModel *noiseModel, const std::vector<int>& keys);

protected:
    /**
    * Constructor - only for subclasses, as this does not set keys.
    */
    NoiseModelFactor(GaussianNoiseModel *noiseModel);

public:
    /** get the dimension of the factor (number of rows on linearization) */
    virtual int dim() const;

    /// Access the factor's involved variable keys
    const std::vector<int>& keys() const;

    /// Access the factor's involved variable keys
    std::vector<int>& keys();

    /** Iterator at beginning of involved variable keys */
    std::vector<int>::const_iterator begin() const;

    /** Iterator at end of involved variable keys */
    std::vector<int>::const_iterator end() const;


    /** Iterator at beginning of involved variable keys */
    std::vector<int>::iterator begin();

    /** Iterator at end of involved variable keys */
    std::vector<int>::iterator end();

    bool empty();


    /**
    * @return the number of variables involved in this factor
    */
    int size() const;

    /// access to the noise model
    GaussianNoiseModel *noiseModel() const;

    /**
    * Error function *without* the NoiseModel, \f$ z-h(x) \f$.
    * Override this method to finish implementing an N-way factor.
    * If the optional arguments is specified, it should compute
    * both the function evaluation and its derivative(s) in H.
    */
    virtual minivector unwhitenedError(const std::map<int, minimatrix*>& x) const;
    virtual minivector unwhitenedError(const std::map<int, minimatrix*>& x,std::vector<minimatrix> &H) const;

    /**
    * Calculate the error of the factor.
    * This is the log-likelihood, e.g. \f$ 0.5(h(x)-z)^2/\sigma^2 \f$ in case of Gaussian.
    * In this class, we take the raw prediction error \f$ h(x)-z \f$, ask the noise model
    * to transform it to \f$ (h(x)-z)^2/\sigma^2 \f$, and then multiply by 0.5.
    */
    virtual double error(const std::map<int, minimatrix*>& c) const;

    /**
    * Linearize a non-linearFactor to get a GaussianFactor,
    * \f$ Ax-b \approx h(x+\delta x)-z = h(x) + A \delta x - z \f$
    * Hence \f$ b = z - h(x) = - \mathtt{error\_vector}(x) \f$
    */
    virtual RealGaussianFactor*  linearize(const std::map<int, minimatrix*>& x,int factorizaton=0) const;
    virtual NoiseModelFactor *clone() const;

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

class NoiseModelFactor1 : public NoiseModelFactor
{

public:
    /** Default constructor for I/O only */
    NoiseModelFactor1() : NoiseModelFactor() {}

    virtual ~NoiseModelFactor1() {}

    inline int key() const
    {
        return keys_.front();
    }

    /**
    *  Constructor
    *  @param noiseModel  pointer to noise model
    *  @param key1 by which to look up X value in Values
    */
    NoiseModelFactor1(GaussianNoiseModel *noiseModel, int key1) : NoiseModelFactor(noiseModel)
    {
        std::vector<int> key;
        keys_.push_back(key1);
    }

    /** Calls the 1-key specific version of evaluateError, which is pure virtual
    *  so must be implemented in the derived class.
    */

    virtual minivector unwhitenedError(const std::map<int, minimatrix*>& x) const
    {
        std::map<int, minimatrix*>::const_iterator xbegin = x.find(keys_[0]);
        return evaluateError(xbegin->second);
    }
    virtual minivector unwhitenedError(const std::map<int, minimatrix*>& x,std::vector<minimatrix> &H) const
    {
        std::map<int, minimatrix*>::const_iterator xbegin = x.find(keys_[0]);
        return evaluateError(xbegin->second, H.front());
    }
    virtual minivector evaluateError(const minimatrix* x) const
    {
        minivector xb;
        return xb;
    }
    virtual minivector evaluateError(const minimatrix* x, minimatrix &H) const
    {
        minivector xb;
        return xb;
    }


}; // \class NoiseModelFactor1

/* ************************************************************************* */
/** A convenient base class for creating your own NoiseModelFactor with 2
 * variables.  To derive from this class, implement evaluateError(). */
class NoiseModelFactor2 : public NoiseModelFactor
{
public:
    /**
    * Default Constructor for I/O
    */
    NoiseModelFactor2( ) : NoiseModelFactor() {}

    /**
    * Constructor
    * @param noiseModel  pointer to noise model
    * @param j1 key of the first variable
    * @param j2 key of the second variable
    */
    NoiseModelFactor2(GaussianNoiseModel *noiseModel, int j1, int j2) : NoiseModelFactor(noiseModel)
    {
        keys_.reserve(2);
        keys_.push_back(j1);
        keys_.push_back(j2);
    }

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

    /**
    *  Override this method to finish implementing a binary factor.
    *  If any of the optional Matrix reference arguments are specified, it should compute
    *  both the function evaluation and its derivative(s) in X1 (and/or X2).
    */

    virtual minivector unwhitenedError(const std::map<int, minimatrix*>& x) const
    {
        std::map<int, minimatrix*>::const_iterator xbegin = x.find(key1());
        std::map<int, minimatrix*>::const_iterator xsecond = x.find(key2());
        return evaluateError(xbegin->second, xsecond->second);
    }
    virtual minivector unwhitenedError(const std::map<int, minimatrix*>& x, std::vector<minimatrix> &H) const
    {
        std::map<int, minimatrix*>::const_iterator xbegin = x.find(key1());
        std::map<int, minimatrix*>::const_iterator xsecond = x.find(key2());

        return evaluateError(xbegin->second, xsecond->second, *(H.begin()), *(H.begin() + 1));
    }
    /**
    *  Override this method to finish implementing a binary factor.
    *  If any of the optional Matrix reference arguments are specified, it should compute
    *  both the function evaluation and its derivative(s) in X1 (and/or X2).
    */
    virtual minivector
    evaluateError(const minimatrix* X1, const minimatrix* X2) const
    {
        minivector VV;
        return VV;
    }

    virtual minivector
    evaluateError(const minimatrix* X1, const minimatrix* X2, minimatrix &H1, minimatrix &H2) const
    {
        minivector VV;
        return VV;
    }


}; // \class NoiseModelFactor2

/** A convenient base class for creating your own NoiseModelFactor with 3
 * variables.  To derive from this class, implement evaluateError(). */
class NoiseModelFactor3: public NoiseModelFactor {

public:

  /**
   * Default Constructor for I/O
   */
  NoiseModelFactor3() : NoiseModelFactor() {}

  /**
   * Constructor
   * @param noiseModel shared pointer to noise model
   * @param j1 key of the first variable
   * @param j2 key of the second variable
   * @param j3 key of the third variable
   */
  NoiseModelFactor3(GaussianNoiseModel* noiseModel, int j1, int j2, int j3): NoiseModelFactor(noiseModel)
    {
      keys_.reserve(3);
      keys_.push_back(j1);keys_.push_back(j2);keys_.push_back(j3);
    }

  virtual ~NoiseModelFactor3() {}

  /** methods to retrieve keys */
  inline int key1() const { return keys_[0]; }
  inline int key2() const { return keys_[1]; }
  inline int key3() const { return keys_[2]; }

  /** Calls the 3-key specific version of evaluateError, which is pure virtual
   * so must be implemented in the derived class. */
    virtual minivector unwhitenedError(const std::map<int, minimatrix*>& x) const
    {
        std::map<int, minimatrix*>::const_iterator xbegin = x.find(key1());
        std::map<int, minimatrix*>::const_iterator xsecond = x.find(key2());
        std::map<int, minimatrix*>::const_iterator xthird = x.find(key3());
        return evaluateError(xbegin->second, xsecond->second,xthird->second);
    }
    virtual minivector unwhitenedError(const std::map<int, minimatrix*>& x, std::vector<minimatrix> &H) const
    {
        std::map<int, minimatrix*>::const_iterator xbegin = x.find(key1());
        std::map<int, minimatrix*>::const_iterator xsecond = x.find(key2());
        std::map<int, minimatrix*>::const_iterator xthird = x.find(key3());

        return evaluateError(xbegin->second, xsecond->second,xthird->second, *(H.begin()), *(H.begin() + 1), *(H.begin() + 2));
    }

  /**
   *  Override this method to finish implementing a trinary factor.
   *  If any of the optional Matrix reference arguments are specified, it should compute
   *  both the function evaluation and its derivative(s) in X1 (and/or X2, X3).
   */
   virtual minivector
    evaluateError(const minimatrix* X1, const minimatrix* X2, const minimatrix* X3) const
    {
        minivector VV;
        return VV;
    }

    virtual minivector
    evaluateError(const minimatrix* X1, const minimatrix* X2, const minimatrix* X3,
                  minimatrix &H1, minimatrix &H2, minimatrix &H3) const
    {
        minivector VV;
        return VV;
    }

}; // \class NoiseModelFactor3


/* ************************************************************************* */
void check(GaussianNoiseModel *noiseModel, int m);

};
#endif // NONLINEARFACTOR_H
