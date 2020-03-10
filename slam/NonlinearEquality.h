/*
 * @file NonlinearEquality.h
 * @brief Factor to handle enforced equality between factors
 */

#pragma once

#include "../nonlinear/NonlinearFactor.h"
#include "../linear/JacobianFactor.h"

#include <limits>
#include <iostream>
#include <cmath>

namespace minisam {

/**
 * An equality factor that forces either one variable to a constant,
 * or a set of variables to be equal to each other.
 *
 * Depending on flag, throws an error at linearization if the constraints are not met.
 *
 * Switchable implementation:
 *   - ALLLOW_ERROR : if we allow that there can be nonzero error, does not throw, and uses gain
 *   - ONLY_EXACT   : throws error at linearization if not at exact feasible point, and infinite error
 *
 * \nosubgrouping
 */
class NonlinearEquality: public NoiseModelFactor1
{

public:

  // feasible value
  minimatrix* feasible_;

  // error handling flag
  bool allow_error_;

  // error gain in allow error case
  double error_gain_;

  /**
   * Function that compares two values
   */
//  bool (*compare_)(const T& a, const T& b);
  double comparetol_;

public:
  /** default constructor - only for serialization */
  NonlinearEquality() {
  }

  virtual ~NonlinearEquality() {
      if(feasible_!=NULL)
      {
          delete feasible_;
          feasible_=NULL;
      }
      if(noiseModel_!=NULL)
      {
          delete noiseModel_;
          noiseModel_=NULL;
      }

  }

  /// @name Standard Constructors
  /// @{

  /**
   * Constructor - allows inexact evaluation
   */
  NonlinearEquality(int j, minimatrix* feasible, double error_gain=0.0,double comparetol=1.0e-9,bool allow_error=true):
      NoiseModelFactor1(ConstrainedNoiseModel::All(feasible->dimension),
          j), feasible_(feasible), allow_error_(allow_error), error_gain_(error_gain), //
      comparetol_(comparetol)
  {

  }

  /// @}

  /// @name Standard Interface
  /// @{

  /** actual error function calculation */
  virtual double error(const std::map<int,minimatrix*>& c) const {
   minimatrix* xj=c.at(this->key());

    minivector e = this->unwhitenedError(c);
    if (allow_error_ || !minimatrix_equal(xj, feasible_,comparetol_)) {
      return error_gain_ *miniblas_vector_ddot(e,e); //dot(e, e);
    } else {
      return 0.0;
    }
  }

  /** error function */

  virtual minivector  evaluateError(const minimatrix* X) const
  {
    const size_t nj = feasible_->dimension;
    if (allow_error_) {
      return minivector(X->LocalCoordinates(feasible_));
    }
    else if (minimatrix_equal(feasible_, X,comparetol_))
    {
      return minivector(nj, 0.0); // set error to infinity if not equal
    } else {
        throw std::invalid_argument(
            "Linearization point not feasible for this key");
      return minivector(nj, std::numeric_limits<double>::infinity()); // set error to infinity if not equal
    }
  }

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


  virtual minivector  evaluateError(const minimatrix* X,minimatrix &H) const
  {
    const size_t nj = feasible_->dimension;
    if (allow_error_) {
        minimatrix_resize(&H,nj,nj);
        minimatrix_set_identity(&H);
      return minivector(X->LocalCoordinates(feasible_));
    }
    else if (minimatrix_equal(feasible_, X,comparetol_))
    {
        minimatrix_resize(&H,nj,nj);
        minimatrix_set_identity(&H);
      return minivector(nj, 0.0); // set error to infinity if not equal
    } else {
        throw std::invalid_argument(
            "Linearization point not feasible for this key. ");
      return minivector(nj, std::numeric_limits<double>::infinity()); // set error to infinity if not equal
    }
  }





  // Linearize is over-written, because base linearization tries to whiten
  virtual RealGaussianFactor*  linearize(const std::map<int, minimatrix*>& x,int factorizaton=0) const
  {

     const minimatrix* xj = x.at(this->key());
     minimatrix A;
    minivector b = evaluateError(xj, A);
    GaussianNoiseModel* model = ConstrainedNoiseModel::All(b.size1);
    return new JacobianFactor(this->key(), A, b, model);

  }


  /// @return a deep copy of this factor
  virtual NoiseModelFactor* clone() const {
    return new NonlinearEquality(key(),feasible_,error_gain_,comparetol_,allow_error_);


  }

  /// @}

};
// \class NonlinearEquality


/* ************************************************************************* */
/**
 * Simple unary equality constraint - fixes a value for a variable
 */
class NonlinearEquality1: public NoiseModelFactor1
{

public:


  /** default constructor to allow for serialization */
  NonlinearEquality1() {
  }

  minimatrix* value_; /// fixed value for variable


public:

  /**
   * Constructor
   * @param value feasible value that values(key) shouild be equal to
   * @param key the key for the unknown variable to be constrained
   * @param mu a parameter which really turns this into a strong prior
   *
   */
  NonlinearEquality1(minimatrix* value, int key, double mu = 1000.0) :
      NoiseModelFactor1( ConstrainedNoiseModel::All(value->dimension,
              std::abs(mu)), key), value_(value)
  {
  }

    NonlinearEquality1(minimatrix* value, int key, GaussianNoiseModel* model) :
      NoiseModelFactor1(noiseModel_, key), value_(value)
  {
  }


  virtual ~NonlinearEquality1() {
  }

  /// @return a deep copy of this factor
  virtual NoiseModelFactor* clone() const {
    return new NonlinearEquality1(value_,key(),noiseModel_);
  }


   virtual minivector evaluateError(const minimatrix* X) const
  {
   return minivector(value_->LocalCoordinates(X));
  }


  virtual minivector  evaluateError(const minimatrix* X,minimatrix &H) const
  {
    int dim=X->dimension;
    minimatrix_resize(&H,dim,dim);
    minimatrix_set_identity(&H);
    // manifold equivalent of h(x)-z -> log(z,h(x))
   return minivector(value_->LocalCoordinates(X));

  }


};
// \NonlinearEquality1

/* ************************************************************************* */
/**
 * Simple binary equality constraint - this constraint forces two factors to
 * be the same.
 */
class NonlinearEquality2: public NoiseModelFactor2
{

public:


  /** default constructor to allow for serialization */
  NonlinearEquality2() {
  }

public:



  NonlinearEquality2(int key1, int key2, int dimension,double mu = 1000.0) :
      NoiseModelFactor2(ConstrainedNoiseModel::All(dimension, std::abs(mu)), key1, key2) {
  }

   NonlinearEquality2(int key1, int key2,GaussianNoiseModel* model) :
      NoiseModelFactor2(model, key1, key2) {
  }
  virtual ~NonlinearEquality2() {
  }

  /// @return a deep copy of this factor
  virtual NoiseModelFactor* clone() const {
    return new NonlinearEquality2(key1(),key2(),noiseModel_);
  }

  /** g(x) with optional derivative2 */
   virtual minivector
    evaluateError(const minimatrix* X1, const minimatrix* X2) const
{

    return minivector(X1->LocalCoordinates(X2));
  }
  virtual minivector
    evaluateError(const minimatrix* X1, const minimatrix* X2, minimatrix &H1, minimatrix &H2) const
 {
    size_t p = X1->dimension;
    minimatrix_resize(&H1,p,p);
    minimatrix_set_neg_identity(&H1);
    minimatrix_resize(&H2,p,p);
    minimatrix_set_identity(&H2);
    return minivector(X1->LocalCoordinates(X2));
  }
};
// \NonlinearEquality2
};// namespace minisam




