#ifndef EXTENDEDKALMANFILTER_H
#define EXTENDEDKALMANFILTER_H

/**
 * @file    ExtendedKalmanFilter.h
 * @brief   Class to perform generic Kalman Filtering using nonlinear factor graphs
 */

// \callgraph
#pragma once

#include "../nonlinear/NonlinearFactorGraph.h"
#include "../nonlinear/NonlinearFactor.h"

namespace minisam
{

/**
 * This is a generic Extended Kalman Filter class implemented using nonlinear factors. GTSAM
 * basically does SRIF with Cholesky to solve the filter problem, making this an efficient,
 *numerically
 * stable Kalman Filter implementation.
 *
 * The Kalman Filter relies on two models: a motion model that predicts the next state using
 * the current state, and a measurement model that predicts the measurement value at a given
 * state. Because these two models are situation-dependent, base classes for each have been
 * provided above, from which the user must derive a class and incorporate the actual modeling
 * equations.
 *
 * The class provides a "predict" and "update" function to perform these steps independently.
 * TODO: a "predictAndUpdate" that combines both steps for some computational savings.
 * \nosubgrouping
 */


class ExtendedKalmanFilter
{


 public:

  minimatrix* x_;                                     // linearization point
  JacobianFactor* priorFactor_;  // Gaussian density on x_
  Factorization factorizatiotype_;
  JacobianFactor* postFactor_;





 public:
  /// @name Standard Constructors
  /// @{

  ExtendedKalmanFilter(int  key_initial, minimatrix* x_initial, GaussianNoiseModel* P_initial,Factorization ftype=CHOLESKY);

  ~ExtendedKalmanFilter()
  {
      if(priorFactor_->model_!=NULL)
  {
   delete priorFactor_->model_;
   priorFactor_->model_=NULL;
  }
  delete priorFactor_;

  if(postFactor_->model_!=NULL)
  {
   delete postFactor_->model_;
   postFactor_->model_=NULL;
  }
  delete postFactor_;

  delete x_;

  }

  minimatrix* solve_(const GaussianFactorGraph& linearFactorGraph, const std::map<int,minimatrix*>& linearizationPoints,
                  int x);
  /// @}

  /// @name Interface
  /// @{

  /**
   * Calculate predictive density P(x_) ~ \int  P(x_min) P(x_min, x_)
   * The motion model should be given as a factor with key1 for x_min and key2_ for x
   */
  void predict(NoiseModelFactor* motionFactor,
  GaussianFactorGraph& linearFactorGraph);
  /**
   * Calculate posterior density P(x_) ~ L(z|x) P(x)
   * The likelihood L(z|x) should be given as a unary factor on x
   */
  void  update(NoiseModelFactor* measurementFactor,
  GaussianFactorGraph& linearFactorGraph);


  /// Return current predictive (if called after predict)/posterior (if called after update)
   JacobianFactor* Density() const {
    return priorFactor_;
  }

  /// @}
};

};  // namespace

#endif // EXTENDEDKALMANFILTER_H
