#ifndef KALMANFILTER_H
#define KALMANFILTER_H


/**
 * @file KalmanFilter.h
 * @brief Simple linear Kalman filter. Implemented using factor graphs, i.e., does Cholesky-based SRIF, really.
 */


#include "../linear/GaussianDensity.h"
#include "../linear/GaussianFactorGraph.h"
#include "../linear/NoiseModel.h"
#include "../inference/Ordering.h"


namespace minisam {

/**
 * Kalman Filter class
 *
 * Knows how to maintain a Gaussian density under linear-Gaussian motion and
 * measurement models. It uses the square-root information form, as usual in GTSAM.
 *
 * The filter is functional, in that it does not have state: you call init() to create
 * an initial state, then predict() and update() that create new states out of an old state.
 */
class  KalmanFilter
{

public:

    /**
     *  This Kalman filter is a Square-root Information filter
     *  The type below allows you to specify the factorization variant.
     */


private:

    int n_; /** dimensionality of state */
    minimatrix I_; /** identity matrix of size n*n */
    Factorization factorizatiotype_;

    GaussianDensity* solve(const GaussianFactorGraph& factorGraph) const;
    GaussianDensity* fuse( GaussianDensity* p, RealGaussianFactor* newFactor) const;

public:

    // Constructor
    KalmanFilter(int n, Factorization ftype=CHOLESKY) :
        n_(n), I_(minimatrix(n,n)), factorizatiotype_(ftype)
    {
        minimatrix_set_identity(&I_);
    }

    /**
     * Create initial state, i.e., prior density at time k=0
     * In Kalman Filter notation, these are x_{0|0} and P_{0|0}
     * @param x0 estimate at time 0
     * @param P0 covariance at time 0, given as a diagonal Gaussian 'model'
     */
    GaussianDensity* init(const minivector& x0, GaussianNoiseModel* P0) const;

    /// version of init with a full covariance matrix
    GaussianDensity* init(const minivector& x0, const minimatrix& P0) const;


    /** Return step index k, starts at 0, incremented at each predict. */
    static int step(const GaussianDensity* p)
    {
        return p->keys().front();
    }

    /**
     * Predict the state P(x_{t+1}|Z^t)
     *   In Kalman Filter notation, this is x_{t+1|t} and P_{t+1|t}
     * Details and parameters:
     *   In a linear Kalman Filter, the motion model is f(x_{t}) = F*x_{t} + B*u_{t} + w
     *   where F is the state transition model/matrix, B is the control input model,
     *   and w is zero-mean, Gaussian white noise with covariance Q.
     */
    GaussianDensity* predict( GaussianDensity* p, const minimatrix& F, const minimatrix& B,
                            const minivector& u, GaussianNoiseModel* modelQ) const;


    /**
     * Predict the state P(x_{t+1}|Z^t)
     *   In Kalman Filter notation, this is x_{t+1|t} and P_{t+1|t}
     * Details and parameters:
     *   In a linear Kalman Filter, the motion model is f(x_{t}) = F*x_{t}+ w
     *   where F is the state transition model/matrix,
     *   and w is zero-mean, Gaussian white noise with covariance Q.
     */
    GaussianDensity* predictNoControl( GaussianDensity* p, const minimatrix& F,
                                        GaussianNoiseModel* modelQ) const;

    /*
     *  Version of predict with full covariance
     *  Q is normally derived as G*w*G^T where w models uncertainty of some
     *  physical property, such as velocity or acceleration, and G is derived from physics.
     *  This version allows more realistic models than a diagonal covariance matrix.
     */
    GaussianDensity* predictQ( GaussianDensity* p, const minimatrix& F, const minimatrix& B,
                             const minivector& u, const minimatrix& Q) const;

    /**
     * Predict the state P(x_{t+1}|Z^t)
     *   In Kalman Filter notation, this is x_{t+1|t} and P_{t+1|t}
     *   After the call, that is the density that can be queried.
     * Details and parameters:
     *   This version of predict takes GaussianFactor motion model [A0 A1 b]
     *   with an optional noise model.
     */
    GaussianDensity* predict2( GaussianDensity* p, const minimatrix& A0, const minimatrix& A1,
                             const minivector& b,  GaussianNoiseModel* model) const;

    /**
     * Update Kalman filter with a measurement
     * For the Kalman Filter, the measurement function, h(x_{t}) = z_{t}
     * will be of the form h(x_{t}) = H*x_{t} + v
     * where H is the observation model/matrix, and v is zero-mean,
     * Gaussian white noise with covariance R.
     * In this version, R is restricted to diagonal Gaussians (model parameter)
     */
    GaussianDensity* update( GaussianDensity* p, const minimatrix& H, const minivector& z,
                            GaussianNoiseModel* model) const;

    /*
     *  Version of update with full covariance
     *  Q is normally derived as G*w*G^T where w models uncertainty of some
     *  physical property, such as velocity or acceleration, and G is derived from physics.
     *  This version allows more realistic models than a diagonal covariance matrix.
     */
    GaussianDensity* updateQ( GaussianDensity* p, const minimatrix& H, const minivector& z,
                            const minimatrix& Q) const;
};

};

#endif // KALMANFILTER_H
