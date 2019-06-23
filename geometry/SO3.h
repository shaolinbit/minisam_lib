#ifndef SO3_H
#define SO3_H


/**
 * @file    SO3.h
 * @brief   3*3 matrix representation of SO(3)
 * @author
 * @date
 */

#pragma once

#include "../base/Matrix.h"
#include "../base/MatCal.h"

#include <cmath>
#include <iosfwd>

namespace minisam
{

/**
 *  True SO(3), i.e., 3*3 matrix subgroup
 *  We guarantee (all but first) constructors only generate from sub-manifold.
 *  However, round-off errors in repeated composition could move off it...
 */
class SO3
{
public:
    Eigen::Matrix3d Matrix3;
public:
    enum
    {
        dimension = 3
    };

    /// @name Constructors
    /// @{

    /// Constructor from AngleAxisd
    SO3() :
        Matrix3(Eigen::MatrixXd::Identity(3,3))
    {
    }

    /// Constructor from Eigen Matrix
    SO3(const Eigen::MatrixXd& R) :
        Matrix3(R.eval())
    {
    }

    /// Constructor from AngleAxisd
    SO3(const Eigen::AngleAxisd& angleAxis) :
        Matrix3(angleAxis)
    {
    }

    /// Static, named constructor TODO think about relation with above
    static SO3 AxisAngle(const Eigen::VectorXd& axis, double theta);

    /// @}

    /// @name Group
    /// @{

    /// identity rotation for group operation
    static SO3 identity();

    /// inverse of a rotation = transpose
    SO3 inverse() const;

    /// @}
    /// @name Lie Group
    /// @{

    /**
     * Exponential map at identity - create a rotation from canonical coordinates
     * \f$ [R_x,R_y,R_z] \f$ using Rodrigues' formula
     */
    static SO3 Expmap(const Eigen::Vector3d& omega);
    static SO3 Expmap(const Eigen::Vector3d& omega,Eigen::MatrixXd* H);
    /// Derivative of Expmap
    static Eigen::Matrix3d ExpmapDerivative(const Eigen::Vector3d& omega);

    /**
     * Log map at identity - returns the canonical coordinates
     * \f$ [R_x,R_y,R_z] \f$ of this rotation
     */
    static Eigen::Vector3d Logmap(const SO3& R);//, ChartJacobian H = boost::none);
    static Eigen::Vector3d Logmap(const SO3& R, Eigen::Matrix3d* H);
    /// Derivative of Logmap
    static Eigen::Matrix3d LogmapDerivative(const Eigen::Vector3d& omega);

    Eigen::Matrix3d AdjointMap() const
    {
        return this->Matrix3;
    }
    Eigen::Vector3d operator*(const Eigen::Vector3d& p) const
    {
        return Matrix3*p;
    }

    Eigen::Matrix3d operator*(const  Eigen::Matrix3d& p) const
    {
        return Matrix3*p;
    }

    Eigen::Matrix3d operator*(double p) const
    {
        return Matrix3*p;
    }


    // Chart at origin
    struct ChartAtOrigin
    {
        static SO3 retract(const Eigen::VectorXd& omega)
        {
            return Expmap(omega);
        }
        static SO3 retract(const Eigen::VectorXd& omega, Eigen::MatrixXd* H)
        {
            return Expmap(omega, H);
        }
        static Eigen::Vector3d Local(const SO3& R)
        {
            return Logmap(R);
        }
        static Eigen::Vector3d Local(const SO3& R, Eigen::Matrix3d* H)
        {
            return Logmap(R, H);
        }
    };

    /// @}
};

// This namespace exposes two functors that allow for saving computation when
// exponential map and its derivatives are needed at the same location in so<3>
// The second functor also implements dedicated methods to apply dexp and/or inv(dexp)

/// Functor implementing Exponential map
class ExpmapFunctor
{
protected:
    const double theta2;
    Eigen::Matrix3d W, K, KK;
    bool nearZero;
    double theta, sin_theta, one_minus_cos;  // only defined if !nearZero

    void init(bool nearZeroApprox = false);

public:
    /// Constructor with element of Lie algebra so(3)
    ExpmapFunctor(const Eigen::VectorXd& omega, bool nearZeroApprox = false);

    /// Constructor with axis-angle
    ExpmapFunctor(const Eigen::VectorXd& axis, double angle, bool nearZeroApprox = false);

    /// Rodrigues formula
    SO3 expmap() const;
};

/// Functor that implements Exponential map *and* its derivatives
class DexpFunctor : public ExpmapFunctor
{
    const Eigen::Vector3d omega;
    double a, b;
    Eigen::Matrix3d dexp_;

public:
    /// Constructor with element of Lie algebra so(3)
    DexpFunctor(const Eigen::VectorXd& omega, bool nearZeroApprox = false);

    // NOTE(luca): Right Jacobian for Exponential map in SO(3) - equation
    // (10.86) and following equations in G.S. Chirikjian, "Stochastic Models,
    // Information Theory, and Lie Groups", Volume 2, 2008.
    //   expmap(omega + v) \approx expmap(omega) * expmap(dexp * v)
    // This maps a perturbation v in the tangent space to
    // a perturbation on the manifold Expmap(dexp * v) */
    const Eigen::Matrix3d& dexp() const
    {
        return dexp_;
    }

    /// Multiplies with dexp(), with optional derivatives
    Eigen::Vector3d applyDexp(const Eigen::Vector3d& v) const;
    Eigen::Vector3d applyDexp(const Eigen::Vector3d& v,Eigen::Matrix3d* H1) const;
    Eigen::Vector3d applyDexp(const Eigen::Vector3d& v,Eigen::Matrix3d* H1,
                              Eigen::Matrix3d* H2) const;

    /// Multiplies with dexp().inverse(), with optional derivatives
    Eigen::Vector3d applyInvDexp(const Eigen::Vector3d& v) const;
    Eigen::Vector3d applyInvDexp(const Eigen::Vector3d& v,
                                 Eigen::Matrix3d* H1,
                                 Eigen::Matrix3d* H2 ) const;
};
};

#endif // SO3_H
