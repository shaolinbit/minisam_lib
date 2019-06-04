#ifndef IMUBIAS_H
#define IMUBIAS_H

/**
 * @file ImuBias.h
 * @date  Feb 2, 2012
 * @author
 */

#pragma once
#include <iosfwd>
#include "../base/Matrix.h"

class ConstantBias
{
private:
    Eigen::Vector3d biasAcc_;
    Eigen::Vector3d biasGyro_;

public:
    /// dimension of the variable - used to autodetect sizes
    static const int dimension = 6;

    ConstantBias() :
        biasAcc_(0.0, 0.0, 0.0), biasGyro_(0.0, 0.0, 0.0)
    {
    }

    ConstantBias(const  Eigen::Vector3d& biasAcc, const  Eigen::Vector3d& biasGyro) :
        biasAcc_(biasAcc), biasGyro_(biasGyro)
    {
    }

    explicit ConstantBias(const  Eigen::VectorXd& v) :
        biasAcc_(v.head<3>()), biasGyro_(v.tail<3>())
    {
    }

    /** return the accelerometer and gyro biases in a single vector */
    Eigen::VectorXd vector() const
    {
        Eigen::VectorXd v(6);
        v << biasAcc_, biasGyro_;
        return v;
    }

    /** get accelerometer bias */
    const Eigen::Vector3d& accelerometer() const
    {
        return biasAcc_;
    }

    /** get gyroscope bias */
    const Eigen::Vector3d& gyroscope() const
    {
        return biasGyro_;
    }

    /** Correct an accelerometer measurement using this bias model, and optionally compute Jacobians */
    Eigen::Vector3d correctAccelerometer(const Eigen::Vector3d& measurement,
                                         Eigen::MatrixXd* H1=NULL,
                                         Eigen::MatrixXd* H2=NULL) const
    {
        if (H1!=NULL)
            (*H1) << - Eigen::MatrixXd::Identity(3,3), Eigen::MatrixXd::Zero(3,3);
        if (H2!=NULL)
            (*H2) << Eigen::MatrixXd::Identity(3,3);
        return measurement - biasAcc_;
    }

    /** Correct a gyroscope measurement using this bias model, and optionally compute Jacobians */
    Eigen::Vector3d correctGyroscope(const Eigen::Vector3d& measurement,
                                     Eigen::MatrixXd* H1=NULL,
                                     Eigen::MatrixXd* H2=NULL) const
    {
        if (H1!=NULL)
            (*H1) << Eigen::MatrixXd::Zero(3,3), -Eigen::MatrixXd::Identity(3,3);
        if (H2!=NULL)
            (*H2) << Eigen::MatrixXd::Identity(3,3);
        return measurement - biasGyro_;
    }

    /// @}
    /// @name Testable
    /// @{

    /// ostream operator
    friend std::ostream& operator<<(std::ostream& os,
                                    const ConstantBias& bias);

    /// print with optional string
    //void print(const std::string& s = "") const;

    /** equality up to tolerance
    inline bool equals(const ConstantBias& expected, double tol = 1e-5) const {
      return equal_with_abs_tol(biasAcc_, expected.biasAcc_, tol)
          && equal_with_abs_tol(biasGyro_, expected.biasGyro_, tol);
    }*/

    /// @}
    /// @name Group
    /// @{

    /** identity for group operation */
    static ConstantBias identity()
    {
        return ConstantBias();
    }

    /** inverse */
    inline ConstantBias operator-() const
    {
        return ConstantBias(-biasAcc_, -biasGyro_);
    }

    /** addition of vector on right */
    ConstantBias operator+(const Eigen::VectorXd& v) const
    {
        return ConstantBias(biasAcc_ + v.head<3>(), biasGyro_ + v.tail<3>());
    }

    /** addition */
    ConstantBias operator+(const ConstantBias& b) const
    {
        return ConstantBias(biasAcc_ + b.biasAcc_, biasGyro_ + b.biasGyro_);
    }

    /** subtraction */
    ConstantBias operator-(const ConstantBias& b) const
    {
        return ConstantBias(biasAcc_ - b.biasAcc_, biasGyro_ - b.biasGyro_);
    }

    /// @}

    /// @name Deprecated
    /// @{
    ConstantBias inverse()
    {
        return -(*this);
    }
    ConstantBias compose(const ConstantBias& q)
    {
        return (*this) + q;
    }
    ConstantBias between(const ConstantBias& q)
    {
        return q - (*this);
    }



    Eigen::VectorXd localCoordinates(const ConstantBias& q)
    {
        return between(q).vector();
    }
    ConstantBias retract(const Eigen::VectorXd& v)
    {
        return compose(ConstantBias(v));
    }
    static Eigen::VectorXd Logmap(const ConstantBias& p)
    {
        return p.vector();
    }
    static ConstantBias Expmap(const Eigen::VectorXd& v)
    {
        return ConstantBias(v);
    }
    /// @}


}; // ConstantBias class
ConstantBias ConstantBiasBetween(const ConstantBias& q1,const ConstantBias& q2,
                                 Eigen::MatrixXd* H1,Eigen::MatrixXd* H2);
#endif // IMUBIAS_H
