#ifndef IMUBIAS_H
#define IMUBIAS_H

/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 * @file ImuBias.h
 * @date  Feb 2, 2012
 * @author Vadim Indelman, Stephen Williams
 */

#pragma once
#include <iosfwd>
#include "../base/Matrix.h"
#include "../base/MatCal.h"
namespace minisam
{

class ConstantBias
{
private:
    Eigen::Vector3d biasAcc_;
    Eigen::Vector3d biasGyro_;

public:
    /// dimension of the variable - used to autodetect sizes
    static const int dimension = 6;

    ConstantBias();

    ConstantBias(const  Eigen::Vector3d& biasAcc, const  Eigen::Vector3d& biasGyro);

    explicit ConstantBias(const  Eigen::VectorXd& v);

    /** return the accelerometer and gyro biases in a single vector */
    Eigen::VectorXd vector() const;

    /** get accelerometer bias */
    const Eigen::Vector3d& accelerometer() const;

    /** get gyroscope bias */
    const Eigen::Vector3d& gyroscope() const;
    /** Correct an accelerometer measurement using this bias model, and optionally compute Jacobians */
    Eigen::Vector3d correctAccelerometer(const Eigen::Vector3d& measurement,
                                         Eigen::MatrixXd* H1=NULL,
                                         Eigen::MatrixXd* H2=NULL) const;

    /** Correct a gyroscope measurement using this bias model, and optionally compute Jacobians */
    Eigen::Vector3d correctGyroscope(const Eigen::Vector3d& measurement,
                                     Eigen::MatrixXd* H1=NULL,
                                     Eigen::MatrixXd* H2=NULL) const;
    /// @}
    /// @name Testable
    /// @{

    /// ostream operator
    friend std::ostream& operator<<(std::ostream& os,
                                    const ConstantBias& bias);

    /// @}
    /// @name Group
    /// @{

    /** identity for group operation */
    static ConstantBias identity();

    /** inverse */
    ConstantBias operator-() const;

    /** addition of vector on right */
    ConstantBias operator+(const Eigen::VectorXd& v) const;

    /** addition */
    ConstantBias operator+(const ConstantBias& b) const;

    /** subtraction */
    ConstantBias operator-(const ConstantBias& b) const;

    /**equal**/
     ConstantBias operator=(const ConstantBias& b) const;

    /// @}

    /// @name Deprecated
    /// @{
    ConstantBias inverse();
    ConstantBias compose(const ConstantBias& q);
    ConstantBias between(const ConstantBias& q);



    Eigen::VectorXd localCoordinates(const ConstantBias& q);
    ConstantBias retract(const Eigen::VectorXd& v);
    static Eigen::VectorXd Logmap(const ConstantBias& p);
    static ConstantBias Expmap(const Eigen::VectorXd& v);
    /// @}


}; // ConstantBias class
ConstantBias ConstantBiasBetween(const ConstantBias& q1,const ConstantBias& q2,
                                 Eigen::MatrixXd* H1,Eigen::MatrixXd* H2);

};
#endif // IMUBIAS_H
