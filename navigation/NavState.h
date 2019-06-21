#ifndef NAVSTATE_H
#define NAVSTATE_H

/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 * @file    NavState.h
 * @brief   Navigation state composing of attitude, position, and velocity
 * @author  Frank Dellaert
 * @date    July 2015
 **/

#pragma once

#include "../geometry/Pose3.h"
/// Velocity is currently typedef'd to Vector3
namespace minisam
{

/**
 * Navigation state: Pose (rotation, translation) + velocity
 */
class NavState
{
private:

    Rot3 R_; ///< Rotation nRb, rotates points/velocities in body to points/velocities in nav
    Eigen::Vector3d t_; ///< position n_t, in nav frame
    Eigen::Vector3d v_; ///< velocity n_v in nav frame

public:

    enum
    {
        dimension = 9
    };

    typedef std::pair<Eigen::Vector3d, Eigen::Vector3d> PositionAndVelocity;

    /// @name Constructors
    /// @{

    /// Default constructor
    NavState() :
        t_(0, 0, 0), v_(Eigen::Vector3d::Zero())
    {
    }
    /// Construct from attitude, position, velocity
    NavState(const Rot3& R, const Eigen::Vector3d & t, const Eigen::Vector3d & v) :
        R_(R), t_(t), v_(v)
    {
    }
    /// Construct from pose and velocity
    NavState(const Pose3& pose, const Eigen::Vector3d & v) :
        R_(pose.rotation()), t_(pose.translation()), v_(v)
    {
    }
    /// Construct from SO(3) and R^6
    NavState(const Eigen::Matrix3d& R, const Eigen::VectorXd tv) :
        R_(R), t_(tv.head<3>()), v_(tv.tail<3>())
    {
    }
    /// Named constructor with derivatives
    static NavState Create(const Rot3& R, const Eigen::Vector3d& t, const Eigen::Vector3d & v,
                           Eigen::MatrixXd* H1, Eigen::MatrixXd* H2,
                           Eigen::MatrixXd* H3);
    /// Named constructor with derivatives
    static NavState FromPoseVelocity(const Pose3& pose, const Eigen::Vector3d & vel,
                                     Eigen::MatrixXd* H1, Eigen::MatrixXd* H2);

    /// @}
    /// @name Component Access
    /// @{

    const Rot3& attitude(Eigen::MatrixXd*  H = NULL) const;
    const Eigen::Vector3d & position(Eigen::MatrixXd* H = NULL) const;
    const Eigen::Vector3d & velocity(Eigen::MatrixXd* H = NULL) const;

    const Pose3 pose() const
    {
        return Pose3(attitude(), position());
    }

    /// @}
    /// @name Derived quantities
    /// @{

    /// Return rotation matrix. Induces computation in quaternion mode
    Eigen::Matrix3d R() const
    {
        return R_.matrix();
    }
    /// Return quaternion. Induces computation in matrix mode
    QQuaternion quaternion() const
    {
        return R_.toQuaternion();
    }
    /// Return position as Vector3
    Eigen::Vector3d  t() const
    {
        return t_;
    }
    /// Return velocity as Vector3. Computation-free.
    const Eigen::Vector3d & v() const
    {
        return v_;
    }
    // Return velocity in body frame
    Eigen::Vector3d  bodyVelocity(Eigen::MatrixXd* H = NULL) const;

    /// Return matrix group representation, in MATLAB notation:
    /// nTb = [nRb 0 n_t; 0 nRb n_v; 0 0 1]
    /// With this embedding in GL(3), matrix product agrees with compose
    Eigen::MatrixXd matrix() const;


    // Tangent space sugar.
    static Eigen::Block<Eigen::VectorXd, 3, 1> dR(Eigen::VectorXd& v)
    {
        return v.segment<3>(0);
    }
    static Eigen::Block<Eigen::VectorXd, 3, 1> dP(Eigen::VectorXd& v)
    {
        return v.segment<3>(3);
    }
    static Eigen::Block<Eigen::VectorXd, 3, 1> dV(Eigen::VectorXd& v)
    {
        return v.segment<3>(6);
    }
    static Eigen::Block<const Eigen::VectorXd, 3, 1> dR(const Eigen::VectorXd& v)
    {
        return v.segment<3>(0);
    }
    static Eigen::Block<const Eigen::VectorXd, 3, 1> dP(const Eigen::VectorXd& v)
    {
        return v.segment<3>(3);
    }
    static Eigen::Block<const Eigen::VectorXd, 3, 1> dV(const Eigen::VectorXd& v)
    {
        return v.segment<3>(6);
    }

    /// retract with optional derivatives
    NavState retract(const Eigen::VectorXd& v, //
                     Eigen::MatrixXd* H1 = NULL,Eigen::MatrixXd* H2 =
                         NULL) const;

    /// localCoordinates with optional derivatives
    Eigen::VectorXd localCoordinates(const NavState& g, //
                                     Eigen::MatrixXd* H1 = NULL, Eigen::MatrixXd* H2 =NULL) const;

    /// @}
    /// @name Dynamics
    /// @{

    /// Integrate forward in time given angular velocity and acceleration in body frame
    /// Uses second order integration for position, returns derivatives except dt.
    NavState update(const Eigen::Vector3d& b_acceleration, const Eigen::Vector3d& b_omega,
                    const double dt, Eigen::MatrixXd* F, Eigen::MatrixXd* G1,
                    Eigen::MatrixXd* G2) const;

    /// Compute tangent space contribution due to Coriolis forces
    Eigen::VectorXd coriolis(double dt, const Eigen::Vector3d& omega, bool secondOrder = false,
                             Eigen::MatrixXd* H = NULL) const;

    /// Correct preintegrated tangent vector with our velocity and rotated gravity,
    /// taking into account Coriolis forces if omegaCoriolis is given.
    Eigen::VectorXd correctPIM(const  Eigen::VectorXd& pim, double dt, const Eigen::Vector3d& n_gravity,
                               Eigen::Vector3d omegaCoriolis, bool use2ndOrderCoriolis =
                                   false, Eigen::MatrixXd* H1 = NULL,
                               Eigen::MatrixXd* H2 =NULL) const;

};
};

#endif
