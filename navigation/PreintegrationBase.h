#ifndef PREINTEGRATIONBASE_H
#define PREINTEGRATIONBASE_H
/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 *  @file  PreintegrationBase.h
 *  @author Luca Carlone
 *  @author Stephen Williams
 *  @author Richard Roberts
 *  @author Vadim Indelman
 *  @author David Jensen
 *  @author Frank Dellaert
 **/


#pragma once

#include "../navigation/PreintegrationParams.h"
#include "../navigation/NavState.h"
#include "../navigation/ImuBias.h"
#include "../linear/NoiseModel.h"

#include <iosfwd>

namespace minisam
{

/**
 * PreintegrationBase is the base class for PreintegratedMeasurements
 * (in ImuFactor) and CombinedPreintegratedMeasurements (in CombinedImuFactor).
 * It includes the definitions of the preintegrated variables and the methods
 * to access, print, and compare them.
 */
class PreintegrationBase
{

protected:

    PreintegrationParams* p_;

    /// Acceleration and gyro bias used for preintegration
    ConstantBias biasHat_;

    /// Time interval from i to j
    double deltaTij_;

    /// Default constructor for serialization
    PreintegrationBase() {}

public:
    /// @name Constructors
    /// @{

    /**
     *  Constructor, initializes the variables in the base class
     *  @param p    Parameters, typically fixed in a single application
     *  @param bias Current estimate of acceleration and rotation rate biases
     */
    PreintegrationBase(PreintegrationParams* p,
                       const ConstantBias& biasHat = ConstantBias());

    PreintegrationBase(const PreintegrationBase& other);

    /// @}

    /// @name Basic utilities
    /// @{
    /// Re-initialize PreintegratedMeasurements
    virtual void resetIntegration()=0;

    /// Re-initialize PreintegratedMeasurements and set new bias
    void resetIntegrationAndSetBias(const ConstantBias& biasHat);

    /// check parameters equality: checks whether  pointer points to same Params object.
    bool matchesParamsWith(const PreintegrationBase& other) const;

    /// pointer to params
    PreintegrationParams* params() const;

    /// const reference to params
    PreintegrationParams* p() const;

    ///@}

    /// @name Instance variables access
    /// @{
    const ConstantBias& biasHat() const;

    void SetbiasHat(const ConstantBias& cb);

    double deltaTij() const;

    void SetdeltaTij(double tij);

    void setParam(PreintegrationParams* p);

    virtual Eigen::Vector3d  deltaPij() const=0;
    virtual Eigen::Vector3d  deltaVij() const=0;
    virtual Rot3     deltaRij() const=0;
    virtual NavState deltaXij() const=0;

    // Exposed for MATLABEigen::MatrixXd::Identity(3,3)
    Eigen::VectorXd biasHatVector() const;
    /// @}

    /// @name Main functionality
    /// @{

    /// Subtract estimate and correct for sensor pose
    /// Compute the derivatives due to non-identity body_P_sensor (rotation and centrifugal acc)
    /// Ignore D_correctedOmega_measuredAcc as it is trivially zero
    std::pair<Eigen::Vector3d, Eigen::Vector3d> correctMeasurementsBySensorPose(
        const Eigen::Vector3d& unbiasedAcc, const Eigen::Vector3d& unbiasedOmega,
        Eigen::Matrix3d* correctedAcc_H_unbiasedAcc = NULL,
        Eigen::Matrix3d* correctedAcc_H_unbiasedOmega = NULL,
        Eigen::Matrix3d* correctedOmega_H_unbiasedOmega = NULL) const;

    /// Update preintegrated measurements and get derivatives
    /// It takes measured quantities in the j frame
    /// Modifies preintegrated quantities in place after correcting for bias and possibly sensor pose
    virtual void update(const Eigen::Vector3d& measuredAcc, const Eigen::Vector3d& measuredOmega,
                        const double dt, Eigen::MatrixXd* A, Eigen::MatrixXd* B, Eigen::MatrixXd* C)=0;

    /// Version without derivatives
    virtual void integrateMeasurement(const Eigen::Vector3d& measuredAcc,
                                      const Eigen::Vector3d& measuredOmega, const double dt);

    /// Given the estimate of the bias, return a NavState tangent vector
    /// summarizing the preintegrated IMU measurements so far
    virtual Eigen::VectorXd biasCorrectedDelta(const Eigen::VectorXd& bias_i,
            Eigen::MatrixXd* H = NULL) const=0;

    /// Predict state at time j
    NavState predict(const NavState& state_i, const Eigen::VectorXd& bias_i,
                     Eigen::MatrixXd& H1,
                     Eigen::MatrixXd& H2) const;

    NavState predict(const NavState& state_i, const Eigen::VectorXd& bias_i) const;

    /// Calculate error given navStates
    Eigen::VectorXd computeError(const NavState& state_i, const NavState& state_j,
                                 const Eigen::VectorXd& bias_i,
                                 Eigen::MatrixXd& H1, Eigen::MatrixXd& H2,
                                 Eigen::MatrixXd& H3) const;
    Eigen::VectorXd computeError(const NavState& state_i, const NavState& state_j,
                                 const Eigen::VectorXd& bias_i) const;

    /// Compute errors w.r.t. preintegrated measurements and jacobians wrt pose_i, vel_i, bias_i, pose_j, bias_j
    Eigen::VectorXd computeErrorAndJacobians(const Pose3& pose_i, const Eigen::VectorXd& vel_i,
            const Pose3& pose_j, const Eigen::VectorXd& vel_j,
            const Eigen::VectorXd& bias_i, Eigen::MatrixXd& H1, Eigen::MatrixXd& H2,
            Eigen::MatrixXd& H3, Eigen::MatrixXd& H4, Eigen::MatrixXd& H5) const;

    Eigen::VectorXd computeErrorAndJacobians(const Pose3& pose_i, const Eigen::VectorXd& vel_i,
            const Pose3& pose_j, const Eigen::VectorXd& vel_j,
            const Eigen::VectorXd& bias_i) const;

};
};

#endif
