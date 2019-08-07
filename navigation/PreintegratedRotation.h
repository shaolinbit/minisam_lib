#ifndef PREINTEGRATEDROTATION_H
#define PREINTEGRATEDROTATION_H

/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 *  @file  PreintegratedRotation.h
 *  @author Luca Carlone
 *  @author Stephen Williams
 *  @author Richard Roberts
 *  @author Vadim Indelman
 *  @author David Jensen
 *  @author Frank Dellaert
 **/

#pragma once

#include "../geometry/Pose3.h"
#include "../base/Matrix.h"
#include "../base/MatCal.h"

namespace minisam
{

/// Parameters for pre-integration:
/// Usage: Create just a single Params and pass a  pointer to the constructor
struct PreintegratedRotationParams
{
    Eigen::MatrixXd gyroscopeCovariance;  ///< continuous-time "Covariance" of gyroscope measurements
    Eigen::Vector3d omegaCoriolis;  ///< Coriolis constant
    Pose3* body_P_sensor;    ///< The pose of the sensor in the body frame

    PreintegratedRotationParams() : gyroscopeCovariance(Eigen::MatrixXd::Identity(3,3)),omegaCoriolis(Eigen::Vector3d::Zero(3)),body_P_sensor(NULL) {}

    virtual ~PreintegratedRotationParams() {}

    void setGyroscopeCovariance(const Eigen::MatrixXd& cov)
    {
        gyroscopeCovariance = cov;
    }
    void setOmegaCoriolis(const Eigen::VectorXd& omega)
    {
        omegaCoriolis=omega;
    }
    void setBodyPSensor(const Pose3& pose)
    {
        *body_P_sensor=pose;
    }

    const Eigen::MatrixXd& getGyroscopeCovariance()     const
    {
        return gyroscopeCovariance;
    }
    Eigen::Vector3d getOmegaCoriolis() const
    {
        return omegaCoriolis;
    }
    Pose3*   getBodyPSensor()   const
    {
        return body_P_sensor;
    }
};

/**
 * PreintegratedRotation is the base class for all PreintegratedMeasurements
 * classes (in AHRSFactor, ImuFactor, and CombinedImuFactor).
 * It includes the definitions of the preintegrated rotation.
 */
class PreintegratedRotation
{
protected:
    /// Parameters
    PreintegratedRotationParams* p_;

    double deltaTij_;           ///< Time interval from i to j
    Rot3 deltaRij_;             ///< Preintegrated relative orientation (in frame i)
    Eigen::MatrixXd delRdelBiasOmega_;  ///< Jacobian of preintegrated rotation w.r.t. angular rate bias

    /// Default constructor for serialization
    PreintegratedRotation() {}

public:
    /// @name Constructors
    /// @{

    /// Default constructor, resets integration to zero
    explicit PreintegratedRotation(PreintegratedRotationParams* p);

    /// Explicit initialization of all class members
    PreintegratedRotation(PreintegratedRotationParams* p,
                          double deltaTij, const Rot3& deltaRij,
                          const Eigen::MatrixXd& delRdelBiasOmega);
    /// @}

    /// @name Basic utilities
    /// @{

    /// Re-initialize PreintegratedMeasurements
    void resetIntegration();

    /// check parameters equality: checks whether  pointer points to same Params object.
    bool matchesParamsWith(const PreintegratedRotation& other) const;
    /// @}

    /// @name Access instance variables
    /// @{
    const PreintegratedRotationParams* params() const;
    const double& deltaTij() const;
    const Rot3& deltaRij() const;
    const Eigen::MatrixXd& delRdelBiasOmega() const;
    /// @}

    /// @name Main functionality
    /// @{

    /// Take the gyro measurement, correct it using the (constant) bias estimate
    /// and possibly the sensor pose, and then integrate it forward in time to yield
    /// an incremental rotation.
    Rot3 incrementalRotation(const Eigen::VectorXd& measuredOmega, const Eigen::VectorXd& biasHat, double deltaT,
                             Eigen::Matrix3d* D_incrR_integratedOmega) const;
    // OptionalJacobian<3, 3> D_incrR_integratedOmega) const;

    /// Calculate an incremental rotation given the gyro measurement and a time interval,
    /// and update both deltaTij_ and deltaRij_.
    void integrateMeasurement(const Eigen::VectorXd& measuredOmega, const Eigen::VectorXd& biasHat, double deltaT,
                              Eigen::MatrixXd* D_incrR_integratedOmega=NULL,
                              Eigen::Matrix3d*  F =NULL);


    /// Return a bias corrected version of the integrated rotation, with optional Jacobian
    Rot3 biascorrectedDeltaRij(const Eigen::VectorXd& biasOmegaIncr,
                               Eigen::MatrixXd& H);
    Rot3 biascorrectedDeltaRij(const Eigen::VectorXd& biasOmegaIncr);

    /// Integrate coriolis correction in body frame rot_i
    Eigen::VectorXd integrateCoriolis(const Rot3& rot_i) const;

    /// @}

};
};
#endif
