#ifndef PREINTEGRATEDROTATION_H
#define PREINTEGRATEDROTATION_H

/**
 *  @file  PreintegratedRotation.h
 *  @author
 **/

#pragma once

#include "../geometry/Pose3.h"
#include "../base/Matrix.h"

/// Parameters for pre-integration:
/// Usage: Create just a single Params and pass a shared pointer to the constructor
struct PreintegratedRotationParams
{
    Eigen::MatrixXd gyroscopeCovariance;  ///< continuous-time "Covariance" of gyroscope measurements
    Eigen::Vector3d omegaCoriolis;  ///< Coriolis constant
    Pose3* body_P_sensor;    ///< The pose of the sensor in the body frame

    PreintegratedRotationParams() : gyroscopeCovariance(Eigen::MatrixXd::Identity(3,3)),omegaCoriolis(Eigen::Vector3d::Zero(3)),body_P_sensor(NULL) {}

    virtual ~PreintegratedRotationParams() {}

    //virtual void print(const std::string& s) const;
    //virtual bool equals(const PreintegratedRotationParams& other, double tol=1e-9) const;

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
    explicit PreintegratedRotation(PreintegratedRotationParams* p) : p_(p)
    {
        resetIntegration();
    }

    /// Explicit initialization of all class members
    PreintegratedRotation(PreintegratedRotationParams* p,
                          double deltaTij, const Rot3& deltaRij,
                          const Eigen::MatrixXd& delRdelBiasOmega)
        : p_(p), deltaTij_(deltaTij), deltaRij_(deltaRij), delRdelBiasOmega_(delRdelBiasOmega) {}

    /// @}

    /// @name Basic utilities
    /// @{

    /// Re-initialize PreintegratedMeasurements
    void resetIntegration();

    /// check parameters equality: checks whether shared pointer points to same Params object.
    bool matchesParamsWith(const PreintegratedRotation& other) const
    {
        return p_ == other.p_;
    }
    /// @}

    /// @name Access instance variables
    /// @{
    const PreintegratedRotationParams* params() const
    {
        return p_;
    }
    const double& deltaTij() const
    {
        return deltaTij_;
    }
    const Rot3& deltaRij() const
    {
        return deltaRij_;
    }
    const Eigen::MatrixXd& delRdelBiasOmega() const
    {
        return delRdelBiasOmega_;
    }
    /// @}

    /// @name Testable
    /// @{
    //void print(const std::string& s) const;
    //bool equals(const PreintegratedRotation& other, double tol) const;
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
    //  OptionalJacobian<3, 3> D_incrR_integratedOmega = boost::none,
    //  OptionalJacobian<3, 3> F = boost::none);

    /// Return a bias corrected version of the integrated rotation, with optional Jacobian
    Rot3 biascorrectedDeltaRij(const Eigen::VectorXd& biasOmegaIncr,
                               Eigen::MatrixXd& H);
    Rot3 biascorrectedDeltaRij(const Eigen::VectorXd& biasOmegaIncr);
    //OptionalJacobian<3, 3> H = boost::none) const;

    /// Integrate coriolis correction in body frame rot_i
    Eigen::VectorXd integrateCoriolis(const Rot3& rot_i) const;

    /// @}

};

#endif
