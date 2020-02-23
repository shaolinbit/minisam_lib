#ifndef PREINTEGRATEDROTATION_H
#define PREINTEGRATEDROTATION_H

/**
 *  @file  PreintegratedRotation.h
 **/

#pragma once

#include "../geometry/Pose3.h"
#include "../mat/Matrix.h"
#include "../mat/MatCal.h"

namespace minisam
{

/// Parameters for pre-integration:
/// Usage: Create just a single Params and pass a  pointer to the constructor
struct PreintegratedRotationParams:minimatrix
{

    PreintegratedRotationParams():minimatrix(4,3)
    {
        minimatrix_set_zero(this);
        data[0]=1.0;
        data[4]=1.0;
        data[8]=1.0;
    }

    virtual ~PreintegratedRotationParams()
    {

    }

    void setGyroscopeCovariance(const minimatrix& cov)
    {
        minimatrix gyroscopeCovariance=minimatrix_blockmatrix_var(this,0,0,3,3);
        minimatrix_memcpy(&gyroscopeCovariance,cov);
    }
    void setOmegaCoriolis(const minivector& omega)
    {
        minivector omegaCoriolis=minimatrix_row(this,3);
        minivector_memcpy(&omegaCoriolis,omega);
    }


     minimatrix getGyroscopeCovariance()     const
    {
        minimatrix gyroscopeCovariance=minimatrix_blockmatrix(*this,0,0,3,3);
        return gyroscopeCovariance;
    }
     minivector getOmegaCoriolis() const
    {
        minivector omegaCoriolis=minimatrix_row(*this,3);
        return omegaCoriolis;
    }
};

/**
 * PreintegratedRotation is the base class for all PreintegratedMeasurements
 * classes (in AHRSFactor, ImuFactor, and CombinedImuFactor).
 * It includes the definitions of the preintegrated rotation.
 */
class PreintegratedRotation//not Ready for Parallel computing.
{
protected:
    /// Parameters
    PreintegratedRotationParams* p_;

    double deltaTij_;           ///< Time interval from i to j
    Rot3 deltaRij_;             ///< Preintegrated relative orientation (in frame i)
    minimatrix delRdelBiasOmega_;  ///< Jacobian of preintegrated rotation w.r.t. angular rate bias

    /// Default constructor for serialization
    PreintegratedRotation() {}

    ~PreintegratedRotation()
    {
    }

public:
    /// @name Constructors
    /// @{

    /// Default constructor, resets integration to zero
    explicit PreintegratedRotation(PreintegratedRotationParams* p);

    /// Explicit initialization of all class members
    PreintegratedRotation(PreintegratedRotationParams* p,
                          double deltaTij, const Rot3& deltaRij,
                          const minimatrix& delRdelBiasOmega);
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
    const minimatrix& delRdelBiasOmega() const;
    /// @}

    /// @name Main functionality
    /// @{

    /// Take the gyro measurement, correct it using the (constant) bias estimate
    /// and possibly the sensor pose, and then integrate it forward in time to yield
    /// an incremental rotation.
    Rot3 incrementalRotation(const minivector& measuredOmega, const minivector& biasHat, double deltaT,
                             minimatrix* D_incrR_integratedOmega) const;

    /// Calculate an incremental rotation given the gyro measurement and a time interval,
    /// and update both deltaTij_ and deltaRij_.
    void integrateMeasurement(const minivector& measuredOmega, const minivector& biasHat, double deltaT,
                              minimatrix* D_incrR_integratedOmega=NULL,
                              minimatrix*  F =NULL);


    /// Return a bias corrected version of the integrated rotation, with optional Jacobian
    Rot3 biascorrectedDeltaRij(const minivector& biasOmegaIncr,
                               minimatrix& H);
    Rot3 biascorrectedDeltaRij(const minivector& biasOmegaIncr);

    /// Integrate coriolis correction in body frame rot_i
    minivector integrateCoriolis(const Rot3& rot_i) const;

    /// @}

};
};
#endif
