#ifndef PREINTEGRATIONPARAMS_H
#define PREINTEGRATIONPARAMS_H

/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 *  @file  PreintegrationParams.h
 *  @author Frank Dellaert
 **/

#pragma once

#include "../navigation/PreintegratedRotation.h"
namespace minisam
{

/// Parameters for pre-integration:
/// Usage: Create just a single Params and pass a  pointer to the constructor
struct PreintegrationParams: PreintegratedRotationParams
{
    Eigen::Matrix3d accelerometerCovariance; ///< continuous-time "Covariance" of accelerometer
    Eigen::Matrix3d integrationCovariance; ///< continuous-time "Covariance" describing integration uncertainty
    bool use2ndOrderCoriolis; ///< Whether Eigen::MatrixXd::Identity(3,3)to use second order Coriolis integration
    Eigen::Vector3d n_gravity; ///< Gravity vector in nav frame

    /// The Params constructor insists on getting the navigation frame gravity vector
    /// For convenience, two commonly used conventions are provided by named constructors below
     PreintegrationParams(const Eigen::Vector3d& n_gravity)
        :PreintegratedRotationParams(PreintegratedRotationParams()), accelerometerCovariance(Eigen::MatrixXd::Identity(3,3)),
         integrationCovariance(Eigen::MatrixXd::Identity(3,3)),
         use2ndOrderCoriolis(false),
         n_gravity(n_gravity) {}

    // Default Params for a Z-down navigation frame, such as NED: gravity points along positive Z-axis
    static PreintegrationParams MakeSharedD(double g = 9.81)
    {
        return PreintegrationParams(Eigen::Vector3d(0, 0, g));
    }

    // Default Params for a Z-up navigation frame, such as ENU: gravity points along negative Z-axis
    static PreintegrationParams MakeSharedU(double g = 9.81)
    {
        return PreintegrationParams(Eigen::Vector3d(0, 0, -g));
    }



     void setAccelerometerCovariance(const Eigen::Matrix3d& cov)
    {
        accelerometerCovariance = cov;
    }
     void setIntegrationCovariance(const Eigen::Matrix3d& cov)
    {
        integrationCovariance = cov;
    }
     void setUse2ndOrderCoriolis(bool flag)
    {
        use2ndOrderCoriolis = flag;
    }

     const Eigen::Matrix3d& getAccelerometerCovariance() const
    {
        return accelerometerCovariance;
    }
    const Eigen::Matrix3d& getIntegrationCovariance()   const
    {
        return integrationCovariance;
    }
     bool           getUse2ndOrderCoriolis()     const
    {
        return use2ndOrderCoriolis;
    }

    bool operator==(const PreintegrationParams other)const
    {
        bool equp=(this->accelerometerCovariance==other.accelerometerCovariance)&&(this->integrationCovariance==other.integrationCovariance)
                  &&(this->use2ndOrderCoriolis&&other.use2ndOrderCoriolis)&&(this->n_gravity==other.n_gravity);
        return equp;
    }

protected:
    /// Default constructor for serialization only: uninitialized!
    PreintegrationParams() {}
};
};
#endif
