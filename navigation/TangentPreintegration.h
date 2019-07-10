#ifndef TANGENTPREINTEGRATION_H
#define TANGENTPREINTEGRATION_H

/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 *  @file   TangentPreintegration.h
 *  @author Frank Dellaert
 *  @author Adam Bry
 **/


#pragma once

#include "../navigation/PreintegrationBase.h"

namespace minisam
{

/**
 * Integrate on the 9D tangent space of the NavState manifold.
 * See extensive discussion in ImuFactor.lyx
 */
class TangentPreintegration : public PreintegrationBase
{
protected:

    /**
     * Preintegrated navigation state, as a 9D vector on tangent space at frame i
     * Order is: theta, position, velocity
     */
    Eigen::VectorXd preintegrated_;
    Eigen::MatrixXd preintegrated_H_biasAcc_;    ///< Jacobian of preintegrated_ w.r.t. acceleration bias
    Eigen::MatrixXd preintegrated_H_biasOmega_;  ///< Jacobian of preintegrated_ w.r.t. angular rate bias

    /// Default constructor for serialization
    TangentPreintegration()
    {
        resetIntegration();
    }

public:
    /// @name Constructors
    /// @{

    /**
     *  Constructor, initializes the variables in the base class
     *  @param p    Parameters, typically fixed in a single application
     *  @param bias Current estimate of acceleration and rotation rate biases
     */
    TangentPreintegration(PreintegrationParams*  p,
                          const ConstantBias& biasHat = ConstantBias());

    /// @}

    /// @name Basic utilities
    /// @{
    /// Re-initialize PreintegratedMeasurements
    void resetIntegration() override;

    /// @}

    /// @name Instance variables access
    /// @{
    Eigen::Vector3d  deltaPij() const override
    {
        return preintegrated_.segment<3>(3);
    }
    Eigen::Vector3d  deltaVij() const override
    {
        return preintegrated_.tail<3>();
    }
    Rot3     deltaRij() const override
    {
        return Rot3::Expmap(theta());
    }
    NavState deltaXij() const override
    {
        return NavState().retract(preintegrated_);
    }

    const Eigen::VectorXd& preintegrated() const
    {
        return preintegrated_;
    }
    Eigen::Vector3d theta() const
    {
        return preintegrated_.head<3>();
    }
    const Eigen::MatrixXd& preintegrated_H_biasAcc() const
    {
        return preintegrated_H_biasAcc_;
    }
    const Eigen::MatrixXd& preintegrated_H_biasOmega() const
    {
        return preintegrated_H_biasOmega_;
    }

    ///@}

    /// @name Main functionality
    /// @{

    // Update integrated vector on tangent manifold preintegrated with acceleration
    // Static, functional version.
    static Eigen::VectorXd UpdatePreintegrated(const Eigen::Vector3d& a_body,
            const Eigen::Vector3d& w_body, const double dt,
            const Eigen::VectorXd& preintegrated,
            Eigen::MatrixXd* A = NULL,
            Eigen::MatrixXd* B = NULL,
            Eigen::MatrixXd* C = NULL);

    /// Update preintegrated measurements and get derivatives
    /// It takes measured quantities in the j frame
    /// Modifies preintegrated quantities in place after correcting for bias and possibly sensor pose
    void update(const Eigen::Vector3d& measuredAcc, const Eigen::Vector3d& measuredOmega,
                const double dt,Eigen::MatrixXd* A, Eigen::MatrixXd* B, Eigen::MatrixXd* C) override;

    /// Given the estimate of the bias, return a NavState tangent vector
    /// summarizing the preintegrated IMU measurements so far
    Eigen::VectorXd biasCorrectedDelta(const Eigen::VectorXd& bias_i,
                                       Eigen::MatrixXd* H = NULL) const override;

    // Compose the two pre-integrated 9D-vectors zeta01 and zeta02, with derivatives
    static Eigen::VectorXd Compose(const Eigen::VectorXd& zeta01, const Eigen::VectorXd& zeta12,
                                   double deltaT12,
                                   Eigen::MatrixXd* H1 = NULL,
                                   Eigen::MatrixXd* H2 = NULL);

    /// Merge in a different set of measurements and update bias derivatives accordingly
    /// The derivatives apply to the preintegrated Vector9
    void mergeWith(const TangentPreintegration& pim, Eigen::MatrixXd* H1, Eigen::MatrixXd* H2);
    /// @}

    /** Dummy clone for MATLAB */
    virtual TangentPreintegration* clone() const
    {
        TangentPreintegration* xb=new TangentPreintegration();
        return xb;
    }
};
};

#endif
