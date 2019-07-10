#ifndef MANIFOLDPREINTEGRATION_H
#define MANIFOLDPREINTEGRATION_H

/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 *  @file  ManifoldPreintegration.h
 *  @author Luca Carlone
 *  @author Stephen Williams
 *  @author Richard Roberts
 *  @author Vadim Indelman
 *  @author David Jensen
 *  @author Frank Dellaert
 **/
#pragma once

#include "../navigation/NavState.h"
#include "../navigation/PreintegrationBase.h"

namespace minisam
{

/**
 * IMU pre-integration on NavSatet manifold.
 * This corresponds to the original RSS paper (with one difference: V is rotated)
 */
class ManifoldPreintegration : public PreintegrationBase
{
protected:

    /**
     * Pre-integrated navigation state, from frame i to frame j
     * Note: relative position does not take into account velocity at time i, see deltap+, in [2]
     * Note: velocity is now also in frame i, as opposed to deltaVij in [2]
     */
    NavState deltaXij_;
    Eigen::Matrix3d delRdelBiasOmega_; ///< Jacobian of preintegrated rotation w.r.t. angular rate bias
    Eigen::Matrix3d delPdelBiasAcc_;   ///< Jacobian of preintegrated position w.r.t. acceleration bias
    Eigen::Matrix3d delPdelBiasOmega_; ///< Jacobian of preintegrated position w.r.t. angular rate bias
    Eigen::Matrix3d delVdelBiasAcc_;   ///< Jacobian of preintegrated velocity w.r.t. acceleration bias
    Eigen::Matrix3d delVdelBiasOmega_; ///< Jacobian of preintegrated velocity w.r.t. angular rate bias

    /// Default constructor for serialization
    ManifoldPreintegration()
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
    ManifoldPreintegration( PreintegrationParams* p,
                            const ConstantBias& biasHat = ConstantBias());

    /// @}

    /// @name Basic utilities
    /// @{
    /// Re-initialize PreintegratedMeasurements
    void resetIntegration() override;

    /// @}

    /// @name Instance variables access
    /// @{
    NavState deltaXij() const override
    {
        return deltaXij_;
    }
    Rot3     deltaRij() const override
    {
        return deltaXij_.attitude();
    }
    Eigen::Vector3d  deltaPij() const override
    {
        return deltaXij_.position();
    }
    Eigen::Vector3d  deltaVij() const override
    {
        return deltaXij_.velocity();
    }

    Eigen::Matrix3d  delRdelBiasOmega() const
    {
        return delRdelBiasOmega_;
    }
    Eigen::Matrix3d  delPdelBiasAcc() const
    {
        return delPdelBiasAcc_;
    }
    Eigen::Matrix3d  delPdelBiasOmega() const
    {
        return delPdelBiasOmega_;
    }
    Eigen::Matrix3d  delVdelBiasAcc() const
    {
        return delVdelBiasAcc_;
    }
    Eigen::Matrix3d  delVdelBiasOmega() const
    {
        return delVdelBiasOmega_;
    }
    void  SetdelRdelBiasOmega(const Eigen::Matrix3d& n3d )
    {
        delRdelBiasOmega_=n3d;
    }
    void    SetdelPdelBiasAcc(const Eigen::Matrix3d& n3d )
    {
        delPdelBiasAcc_=n3d;
    }
    void  SetdelPdelBiasOmega(const Eigen::Matrix3d& n3d )
    {
        delPdelBiasOmega_=n3d;
    }
    void  SetdelVdelBiasAcc(const Eigen::Matrix3d& n3d )
    {
        delVdelBiasAcc_=n3d;
    }
    void   SetdelVdelBiasOmega(const Eigen::Matrix3d& n3d )
    {
        delVdelBiasOmega_=n3d;
    }
    /// @}

    /// @name Main functionality
    /// @{

    /// Update preintegrated measurements and get derivatives
    /// It takes measured quantities in the j frame
    /// Modifies preintegrated quantities in place after correcting for bias and possibly sensor pose
    void update(const Eigen::Vector3d& measuredAcc, const Eigen::Vector3d& measuredOmega, const double dt,
                Eigen::MatrixXd* A, Eigen::MatrixXd* B, Eigen::MatrixXd* C) override;

    /// Given the estimate of the bias, return a NavState tangent vector
    /// summarizing the preintegrated IMU measurements so far
    Eigen::VectorXd biasCorrectedDelta(const Eigen::VectorXd& bias_i,
                                       Eigen::MatrixXd* H = NULL) const override;

    /** Dummy clone for MATLAB */
    virtual ManifoldPreintegration* clone() const
    {
        ManifoldPreintegration* Mfp=new ManifoldPreintegration(p(),biasHat());
        Mfp->SetdelRdelBiasOmega(delRdelBiasOmega());
        Mfp->SetdelPdelBiasAcc(delPdelBiasAcc());
        Mfp->SetdelPdelBiasOmega(delPdelBiasOmega());
        Mfp->SetdelVdelBiasAcc(delVdelBiasAcc());
        Mfp->SetdelVdelBiasOmega(delVdelBiasOmega());
        return Mfp;
    }

    /// @}
};
};
#endif
