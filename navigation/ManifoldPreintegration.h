#ifndef MANIFOLDPREINTEGRATION_H
#define MANIFOLDPREINTEGRATION_H

/**
 *  @file  ManifoldPreintegration.h
 **/
#pragma once

#include "../navigation/PreintegrationParams.h"
#include "../navigation/ImuBias.h"
#include "../navigation/NavState.h"
#include "../mat/Matrix.h"
#include "../mat/MatCal.h"

namespace minisam
{

/**
 * IMU pre-integration on NavSatet manifold.
 * This corresponds to the original RSS paper (with one difference: V is rotated)
 */
class ManifoldPreintegration :public minimatrix// public PreintegrationBase
{
public:

    /**
     * Pre-integrated navigation state, from frame i to frame j
     * Note: relative position does not take into account velocity at time i, see deltap+, in [2]
     * Note: velocity is now also in frame i, as opposed to deltaVij in [2]
     */
    /// Default constructor for serialization
    ManifoldPreintegration():
#ifdef USE_QUATERNIONS
        minimatrix(35,3)
#else
        minimatrix(34,3)
#endif
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
    ManifoldPreintegration(const PreintegrationParams& p,
                           const ConstantBias& biasHat = ConstantBias());
    ManifoldPreintegration(const ManifoldPreintegration& other);

    ManifoldPreintegration(const minimatrix& other);
    /// @}

    /// @name Basic utilities
    /// @{
    /// Re-initialize PreintegratedMeasurements
    void resetIntegration();

    /// @}

    /// @name Instance variables access
    /// @{
    NavState deltaXij() const;

    Rot3   deltaRij() const
    {
        return deltaXij().attitude();
    }
    minivector  deltaPij() const
    {
#ifdef USE_QUATERNIONS
        minivector t_=minimatrix_row(*this,33);
#else
        minivector t_=minimatrix_row(*this,32);
#endif
        return t_;
    }
    minivector  deltaVij() const
    {
#ifdef USE_QUATERNIONS
        minivector v_=minimatrix_row(*this,34);
#else
        minivector v_=minimatrix_row(*this,33);
#endif
        return v_;
    }
    double deltaTij() const
    {
        return data[39];
    }
    PreintegrationParams p() const
    {
        minimatrix PreintegrationParams_=minimatrix_blockmatrix(*this,0,0,11,3);
        return PreintegrationParams(&PreintegrationParams_);
    }

    minimatrix  delRdelBiasOmega() const
    {
        minimatrix delRdelBiasOmega_=minimatrix_blockmatrix(*this,14,0,3,3);
        return delRdelBiasOmega_;
    }
    minimatrix  delPdelBiasAcc() const
    {
        minimatrix delPdelBiasAcc_=minimatrix_blockmatrix(*this,17,0,3,3);
        return delPdelBiasAcc_;
    }
    minimatrix  delPdelBiasOmega() const
    {
        minimatrix delPdelBiasOmega_=minimatrix_blockmatrix(*this,20,0,3,3);
        return delPdelBiasOmega_;
    }
    minimatrix  delVdelBiasAcc() const
    {
        minimatrix delVdelBiasAcc_=minimatrix_blockmatrix(*this,23,0,3,3);
        return delVdelBiasAcc_;
    }
    minimatrix  delVdelBiasOmega() const
    {
        minimatrix delVdelBiasOmega_=minimatrix_blockmatrix(*this,26,0,3,3);
        return delVdelBiasOmega_;
    }
    void  SetdelRdelBiasOmega(const minimatrix& n3d )
    {
        minimatrix delRdelBiasOmega_=minimatrix_blockmatrix_var(this,14,0,3,3);
        minimatrix_memcpy(&delRdelBiasOmega_,n3d);
    }
    void    SetdelPdelBiasAcc(const minimatrix& n3d )
    {
        minimatrix delPdelBiasAcc_=minimatrix_blockmatrix_var(this,17,0,3,3);
        minimatrix_memcpy(&delPdelBiasAcc_,n3d);
    }
    void  SetdelPdelBiasOmega(const minimatrix& n3d )
    {
        minimatrix delPdelBiasOmega_=minimatrix_blockmatrix_var(this,20,0,3,3);
        minimatrix_memcpy(&delPdelBiasOmega_,n3d);

    }
    void  SetdelVdelBiasAcc(const minimatrix& n3d )
    {
        minimatrix delVdelBiasAcc_=minimatrix_blockmatrix_var(this,23,0,3,3);
        minimatrix_memcpy(&delVdelBiasAcc_,n3d);
    }
    void   SetdelVdelBiasOmega(const minimatrix& n3d )
    {
        minimatrix delVdelBiasOmega_=minimatrix_blockmatrix_var(this,26,0,3,3);
        minimatrix_memcpy(&delVdelBiasOmega_,n3d);
    }
    /// @}

    /// @name Main functionality
    /// @{
    void integrateMeasurement(const minivector& measuredAcc,
                              const minivector& measuredOmega, const double dt);


    /// Predict state at time j
    NavState predict(const NavState& state_i, const ConstantBias& bias_i,
                     minimatrix& H1,
                     minimatrix& H2) const;


    NavState predict(const NavState& state_i, const ConstantBias& bias_i) const;

    /// Calculate error given navStates
    minivector computeError(const NavState& state_i, const NavState& state_j,
                            const ConstantBias& bias_i,
                            minimatrix& H1, minimatrix& H2,
                            minimatrix& H3) const;
    minivector computeError(const NavState& state_i, const NavState& state_j,
                            const ConstantBias& bias_i) const;

    /// Compute errors w.r.t. preintegrated measurements and jacobians wrt pose_i, vel_i, bias_i, pose_j, bias_j
    minivector computeErrorAndJacobians(const minimatrix* pose_i, const minimatrix* vel_i,
                                        const minimatrix* pose_j, const minimatrix* vel_j,
                                        const ConstantBias& bias_i, minimatrix& H1, minimatrix& H2,
                                        minimatrix& H3, minimatrix& H4, minimatrix& H5) const;

    minivector computeErrorAndJacobians(const minimatrix* pose_i, const minimatrix* vel_i,
                                        const minimatrix* pose_j, const minimatrix* vel_j,
                                        const ConstantBias& bias_i) const;
    /// Update preintegrated measurements and get derivatives
    /// It takes measured quantities in the j frame
    /// Modifies preintegrated quantities in place after correcting for bias and possibly sensor pose
    void update(const minivector& measuredAcc, const minivector& measuredOmega, const double dt,
                minimatrix* A, minimatrix* B, minimatrix* C) ;

    /// Given the estimate of the bias, return a NavState tangent vector
    /// summarizing the preintegrated IMU measurements so far
    minivector biasCorrectedDelta(const ConstantBias& bias_i,
                                  minimatrix* H = NULL) const ;

    ManifoldPreintegration* clone() const
    {
        ManifoldPreintegration* Mfp=new ManifoldPreintegration(*this);
        return Mfp;
    }

    /// @}
};
};
#endif
