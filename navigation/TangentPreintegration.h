#ifndef TANGENTPREINTEGRATION_H
#define TANGENTPREINTEGRATION_H


/**
 *  @file   TangentPreintegration.h
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
 * Integrate on the 9D tangent space of the NavState manifold.
 * See extensive discussion in ImuFactor.lyx
 */
class TangentPreintegration : public minimatrix
{

    /**
     * Preintegrated navigation state, as a 9D vector on tangent space at frame i
     * Order is: theta, position, velocity
     */
    /// Acceleration and gyro bias used for preintegration
    /// Time interval from i to j
public:
    /// Default constructor for serialization
    TangentPreintegration():minimatrix(35,3)
    {
        resetIntegration();
    }

    ~TangentPreintegration()
    {
    }

public:
    /// @name Constructors
    /// @{

    /**
     *  Constructor, initializes the variables in the base class
     *  @param p    Parameters, typically fixed in a single application
     *  @param bias Current estimate of acceleration and rotation rate biases
     */
    TangentPreintegration(const PreintegrationParams&  p,
                          const ConstantBias& biasHat = ConstantBias());

    TangentPreintegration(const TangentPreintegration& other);

    TangentPreintegration(const minimatrix& other);

    TangentPreintegration(const minimatrix* other);

    /// @}

    /// @name Basic utilities
    /// @{
    /// Re-initialize PreintegratedMeasurements
    void resetIntegration();


    /// @}

    /// @name Instance variables access
    /// @{
    minivector theta() const
    {
        return minimatrix_row(*this,14);
    }

    minivector  deltaPij() const
    {
        return minimatrix_row(*this,15);

    }
    minivector   deltaVij() const
    {
        return minimatrix_row(*this,16);
    }
    Rot3   deltaRij() const
    {
        return Rot3::Expmap(theta());
    }
    NavState deltaXij() const
    {
        return NavState().retract(preintegrated());
    }

    minivector preintegrated() const
    {
        minivector preintegrated_(9);
        preintegrated_.data[0]=minimatrix_get(*this,14,0);
        preintegrated_.data[1]=minimatrix_get(*this,14,1);
        preintegrated_.data[2]=minimatrix_get(*this,14,2);
        preintegrated_.data[3]=minimatrix_get(*this,15,0);
        preintegrated_.data[4]=minimatrix_get(*this,15,1);
        preintegrated_.data[5]=minimatrix_get(*this,15,2);
        preintegrated_.data[6]=minimatrix_get(*this,16,0);
        preintegrated_.data[7]=minimatrix_get(*this,16,1);
        preintegrated_.data[8]=minimatrix_get(*this,16,2);

        return preintegrated_;
    }

    minimatrix preintegrated_H_biasAcc() const
    {
        minimatrix preintegrated_H_biasAcc_=minimatrix_blockmatrix(*this,17,0,9,3);
        return preintegrated_H_biasAcc_;
    }
    minimatrix preintegrated_H_biasOmega() const
    {
        minimatrix preintegrated_H_biasOmega_=minimatrix_blockmatrix(*this,26,0,9,3);
        return preintegrated_H_biasOmega_;
    }

    double deltaTij() const
    {
        return data[39];
    }
    const ConstantBias& biasHat() const
    {
        minimatrix mb=minimatrix_blockmatrix(*this,11,0,2,3);
        return ConstantBias(&mb);
    }

    const PreintegrationParams& p() const
    {
        minimatrix PreintegrationParams_=minimatrix_blockmatrix(*this,0,0,11,3);
        return PreintegrationParams(&PreintegrationParams_);
    }
    ///@}
    /// @name test
    /// @{
    void print();
    /// @}
    /// @name Main functionality
    /// @{
    void integrateMeasurement(const minivector& measuredAcc,
                              const minivector& measuredOmega, const double dt);

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
    // Update integrated vector on tangent manifold preintegrated with acceleration
    // Static, functional version.
    minivector UpdatePreintegrated(const minivector& a_body,
                                   const minivector& w_body, const double dt,
                                   const minivector& preintegrated,
                                   minimatrix* A = NULL,
                                   minimatrix* B = NULL,
                                   minimatrix* C = NULL);

    /// Update preintegrated measurements and get derivatives
    /// It takes measured quantities in the j frame
    /// Modifies preintegrated quantities in place after correcting for bias and possibly sensor pose
    void update(const minivector& measuredAcc, const minivector& measuredOmega,
                const double dt,minimatrix* A, minimatrix* B, minimatrix* C) ;

    /// Given the estimate of the bias, return a NavState tangent vector
    /// summarizing the preintegrated IMU measurements so far
    minivector biasCorrectedDelta(const ConstantBias& bias_i,
                                  minimatrix* H = NULL) const;

    // Compose the two pre-integrated 9D-vectors zeta01 and zeta02, with derivatives
    static minivector Compose(const minivector& zeta01, const minivector& zeta12,
                              double deltaT12,
                              minimatrix* H1 = NULL,
                              minimatrix* H2 = NULL);

    /// Merge in a different set of measurements and update bias derivatives accordingly
    /// The derivatives apply to the preintegrated Vector9
    void mergeWith(const TangentPreintegration& pim, minimatrix* H1, minimatrix* H2);
    /// @}

    TangentPreintegration* clone() const
    {
        TangentPreintegration* xb=new TangentPreintegration(*this);
        return xb;
    }
};
};

#endif
