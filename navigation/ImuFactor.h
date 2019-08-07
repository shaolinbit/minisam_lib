#ifndef IMUFACTOR_H
#define IMUFACTOR_H

/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 *  @file  ImuFactor.h
 *  @author Luca Carlone
 *  @author Stephen Williams
 *  @author Richard Roberts
 *  @author Vadim Indelman
 *  @author David Jensen
 *  @author Frank Dellaert
 **/

#pragma once

#include "../nonlinear/NonlinearFactor.h"
#include "../navigation/ManifoldPreintegration.h"
#include "../navigation/TangentPreintegration.h"


namespace minisam
{

#ifdef TANGENT_PREINTEGRATION
typedef TangentPreintegration PreintegrationType;
#else
typedef ManifoldPreintegration PreintegrationType;
#endif

/*
 * If you are using the factor, please cite:
 * L. Carlone, Z. Kira, C. Beall, V. Indelman, F. Dellaert, "Eliminating
 * conditionally independent sets in factor graphs: a unifying perspective based
 * on smart factors", Int. Conf. on Robotics and Automation (ICRA), 2014.
 *
 * C. Forster, L. Carlone, F. Dellaert, D. Scaramuzza, "IMU Preintegration on
 * Manifold for Efficient Visual-Inertial Maximum-a-Posteriori Estimation",
 * Robotics: Science and Systems (RSS), 2015.
 *
 * REFERENCES:
 * [1] G.S. Chirikjian, "Stochastic Models, Information Theory, and Lie Groups",
 *     Volume 2, 2008.
 * [2] T. Lupton and S.Sukkarieh, "Visual-Inertial-Aided Navigation for
 *     High-Dynamic Motion in Built Environments Without Initial Conditions",
 *     TRO, 28(1):61-76, 2012.
 * [3] L. Carlone, S. Williams, R. Roberts, "Preintegrated IMU factor:
 *     Computation of the Jacobian Matrices", Tech. Report, 2013.
 * [4] C. Forster, L. Carlone, F. Dellaert, D. Scaramuzza, "IMU Preintegration on
 *     Manifold for Efficient Visual-Inertial Maximum-a-Posteriori Estimation",
 *     Robotics: Science and Systems (RSS), 2015.
 */

/**
 * PreintegratedIMUMeasurements accumulates (integrates) the IMU measurements
 * (rotation rates and accelerations) and the corresponding covariance matrix.
 * The measurements are then used to build the Preintegrated IMU factor.
 * Integration is done incrementally (ideally, one integrates the measurement
 * as soon as it is received from the IMU) so as to avoid costly integration
 * at time of factor construction.
 *
 * @addtogroup SLAM
 */
class PreintegratedImuMeasurements: public PreintegrationType
{

    friend class ImuFactor;
    friend class ImuFactor2;

protected:

    Eigen::MatrixXd preintMeasCov_; ///< COVARIANCE OF: [PreintPOSITION PreintVELOCITY PreintROTATION]
    ///< (first-order propagation from *measurementCovariance*).
    PreintegrationParams* pp_;
    ConstantBias ConstantBias_;
public:

    /// Default constructor for serialization and Cython wrapper
    PreintegratedImuMeasurements();

    /**
      *  Constructor, initializes the class with no measurements
      *  @param bias Current estimate of acceleration and rotation rate biases
      *  @param p    Parameters, typically fixed in a single application
      */
    PreintegratedImuMeasurements(PreintegrationParams* p,
                                 const ConstantBias& biasHat = ConstantBias());

    /**
      *  Construct preintegrated directly from members: base class and preintMeasCov
      *  @param base               PreintegrationType instance
      *  @param preintMeasCov      Covariance matrix used in noise model.
      */
    PreintegratedImuMeasurements(const PreintegrationType& base, const Eigen::MatrixXd& preintMeasCov);

    /// Re-initialize PreintegratedIMUMeasurements
    void resetIntegration() override;

    /**
     * Add a single IMU measurement to the preintegration.
     * @param measuredAcc Measured acceleration (in body frame, as given by the sensor)
     * @param measuredOmega Measured angular velocity (as given by the sensor)
     * @param dt Time interval between this and the last IMU measurement
     */
    void integrateMeasurement(const Eigen::Vector3d& measuredAcc,
                              const Eigen::Vector3d& measuredOmega, const double dt) override;

    /// Add multiple measurements, in matrix columns
    void integrateMeasurements(const Eigen::MatrixXd& measuredAccs, const Eigen::MatrixXd& measuredOmegas,
                               const Eigen::MatrixXd& dts);

    /// Return pre-integrated measurement covariance
    Eigen::MatrixXd preintMeasCov() const;

#ifdef TANGENT_PREINTEGRATION
    /// Merge in a different set of measurements and update bias derivatives accordingly
    void mergeWith(const PreintegratedImuMeasurements& pim, Eigen::MatrixXd* H1, Eigen::MatrixXd* H2);
#endif
};

/**
 * ImuFactor is a 5-ways factor involving previous state (pose and velocity of
 * the vehicle at previous time step), current state (pose and velocity at
 * current time step), and the bias estimate. Following the preintegration
 * scheme proposed in [2], the ImuFactor includes many IMU measurements, which
 * are "summarized" using the PreintegratedIMUMeasurements class.
 * Note that this factor does not model "temporal consistency" of the biases
 * (which are usually slowly varying quantities), which is up to the caller.
 * See also CombinedImuFactor for a class that does this for you.
 *
 * @addtogroup SLAM
 */
class ImuFactor: public NoiseModelFactor
{
private:
    PreintegratedImuMeasurements _PIM_;
public:

    /** Default constructor - only use for serialization */
    ImuFactor();
    /**
     * Constructor
     * @param pose_i Previous pose key
     * @param vel_i  Previous velocity key
     * @param pose_j Current pose key
     * @param vel_j  Current velocity key
     * @param bias   Previous bias key
     */
    ImuFactor(int pose_i, int vel_i, int pose_j, int vel_j, int bias,
              const PreintegratedImuMeasurements& preintegratedMeasurements);

    virtual ~ImuFactor()
    {
    }

    inline int key1() const
    {
        return keys_[0];
    }
    inline int key2() const
    {
        return keys_[1];
    }
    inline int key3() const
    {
        return keys_[2];
    }
    inline int key4() const
    {
        return keys_[3];
    }
    inline int key5() const
    {
        return keys_[4];
    }

    /// @return a deep copy of this factor
    virtual NonlinearFactor* clone() const;

    /** Access the preintegrated measurements. */

    const PreintegratedImuMeasurements& preintegratedMeasurements() const
    {
        return _PIM_;
    }

    virtual Eigen::VectorXd unwhitenedError(const std::map<int, Pose3>& x1,const std::map<int, Eigen::VectorXd>& x2, std::vector<Eigen::MatrixXd> &H) const;

    virtual Eigen::VectorXd unwhitenedError(const std::map<int, Pose3>& x1,const std::map<int, Eigen::VectorXd>& x2) const;
    /** implement functions needed to derive from Factor */

    /// vector of errors
    Eigen::VectorXd evaluateError(const Pose3& pose_i, const Eigen::VectorXd& vel_i,
                                  const Pose3& pose_j, const Eigen::VectorXd& vel_j,
                                  const Eigen::VectorXd& bias_i, Eigen::MatrixXd& H1, Eigen::MatrixXd& H2,
                                  Eigen::MatrixXd& H3, Eigen::MatrixXd& H4, Eigen::MatrixXd& H5) const;
    Eigen::VectorXd evaluateError(const Pose3& pose_i, const Eigen::VectorXd& vel_i,
                                  const Pose3& pose_j, const Eigen::VectorXd& vel_j,
                                  const Eigen::VectorXd& bias_i) const;

#ifdef TANGENT_PREINTEGRATION
    /// Merge two pre-integrated measurement classes
    static PreintegratedImuMeasurements Merge(
        const PreintegratedImuMeasurements& pim01,
        const PreintegratedImuMeasurements& pim12);

    /// Merge two factors
    static ImuFactor* Merge(ImuFactor* f01,  ImuFactor*  f12);
#endif

};
// class ImuFactor

/**
 * ImuFactor2 is a ternary factor that uses NavStates rather than Pose/Velocity.
 * @addtogroup SLAM
 */
class ImuFactor2 : public NoiseModelFactor
{
private:
    PreintegratedImuMeasurements _PIM_;

public:

    /** Default constructor - only use for serialization */
    ImuFactor2() {}

    inline int key1() const
    {
        return keys_[0];
    }
    inline int key2() const
    {
        return keys_[1];
    }
    inline int key3() const
    {
        return keys_[2];
    }

    /**
     * Constructor
     * @param state_i Previous state key
     * @param state_j Current state key
     * @param bias    Previous bias key
     */
    ImuFactor2(int state_i, int state_j, int bias,
               const PreintegratedImuMeasurements& preintegratedMeasurements);

    virtual ~ImuFactor2()
    {
    }

    /// @return a deep copy of this factor
    virtual NonlinearFactor* clone() const;

    /** Access the preintegrated measurements. */

    const PreintegratedImuMeasurements& preintegratedMeasurements() const
    {
        return _PIM_;
    }

    /** implement functions needed to derive from Factor */


    /// vector of errors
    Eigen::VectorXd evaluateError(const NavState& state_i, const NavState& state_j,
                                  const Eigen::VectorXd& bias_i,  //
                                  Eigen::MatrixXd& H1,
                                  Eigen::MatrixXd& H2,
                                  Eigen::MatrixXd& H3) const;

};
// class ImuFactor2
};
#endif
