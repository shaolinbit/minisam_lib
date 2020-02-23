#ifndef IMUFACTOR_H
#define IMUFACTOR_H

/**
 *  @file  ImuFactor.h
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
class PreintegratedImuMeasurements:public minimatrix// public PreintegrationType
{

    friend class ImuFactor;
    friend class ImuFactor2;
public:

    /// Default constructor for serialization and Cython wrapper
    PreintegratedImuMeasurements();
    ~PreintegratedImuMeasurements() {}

    /**
      *  Constructor, initializes the class with no measurements
      *  @param bias Current estimate of acceleration and rotation rate biases
      *  @param p    Parameters, typically fixed in a single application
      */
    PreintegratedImuMeasurements(const PreintegrationParams& p,
                                 const ConstantBias& biasHat = ConstantBias());

    PreintegratedImuMeasurements(const PreintegratedImuMeasurements& other);
    /**
      *  Construct preintegrated directly from members: base class and preintMeasCov
      *  @param base               PreintegrationType instance
      *  @param preintMeasCov      Covariance matrix used in noise model.
      */
    PreintegratedImuMeasurements(const PreintegrationType& base, const minimatrix& preintMeasCov);


    /**
     * Add a single IMU measurement to the preintegration.
     * @param measuredAcc Measured acceleration (in body frame, as given by the sensor)
     * @param measuredOmega Measured angular velocity (as given by the sensor)
     * @param dt Time interval between this and the last IMU measurement
     */
    void integrateMeasurement(const minivector& measuredAcc,
                              const minivector& measuredOmega, const double dt);

    /// Add multiple measurements, in matrix columns
    void integrateMeasurements(const minimatrix& measuredAccs, const minimatrix& measuredOmegas,
                               const minimatrix& dts);

    /// Return pre-integrated measurement covariance
     minimatrix preintMeasCov() const;

     minimatrix gyroscopeCovariance() const;
     minimatrix  accelerometerCovariance() const;
     minimatrix integrationCovariance() const;



    void print() const;

#ifdef TANGENT_PREINTEGRATION
    /// Merge in a different set of measurements and update bias derivatives accordingly
    void mergeWith(const PreintegratedImuMeasurements& pim, minimatrix* H1, minimatrix* H2);
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
public:
    minimatrix* _PIM_;
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
        delete noiseModel_;
        delete _PIM_;
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
    virtual NoiseModelFactor* clone() const;

    /** Access the preintegrated measurements. */

    const PreintegratedImuMeasurements& preintegratedMeasurements() const
    {
        return PreintegratedImuMeasurements(*_PIM_);
    }

    virtual minivector unwhitenedError(const std::map<int, minimatrix*>& x, std::vector<minimatrix> &H) const;

    virtual minivector unwhitenedError(const std::map<int, minimatrix*>& x) const;
    /** implement functions needed to derive from Factor */

    /// vector of errors
    minivector evaluateError(const minimatrix* pose_i, const minimatrix* vel_i,
                             const minimatrix* pose_j, const minimatrix* vel_j,
                             const minimatrix* bias_i, minimatrix& H1, minimatrix& H2,
                             minimatrix& H3, minimatrix& H4, minimatrix& H5) const;
    minivector evaluateError(const minimatrix* pose_i, const minimatrix* vel_i,
                             const minimatrix* pose_j, const minimatrix* vel_j,
                             const minimatrix* bias_i) const;

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
public:
    minimatrix* _PIM_;//PreintegratedImuMeasurements _PIM_;

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
    virtual NoiseModelFactor* clone() const;

    /** Access the preintegrated measurements. */

    const PreintegratedImuMeasurements& preintegratedMeasurements() const
    {
        return PreintegratedImuMeasurements(*_PIM_);
    }

    /** implement functions needed to derive from Factor */
    virtual minivector unwhitenedError(const std::map<int, minimatrix*>& x, std::vector<minimatrix> &H) const;

    virtual minivector unwhitenedError(const std::map<int, minimatrix*>& x) const;

    /// vector of errors
    minivector evaluateError(const minimatrix* state_i, const minimatrix* state_j,
                             const minimatrix* bias_i) const;
    minivector evaluateError(const minimatrix* state_i, const minimatrix* state_j,
                             const minimatrix* bias_i,  //
                             minimatrix& H1,
                             minimatrix& H2,
                             minimatrix& H3) const;

};
// class ImuFactor2
};
#endif
