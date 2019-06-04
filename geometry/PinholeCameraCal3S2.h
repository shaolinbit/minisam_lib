#ifndef PINHOLECAMERACAL3S2_H
#define PINHOLECAMERACAL3S2_H

/**
 * @file PinholeCameraCal3S2.h
 * @brief Base class for all pinhole cameras
 * @author
 * @date
 */

#pragma once

#include "../geometry/PinholePoseCal3S2.h"

/**
 * A pinhole camera class that has a Pose3 and a Calibration.
 * Use PinholePose if you will not be optimizing for Calibration
 * @addtogroup geometry
 * \nosubgrouping
 */
class PinholeCameraCal3S2: public PinholeBaseKCal3S2
{
private:

    //typedef PinholeBaseK<Calibration> Base; ///< base class has 3D pose as private member
    Cal3_S2 K_; ///< Calibration, part of class now

    // Get dimensions of calibration type at compile time
    static const int DimK = Cal3_S2::dimension;

public:

    enum
    {
        dimension = 6 + DimK
    }; ///< Dimension depends on calibration

    /// @name Standard Constructors
    /// @{

    /** default constructor */
    PinholeCameraCal3S2()
    {
    }

    /** constructor with pose */
    explicit PinholeCameraCal3S2(const Pose3& pose) :
        PinholeBaseKCal3S2(pose)
    {
    }

    /** constructor with pose and calibration */
    PinholeCameraCal3S2(const Pose3& pose, const Cal3_S2& K) :
        PinholeBaseKCal3S2(pose), K_(K)
    {
    }

    /// @}
    /// @name Named Constructors
    /// @{

    /**
     * Create a level camera at the given 2D pose and height
     * @param K the calibration
     * @param pose2 specifies the location and viewing direction
     * (theta 0 = looking in direction of positive X axis)
     * @param height camera height
     */
    static PinholeCameraCal3S2 Level(const Cal3_S2 &K, const Pose2& pose2,
                                     double height)
    {
        return PinholeCameraCal3S2(PinholeBaseKCal3S2::LevelPose(pose2, height), K);
    }

    /// PinholeCamera::level with default calibration
    static PinholeCameraCal3S2 Level(const Pose2& pose2, double height)
    {
        return PinholeCameraCal3S2::Level(Cal3_S2(), pose2, height);
    }

    /**
     * Create a camera at the given eye position looking at a target point in the scene
     * with the specified up direction vector.
     * @param eye specifies the camera position
     * @param target the point to look at
     * @param upVector specifies the camera up direction vector,
     *        doesn't need to be on the image plane nor orthogonal to the viewing axis
     * @param K optional calibration parameter
     */
    static PinholeCameraCal3S2 Lookat(const Eigen::Vector3d& eye, const Eigen::Vector3d& target,
                                      const Eigen::Vector3d& upVector, const Cal3_S2& K = Cal3_S2())
    {
        return PinholeCameraCal3S2(PinholeBaseKCal3S2::LookatPose(eye, target, upVector), K);
    }

    // Create PinholeCamera, with derivatives
    static PinholeCameraCal3S2 Create(const Pose3& pose, const Cal3_S2 &K,
                                      Eigen::MatrixXd* H1, //
                                      Eigen::MatrixXd* H2)
    {
        // typedef Eigen::Matrix<double, DimK, 6> MatrixK6;
        if (H1!=NULL)
        {
            *H1<<Eigen::MatrixXd::Identity(6,6),Eigen::MatrixXd::Zero(DimK,6);
        }
        //   *H1 << I_6x6, MatrixK6::Zero();
        // typedef Eigen::Matrix<double, 6, DimK> Matrix6K;
        // typedef Eigen::Matrix<double, DimK, DimK> MatrixK;
        if (H2!=NULL)
            *H2 << Eigen::MatrixXd::Zero(6,DimK),Eigen::MatrixXd::Identity(DimK,DimK);
        return PinholeCameraCal3S2(pose,K);
    }

    /// @}
    /// @name Advanced Constructors
    /// @{

    /// Init from vector, can be 6D (default calibration) or dim
    explicit PinholeCameraCal3S2(const Eigen::VectorXd &v) :
        PinholeBaseKCal3S2(v.head<6>())
    {
        if (v.size() > 6)
            K_ = Cal3_S2(v.tail<DimK>());
    }

    /// Init from Vector and calibration
    PinholeCameraCal3S2(const Eigen::VectorXd &v, const Eigen::VectorXd &K) :
        PinholeBaseKCal3S2(v), K_(K)
    {
    }

    /// @}
    /// @name Standard Interface
    /// @{

    ~PinholeCameraCal3S2()
    {
    }

    /// return pose
    const Pose3& pose() const
    {
        return PinholeBaseKCal3S2::pose();
    }

    /// return pose, with derivative
    const Pose3& getPose(Eigen::MatrixXd* H) const
    {
        if (H!=NULL)
        {
            H->setZero();
            H->template block<6, 6>(0, 0) = Eigen::MatrixXd::Identity(6,6);
        }
        return PinholeBaseKCal3S2::pose();
    }

    /// return calibration
    Cal3_S2& calibration()
    {
        return K_;
    }

    /// @}
    /// @name Manifold
    /// @{

    /// @deprecated
    int dim() const
    {
        return dimension;
    }

    /// @deprecated
    static int Dim()
    {
        return dimension;
    }

// typedef Eigen::Matrix<double, dimension, 1> VectorK6;

    /// move a cameras according to d
    PinholeCameraCal3S2 retract(const Eigen::VectorXd& d)
    {
        if ( d.size() == 6)
        {
            // return PinholeCamera(this->pose().retract(d), calibration());
            return PinholeCameraCal3S2(Pose3::ChartAtOrigin::Retract(d), calibration());
        }
        else
        {
            //return PinholeCameraCal3S2(this->pose().retract(d.head<6>()),
            //calibration().retract(d.tail(calibration().dim())));
            return PinholeCameraCal3S2(Pose3::ChartAtOrigin::Retract(d.head<6>()),
                                       calibration().retract(d.tail(calibration().dim())));

        }

    }

    /// return canonical coordinate
    Eigen::VectorXd localCoordinates(PinholeCameraCal3S2& T2)
    {
        Eigen::VectorXd d;
        // d.template head<6>() = this->pose().localCoordinates(T2.pose());
        d.template head<6>() = Pose3::ChartAtOrigin::Local(T2.pose());
        d.template tail<DimK>() = calibration().localCoordinates(T2.calibration());
        return d;
    }

    /// for Canonical
    static PinholeCameraCal3S2 identity()
    {
        return PinholeCameraCal3S2(); // assumes that the default constructor is valid
    }

    /// @}
    /// @name Transformations and measurement functions
    /// @{

    //typedef Eigen::Matrix<double, 2, DimK> Matrix2K;

    /** Templated projection of a 3D point or a point at infinity into the image
     *  @param pw either a Point3 or a Unit3, in world coordinates
     */
    //template<class POINT>
    Eigen::Vector2d _project2Point(const Eigen::Vector3d& pw, Eigen::MatrixXd*  Dcamera,
                                   Eigen::MatrixXd* Dpoint)
    {
        // We just call 3-derivative version in Base
        Eigen::MatrixXd Dpose(2,6);
        Eigen::MatrixXd Dcal(2,DimK);
        // Eigen::Matrix<double, 2, DimK> Dcal;
        Eigen::Vector2d pi = PinholeBaseKCal3S2::projectPoint(pw, Dcamera ? &Dpose : NULL, Dpoint,
                             Dcamera ? &Dcal : NULL);
        if (Dcamera!=NULL)
            *Dcamera << Dpose, Dcal;
        return pi;
    }
    // template<class POINT>
    Eigen::Vector2d _project2Unit(const Unit3& pw, Eigen::MatrixXd*  Dcamera,
                                  Eigen::MatrixXd* Dpoint)
    {
        // We just call 3-derivative version in Base
        Eigen::MatrixXd Dpose(2,6);
        Eigen::MatrixXd Dcal(2,DimK);
        //Eigen::Matrix<double, 2, DimK> Dcal;
        Eigen::Vector2d pi = PinholeBaseKCal3S2::projectUnit(pw, Dcamera ? &Dpose : NULL, Dpoint,
                             Dcamera ? &Dcal : NULL);
        if (Dcamera!=NULL)
            *Dcamera << Dpose, Dcal;
        return pi;
    }

    /// project a 3D point from world coordinates into the image
    Eigen::Vector2d project2Point(const Eigen::Vector3d& pw,  Eigen::MatrixXd*  Dcamera,
                                  Eigen::MatrixXd*  Dpoint)
    {
        return _project2Point(pw, Dcamera, Dpoint);
    }

    /// project a point at infinity from world coordinates into the image
    Eigen::Vector2d project2Unit(const Unit3& pw, Eigen::MatrixXd*  Dcamera,
                                 Eigen::MatrixXd*  Dpoint)
    {
        return _project2Unit(pw, Dcamera, Dpoint);
    }

    /**
     * Calculate range to a landmark
     * @param point 3D location of landmark
     * @return range (double)
     */
    double range(const Eigen::Vector3d& point, Eigen::MatrixXd* Dcamera,
                 Eigen::MatrixXd* Dpoint) const
    {
        Eigen::MatrixXd Dpose_(1,6);
        double result = this->pose().range(point, Dcamera ? &Dpose_ : NULL, Dpoint);
        if (Dcamera!=NULL)
            *Dcamera << Dpose_, Eigen::Matrix<double, 1, DimK>::Zero();
        return result;
    }

    /**
     * Calculate range to another pose
     * @param pose Other SO(3) pose
     * @return range (double)
     */
    double range(const Pose3& pose, Eigen::MatrixXd* Dcamera,
                 Eigen::MatrixXd* Dpose) const
    {
        Eigen::MatrixXd Dpose_(1,6);
        double result = this->pose().range(pose, Dcamera ? &Dpose_ : NULL, Dpose);
        if (Dcamera!=NULL)
            *Dcamera << Dpose_, Eigen::MatrixXd::Zero(1,DimK);
        return result;
    }

    /**
     * Calculate range to another camera
     * @param camera Other camera
     * @return range (double)
     */
// template<class CalibrationB>
    double range(const PinholeCameraCal3S2& camera,
                 Eigen::MatrixXd*  Dcamera,
                 Eigen::MatrixXd*  Dother) const
    {
        Eigen::MatrixXd Dcamera_(1,6), Dother_(1,6);
        double result = this->pose().range(camera.pose(), Dcamera ? &Dcamera_ : NULL,
                                           Dother ? &Dother_ : NULL);
        if (Dcamera!=NULL)
        {
            *Dcamera << Dcamera_, Eigen::MatrixXd::Zero(1,DimK);
        }
        if (Dother!=NULL)
        {
            Dother->setZero();
            Dother->template block<1, 6>(0, 0) = Dother_;
        }
        return result;
    }

    /**
     * Calculate range to a calibrated camera
     * @param camera Other camera
     * @return range (double)
     */
    double range(const CalibratedCamera& camera,
                 Eigen::MatrixXd*  Dcamera,
                 Eigen::MatrixXd*   Dother) const
    {
        return range(camera.pose(), Dcamera, Dother);
    }

};

#endif // PINHOLECAMERA_H
