#ifndef PINHOLEPOSECAL3S2_H
#define PINHOLEPOSECAL3S2_H

/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 * @file   PinholePose.h
 * @brief  Pinhole camera with known calibration
 * @author Yong-Dian Jian
 * @author Frank Dellaert
 * @date   Feb 20, 2015
 */
#pragma once

#include "../geometry/CalibratedCamera.h"
#include "../geometry/Cal3_S2.h"

namespace minisam
{

/**
 * A pinhole camera class that has a Pose3 and a *fixed* Calibration.
 * @addtogroup geometry
 * \nosubgrouping
 */
class PinholeBaseKCal3S2: public PinholeBase
{

private:

    // Get dimensions of calibration type at compile time
    static const int DimK = Cal3_S2::dimension;

public:

    /// @name Standard Constructors
    /// @{

    /** default constructor */
    PinholeBaseKCal3S2()
    {
    }

    /** constructor with pose */
    explicit PinholeBaseKCal3S2(const Pose3& pose) :
        PinholeBase(pose)
    {
    }

    /// @}
    /// @name Advanced Constructors
    /// @{

    explicit PinholeBaseKCal3S2(const Eigen::VectorXd &v) :
        PinholeBase(v)
    {
    }

    /// @}
    /// @name Standard Interface
    /// @{

    virtual ~PinholeBaseKCal3S2()
    {
    }

    /// return calibration
    virtual  Cal3_S2& calibration()  = 0;

    /// @}
    /// @name Transformations and measurement functions
    /// @{

    /// Project a point into the image and check depth
    std::pair<Eigen::Vector2d, bool> projectSafe(const Eigen::Vector3d& pw)
    {
        std::pair<Eigen::Vector2d, bool> pn = PinholeBase::projectSafe(pw);
        pn.first = calibration().uncalibrate(pn.first);
        return pn;
    }

    /** project a point from world coordinate to the image
     *  @param pw is a point in the world coordinates
     */
    Eigen::Vector2d projectPoint(const Eigen::Vector3d& pw)
    {
        const Eigen::Vector2d pn = PinholeBase::project2Point(pw); // project to normalized coordinates
        return calibration().uncalibrate(pn); // uncalibrate to pixel coordinates
    }

    /** project a point from world coordinate to the image
     *  @param pw is a point at infinity in the world coordinates
     */
    Eigen::Vector2d projectUnit(const Unit3& pw)
    {
        const Unit3 pc = pose().rotation().unrotateUnit(pw); // convert to camera frame
        const Eigen::Vector2d pn = PinholeBase::ProjectUnit(pc); // project to normalized coordinates
        return calibration().uncalibrate(pn);  // uncalibrate to pixel coordinates
    }

    /** Templated projection of a point (possibly at infinity) from world coordinate to the image
     *  @param pw is a 3D point or aUnit3 (point at infinity) in world coordinates
     *  @param Dpose is the Jacobian w.r.t. pose3
     *  @param Dpoint is the Jacobian w.r.t. point3
     *  @param Dcal is the Jacobian w.r.t. calibration
     */
    Eigen::Vector2d _projectPoint(const Eigen::Vector3d& pw,
                                  Eigen::MatrixXd* Dpose,
                                  Eigen::MatrixXd* Dpoint,
                                  Eigen::MatrixXd* Dcal)
    {

        // project to normalized coordinates
        const Eigen::Vector2d pn = PinholeBase::project2Point(pw, Dpose, Dpoint);

        // uncalibrate to pixel coordinates
        Eigen::MatrixXd Dpi_pn(2,2);
        const Eigen::Vector2d pi = calibration().uncalibrate(pn, Dcal,
                                   Dpose || Dpoint ? &Dpi_pn : NULL);

        // If needed, apply chain rule
        if (Dpose!=NULL)
            *Dpose = Dpi_pn * *Dpose;
        if (Dpoint!=NULL)
            *Dpoint = Dpi_pn * *Dpoint;

        return pi;
    }

    Eigen::Vector2d _projectUnit(const Unit3& pw,Eigen::MatrixXd* Dpose,
                                 Eigen::MatrixXd* Dpoint,
                                 Eigen::MatrixXd* Dcal)
    {

        // project to normalized coordinates
        const Eigen::Vector2d pn = PinholeBase::project2Unit(pw, Dpose, Dpoint);

        // uncalibrate to pixel coordinates
        Eigen::MatrixXd Dpi_pn(2,2);
        const Eigen::Vector2d pi = calibration().uncalibrate(pn, Dcal,
                                   Dpose || Dpoint ? &Dpi_pn : NULL);

        // If needed, apply chain rule
        if (Dpose!=NULL)
            *Dpose = Dpi_pn * *Dpose;
        if (Dpoint!=NULL)
            *Dpoint = Dpi_pn * *Dpoint;

        return pi;
    }

    /// project a 3D point from world coordinates into the image
    Eigen::Vector2d projectPoint(const Eigen::Vector3d& pw,Eigen::MatrixXd* Dpose,
                                 Eigen::MatrixXd* Dpoint,
                                 Eigen::MatrixXd* Dcal)
    {
        return _projectPoint(pw, Dpose, Dpoint, Dcal);
    }

    /// project a point at infinity from world coordinates into the image
    Eigen::Vector2d projectUnit(const Unit3& pw, Eigen::MatrixXd* Dpose,
                                Eigen::MatrixXd* Dpoint,
                                Eigen::MatrixXd* Dcal)
    {
        return _projectUnit(pw, Dpose, Dpoint, Dcal);
    }

    /// backproject a 2-dimensional point to a 3-dimensional point at given depth
    Eigen::Vector3d backproject(const Eigen::Vector2d& p, double depth)
    {
        const Eigen::Vector2d pn = calibration().calibrate(p);
        return pose().transform_from(backproject_from_camera(pn, depth));
    }

    /// backproject a 2-dimensional point to a 3-dimensional point at infinity
    Unit3 backprojectPointAtInfinity(const Eigen::Vector2d& p)
    {
        const Eigen::Vector2d pn = calibration().calibrate(p);
        const Unit3 pc(pn(0), pn(1), 1.0); //by convention the last element is 1
        return pose().rotation().rotateUnit(pc);
    }

    /**
     * Calculate range to a landmark
     * @param point 3D location of landmark
     * @return range (double)
     */
    double range(const Eigen::Vector3d& point,
                 Eigen::MatrixXd* Dcamera,
                 Eigen::MatrixXd* Dpoint) const
    {
        return pose().range(point, Dcamera, Dpoint);
    }

    /**
     * Calculate range to another pose
     * @param pose Other SO(3) pose
     * @return range (double)
     */
    double range(const Pose3& pose, Eigen::MatrixXd* Dcamera,
                 Eigen::MatrixXd*  Dpose) const
    {
        return this->pose().range(pose, Dcamera, Dpose);
    }

    /**
     * Calculate range to a CalibratedCamera
     * @param camera Other camera
     * @return range (double)
     */
    double range(const CalibratedCamera& camera, Eigen::MatrixXd*  Dcamera,
                 Eigen::MatrixXd*  Dother) const
    {
        return pose().range(camera.pose(), Dcamera, Dother);
    }

    /**
     * Calculate range to a PinholePoseK derived class
     * @param camera Other camera
     * @return range (double)
     */
    double range(const PinholeBaseKCal3S2& camera,
                 Eigen::MatrixXd* Dcamera,
                 Eigen::MatrixXd* Dother) const
    {
        return pose().range(camera.pose(), Dcamera, Dother);
    }

    ///@}


};
// end of class PinholeBaseK

/**
 * A pinhole camera class that has a Pose3 and a *fixed* Calibration.
 * Instead of using this class, one might consider calibrating the measurements
 * and using CalibratedCamera, which would then be faster.
 * @addtogroup geometry
 * \nosubgrouping
 */
class PinholePoseCal3S2: public PinholeBaseKCal3S2
{

private:

    Cal3_S2* K_; ///< pointer to fixed calibration

public:

    enum
    {
        dimension = 6
    }; ///< There are 6 DOF to optimize for

    /// @name Standard Constructors
    /// @{

    /** default constructor */
    PinholePoseCal3S2()
    {
    }

    /** constructor with pose, uses default calibration */
    explicit PinholePoseCal3S2(const Pose3& pose) :
        PinholeBaseKCal3S2(pose), K_(new Cal3_S2())
    {
    }

    /** constructor with pose and calibration */
    PinholePoseCal3S2(const Pose3& pose, Cal3_S2* K) :
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
    static PinholePoseCal3S2 Level( Cal3_S2* K,
                                    const Pose2& pose2, double height)
    {
        return PinholePoseCal3S2(PinholeBaseKCal3S2::LevelPose(pose2, height), K);
    }

    /// PinholePose::level with default calibration
    static PinholePoseCal3S2 Level(const Pose2& pose2, double height)
    {
        return PinholePoseCal3S2::Level(new Cal3_S2(), pose2, height);
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
    static PinholePoseCal3S2 Lookat(const Eigen::Vector3d& eye, const Eigen::Vector3d& target,
                                    const Eigen::Vector3d& upVector,  Cal3_S2* K =
                                        new Cal3_S2())
    {
        return PinholePoseCal3S2(PinholeBaseKCal3S2::LookatPose(eye, target, upVector), K);
    }

    /// @}
    /// @name Advanced Constructors
    /// @{

    /// Init from 6D vector
    explicit PinholePoseCal3S2(const Eigen::VectorXd  &v) :
        PinholeBaseKCal3S2(v), K_(new Cal3_S2())
    {
    }

    /// Init from Vector and calibration
    PinholePoseCal3S2(const Eigen::VectorXd  &v, const Eigen::VectorXd &K) :
        PinholeBaseKCal3S2(v), K_(new Cal3_S2(K))
    {
    }

    /// stream operator
    friend std::ostream& operator<<(std::ostream &os, const PinholePoseCal3S2& camera)
    {
        os << "{R: " << camera.pose().rotation().rpy().transpose();
        os << ", t: " << camera.pose().translation().transpose();
        if (!camera.K_)
            os << ", K: none";
        else
            os << ", K: " << *camera.K_;
        os << "}";
        return os;
    }
    /// @}
    /// @name Standard Interface
    /// @{

    ~PinholePoseCal3S2()
    {
    }

    /// return  pointer to calibration
    const Cal3_S2* sharedCalibration() const
    {
        return K_;
    }

    /// return calibration
    virtual  Cal3_S2& calibration()
    {
        return *K_;
    }

    /** project a point from world coordinate to the image, 2 derivatives only
     *  @param pw is a point in world coordinates
     *  @param Dpose is the Jacobian w.r.t. the whole camera (really only the pose)
     *  @param Dpoint is the Jacobian w.r.t. point3
     */
    Eigen::Vector2d project2Point(const Eigen::Vector3d& pw, Eigen::MatrixXd* Dpose,
                                  Eigen::MatrixXd* Dpoint)
    {
        return PinholeBaseKCal3S2::projectPoint(pw, Dpose, Dpoint,NULL);
    }

    /// project2 version for point at infinity
    Eigen::Vector2d project2Unit(const Unit3& pw, Eigen::MatrixXd* Dpose,
                                 Eigen::MatrixXd* Dpoint)
    {
        return PinholeBaseKCal3S2::projectUnit(pw, Dpose, Dpoint,NULL);
    }

    /// @}
    /// @name Manifold
    /// @{

    /// @deprecated
    int dim() const
    {
        return 6;
    }

    /// @deprecated
    static int Dim()
    {
        return 6;
    }

    /// move a cameras according to d
    PinholePoseCal3S2 retract(const Eigen::VectorXd& d) const
    {
        return PinholePoseCal3S2(Pose3::ChartAtOrigin::retract(d), K_);
    }

    /// return canonical coordinate
    Eigen::VectorXd localCoordinates(const PinholePoseCal3S2& p) const
    {
        return Pose3::ChartAtOrigin::Local(p.PinholeBaseKCal3S2::pose());
    }

    /// for Canonical
    static PinholePoseCal3S2 identity()
    {
        return PinholePoseCal3S2(); // assumes that the default constructor is valid
    }

    /// @}
};

};
// end of class PinholePose
#endif // PINHOLEPOSE_H
