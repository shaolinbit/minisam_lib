#ifndef  CALIBRATEDCAMERA_H
#define  CALIBRATEDCAMERA_H

/**
 * @file CalibratedCamera.h
 * @brief Calibrated camera for which only pose is unknown
 * @date
 * @author
 */

#pragma once

#include "../geometry/Pose3.h"


/**
 * A pinhole camera class that has a Pose3, functions as base class for all pinhole cameras
 * @addtogroup geometry
 * \nosubgrouping
 */
class PinholeBase
{

public:

    /** Pose Concept requirements */
    //typedef Rot3 Rotation;
// typedef Point3 Translation;

    /**
     *  Some classes template on either PinholeCamera or StereoCamera,
     *  and this typedef informs those classes what "project" returns.
     */
    //typedef Point2 Measurement;

private:

    Pose3 pose_; ///< 3D pose of camera

protected:

    /// @name Derivatives
    /// @{

    /**
     * Calculate Jacobian with respect to pose
     * @param pn projection in normalized coordinates
     * @param d disparity (inverse depth)
     */
    static Eigen::MatrixXd Dpose(const Eigen::Vector2d& pn, double d);

    /**
     * Calculate Jacobian with respect to point
     * @param pn projection in normalized coordinates
     * @param d disparity (inverse depth)
     * @param Rt transposed rotation matrix
     */
    static Eigen::MatrixXd Dpoint(const Eigen::Vector2d& pn, double d, const Eigen::Matrix3d& Rt);

    /// @}

public:

    /// @name Static functions
    /// @{

    /**
     * Create a level pose at the given 2D pose and height
     * @param K the calibration
     * @param pose2 specifies the location and viewing direction
     * (theta 0 = looking in direction of positive X axis)
     * @param height camera height
     */
    static Pose3 LevelPose(const Pose2& pose2, double height);

    /**
     * Create a camera pose at the given eye position looking at a target point in the scene
     * with the specified up direction vector.
     * @param eye specifies the camera position
     * @param target the point to look at
     * @param upVector specifies the camera up direction vector,
     *        doesn't need to be on the image plane nor orthogonal to the viewing axis
     */
    static Pose3 LookatPose(const Eigen::Vector3d& eye, const Eigen::Vector3d& target,
                            const Eigen::Vector3d& upVector);

    /// @}
    /// @name Standard Constructors
    /// @{

    /** default constructor */
    PinholeBase()
    {
    }

    /** constructor with pose */
    explicit PinholeBase(const Pose3& pose) :
        pose_(pose)
    {
    }

    /// @}
    /// @name Advanced Constructors
    /// @{

    explicit PinholeBase(const Eigen::VectorXd &v) :
        pose_(Pose3::Expmap(v))
    {
    }

    /// @}
    /// @name Testable
    /// @{

    /// assert equality up to a tolerance
    //bool equals(const PinholeBase &camera, double tol = 1e-9) const;

    /// print
//void print(const std::string& s = "PinholeBase") const;

    /// @}
    /// @name Standard Interface
    /// @{

    /// return pose, constant version
    const Pose3& pose() const
    {
        return pose_;
    }

    /// get rotation
    const Rot3& rotation() const
    {
        return pose_.rotation();
    }

    /// get translation
    const Eigen::Vector3d& translation() const
    {
        return pose_.translation();
    }

    /// return pose, with derivative
    const Pose3& getPose() const;
    const Pose3& getPose(Eigen::MatrixXd* H) const;
    /// @}
    /// @name Transformations and measurement functions
    /// @{

    /**
     * Project from 3D point in camera coordinates into image
     * Does *not* throw a CheiralityException, even if pc behind image plane
     * @param pc point in camera coordinates
     */
    static Eigen::Vector2d ProjectPoint(const Eigen::Vector3d& pc);
    static Eigen::Vector2d ProjectPoint(const Eigen::Vector3d& pc, //
                                        Eigen::MatrixXd* Dpoint);
    /**
     * Project from 3D point at infinity in camera coordinates into image
     * Does *not* throw a CheiralityException, even if pc behind image plane
     * @param pc point in camera coordinates
     */
    static Eigen::Vector2d ProjectUnit(const Unit3& pc);
    static Eigen::Vector2d ProjectUnit(const Unit3& pc, //
                                       Eigen::Matrix2d* Dpoint);

    /// Project a point into the image and check depth
    std::pair<Eigen::Vector2d, bool> projectSafe(const Eigen::Vector3d& pw) const;

    /** Project point into the image
     * Throws a CheiralityException if point behind image plane iff THROW_CHEIRALITY_EXCEPTION
     * @param point 3D point in world coordinates
     * @return the intrinsic coordinates of the projected point
     */
    Eigen::Vector2d project2Point(const Eigen::Vector3d& point) const;
    Eigen::Vector2d project2Point(const Eigen::Vector3d& point, Eigen::MatrixXd* Dpose,
                                  Eigen::MatrixXd* Dpoint) const;

    /** Project point at infinity into the image
     * Throws a CheiralityException if point behind image plane iff THROW_CHEIRALITY_EXCEPTION
     * @param point 3D point in world coordinates
     * @return the intrinsic coordinates of the projected point
     */
    Eigen::Vector2d  project2Unit(const Unit3& point) const;
    Eigen::Vector2d  project2Unit(const Unit3& point,
                                  Eigen::MatrixXd* Dpose,
                                  Eigen::MatrixXd* Dpoint) const;

    /// backproject a 2-dimensional point to a 3-dimensional point at given depth
    static Eigen::Vector3d backproject_from_camera(const Eigen::Vector2d& p, const double depth);

    /// @}
    /// @name Advanced interface
    /// @{

    /**
     * Return the start and end indices (inclusive) of the translation component of the
     * exponential map parameterization
     * @return a pair of [start, end] indices into the tangent space vector
     */
    inline static std::pair<int, int> translationInterval()
    {
        return std::make_pair(3, 5);
    }

    /// @}

//private:

    /** Serialization function */
    //friend class boost::serialization::access;
    //template<class Archive>
    //void serialize(Archive & ar, const unsigned int /*version*/) {
    //  ar & BOOST_SERIALIZATION_NVP(pose_);
// }

};
// end of class PinholeBase

/**
 * A Calibrated camera class [R|-R't], calibration K=I.
 * If calibration is known, it is more computationally efficient
 * to calibrate the measurements rather than try to predict in pixels.
 * @addtogroup geometry
 * \nosubgrouping
 */
class  CalibratedCamera: public PinholeBase
{

public:

    enum
    {
        dimension = 6
    };

    /// @name Standard Constructors
    /// @{

    /// default constructor
    CalibratedCamera()
    {
    }

    /// construct with pose
    explicit CalibratedCamera(const Pose3& pose) :
        PinholeBase(pose)
    {
    }

    /// @}
    /// @name Named Constructors
    /// @{

    // Create CalibratedCamera, with derivatives
    static CalibratedCamera Create(const Pose3& pose)
    {
        //if (H1)
        //    *H1 << I_6x6;
        return CalibratedCamera(pose);
    }

    static CalibratedCamera Create(const Pose3& pose,
                                   Eigen::MatrixXd* H1)
    {
        // if (H1)
        *H1 << Eigen::MatrixXd::Identity(6,6);
        return CalibratedCamera(pose);
    }
    /**
     * Create a level camera at the given 2D pose and height
     * @param pose2 specifies the location and viewing direction
     * @param height specifies the height of the camera (along the positive Z-axis)
     * (theta 0 = looking in direction of positive X axis)
     */
    static CalibratedCamera Level(const Pose2& pose2, double height);

    /**
     * Create a camera at the given eye position looking at a target point in the scene
     * with the specified up direction vector.
     * @param eye specifies the camera position
     * @param target the point to look at
     * @param upVector specifies the camera up direction vector,
     *        doesn't need to be on the image plane nor orthogonal to the viewing axis
     */
    static CalibratedCamera Lookat(const Eigen::Vector3d& eye, const Eigen::Vector3d& target,
                                   const Eigen::Vector3d& upVector);

    /// @}
    /// @name Advanced Constructors
    /// @{

    /// construct from vector
    explicit CalibratedCamera(const Eigen::VectorXd &v) :
        PinholeBase(v)
    {
    }

    /// @}
    /// @name Standard Interface
    /// @{

    /// destructor
    virtual ~CalibratedCamera()
    {
    }

    /// @}
    /// @name Manifold
    /// @{

    /// move a cameras pose according to d
    CalibratedCamera retract(const Eigen::VectorXd& d) const;

    /// Return canonical coordinate
    Eigen::VectorXd localCoordinates(const CalibratedCamera& T2) const;

    /// @deprecated
    inline int dim() const
    {
        return 6;
    }

    /// @deprecated
    inline static int Dim()
    {
        return 6;
    }

    /// @}
    /// @name Transformations and measurement functions
    /// @{

    /**
     * @deprecated
     * Use project2, which is more consistently named across Pinhole cameras
     */
    Eigen::Vector2d project(const Eigen::Vector3d& point) const;

    Eigen::Vector2d project(const Eigen::Vector3d& point, Eigen::MatrixXd* Dcamera,
                            Eigen::MatrixXd* Dpoint) const;

    /// backproject a 2-dimensional point to a 3-dimensional point at given depth
    Eigen::Vector3d backproject(const Eigen::Vector2d& pn, double depth) const
    {
        return pose().transform_from(backproject_from_camera(pn, depth));
    }

    /**
     * Calculate range to a landmark
     * @param point 3D location of landmark
     * @return range (double)
     */
    double range(const Eigen::Vector3d& point) const
    {
        return pose().range(point);
    }

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
    double range(const Pose3& pose) const
    {
        return this->pose().range(pose);
    }
    double range(const Pose3& pose, Eigen::MatrixXd* Dcamera,
                 Eigen::MatrixXd*Dpose) const
    {
        return this->pose().range(pose, Dcamera, Dpose);
    }

    /**
     * Calculate range to another camera
     * @param camera Other camera
     * @return range (double)
     */
    double range(const CalibratedCamera& camera) const
    {
        return pose().range(camera.pose());
    }
    double range(const CalibratedCamera& camera, //
                 Eigen::MatrixXd* H1, //
                 Eigen::MatrixXd* H2) const
    {
        return pose().range(camera.pose(), H1, H2);
    }
};


#endif // CALIBRATEDCAMERA_H
