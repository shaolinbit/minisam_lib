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

namespace minisam
{

/**
 * A pinhole camera class that has a Pose3, functions as base class for all pinhole cameras
 * @addtogroup geometry
 * \nosubgrouping
 */
class PinholeBase:public Pose3
{

public:
    /// @name Standard Constructors
    /// @{

    /** default constructor */
    PinholeBase():Pose3()
    {
    }

    /** constructor with pose */
    explicit PinholeBase(const Pose3& pose) :
        Pose3(pose)
    {
    }

    /// @}
    /// @name Advanced Constructors
    /// @{

    explicit PinholeBase(const minivector &v) :
        Pose3(Pose3::Expmap(v))
    {
    }

    /// @}
    /// @name Standard Interface
    /// @{

    /// return pose, constant version
    const Pose3& pose() const
    {
        return *this;
    }

    /// get rotation
     Rot3 rotation() const
    {
        return this->rotation();
    }

    /// get translation
     minivector translation() const
    {
        return this->translation();
    }

    /// return pose, with derivative
     Pose3 getPose() const;
     Pose3 getPose(minimatrix* H) const;
    /// @}
    /// @name Transformations and measurement functions
    /// @{

    /**
     * Project from 3D point in camera coordinates into image
     * Does *not* throw a CheiralityException, even if pc behind image plane
     * @param pc point in camera coordinates
     */
    static minivector ProjectPoint(const minivector& pc);
    static minivector ProjectPoint(const minivector& pc, //
                                   minimatrix* Dpoint);
    /**
     * Project from 3D point at infinity in camera coordinates into image
     * Does *not* throw a CheiralityException, even if pc behind image plane
     * @param pc point in camera coordinates
     */
    static minivector ProjectUnit(const Unit3& pc);
    static minivector ProjectUnit(Unit3& pc, //
                                  minimatrix* Dpoint);

    /// Project a point into the image and check depth
    std::pair<minivector, bool> projectSafe(const minivector& pw) const;

    /** Project point into the image
     * Throws a CheiralityException if point behind image plane iff THROW_CHEIRALITY_EXCEPTION
     * @param point 3D point in world coordinates
     * @return the intrinsic coordinates of the projected point
     */
    minivector project2Point(const minivector& point) const;
    minivector project2Point(const minivector& point, minimatrix* Dpose,
                             minimatrix* Dpoint) const;

    /** Project point at infinity into the image
     * Throws a CheiralityException if point behind image plane iff THROW_CHEIRALITY_EXCEPTION
     * @param point 3D point in world coordinates
     * @return the intrinsic coordinates of the projected point
     */
    minivector  project2Unit( Unit3& point) const;
    minivector  project2Unit( Unit3& point,
                              minimatrix* Dpose,
                              minimatrix* Dpoint) const;

    /// backproject a 2-dimensional point to a 3-dimensional point at given depth
    static minivector backproject_from_camera(const minivector& p, const double depth);

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

};

 /**
     * Create a level pose at the given 2D pose and height
     * @param K the calibration
     * @param pose2 specifies the location and viewing direction
     * (theta 0 = looking in direction of positive X axis)
     * @param height camera height
     */
   Pose3 PinholeBaseLevelPose(const Pose2& pose2, double height);

    /**
     * Create a camera pose at the given eye position looking at a target point in the scene
     * with the specified up direction vector.
     * @param eye specifies the camera position
     * @param target the point to look at
     * @param upVector specifies the camera up direction vector,
     *        doesn't need to be on the image plane nor orthogonal to the viewing axis
     */
     Pose3 PinholeBaseLookatPose(const minivector& eye, const minivector& target,
                            const minivector& upVector);
    /// @name Derivatives
    /// @{

    /**
     * Calculate Jacobian with respect to pose
     * @param pn projection in normalized coordinates
     * @param d disparity (inverse depth)
     */
     minimatrix PinholeBaseDpose(const minivector& pn, double d);

    /**
     * Calculate Jacobian with respect to point
     * @param pn projection in normalized coordinates
     * @param d disparity (inverse depth)
     * @param Rt transposed rotation matrix
     */
     minimatrix PinholeBaseDpoint(const minivector& pn, double d, const minimatrix& Rt);

    /// @}
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

    /// @name Standard Constructors
    /// @{

    /// default constructor
    CalibratedCamera()
    {
    }

    CalibratedCamera(const minimatrix* m_memory)
    {
    size1=m_memory->size1;
    size2=m_memory->size2;
    prd=m_memory->prd;
    data=m_memory->data;
    owner=0;
    dimension=m_memory->dimension;

    }

    /// construct with pose
    explicit CalibratedCamera(const Pose3& pose) :
        PinholeBase(pose)
    {
    }

    /// @}
    /// @name Named Constructors
    /// @{
    static CalibratedCamera Create(const Pose3& pose)
    {
        return CalibratedCamera(pose);
    }

    static CalibratedCamera Create(const Pose3& pose,
                                   minimatrix* H1)
    {
        minimatrix_set_identity(H1);
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
    static CalibratedCamera Lookat(const minivector& eye, const minivector& target,
                                   const minivector& upVector);

    /// @}
    /// @name Advanced Constructors
    /// @{

    /// construct from vector
    explicit CalibratedCamera(const minivector &v) :
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
    virtual minimatrix* Retract(const minimatrix* d) const;
    virtual minimatrix LocalCoordinates(const minimatrix* T2) const;



    inline int dim() const
    {
        return 6;
    }

    /// @name Transformations and measurement functions
    /// @{

    /**
     * @deprecated
     * Use project2, which is more consistently named across Pinhole cameras
     */
    minivector project(const minivector& point) const;

    minivector project(const minivector& point, minimatrix* Dcamera,
                       minimatrix* Dpoint) const;

    /// backproject a 2-dimensional point to a 3-dimensional point at given depth
    minivector backproject(const minivector& pn, double depth) const;
    /**
     * Calculate range to a landmark
     * @param point 3D location of landmark
     * @return range (double)
     */
    double range(const minivector& point) const;

    double range(const minivector& point,
                 minimatrix* Dcamera,
                 minimatrix* Dpoint) const;

    /**
     * Calculate range to another pose
     * @param pose Other SO(3) pose
     * @return range (double)
     */
    double range(const Pose3& pose) const;
    double range(const Pose3& pose, minimatrix* Dcamera,
                 minimatrix*Dpose) const;

    /**
     * Calculate range to another camera
     * @param camera Other camera
     * @return range (double)
     */
    double range(const CalibratedCamera& camera) const;
    double range(const CalibratedCamera& camera, //
                 minimatrix* H1, //
                 minimatrix* H2) const;
};

};
#endif // CALIBRATEDCAMERA_H
