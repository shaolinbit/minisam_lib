#ifndef SIMPLECAMERA_H
#define SIMPLECAMERA_H



/**
 * @file SimpleCamera.h
 * @brief A simple camera class with a Cal3_S2 calibration
 */

#pragma once
#include "../geometry/PinholeCameraCal3S2.h"
#include "../geometry/Cal3_S2.h"

namespace minisam
{

class SimpleCamera : public PinholeCameraCal3S2
{

public:

    /// @name Standard Constructors
    /// @{

    /** default constructor */
    SimpleCamera() :
        PinholeCameraCal3S2()
    {
    }

    /** constructor with pose */
    explicit SimpleCamera(const Pose3& pose) :
        PinholeCameraCal3S2(pose)
    {
    }

    /** constructor with pose and calibration */
    SimpleCamera(const Pose3& pose, const Cal3_S2& K) :
        PinholeCameraCal3S2(pose, K)
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
    SimpleCamera Level(const Cal3_S2 &K, const Pose2& pose2,
                       double height);

    /// PinholeCamera::level with default calibration
    SimpleCamera Level(const Pose2& pose2, double height);
    /**
     * Create a camera at the given eye position looking at a target point in the scene
     * with the specified up direction vector.
     * @param eye specifies the camera position
     * @param target the point to look at
     * @param upVector specifies the camera up direction vector,
     *        doesn't need to be on the image plane nor orthogonal to the viewing axis
     * @param K optional calibration parameter
     */
    SimpleCamera Lookat(const minivector& eye, const minivector& target,
                        const minivector& upVector, const Cal3_S2& K = Cal3_S2());

    /// @}
    /// @name Advanced Constructors
    /// @{

    /// Init from vector, can be 6D (default calibration) or dim
    explicit SimpleCamera(const minivector &v) :
        PinholeCameraCal3S2(v)
    {
    }

    /// Init from Vector and calibration
    SimpleCamera(const minivector &v, const minivector &K) :
        PinholeCameraCal3S2(v, K)
    {
    }

    /// Copy this object as its actual derived type.
    SimpleCamera* clone() const;
    /// @}
    /// @name Manifold
    /// @{

    /// move a cameras according to d
    virtual minimatrix* Retract(const minimatrix* mpose);
    /// @}

};

/// Recover camera from 3*4 camera matrix
SimpleCamera simpleCamera(const minimatrix& P);
};

#endif // SIMPLECAMERA_H
