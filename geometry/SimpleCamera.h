#ifndef SIMPLECAMERA_H
#define SIMPLECAMERA_H
#ifndef DLL_PUBLIC
#define DLL_PUBLIC __attribute__ ((visibility("default")))
#endif

/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 * @file SimpleCamera.h
 * @brief A simple camera class with a Cal3_S2 calibration
 * @date Aug 16, 2009
 * @author Frank Dellaert
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
    static SimpleCamera Level(const Cal3_S2 &K, const Pose2& pose2,
                              double height);

    /// PinholeCamera::level with default calibration
    static SimpleCamera Level(const Pose2& pose2, double height);
    /**
     * Create a camera at the given eye position looking at a target point in the scene
     * with the specified up direction vector.
     * @param eye specifies the camera position
     * @param target the point to look at
     * @param upVector specifies the camera up direction vector,
     *        doesn't need to be on the image plane nor orthogonal to the viewing axis
     * @param K optional calibration parameter
     */
    static SimpleCamera Lookat(const Eigen::Vector3d& eye, const Eigen::Vector3d& target,
                               const Eigen::Vector3d& upVector, const Cal3_S2& K = Cal3_S2());

    /// @}
    /// @name Advanced Constructors
    /// @{

    /// Init from vector, can be 6D (default calibration) or dim
    explicit SimpleCamera(const Eigen::VectorXd &v) :
        PinholeCameraCal3S2(v)
    {
    }

    /// Init from Vector and calibration
    SimpleCamera(const Eigen::VectorXd &v, const Eigen::VectorXd &K) :
        PinholeCameraCal3S2(v, K)
    {
    }

    /// Copy this object as its actual derived type.
    SimpleCamera* clone() const;
    /// @}
    /// @name Manifold
    /// @{

    /// move a cameras according to d
    SimpleCamera retract(const Eigen::VectorXd& d);
    /// @}

};

/// Recover camera from 3*4 camera matrix
SimpleCamera simpleCamera(const Eigen::MatrixXd& P);
};

#endif // SIMPLECAMERA_H
